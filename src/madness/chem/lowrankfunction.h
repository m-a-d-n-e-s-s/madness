//
// Created by Florian Bischoff on 8/10/23.
//

#ifndef MADNESS_LOWRANKFUNCTION_H
#define MADNESS_LOWRANKFUNCTION_H


#include<madness/mra/mra.h>
#include<madness/mra/vmra.h>
#include<madness/world/timing_utilities.h>
#include<madness/chem/electronic_correlation_factor.h>
#include<madness/chem/IntegratorXX.h>
#include <random>
#include <algorithm>



namespace madness {

    struct LowRankFunctionParameters : QCCalculationParametersBase {

        LowRankFunctionParameters() : QCCalculationParametersBase() {

            // initialize with: key, value, comment (optional), allowed values (optional)
            initialize<double>("radius",2.0,"the radius");
            initialize<double>("gamma",1.0,"the exponent of the correlation factor");
            initialize<double>("volume_element",0.1,"volume covered by each grid point");
            initialize<double>("tol",1.e-8,"rank-reduced cholesky tolerance");
            initialize<std::string>("f12type","Slater","correlation factor",{"Slater","SlaterF12"});
            initialize<std::string>("orthomethod","cholesky","orthonormalization",{"cholesky","canonical","symmetric"});
            initialize<std::string>("transpose","slater2","transpose of the matrix",{"slater1","slater2"});
            initialize<std::string>("gridtype","random","the grid type",{"random","adaptive","twostage","harmonics"});
            initialize<int>("optimize",1,"number of optimization iterations");
            initialize<int>("lmax",2,"max angular momentum for the RI harmonics");
            initialize<bool>("canonicalize",false,"canonicalize the rep, i.e. metric is the identity");
            initialize<std::vector<double>>("tempered",{0.05,2.0,3.0},"zeta_min,zeta_max,factor");
            initialize<double>("adaptive_coarse_factor",8.0,"coarse grid uses volume_element*factor");
            initialize<double>("adaptive_refine_radius",0.5,"local refinement radius as fraction of radius");
            initialize<double>("adaptive_significance_ratio",0.2,"relative threshold for significant coarse Y norms");
            initialize<int>("adaptive_max_centers",16,"max number of significant coarse centers to refine");
            initialize<int>("adaptive_min_centers",2,"minimum number of coarse centers to keep");
        }
        [[nodiscard]] std::string get_tag() const override {
            return {"lrf"};
        }


        void read_and_set_derived_values(World& world, const commandlineparser& parser, std::string tag) {
            read_input_and_commandline_options(world,parser,tag);
        }

        [[nodiscard]] double radius() const {return get<double>("radius");}
        [[nodiscard]] double gamma() const {return get<double>("gamma");}
        [[nodiscard]] double volume_element() const {return get<double>("volume_element");}
        [[nodiscard]] double tol() const {return get<double>("tol");}
        [[nodiscard]] int optimize() const {return get<int>("optimize");}
        [[nodiscard]] int lmax() const {return get<int>("lmax");}
        [[nodiscard]] bool canonicalize() const {return get<bool>("canonicalize");}
        [[nodiscard]] std::string gridtype() const {return get<std::string>("gridtype");}
        [[nodiscard]] std::string orthomethod() const {return get<std::string>("orthomethod");}
        [[nodiscard]] std::string f12type() const {return get<std::string>("f12type");}
        [[nodiscard]] std::vector<double> tempered() const {return get<std::vector<double>>("tempered");}
        [[nodiscard]] double adaptive_coarse_factor() const {return get<double>("adaptive_coarse_factor");}
        [[nodiscard]] double adaptive_refine_radius() const {return get<double>("adaptive_refine_radius");}
        [[nodiscard]] double adaptive_significance_ratio() const {return get<double>("adaptive_significance_ratio");}
        [[nodiscard]] int adaptive_max_centers() const {return get<int>("adaptive_max_centers");}
        [[nodiscard]] int adaptive_min_centers() const {return get<int>("adaptive_min_centers");}
    };


    class gridbase {
    public:
        [[nodiscard]] double get_volume_element() const {return volume_element;}
        [[nodiscard]] double get_radius() const {return radius;}

        virtual ~gridbase() = default;
        // visualize the grid in xyz format
        template<std::size_t NDIM>
        void visualize(const std::string filename, const std::vector<Vector<double,NDIM>>& grid) const {
            print("visualizing grid to file",filename);
            print("a total of",grid.size(),"grid points");
            std::ofstream file(filename);
            for (const auto& r : grid) {
                // formatted output
                file << std::fixed << std::setprecision(6);
                for (int i=0; i<NDIM; ++i) file << r[i] << " ";
                file << std::endl;
            }
            file.close();
        }

    protected:
        double volume_element=0.1;
        double radius=3;
        bool do_print=false;
    };

    template<std::size_t NDIM>
    class dftgrid : public gridbase {
    public:
        GridBuilder builder;
        explicit dftgrid(const double volume_element, const double radius) {
            // increase number of radial grid points until the volume element is below the threshold
            double current_ve=1.0;
            std::size_t nradial=10;
            while (current_ve>volume_element) {
                nradial+=10;
                print("trying nradial",nradial);
                GridBuilder tmp;
                tmp.set_angular_order(7);
                tmp.set_nradial(nradial);
                tmp.make_grid();
                double volume=4./3. *M_PI * std::pow(radius,3.0);
                auto npoints=tmp.get_points().size();
                current_ve=volume/npoints;
                print("volume, npoints, volume element",volume,npoints,current_ve);
            }
            builder.set_nradial(nradial);
            builder.set_angular_order(7);
        }


        explicit dftgrid(const std::size_t nradial, const std::size_t angular_order, const coord_3d origin=coord_3d()) {
            static_assert(NDIM==3,"DFT Grids only defined for NDIM=3");
            builder.set_nradial(nradial);
            builder.set_angular_order(angular_order);
            builder.set_origin(origin);
            builder.make_grid();
        }

        std::vector<Vector<double,NDIM>> get_grid() const {
            return builder.get_points();
        }

    };

    /// grid with random points around the origin, with a Gaussian distribution
    template<std::size_t NDIM>
    class randomgrid : public gridbase {
    public:
        randomgrid(const double volume_element, const double radius, const Vector<double,NDIM> origin=Vector<double,NDIM>(0.0))
            : gridbase(), origin(origin) {
            this->volume_element=volume_element;
            this->radius=radius;
        }

        std::vector<Vector<double,NDIM>> get_grid() const {
            std::vector<Vector<double, NDIM>> grid;
            long npoint_within_volume=volume()/volume_element;
            if (do_print) print("npoint_within_volume",npoint_within_volume);

            auto cell = FunctionDefaults<NDIM>::get_cell();
            auto is_in_cell = [&cell](const Vector<double, NDIM>& r) {
                for (size_t d = 0; d < NDIM; ++d) if (r[d] < cell(d, 0) or r[d] > cell(d, 1)) return false;
                return true;
            };
            double rad=radius;
            auto o=origin;
            auto is_in_sphere = [&rad,&o](const Vector<double, NDIM>& r) {
                return ((r-o).normf()<rad);
            };

            // set variance such that about 70% of all grid points sits within the radius
            double variance=radius;
            if (NDIM==2) variance=0.6*radius;
            if (NDIM==3) variance=0.5*radius;
            long maxrank=10*npoint_within_volume;
            long rank=0;
            for (int r=0; r<maxrank; ++r) {
                auto tmp = gaussian_random_distribution(origin, variance);
                if (not is_in_cell(tmp)) continue;
                if (is_in_sphere(tmp)) ++rank;
                grid.push_back(tmp);
                if (rank==npoint_within_volume) break;
            }
            if (do_print) {
                print("origin                ",origin);
                print("radius                ",radius);
                print("grid points in volume ",rank);
                print("total grid points     ",grid.size());
                print("ratio                 ",rank/double(grid.size()));
                print("volume element        ",volume()/rank);
            }
            return grid;
        }

        Vector<double,NDIM> get_origin() const {
            return origin;
        }

    private:

        double volume() const {
            MADNESS_CHECK(NDIM>0 and NDIM<4);
            if (NDIM==1) return 2.0*radius;
            else if (NDIM==2) return constants::pi*radius*radius;
            else if (NDIM==3) return 4.0 / 3.0 * constants::pi * std::pow(radius, 3.0);
            else {
                MADNESS_EXCEPTION("invalid NDIM",1);
                return 0.0;
            }
        }

        static Vector<double,NDIM> gaussian_random_distribution(const Vector<double,NDIM> origin, double variance) {
            std::random_device rd{};
            std::mt19937 gen{rd()};
            Vector<double,NDIM> result;
            for (size_t i = 0; i < NDIM; ++i) {
                std::normal_distribution<> d{origin[i], variance};
                result[i]=d(gen);
            }

            return result;
        }

        Vector<double,NDIM> origin;

    };

    template<std::size_t NDIM>
    struct cartesian_grid {
        Vector<double,NDIM> lovec,hivec;
        std::vector<long> stride;
        long index=0;
        long n_per_dim=0;
        long total_n=0;
        Vector<double,NDIM> increment;

        cartesian_grid(const double volume_per_gridpoint, const double radius) {
            double length_per_gridpoint=std::pow(volume_per_gridpoint,1.0/NDIM);
            n_per_dim=ceil(2.0*radius/length_per_gridpoint);
            print("length per gridpoint, n_per_dim",length_per_gridpoint,n_per_dim);
            print("volume_per_gridpoint",std::pow(length_per_gridpoint,NDIM));
            initialize(-radius,radius);
            print("increment",increment);
        }


        cartesian_grid(const long n_per_dim, const double lo, const double hi)
                : n_per_dim(n_per_dim) {
            initialize(lo,hi);
        }

        cartesian_grid(const cartesian_grid<NDIM>& other) : lovec(other.lovec),
                                                            hivec(other.hivec), stride(other.stride), index(0), n_per_dim(other.n_per_dim),
                                                            total_n(other.total_n), increment(other.increment) {
        }

        cartesian_grid& operator=(const cartesian_grid<NDIM>& other) {
            cartesian_grid<NDIM> tmp(other);
            std::swap(*this,other);
            return *this;
        }

        void initialize(const double lo, const double hi) {
            lovec.fill(lo);
            hivec.fill(hi);
            increment=(hivec-lovec)*(1.0/double(n_per_dim-1));
            stride=std::vector<long>(NDIM,1l);
            total_n=std::pow(n_per_dim,NDIM);
            for (long i=NDIM-2; i>=0; --i) stride[i]=n_per_dim*stride[i+1];

        }

        double volume_per_gridpoint() const{
            double volume=1.0;
            for (int i=0; i<NDIM; ++i) volume*=(hivec[i]-lovec[i]);
            return volume/total_n;
        }

        void operator++() {
            index++;
        }

        bool operator()() const {
            return index < total_n;
        }

        Vector<double,NDIM> get_coordinates() const {
            Vector<double,NDIM> tmp(NDIM);
            for (int idim=0; idim<NDIM; ++idim) {
                tmp[idim]=(index/stride[idim])%n_per_dim;
            }
            return lovec+tmp*increment;
        }

    };

    /// given a molecule, return a suitable grid
    template<std::size_t NDIM>
    class molecular_grid : public gridbase {

    public:
        /// ctor takes centers of the grid and the grid parameters
        molecular_grid(const std::vector<Vector<double,NDIM>> origins, const LowRankFunctionParameters& params)
            : centers(origins)
        {
            if (centers.size()==0) centers.push_back(Vector<double,NDIM>(0) );
            // treat harmonics and twostage like random sampling for constructing the molecular grid
            if (params.gridtype()=="random" or params.gridtype()=="adaptive" or params.gridtype()=="harmonics" or params.gridtype()=="twostage") grid_builder=std::make_shared<randomgrid<NDIM>>(params.volume_element(),params.radius());
            // else if (params.gridtype()=="cartesian") grid_builder=std::make_shared<cartesian_grid<NDIM>>(params.volume_element(),params.radius());
            else if (params.gridtype()=="dftgrid") {
                if constexpr (NDIM==3) {
                    grid_builder=std::make_shared<dftgrid<NDIM>>(params.volume_element(),params.radius());
                } else {
                    MADNESS_EXCEPTION("no dft grid with NDIM != 3",1);
                }
            } else {
                MADNESS_EXCEPTION("no such grid type",1);
            }
        }


        /// ctor takes centers of the grid and the grid builder
        molecular_grid(const std::vector<Vector<double,NDIM>> origins, std::shared_ptr<gridbase> grid)
            : centers(origins), grid_builder(grid) {
            if (centers.size()==0) centers.push_back({0,0,0});
        }

        /// ctor takes molecule and grid builder
        molecular_grid(const Molecule& molecule, std::shared_ptr<gridbase> grid)
        : molecular_grid(molecule.get_all_coords_vec(),grid) {}

        std::vector<Vector<double,NDIM>> get_grid() const {
            MADNESS_CHECK_THROW(grid_builder,"no grid builder given in molecular_grid");
            MADNESS_CHECK_THROW(centers.size()>0,"no centers given in molecular_grid");
            std::vector<Vector<double,NDIM>> grid;
            for (const auto& coords : centers) {
                print("atom sites",coords);
                if (auto builder=dynamic_cast<dftgrid<NDIM>*>(grid_builder.get())) {
                    if constexpr (NDIM==3) {
                        dftgrid<NDIM> b1(builder->builder.get_nradial(),builder->builder.get_angular_order(),coords);
                        auto atomgrid=b1.get_grid();
                        grid.insert(grid.end(),atomgrid.begin(),atomgrid.end());
                    } else {
                        MADNESS_EXCEPTION("no DFT grid for NDIM /= 3",1);
                    }
                } else if (auto builder=dynamic_cast<randomgrid<NDIM>*>(grid_builder.get())) {
                    randomgrid<NDIM> b1(builder->get_volume_element(),builder->get_radius(),coords);
                    auto atomgrid=b1.get_grid();
                    grid.insert(grid.end(),atomgrid.begin(),atomgrid.end());
                } else {
                    MADNESS_EXCEPTION("no such grid builder",1);
                }
            }
            return grid;
        }

    private:
        std::vector<Vector<double,NDIM>> centers;
        std::shared_ptr<gridbase> grid_builder;

    };

template<std::size_t PDIM>
struct particle {
    std::array<int,PDIM> dims;

    /// default constructor
    particle() = default;

    /// convenience for particle 1 (the left/first particle)
    static particle particle1() {
        particle p;
        for (size_t i=0; i<PDIM; ++i) p.dims[i]=i;
        return p;
    }
    /// convenience for particle 2 (the right/second particle)
    static particle particle2() {
        particle p;
        for (size_t i=0; i<PDIM; ++i) p.dims[i]=i+PDIM;
        return p;
    }

    /// return the other particle
    particle complement() const {
        MADNESS_CHECK(is_first() or is_last());
        if (is_first()) return particle2();
        return particle1();
    }

    explicit particle(const int p) : particle(std::vector<int>(1,p)) {}
    explicit particle(const int p1, const int p2) : particle(std::vector<int>({p1,p2})) {}
    explicit particle(const int p1, const int p2,const int p3) : particle(std::vector<int>({p1,p2,p3})) {}
    explicit particle(const std::vector<int> p) {
        for (int i=0; i<PDIM; ++i) dims[i]=p[i];
    }

    [[nodiscard]] std::string str() const {
        std::stringstream ss;
        ss << *this;
        return ss.str();
    }


    /// type conversion to std::array
    std::array<int,PDIM> get_array() const {
        return dims;
    }


    /// assuming two particles only
    bool is_first() const {return dims[0]==0;}
    /// assuming two particles only
    bool is_last() const {return dims[0]==(PDIM);}

    template<std::size_t DUMMYDIM=PDIM>
    typename std::enable_if_t<DUMMYDIM==1, std::tuple<int>>
    get_tuple() const {return std::tuple<int>(dims[0]);}

    template<std::size_t DUMMYDIM=PDIM>
    typename std::enable_if_t<DUMMYDIM==2, std::tuple<int,int>>
    get_tuple() const {return std::tuple<int,int>(dims[0],dims[1]);}

    template<std::size_t DUMMYDIM=PDIM>
    typename std::enable_if_t<DUMMYDIM==3, std::tuple<int,int,int>>
    get_tuple() const {return std::tuple<int,int,int>(dims[0],dims[1],dims[2]);}
};

template<std::size_t PDIM>
std::ostream& operator<<(std::ostream& os, const particle<PDIM>& p) {
    os << "(";
    for (auto i=0; i<PDIM-1; ++i) os << p.dims[i] << ";";
    os << p.dims[PDIM-1] << ")";
    return os;
}

/// the low-rank functor is what the LowRankFunction will represent

/// derive from this class :
/// must implement in inner product
/// may implement an operator()(const coord_nd&)
template<typename T, std::size_t NDIM, std::size_t LDIM=NDIM/2>
struct LRFunctorBase : public FunctionFunctorInterface<T,NDIM> {

    virtual ~LRFunctorBase() {};
    virtual std::vector<Function<T,LDIM>> inner(const std::vector<Function<T,LDIM>>& rhs,
                                        const particle<LDIM> p1, const particle<LDIM> p2) const =0;

    virtual Function<T,LDIM> inner(const Function<T,LDIM>& rhs, const particle<LDIM> p1, const particle<LDIM> p2) const {
        return inner(std::vector<Function<T,LDIM>>({rhs}),p1,p2)[0];
    }

    /// evaluate the functor at a given point, e.g. for plotting
    virtual T operator()(const Vector<T,NDIM>& r) const =0;

    virtual typename Tensor<T>::scalar_type norm2() const {
        MADNESS_EXCEPTION("L2 norm not implemented",1);
    }

    /// introspection
    virtual std::string type() const = 0;

    virtual World& world() const =0;
    friend std::vector<Function<T,LDIM>> inner(const LRFunctorBase& functor, const std::vector<Function<T,LDIM>>& rhs,
                                               const particle<LDIM> p1, const particle<LDIM> p2) {
        return functor.inner(rhs,p1,p2);
    }
    friend Function<T,LDIM> inner(const LRFunctorBase& functor, const Function<T,LDIM>& rhs,
                                               const particle<LDIM> p1, const particle<LDIM> p2) {
        return functor.inner(rhs,p1,p2);
    }

};

template<typename T, std::size_t NDIM, std::size_t LDIM=NDIM/2>
struct LRFunctorF12 : public LRFunctorBase<T,NDIM> {
//    LRFunctorF12() = default;
    LRFunctorF12(const std::shared_ptr<SeparatedConvolution<T,LDIM>> f12, const std::vector<Function<T,LDIM>>& a,
                 const std::vector<Function<T,LDIM>>& b) : f12(f12), a(a), b(b) {

        // if a or b are missing, they are assumed to be 1
        // you may not provide a or b, but if you do they have to have the same size because they are summed up
        if (a.size()>0 and b.size()>0)
            MADNESS_CHECK_THROW(a.size()==b.size(), "a and b must have the same size");
        if (a.size()==0) this->a.resize(b.size());
        if (b.size()==0) this->b.resize(a.size());
        MADNESS_CHECK_THROW(this->a.size()==this->b.size(), "a and b must have the same size");
    }

    /// delegate to the other ctor with vector arguments
    LRFunctorF12(const std::shared_ptr<SeparatedConvolution<T,LDIM>> f12,
                 const Function<T,LDIM>& a, const Function<T,LDIM>& b)
                 : LRFunctorF12(f12,std::vector<Function<T,LDIM>>({a}), std::vector<Function<T,LDIM>>({b})) {}


    std::string type() const override {return "LRFunctorF12";}
private:
    std::shared_ptr<SeparatedConvolution<T,LDIM>> f12;  ///< a two-particle function
    std::vector<Function<T,LDIM>> a,b;   ///< the lo-dim functions
public:

    World& world() const override {return f12->get_world();}
    std::vector<Function<T,LDIM>> inner(const std::vector<Function<T,LDIM>>& rhs,
                                        const particle<LDIM> p1, const particle<LDIM> p2) const override {

        std::vector<Function<T,LDIM>> result;
        // functor is now \sum_i a_i(1) b_i(2) f12
        // result(1) = \sum_i \int a_i(1) f(1,2) b_i(2) rhs(2) d2
        //            = \sum_i a_i(1) \int f(1,2) b_i(2) rhs(2) d2
        World& world=rhs.front().world();

        const int nbatch=30;
        for (size_t i=0; i<rhs.size(); i+=nbatch) {
            std::vector<Function<T,LDIM>> rhs_batch;
            auto begin= rhs.begin()+i;
            auto end= size_t(i+nbatch)<rhs.size() ? rhs.begin()+i+nbatch : rhs.end();
            std::copy(begin,end, std::back_inserter(rhs_batch));
            auto tmp2= zero_functions_compressed<T,LDIM>(world,rhs_batch.size());

            if (a.size()==0) tmp2=apply(world,*(f12),rhs_batch);

            for (size_t ia=0; ia<a.size(); ia++) {
                auto premultiply= p1.is_first() ? a[ia] : b[ia];
                auto postmultiply= p1.is_first() ? b[ia] : a[ia];

                auto tmp=copy(world,rhs_batch);
                if (premultiply.is_initialized()) tmp=rhs_batch*premultiply;
                auto tmp1=apply(world,*(f12),tmp);
                if (postmultiply.is_initialized()) tmp1=tmp1*postmultiply;
                tmp2+=tmp1;
            }

            for (auto& t : tmp2) result.push_back(t);
        }
        return result;
    }

    typename Tensor<T>::scalar_type norm2() const override {
        const Function<T, LDIM> one = FunctionFactory<T, LDIM>(world()).f(
                [](const Vector<double, LDIM>& r) { return 1.0; });
        std::vector<Function<T, LDIM>> pre, post;
        std::size_t sz = a.size();
        if (sz == 0) {
            pre = std::vector<Function<T, LDIM>>(1, one);
            post = std::vector<Function<T, LDIM>>(1, one);
        } else {
            pre = (a.front().is_initialized()) ? a : std::vector<Function<T, LDIM>>(sz, one);
            post = (b.front().is_initialized()) ? b : std::vector<Function<T, LDIM>>(sz, one);
        }

        const SeparatedConvolution<T,LDIM>& f12a=*(f12);
        const SeparatedConvolution<T,LDIM> f12sq= SeparatedConvolution<T,LDIM>::combine(f12a,f12a);

        // \int f(1,2)^2 d1d2 = \int f(1,2)^2 pre(1)^2 post(2)^2 d1 d2
        // || \sum_i f(1,2) a_i(1) b_i(2) || = \int ( \sum_{ij} a_i(1) a_j(1) f(1,2)^2 b_i(1) b_j(2) ) d1d2
        std::vector<Function<T,LDIM>> aa,bb;
        for (std::size_t i=0; i<pre.size(); ++i) {
            aa=append(aa,(pre[i]*pre));
            bb=append(bb,(post[i]*post));
        }
//        typename Tensor<T>::scalar_type term1 =madness::inner(mul(world(),post,post),f12sq(mul(world(),pre,pre)));
        typename Tensor<T>::scalar_type term1 =madness::inner(bb,f12sq(aa));
        return sqrt(term1);

    }

    T operator()(const Vector<double,NDIM>& r) const override {

        if (a.size()==0) return 0.0;
        auto split = [](const Vector<double,NDIM>& r) {
            Vector<double,LDIM> first, second;
            for (size_t i=0; i<LDIM; ++i) {
                first[i]=r[i];
                second[i]=r[i+LDIM];
            }
            return std::make_pair(first,second);
        };

        double gamma=f12->info.mu;
        auto [first,second]=split(r);


        double result=0.0;
        for (std::size_t ia=0; ia<a.size(); ++ia) {
            double result1=1.0;
            if (a[ia].is_initialized()) result1*=a[ia](first);
            if (b[ia].is_initialized()) result1*=b[ia](second);
            if (f12->info.type==OT_SLATER) result1*=exp(-gamma*(first-second).normf());
            else if (f12->info.type==OT_F12) result1*=(1.0-exp(-gamma* madness::inner(first-second,first-second)));
            else if (f12->info.type==OT_GAUSS) result1*=exp(-gamma* madness::inner(first-second,first-second));
            else {
                MADNESS_EXCEPTION("no such operator_type",1);
            }
            result+=result1;
        }
        return result;

    }
};

template<typename T, std::size_t NDIM, std::size_t LDIM=NDIM/2>
struct LRFunctorPure : public LRFunctorBase<T,NDIM> {
    LRFunctorPure() = default;
    LRFunctorPure(const Function<T,NDIM>& f) : f(f) {}
    World& world() const override {return f.world();}

    Function<T, NDIM> f;    ///< a hi-dim function

    std::vector<Function<T,LDIM>> inner(const std::vector<Function<T,LDIM>>& rhs,
                                        const particle<LDIM> p1, const particle<LDIM> p2) const override {
        return madness::innerXX<LDIM>(f,rhs,p1.get_array(),p2.get_array());
    }

    std::string type() const override {return "LRFunctorPure";}

    /// evaluate the functor at a given point, e.g. for plotting
    T operator()(const Vector<double,NDIM>& r) const override {
        return f(r);
    }

    typename Tensor<T>::scalar_type norm2() const override {
        return f.norm2();
    }
};


    /// LowRankFunction represents a hi-dimensional (NDIM) function as a sum of products of low-dimensional (LDIM) functions

    /// a LowRankFunction can be created from a hi-dim function directly, or from a composite like f(1,2) phi(1) psi(2),
    /// where f(1,2) is a two-particle function (e.g. a Slater function)
    /// there are two possible representation
    ///  canonical:     f(1,2) = \sum_i g_i(1) h_i(2)
    ///  general:       f(1,2) = \sum_{ij} g_i(1) M_{ij} h_j(2)
    /// for the time being we don't require g or h to be orthogonal or normalized
    template<typename T, std::size_t NDIM, std::size_t LDIM=NDIM/2>
    class LowRankFunction {
    public:

        World& world;
        double rank_revealing_tol=1.e-8;     // rrcd tol
        std::string orthomethod="canonical";
        bool do_print=false;
        std::vector<Function<T,LDIM>> g,h;
        Tensor<T> metric;   ///< the coupling matrix in the general representation, empty in the canonical representation
        const particle<LDIM> p1=particle<LDIM>::particle1();
        const particle<LDIM> p2=particle<LDIM>::particle2();

        LowRankFunction(World& world) : world(world) {}

        LowRankFunction(std::vector<Function<T,LDIM>> g, std::vector<Function<T,LDIM>> h,
                        double tol, std::string orthomethod, const Tensor<T> metric=Tensor<T>()) : world(g.front().world()),
                        rank_revealing_tol(tol), orthomethod(orthomethod), g(g), h(h), metric(metric) {
        }

        /// shallow copy ctor
        LowRankFunction(const LowRankFunction& other) : world(other.world),
            rank_revealing_tol(other.rank_revealing_tol), orthomethod(other.orthomethod),
            g(other.g), h(other.h), metric(other.metric) {
        }

        /// deep copy
        friend LowRankFunction copy(const LowRankFunction& other) {
            return LowRankFunction<T,NDIM>(madness::copy(other.g),madness::copy(other.h),
                other.rank_revealing_tol,other.orthomethod, madness::copy(other.metric));
        }

        LowRankFunction& operator=(const LowRankFunction& f) { // Assignment required for storage in vector
            if (this == &f) return *this;
            rank_revealing_tol = f.rank_revealing_tol;
            orthomethod = f.orthomethod;
            do_print = f.do_print;
            g = f.g;
            h = f.h;
            metric = f.metric;
            return *this;
        }

        /// function evaluation
        T operator()(const Vector<double,NDIM>& r) const {
            Vector<double,LDIM> first, second;
            for (int i=0; i<LDIM; ++i) {
                first[i]=r[i];
                second[i]=r[i+LDIM];
            }
            double result=0.0;
            Tensor<T> gvec(g.size()), hvec(h.size());
            for (int i=0; i<g.size(); ++i) gvec(i)=g[i](first);
            for (int i=0; i<h.size(); ++i) hvec(i)=h[i](second);
            if (is_canonical()) result=gvec.trace(hvec);
            else result=gvec.trace(inner(metric,hvec,1,0));
            return result;
        }

        /// the canonical representation is a special case of the general representation where the coupling matrix
        /// is the identity, so we can check for that
        bool is_canonical() const {
            // some sanity check
            if (metric.size()==0) MADNESS_CHECK_THROW(g.size()==h.size(),"inconsistent sizes in LRF");
            if (metric.size()>0) {
                MADNESS_CHECK_THROW(g.size()==metric.dim(0),"inconsistent g sizes in LRF");
                MADNESS_CHECK_THROW(h.size()==metric.dim(1),"inconsistent h sizes in LRF");
            }
            return metric.size()==0;
        }

        /*
         * arithmetic section
         */

        /// addition
        LowRankFunction operator+(const LowRankFunction& b) const {
            LowRankFunction<T,NDIM> result=copy(*this);
            result+=b;
            return result;
        }
        /// subtraction
        LowRankFunction operator-(const LowRankFunction& b) const {
            LowRankFunction<T,NDIM> result=copy(*this);
            result-=b;
            return result;
        }

        /// in-place addition
        LowRankFunction& operator+=(const LowRankFunction& b) {

            // result metric is (m1  0 \\ 0  m2)
            if (metric or b.metric) {
                Tensor<T> tmp_metric=Tensor<T>(g.size() + b.g.size(), h.size() + b.h.size());
                Tensor<T> ametric= (metric) ? metric : identity_matrix<T>(g.size());
                tmp_metric(Slice(0,g.size()-1),Slice(0,h.size()-1))=ametric;

                Tensor<T> bmetric= (b.metric) ? b.metric : identity_matrix<T>(b.g.size());
                tmp_metric(Slice(g.size(),g.size()+b.g.size()-1),Slice(h.size(),h.size()+b.h.size()-1))=bmetric;
                metric=tmp_metric;
            }

            g=append(g,copy(b.g));
            h=append(h,copy(b.h));

            return *this;
        }

        /// in-place subtraction
        LowRankFunction& operator-=(const LowRankFunction& b) {
            LowRankFunction<T,NDIM> tmp(b);     // shallow
            if (b.metric) tmp.metric=-1.0*b.metric;           // deep
            else tmp.metric=-1.0*identity_matrix<T>(b.g.size());
            return (*this)+=tmp;
        }

        /// scale by a scalar
        template<typename Q>
        LowRankFunction operator*(const Q a) const {
            return LowRankFunction<TensorResultType<T,Q>,NDIM>(g * a, Q(h),rank_revealing_tol,orthomethod,metric);
        }

        /// out-of-place scale by a scalar (no type conversion)
        LowRankFunction operator*(const T a) const {
            return LowRankFunction(g * a, h,rank_revealing_tol,orthomethod,metric);
        }

        /// multiplication with a scalar
        friend LowRankFunction operator*(const T a, const LowRankFunction& other) {
            return other*a;
        }

        /// in-place scale by a scalar (no type conversion)
        LowRankFunction& operator*=(const T a) {
            g=g*a;
            return *this;
        }

        /// l2 norm
        typename TensorTypeData<T>::scalar_type norm2() const {
            auto g_ij=matrix_inner(world,g,g);
            auto h_ij=matrix_inner(world,h,h);
            // ||this||^2 = <g_i| g_k> M_ij <h_j| h_l> M_kl
            // in the canonical case, M_ij = delta_ij, so we just need to contract g_ij and h_ij
            if (is_canonical()) return sqrt(g_ij.trace(h_ij));

            /// in the general case we have to contract with the coupling matrix metric
            auto SM1=inner(g_ij,metric,1,0);     // S_il = <g_i| g_k> M_kl
            auto SM2=inner(metric,h_ij,1,0);     // S_il = M_ij <h_j| h_l>
            return sqrt(SM1.trace(SM2));

        }

        std::vector<Function<T,LDIM>> get_functions(const particle<LDIM>& p) const {
            MADNESS_CHECK(p.is_first() or p.is_last());
            if (p.is_first()) return g;
            return h;
        }

        std::vector<Function<T,LDIM>> get_g() const {return g;}
        std::vector<Function<T,LDIM>> get_h() const {return h;}

        Tensor<long> rank() const {
            Tensor<long> r(2);
            r(0l)=g.size();
            r(1l)=h.size();
            if (metric.size()>0) {
                MADNESS_ASSERT(r(0l)==metric.dim(0));
                MADNESS_ASSERT(r(1l)==metric.dim(1));
            } else {
                MADNESS_ASSERT(r(0l)==r(1l));
            }
            return r;
        }

        /// return the size in GByte
        double size() const {
            double sz=get_size(world,g);
            sz+=get_size(world,h);
            return sz;
        }

        /// f(1,2) = \sum_{pq} g_p(1) M_{pq} h_q(2)
        Function<T,NDIM> reconstruct() const {
            std::vector<Function<T,LDIM>> gtilde=g;
            if (not is_canonical()) gtilde=transform(world,g,metric);        // gtilde_p(1) = \sum_q g_q(1) M_{qp}
            MADNESS_ASSERT(gtilde.size()==h.size());
            auto fapprox=hartree_product(gtilde[0],h[0]);
            for (int i=1; i<gtilde.size(); ++i) fapprox+=hartree_product(gtilde[i],h[i]);
            return fapprox;
        }

        /// remove linear dependencies in g and h; result will have the form
        /// f(1,2) = \sum_{pq} g_p(1) M_{pq} h_q(2)
        /// {g} and {h} are not necessarily orthogonal, but they are linearly independent
        /// the coupling matrix M is not necessarily quadratic
        void remove_linear_dependencies(double tol=-1.0) {
            if (tol<0.0) tol=rank_revealing_tol;

            Tensor<T> ovlp_g=matrix_inner(world,g,g);
            Tensor<T> ovlp_h=matrix_inner(world,h,h);
            auto [pg,Xg,Xg_inv]=rr_cholesky_matrix_and_reorder(g,ovlp_g,tol);
            auto [ph,Xh,Xh_inv]=rr_cholesky_matrix_and_reorder(h,ovlp_h,tol);
            // cite from rr_cholesky_matrix_and_reorder
            //             pv: reordered original, linearly-independent functions,
            //             X: the (r,r) transformation matrix from pv to orthonormalized functions:
            //                 v_ortho_i = v_p X_ip  =  pv_j X_ji
            //             X_inv: its (r,p) "inverse" from original functions to orthonormalized functions:
            //                 v_p = v_ortho_i X_inv_ip
            // note that X is not the inverse of X_inv, that would be X_full
            g=pg;
            h=ph;
            // g M h = pg X X^(-1) M X^(-T) X^T ph
            auto XXinv_g=inner(Xg,Xg_inv);
            auto XXinv_h=inner(Xh,Xh_inv);
            if (metric.size()>0) {
                metric=inner(XXinv_g,inner(metric,XXinv_h,1,1),1,0);      // Xg . M . Xh^T
            } else {
                // if it was canonical keep it canonical
                metric = inner(XXinv_g,XXinv_h,1,1);      // Xg . Xh^T
                canonicalize();
            }
        }

        void reorthonormalize() {
            remove_linear_dependencies();
        }


        /// perform RRCD on the overlap matrix of the functions, reorder the functions according to the pivoting,
        /// and return the upper triangular matrix from the Cholesky decomposition

        /// note: the rr_cholesky is more efficient than a full SVD for removing linear dependencies, but it does not
        /// give an orthonormal basis, so an additional orthonormalization step is needed if orthogonality is desired.
        /// To orthonormalize:
        /// auto [pv, X] = rr_cholesky_matrix_and_reorder(v,ovlp,tol);
        /// v_ortho=transform(world,pv,X); // v_ortho_i = sum_j v_j U_{ji}^(-1)
        /// @param[in] v the input functions
        /// @param[in] ovlp the overlap matrix of the input functions, will be destroyed in-place by the rr_cholesky
        /// @param[in] tol the tolerance for linear dependence in the rr_cholesky
        /// @return a tuple containing
        ///             pv: reordered original, linearly-independent functions,
        ///             X: the (r,r) transformation matrix from pv to orthonormalized functions:
        ///                 v_ortho_i = v_p X_pi  =  pv_j X_ji
        ///             X_inv: its (r,p) "inverse" from original functions to orthonormalized functions:
        ///                 v_p = v_ortho_i X_inv_ip
        /// note that X is not the inverse of X_inv, that would be X_full
        /// but v_p(x) = pv_i(x) X_ij X_inv_jp  is numerically the same
        static std::tuple<std::vector<Function<T,LDIM>>, Tensor<T>, Tensor<T>>
        rr_cholesky_matrix_and_reorder(const std::vector<Function<T,LDIM>>& v, Tensor<T> ovlp, const double tol) {
            int rank;
            Tensor<int> piv;

            rr_cholesky(ovlp,tol,piv,rank); // destroys ovlp and gives back Upper ∆ Matrix from CCD

            // rearrange and truncate the functions according to the pivoting of the rr_cholesky
            std::vector<Function<T,LDIM> > pv(rank);
            for(integer i=0;i<rank;++i) pv[i]=v[piv[i]];

            // compute transformation matrix v -> w = v_ortho
            // w_i = v_p X_pi
            //     = v_piv(i) X_ij
            // v_p = w_i X^(-1)_ip = w_i Xinv_ip
            // pv_i = piv(v_p)      // pv is not orthogonal or normalized, just linearly independent
            // this returns: pv, X_ij, X_inv_ip

            // no need to invert all of ovlp, only the upper left rank x rank block
            Tensor<T> U=ovlp(Slice(0,rank-1),Slice(0,rank-1));
            Tensor<T> Xinv=ovlp(Slice(0,rank-1),_); // has dimensions (rank,original_size)
            Tensor<T> X=inverse(U);             // has dimensions (rank,rank)

            // sanity check
            // result = |i~> = sum_r |i> U_{ir}^(-1)
            // with <i~|j~> = sum_{r} U_{ir}^(-1) <i|j> U_{jp}^(-1)
            //              = sum_{r} U_{ir}^(-1) (U^T U)_{ij} U_{jr}^(-1)
            //              = sum_{r} U_{ir}^(-1) U_{ip} U_{jp} U_{jp}^(-1)
            //              = delta_{ij}        // has matrix dimension (r,r)

            // permute X_inv, so that the following identity holds:
            // also the same: v_p(x) = pv_i(x) X_ij X_inv_jq P_qp
            Tensor<T> tmp=copy(Xinv);
            for (int i=0; i<tmp.dim(1); ++i) Xinv(_,piv(i))=tmp(_,i);

            return std::make_tuple(pv,X,Xinv);
        }

    private:
        /// orthonormalize the argument vector
        std::vector<Function<T,LDIM>> orthonormalize(const std::vector<Function<T,LDIM>>& g) const {

            MADNESS_CHECK_THROW(is_canonical(),"no orthonormalization unless canonicalized");
            double tol=rank_revealing_tol;
            std::vector<Function<T,LDIM>> g2;
            auto ovlp=matrix_inner(world,g,g);
            if (orthomethod=="canonical") {
                tol*=0.01;
                print("orthonormalizing with method/tol",orthomethod,tol);
                g2=orthonormalize_canonical(g,ovlp,tol);
            } else if (orthomethod=="cholesky") {
                print("orthonormalizing with method/tol",orthomethod,tol);
                g2=orthonormalize_rrcd(g,ovlp,tol);
            }
            else {
                MADNESS_EXCEPTION("no such orthomethod",1);
            }
            double tight_thresh=FunctionDefaults<3>::get_thresh()*0.1;
            return truncate(g2,tight_thresh);
        }

    public:

        /// optimize the lrf using the lrfunctor

        /// @param[in]  nopt       number of iterations (wrt to Alg. 4.3 in Halko)
        void optimize(const LRFunctorBase<T,NDIM>& lrfunctor1, const long nopt=1) {
            timer t(world);
            MADNESS_CHECK_THROW(is_canonical(),"currently only optimization of canonical LRFs supported");
            t.do_print=do_print;
            for (int i=0; i<nopt; ++i) {
                // orthonormalize h
                h=truncate(orthonormalize(h));
                t.tag("ortho1");
                g=truncate(inner(lrfunctor1,h,p2,p1));
                t.tag("inner1");
                g=truncate(orthonormalize(g));
                t.tag("ortho2");
                h=truncate(inner(lrfunctor1,g,p1,p1));
                t.tag("inner2");
            }
        }
        void canonicalize() {
            if (is_canonical()) return;
            g=transform(world,g,metric);
            metric=Tensor<T>();
        }

        /// compute the l2 error |functor - \sum_i g_ih_i|_2

        /// \int (f(1,2) - gh(1,2))^2 = \int f(1,2)^2 - 2\int f(1,2) gh(1,2) + \int gh(1,2)^2
        /// since we are subtracting large numbers the numerics are sensitive, and NaN may be returned..
        double l2error(const LRFunctorBase<T,NDIM>& lrfunctor1) const {

            timer t(world);
            t.do_print=do_print;

            // \int f(1,2)^2 d1d2
            double term1 = lrfunctor1.norm2();
            term1=term1*term1;
            t.tag("computing term1");

            // \int f(1,2) \sum_{ij} g_i(1) m_{ij} h_j(2) d1d2
            // = \sum_{ij} m_{ij} \int (\int f(1,2) g_i(1) d1) h_j(2) d2
            double term2=0.0;
            if (is_canonical()) {
                term2=madness::inner(g,inner(lrfunctor1,h,p2,p1));
            } else {
                std::vector<Function<T,LDIM>> fh=inner(lrfunctor1,h,p2,p1);
                Tensor<T> fgh=matrix_inner(world,g,fh);
                term2=fgh.trace(metric);
            }
            t.tag("computing term2");

            double term3=0.0;
            if (is_canonical()) {
                auto tmp1=matrix_inner(world,h,h);
                auto tmp2=matrix_inner(world,g,g);
                term3=tmp1.trace(tmp2);
            } else {
                // general case: no orthogonality, so we have to contract with the coupling matrix metric
                // \int gh(1,2)^2 d1d2 = \int \sum_{ijkl} g_i(1) m_{ij} g_k(1) m_{kl} h_j(2) h_l(2) d1d2
                //   = \sum_{ijkl} m_{ij} m_{kl} <g_i | g_k> <h_j | h_l>
                auto gmat=matrix_inner(world,g,g);
                auto hmat=matrix_inner(world,h,h);
                auto SM1 = inner(gmat,metric);          // gmat_ik m_kl = SM1_il
                auto SM2 = inner(metric, hmat);         // m_ij hmat_jl = SM2_il
                term3=SM1.trace(SM2);
            }
//            print("term3/a/diff",term3a,term3,term3-term3a);
            t.tag("computing term3");

            double arg=term1-2.0*term2+term3;
            if (arg<0.0) {
                print("negative l2 error");
                arg*=-1.0;
//                throw std::runtime_error("negative argument in l2error");
            }
            double error=sqrt(arg)/sqrt(term1);
            if (world.rank()==0 and do_print) {
                print("term1,2,3, error",term1, term2, term3, "  --",error);
            }

            return error;
        }

    };

    /// compute the inner product to 2 LowRankFunctions
    template<typename T, std::size_t NDIM>
    T inner(const LowRankFunction<T,NDIM>& a, const LowRankFunction<T,NDIM>& b) {
        World& world=a.world;

        // result = <a.g_i| b.g_k> a.M_ij <a.h_j| b.h_l> b.M_kl
        auto g_ik=matrix_inner(world,a.g,b.g);
        auto h_jl=matrix_inner(world,a.h,b.h);

        /// in the general case we have to contract with the coupling matrix metric
        auto SM1 = (b.metric.size()) ? inner(g_ik,b.metric,1,0) : g_ik;     // S_il = <g_i| g_k> M_kl
        auto SM2 = (a.metric.size()) ? inner(a.metric,h_jl,1,0) : h_jl;     // S_il = M_ij <h_j| h_l>
        return SM1.trace(SM2);
    }

    /**
     * inner product: LowRankFunction lrf; Function f, g; double d
     *  lrf(1,3) = inner(lrf(1,2), lrf(2,3))
     *  lrf(1,3) = inner(lrf(1,2), f(2,3))
     *  g(1) = inner(lrf(1,2), f(2))
     *  d = inner(lrf(1,2), f(1,2))
     *  d = inner(lrf(1,2), lrf(1,2))
     */

    ///  lrf(1,3) = inner(full(1,2), lrf(2,3))

    /// @param[in] f1 the first function
    /// @param[in] f2 the second function
    /// @param[in] p1 the integration variable of the first function
    /// @param[in] p2 the integration variable of the second function
    template<typename T, std::size_t NDIM, std::size_t PDIM>
    LowRankFunction<T,NDIM> inner(const Function<T,NDIM>& f1, const LowRankFunction<T,NDIM>& f2,
                                  const particle<PDIM> p1, const particle<PDIM> p2) {
        // MADNESS_CHECK_THROW(f2.is_canonical(),"need canonical representation for inner product of two low-rank functions");
        auto result=inner(f2,f1,p2,p1);
        std::swap(result.g,result.h);
        if (result.metric) result.metric=transpose(result.metric);
        return result;
    }

    ///  lrf(1,3) = inner(lrf(1,2), full(2,3))

    /// @param[in] f1 the first function
    /// @param[in] f2 the second function
    /// @param[in] p1 the integration variable of the first function
    /// @param[in] p2 the integration variable of the second function
    template<typename T, std::size_t NDIM, std::size_t PDIM>
    LowRankFunction<T,NDIM> inner(const LowRankFunction<T,NDIM>& f1, const Function<T,NDIM>& f2,
                                  const particle<PDIM> p1, const particle<PDIM> p2) {
        static_assert(TensorTypeData<T>::iscomplex==false, "complex inner in LowRankFunction not implemented");
        World& world=f1.world;
        static_assert(2*PDIM==NDIM);
        // int f(1,2) k(2,3) d2 = \sum \int g_i(1) M_ij h_j(2) k(2,3) d2
        //                      = \sum g_i(1) M_ij \int h_j(2) k(2,3) d2
        //                      = \sum g_i(1) M_ij k_j(3)
        LowRankFunction<T, NDIM> result(world);
        if (p1.is_last()) { // integrate over 2: result(1,3) = lrf(1,2) f(2,3)
            result.g = f1.g;
            change_tree_state(f1.h,reconstructed);
            result.h=innerXX<PDIM>(f2,f1.h,p2.get_array(),particle<PDIM>::particle1().get_array());
            result.metric=f1.metric;
        } else if (p1.is_first()) { // integrate over 1: result(2,3) = lrf(1,2) f(1,3)
            // int f(1,2) k(2,3) d2 = \sum \int g_i(1) M_ij h_j(2) k(1,3) d1
            //                      = \sum h_j(2) M_ij \int g_i(1) k(1,3) d1
            //                      = \sum h_j(2) M_ji k_i(3)
            result.g = f1.h;        // correct! second variable of f1 becomes first variable of result
            change_tree_state(f1.g,reconstructed);
            result.h=innerXX<PDIM>(f2,f1.g,p2.get_array(),particle<PDIM>::particle1().get_array());
            if (f1.metric.size()>0) result.metric=transpose(f1.metric);
        }
        return result;
    }

    ///  lrf(1,3) = inner(lrf(1,2), lrf(2,3))

    /// @param[in] f1 the first function
    /// @param[in] f2 the second function
    /// @param[in] p1 the integration variable of the first function
    /// @param[in] p2 the integration variable of the second function
    template<typename T, std::size_t NDIM, std::size_t PDIM>
    LowRankFunction<T,NDIM> inner(const LowRankFunction<T,NDIM>& f1, const LowRankFunction<T,NDIM>& f2,
                                  const particle<PDIM> p1, const particle<PDIM> p2) {
        World& world=f1.world;
        static_assert(2*PDIM==NDIM);

        // integrate over 2: result(1,3) =
        // p1,p2 = 1,0: inner(lrf(1,2) ,lrf(2,3) ) = \sum_ij g1_i(1) M_ij <h1_j(2) g2_k(2)> h2_l(3) M_kl
        // p1,p2 = 0,0: inner(lrf(2,1) ,lrf(2,3) ) = \sum_ij h1_j(1) M_ij <g1_i(2) g2_k(2)> h2_l(3) M_kl
        // p1,p2 = 0,1: inner(lrf(2,1) ,lrf(3,2) ) = \sum_ij h1_j(1) M_ij <g1_i(2) h2_l(2)> g2_k(3) M_kl
        // p1,p2 = 1,1: inner(lrf(1,2) ,lrf(3,2) ) = \sum_ij g1_i(1) M_ij <h1_j(2) h2_l(2)> g2_k(3) M_kl
        auto matrix=matrix_inner(world,f1.get_functions(p1),f2.get_functions(p2));
        int index1 = (p1.is_first()) ? 0 : 1;
        int index2 = (p2.is_first()) ? 0 : 1;
        auto tmp=inner(f1.metric,matrix,index1,0);
        auto metric=inner(tmp,f2.metric,1,index2);
        auto gg=copy(world,f1.get_functions(p1.complement()));
        auto hh=copy(world,f2.get_functions(p2.complement()));
        return LowRankFunction<T,NDIM>(gg,hh,f1.rank_revealing_tol,f1.orthomethod,metric);
    }

    ///  f(1) = inner(lrf(1,2), f(2))

    /// @param[in] f1 the first function
    /// @param[in] vf vector of the second functions
    /// @param[in] p1 the integration variable of the first function
    /// @param[in] p2 the integration variable of the second function, dummy variable for consistent notation
    template<typename T, std::size_t NDIM, std::size_t PDIM>
    std::vector<Function<T,NDIM-PDIM>> inner(const LowRankFunction<T,NDIM>& f1, const std::vector<Function<T,PDIM>>& vf,
                                  const particle<PDIM> p1, const particle<PDIM> p2=particle<PDIM>::particle1()) {
        World& world=f1.world;
        static_assert(2*PDIM==NDIM);
        MADNESS_CHECK(p2.is_first());

        // inner(lrf(1,2), f_k(2) ) = \sum_i g1_i(1) M_{ij} <h1_j(2) f_k(2)>
        // inner(lrf(2,1), f_k(2) ) = \sum_i h1_j(1) M_{ij} <g1_i(2) f_k(2)>
        auto matrix_jk=matrix_inner(world,f1.get_functions(p1),vf);
        int index1 = (p1.is_first()) ? 0 : 1;
        auto matrix = inner(f1.metric,matrix_jk,index1,0);
        return transform(world,f1.get_functions(p1.complement()),matrix);
    }

    ///  f(1) = inner(lrf(1,2), f(2))

    /// @param[in] f1 the first function
    /// @param[in] vf the second function
    /// @param[in] p1 the integration variable of the first function
    /// @param[in] p2 the integration variable of the second function, dummy variable for consistent notation
    template<typename T, std::size_t NDIM, std::size_t PDIM>
    Function<T,NDIM-PDIM> inner(const LowRankFunction<T,NDIM>& f1, const Function<T,PDIM>& f2,
                                        const particle<PDIM> p1, const particle<PDIM> p2=particle<PDIM>::particle1()) {
        return inner(f1,std::vector<Function<T,PDIM>>({f2}),p1,p2)[0];
    }

    /// Factory class to compute a low-rank approximation of a given hi-dimensional function using randomized
    /// projection and rank-revealing QR decomposition (RRCD)
    template<typename T, std::size_t NDIM, std::size_t LDIM=NDIM/2>
    class LowRankFunctionFactory {
    public:

        const particle<LDIM> p1=particle<LDIM>::particle1();
        const particle<LDIM> p2=particle<LDIM>::particle2();

        LowRankFunctionParameters parameters;
        std::vector<Vector<double,LDIM>> origins;  ///< origins of the molecular grid

        LowRankFunctionFactory() = default;
        LowRankFunctionFactory(const LowRankFunctionParameters param, const std::vector<Vector<double,LDIM>> origins={})
                : parameters(param), origins(origins) {}

        LowRankFunctionFactory(const LowRankFunctionParameters param, const Molecule& molecule)
                : LowRankFunctionFactory(param,molecule.get_all_coords_vec()){}

        LowRankFunctionFactory(const LowRankFunctionFactory& other) = default;

        LowRankFunctionFactory& set_centers(const std::vector<Vector<double,LDIM>> centers) {
            origins=centers;
            return *this;
        }
        LowRankFunctionFactory& set_radius(const double radius) {
            parameters.set_user_defined_value("radius",radius);
            return *this;
        }
        LowRankFunctionFactory& set_volume_element(const double volume_element) {
            parameters.set_user_defined_value("volume_element",volume_element);
            return *this;
        }
        LowRankFunctionFactory& set_rank_revealing_tol(const double rrtol) {
            parameters.set_user_defined_value("tol",rrtol);
            return *this;
        }
        LowRankFunctionFactory& set_orthomethod(const std::string orthomethod) {
            parameters.set_user_defined_value("orthomethod",orthomethod);
            return *this;
        }
        LowRankFunctionFactory& set_canonicalize(const bool canonicalize) {
            parameters.set_user_defined_value("canonicalize",canonicalize);
            return *this;
        }
        LowRankFunctionFactory& set_gridtype(const std::string& gridtype) {
            parameters.set_user_defined_value("gridtype",gridtype);
            return *this;
        }
        LowRankFunctionFactory& set_adaptive_coarse_factor(const double factor) {
            parameters.set_user_defined_value("adaptive_coarse_factor",factor);
            return *this;
        }
        LowRankFunctionFactory& set_adaptive_refine_radius(const double ratio) {
            parameters.set_user_defined_value("adaptive_refine_radius",ratio);
            return *this;
        }
        LowRankFunctionFactory& set_adaptive_significance_ratio(const double ratio) {
            parameters.set_user_defined_value("adaptive_significance_ratio",ratio);
            return *this;
        }
        LowRankFunctionFactory& set_adaptive_max_centers(const int max_centers) {
            parameters.set_user_defined_value("adaptive_max_centers",max_centers);
            return *this;
        }
        LowRankFunctionFactory& set_adaptive_min_centers(const int min_centers) {
            parameters.set_user_defined_value("adaptive_min_centers",min_centers);
            return *this;
        }

        struct YFormationResult {
            std::vector<Function<T,LDIM>> Y;
            std::vector<double> norms;
            std::vector<std::size_t> significant_indices;
        };

    private:
        std::vector<Vector<double,LDIM>> make_uniform_random_grid_in_cell(const double volume_element) const {
            auto cell = FunctionDefaults<LDIM>::get_cell();
            double volume = 1.0;
            for (size_t d = 0; d < LDIM; ++d) volume *= (cell(d,1) - cell(d,0));
            long npoint = std::max(1l, long(volume / volume_element));

            std::random_device rd{};
            std::mt19937 gen{rd()};
            std::vector<std::uniform_real_distribution<double>> dist;
            for (size_t d = 0; d < LDIM; ++d) dist.emplace_back(cell(d,0), cell(d,1));

            std::vector<Vector<double,LDIM>> grid;
            grid.reserve(npoint);
            for (long i = 0; i < npoint; ++i) {
                Vector<double,LDIM> r;
                for (size_t d = 0; d < LDIM; ++d) r[d] = dist[d](gen);
                grid.push_back(r);
            }
            return grid;
        }

        std::vector<std::size_t> pick_significant_indices(const std::vector<double>& norms) const {
            if (norms.empty()) return {};
            std::vector<std::size_t> idx(norms.size());
            for (std::size_t i = 0; i < idx.size(); ++i) idx[i] = i;
            std::sort(idx.begin(), idx.end(), [&](const std::size_t a, const std::size_t b) {
                return norms[a] > norms[b];
            });

            const double maxnorm = norms[idx.front()];
            const double cutoff = std::max(parameters.tol(), maxnorm * parameters.adaptive_significance_ratio());
            const std::size_t max_keep = std::max(1, parameters.adaptive_max_centers());
            const std::size_t min_keep = std::max(1, parameters.adaptive_min_centers());

            std::vector<std::size_t> keep;
            for (auto i : idx) {
                if (norms[i] >= cutoff) keep.push_back(i);
                if (keep.size() >= max_keep) break;
            }
            if (keep.size() < min_keep) {
                keep.clear();
                const std::size_t n = std::min<std::size_t>(max_keep, idx.size());
                for (std::size_t i = 0; i < n; ++i) keep.push_back(idx[i]);
            }
            return keep;
        }

    public:
        LowRankFunction<T,NDIM> project(const LRFunctorBase<T,NDIM>& lrfunctor) const {
            World& world=lrfunctor.world();
            bool do_print=true;
            timer t1(world);
            t1.do_print=do_print;

            // get sampling grid
            std::vector<Vector<double,LDIM>> grid;
            if (parameters.gridtype()=="adaptive") {
                const double coarse_ve = parameters.volume_element() * std::max(1.0, parameters.adaptive_coarse_factor());
                auto coarse_grid = make_uniform_random_grid_in_cell(coarse_ve);
                auto coarse = Yformer(lrfunctor, coarse_grid, parameters, 30.0, 0.0);
                auto centers = pick_significant_indices(coarse.norms);


                // Baseline global coverage comes from the coarse probe itself.
                grid = coarse_grid;

                const double local_radius = std::max(parameters.volume_element(),
                                                     parameters.radius() * parameters.adaptive_refine_radius());
                for (auto icenter : centers) {
                    randomgrid<LDIM> rg(parameters.volume_element(), local_radius, coarse_grid[icenter]);
                    auto local = rg.get_grid();
                    grid.insert(grid.end(), local.begin(), local.end());
                }

                // Robust fallback if coarse probing produced no points (degenerate cell/VE settings).
                if (grid.empty()) {
                    grid = make_uniform_random_grid_in_cell(parameters.volume_element());
                }
            } else if (parameters.gridtype()=="twostage") {
                // go through all centers and place a probe grid point in each octand around the center, at distance 0.2
                // if at least on of the probe points has a significant response, keep the center and continue
                // with the random grid around the center, otherwise discard the center
                const double probe_distance = 0.2;
                const double probe_ve = parameters.volume_element() * std::max(1.0, parameters.adaptive_coarse_factor());
                std::vector<Vector<double,LDIM>> probe_points;
                for (const auto& origin : origins) {
                    for (int octant = 0; octant < (1 << LDIM); ++octant) {
                        Vector<double,LDIM> probe_point = origin;
                        for (size_t d = 0; d < LDIM; ++d) {
                            probe_point[d] += ((octant & (1 << d)) ? 1 : -1) * probe_distance;
                        }
                        probe_points.push_back(probe_point);
                    }
                    auto probe = Yformer(lrfunctor, probe_points, parameters, 30.0, 0.0);
                    if (std::any_of(probe.norms.begin(), probe.norms.end(), [&](double norm) {
                        return norm >= parameters.tol();
                    })) {
                        randomgrid<LDIM> rg(parameters.volume_element(), parameters.radius(), origin);
                        auto local = rg.get_grid();
                        grid.insert(grid.end(), local.begin(), local.end());
                    }
                }
                MADNESS_CHECK_THROW(not grid.empty(),"grid is empty");
            } else {
                molecular_grid<LDIM> mgrid(origins,parameters);
                grid=mgrid.get_grid();
            }
            print("initial grid size",grid.size());

            auto yformed = Yformer(lrfunctor,grid,parameters);
            auto Y = yformed.Y;
            if (Y.empty()) {
                auto retry = Yformer(lrfunctor,grid,parameters,30.0,0.0);
                Y = retry.Y;
            }
            MADNESS_CHECK_THROW(!Y.empty(),"Yformer generated no basis functions for projection");
            t1.tag("Yforming");
            print("y.size()",Y.size());

            double tol=parameters.tol();
            Tensor<T> X;

            auto ovlp=matrix_inner(world,Y,Y);  // error in symmetric matrix_inner, use non-symmetric form here!
            auto [pY, t, tinv]=LowRankFunction<T,NDIM>::rr_cholesky_matrix_and_reorder(Y,ovlp,tol);
            t1.tag("remove linear dependence");
            // pY is now a set of linearly independent basis functions

            if (parameters.orthomethod()=="cholesky") {
                X=t;
            } else if (parameters.orthomethod()=="symmetric") {
                ovlp=matrix_inner(world,pY,pY);
                X=orthonormalize_symmetric_matrix(ovlp);
            } else if (parameters.orthomethod()=="canonical") {
                ovlp=matrix_inner(world,pY,pY);
                // since we have already removed the lindep functions with the rr_cholesky, set tol=0.0
                X=canonical_orthonormalization_matrix<T,NDIM>(world,ovlp,0.0);
            } else {
                print("unknown orthogonalization method",parameters.orthomethod());
                MADNESS_EXCEPTION("no such orthomethod",1);
            }

            // some diagnosis on the numerics of the metric
            Tensor<T> metric=inner(X,X,1,1); // metric = t t^T
            auto condition=condition_number(metric);
            if (condition.front()>1.e5) {
                auto conditionX = condition_number(X);
                print("warning: ill-conditioned half-metric X in low-rank function projection",conditionX);
                print("warning: ill-conditioned metric in low-rank function projection",condition);
            }

            // default: g is just Y without linear dependencies, no orthonormalization
            auto g=pY;

            // canonicalization: orthonormalize g
            if (parameters.canonicalize()) {
                g=truncate(transform(world,g,X));
                metric.clear();
                t1.tag("Y orthonormalization");
            }
            auto h=truncate(inner(lrfunctor,g,p1,p1));
            t1.tag("Y backprojection");

            LowRankFunction<T,NDIM> result(g,h,parameters.tol(),parameters.orthomethod(),metric);
            result.remove_linear_dependencies(); // improves numerical stability
            t1.tag("removing lindep");

            return result;
        }


        /// apply a rhs on grid points to the hi-dim function and form Y = A_ij w_j (in Halko's language).
        /// The RHS selection (localized Gaussian/exponential vs harmonics) is controlled by
        /// the single `gridtype` parameter: values "random"/"adaptive"/"twostage" use
        /// localized Gaussian RHS (the previous "exponential" behavior), whereas
        /// "harmonics" uses the RI harmonics-based RHS.
        YFormationResult Yformer(const LRFunctorBase<T,NDIM>& lrfunctor1,
            const std::vector<Vector<double,LDIM>>& grid,
            const LowRankFunctionParameters& parameters,
            const double exponent=30.0,
            const double significance_tol=-1.0) const {

            World& world=lrfunctor1.world();
            std::vector<Function<double,LDIM>> Y;
            if (parameters.gridtype()=="harmonics") { // use harmonics-based LHS
                Y = harmonic_basis(world, {parameters.radius(), parameters.radius() * 10.0, 2.0}, parameters.lmax(), origins);
            } else {
                // default: use localized Gaussian RHS (for gridtype values other than "harmonics")
                std::vector<Function<double,LDIM>> omega;
                double coeff=std::pow(2.0*exponent/constants::pi,0.25*LDIM);
                for (const auto& point : grid) {
                    omega.push_back(FunctionFactory<double,LDIM>(world)
                        // .thresh(1.e-6)      // make sure this is not undersampled
                        .special_points({point})
                                            .functor([&point,&exponent,&coeff](const Vector<double,LDIM>& r)
                                                     {
                                                         auto r_rel=r-point;
                                                         return coeff*exp(-exponent*madness::inner(r_rel,r_rel));
                                                     }));
                }
                Y=inner(lrfunctor1,omega,p2,p1);
            }
            auto norms=norm2s(world,Y);
            std::vector<Function<double,LDIM>> Ynormalized;
            std::vector<std::size_t> significant_indices;
            const double ytol = (significance_tol<0.0) ? parameters.tol() : significance_tol;

            for (size_t i=0; i<Y.size(); ++i) {
                if (norms[i]>ytol) {
                    Ynormalized.push_back(Y[i]);
                    significant_indices.push_back(i);
                }
            }
            if (not Ynormalized.empty()) normalize(world,Ynormalized);
            return {Ynormalized,norms,significant_indices};
        }

        /// return a set of solid harmonic functions up to lmax with a zeta range on the given (atomic) centers

        /// @param[in] World
        /// @param[in] zeta_range range of the Gaussian exponent zeta, given as {zeta_min, zeta_max, zeta_step}
        ///             for even-tempered progression
        /// @param[in] lmax maximum angular momentum of the solid harmonics
        /// @param[in] centers the centers of the solid harmonics, typically the atomic positions
        static std::vector<Function<T, LDIM>> harmonic_basis(World& world,
                                                             const std::vector<double> zeta_range, const int lmax,
                                                             const std::vector<Vector<double, LDIM>> centers)
        {
            std::vector<Function<T, LDIM>> harmonics;
            // harmonics-based Y
            struct GaussianFunction : public FunctionFunctorInterface<double, LDIM>
            {
                std::vector<long> ijk; // angular momentum
                Vector<double, LDIM> center;
                double zeta = 1.0;
                GaussianFunction() = default;

                GaussianFunction(const std::vector<long>& ijk, const Vector<double, LDIM>& center, const double zeta)
                    : ijk(ijk), center(center), zeta(zeta)
                {
                }

                double operator()(const Vector<double, LDIM>& r) const override
                {
                    auto r_rel = r - center;
                    double val = exp(-zeta * madness::inner(r_rel, r_rel));
                    for (int d = 0; d < LDIM; ++d) val *= std::pow(r_rel[d], ijk[d]);
                    return val;
                }
            };

            // compute the monomial exponents for the Cartesian Gaussian functions up to a certain angular momentum
            std::vector<std::vector<long>> types;
            for (int l = 0; l <= lmax; ++l)
            {
                std::vector<long> ijk(LDIM, 0);
                std::function<void(int, int)> generate_ijk = [&](int pos, int remaining_l)
                {
                    if (pos == LDIM - 1)
                    {
                        ijk[pos] = remaining_l;
                        types.push_back(ijk);
                        return;
                    }
                    for (int i = 0; i <= remaining_l; ++i)
                    {
                        ijk[pos] = i;
                        generate_ijk(pos + 1, remaining_l - i);
                    }
                };
                generate_ijk(0, l);
            }

            {
                std::vector<double> zetas;
                double z = zeta_range[0];
                while (z < zeta_range[1])
                {
                    zetas.push_back(z);
                    z *= zeta_range[2];
                    if (zetas.size() > 100)
                    {
                        print("too many zetas in even-tempered basis, check your parameters.tempered", zeta_range);
                        break;
                    }
                }

                for (auto type : types)
                {
                    for (auto zeta : zetas)
                    {
                        for (auto& center : centers)
                        {
                            GaussianFunction gf(type, center, zeta);
                            Function<double, LDIM> f = FunctionFactory<double, LDIM>(world).functor(gf);
                            harmonics.push_back(f);
                        }
                    }
                }
            }
            return harmonics;
        }


    };


} // namespace madness

#endif //MADNESS_LOWRANKFUNCTION_H
