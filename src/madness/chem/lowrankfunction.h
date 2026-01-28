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
            initialize<std::string>("gridtype","random","the grid type",{"random","cartesian","dftgrid"});
            initialize<std::string>("rhsfunctiontype","exponential","the type of function",{"exponential"});
            initialize<int>("optimize",1,"number of optimization iterations");
        }
        std::string get_tag() const override {
            return std::string("lrf");
        }


        void read_and_set_derived_values(World& world, const commandlineparser& parser, std::string tag) {
            read_input_and_commandline_options(world,parser,tag);
        }

        double radius() const {return get<double>("radius");}
        double gamma() const {return get<double>("gamma");}
        double volume_element() const {return get<double>("volume_element");}
        double tol() const {return get<double>("tol");}
        int optimize() const {return get<int>("optimize");}
        std::string gridtype() const {return get<std::string>("gridtype");}
        std::string orthomethod() const {return get<std::string>("orthomethod");}
        std::string rhsfunctiontype() const {return get<std::string>("rhsfunctiontype");}
        std::string f12type() const {return get<std::string>("f12type");}
    };


    class gridbase {
    public:
        double get_volume_element() const {return volume_element;}
        double get_radius() const {return radius;}

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
            if (NDIM==2) return constants::pi*radius*radius;
            if (NDIM==3) return 4.0 / 3.0 * constants::pi * std::pow(radius, 3.0);
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
        long n_per_dim;
        long total_n;
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
            if (params.gridtype()=="random") grid_builder=std::make_shared<randomgrid<NDIM>>(params.volume_element(),params.radius());
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
        molecular_grid(const Molecule& molecule, std::shared_ptr<gridbase> grid) : molecular_grid(molecule.get_all_coords_vec(),grid) {}

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

    particle(const int p) : particle(std::vector<int>(1,p)) {}
    particle(const int p1, const int p2) : particle(std::vector<int>({p1,p2})) {}
    particle(const int p1, const int p2,const int p3) : particle(std::vector<int>({p1,p2,p3})) {}
    particle(const std::vector<int> p) {
        for (int i=0; i<PDIM; ++i) dims[i]=p[i];
    }

    std::string str() const {
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
struct LRFunctorBase {

    virtual ~LRFunctorBase() {};
    virtual std::vector<Function<T,LDIM>> inner(const std::vector<Function<T,LDIM>>& rhs,
                                        const particle<LDIM> p1, const particle<LDIM> p2) const =0;

    virtual Function<T,LDIM> inner(const Function<T,LDIM>& rhs, const particle<LDIM> p1, const particle<LDIM> p2) const {
        return inner(std::vector<Function<T,LDIM>>({rhs}),p1,p2)[0];
    }

    virtual T operator()(const Vector<T,NDIM>& r) const =0;
    virtual typename Tensor<T>::scalar_type norm2() const {
        MADNESS_EXCEPTION("L2 norm not implemented",1);
    }

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


private:
    std::shared_ptr<SeparatedConvolution<T,LDIM>> f12;  ///< a two-particle function
    std::vector<Function<T,LDIM>> a,b;   ///< the lo-dim functions
public:

    World& world() const {return f12->get_world();}
    std::vector<Function<T,LDIM>> inner(const std::vector<Function<T,LDIM>>& rhs,
                                        const particle<LDIM> p1, const particle<LDIM> p2) const {

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

    typename Tensor<T>::scalar_type norm2() const {
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

    T operator()(const Vector<double,NDIM>& r) const {

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
            if (b[ia].is_initialized()) result1*=b[ia](first);
            if (f12->info.type==OT_SLATER) result1*=exp(-gamma*(first-second).normf());
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
    World& world() const {return f.world();}

    Function<T, NDIM> f;    ///< a hi-dim function

    std::vector<Function<T,LDIM>> inner(const std::vector<Function<T,LDIM>>& rhs,
                                        const particle<LDIM> p1, const particle<LDIM> p2) const {
        return madness::innerXX<LDIM>(f,rhs,p1.get_array(),p2.get_array());
//        std::vector<Function<T,LDIM>> result;
//        for (const auto& r : rhs) result.push_back(madness::inner(f,r,p1.get_tuple(),p2.get_tuple()));
//        return result;
    }

    T operator()(const Vector<double,NDIM>& r) const {
        return f(r);
    }

    typename Tensor<T>::scalar_type norm2() const {
        return f.norm2();
    }
};


    /// LowRankFunction represents a hi-dimensional (NDIM) function as a sum of products of low-dimensional (LDIM) functions

    /// f(1,2) = \sum_i g_i(1) h_i(2)
    /// a LowRankFunction can be created from a hi-dim function directly, or from a composite like f(1,2) phi(1) psi(2),
    /// where f(1,2) is a two-particle function (e.g. a Slater function)
    template<typename T, std::size_t NDIM, std::size_t LDIM=NDIM/2>
    class LowRankFunction {
    public:

        World& world;
        double rank_revealing_tol=1.e-8;     // rrcd tol
        std::string orthomethod="canonical";
        bool do_print=false;
        std::vector<Function<T,LDIM>> g,h;
        const particle<LDIM> p1=particle<LDIM>::particle1();
        const particle<LDIM> p2=particle<LDIM>::particle2();

        LowRankFunction(World& world) : world(world) {}

        LowRankFunction(std::vector<Function<T,LDIM>> g, std::vector<Function<T,LDIM>> h,
                        double tol, std::string orthomethod) : world(g.front().world()),
                        rank_revealing_tol(tol), orthomethod(orthomethod), g(g), h(h) {

        }

        /// shallow copy ctor
        LowRankFunction(const LowRankFunction& other) : world(other.world),
            rank_revealing_tol(other.rank_revealing_tol), orthomethod(other.orthomethod),
            g(other.g), h(other.h) {
        }

        /// deep copy
        friend LowRankFunction copy(const LowRankFunction& other) {
            return LowRankFunction<T,NDIM>(madness::copy(other.g),madness::copy(other.h),other.rank_revealing_tol,other.orthomethod);
        }

        LowRankFunction& operator=(const LowRankFunction& f) { // Assignment required for storage in vector
            LowRankFunction ff(f);
            std::swap(ff.g,g);
            std::swap(ff.h,h);
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
            for (int i=0; i<rank(); ++i) result+=g[i](first)*h[i](second);
            return result;
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

            g=append(g,copy(b.g));
            h=append(h,copy(b.h));
            return *this;
        }

        /// in-place subtraction
        LowRankFunction& operator-=(const LowRankFunction& b) {
            g=append(g,-1.0*b.g);   // operator* implies deep copy of b.g
            h=append(h,copy(b.h));
            return *this;
        }

        /// scale by a scalar
        template<typename Q>
        LowRankFunction operator*(const Q a) const {
            return LowRankFunction<TensorResultType<T,Q>,NDIM>(g * a, Q(h),rank_revealing_tol,orthomethod);
        }

        /// out-of-place scale by a scalar (no type conversion)
        LowRankFunction operator*(const T a) const {
            return LowRankFunction(g * a, h,rank_revealing_tol,orthomethod);
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
            auto tmp1=matrix_inner(world,h,h);
            auto tmp2=matrix_inner(world,g,g);
            return sqrt(tmp1.trace(tmp2));
        }

        std::vector<Function<T,LDIM>> get_functions(const particle<LDIM>& p) const {
            MADNESS_CHECK(p.is_first() or p.is_last());
            if (p.is_first()) return g;
            return h;
        }

        std::vector<Function<T,LDIM>> get_g() const {return g;}
        std::vector<Function<T,LDIM>> get_h() const {return h;}

        long rank() const {return g.size();}

        /// return the size in GByte
        double size() const {
            double sz=get_size(world,g);
            sz+=get_size(world,h);
            return sz;
        }

        Function<T,NDIM> reconstruct() const {
            auto fapprox=hartree_product(g[0],h[0]);
            for (int i=1; i<g.size(); ++i) fapprox+=hartree_product(g[i],h[i]);
            return fapprox;
        }

        /// orthonormalize the argument vector
        std::vector<Function<T,LDIM>> orthonormalize(const std::vector<Function<T,LDIM>>& g) const {

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


        /// optimize the lrf using the lrfunctor

        /// @param[in]  nopt       number of iterations (wrt to Alg. 4.3 in Halko)
        void optimize(const LRFunctorBase<T,NDIM>& lrfunctor1, const long nopt=1) {
            timer t(world);
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

        /// remove linear dependencies using rank-revealing cholesky decomposition
        ///
        /// @return this with g orthonormal and h orthogonal
        void reorthonormalize_rrcd(double thresh=-1.0) {

            // using the following terms from the CD of the overlap matrix (see rr_cholesky)
            //    PT A P = U^T U    // P permutation matrix, U upper triangular
            // for this case:
            //    ovlp = <g_i | g_j> = P UT U PT
            //    (UT)^-1 = U^(-1)T         // pseudo inversion and transposition commute
            //    PT = P^(-1)               // P is a unitary permutation matrix
            //    U^(-1) = UT (U UT)^(-1)
            // the last line assumes U is of the form (r,n) with r the rank, and n the original number of functions
            //
            //    f(1,2) = sum_i  g_i(1) h_i(2)
            //
            // start with h
            //   Sh = <h_i | h_j> = Ph UhT Uh PhT
            //   h_ortho_i = h_j Ph_jk Uh^(-1)_ki
            // as <h_ortho_i | h_ortho_j> = UhT^(-1) PhT <h_i | h_j> Ph (Uh^(-1)) = UhT^(-1) PhT Ph UhT Uh PhT Ph Uh^(-1) = 1
            //   thus
            //   f(1,2) = g_i(1) h_i(2)
            //          = g_i(1) h_j(2) Ph_jk Uh^(-1)_kl Uh_lm PhT_mi
            //          = ( g_i(1) Ph_im Uh_lm) ( h_j(2) Ph_jk Uh^(-1)_kl)
            //          = g'_l(1) h'_l(2)
            // yields orthogonalized h'(2).
            // now orthonormalize g'
            //   g'_i = g'_j Ph_jk Uh_ki
            //   Sg' = <g'_i | g'_j> =  UhT_ik Ph_lk <g_l | g_m> Ph_mn Uh_nj
            // as <g'_ortho_i | g'_ortho_j> = .. see above .. = 1
            // thus
            //   f(1,2) = g'_i(1) h'_i(2)
            //          = g'_i(1) h'_j(2) Pg'_jk Ug'^(-1)_kl Ug'_lm Pg'T_mi
            //          = ( g'_i(1) Pg'T_mi Ug'lm ) ( h'_j(2) Pg'_jk Ug'^(-1)_kl )
            //          = g''_l(1) h''_l(2)
            // yields orthonormal g''(1) and orthogonal h''(2)
            //
            // collecting everything:
            //   Sh = <h_i | h_j> = Ph UhT Uh PhT
            //   Uh^(-1) = UhT (Uh UhT)^(-1)
            //   Sg = <g_i | g_j> = Pg UgT Ug PgT
            //   g'_i = g'_j Ph_jk Uh_ik
            //   h'_i = h_j Ph_jk Uh^(-1)_ki
            //   Sg' = <g'_i | g'_j> =  UhT_ki Ph_kl <g_l | g_m> Ph_mn Uh_jn = Pg' Ug'T Ug' Pg'T
            //   Ug'^(-1) = Ug'T (Ug' Ug'T)^(-1)
            //   g''_i(1) = g'_j(1) Pg'T_mj Ug'im
            //            = g_k(1) Ph_kj Uh_mj Pg'T_nm Ug'in
            //            = g_k(1) Ph_kj Tg_ji
            //   h''_i(2) = h'_j(2) Pg'_jk Ug'^(-1)_ki )
            //            = h_l(2) Ph_lj Uh^(-1)_jk Pg'_kn Ug'^(-1)_ni
            //            = h_l(2) Ph_lj Th_ji
            double tol=rank_revealing_tol;

            // no pivot_inverse is necessary..
            auto pivot_vec = [&] (const std::vector<Function<T,LDIM>>& v, const Tensor<integer>& ipiv) {
                std::size_t rank=ipiv.size();
                std::vector<Function<T,LDIM> > pv(rank);
                for(int i=0;i<ipiv.size();++i) pv[i]=v[ipiv[i]];       // all elements in [rank,v.size()) = 0 due to lindep
                return pv;
            };

//            auto pivot_mat = [&] (const Tensor<T>& t, const Tensor<integer>& ipiv) {
//                // std::size_t rank=ipiv.size();
//                Tensor<T> pt(t.dim(0),t.dim(1));
//                // for(int i=0;i<ipiv.size();++i) pt(i,_)=t(ipiv[i],_);       // all elements in [rank,v.size()) = 0 due to lindep
//                for(int i=0;i<ipiv.size();++i) pt(ipiv[i],_)=t(i,_);       // all elements in [rank,v.size()) = 0 due to lindep
//                return pt;
//            };

            auto pivot = [](const Tensor<integer>& ipiv) {
                Tensor<T> piv(ipiv.size(),ipiv.size());
                for (int i=0; i<ipiv.size(); ++i) piv(ipiv[i],i)=1.0;
                return piv;
            };

            auto U = [&](Tensor<T> ovlp) {
                Tensor<integer> piv;
                auto A=copy(ovlp);
                int rank;
                rr_cholesky(A,tol,piv,rank); // destroys ovlp and gives back Upper ∆ Matrix from CD
                auto U1=copy(A(Slice(0,rank-1),_));
                // test
                {
                    auto P = pivot(piv);
                    auto PAP = inner(P,inner(ovlp,P),0,0);; // Pt A P
                    auto test=madness::inner(U1,U1,0,0);    // UT U
                    double err=(PAP-test).normf()/PAP.normf();
                    // print("A",ovlp);
                    // print("P",P);
                    // print("PAP",PAP);
                    // print("test",test);
                    print("error in rr_cholesky",err,"rank/tot",rank,ovlp.dim(0));
                    print("dim(U1)",U1.dim(0),U1.dim(1));
                    MADNESS_CHECK_THROW(err<1.e-10,"rr_cholesky failed");
                }
                return std::make_pair(U1,piv);
            };

            // inverse is U^(-1) = U^T (U U^T)^(-1)
            auto pseudo_invert = [&](Tensor<T>& U) {
                // print("U",U);
                Tensor<T> tmp=madness::inner(U,U,1,1);      // U U^T
                Tensor<T> tmp_inv=inverse(tmp);
                auto Uinv=madness::inner(U,tmp_inv,0,0);      // U^T (U U^T)^(-1)
                // test
                if (0) {            // this is a left inverse, and not applicable here
                    auto test1=madness::inner(Uinv,U,1,0);   // U^(-1) U
                    for (int i=0; i<test1.dim(0); ++i) test1(i,i)-=1.0;
                    print("dimensions in pseudo_invert U-1 U",U.dim(0),U.dim(1),test1.dim(0),test1.dim(1));
                    print("err in pseudo_invert",test1.normf());
                }
                {               // this is the right inverse, and applicable here
                    auto test1=madness::inner(U,Uinv,1,0);   // U U^(-1)
                    for (int i=0; i<test1.dim(0); ++i) test1(i,i)-=1.0;
                    print("dimensions in pseudo_invert U U-1",U.dim(0),U.dim(1),test1.dim(0),test1.dim(1));
                    print("err in pseudo_invert",test1.normf());
                }

                return Uinv;
            };

            auto Sh=matrix_inner(world,h,h);
            auto [Uh,ph]=U(Sh);              // Uh(r,n), Ph(n,n);
            auto Uh_inv=pseudo_invert(Uh);

            auto Sg=matrix_inner(world,g,g);
            auto PhUhT=madness::inner(pivot( ph),Uh,1,1);   // Ph_jk Uh_ik
            auto Sgprime=madness::inner(PhUhT,madness::inner(Sg,PhUhT,1,0),0,0); // UP Sg PU
            auto [Ugprime,pgprime]=U(Sgprime);
            auto Ugprime_inv=pseudo_invert(Ugprime);
            auto Tg=madness::inner(Uh, inner(pivot(pgprime),Ugprime,0,1),0,0);  // UhT Pg'T Ug'^(-1)T
            auto Th=madness::inner(Uh_inv, inner(pivot(pgprime),Ugprime_inv,1,0),1,0);  // Uh Pg' Ug'^(-1)
            print("dim(Tg)",Tg.dim(0),Tg.dim(1));
            print("dim(Th)",Th.dim(0),Th.dim(1));
            print("g.size()",g.size());
            print("h.size()",h.size());

            g=pivot_vec(g,ph);
            h=pivot_vec(h,ph);
            g=transform(world,g,Tg);
            h=transform(world,h,Th);


        }

        void reorthonormalize(double thresh=-1.0) {
            if (orthomethod=="canonical") {
                reorthonormalize_canonical(thresh);
            } else if (orthomethod=="cholesky") {
                reorthonormalize_rrcd(thresh);
            } else {
                MADNESS_EXCEPTION("no such orthomethod",1);
            }
        }
        /// after external operations g might not be orthonormal and/or optimal -- reorthonormalize

        /// orthonormalization similar to Bischoff, Harrison, Valeev, JCP 137 104103 (2012), Sec II C 3
        /// f   =\sum_i g_i h_i
        ///     = g X- (X+)^T (Y+)^T Y- h
        ///     = g X-  U S V^T  Y- h
        ///     = g (X- U) (S V^T Y-) h
        /// requires 2 matrix_inner and 2 transforms. g and h are optimal
        /// @param[in]  thresh        SVD threshold
        void reorthonormalize_canonical(double thresh=-1.0) {
            if (thresh<0.0) thresh=rank_revealing_tol;
            Tensor<T> ovlp_g = matrix_inner(world, g, g);
            Tensor<T> ovlp_h = matrix_inner(world, h, h);

            ovlp_g=0.5*(ovlp_g+transpose(ovlp_g));
            ovlp_h=0.5*(ovlp_h+transpose(ovlp_h));

            auto [eval_g, evec_g] = syev(ovlp_g);
            auto [eval_h, evec_h] = syev(ovlp_h);

            // get relevant part of the eigenvalues and eigenvectors
            // eigenvalues are sorted in ascending order
            auto get_slice = [](auto eval, double thresh) {
                // remove small/negative eigenvalues
                eval.screen(thresh);
                Slice s;
                for (int i=0; i<eval.size(); ++i) {
                    MADNESS_CHECK_THROW(eval[i]>=0.0,"negative eigenvalues in reorthonormalize");
                    if (eval[i]>thresh) {
                        return s=Slice(i,-1);       // from i to the end
                        break;
                    }
                }
                return s;
            };

            Slice gslice=get_slice(eval_g,1.e-12);
            Slice hslice=get_slice(eval_h,1.e-12);

            Tensor<T> Xplus=copy(evec_g(_,gslice));
            Tensor<T> Xminus=copy(evec_g(_,gslice));
            Tensor<T> Yplus=copy(evec_h(_,hslice));
            Tensor<T> Yminus=copy(evec_h(_,hslice));
            eval_g=copy(eval_g(gslice));
            eval_h=copy(eval_h(hslice));

            for (int i=0; i<eval_g.size(); ++i) Xplus(_,i)*=std::pow(eval_g(i),0.5);
            for (int i=0; i<eval_g.size(); ++i) Xminus(_,i)*=std::pow(eval_g(i),-0.5);
            for (int i=0; i<eval_h.size(); ++i) Yplus(_,i)*=std::pow(eval_h(i),0.5);
            for (int i=0; i<eval_h.size(); ++i) Yminus(_,i)*=std::pow(eval_h(i),-0.5);

            Tensor<T> M=madness::inner(Xplus,Yplus,0,0);    // (X+)^T Y+
            auto [U,s,VT]=svd(M);

            // truncate
            typename Tensor<T>::scalar_type s_accumulated=0.0;
            int i=s.size()-1;
            for (;i>=0; i--) {
                s_accumulated+=s[i];
                if (s_accumulated>thresh) {
                    i++;
                    break;
                }
            }
            for (int j=0; j<s.size(); ++j) VT(j,_)*=s[j];
            Tensor<T> XX=madness::inner(Xminus,U,1,0);
            Tensor<T> YY=madness::inner(Yminus,VT,1,1);

            g=truncate(transform(world,g,XX));
            h=truncate(transform(world,h,YY));
        }


        double check_orthonormality(const std::vector<Function<T,LDIM>>& v) const {
            Tensor<T> ovlp=matrix_inner(world,v,v);
            return check_orthonormality(ovlp);
        }

        double check_orthonormality(const Tensor<T>& ovlp) const {
            timer t(world);
            t.do_print=do_print;
            Tensor<T> ovlp2=ovlp;
            for (int i=0; i<ovlp2.dim(0); ++i) ovlp2(i,i)-=1.0;
            if (world.rank()==0 and do_print) {
                print("absmax",ovlp2.absmax());
                print("l2",ovlp2.normf()/ovlp2.size());
            }
            t.tag("check_orthonoramality");
            return ovlp.absmax();
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

            // \int f(1,2) pre(1) post(2) \sum_i g(1) h(2) d1d2
//            double term2=madness::inner(pre*g,f12(post*h));
            double term2=madness::inner(g,inner(lrfunctor1,h,p2,p1));
            t.tag("computing term2");

            // g functions are orthonormal
            // \int gh(1,2)^2 d1d2 = \int \sum_{ij} g_i(1) g_j(1) h_i(2) h_j(2) d1d2
            //   = \sum_{ij} \int g_i(1) g_j(1) d1 \int h_i(2) h_j(2) d2
            //   = \sum_{ij} delta_{ij} \int h_i(2) h_j(2) d2
            //   = \sum_{i} \int h_i(2) h_i(2) d2
            double zero=check_orthonormality(g);
            if (zero>1.e-10) print("g is not orthonormal",zero);
            // double term3a=madness::inner(h,h);
            auto tmp1=matrix_inner(world,h,h);
            auto tmp2=matrix_inner(world,g,g);
            double term3=tmp1.trace(tmp2);
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

    // This interface is necessary to compute inner products
    template<typename T, std::size_t NDIM>
    double inner(const LowRankFunction<T,NDIM>& a, const LowRankFunction<T,NDIM>& b) {
        World& world=a.world;
        return (matrix_inner(world,a.g,b.g).emul(matrix_inner(world,a.h,b.h))).sum();
    }



//    template<typename T, std::size_t NDIM>
//    LowRankFunction<T,NDIM> inner(const Function<T,NDIM>& lhs, const LowRankFunction<T,NDIM>& rhs,
//                                  const std::tuple<int> v1, const std::tuple<int> v2) {
//        World& world=rhs.world;
//        // int lhs(1,2) rhs(2,3) d2 = \sum \int lhs(1,2) g_i(2) h_i(3) d2
//        //                      = \sum \int lhs(1,2) g_i(2) d2 h_i(3)
//        LowRankFunction<T, NDIM + NDIM - 2> result(world);
//        result.h=rhs.h;
//        decltype(rhs.g) g;
//        for (int i=0; i<rhs.rank(); ++i) {
//            g.push_back(inner(lhs,rhs.g[i],{v1},{0}));
//        }
//        result.g=g;
//        return result;
//    }

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
        auto result=inner(f2,f1,p2,p1);
        std::swap(result.g,result.h);
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
        // int f(1,2) k(2,3) d2 = \sum \int g_i(1) h_i(2) k(2,3) d2
        //                      = \sum g_i(1) \int h_i(2) k(2,3) d2
        LowRankFunction<T, NDIM> result(world);
        if (p1.is_last()) { // integrate over 2: result(1,3) = lrf(1,2) f(2,3)
            result.g = f1.g;
            change_tree_state(f1.h,reconstructed);
            result.h=innerXX<PDIM>(f2,f1.h,p2.get_array(),particle<PDIM>::particle1().get_array());
        } else if (p1.is_first()) { // integrate over 1: result(2,3) = lrf(1,2) f(1,3)
            result.g = f1.h;        // correct! second variable of f1 becomes first variable of result
            change_tree_state(f1.g,reconstructed);
            result.h=innerXX<PDIM>(f2,f1.g,p2.get_array(),particle<PDIM>::particle1().get_array());
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

        // inner(lrf(1,2) ,lrf(2,3) ) = \sum_ij g1_i(1) <h1_i(2) g2_j(2)> h2_j(3)
        auto matrix=matrix_inner(world,f2.get_functions(p2),f1.get_functions(p1));
        auto htilde=transform(world,f2.get_functions(p2.complement()),matrix);
        auto gg=copy(world,f1.get_functions(p1.complement()));
        return LowRankFunction<T,NDIM>(gg,htilde,f1.rank_revealing_tol,f1.orthomethod);
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

        // inner(lrf(1,2), f_k(2) ) = \sum_i g1_i(1) <h1_i(2) f_k(2)>
        auto matrix=matrix_inner(world,f1.get_functions(p1),vf);
        return transform(world,f1.get_functions(p1.complement()),matrix);
    }

    ///  f(1) = inner(lrf(1,2), f(2))

    /// @param[in] f1 the first function
    /// @param[in] vf the second function
    /// @param[in] p1 the integration variable of the first function
    /// @param[in] p2 the integration variable of the second function, dummy variable for consistent notation
    template<typename T, std::size_t NDIM, std::size_t PDIM>
    Function<T,NDIM> inner(const LowRankFunction<T,NDIM>& f1, const Function<T,PDIM>& f2,
                                        const particle<PDIM> p1, const particle<PDIM> p2=particle<PDIM>::particle1()) {
        return inner(f1,std::vector<Function<T,PDIM>>({f2}),p1,p2)[0];
    }

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

        LowRankFunction<T,NDIM> project(const LRFunctorBase<T,NDIM>& lrfunctor) const {
            World& world=lrfunctor.world();
            bool do_print=true;
            timer t1(world);
            t1.do_print=do_print;
            auto orthomethod=parameters.orthomethod();
            auto rank_revealing_tol=parameters.tol();

            // get sampling grid
            molecular_grid<LDIM> mgrid(origins,parameters);
            auto grid=mgrid.get_grid();
            if (world.rank()==0) print("grid size",grid.size());

            auto Y=Yformer(lrfunctor,grid,parameters.rhsfunctiontype());
            t1.tag("Yforming");

            auto ovlp=matrix_inner(world,Y,Y);  // error in symmetric matrix_inner, use non-symmetric form here!
            t1.tag("compute ovlp");
            auto g=truncate(orthonormalize_rrcd(Y,ovlp,rank_revealing_tol));
            t1.tag("rrcd/truncate/thresh");
            auto sz=get_size(world,g);
            if (world.rank()==0 and do_print) print("gsize",sz);
//            check_orthonormality(g);

            if (world.rank()==0 and do_print) {
                print("Y.size()",Y.size());
                print("g.size()",g.size());
            }

            auto h=truncate(inner(lrfunctor,g,p1,p1));
            t1.tag("Y backprojection with truncation");
            return LowRankFunction<T,NDIM>(g,h,parameters.tol(),parameters.orthomethod());

        }

        /// apply a rhs (delta or exponential) on grid points to the hi-dim function and form Y = A_ij w_j (in Halko's language)
        std::vector<Function<T,LDIM>> Yformer(const LRFunctorBase<T,NDIM>& lrfunctor1, const std::vector<Vector<double,LDIM>>& grid,
                                              const std::string rhsfunctiontype, const double exponent=30.0) const {

            World& world=lrfunctor1.world();
            std::vector<Function<double,LDIM>> Y;
            if (rhsfunctiontype=="exponential") {
                std::vector<Function<double,LDIM>> omega;
                double coeff=std::pow(2.0*exponent/constants::pi,0.25*LDIM);
                for (const auto& point : grid) {
                    omega.push_back(FunctionFactory<double,LDIM>(world)
                                            .functor([&point,&exponent,&coeff](const Vector<double,LDIM>& r)
                                                     {
                                                         auto r_rel=r-point;
                                                         return coeff*exp(-exponent*madness::inner(r_rel,r_rel));
                                                     }));
                }
                Y=inner(lrfunctor1,omega,p2,p1);
            } else {
                MADNESS_EXCEPTION("confused rhsfunctiontype",1);
            }
            auto norms=norm2s(world,Y);
            std::vector<Function<double,LDIM>> Ynormalized;

            for (size_t i=0; i<Y.size(); ++i) if (norms[i]>parameters.tol()) Ynormalized.push_back(Y[i]);
            normalize(world,Ynormalized);
            return Ynormalized;
        }

    };


} // namespace madness

#endif //MADNESS_LOWRANKFUNCTION_H
