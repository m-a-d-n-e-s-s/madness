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
        static constexpr char const *tag = "lrf";

        LowRankFunctionParameters() : QCCalculationParametersBase()
        {
            set_defaults();
        }

        LowRankFunctionParameters(World& world, const commandlineparser& parser) : QCCalculationParametersBase()
        {
            set_defaults();
            read_and_set_derived_values(world,parser,"lrf");
        }


        void set_defaults() {
            // initialize with: key, value, comment (optional), allowed values (optional)
            initialize<double>("radius",2.0,"the radius");
            initialize<double>("gamma",1.0,"the exponent of the correlation factor");
            initialize<double>("volume_element",0.1,"initial volume covered by each grid point");
            initialize<double>("tol",1.e-8,"rank-reduced cholesky tolerance");
            initialize<std::string>("f12type","Slater","correlation factor",{"Slater","SlaterF12"});
            // orthomethod removed — always use cholesky (rr_cholesky) orthonormalization
            initialize<std::string>("transpose","slater2","transpose of the matrix",{"slater1","slater2"});
            initialize<std::string>("gridtype","random","the grid type",{"random","adaptive","twostage","harmonics"});
            initialize<int>("optimize",1,"number of optimization iterations");
            initialize<int>("lmax",2,"max angular momentum for the RI harmonics");
            initialize<bool>("canonicalize",false,"canonicalize the rep, i.e. metric is the identity (non-canon. sensitive to thresh!)");
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
        // orthomethod removed — always cholesky
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

            // A point generated around centers[ic] is kept iff centers[ic] is its
            // nearest center — i.e. it lies in the Voronoi cell of its generator.
            // This eliminates the over-density that arises when neighbouring
            // per-center Gaussian clouds overlap. With a single center the test
            // is vacuous.
            const bool do_voronoi = (centers.size() > 1);
            auto in_own_voronoi_cell = [&](const Vector<double,NDIM>& p, const std::size_t ic) {
                double d_own = (p - centers[ic]).normf();
                for (std::size_t jc = 0; jc < centers.size(); ++jc) {
                    if (jc == ic) continue;
                    if ((p - centers[jc]).normf() < d_own) return false;
                }
                return true;
            };

            std::vector<Vector<double,NDIM>> grid;
            for (std::size_t ic = 0; ic < centers.size(); ++ic) {
                const auto& coords = centers[ic];
                print("atom sites",coords);
                std::vector<Vector<double,NDIM>> atomgrid;
                if (auto builder=dynamic_cast<dftgrid<NDIM>*>(grid_builder.get())) {
                    if constexpr (NDIM==3) {
                        dftgrid<NDIM> b1(builder->builder.get_nradial(),builder->builder.get_angular_order(),coords);
                        atomgrid=b1.get_grid();
                    } else {
                        MADNESS_EXCEPTION("no DFT grid for NDIM /= 3",1);
                    }
                } else if (auto builder=dynamic_cast<randomgrid<NDIM>*>(grid_builder.get())) {
                    randomgrid<NDIM> b1(builder->get_volume_element(),builder->get_radius(),coords);
                    atomgrid=b1.get_grid();
                } else {
                    MADNESS_EXCEPTION("no such grid builder",1);
                }

                if (do_voronoi) {
                    std::size_t kept = 0;
                    for (const auto& p : atomgrid) {
                        if (in_own_voronoi_cell(p, ic)) {
                            grid.push_back(p);
                            ++kept;
                        }
                    }
                    print("center",ic,"generated",atomgrid.size(),"kept",kept,"after Voronoi filter");
                } else {
                    grid.insert(grid.end(),atomgrid.begin(),atomgrid.end());
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

    /// set the MRA threshold on all internal functions and reproject to the new precision
    virtual void set_thresh(double thresh) {}

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

    /// access the F12 operator
    const std::shared_ptr<SeparatedConvolution<T,LDIM>>& get_f12() const { return f12; }
    /// access the a functions (particle 1 pre-multipliers)
    const std::vector<Function<T,LDIM>>& get_a() const { return a; }
    /// access the b functions (particle 2 pre-multipliers)
    const std::vector<Function<T,LDIM>>& get_b() const { return b; }

    void set_thresh(double thresh) override {
        for (auto& f : a) if (f.is_initialized()) f.set_thresh(thresh);
        for (auto& f : b) if (f.is_initialized()) f.set_thresh(thresh);
    }
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
        // Prefer the specialized combine rule when available; fall back to
        // the generic outer-product construction for SepConvs built from
        // raw (c, α) tensors whose OperatorInfo type is UNDEFINED.
        auto build_f12sq = [&]() -> SeparatedConvolution<T,LDIM> {
            try {
                return SeparatedConvolution<T,LDIM>::combine(f12a, f12a);
            } catch (const madness::MadnessException&) {
                MADNESS_CHECK_THROW(f12a.stored_coeffs.size() > 0,
                    "LRFunctorF12::norm2: cannot combine f12 with itself and no "
                    "stored (c, α) available for generic fallback");
                return SeparatedConvolution<T,LDIM>::combine_generic(
                    world(), f12a.stored_coeffs, f12a.stored_exponents,
                             f12a.stored_coeffs, f12a.stored_exponents,
                    f12a.info.lo, f12a.info.thresh, /*prune=*/true);
            }
        };
        const SeparatedConvolution<T,LDIM> f12sq = build_f12sq();

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

    void set_thresh(double thresh) override {
        if (f.is_initialized()) f.set_thresh(thresh);
    }

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
        bool do_print=false;
        std::vector<Function<T,LDIM>> g,h;
        Tensor<T> metric;   ///< the coupling matrix in the general representation, empty in the canonical representation
        const particle<LDIM> p1=particle<LDIM>::particle1();
        const particle<LDIM> p2=particle<LDIM>::particle2();

        LowRankFunction(World& world) : world(world) {}

        LowRankFunction(std::vector<Function<T,LDIM>> g, std::vector<Function<T,LDIM>> h,
                        double tol, const Tensor<T> metric=Tensor<T>()) : world(g.front().world()),
                        rank_revealing_tol(tol), g(g), h(h), metric(metric) {
        }

        /// shallow copy ctor
        LowRankFunction(const LowRankFunction& other) : world(other.world),
            rank_revealing_tol(other.rank_revealing_tol),
            g(other.g), h(other.h), metric(other.metric) {
        }

        /// deep copy
        friend LowRankFunction copy(const LowRankFunction& other) {
            return LowRankFunction<T,NDIM>(madness::copy(other.g),madness::copy(other.h),
                other.rank_revealing_tol, madness::copy(other.metric));
        }

        LowRankFunction& operator=(const LowRankFunction& f) { // Assignment required for storage in vector
            if (this == &f) return *this;
            rank_revealing_tol = f.rank_revealing_tol;
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
                if (ametric) tmp_metric(Slice(0,g.size()-1),Slice(0,h.size()-1))=ametric;

                Tensor<T> bmetric= (b.metric) ? b.metric : identity_matrix<T>(b.g.size());
                if (bmetric) tmp_metric(Slice(g.size(),g.size()+b.g.size()-1),Slice(h.size(),h.size()+b.h.size()-1))=bmetric;
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
            return LowRankFunction<TensorResultType<T,Q>,NDIM>(g * a, Q(h),rank_revealing_tol,metric);
        }

        /// out-of-place scale by a scalar (no type conversion)
        LowRankFunction operator*(const T a) const {
            return LowRankFunction(g * a, h,rank_revealing_tol,metric);
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
        /// orthonormalize the argument vector via rank-revealing Cholesky
        std::vector<Function<T,LDIM>> orthonormalize(const std::vector<Function<T,LDIM>>& g) const {
            MADNESS_CHECK_THROW(is_canonical(),"no orthonormalization unless canonicalized");
            auto ovlp=matrix_inner(world,g,g);
            auto g2=orthonormalize_rrcd(g,ovlp,rank_revealing_tol);
            double tight_thresh=FunctionDefaults<LDIM>::get_thresh()*0.1;
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

        /// Cheap relative L2 error via Pythagoras (canonical form only).
        ///
        /// When h is orthonormal AND g is the reprojection of h
        /// (g_k = <f(r1,r2) | h_k(r2)>_{r2}), the following identity holds:
        ///     ||f - sum_k g_k h_k||^2 = ||f||^2 - sum_k ||g_k||^2
        /// (cross and self terms collapse to the same sum_k ||g_k||^2).
        /// No functor applications are needed — only norms of the stored g.
        ///
        /// @param[in] f_norm_sq  precomputed ||f||^2 (typically lrfunctor.norm2()^2)
        /// @pre the LRF is in canonical form and g holds the reprojected factors.
        double l2error_pythagoras(double f_norm_sq) const {
            MADNESS_CHECK_THROW(is_canonical(),
                "l2error_pythagoras requires canonical form (h orthonormal)");
            MADNESS_CHECK_THROW(f_norm_sq > 0.0,
                "l2error_pythagoras requires a positive ||f||^2");
            auto gnorms = norm2s(world, g);
            double sum_gk_sq = 0.0;
            for (double n : gnorms) sum_gk_sq += n * n;
            double arg = f_norm_sq - sum_gk_sq;
            if (arg < 0.0) arg = -arg;  // numerical noise near convergence
            print("l2 error via Pythagoras: (abs., rel.)", sqrt(arg),sqrt(arg) / sqrt(f_norm_sq));
            return sqrt(arg) / sqrt(f_norm_sq);
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
        auto tmp = (f1.metric.size()>0) ? inner(f1.metric,matrix,index1,0) : matrix;
        auto metric = (f2.metric.size()>0) ? inner(tmp,f2.metric,1,index2) : tmp;
        auto gg=copy(world,f1.get_functions(p1.complement()));
        auto hh=copy(world,f2.get_functions(p2.complement()));
        return LowRankFunction<T,NDIM>(gg,hh,f1.rank_revealing_tol,metric);
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
        auto matrix = (f1.metric.size()>0) ? inner(f1.metric,matrix_jk,index1,0) : matrix_jk;
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

    /// Factory class to compute a low-rank approximation of a given hi-dimensional function
    /// f(r1,r2) using randomized projection and rank-revealing Cholesky decomposition (RRCD).
    ///
    /// The approximation takes the form
    ///   f(r1,r2) ≈ Σ_{ij} g_i(r1) X_{ij} h_j(r2)
    /// where X is the "half-metric", g are (non-orthogonal) basis functions, and h are
    /// backprojected functions. In the canonical form (canonicalize=true), X is absorbed
    /// into g and the representation simplifies to Σ_i g_i(r1) h_i(r2).
    ///
    /// ## Error analysis
    ///
    /// The L2 error ||f - f_approx|| has three independent sources that add in quadrature:
    ///
    ///   ε_total² ≈ ε_sampling² + ε_truncation² + ε_conditioning²
    ///
    /// ### 1. Sampling error (grid coverage)
    ///
    /// The grid probes f's column space via localized Gaussians. If the grid misses
    /// regions where f has significant weight, those contributions cannot be represented.
    /// Controlled by: volume_element (ve), radius, center placement.
    ///
    /// Heuristic: for well-localized functions (Gaussians on atomic centers), ve ≈ 0.1
    /// with radius ≈ 2-3 suffices for errors down to ~1e-3. Finer grids increase
    /// Y-formation cost linearly (O(grid_size) functor applications) but the rank
    /// saturates quickly — most extra grid points are linearly dependent.
    ///
    /// ### 2. Truncation error (rank reduction)
    ///
    /// The RRCD with tolerance tol discards Y components with overlap eigenvalues < tol.
    /// For Gaussian kernels the eigenvalues decay rapidly, giving truncation error that
    /// scales roughly as tol^(0.3-0.5) (empirically). Setting tol = ε² gives truncation
    /// error ~ ε.
    ///
    /// ### 3. Conditioning (eliminated by power-iteration algorithm)
    ///
    /// The current algorithm orthonormalizes h = f·Y (not Y), so the overlap
    /// <h|h> has eigenvalues weighted by σ_k² (singular values of f).  Rank
    /// reduction on <h|h> is always well-conditioned — no half-metric or
    /// full-metric amplification.  The result is canonical by construction
    /// (no coupling metric), and conditioning is not an error source.
    ///
    /// ### Parameter heuristics for a target error ε
    ///
    ///   tol    = ε²           (truncation error ~ √tol ~ ε)
    ///   ve     ~ C · ε^(2/LDIM)  with C = 2/1/0.5 for LDIM=1/2/3
    ///
    /// ### Algorithm (project())
    ///
    /// 1. Build grid around origins (Gaussian-distributed, σ² = 0.5·radius)
    /// 2. Form Y basis: multi-scale Gaussian RHS at grid points
    /// 3. Power-iteration projection: Y → h_raw → rr_cholesky(<h|h>) → h_ortho → g
    /// 4. Tol sweep: incrementally tighten rr_cholesky tolerance, reuse cached <h|h>
    /// 5. If grid-limited: halve ve, rebuild Y, re-project from scratch
    /// 6. If thresh-limited: tighten MRA thresh, rebuild Y, re-project
    /// 7. L2 error via Pythagoras: ||f||² − Σ||g_k||² (no functor applications)
    ///
    /// The project() method automates parameter selection from a single
    /// target L2 error.
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
        // set_orthomethod removed — always use cholesky
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


        /// Numerically stable projection using power-iteration.
        ///
        /// this method orthonormalizes h = f·Y. The overlap <h|h> has eigenvalues
        /// weighted by σ_k² (singular values of f), making the rank reduction well-conditioned.
        /// The result is canonical by construction (no metric).
        ///
        /// Algorithm:
        ///   1. rr_cholesky(<Y|Y>, tol_coarse) → pY              (coarse rank reduction)
        ///   2. h_raw = inner(f, pY, p1, p1)                     (backproject to r2)
        ///   3. rr_cholesky(<h_raw|h_raw>, tol) → ph              (well-conditioned rank)
        ///   4. orthonormalize ph → h_ortho                       (via Cholesky factor)
        ///   5. g = inner(f, h_ortho, p2, p1)                     (reproject to r1)
        ///   6. return LowRankFunction(g, h_ortho)                (canonical)
        ///
        /// @param[in] lrfunctor  the functor
        /// @param[in] Y          the probe basis from Yformer
        /// @param[in] tol        the target Cholesky tolerance
        /// @return               canonical LowRankFunction
        LowRankFunction<T,NDIM> project_from_Y_stable(
            const LRFunctorBase<T,NDIM>& lrfunctor,
            const std::vector<Function<T,LDIM>>& Y,
            const double tol) const
        {
            World& world = lrfunctor.world();
            timer t1(world);
            t1.do_print = true;

            double tight_thresh = FunctionDefaults<LDIM>::get_thresh() * 0.1;

            // Step 1: coarse rank reduction on Y (cheap, removes exact lindep only)
            // Use TIGHT tolerance so we don't prematurely discard info; the
            // authoritative rank reduction happens on <h|h> below.
            double tol_coarse = std::max(1.e-14, tol * 0.01);
            auto ovlp_Y = matrix_inner(world, Y, Y);
            auto [pY, Ut_Y, Xinv_Y] = LowRankFunction<T,NDIM>::rr_cholesky_matrix_and_reorder(
                Y, ovlp_Y, tol_coarse);
            t1.tag("stable: rr_cholesky on Y (coarse)");
            print("stable: Y.size =", Y.size(), "pY.size =", pY.size());

            // Step 2: backproject pY to r2 space — h_raw inherits σ_k² conditioning
            auto h_raw = truncate(inner(lrfunctor, pY, p1, p1), tight_thresh);
            t1.tag("stable: backprojection (h_raw)");

            // Step 3: rank reduction on <h|h> — well-conditioned
            auto ovlp_h = matrix_inner(world, h_raw, h_raw);
            auto [ph, Ut_h, Xinv_h] = LowRankFunction<T,NDIM>::rr_cholesky_matrix_and_reorder(
                h_raw, ovlp_h, tol);
            t1.tag("stable: rr_cholesky on h");
            print("stable: h_raw.size =", h_raw.size(), "ph.size =", ph.size());

            // Step 4: orthonormalize ph → h_ortho via Cholesky factor X = U^{-1}
            // Since <h|h> is well-conditioned, this orthonormalization is stable
            auto h_ortho = truncate(transform(world, ph, Ut_h), tight_thresh);
            t1.tag("stable: orthonormalization");

            // Step 5: reproject to r1 space — g lives in r1 and is well-behaved
            auto g = truncate(inner(lrfunctor, h_ortho, p2, p1), tight_thresh);
            t1.tag("stable: reprojection (g)");

            // Step 6: canonical form, no metric
            LowRankFunction<T,NDIM> result(g, h_ortho, tol);
            print("stable: final rank =", g.size());

            return result;
        }

        /// State for the stable power-iteration projection.
        /// Stores h_raw and its overlap for incremental tol-tightening on <h|h>.
        struct ProjectionState {
            std::vector<Function<T,LDIM>> h_raw;    ///< backprojected h (before orthonorm)
            std::vector<Function<T,LDIM>> h_ortho;  ///< orthonormalized h
            std::vector<Function<T,LDIM>> g;        ///< reprojected g
            Tensor<T> ovlp_h;                       ///< cached <h_raw|h_raw>
            double tol = 0.0;
            double f_norm_sq = -1.0;                ///< cached ||f||^2 for Pythagoras l2error

            long rank() const { return g.size(); }
        };

        /// Build a grid around the given origins (or uniform random if no origins).
        std::vector<Vector<double,LDIM>> build_grid(double ve, long max_points = 0) const {
            std::vector<Vector<double,LDIM>> grid;
            if (!origins.empty()) {
                for (const auto& origin : origins) {
                    randomgrid<LDIM> rg(ve, parameters.radius(), origin);
                    auto local = rg.get_grid();
                    grid.insert(grid.end(), local.begin(), local.end());
                }
            } else {
                grid = make_uniform_random_grid_in_cell(ve);
            }
            if (max_points > 0 && long(grid.size()) > max_points)
                grid.resize(max_points);
            return grid;
        }

        /// Form Y basis from grid, with fallback to loose significance threshold.
        std::vector<Function<T,LDIM>> form_Y(
            const LRFunctorBase<T,NDIM>& lrfunctor,
            const std::vector<Vector<double,LDIM>>& grid,
            const LowRankFunctionParameters& params) const
        {
            auto yresult = Yformer(lrfunctor, grid, params);
            auto Y = yresult.Y;
            if (Y.empty()) {
                yresult = Yformer(lrfunctor, grid, params, 30.0, 0.0);
                Y = yresult.Y;
            }
            return Y;
        }

        /// RAII guard for temporarily changing FunctionDefaults and functor thresh.
        struct ThreshGuard {
            double orig_ld, orig_nd;
            LRFunctorBase<T,NDIM>& functor;
            ThreshGuard(double tight, LRFunctorBase<T,NDIM>& f)
                : orig_ld(FunctionDefaults<LDIM>::get_thresh()),
                  orig_nd(FunctionDefaults<NDIM>::get_thresh()), functor(f) {
                FunctionDefaults<LDIM>::set_thresh(tight);
                FunctionDefaults<NDIM>::set_thresh(tight);
                functor.set_thresh(tight);
            }
            ~ThreshGuard() {
                FunctionDefaults<LDIM>::set_thresh(orig_ld);
                FunctionDefaults<NDIM>::set_thresh(orig_nd);
                functor.set_thresh(orig_ld);
            }
        };

        /// Stable projection using power iteration.
        ///
        /// If prev state is empty, does full projection (Y → h → ortho → g).
        /// If prev has content, incrementally extends: rr_cholesky on ovlp_h with
        /// tighter tol, orthonormalize and reproject only the new functions.
        ///
        /// The overlap matrix ovlp_h is cached and reused across calls for pivot
        /// consistency (parallel matrix_inner is non-deterministic).
        ProjectionState project_stable(
            const LRFunctorBase<T,NDIM>& lrfunctor,
            const std::vector<Function<T,LDIM>>& Y,
            const double tol,
            const ProjectionState& prev) const
        {
            World& world = lrfunctor.world();
            timer t1(world);
            t1.do_print = true;
            double tight_thresh = FunctionDefaults<LDIM>::get_thresh() * 0.1;

            ProjectionState state;
            state.tol = tol;
            // cache ||f||^2 for cheap Pythagoras l2error (canonical form)
            if (prev.f_norm_sq > 0.0) {
                state.f_norm_sq = prev.f_norm_sq;
            } else {
                double fn = lrfunctor.norm2();
                state.f_norm_sq = fn * fn;
                t1.tag("stable: ||f||^2");
            }

            if (prev.rank() == 0) {
                // --- full projection ---
                // coarse rank reduction on Y
                double tol_coarse = std::max(1.e-14, tol * 0.01);
                auto ovlp_Y = matrix_inner(world, Y, Y);
                auto [pY, Ut_Y, Xinv_Y] = LowRankFunction<T,NDIM>::rr_cholesky_matrix_and_reorder(
                    Y, ovlp_Y, tol_coarse);
                t1.tag("stable: rr_cholesky on Y");

                // backproject to h_raw
                state.h_raw = truncate(inner(lrfunctor, pY, p1, p1), tight_thresh);
                t1.tag("stable: backprojection (h_raw)");

                // cache <h|h> and do authoritative rank reduction
                state.ovlp_h = matrix_inner(world, state.h_raw, state.h_raw);
                auto [ph, Ut_h, Xinv_h] = LowRankFunction<T,NDIM>::rr_cholesky_matrix_and_reorder(
                    state.h_raw, copy(state.ovlp_h), tol);
                t1.tag("stable: rr_cholesky on h");

                state.h_ortho = truncate(transform(world, ph, Ut_h), tight_thresh);
                t1.tag("stable: orthonormalization");

                state.g = truncate(inner(lrfunctor, state.h_ortho, p2, p1), tight_thresh);
                t1.tag("stable: reprojection (g)");
                LowRankFunction<T,NDIM> lrf(state.g, state.h_ortho, tol);
                double current_error = lrf.l2error_pythagoras(state.f_norm_sq);
                t1.tag("stable: compute error in optimization");


                print("stable: pY =", pY.size(), "h_raw =", state.h_raw.size(),
                      "rank =", state.g.size());

            } else {
                // --- incremental tol-tightening on <h|h> ---
                // reuse prev.h_raw and prev.ovlp_h (must be same — cached)
                state.h_raw = prev.h_raw;
                state.ovlp_h = prev.ovlp_h;

                Tensor<T> ovlp_work = copy(state.ovlp_h);
                auto [ph, Ut_h, Xinv_h] = LowRankFunction<T,NDIM>::rr_cholesky_matrix_and_reorder(
                    state.h_raw, ovlp_work, tol);
                t1.tag("stable: rr_cholesky on h (incremental)");

                long r_prev = prev.rank();
                long r_new = ph.size();

                if (r_new <= r_prev) {
                    state.h_ortho = prev.h_ortho;
                    state.g = prev.g;
                } else {
                    // only orthonormalize and reproject the new functions
                    Tensor<T> X_new_cols = Ut_h(_, Slice(r_prev, r_new - 1));
                    auto h_ortho_new = truncate(transform(world, ph, X_new_cols), tight_thresh);
                    t1.tag("stable: orthonormalization (new only)");
                    auto g_new = truncate(inner(lrfunctor, h_ortho_new, p2, p1), tight_thresh);
                    t1.tag("stable: reprojection (new only)");

                    state.h_ortho = append(prev.h_ortho, h_ortho_new);
                    state.g = append(prev.g, g_new);
                }
                print("stable incremental: r_prev =", r_prev, "r_new =", state.g.size());
            }

            return state;
            }

        /// Diagnose the dominant error source.
        /// With the stable (canonical) algorithm, conditioning is not an issue.
        /// Thresh-limited only if:
        ///   - the last tol-tightening did not improve the error by at least 20%
        ///     (error >= 0.8 * prev_error), AND
        ///   - sqrt(rank) * thresh > 0.3 * error (realistic noise combination), AND
        ///   - tol is at the floor (~1e-14), i.e. tol-tightening has no more room.
        /// The sqrt(rank) factor replaces the worst-case rank factor: for canonical
        /// form with uncorrelated truncation noise, combined error scales like
        /// sqrt(rank) * thresh, not rank * thresh. The rank-linear bound triggers
        /// spuriously at moderate rank. The tol-floor guard prevents switching to
        /// expensive thresh work when grid augmentation hasn't been tried yet.
        /// Tol-limited: tol probe on ovlp_h shows rank increase (preferred when available).
        /// Grid-limited: neither of the above — try more grid points.
        std::string diagnose_error(double error, double prev_error, long rank,
                                   double tol, double thresh,
                                   const Tensor<T>& ovlp_h) const {
            double eps_trunc = std::sqrt(double(rank)) * thresh;
            bool error_stalled = (prev_error > 0 && error >= prev_error * 0.8);
            bool tol_at_floor = (tol <= 1.01e-14);
            print("  diagnosis: sqrt(rank)*thresh =", eps_trunc, "error =", error,
                  "rank =", rank, "tol =", tol,
                  "stalled =", error_stalled, "tol_at_floor =", tol_at_floor);

            // tol probe on <h|h> — check first, prefer tol if available
            double tol_probe = std::max(1.e-14, tol * 0.01);
            if (tol_probe < tol && ovlp_h.size() > 0) {
                Tensor<T> ovlp_copy = copy(ovlp_h);
                int rank_probe;
                Tensor<integer> piv_probe;
                rr_cholesky(ovlp_copy, tol_probe, piv_probe, rank_probe);
                print("  tol probe on <h|h>: rank", rank, "->", rank_probe);
                if (rank_probe > rank * 1.1) {
                    print("  -> tol-limited");
                    return "tol";
                }
            }

            // thresh: only when error stalled AND realistic noise bound is significant
            // AND we've already exhausted tol headroom (else grid should be tried first).
            if (error_stalled && tol_at_floor && eps_trunc > 0.3 * error) {
                print("  -> thresh-limited");
                return "thresh";
            }
            print("  -> grid-limited");
            return "grid";
        }

        /// Adaptively project a high-dimensional functor to a low-rank representation.
        ///
        /// Uses the numerically stable power-iteration algorithm:
        ///   Y → h = f·Y → orthonormalize h → g = f·h_ortho → canonical LRF
        ///
        /// Single refinement loop with three actions:
        ///   - tol-limited: tighten tol on <h|h>, incrementally extend h_ortho and g
        ///   - grid-limited: halve ve and re-project from scratch (denser Y sampling)
        ///   - thresh-limited: tighten MRA thresh, reform everything
        ///   - tol-degraded: revert to best state, route to thresh or grid
        ///
        /// @param[in] lrfunctor       the high-dimensional functor to approximate
        /// @param[in] target_l2error  the desired L2 error (default: FunctionDefaults thresh)
        /// @param[in] max_iter        maximum refinement iterations (default: 5)
        /// @return                    canonical low-rank approximation
        LowRankFunction<T,NDIM> project(
            const LRFunctorBase<T,NDIM>& lrfunctor,
            const double target_l2error = FunctionDefaults<NDIM>::get_thresh(),
            const int max_iter = 5) const
        {
            World& world = lrfunctor.world();
            timer t1(world);
            t1.do_print = true;
            const double eps = target_l2error;
            const long max_Y_size = 5000;

            double tol = std::max(1.e-14, std::min(1.e-3, eps * eps));
            double current_thresh = FunctionDefaults<LDIM>::get_thresh();
            double ve = parameters.volume_element();

            print("project: eps =", eps, "tol =", tol,
                  "ve =", ve, "thresh =", current_thresh);
            t1.tag("project: parameter setup");

            LowRankFunctionParameters local_params = parameters;
            local_params.set_derived_value("volume_element", ve);
            local_params.set_derived_value("tol", tol);

            auto grid = build_grid(ve);
            print("project: initial grid size", grid.size());
            t1.tag("project: grid construction");

            auto Y = form_Y(lrfunctor, grid, local_params);
            MADNESS_CHECK_THROW(!Y.empty(), "project: no basis functions");
            print("project: initial Y size", Y.size());
            t1.tag("project: Yforming");

            // initial stable projection
            ProjectionState empty_state;
            auto state = project_stable(lrfunctor, Y, tol, empty_state);
            auto lrf = LowRankFunction<T,NDIM>(state.g, state.h_ortho,
                tol);
            t1.tag("project: initial projection");

            double current_error = lrf.l2error_pythagoras(state.f_norm_sq);
            long current_rank = state.rank();
            print("project: iter 0, l2error =", current_error,
                  "rank =", current_rank, "target =", eps);
            t1.tag("project: l2error");

            if (current_error <= eps) {
                print("project: converged at iteration 0");
                return lrf;
            }

            // ThreshGuard for optional thresh-tightening
            std::unique_ptr<ThreshGuard> thresh_guard;
            double prev_error = -1.0;
            int iter = 0;

            // track best state seen across all iterations
            auto best_lrf = lrf;
            double best_error = current_error;

            while (iter <= max_iter) {
                if (state.ovlp_h.size() == 0) state.ovlp_h=matrix_inner(world,state.h_ortho,state.h_ortho);
                auto diagnosis = diagnose_error(current_error, prev_error, current_rank,
                                                tol, current_thresh, state.ovlp_h);
                prev_error = current_error;

                // --- tol-limited: try incremental tol-tightening ---
                if (diagnosis == "tol") {
                    auto prev_state = state;
                    auto prev_lrf = lrf;
                    double prev_tol = tol;

                    tol = std::max(1.e-14, tol * 0.01);
                    state = project_stable(lrfunctor, Y, tol, prev_state);
                    lrf = LowRankFunction<T,NDIM>(state.g, state.h_ortho,
                        tol);
                    t1.tag("project: tol refinement iter " + std::to_string(iter));

                    current_error = lrf.l2error_pythagoras(state.f_norm_sq);
                    current_rank = state.rank();
                    print("project: iter", iter, "(tol), l2error =",
                          current_error, "rank =", current_rank);
                    t1.tag("project: l2error iter " + std::to_string(iter));
                    if (current_error <= eps) return lrf;
                    if (current_error < best_error) {
                        best_error = current_error;
                        best_lrf = lrf;
                    }

                    if (current_error > prev_error * 1.1) {
                        // tol made it worse — revert. Tol-degradation is a strong
                        // thresh-limit signal: we're trying to resolve singular values
                        // smaller than the current MRA thresh can represent, so the
                        // new directions are numerical noise. Route to thresh-tightening
                        // (if not already done); otherwise grid augmentation is the
                        // only remaining avenue.
                        print("project: tol degraded, reverting");
                        state = prev_state;
                        lrf = prev_lrf;
                        tol = prev_tol;
                        current_error = prev_error;
                        current_rank = state.rank();
                        diagnosis = thresh_guard ? "grid" : "thresh";
                        // NB: do NOT ++iter here — the branch below does so
                    } else {
                        continue;  // tol succeeded — don't count as iteration
                    }
                }

                // --- grid-limited (or tol-degraded fallthrough): fresh Y at finer ve ---
                if (diagnosis == "grid") {
                    ++iter;
                    if (iter>max_iter) break;
                    ve *= 0.5;
                    print("project: grid-limited, halving ve to", ve, "and re-projecting from scratch");
                    local_params.set_derived_value("volume_element", ve);

                    auto new_grid = build_grid(ve, max_Y_size);
                    if (new_grid.empty()) break;
                    print("project: new grid size", new_grid.size());

                    Y = form_Y(lrfunctor, new_grid, local_params);
                    t1.tag("project: grid Yforming iter " + std::to_string(iter));
                    if (Y.empty()) break;
                    print("project: new Y size", Y.size());

                    // reset tol and prev_error for a fresh tol sweep
                    tol = std::max(1.e-14, std::min(1.e-3, eps * eps));
                    prev_error = -1.0;

                    ProjectionState fresh;
                    state = project_stable(lrfunctor, Y, tol, fresh);
                    lrf = LowRankFunction<T,NDIM>(state.g, state.h_ortho,
                        tol);
                    t1.tag("project: grid projection iter " + std::to_string(iter));

                    current_error = lrf.l2error_pythagoras(state.f_norm_sq);
                    current_rank = state.rank();
                    print("project: iter", iter, "(grid ve=", ve, "), l2error =",
                          current_error, "rank =", current_rank);
                    t1.tag("project: l2error iter " + std::to_string(iter));
                    if (current_error <= eps) return lrf;

                    // stagnation: grid halving didn't beat the previous best
                    if (current_error >= best_error * 0.9) {
                        print("project: grid halving stagnated (best=", best_error,
                              "), stopping");
                        break;
                    }
                    if (current_error < best_error) {
                        best_error = current_error;
                        best_lrf = lrf;
                    }
                    continue;
                }

                // --- thresh-limited: tighten thresh, reform everything ---
                if (diagnosis == "thresh") {
                    ++iter;
                    if (iter>max_iter) break;
                    if (thresh_guard) {
                        print("project: thresh already tightened, giving up");
                        break;
                    }
                    current_thresh *= 0.1;
                    print("project: tightening thresh to", current_thresh);
                    thresh_guard = std::make_unique<ThreshGuard>(current_thresh,
                        const_cast<LRFunctorBase<T,NDIM>&>(lrfunctor));
                    t1.tag("project: thresh tightening");

                    Y = form_Y(lrfunctor, grid, local_params);
                    if (Y.empty()) break;
                    t1.tag("project: tight Yforming");

                    tol = std::max(1.e-14, std::min(1.e-3, eps * eps));
                    ProjectionState empty;
                    state = project_stable(lrfunctor, Y, tol, empty);
                    lrf = LowRankFunction<T,NDIM>(state.g, state.h_ortho,
                        tol);
                    t1.tag("project: tight projection");

                    current_error = lrf.l2error_pythagoras(state.f_norm_sq);
                    current_rank = state.rank();
                    print("project: iter", iter, "(thresh), l2error =",
                          current_error, "rank =", current_rank);
                    t1.tag("project: l2error iter " + std::to_string(iter));
                    if (current_error <= eps) return lrf;
                }
            }

            print("WARNING project: did NOT converge to target", eps,
                  "; achieved", best_error);
            return best_lrf;
        }

    private:
        /// Build a Taylor-expansion-based Y basis spanning the range of an F12 operator.
        ///
        /// For each Gaussian term c_m exp(-alpha_m r^2) in the operator's GFit
        /// expansion and each multi-index k with |k| <= max_order, emits
        ///     phi_{m,k}(r) = r^k * exp(-alpha_m |r|^2)
        /// when the Taylor coefficient |c_m * prod_d (2 alpha_m)^{k_d}/k_d!|
        /// exceeds coeff_cutoff. These monomial-Gaussians probe the operator's
        /// range analytically rather than by random sampling.
        std::vector<Function<T,LDIM>> build_taylor_y(
            const LRFunctorF12<T,NDIM>& functor,
            const int max_order,
            const double coeff_cutoff) const
        {
            World& world = functor.world();
            auto f12 = functor.get_f12();
            GFit<T, LDIM> fit(f12->info);
            Tensor<T> coeffs = fit.coeffs();
            Tensor<T> expnts = fit.exponents();
            const long M = coeffs.dim(0);

            auto multi_indices = generate_multi_indices(max_order);

            std::vector<Function<T,LDIM>> Y;
            long skipped = 0;
            for (long m = 0; m < M; ++m) {
                double cm = coeffs(m);
                double am = expnts(m);
                for (const auto& kidx : multi_indices) {
                    double coeff_mk = cm;
                    for (std::size_t d = 0; d < LDIM; ++d) {
                        for (int j = 0; j < kidx[d]; ++j) {
                            coeff_mk *= (2.0 * am) / double(j + 1);
                        }
                    }
                    if (std::abs(coeff_mk) < coeff_cutoff) { ++skipped; continue; }

                    std::array<int, LDIM> k_copy = kidx;
                    Function<T,LDIM> phi = FunctionFactory<T,LDIM>(world)
                        .special_points({Vector<double,LDIM>(0.0)})
                        .functor([k_copy, am](const Vector<double,LDIM>& r) -> T {
                            double val = std::exp(-am * madness::inner(r, r));
                            for (std::size_t d = 0; d < LDIM; ++d) {
                                for (int p = 0; p < k_copy[d]; ++p) val *= r[d];
                            }
                            return T(val);
                        });
                    Y.push_back(phi);
                }
            }
            print("build_taylor_y: M =", M, "max_order =", max_order,
                  "Y =", Y.size(), "(skipped", skipped, "below cutoff", coeff_cutoff, ")");
            return Y;
        }

        /// Generate all multi-indices k = (k_1, ..., k_LDIM) with total degree |k| <= max_order.
        static std::vector<std::array<int, LDIM>> generate_multi_indices(int max_order) {
            std::vector<std::array<int, LDIM>> result;
            std::array<int, LDIM> k{};
            std::function<void(std::size_t, int)> recurse = [&](std::size_t dim, int remaining) {
                if (dim == LDIM) {
                    result.push_back(k);
                    return;
                }
                for (int ki = 0; ki <= remaining; ++ki) {
                    k[dim] = ki;
                    recurse(dim + 1, remaining - ki);
                }
            };
            recurse(0, max_order);
            return result;
        }

    public:
        /// Construct LRF directly from the Gaussian expansion of an F12 operator.
        ///
        /// Instead of random-Y probing, uses the known separable structure:
        /// each Gaussian term c_m exp(-alpha_m |r1-r2|^2) is Taylor-expanded via
        ///     exp(2 alpha_m r1.r2) = sum_k (2 alpha)^|k| r1^k r2^k / k!
        /// giving explicit g and h basis functions. After rr_cholesky rank
        /// reduction, yields a canonical LRF accurate to the operator's Gaussian
        /// fit precision (~1e-8), bypassing the O(1/r^{3/2}) error floor of
        /// random probing for Slater-type kernels.
        ///
        /// @param[in] functor          the F12 functor with operator and a,b functions
        /// @param[in] target_error     controls coefficient cutoff and rank reduction tol
        /// @param[in] max_taylor_order max total degree of Taylor expansion per Gaussian term
        /// @return                     canonical LowRankFunction
        LowRankFunction<T,NDIM> project_from_operator(
            const LRFunctorF12<T,NDIM>& functor,
            const double target_error = FunctionDefaults<NDIM>::get_thresh(),
            const int max_taylor_order = 15) const
        {
            World& world = functor.world();
            timer t1(world);
            t1.do_print = true;

            double coeff_cutoff = target_error * target_error * 1.e-4;
            auto Y = build_taylor_y(functor, max_taylor_order, coeff_cutoff);
            t1.tag("direct: Y construction");

            if (Y.empty()) {
                print("WARNING project_from_operator: no significant Y terms");
                return LowRankFunction<T,NDIM>(world);
            }

            // 4. Use project_stable with the Taylor-constructed Y basis
            //    This lets the SeparatedConvolution handle coefficient summation
            //    correctly: h_raw = K*(a*Y) sums over ALL Gaussian terms, preserving
            //    cancellations that the pure Taylor construction loses.
            double tol = std::max(1.e-14, std::min(1.e-3, target_error * target_error));
            ProjectionState empty;
            auto state = project_stable(functor, Y, tol, empty);
            auto lrf = LowRankFunction<T,NDIM>(state.g, state.h_ortho, tol);

            double current_error = lrf.l2error_pythagoras(state.f_norm_sq);
            print("direct: Pythagoras error =", current_error, "rank =", state.rank());
            t1.tag("direct: projection");

            // incremental tol-tightening if needed
            while (current_error > target_error && state.ovlp_h.size() > 0) {
                double tighter = std::max(1.e-14, tol * 0.01);
                if (tighter >= tol) break;
                auto state2 = project_stable(functor, Y, tighter, state);
                auto lrf2 = LowRankFunction<T,NDIM>(state2.g, state2.h_ortho, tighter);
                double error2 = lrf2.l2error_pythagoras(state2.f_norm_sq);
                print("direct: tol", tol, "->", tighter, ": error =", error2, "rank =", state2.rank());
                if (error2 >= current_error * 0.9) break;  // stagnated
                lrf = lrf2;
                state = state2;
                current_error = error2;
                tol = tighter;
                t1.tag("direct: tol tightening");
            }

            print("direct: final error =", current_error, "rank =", lrf.rank());
            return lrf;
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
                Y = harmonic_basis(world, parameters.tempered(), parameters.lmax(), origins);
            } else {
                // default: use localized Gaussian RHS (for gridtype values other than "harmonics")
                std::vector<Function<double,LDIM>> omega;
                for (const auto& point : grid) {
                    double coeff=std::pow(2.0*exponent/constants::pi,0.25*LDIM);
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
