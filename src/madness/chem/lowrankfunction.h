//
// Created by Florian Bischoff on 8/10/23.
//

#ifndef MADNESS_LOWRANKFUNCTION_H
#define MADNESS_LOWRANKFUNCTION_H


#include<madness/mra/mra.h>
#include<madness/mra/vmra.h>
#include<madness/world/timing_utilities.h>
#include<madness/chem/electronic_correlation_factor.h>
#include <random>


namespace madness {

    struct LowRankFunctionParameters : QCCalculationParametersBase {

        LowRankFunctionParameters() : QCCalculationParametersBase() {

            // initialize with: key, value, comment (optional), allowed values (optional)
            initialize<double>("radius",2.0,"the radius");
            initialize<double>("gamma",1.0,"the exponent of the correlation factor");
            initialize<double>("volume_element",0.1,"volume covered by each grid point");
            initialize<bool>("hard_shell",true,"radius is hard");
            initialize<double>("tol",1.e-8,"rank-reduced cholesky tolerance");
            initialize<std::string>("f12type","Slater","correlation factor",{"Slater","SlaterF12"});
            initialize<std::string>("orthomethod","cholesky","orthonormalization",{"cholesky","canonical","symmetric"});
            initialize<std::string>("transpose","slater2","transpose of the matrix",{"slater1","slater2"});
            initialize<std::string>("gridtype","random","the grid type",{"random","cartesian","spherical"});
            initialize<std::string>("rhsfunctiontype","exponential","the type of function",{"delta","exponential"});
            initialize<int>("optimize",1,"number of optimization iterations");
        }

        void read_and_set_derived_values(World& world, const commandlineparser& parser, std::string tag) {
            read_input_and_commandline_options(world,parser,tag);
        }

        double radius() const {return get<double>("radius");}
        double gamma() const {return get<double>("gamma");}
        double volume_element() const {return get<double>("volume_element");}
        double tol() const {return get<double>("tol");}
        bool hard_shell() const {return get<bool>("hard_shell");}
        int optimize() const {return get<int>("optimize");}
        std::string gridtype() const {return get<std::string>("gridtype");}
        std::string orthomethod() const {return get<std::string>("orthomethod");}
        std::string rhsfunctiontype() const {return get<std::string>("rhsfunctiontype");}
        std::string f12type() const {return get<std::string>("f12type");}
    };



    template<std::size_t NDIM>
    class cgrid {
    public:
        double volume_element=0.1;
        double radius=3;
        bool hard_shell=true;

        cgrid(const double volume_element, const double radius, bool hard_shell)
                :volume_element(volume_element), radius(radius), hard_shell(hard_shell) {
        };

        std::vector<Vector<double,NDIM>> get_coordinates() const {
            // 1D grid
            double volume_element_1d=std::pow(volume_element,1./NDIM);
            long ngrid=std::ceil(radius/volume_element_1d);
            double stepsize=radius/ngrid;
            double scale=1.0;
            if (not hard_shell) scale=std::pow(2.0,1.0/(ngrid+1));
            print("scale",scale);

            std::vector<double> coord1d;
            print("volume element, stepsize, ngrid" ,volume_element, std::pow(stepsize,NDIM),stepsize,ngrid);
            for (int i=0; i<ngrid+1; ++i) {
                if (not hard_shell) stepsize*=scale;
                double c=i*stepsize;
                coord1d.push_back(c);
                if (i!=0) coord1d.push_back(-c);
                print("coord",c);
            }
            print("coord1d",coord1d.size());
            std::vector<Vector<double,NDIM>> result;
            for (int i=0; i<coord1d.size(); ++i) {
                for (int j=0; j<coord1d.size(); ++j) {
                    for (int k=0; k<coord1d.size(); ++k) {
                        Vector<double,NDIM> c{coord1d[i],coord1d[j],coord1d[k]};
                        double cutoff = hard_shell ? radius : 2.0*radius;
                        if (c.normf()<cutoff) result.push_back(c);
                    }
                }
            }
            print("result.size()",result.size());
//        for (auto& r: result) print(r);
            return  result;

        }
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


    template<std::size_t NDIM>
    struct randomgaussian {
        Vector<double,NDIM> random_origin;
        double exponent;
        double radius=2;
        randomgaussian(double exponent, double radius) : exponent(exponent), radius(radius) {
//            Vector<double,NDIM> ran; // [0,1]
//            RandomVector(NDIM,ran.data());
            Vector<double,NDIM> ran= this->gaussian_random_distribution(0,radius);
            random_origin=2.0*radius*ran-Vector<double,NDIM>(radius);
//            print("origin at ",random_origin, ", exponent",exponent);
        }
        double operator()(const Vector<double,NDIM>& r) const {
    //        return exp(-exponent*inner(r-random_origin,r-random_origin));
            double r2=inner(r-random_origin,r-random_origin);
            return exp(-exponent*r2);
        }

        static Vector<double,NDIM> gaussian_random_distribution(double mean, double variance) {
            std::random_device rd{};
            std::mt19937 gen{rd()};
            std::normal_distribution<> d{mean, variance};
            Vector<double,NDIM> result;
            for (int i = 0; i < NDIM; ++i) result[i]=d(gen);
            return result;
        }
    };


    /// LowRankFunction represents a hi-dimensional (NDIM) function as a sum of products of low-dimensional (LDIM) functions

    /// f(1,2) = \sum_i g_i(1) h_i(2)
    /// a LowRankFunction can be created from a hi-dim function directly, or from a composite like f(1,2) phi(1) psi(2),
    /// where f(1,2) is a two-particle function (e.g. a Slater function)
    template<typename T, std::size_t NDIM, std::size_t LDIM=NDIM/2>
    class LowRankFunction {
    public:

        /// what the LowRankFunction will represent
        struct LRFunctor {
            LRFunctor() = default;

            Function<T, NDIM> f;    ///< a hi-dim function
            std::shared_ptr<SeparatedConvolution<T,LDIM>> f12;  ///< a two-particle function
            Function<T,LDIM> a,b;   ///< the lo-dim functions
            bool has_f() const {
                return f.is_initialized();
            }
            bool has_f12() const {
                return (f12.get());
            }
            T operator()(const Vector<double,NDIM>& r) const {

                if (f12->info.type==OT_SLATER) {
                    double gamma=f12->info.mu;
                    Vector<double,LDIM> first, second;
                    for (int i=0; i<LDIM; ++i) {
                        first[i]=r[i];
                        second[i]=r[i+LDIM];
                    }
                    MADNESS_CHECK(has_f12());
                    double result=a(first)*b(second)*exp(-gamma*(first-second).normf());
                    return result;
                } else {
                    return 1.0;
                }
            }
        };

        template<std::size_t PDIM>
        struct particle {
            std::array<int,PDIM> dims;

            /// default constructor
            particle() = default;

            /// convenience for particle 1 (the left/first particle)
            static particle particle1() {
                particle p;
                for (int i=0; i<PDIM; ++i) p.dims[i]=i;
                return p;
            }
            /// convenience for particle 2 (the right/second particle)
            static particle particle2() {
                particle p;
                for (int i=0; i<PDIM; ++i) p.dims[i]=i+PDIM;
                return p;
            }


            particle(const int p) : particle(std::vector<int>(1,p)) {}
            particle(const int p1, const int p2) : particle(std::vector<int>({p1,p2})) {}
            particle(const int p1, const int p2,const int p3) : particle(std::vector<int>({p1,p2,p3})) {}
            particle(const std::vector<int> p) {
                for (int i=0; i<PDIM; ++i) dims[i]=p[i];
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


        /// inner product: result(1) = \int f(1,2) rhs(2) d2

        /// @param[in]  functor the hidim function
        /// @param[in]  rhs     the rhs
        /// @param[in]  p1      the variable in f(1,2) to be integrated over
        /// @param[in]  p2      the variable in rhs to be integrated over (usually all of them)
        static std::vector<Function<T,LDIM>> inner(LRFunctor& functor, const std::vector<Function<T,LDIM>>& rhs,
                                                   const particle<LDIM> p1, const particle<LDIM> p2) {
            std::vector<Function<T,LDIM>> result;
            MADNESS_CHECK(functor.has_f() xor functor.has_f12());
            MADNESS_CHECK(p1.is_first() xor p1.is_last());
            MADNESS_CHECK(p2.is_first());

            if (functor.has_f()) {
                for (const auto& r : rhs) result.push_back(madness::inner(functor.f,r,p1.get_tuple(),p2.get_tuple()));

            } else if (functor.has_f12()) {
                // functor is now a(1) b(2) f12
                // result(1) = \int a(1) f(1,2) b(2) rhs(2) d2
                World& world=rhs.front().world();
                auto premultiply= p1.is_first() ? functor.a : functor.b;
                auto postmultiply= p1.is_first() ? functor.b : functor.a;

                const int nbatch=30;
                for (int i=0; i<rhs.size(); i+=nbatch) {
                    std::vector<Function<T,LDIM>> tmp;
                    auto begin= rhs.begin()+i;
                    auto end= (i+nbatch)<rhs.size() ? rhs.begin()+i+nbatch : rhs.end();
                    std::copy(begin,end, std::back_inserter(tmp));

                    if (premultiply.is_initialized()) tmp=tmp*premultiply;
                    auto tmp1=apply(world,*(functor.f12),tmp);
                    if (postmultiply.is_initialized()) tmp1=tmp1*postmultiply;
                    for (auto& t : tmp1) result.push_back(t);

                }

            } else {
                MADNESS_EXCEPTION("confused functor in LowRankFunction",1);
            }
            return result;

        }

        /// inner product: result(1) = \int f(1,2) delta(2) d2

        /// @param[in]  functor the hidim function
        /// @param[in]  grid     grid points with delta functions
        /// @param[in]  p1      the variable in f(1,2) to be integrated over
        /// @param[in]  p2      the variable in rhs to be integrated over (usually all of them)
        static std::vector<Function<T,LDIM>> inner(LRFunctor& functor, const std::vector<Vector<double,LDIM>>& grid,
                                                   const particle<LDIM> p1, const particle<LDIM> p2) {
            std::vector<Function<T,LDIM>> result;
            MADNESS_CHECK(functor.has_f() xor functor.has_f12());
            MADNESS_CHECK(p1.is_first() xor p1.is_last());
            MADNESS_CHECK(p2.is_first());

            if (functor.has_f()) {
                MADNESS_EXCEPTION("no grid points with an explicit hi-dim function",1);

            } else if (functor.has_f12()) {
                MADNESS_CHECK(functor.f12->info.type==OT_SLATER);
                // functor is now a(1) b(2) f12
                // result(1) = \int a(1) f(1,2) b(2) delta(R-2) d2
                //           = a(1) f(1,R) b(R)
                World& world=functor.f12->get_world();
                auto premultiply= p1.is_first() ? functor.a : functor.b;
                auto postmultiply= p1.is_first() ? functor.b : functor.a;

                std::vector<Function<T,LDIM>> f1R= slater_functions_on_grid(world,grid);
                print("bla1");
                if (premultiply.is_initialized()) {
                    for (int i=0; i<grid.size(); ++i) f1R[i] = f1R[i]*premultiply(grid[i]);
                }
                print("bla2");
                if (postmultiply.is_initialized()) f1R=f1R*postmultiply;
                print("bla3");
                result=f1R;

            } else {
                MADNESS_EXCEPTION("confused functor in LowRankFunction",1);
            }
            return result;

        }


        static std::vector<Function<T,LDIM>> slater_functions_on_grid(World& world, const std::vector<Vector<double,LDIM>>& grid) {
            std::vector<Function<T,LDIM>> result;
            for (const auto& point : grid) {
                auto sl=[&point](const Vector<double,LDIM>& r) {
                    return exp(-sqrt(madness::inner(r-point,r-point)+1.e-8));
                };
                result.push_back(FunctionFactory<T,LDIM>(world).functor(sl));
            }
            return result;
        }


        World& world;
        double rank_revealing_tol=1.e-8;     // rrcd tol
        std::string orthomethod="symmetric";
        bool do_print=true;
        std::vector<Function<T,LDIM>> g,h;
        LRFunctor lrfunctor;
        particle<LDIM> p1=particle<LDIM>::particle1();
        particle<LDIM> p2=particle<LDIM>::particle2();

        LowRankFunction(World& world) : world(world) {}

        /// construct from the hi-dim function f
        LowRankFunction(const Function<T,NDIM>& f) : LowRankFunction(f.world()) {
            lrfunctor.f=f;
        }

        /// construct from the hi-dim function  f12*a(1)(b(2)
        LowRankFunction(const std::shared_ptr<SeparatedConvolution<T,LDIM>> f12, const Function<T,LDIM>& a,
                        const Function<T,LDIM>& b) : LowRankFunction(a.world()) {
            lrfunctor.a=a;
            lrfunctor.b=b;
            lrfunctor.f12=f12;
        }

        LowRankFunction(std::vector<Function<T,LDIM>> g, std::vector<Function<T,LDIM>> h)
                : world(g.front().world()), g(g), h(h) {}

        LowRankFunction(const LowRankFunction& a) : world(a.world), g(copy(world, a.g)), h(copy(world, a.h)) {} // Copy constructor necessary

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
            return LowRankFunction(g + b.g, h + b.h);
        }
        /// subtraction
        LowRankFunction operator-(const LowRankFunction& b) const {
            return LowRankFunction(g - b.g, h - b.h);
        }

        /// in-place addition
        LowRankFunction& operator+=(const LowRankFunction& b) {
            g+=b.g;
            h+=b.h;
            return *this;
        }

        /// in-place subtraction
        LowRankFunction& operator-=(const LowRankFunction& b) {
            g-=b.g;
            h-=b.h;
            return *this;
        }

        /// scale by a scalar
        template<typename Q>
        LowRankFunction operator*(const Q a) const {
            return LowRankFunction<TensorResultType<T,Q>,NDIM>(g * a, Q(h));
        }

        /// in-place scale by a scalar (no type conversion)
        LowRankFunction operator*(const T a) const {
            return LowRankFunction(g * a, h);
        }

        /// following Halko

        /// ||A - Q Q^T A||< epsilon
        /// Y = A Omega && Q = QR(Y)
        /// || f(1,2) - \sum_i g_i(1) h_i(2) || < epsilon
        /// Y_i(1) = \int f(1,2) Omega_i(2) d2 && g_i(1) = QR(Y_i(1)) && h_i(2) = \int g_i^*(1) f(1,2) d1
        void project(const LowRankFunctionParameters& params) {
            timer t1(world);
            t1.do_print=do_print;
            orthomethod=params.orthomethod();
            rank_revealing_tol=params.tol();


            // make grid
            std::vector<Vector<double,LDIM>> grid;
            if (params.gridtype()=="random") {
                double volume=4.0/3.0*constants::pi *std::pow(params.radius(),3.0);
                long rank = lround(volume/(0.67*params.volume_element() ));
                for (int i=0; i<rank; ++i) {
                    auto tmp=randomgaussian<LDIM>::gaussian_random_distribution(0,params.radius());
                    auto cell=FunctionDefaults<LDIM>::get_cell();
                    auto is_in_cell = [&cell](const Vector<double,LDIM>& r) {
                        for (int d=0; d<LDIM; ++d) if (r[d]<cell(d,0) or r[d]>cell(d,1)) return false;
                        return true;
                    };
                    if (not is_in_cell(tmp)) continue;
                    grid.push_back(tmp);
                }
                if (world.rank()==0 and do_print) {
                    double volume = 4.0/3.0*constants::pi * std::pow(params.radius(),3.0);
                    print("volume element in random grid",volume/(0.67*rank));
                }


            } else if (params.gridtype()=="spherical") {
                cgrid<LDIM> cg(params.volume_element(),params.radius(),params.hard_shell());
                grid=cg.get_coordinates();
            } else {
                MADNESS_EXCEPTION("unknown grid type in project",1);
            }


            auto Y=Yformer(grid,params.rhsfunctiontype());
            t1.tag("Yforming");

            auto ovlp=matrix_inner(world,Y,Y);  // error in symmetric matrix_inner, use non-symmetric form here!
            t1.tag("compute ovlp");
            g=truncate(orthonormalize_rrcd(Y,ovlp,rank_revealing_tol));
            t1.tag("rrcd/truncate/thresh");
            auto sz=get_size(world,g);
            print("gsize",sz);
            check_orthonormality(g);

            if (world.rank()==0 and do_print) {
                print("Y.size()",Y.size());
                print("g.size()",g.size());
            }

            h=truncate(inner(lrfunctor,g,p1,p1));
            t1.tag("Y backprojection with truncation");

        }

        /// apply a rhs (delta or exponential) on grid points to the hi-dim function and form Y = A_ij w_j (in Halko's language)
        std::vector<Function<T,LDIM>> Yformer(const std::vector<Vector<double,LDIM>>& grid, const std::string rhsfunctiontype,
                                              const double exponent=30.0) {

            std::vector<Function<double,LDIM>> Y;
            if (rhsfunctiontype=="delta") {
                Y=inner(lrfunctor,grid,p2,p1);

            } else if (rhsfunctiontype=="exponential") {
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
                Y=inner(lrfunctor,omega,p2,p1);
            } else {
                MADNESS_EXCEPTION("confused rhsfunctiontype",1);
            }
            auto norms=norm2s(world,Y);
            std::vector<Function<double,LDIM>> Ynormalized;

            for (int i=0; i<Y.size(); ++i) if (norms[i]>rank_revealing_tol) Ynormalized.push_back(Y[i]);
            normalize(world,Ynormalized);
            return Ynormalized;
        }

        long rank() const {return g.size();}

        Function<T,NDIM> reconstruct() const {
            auto fapprox=hartree_product(g[0],h[0]);
            for (int i=1; i<g.size(); ++i) fapprox+=hartree_product(g[i],h[i]);
            return fapprox;
        }


        std::vector<Function<T,LDIM>> orthonormalize(const std::vector<Function<T,LDIM>>& g) const {

            double tol=rank_revealing_tol;
            std::vector<Function<T,LDIM>> g2;
            auto ovlp=matrix_inner(world,g,g);
            if (orthomethod=="symmetric") {
                print("orthonormalizing with method/tol",orthomethod,tol);
                g2=orthonormalize_symmetric(g,ovlp);
            } else if (orthomethod=="canonical") {
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


        /// optimize using Cholesky decomposition

        /// @param[in]  nopt       number of iterations (wrt to Alg. 4.3 in Halko)
        void optimize(const long nopt=1) {
            timer t(world);
            t.do_print=do_print;
            for (int i=0; i<nopt; ++i) {
                // orthonormalize h
                h=truncate(orthonormalize(h));
                t.tag("ortho1");
                g=truncate(inner(lrfunctor,h,p2,p1));
                t.tag("inner1");
                g=truncate(orthonormalize(g));
                t.tag("ortho2");
                h=truncate(inner(lrfunctor,g,p1,p1));
                t.tag("inner2");
            }
        }

        /// after external operations g might not be orthonormal and/or optimal

        /// orthonormalization similar to Bischoff, Harrison, Valeev, JCP 137 104103 (2012), Sec II C 3
        /// f   =\sum_i g_i h_i
        ///     = g X- (X+)^T (Y+)^T Y- h
        ///     = g X-  U S V^T  Y- h
        ///     = g (X- U) (S V^T Y-) h
        /// requires 2 matrix_inner and 2 transforms. g and h are optimal, but contain all cusps etc..
        /// @param[in]  thresh        SVD threshold
        void reorthonormalize(double thresh=-1.0) {
            if (thresh<0.0) thresh=rank_revealing_tol;
            Tensor<T> ovlp_g = matrix_inner(world, g, g);
            Tensor<T> ovlp_h = matrix_inner(world, h, h);
            auto [eval_g, evec_g] = syev(ovlp_g);
            auto [eval_h, evec_h] = syev(ovlp_h);
            Tensor<T> Xplus=copy(evec_g);
            Tensor<T> Xminus=copy(evec_g);
            Tensor<T> Yplus=copy(evec_h);
            Tensor<T> Yminus=copy(evec_h);
//            print("eval_g",eval_g);
//            print("eval_h",eval_h);
            for (int i=0; i<eval_g.size(); ++i) Xplus(_,i)*=std::pow(eval_g(i),0.5);
            for (int i=0; i<eval_g.size(); ++i) Xminus(_,i)*=std::pow(eval_g(i),-0.5);
            for (int i=0; i<eval_h.size(); ++i) Yplus(_,i)*=std::pow(eval_h(i),0.5);
            for (int i=0; i<eval_h.size(); ++i) Yminus(_,i)*=std::pow(eval_h(i),-0.5);

            // test
            //{
            //    auto gortho = transform(world, h, Yminus);
            //    auto one=matrix_inner(world,gortho,gortho);
            //    check_orthonormality(one);
            //}
            //{
            //    auto gortho = transform(world, g, Xminus);
            //    auto one=matrix_inner(world,gortho,gortho);
            //    check_orthonormality(one);
            //}
            //{
            //    auto one=madness::inner(Xminus,Xplus,1,1);
            //    check_orthonormality(one);
            //    one=madness::inner(Yminus,Yplus,1,1);
            //    check_orthonormality(one);
            //}

            Tensor<T> M=madness::inner(Xplus,Yplus,0,0);    // (X+)^T Y+
            auto [U,s,VT]=svd(M);
//            print("s",s);

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
            print("accumulated s",s_accumulated,thresh,s.size(),i);
            for (int j=0; j<s.size(); ++j) VT(j,_)*=s[j];
            // test
//            auto mm=madness::inner(U,VT,1,0);
//            double none=(mm-M).normf();
//            print("U S V^T - M",none);


            Tensor<T> XX=madness::inner(Xminus,U,1,0);
            Tensor<T> YY=madness::inner(Yminus,VT,1,1);

            g=truncate(transform(world,g,XX));
            h=truncate(transform(world,h,YY));
//            {
//                auto one=matrix_inner(world,g,g);
//                check_orthonormality(one);
//            }
//            {
//                auto one=matrix_inner(world,h,h);
//                check_orthonormality(one);
//            }



        }


        double check_orthonormality(const std::vector<Function<T,LDIM>>& v) const {
            Tensor<T> ovlp=matrix_inner(world,v,v);
            return check_orthonormality(ovlp);
        }

        double check_orthonormality(const Tensor<T>& ovlp) const {
            timer t(world);
            Tensor<T> ovlp2=ovlp;
            for (int i=0; i<ovlp2.dim(0); ++i) ovlp2(i,i)-=1.0;
            print("absmax",ovlp2.absmax());
            print("l2",ovlp2.normf()/ovlp2.size());
            t.tag("check_orthonoramality");
            return ovlp.absmax();
        }

        double explicit_error() const {
            auto fapprox=reconstruct();
            return (lrfunctor.f-fapprox).norm2();
        }

        double randomized_error() const {
            return 1.e9;
        }

        double error() const {
            if (LDIM<3) return explicit_error();
            else return randomized_error();
        }

        /// compute the l2 error |functor - \sum_i g_ih_i|_2

        /// \int (f(1,2) - gh(1,2))^2 = \int f(1,2)^2 - 2\int f(1,2) gh(1,2) + \int gh(1,2)^2
        double l2error() const {
            MADNESS_CHECK(lrfunctor.has_f12());

            timer t(world);
            t.do_print=do_print;
            const Function<T,LDIM> one=FunctionFactory<T,LDIM>(world).f([](const Vector<double,LDIM>& r){return 1.0;});
            const Function<T,LDIM> pre=(lrfunctor.a.is_initialized()) ? lrfunctor.a : one;
            const Function<T,LDIM> post=(lrfunctor.b.is_initialized()) ? lrfunctor.b : one;
            const SeparatedConvolution<T,LDIM>& f12=*lrfunctor.f12;
            const SeparatedConvolution<T,LDIM> f12sq= SeparatedConvolution<T,LDIM>::combine(f12,f12);

            // \int f(1,2)^2 d1d2 = \int f(1,2)^2 pre(1)^2 post(2)^2 d1 d2
            double term1 =madness::inner(post*post,f12sq(pre*pre));
            t.tag("computing term1");

            // \int f(1,2) pre(1) post(2) \sum_i g(1) h(2) d1d2
            double term2=madness::inner(pre*g,f12(post*h));
            t.tag("computing term2");

            // g functions are orthonormal
            // \int gh(1,2)^2 d1d2 = \int \sum_{ij} g_i(1) g_j(1) h_i(2) h_j(2) d1d2
            //   = \sum_{ij} \int g_i(1) g_j(1) d1 \int h_i(2) h_j(2) d2
            //   = \sum_{ij} delta_{ij} \int h_i(2) h_j(2) d2
            //   = \sum_{i} \int h_i(2) h_i(2) d2
            double zero=check_orthonormality(g);
            if (zero>1.e-10) print("g is not orthonormal",zero);
            double term3=madness::inner(h,h);
            t.tag("computing term3");

            double error=sqrt(term1-2.0*term2+term3)/sqrt(term1);
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



    template<typename T, std::size_t NDIM>
    LowRankFunction<T,NDIM> inner(const Function<T,NDIM>& lhs, const LowRankFunction<T,NDIM>& rhs, const std::tuple<int> v1, const std::tuple<int> v2) {
        World& world=rhs.world;
        // int lhs(1,2) rhs(2,3) d2 = \sum \int lhs(1,2) g_i(2) h_i(3) d2
        //                      = \sum \int lhs(1,2) g_i(2) d2 h_i(3)
        LowRankFunction<T, NDIM + NDIM - 2> result(world);
        result.h=rhs.h;
        decltype(rhs.g) g;
        for (int i=0; i<rhs.rank(); ++i) {
            g.push_back(inner(lhs,rhs.g[i],{v1},{0}));
        }
        result.g=g;
        return result;
    }

    template<typename T, std::size_t NDIM>
    LowRankFunction<T,NDIM> inner(const LowRankFunction<T,NDIM>& f, const Function<T,NDIM>& g, const std::tuple<int> v1, const std::tuple<int> v2) {
        World& world=f.world;
        // int f(1,2) k(2,3) d2 = \sum \int g_i(1) h_i(2) k(2,3) d2
        //                      = \sum g_i(1) \int h_i(2) k(2,3) d2
        LowRankFunction<T, NDIM + NDIM - 2> result(world);
        result.g=f.g;
        decltype(f.h) h;
        for (int i=0; i<f.rank(); ++i) {
            h.push_back(inner(f.h[i],g,{0},{v2}));
        }
        result.h=h;
        return result;
    }



} // namespace madness

#endif //MADNESS_LOWRANKFUNCTION_H
