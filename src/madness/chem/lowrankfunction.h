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


    template<std::size_t NDIM>
    struct cartesian_grid {
        Vector<double,NDIM> lovec,hivec;
        std::vector<long> stride;
        long index=0;
        long n_per_dim;
        long total_n;
        Vector<double,NDIM> increment;

        cartesian_grid(const long n_per_dim, const double lo, const double hi)
                : n_per_dim(n_per_dim) {
            lovec.fill(lo);
            hivec.fill(hi);
            increment=(hivec-lovec)*(1.0/double(n_per_dim-1));
            stride=std::vector<long>(NDIM,1l);
            total_n=std::pow(n_per_dim,NDIM);
            for (long i=NDIM-2; i>=0; --i) stride[i]=n_per_dim*stride[i+1];
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
        double rank_revealing_tol=1.e-12;     // rrcd tol
        bool do_print=true;
        bool stable_power_iteration=true;
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

        void project(const double volume_per_point, const double radius, const std::string gridtype, std::string rhsfunctiontype,
                     double tol1) {
            rank_revealing_tol=tol1;
            long rank=0;
            if (gridtype=="random") {
                // number of points within radius = variance: 0.67 * #total points = 0.67*rank
                // volume of sphere 4/3 pi *r^3
                // volume per point=volume/#points in radius = volume / (0.67 * rank)
                // rank= volume/(0.67/volume_per_point)
                double volume=4.0/3.0*constants::pi *std::pow(radius,3.0);
                rank = lround(volume/(0.67*volume_per_point ));
            }
            project(rank,radius,gridtype,rhsfunctiontype);
        }

        /// following Halko

        /// ||A - Q Q^T A||< epsilon
        /// Y = A Omega && Q = QR(Y)
        /// || f(1,2) - \sum_i g_i(1) h_i(2) || < epsilon
        /// Y_i(1) = \int f(1,2) Omega_i(2) d2 && g_i(1) = QR(Y_i(1)) && h_i(2) = \int g_i^*(1) f(1,2) d1
        void project(const long rank, const double radius, const std::string gridtype, std::string rhsfunctiontype) {
            timer t1(world);
            t1.do_print=do_print;

            // make grid
            std::vector<Vector<double,LDIM>> grid;
            if (gridtype=="random") {
                for (int i=0; i<rank; ++i) {
                    auto tmp=randomgaussian<LDIM>::gaussian_random_distribution(0,radius);
                    auto cell=FunctionDefaults<LDIM>::get_cell();
                    auto is_in_cell = [&cell](const Vector<double,LDIM>& r) {
                        for (int d=0; d<LDIM; ++d) if (r[d]<cell(d,0) or r[d]>cell(d,1)) return false;
                        return true;
                    };
                    if (not is_in_cell(tmp)) continue;
                    grid.push_back(tmp);
                }
                if (world.rank()==0 and do_print) {
                    double volume = 4.0/3.0*constants::pi * std::pow(radius,3.0);
                    print("volume element in random grid",volume/(0.67*rank));
                }


            } else if (gridtype=="cartesian") {
                long nperdim=std::lround(std::pow(double(rank),1.0/3.0));
                cartesian_grid<LDIM> cg(nperdim,-radius,radius);
                for (cg.index=0; cg(); ++cg) grid.push_back(cg.get_coordinates());
                if (world.rank()==0 and do_print) print("volume element in cartesian grid",cg.volume_per_gridpoint());
            } else {
                MADNESS_EXCEPTION("unknown grid type in project",1);
            }

            if (world.rank()==0 and do_print) {
                print("grid is",gridtype,"with radius",radius,"and",grid.size(),"gridpoints");
                print("rhs functions are",rhsfunctiontype);
            }

            auto Y=Yformer(grid,rhsfunctiontype);
            t1.tag("Yforming");

            std::ostringstream oss;
            oss << std::scientific << std::setprecision(1) << rank_revealing_tol;
            std::string scientificString = oss.str();
            double err=1.0;
            double tight_thresh=FunctionDefaults<3>::get_thresh()*0.1;
            g=Y;
            while (err>1.e-10) {
                g=truncate(orthonormalize_rrcd(g,rank_revealing_tol),std::min(1.e-3,tight_thresh*100.0));
                err=check_orthonormality(g);
                print("error of non-orthonormality",err);
                t1.tag("Y orthonormalizing with rank_revealing_tol loose_thresh "+scientificString);
            }
            g=truncate(g);

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
                for (const auto& point : grid) {
                    omega.push_back(FunctionFactory<double,LDIM>(world)
                            .functor([&point,&exponent](const Vector<double,LDIM>& r)
                                     {
                                         auto r_rel=r-point;
                                         return exp(-exponent*madness::inner(r_rel,r_rel));
                                     }));
                }
                Y=inner(lrfunctor,omega,p2,p1);
            } else {
                MADNESS_EXCEPTION("confused rhsfunctiontype",1);
            }
            return Y;
        }

        long rank() const {return g.size();}

        Function<T,NDIM> reconstruct() const {
            auto fapprox=hartree_product(g[0],h[0]);
            for (int i=1; i<g.size(); ++i) fapprox+=hartree_product(g[i],h[i]);
            return fapprox;
        }

        /// optimize using Cholesky decomposition

        /// if stable_power_iteration is true, orthonormalize in between applications of the kernel (Alg. 4.4 in Halko)
        /// @param[in]  nopt       number of iterations (wrt to Alg. 4.3 in Halko)
        void optimize(const long nopt=1) {
            timer t(world);
            t.do_print=do_print;
            double tight_thresh=FunctionDefaults<3>::get_thresh()*0.1;
            for (int i=0; i<nopt; ++i) {
                // orthonormalize h
                if (stable_power_iteration) h=truncate(orthonormalize_rrcd(h,rank_revealing_tol),tight_thresh);
                t.tag("ortho1 with rrcd/truncate/tight");
                g=truncate(inner(lrfunctor,h,p2,p1),tight_thresh);
                t.tag("inner1/truncate/tight");
                g=truncate(orthonormalize_rrcd(g,rank_revealing_tol),tight_thresh);
                t.tag("ortho2/truncate/tight");
                h=truncate(inner(lrfunctor,g,p1,p1),tight_thresh);
                t.tag("inner2/truncate/tight");
            }
            t.tag("optimize_fast");
        }

        double check_orthonormality(const std::vector<Function<T,LDIM>>& v) const {
            Tensor<T> ovlp=matrix_inner(world,v,v);
            for (int i=0; i<ovlp.dim(0); ++i) ovlp(i,i)-=1.0;
            return ovlp.normf()/ovlp.size();
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
