//
// Created by Florian Bischoff on 8/10/23.
//

#ifndef MADNESS_LOWRANKFUNCTION_H
#define MADNESS_LOWRANKFUNCTION_H


#include<madness/mra/mra.h>
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


    template<typename T, std::size_t NDIM, std::size_t LDIM=NDIM/2>
    class LowRank {
    public:

        struct LRFunctor {
            LRFunctor() = default;

            Function<T, NDIM> f;
            std::shared_ptr<SeparatedConvolution<T,LDIM>> f12;
            Function<T,LDIM> a,b;
            bool has_f() const {
                return f.is_initialized();
            }
            bool has_f12() const {
                return (f12.get());
            }
            T operator()(const Vector<double,NDIM>& r) const {
                Vector<double,LDIM> first, second;
                for (int i=0; i<LDIM; ++i) {
                    first[i]=r[i];
                    second[i]=r[i+LDIM];
                }
                MADNESS_CHECK(has_f12());
                double result=a(first)*b(second)*exp(-(first-second).normf());
                return result;
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
        std::vector<Function<T,LDIM>> g,h;
        LRFunctor lrfunctor;
        particle<LDIM> p1=particle<LDIM>::particle1();
        particle<LDIM> p2=particle<LDIM>::particle2();

        LowRank(World& world) : world(world) {}

        /// construct from the hi-dim function f
        LowRank(const Function<T,NDIM>& f) : LowRank(f.world()) {
            lrfunctor.f=f;
        }

        /// construct from the hi-dim function  f12*a(1)(b(2)
        LowRank(const std::shared_ptr<SeparatedConvolution<T,LDIM>> f12, const Function<T,LDIM>& a,
                const Function<T,LDIM>& b) : LowRank(a.world()) {
            lrfunctor.a=a;
            lrfunctor.b=b;
            lrfunctor.f12=f12;
        }

        LowRank(std::vector<Function<T,LDIM>> g, std::vector<Function<T,LDIM>> h)
                : world(g.front().world()), g(g), h(h) {}

        LowRank(const LowRank& a) : world(a.world), g(copy(world,a.g)), h(copy(world,a.h)) {} // Copy constructor necessary

        LowRank& operator=(const LowRank& f) { // Assignment required for storage in vector
            LowRank ff(f);
            std::swap(ff.g,g);
            std::swap(ff.h,h);
            return *this;
        }

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

        LowRank operator-(const LowRank& b) const { // Operator- necessary
            return LowRank(g-b.g,h-b.h);
        }

        LowRank& operator+=(const LowRank& b) { // Operator+= necessary
            g+=b.g;
            h+=b.h;
            return *this;
        }

        LowRank operator*(double a) const { // Scale by a constant necessary
            return LowRank(g*a,h);
        }

        LowRank multiply(const std::vector<Function<T,LDIM>>& vec, const long particle) {
            auto result=*this;  // deep copy
            if (particle==0) result.g=g*vec;
            if (particle==1) result.h=h*vec;
            return *this;
        }

        void project(const double volume_per_point, const double radius, const std::string gridtype, std::string rhsfunctiontype) {
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
                double volume = 4.0/3.0*constants::pi * std::pow(radius,3.0);
                print("volume element in random grid",volume/(0.67*rank));


            } else if (gridtype=="cartesian") {
                long nperdim=std::lround(std::pow(double(rank),1.0/3.0));
                cartesian_grid<LDIM> cg(nperdim,-radius,radius);
                for (cg.index=0; cg(); ++cg) grid.push_back(cg.get_coordinates());
                print("volume element in cartesian grid",cg.volume_per_gridpoint());
            } else {
                MADNESS_EXCEPTION("unknown grid type in project",1);
            }

            print("grid is",gridtype,"with radius",radius,"and",grid.size(),"gridpoints");
            print("rhs functions are",rhsfunctiontype);

            auto Y=Yformer(grid,rhsfunctiontype);
            t1.tag("Yforming");

            double tol=1.e-12;
            g=orthonormalize_rrcd(Y,tol);
            t1.tag("Y orthonormalizing with tol"+std::to_string(tol));

            print("Y.size()",Y.size());
            print("g.size()",g.size());

            h=inner(lrfunctor,g,p1,p1);
            t1.tag("Y backprojection");

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

        Tensor<double> orthonormalize(const bool s_in_h) {
            return orthonormalize_svd(s_in_h);
        }


        /// @return     the singular values s
        Tensor<double> orthonormalize_svd(const bool s_in_h) {
            timer t(world);
            /**
             *  |g >< h| = |g_ortho><g_ortho | g> < h | h_ortho ><h_ortho |
             *           = |g_ortho> gg hh <h_ortho |
             *           = |g_ortho> U s VT <h_ortho |
             */
            std::vector<Function<T,LDIM>> g_ortho=orthonormalize_rrcd(g,1.e-8);
            t.tag("ortho1");
            std::vector<Function<T,LDIM>> h_ortho=orthonormalize_rrcd(h,1.e-8);
            t.tag("ortho2");
            auto gg=matrix_inner(world,g_ortho,g);
            auto hh=matrix_inner(world,h,h_ortho);
            auto ovlp=madness::inner(gg,hh);
            Tensor<T> U,VT;
            Tensor<double> s;
            svd(ovlp,U,s,VT);
            auto V=transpose(VT);


            t.tag("svd");
            g=transform(world,g_ortho,U);
            t.tag("transform1");
            h=transform(world,h_ortho,V);
            t.tag("transform2");


            // include singular values into h
            if (s_in_h) for (int i=0; i<h.size(); ++i) h[i]*=s[i];
            t.end("orthonormalization");
            return s;

        }
        void orthonormalize_cd() {
            /**
             *  |g >s< h| = |g_ortho><g_ortho | g> s < h |
             *           = |g_ortho> gg s < h |
             *           = |g_ortho> <h_non_ortho |
             */
            World& world=g.front().world();

            auto g_ortho= orthonormalize_rrcd(g,1.e-10);
//        auto g_ortho= orthonormalize_canonical(g,1.e-8);      // similar to SVD
//        auto g_ortho= orthonormalize_symmetric(g);

            auto ovlp=matrix_inner(world,g_ortho,g_ortho);
            for (int i=0; i<ovlp.dim(0); ++i) ovlp(i,i)-=1.0;
            MADNESS_CHECK(fabs(ovlp.normf()/ovlp.size()<1.e-10));

            Tensor<double> gg=matrix_inner(world,g_ortho,g);
            /// Transforms a vector of functions according to new[i] = sum[j] old[j]*c[j,i]
//            h=truncate(transform(world,h,transpose(gg)));
//            g=truncate(g_ortho);
            h=(transform(world,h,transpose(gg)));
            g=(g_ortho);

        }

        void optimize(const long nopt=2) {
            optimize_cd(nopt);
        }

        /// optimize using Cholesky decomposition
        void optimize_cd(const long nopt) {

            timer t(world);
            for (int iopt=0; iopt<nopt; ++iopt) {
//                this->orthonormalize_cd();      // g orthonormal, h not
                t.tag("ortho1");
                auto gtmp=inner(lrfunctor,h,p2,p1);
                t.tag("inner1");
                g=orthonormalize_rrcd(gtmp,1.e-12);
                t.tag("ortho2");
                h=inner(lrfunctor,g,p1,p1);
                t.tag("inner2");

                if (iopt%2==1) {
                    double err=error();
                    print("optimization iteration",iopt,", error in f_approx_opt",err);
                }
            }
            t.tag("optimize");
        }

        /// assumes g and h to be orthonormal -> simple projection
        void optimize_svd(const long nopt=2) {

            timer t(world);
            for (int iopt=0; iopt<nopt; ++iopt) {
                auto s=orthonormalize(true);    // includes singular values in h -> not normalized any more
                std::vector<Function<double,LDIM>> gtmp(g.size());
                /// remove singular values from h again -> normalized
                if constexpr (LDIM == 1) gtmp= inner(lrfunctor, h, {1}, {0});
                if constexpr (LDIM == 2) gtmp= inner(lrfunctor, h, {2, 3}, {0, 1});
                if constexpr (LDIM == 3) gtmp= inner(lrfunctor, h, {3, 4, 5}, {0, 1, 2});
                std::vector<double> sinv(s.size());
                for (int i=0; i<s.size(); ++i) sinv[i]=1.0/s[i];
                scale(world,gtmp,sinv);
//                g=orthonormalize_canonical(gtmp,1.e-12);
                g=orthonormalize_rrcd(gtmp,1.e-12);
                std::vector<Function<double,LDIM>> htmp(g.size());
                if constexpr (LDIM==1) htmp=inner(lrfunctor,g,{0},{0});
                if constexpr (LDIM==2) htmp=inner(lrfunctor,g,{0,1},{0,1});
                if constexpr (LDIM==3) htmp=inner(lrfunctor,g,{0,1,2},{0,1,2});
                h=htmp;

//                if (g.size()>1) s=orthonormalize(true);

                if (iopt%2==1) {
                    double err=error();
                    print("optimization iteration",iopt,", error in f_approx_opt",err);
                }
            }
            t.tag("optimize");
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

    //    double get() const {return x;}
    };

    // This interface is necessary to compute inner products
    template<typename T, std::size_t NDIM>
    double inner(const LowRank<T,NDIM>& a, const LowRank<T,NDIM>& b) {
        World& world=a.world;
        return (matrix_inner(world,a.g,b.g).emul(matrix_inner(world,a.h,b.h))).sum();
    }



    template<typename T, std::size_t NDIM>
    LowRank<T,NDIM> inner(const Function<T,NDIM>& lhs, const LowRank<T,NDIM>& rhs, const std::tuple<int> v1, const std::tuple<int> v2) {
        World& world=rhs.world;
        // int lhs(1,2) rhs(2,3) d2 = \sum \int lhs(1,2) g_i(2) h_i(3) d2
        //                      = \sum \int lhs(1,2) g_i(2) d2 h_i(3)
        LowRank<T,NDIM+NDIM-2> result(world);
        result.h=rhs.h;
        decltype(rhs.g) g;
        for (int i=0; i<rhs.rank(); ++i) {
            g.push_back(inner(lhs,rhs.g[i],{v1},{0}));
        }
        result.g=g;
        return result;
    }

    template<typename T, std::size_t NDIM>
    LowRank<T,NDIM> inner(const LowRank<T,NDIM>& f, const Function<T,NDIM>& g, const std::tuple<int> v1, const std::tuple<int> v2) {
        World& world=f.world;
        // int f(1,2) k(2,3) d2 = \sum \int g_i(1) h_i(2) k(2,3) d2
        //                      = \sum g_i(1) \int h_i(2) k(2,3) d2
        LowRank<T,NDIM+NDIM-2> result(world);
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
