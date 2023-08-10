//
// Created by Florian Bischoff on 8/10/23.
//

#ifndef MADNESS_LOWRANKFUNCTION_H
#define MADNESS_LOWRANKFUNCTION_H


#include<madness/mra/mra.h>
#include<madness/chem/electronic_correlation_factor.h>

namespace madness {



    template<std::size_t NDIM>
    struct randomgaussian {
        Vector<double,NDIM> random_origin;
        double exponent;
        double radius=2;
        randomgaussian(double exponent, double radius) : exponent(exponent), radius(radius) {
            Vector<double,NDIM> ran; // [0,1]
            RandomVector(NDIM,ran.data());
            random_origin=2.0*radius*ran-Vector<double,NDIM>(radius);
    //        print("origin at ",random_origin, ", exponent",exponent);
        }
        double operator()(const Vector<double,NDIM>& r) const {
    //        return exp(-exponent*inner(r-random_origin,r-random_origin));
            return exp(-exponent*(r-random_origin).normf());
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
        };
        template<std::size_t PDIM>
        struct particle {
            std::array<int,PDIM> dims;
            particle(const int p) : particle(std::vector<int>(1,p)) {}
            particle(const std::vector<int> p) {
                for (int i=0; i<PDIM; ++i) dims[i]=p[i];
            }

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

        static std::vector<Function<T,LDIM>> inner(LRFunctor& functor, const std::vector<Function<T,LDIM>>& rhs,
                                                   const particle<LDIM> p1, const particle<LDIM> p2) {
            std::vector<Function<T,LDIM>> result;
            if (functor.has_f()) {
                for (const auto& r : rhs) result.push_back(madness::inner(functor.f,r,p1.get_tuple(),p2.get_tuple()));
                return result;
            }

        }

        World& world;
        std::vector<Function<T,LDIM>> g,h;
        LRFunctor lrfunctor;

        LowRank(World& world) : world(world) {};

        /// construct from the hi-dim function f
        LowRank(const Function<T,NDIM>& f) : world(f.world()) {
            lrfunctor.f=f;
        }

        /// construct from the hi-dim function  f12*a(1)(b(2)
        LowRank(const std::shared_ptr<SeparatedConvolution<T,LDIM>> f12, const Function<T,LDIM>& a,
                const Function<T,LDIM>& b) : world(a.world) {
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

        /// following Halko

        /// ||A - Q Q^T A||< epsilon
        /// Y = A Omega && Q = QR(Y)
        /// || f(1,2) - \sum_i g_i(1) h_i(2) || < epsilon
        /// Y_i(1) = \int f(1,2) Omega_i(2) d2 && g_i(1) = QR(Y_i(1)) && h_i(2) = \int g_i^*(1) f(1,2) d1
        void project(const long rank, const double radius=3.0) {
            timer t1(world);
            std::vector<Function<double,LDIM>> omega2(rank);
            for (long i=0; i<rank; ++i) {
                omega2[i]=FunctionFactory<double,LDIM>(world).functor(randomgaussian<LDIM>(RandomValue<double>()+3,radius));
            }
            t1.tag("projection 1D functions");

            std::vector<Function<double,LDIM>> Y(rank);

            if constexpr (LDIM==1)  Y=inner(lrfunctor,omega2,{1},{0});
            if constexpr (LDIM==2)  Y=inner(lrfunctor,omega2,{2,3},{0,1});
            if constexpr (LDIM==3)  Y=inner(lrfunctor,omega2,{3,4,5},{0,1,2});
            t1.tag("Yforming");
            print("Y.size()",Y.size());

            g=orthonormalize_canonical(Y,1.e-12);
            print("g.size()",g.size());
            t1.tag("Y orthonormalizing");
            h.resize(g.size());
            if constexpr (LDIM==1) h=inner(lrfunctor,g,{0},{0});
            if constexpr (LDIM==2) h=inner(lrfunctor,g,{0,1},{0,1});
            if constexpr (LDIM==3) h=inner(lrfunctor,g,{0,1,2},{0,1,2});

            t1.tag("Y backprojection");

        }

        long rank() const {return g.size();}

        Function<T,NDIM> reconstruct() const {
            auto fapprox=hartree_product(g[0],h[0]);
            for (int i=1; i<g.size(); ++i) fapprox+=hartree_product(g[i],h[i]);
            return fapprox;
        }

        /// @return     the singular values s
        Tensor<double> orthonormalize(const bool s_in_h) {
            timer t(world);
            /**
             *  |g >< h| = |g_ortho><g_ortho | g> < h | h_ortho ><h_ortho |
             *           = |g_ortho> gg hh <h_ortho |
             *           = |g_ortho> U s VT <h_ortho |
             */
            std::vector<Function<T,LDIM>> g_ortho=orthonormalize_canonical(g,1.e-8);
            std::vector<Function<T,LDIM>> h_ortho=orthonormalize_canonical(h,1.e-8);
            auto gg=matrix_inner(world,g_ortho,g);
            auto hh=matrix_inner(world,h,h_ortho);
            auto ovlp=madness::inner(gg,hh);
            Tensor<T> U,VT;
            Tensor<double> s;
            svd(ovlp,U,s,VT);
            auto V=transpose(VT);

            g=transform(world,g_ortho,U);
            h=transform(world,h_ortho,V);


            // include singular values into h
            if (s_in_h) for (int i=0; i<h.size(); ++i) h[i]*=s[i];
            t.tag("orthonormalization");
            return s;

        }

        /// assumes g and h to be orthonormal -> simple projection
        void optimize(const long nopt=2) {

            timer t(world);
            auto s=orthonormalize(true);    // includes singular values in h -> not normalized any more
            for (int iopt=0; iopt<nopt; ++iopt) {
                std::vector<Function<double,LDIM>> gtmp(g.size());
                /// remove singular values from h again -> normalized
                if constexpr (LDIM == 1) gtmp= inner(lrfunctor, h, {1}, {0});
                if constexpr (LDIM == 2) gtmp= inner(lrfunctor, h, {2, 3}, {0, 1});
                if constexpr (LDIM == 3) gtmp= inner(lrfunctor, h, {3, 4, 5}, {0, 1, 2});
                std::vector<double> sinv(s.size());
                for (int i=0; i<s.size(); ++i) sinv[i]=1.0/s[i];
                scale(world,gtmp,sinv);
                g=orthonormalize_canonical(gtmp,1.e-12);
                std::vector<Function<double,LDIM>> htmp(g.size());
                if constexpr (LDIM==1) htmp=inner(lrfunctor,g,{0},{0});
                if constexpr (LDIM==2) htmp=inner(lrfunctor,g,{0,1},{0,1});
                if constexpr (LDIM==3) htmp=inner(lrfunctor,g,{0,1,2},{0,1,2});
                h=htmp;

                if (g.size()>1) s=orthonormalize(true);

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
