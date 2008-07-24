#ifndef MAD_OPERATOR_H
#define MAD_OPERATOR_H

#include <mra/mra.h>
#include <limits.h>
#include <mra/adquad.h>
#include <tensor/mtxmq.h>
#include <tensor/aligned.h>
#include <linalg/tensor_lapack.h>

namespace madness {

    void aligned_add(long n, double* RESTRICT a, const double* RESTRICT b);
    void aligned_sub(long n, double* RESTRICT a, const double* RESTRICT b);
    void aligned_add(long n, double_complex* RESTRICT a, const double_complex* RESTRICT b);
    void aligned_sub(long n, double_complex* RESTRICT a, const double_complex* RESTRICT b);

    /// a(n,m) --> b(m,n) ... optimized for smallish matrices
    template <typename T>
    static void fast_transpose(long n, long m, const T* RESTRICT a, T* RESTRICT b) {
        // n will always be k or 2k (k=wavelet order) and m will be anywhere
        // from 2^(NDIM-1) to (2k)^(NDIM-1).

//         for (long i=0; i<n; i++)
//             for (long j=0; j<m; j++) 
//                 b[j*n+i] = a[i*m+j];

        if (n==1 || m==1) {
            long nm=n*m;
            for (long i=0; i<nm; i++) b[i] = a[i];
            return;
        }

        long n4 = (n>>2)<<2;
        long m4 = m<<2;
        const T* RESTRICT a0 = a;
        for (long i=0; i<n4; i+=4, a0+=m4) {
            const T* RESTRICT a1 = a0+m;
            const T* RESTRICT a2 = a1+m;
            const T* RESTRICT a3 = a2+m;
            T* RESTRICT bi = b+i;
            for (long j=0; j<m; j++, bi+=n) {
                T tmp0 = a0[j];
                T tmp1 = a1[j];
                T tmp2 = a2[j];
                T tmp3 = a3[j];

                bi[0] = tmp0;
                bi[1] = tmp1;
                bi[2] = tmp2;
                bi[3] = tmp3;
            }
        }

        for (long i=n4; i<n; i++)
            for (long j=0; j<m; j++) 
                b[j*n+i] = a[i*m+j];

    }

    /// a(i,j) --> b(i,j) for i=0..n-1 and j=0..r-1 noting dimensions are a(n,m) and b(n,r).

    /// returns b
    template <typename T>
    static T* shrink(long n, long m, long r, const T* a, T* RESTRICT b) {
        T* result = b;
        for (long i=0; i<n; i++, a+=m, b+=r) {
            for (long j=0; j<r; j++) {
                b[j] = a[j];
            }
        }
        return result;
    }

    extern void bsh_fit(double mu, double lo, double hi, double eps, 
                        Tensor<double> *pcoeff, Tensor<double> *pexpnt, bool prnt=false);

//     extern void bsh_fit_mod(double mu, double lo, double hi, double eps, 
//                             Tensor<double> *pcoeff, Tensor<double> *pexpnt, bool prnt=false);

    /// Simplified interface around hash_map to cache stuff

    /// Since insertions into STL containers have the nasty habit of
    /// invalidating iterators we actually store shared pointers
    template <typename Q>
    class SimpleCache {
    private:
        typedef HASH_MAP_NAMESPACE::hash_map< unsigned long, SharedPtr<Q> > mapT;
        typedef std::pair< unsigned long, SharedPtr<Q> > pairT;
        mapT cache;
        
        // Turns (n,lx) into key
        inline unsigned long key(Level n, long lx) const {
            return (n<<25) | (lx+16777216);
        }
        
        // Turns (n,displacement) into key
        template <int NDIM>
        inline unsigned long key(Level n, const Key<NDIM>& d) const {
            MADNESS_ASSERT((6+NDIM*4) <= sizeof(unsigned long)*8);
            unsigned long k = n<<2;
            for (int i=0; i<NDIM; i++) k = (k<<4) | (d.translation()[i]+7);
            return k;
        }
        
    public:
        SimpleCache() : cache() {};
        
        SimpleCache(const SimpleCache& c) : cache(c.cache) {};
        SimpleCache& operator=(const SimpleCache& c) {
            if (this != &c) {
                cache.clear();
                cache = c.cache;
            }
            return *this;
        }
        
        /// If (n,index) is present return pointer to cached value, otherwise return NULL
        template <typename indexT>
        inline const Q* getptr(Level n,  const indexT& index) const {
            typename mapT::const_iterator test = cache.find(key(n,index));
            if (test == cache.end()) return 0;
            return test->second.get();
        }
        

        /// Set value associated with (n,index)
        template <typename indexT>
        inline void set(Level n, const indexT& index, const Q& val) {
            cache.insert(pairT(key(n,index),SharedPtr<Q>(new Q(val))));
        }
    };

    template <typename Q>
    struct ConvolutionData1D {
        Tensor<Q> R, T;  ///< R=ns, T=T part of ns
        Tensor<Q> RU, RVT, TU, TVT; ///< SVD approximations to R and T
        Tensor<typename Tensor<Q>::scalar_type> Rs, Ts;
        double Rnorm, Tnorm, Rnormf, Tnormf, NSnormf;

        ConvolutionData1D(const Tensor<Q>& R, const Tensor<Q>& T) : R(R), T(T) {
            make_approx(R, RU, Rs, RVT, Rnorm);
            make_approx(T, TU, Ts, TVT, Tnorm);
            Rnormf = R.normf();
            Tnormf = T.normf();
            int k = T.dim[0];
            Tensor<Q> NS = copy(R);
            for (int i=0; i<k; i++) 
                for (int j=0; j<k; j++)
                    NS(i,j) = 0.0;
            NSnormf = NS.normf();
        }

        void make_approx(const Tensor<Q>& R, 
                         Tensor<Q>& RU, Tensor<typename Tensor<Q>::scalar_type>& Rs, Tensor<Q>& RVT, double& norm) {
            int n = R.dim[0];
            svd(R, &RU, &Rs, &RVT);
            for (int i=0; i<n; i++) {
                for (int j=0; j<n; j++) {
                    RVT(i,j) *= Rs[i];
                }
            }
            for (int i=n-1; i>1; i--) { // Form cumulative sum of norms
                Rs[i-1] += Rs[i];
            }
            
            norm = Rs[0];
            if (Rs[0]>0.0) { // Turn into relative errors
                double rnorm = 1.0/norm;
                for (int i=0; i<n; i++) {
                    Rs[i] *= rnorm;
                }
            }
        }
    };


    /// Provides the common functionality/interface of all 1D convolutions

    /// Derived classes must implement rnlp, issmall, natural_level
    template <typename Q>
    class Convolution1D {
    public:
        int k;
        int npt;
        double sign;
        Tensor<double> quad_x;
        Tensor<double> quad_w;
        Tensor<double> c;
        Tensor<double> hgT;
        Tensor<double> hgT2k;

        mutable SimpleCache< Tensor<Q> > rnlp_cache;
        mutable SimpleCache< Tensor<Q> > rnlij_cache;
        mutable SimpleCache< ConvolutionData1D<Q> > ns_cache;
        
        Convolution1D() : k(-1), npt(0), sign(1.0) {}; 

        virtual ~Convolution1D() {};

        Convolution1D(int k, int npt, double sign=1.0) 
            : k(k)
            , npt(npt)
            , sign(sign)
            , quad_x(npt)
            , quad_w(npt)
        {
            MADNESS_ASSERT(autoc(k,&c));

            gauss_legendre(npt,0.0,1.0,quad_x.ptr(),quad_w.ptr());
            MADNESS_ASSERT(two_scale_hg(k,&hgT));
            hgT = transpose(hgT);
            MADNESS_ASSERT(two_scale_hg(2*k,&hgT2k));
            hgT2k = transpose(hgT2k);

            // Cannot construct the coefficients here since the
            // derived class is not yet constructed so cannot call
            // (even indirectly) a virtual method
        }

        /// Compute the projection of the operator onto the double order polynomials
        virtual Tensor<Q> rnlp(long n, long lx) const = 0;

        /// Returns true if the block is expected to be small
        virtual bool issmall(long n, long lx) const = 0;

        /// Returns the level for projection
        virtual Level natural_level() const {return 13;}

        /// Computes the transition matrix elements for the convolution for n,l
        
        /// Returns the tensor 
        /// \code
        ///   r(i,j) = int(K(x-y) phi[n0](x) phi[nl](y), x=0..1, y=0..1)
        /// \endcode
        /// This is computed from the matrix elements over the correlation
        /// function which in turn are computed from the matrix elements
        /// over the double order legendre polynomials.
        const Tensor<Q>& rnlij(long n, long lx) const {
            const Tensor<Q>* p=rnlij_cache.getptr(n,lx);
            if (p) return *p;
    
            long twok = 2*k;
            Tensor<Q>  R(2*twok);
            R(Slice(0,twok-1)) = get_rnlp(n,lx-1);
            R(Slice(twok,2*twok-1)) = get_rnlp(n,lx);
            R.scale(pow(0.5,0.5*n));        
            R = inner(c,R);
            // Enforce symmetry because it seems important ... is it?
            // What about a complex exponents?
//             if (lx == 0) 
//                 for (int i=0; i<k; i++) 
//                     for (int j=0; j<i; j++) 
//                         R(i,j) = R(j,i) = ((i+j)&1) ? 0.0 : 0.5*(R(i,j)+R(j,i));

            rnlij_cache.set(n,lx,R);
            return *rnlij_cache.getptr(n,lx);
        };

        /// Returns a pointer to the cached nonstandard form of the operator
        const ConvolutionData1D<Q>* nonstandard(long n, long lx) const {
            const ConvolutionData1D<Q>* p = ns_cache.getptr(n,lx);
            if (p) return p;

            PROFILE_MEMBER_FUNC(ConvolutionData1D);
            
            Tensor<Q> R(2*k,2*k), T;
            if (issmall(n, lx)) {
                T = Tensor<Q>(k,k);
            }
            else {
                long lx2 = lx*2;
                Slice s0(0,k-1), s1(k,2*k-1);
                
                R(s0,s0) = R(s1,s1) = rnlij(n+1,lx2);
                R(s1,s0) = rnlij(n+1,lx2+1);
                R(s0,s1) = rnlij(n+1,lx2-1);
                
                R = transform(R,hgT);
                // Enforce symmetry because it seems important ... what about complex?????????
//                 if (lx == 0) 
//                     for (int i=0; i<2*k; i++) 
//                         for (int j=0; j<i; j++) 
//                             R(i,j) = R(j,i) = ((i+j)&1) ? 0.0 : 0.5*(R(i,j)+R(j,i));
                
                R = transpose(R);
                T = copy(R(s0,s0));
            }
            
            ns_cache.set(n,lx,ConvolutionData1D<Q>(R,T));

            return ns_cache.getptr(n,lx);
        };

        const Tensor<Q>& get_rnlp(long n, long lx) const 
        {
            const Tensor<Q>* p=rnlp_cache.getptr(n,lx);
            if (p) return *p;

            long twok = 2*k;
            Tensor<Q> r;
            
            if (issmall(n, lx)) {
                r = Tensor<Q>(twok);
            }
            else if (n < natural_level()) {
                Tensor<Q>  R(2*twok);
                R(Slice(0,twok-1)) = get_rnlp(n+1,2*lx);
                R(Slice(twok,2*twok-1)) = get_rnlp(n+1,2*lx+1);
                
                R = transform(R, hgT2k);
                r = copy(R(Slice(0,twok-1)));
            }
            else {
                r = rnlp(n, lx);
            }
            
            rnlp_cache.set(n, lx, r);
            return *rnlp_cache.getptr(n,lx);
        }
    };


    // To test generic convolution by comparing with GaussianConvolution1D
    template <typename Q> 
    class GaussianGenericFunctor {
    private:
        Q coeff;
        double exponent;
    public:
        GaussianGenericFunctor(Q coeff, double exponent)
            : coeff(coeff), exponent(exponent) {}

        Q operator()(double x) const {
            return coeff*exp(-exponent*x*x);
        }
    };


    /// Generic 1D convolution using brute force (i.e., slow) adaptive quadrature for rnlp

    /// Calls op(x) with x in *simulation coordinates* to evaluate the function.
    template <typename Q, typename opT>
    class GenericConvolution1D : public Convolution1D<Q> {
    private:
        opT op;
        long maxl; //< At natural level is l beyond which operator is zero
    public:

        GenericConvolution1D() {}

        GenericConvolution1D(int k, const opT& op) 
            : Convolution1D<Q>(k, 20), op(op), maxl(LONG_MAX-1)
        {
            PROFILE_MEMBER_FUNC(GenericConvolution1D);
            
            // For efficiency carefully compute outwards at the "natural" level
            // until several successive boxes are determined to be zero.  This 
            // then defines the future range of the operator and also serves 
            // to precompute the values used in the rnlp cache.

            Level natl = this->natural_level();
            int nzero = 0;
            for (long lx=0; lx<(1L<<natl); lx++) {
                const Tensor<Q>& rp = this->get_rnlp(natl, lx);
                const Tensor<Q>& rm = this->get_rnlp(natl,-lx);
                if (rp.normf()<1e-12 && rm.normf()<1e-12) nzero++;
                if (nzero == 3) {
                    maxl = lx-2;
                    break;
                }
            }
        }

        struct Shmoo {
            typedef Tensor<Q> returnT;
            int n;
            long lx;
            const GenericConvolution1D<Q,opT>& q;

            Shmoo(int n, long lx, const GenericConvolution1D<Q,opT>* q)
                : n(n), lx(lx), q(*q) {}

            returnT operator()(double x) const {
                int twok = q.k*2;
                double fac = std::pow(0.5,n);
                double phix[twok];
                legendre_scaling_functions(x-lx,twok,phix);
                Q f = q.op(fac*x)*sqrt(fac);
                returnT v(twok);
                for (long p=0; p<twok; p++) v(p) += f*phix[p];
                return v;
            }
        };

        Tensor<Q> rnlp(long n, long lx) const {
            return adq1(lx, lx+1, Shmoo(n, lx, this), 1e-12, 
                        this->npt, this->quad_x.ptr(), this->quad_w.ptr(), 0);
        }

        bool issmall(long n, long lx) const {
            if (lx < 0) lx = 1 - lx;
            // Always compute contributions to nearest neighbor coupling
            // ... we are two levels below so 0,1 --> 0,1,2,3 --> 0,...,7
            if (lx <= 7) return false;
            
            n = this->natural_level()-n;
            if (n >= 0) lx = lx << n;
            else lx = lx >> n;

            return lx >= maxl;
        }
    };
        

    // For complex types return +1 as the sign and leave coeff unchanged
    template <typename Q, bool iscomplex> 
    struct munge_sign_struct {
        static double op(Q& coeff) {
            return 1.0;
        }
    };
        
    // For real types return actual sign and make coeff positive
    template <typename Q>
    struct munge_sign_struct<Q,false> {
        static typename Tensor<Q>::scalar_type op(Q& coeff) {
            if (coeff < 0.0) {
                coeff = -coeff;
                return -1.0;
            }
            else {
                return 1.0;
            }
        }
    };

    template <typename Q>
    typename Tensor<Q>::scalar_type munge_sign(Q& coeff) {
        return munge_sign_struct<Q, TensorTypeData<Q>::iscomplex>::op(coeff);
    }
            

    /// 1D Gaussian convolution with coeff and expnt given in *simulation* coordinates [0,1]
    template <typename Q>
    class GaussianConvolution1D : public Convolution1D<Q> {
    public:
        Q coeff;
        double expnt;

        GaussianConvolution1D() : Convolution1D<Q>(), coeff(0.0), expnt(0.0) {}; 

        GaussianConvolution1D(int k, Q coeff, double expnt, double sign=1.0)
            : Convolution1D<Q>(k,k+11,sign), coeff(coeff), expnt(expnt)
        {}

        /// Compute the projection of the operator onto the double order polynomials
        
        /// The returned reference is to a cached tensor ... if you want to
        /// modify it, take a copy first.
        ///
        /// Return in \c v[p] \c p=0..2*k-1
        /// \code
        /// r(n,l,p) = 2^(-n) * int(K(2^(-n)*(z+l)) * phi(p,z), z=0..1)
        /// \endcode
        /// The kernel is coeff*exp(-expnt*z^2).  This is equivalent to
        /// \code
        /// r(n,l,p) = 2^(-n)*coeff * int( exp(-beta*z^2) * phi(p,z-l), z=l..l+1)
        /// \endcode
        /// where 
        /// \code
        /// beta = alpha * 2^(-2*n) 
        /// \endcode
        Tensor<Q> rnlp(long n, long lx) const {
            PROFILE_MEMBER_FUNC(GaussianConvolution1D);
            int twok = 2*this->k;
            Tensor<Q> v(twok);       // Can optimize this away by passing in
            
            long lkeep = lx;
            if (lx<0) lx = -lx-1;
            
            /* Apply high-order Gauss Legendre onto subintervals
       
               coeff*int(exp(-beta(x+l)**2) * phi[p](x),x=0..1);
               
               The translations internally considered are all +ve, so
               signficant pieces will be on the left.  Finish after things
               become insignificant.  
               
               The resulting coefficients are accurate to about 1e-20.
            */
            
            // Rescale expnt & coeff onto level n so integration range
            // is [l,l+1]
            Q scaledcoeff = coeff*pow(sqrt(0.5),double(n));
            double beta = expnt * pow(0.25,double(n));
            
            // Subdivide interval into nbox boxes of length h
            // ... estimate appropriate size from the exponent.  A
            // Gaussian with real-part of the exponent beta falls in
            // magnitude by a factor of 1/e at x=1/sqrt(beta), and by
            // a factor of e^-49 ~ 5e-22 at x=7/sqrt(beta).  So, if we
            // use a box of size 1/sqrt(beta) we will need at most 7
            // boxes.  Incorporate the coefficient into the screening
            // since it may be large.  We can represent exp(-x^2) over
            // [l,l+1] with a polynomial of order 21 over [l,l+1] to a
            // relative precision of better than machine precision for
            // l=0,1,2,3 and for l>3 the absolute error is less than
            // 1e-23.  We want to compute matrix elements with
            // polynomials of order 2*k-1, so the total order is
            // 2*k+20, which can be integrated with a quadrature rule
            // of npt=k+11.  npt is set in the constructor.
            
            double h = 1.0/sqrt(beta);  // 2.0*sqrt(0.5/beta);
            long nbox = long(1.0/h);
            if (nbox < 1) nbox = 1;
            h = 1.0/nbox;
            
            // Find argmax such that h*scaledcoeff*exp(-argmax)=1e-22 ... if
            // beta*xlo*xlo is already greater than argmax we can neglect this
            // and subsequent boxes
            double argmax = std::abs(log(1e-22/std::abs(scaledcoeff*h)));
            
            for (long box=0; box<nbox; box++) {
                double xlo = box*h + lx;
                if (beta*xlo*xlo > argmax) break;
                for (long i=0; i<this->npt; i++) {
#ifdef IBMXLC
                    double phix[70];
#else
                    double phix[twok];
#endif
                    double xx = xlo + h*this->quad_x(i);
                    Q ee = scaledcoeff*exp(-beta*xx*xx)*this->quad_w(i)*h;
                    legendre_scaling_functions(xx-lx,twok,phix);
                    for (long p=0; p<twok; p++) v(p) += ee*phix[p];
                }
            }
            
            if (lkeep < 0) {
                /* phi[p](1-z) = (-1)^p phi[p](z) */
                for (long p=1; p<twok; p+=2) v(p) = -v(p);
            }

            return v;
        };

        /// Returns true if the block is expected to be small
        bool issmall(long n, long lx) const { 
            double beta = expnt*(pow(0.25,double(n)));
            long ll;
            if (lx > 0)
                ll = lx - 1;
            else if (lx < 0)
                ll = -1 - lx;
            else
                ll = 0;
            
            return (beta*ll*ll > 49.0);      // 49 -> 5e-22     69 -> 1e-30
        };
    };
    

    template <typename Q, int NDIM>
    struct SeparatedConvolutionInternal {
        double norm;
        const ConvolutionData1D<Q>* ops[NDIM];
    };

    template <typename Q, int NDIM>
    struct SeparatedConvolutionData {
        std::vector< SeparatedConvolutionInternal<Q,NDIM> > muops;
        double norm;

        SeparatedConvolutionData(int rank) : muops(rank) {}
        SeparatedConvolutionData(const SeparatedConvolutionData<Q,NDIM>& q) {
            muops = q.muops;
            norm = q.norm;
        }
    };


    /// Convolutions in separated form (including Gaussian)
    template <typename Q, int NDIM>
    class SeparatedConvolution : public WorldObject< SeparatedConvolution<Q,NDIM> > {
    public:
        typedef Q opT;  ///< The apply function uses this to infer resultT=opT*inputT
        bool doleaves;  ///< If should be applied to leaf coefficients ... false by default
        bool dowiden0;  ///< If true widen operator if diagonal term significant ... false by default
        bool dowiden1;  ///< If true widen operator if make off-diaginal contrib ... false by default
    private:
        const int k;
        const int rank;
        const std::vector<long> vk;
        const std::vector<long> v2k;
        mutable Tensor<Q> work5;
        const std::vector<Slice> s0;
        mutable std::vector< SharedPtr< Convolution1D<Q> > > ops;
        mutable SimpleCache< SeparatedConvolutionData<Q,NDIM> > data;
        mutable double nflop[64];

        struct Transformation {
            long r;     // Effective rank of transformation
            const Q* U; // Ptr to matrix
            const Q* VT; 
        };

        template <typename T, typename R>
        void apply_transformation(Level n, long dimk,
                                  const Transformation trans[NDIM],
                                  const Tensor<T>& f,
                                  Tensor<R>& work1,
                                  Tensor<R>& work2,
                                  Tensor<Q>& work3,
                                  const double musign,
                                  Tensor<R>& result) const {
            
            PROFILE_MEMBER_FUNC(SeparatedConvolution);
            long size = 1;
            for (int i=0; i<NDIM; i++) size *= dimk;
            long dimi = size/dimk;

            R* RESTRICT w1=work1.ptr();
            R* RESTRICT w2=work2.ptr();
            Q* RESTRICT w3=work3.ptr();

            const Q* U;

            U = (trans[0].r == dimk) ? trans[0].U : shrink(dimk,dimk,trans[0].r,trans[0].U,w3);
            mTxmq(dimi, trans[0].r, dimk, w1, f.ptr(), U);
            size = trans[0].r * size / dimk;
            dimi = size/dimk;
            nflop[n] += dimi*trans[0].r*dimk*2.0;
            for (int d=1; d<NDIM; d++) {
                U = (trans[d].r == dimk) ? trans[d].U : shrink(dimk,dimk,trans[d].r,trans[d].U,w3);
                mTxmq(dimi, trans[d].r, dimk, w2, w1, U);
                size = trans[d].r * size / dimk;
                dimi = size/dimk;
                std::swap(w1,w2);
                nflop[n] += dimi*trans[d].r*dimk*2.0;
            }

            // If all blocks are full rank we can skip the transposes
            bool doit = false;
            for (int d=0; d<NDIM; d++) doit = doit || trans[d].VT;

            if (doit) {
                for (int d=0; d<NDIM; d++) {
                    if (trans[d].VT) {
                        dimi = size/trans[d].r;
                        mTxmq(dimi, dimk, trans[d].r, w2, w1, trans[d].VT);
                        nflop[n] += dimi*trans[d].r*dimk*2.0;
                        size = dimk*size/trans[d].r;
                    }
                    else {
                        fast_transpose(dimk, dimi, w1, w2);
                    }
                    std::swap(w1,w2);
                }
            }

            // Assuming here that result is contiguous
            MADNESS_ASSERT(size == result.size);
            R* RESTRICT p = result.ptr();
            if (musign > 0) {
                //for (long i=0; i<size; i++) p[i] += w1[i];
                aligned_add(size, p, w1);
            }
            else {
                //for (long i=0; i<size; i++) p[i] -= w1[i];
                aligned_sub(size, p, w1);
            }
        }

        /// Apply one of the separated terms, accumulating into the result
        template <typename T>
        void muopxv_fast(Level n,
                         const ConvolutionData1D<Q>* const ops[NDIM],
                         const Tensor<T>& f, const Tensor<T>& f0,
                         Tensor<TENSOR_RESULT_TYPE(T,Q)>& result,
                         Tensor<TENSOR_RESULT_TYPE(T,Q)>& result0,
                         double tol,
                         const double musign, 
                         Tensor<TENSOR_RESULT_TYPE(T,Q)>& work1,
                         Tensor<TENSOR_RESULT_TYPE(T,Q)>& work2) const {

            PROFILE_MEMBER_FUNC(SeparatedConvolution);
            Transformation trans[NDIM];

            const long twok = 2*k;

            double Rnorm = 1.0;
            for (int d=0; d<NDIM; d++) Rnorm *= ops[d]->Rnorm;
               
            tol = tol/(Rnorm*NDIM);  // Errors are relative within here

            // Determine rank of SVD to use or if to use the full matrix
            long break_even;
            if (NDIM==1) break_even = long(0.5*twok);
            else if (NDIM==2) break_even = long(0.6*twok);
            else if (NDIM==3) break_even=long(0.65*twok);
            else break_even=long(0.7*twok);
            for (int d=0; d<NDIM; d++) {
                long r;
                for (r=0; r<twok; r++) {
                    if (ops[d]->Rs[r] < tol) break;
                }
                if (r >= break_even) {
                    trans[d].r = twok;
                    trans[d].U = ops[d]->R.ptr();
                    trans[d].VT = 0;
                }
                else {
                    r += (r&1L);
                    trans[d].r = std::max(2L,r);
                    trans[d].U = ops[d]->RU.ptr();
                    trans[d].VT = ops[d]->RVT.ptr();
                }
            }
            apply_transformation(n, twok, trans, f, work1, work2, work5, musign, result);

            if (n > 0) {
                if (NDIM==1) break_even = long(0.5*k);
                else if (NDIM==2) break_even = long(0.6*k);
                else if (NDIM==3) break_even=long(0.65*k);
                else break_even=long(0.7*k);
                for (int d=0; d<NDIM; d++) {
                    long r;
                    for (r=0; r<k; r++) {
                        if (ops[d]->Ts[r] < tol) break;
                    }
                    if (r >= break_even) {
                        trans[d].r = k;
                        trans[d].U = ops[d]->T.ptr();
                        trans[d].VT = 0;
                    }
                    else {
                        r += (r&1L);
                        trans[d].r = std::max(2L,r);
                        trans[d].U = ops[d]->TU.ptr();
                        trans[d].VT = ops[d]->TVT.ptr();
                    }
                }
                apply_transformation(n, k, trans, f0, work1, work2, work5, -musign, result0);
            }
        }

        /// Apply one of the separated terms, accumulating into the result

        /// !!! Keep this routine exactly consistent with muopxvT so that
        /// munorm converges correctly
        template <typename T>
        void muopxv(Level n,
                    const ConvolutionData1D<Q>* const ops[NDIM],
                    const Tensor<T>& f, const Tensor<T>& f0,
                    Tensor<TENSOR_RESULT_TYPE(T,Q)>& result,
                    const double tol,
                    const double musign) const {
            PROFILE_MEMBER_FUNC(SeparatedConvolution);
            // Temporaries can be optimized away
            Tensor<TENSOR_RESULT_TYPE(T,Q)> tmp = inner(f,ops[0]->R,0,0);
            for (int d=1; d<NDIM; d++) {
                tmp = inner(tmp,ops[d]->R,0,0);
            }
            result.gaxpy(1.0,tmp,musign);
            
            if (n > 0) {
                tmp = inner(f0,ops[0]->T,0,0);
                for (int d=1; d<NDIM; d++) {
                    tmp = inner(tmp,ops[d]->T,0,0);
                }
                result(s0).gaxpy(1.0,tmp,-musign);
            }
        }

        /// Apply transpose of one of the separated terms, accumulating into the result

        /// This is only used when computing the actual 2-norm by the power method
        /// !!! Keep this routine exactly consistent with muopxv so that
        /// munorm converges correctly
        template <typename T>
        void muopxvT(Level n,
                     const ConvolutionData1D<Q>* ops[],
                     const Tensor<T>& f, const Tensor<T>& f0,
                     Tensor<TENSOR_RESULT_TYPE(T,Q)>& result,
                     const double tol,
                     const double musign) const {
            PROFILE_MEMBER_FUNC(SeparatedConvolution);
            // Temporaries can be optimized away
            Tensor<TENSOR_RESULT_TYPE(T,Q)> tmp = inner(f,ops[0]->R,0,1);
            for (int d=1; d<NDIM; d++) {
                tmp = inner(tmp,ops[d]->R,0,1);
            }
            result.gaxpy(1.0,tmp,musign);
            
            if (n > 0) {
                tmp = inner(f0,ops[0]->T,0,1); // Slice+copy can be optimized away
                for (int d=1; d<NDIM; d++) {
                    tmp = inner(tmp,ops[d]->T,0,1);
                }
                result(s0).gaxpy(1.0,tmp,-musign);
            }
        }


        /// Computes the 2-norm of one of the separated terms
        double munorm(Level n, const ConvolutionData1D<Q>* ops[]) const {
            PROFILE_MEMBER_FUNC(SeparatedConvolution);
            Tensor<Q> f(v2k), f0, ff(v2k);
            
            double tol = 1e-20;
            
            f.fillrandom();
            f.scale(1.0/f.normf());
            double evalp = 1.0, eval, ratio=99.0;
            for (int iter=0; iter<100; iter++) {
                ff.fill(0.0); 
                f0 = copy(f(s0));
                muopxv(n,ops,f,f0,ff,tol,1.0);
                f.fill(0.0);  
                f0 = copy(ff(s0));
                muopxvT(n,ops,ff,f0,f,tol,1.0);
                
                eval = f.normf();
                if (eval == 0.0) break;
                f.scale(1.0/eval);
                eval = sqrt(eval);
                ratio = eval/evalp;
                //std::printf("munorm: %d %10.2e %10.2e %10.2e \n", iter, eval, evalp, ratio);
                if (iter>0 && ratio<1.2 && ratio>0.9999) break; // 1.2 was 1.02;  >0.9999 was >=1.0
                if (iter>10 && eval<tol) break;
                evalp = eval;
                if (iter == 99) throw "munorm failed";
            }
            return eval*ratio;
        }

        /// Computes the Frobenius norm of one of the separated terms
        double munorm2(Level n, const ConvolutionData1D<Q>* ops[]) const {
            PROFILE_MEMBER_FUNC(SeparatedConvolution);
            double prodR=1.0, prodT=1.0;
            for (int d=0; d<NDIM; d++) {
                prodR *= ops[d]->Rnormf;
                prodT *= ops[d]->Tnormf;

            }
            if (n) prodR = sqrt(prodR*prodR - prodT*prodT);

            if (prodR < 1e-8*prodT) {
                double prod=1.0, sum=0.0;
                for (int d=0; d<NDIM; d++) {
                    double a = ops[d]->NSnormf;
                    double b = ops[d]->Tnormf;
                    double aa = std::min(a,b);
                    double bb = std::max(a,b);
                    prod *= bb;
                    if (bb > 0.0) sum +=(aa/bb);
                }
                if (n) prod *= sum;
                prodR = prod;
            }

            return prodR;
        }


        const SeparatedConvolutionInternal<Q,NDIM> getmuop(int mu, Level n, const Key<NDIM>& disp) const {
            PROFILE_MEMBER_FUNC(SeparatedConvolution);
            SeparatedConvolutionInternal<Q,NDIM> op;
            for (int d=0; d<NDIM; d++) {
                op.ops[d] = ops[mu]->nonstandard(n, disp.translation()[d]);
            }

            op.norm = munorm2(n, op.ops);
            //op.norm = munorm(n, op.ops);

//             double newnorm = munorm2(n, op.ops);
//             // This rescaling empirically based upon BSH separated expansion
//             // ... needs more testing.  OK also for TDSE.
//             // All is good except for some 000 blocks which are up to sqrt(k^d) off.
//             for (int d=0; d<NDIM; d++)  {
//                 if (disp[d] == 0) newnorm *= 0.5;
//                 else if (std::abs(disp[d]) == 1) newnorm *= 0.8;
//             }
//            double oldnorm = munorm(n, op.ops);
//             if (oldnorm > 1e-13 && (newnorm < 0.5*oldnorm || newnorm > 2.0*oldnorm) )
//                 print("munorm", n, disp, mu, newnorm, oldnorm, newnorm/oldnorm);

            return op;
        }

        
        /// Returns pointer to cached operator
        const SeparatedConvolutionData<Q,NDIM>* getop(Level n, const Key<NDIM>& d) const {
            PROFILE_MEMBER_FUNC(SeparatedConvolution);
            const SeparatedConvolutionData<Q,NDIM>* p = data.getptr(n,d);
            if (p) return p;

            SeparatedConvolutionData<Q,NDIM> op(rank);
            for (int mu=0; mu<rank; mu++) {
                op.muops[mu] = getmuop(mu, n, d);
            }

            double norm = 0.0;
            for (int mu=0; mu<rank; mu++) {
                const double munorm = op.muops[mu].norm;
                norm += munorm*munorm;
            }
            op.norm = sqrt(norm);
            data.set(n, d, op);
            return data.getptr(n,d);
        }

        void check_cubic() {
            // !!! NB ... cell volume obtained from global defaults
            const Tensor<double>& cell_width = FunctionDefaults<NDIM>::get_cell_width();
            // Check that the cell is cubic since currently is assumed
            for (long d=1; d<NDIM; d++) {
                MADNESS_ASSERT(fabs(cell_width(d)-cell_width(0L)) < 1e-14*cell_width(0L));
            }
        }


    public:

        // For general convolutions
        SeparatedConvolution(World& world,
                             long k, 
                             std::vector< SharedPtr< Convolution1D<Q> > > ops,
                             bool doleaves=false)
            : WorldObject< SeparatedConvolution<Q,NDIM> >(world)
            , doleaves(doleaves)
            , dowiden0(false)
            , dowiden1(false)
            , k(k)
            , rank(ops.size())
            , vk(NDIM,k)
            , v2k(NDIM,2*k)
            , work5(2*k,2*k)
            , s0(std::max(2,NDIM),Slice(0,k-1))
            , ops(ops)
        {

            check_cubic();

            for (int i=0; i<64; i++) nflop[i] = 0.0;

            this->process_pending();
        }

        /// Contstructor for Gaussian Convolutions (mostly for backward compatability)
        SeparatedConvolution(World& world, 
                             int k,
                             const Tensor<Q>& coeff, const Tensor<double>& expnt,
                             bool doleaves = false) 
            : WorldObject< SeparatedConvolution<Q,NDIM> >(world)
            , doleaves(doleaves)
            , dowiden0(false)
            , dowiden1(false)
            , k(k)
            , rank(coeff.dim[0])
            , vk(NDIM,k)
            , v2k(NDIM,2*k)
            , work5(2*k,2*k)
            , s0(std::max(2,NDIM),Slice(0,k-1))
            , ops(coeff.dim[0]) 
        {
            check_cubic();
            double width = FunctionDefaults<NDIM>::get_cell_width()(0L);
            for (int i=0; i<64; i++) nflop[i] = 0.0;

            for (int i=0; i<rank; i++) {
                Q c = coeff(i);
                double sign = munge_sign(c);
                c = std::pow(c, 1.0/NDIM);
                ops[i] = SharedPtr< Convolution1D<Q> >(new GaussianConvolution1D<Q>(k, 
                                                                                    c*width, 
                                                                                    expnt(i)*width*width,
                                                                                    sign));
            }
        }

        const double* get_nflop() const {
            return nflop;
        }

        double norm(Level n, const Key<NDIM>& d) const {
            return getop(n, d)->norm;
        }

        template <typename T>
        Tensor<TENSOR_RESULT_TYPE(T,Q)> apply(const Key<NDIM>& source,
                                              const Key<NDIM>& shift,
                                              const Tensor<T>& coeff,
                                              double tol) const {
            PROFILE_MEMBER_FUNC(SeparatedConvolution);
            typedef TENSOR_RESULT_TYPE(T,Q) resultT;
            const Tensor<T>* input = &coeff;
            Tensor<T> dummy;

            if (coeff.dim[0] == k) {
                // This processes leaf nodes with only scaling
                // coefficients ... FuncImpl::apply by default does not
                // apply the operator to these since for smoothing operators
                // it is not necessary.  It is necessary for operators such
                // as differentiation and time evolution and will also occur
                // if the application of the operator widens the tree.
                dummy = Tensor<T>(v2k);
                dummy(s0) = coeff;
                input = &dummy;
            }
            else {
                MADNESS_ASSERT(coeff.dim[0]==2*k);
            }

            tol = tol/rank; // Error is per separated term

            const SeparatedConvolutionData<Q,NDIM>* op = getop(source.level(), shift);
            //print("sepop",source,shift,op->norm,tol);
            Tensor<resultT> r(v2k), r0(vk);
            Tensor<resultT> work1(v2k,false), work2(v2k,false);

            const Tensor<T> f0 = copy(coeff(s0));
            for (int mu=0; mu<rank; mu++) {
                const SeparatedConvolutionInternal<Q,NDIM>& muop =  op->muops[mu];
                //print(source, shift, mu, muop.norm);
                if (muop.norm > tol) {
                    muopxv_fast(source.level(), muop.ops, *input, f0, r, r0, tol, ops[mu]->sign, 
                                work1, work2);
                    //muopxv(source.level(), muop.ops, *input, f0, r, tol, ops[mu]->sign);
                }
            }
            r(s0).gaxpy(1.0,r0,1.0);
            return r;
        }
        
    };


    /// Factory function generating separated kernel for convolution with 1/r.
    template <typename Q, int NDIM>
    SeparatedConvolution<Q,NDIM> CoulombOperator(World& world,
                                                 long k,
                                                 double lo,
                                                 double eps) {
        const Tensor<double>& cell_width = FunctionDefaults<NDIM>::get_cell_width();
        double hi = cell_width.normf(); // Diagonal width of cell
        const double pi = 3.14159265358979323846264338328;

        // bsh_fit generates representation for 1/4Pir but we want 1/r
        // so have to scale eps by 1/4Pi
        Tensor<double> coeff, expnt;
        bsh_fit(0.0, lo, hi, eps/(4.0*pi), &coeff, &expnt, false);
        coeff.scale(4.0*pi);
        return SeparatedConvolution<Q,NDIM>(world, k, coeff, expnt);
    }

    /// Factory function generating separated kernel for convolution with 1/r.
    template <typename Q, int NDIM>
    SeparatedConvolution<Q,NDIM>* CoulombOperatorPtr(World& world,
                                                     long k,
                                                     double lo,
                                                     double eps) {
        const Tensor<double>& cell_width = FunctionDefaults<NDIM>::get_cell_width();
        double hi = cell_width.normf(); // Diagonal width of cell
        const double pi = 3.14159265358979323846264338328;
        
        // bsh_fit generates representation for 1/4Pir but we want 1/r
        // so have to scale eps by 1/4Pi
        Tensor<double> coeff, expnt;
        bsh_fit(0.0, lo, hi, eps/(4.0*pi), &coeff, &expnt, false);
        coeff.scale(4.0*pi);
        return new SeparatedConvolution<Q,NDIM>(world, k, coeff, expnt);
    }

    /// Factory function generating separated kernel for convolution with exp(-mu*r)/(4*pi*r)
    template <typename Q, int NDIM>
    SeparatedConvolution<Q,NDIM> BSHOperator(World& world,
                                             double mu,
                                             long k,
                                             double lo,
                                             double eps) {
        const Tensor<double>& cell_width = FunctionDefaults<NDIM>::get_cell_width();
        double hi = cell_width.normf(); // Diagonal width of cell
        Tensor<double> coeff, expnt;
        bsh_fit(mu, lo, hi, eps, &coeff, &expnt, false);
        return SeparatedConvolution<Q,NDIM>(world, k, coeff, expnt);
    }

    /// Factory function generating separated kernel for convolution with exp(-mu*r)/(4*pi*r)
    template <typename Q, int NDIM>
    SeparatedConvolution<Q,NDIM>* BSHOperatorPtr(World& world,
                                                 double mu,
                                                 long k,
                                                 double lo,
                                                 double eps) {
        const Tensor<double>& cell_width = FunctionDefaults<NDIM>::get_cell_width();
        double hi = cell_width.normf(); // Diagonal width of cell
        Tensor<double> coeff, expnt;
        bsh_fit(mu, lo, hi, eps, &coeff, &expnt, false);
        return new SeparatedConvolution<Q,NDIM>(world, k, coeff, expnt);
    }

    namespace archive {
        template <class Archive, class T, int NDIM>
        struct ArchiveLoadImpl<Archive,const SeparatedConvolution<T,NDIM>*> {
            static inline void load(const Archive& ar, const SeparatedConvolution<T,NDIM>*& ptr) {
                WorldObject< SeparatedConvolution<T,NDIM> >* p;
                ar & p;
                ptr = static_cast< const SeparatedConvolution<T,NDIM>* >(p);
            }
        };
        
        template <class Archive, class T, int NDIM>
        struct ArchiveStoreImpl<Archive,const SeparatedConvolution<T,NDIM>*> {
            static inline void store(const Archive& ar, const SeparatedConvolution<T,NDIM>*const& ptr) {
                ar & static_cast< const WorldObject< SeparatedConvolution<T,NDIM> >* > (ptr);
            }
        };
    }

}




#endif
