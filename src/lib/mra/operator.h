#ifndef MAD_OPERATOR_H
#define MAD_OPERATOR_H

#include <mra/mra.h>

namespace madness {

    /// Simplified interface around hash_map to cache coeffs
    template <typename Q>
    class CoeffCache {
    private:
        typedef HASH_MAP_NAMESPACE::hash_map< unsigned long, Tensor<Q> > mapT;
        typedef std::pair< unsigned long, Tensor<Q> > pairT;
        mapT cache;
        const typename mapT::iterator end;
        
        // Turns (n,lx) into key
        inline unsigned long key(long n, long lx) const {
            return unsigned((n<<16) + (lx+(1L<<15)));
        };
        
    public:
        CoeffCache() : cache(), end(cache.end()) {};
        
        CoeffCache(const CoeffCache& c) : cache(c.cache), end(cache.end()) {};
        CoeffCache& operator=(const CoeffCache& c) {
            if (this != &c) {
                cache.clear();
                cache = c.cache;
            }
            return *this;
        };
        
        /// If (n,lx) is present return pointer to cached value, otherwise return NULL
        inline Tensor<Q>* getptr(long n,  long lx) {
            typename mapT::iterator test = cache.find(key(n,lx));
            if (test == end) return 0;
            return &((*test).second);
        };
        

        /// Set value associated with (n,lx)
        inline void set(long n, long lx, const Tensor<double>& val) {
            cache.insert(pairT(key(n,lx),copy(val)));
        };
    };

    /// Stores info about 1D Gaussian convolution
    template <typename Q>
    class GaussianConvolution {
    public:
        int k;
        int npt;
        Q coeff;
        double expnt;
        Tensor<Q> c;
        Tensor<double> hgT;
        Tensor<double> quad_x;
        Tensor<double> quad_w;
        CoeffCache<Q> rnlij_cache;
        CoeffCache<Q> ns_cache;
        CoeffCache<Q> ns_T_cache;
        
        GaussianConvolution() : k(-1), npt(-1), coeff(0.0), expnt(0.0) {}; 
        GaussianConvolution(int k, Q coeff, double expnt)
            : k(k), npt(k+11), coeff(coeff), expnt(expnt)
        {
            MADNESS_ASSERT(autoc(k,&c));
            quad_x = Tensor<double>(npt);
            quad_w = Tensor<double>(npt);
            gauss_legendre(npt,0.0,1.0,quad_x.ptr(),quad_w.ptr());
            if (!two_scale_hg(k,&hgT)) throw "need to make a consistent exception framework2!";
            hgT = transpose(hgT);
        }

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
        Tensor<Q>& rnlp(long n, long lx) {
            int twok = 2*k;
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
            if (nbox < 1) nbox = 1;       // If the exponent is ??
            h = 1.0/nbox;
            
            // Find argmax such that h*scaledcoeff*exp(-argmax)=1e-22 ... if
            // beta*xlo*xlo is already greater than argmax we can neglect this
            // and subsequent boxes
            double argmax = fabs(log(1e-22/fabs(scaledcoeff*h)));
            
            for (long box=0; box<nbox; box++) {
                double xlo = box*h + lx;
                if (beta*xlo*xlo > argmax) break;
                for (long i=0; i<npt; i++) {
#ifdef IBMXLC
                    double phix[70];
#else
                    double phix[twok];
#endif
                    double xx = xlo + h*quad_x(i);
                    Q ee = scaledcoeff*exp(-beta*xx*xx)*quad_w(i)*h;
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
        
        const Tensor<Q>& rnlij(long n, long lx);
        const Tensor<Q>& nonstandard(long n, long lx);
        Tensor<Q>& nonstandard_T(long n, long lx);
        bool issmall(long n, long lx);
    };
    


    // For complex types return +1 as the sign and leave coeff unchanged
    template <typename Q, bool iscomplex> 
    struct munge_sign_struct {
        static Q op(Q& coeff) {
            return 1.0;
        }
    };
        
    // For real types return actual sign and make coeff positive
    template <typename Q>
    struct munge_sign_struct<Q,false> {
        static Q op(Q& coeff) {
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
    Q munge_sign(Q& coeff) {
        return munge_sign_struct<Q, TensorTypeData<Q>::iscomplex>::op(coeff);
    }
            
    /// Isotropic convolutions in separated form (presently Gaussian)
    template <typename Q, int NDIM>
    class SeparatedConvolution : public WorldObject< SeparatedConvolution<Q,NDIM> > {
    public:
        typedef Q opT;  ///< The apply function uses this to infer resultT=opT*inputT
    private:
        const int k;
        const double tol;
        const int rank;
        std::vector< GaussianConvolution<Q> > ops;
        std::vector<Q> signs;

    public:

        SeparatedConvolution(World& world,
                             long k, double tol,
                             const Tensor<Q>& coeffs, 
                             const Tensor<double>& expnts)    // restriction to double must be relaxed
            : WorldObject< SeparatedConvolution<Q,NDIM> >(world)
            , k(k)
            , tol(tol)
            , rank(coeffs.dim[0])
            , ops(rank)
            , signs(rank) 
        {
            for (int i=0; i<rank; i++) {
                Q coeff = coeffs(i);
                signs[i] = munge_sign(coeff); 
                coeff = std::pow(coeff,1.0/NDIM);
                ops[i] = GaussianConvolution<Q>(k, coeff, expnts(i));
            }
            this->process_pending();
        }

        
        double norm(const Key<NDIM>& key, const Displacement<NDIM>& d) const {
            if (d.distsq > 2) return 0.0;
            else return 1.0;
        }

        template <typename T>
        Tensor<TENSOR_RESULT_TYPE(T,Q)> apply(const Key<NDIM>& source,
                                              const Displacement<NDIM>& shift,
                                              const Tensor<T>& coeff) const {
            typedef TENSOR_RESULT_TYPE(T,Q) resultT;
            print("sepop",source,shift);
            std::vector<long> v2k(NDIM);
            for (int i=0; i<NDIM; i++) v2k[i] = 2*k;
            return Tensor<resultT>(v2k);
        }
        
//         double munorm(long mu, long n, long x, long y, long z);
//         double norm(long n, long x, long y, long z);
//         void opxv(long n, long x, long y, long z,
//                   const Tensor<double>& f, Tensor<double>& result,
//                   double tol);
//         void opxvt(long n, long x, long y, long z,
//                    const Tensor<double>& f, Tensor<double>& result,
//                    double tol);
//         void muopxv(const long mu, 
//                     const long n, const long x, const long y, const long z,
//                     const Tensor<double>& f, Tensor<double>& result,
//                     const double tol);
//         void muopxvt(long mu, long n, long x, long y, long z,
//                      const Tensor<double>& f, Tensor<double>& result,
//                      double tol);
    };


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
