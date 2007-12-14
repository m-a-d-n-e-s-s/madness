#ifndef MAD_OPERATOR_H
#define MAD_OPERATOR_H

#include <mra/mra.h>

namespace madness {

    /// Simplified interface around hash_map to cache stuff
    template <typename Q>
    class SimpleCache {
    private:
        typedef HASH_MAP_NAMESPACE::hash_map< unsigned long, Q> mapT;
        typedef std::pair< unsigned long, Q> pairT;
        mapT cache;
        const typename mapT::const_iterator end;
        
        // Turns (n,lx) into key
        inline unsigned long key(Level n, long lx) const {
            return (n<<6) | (lx+8);
        }
        
        // Turns (n,displacement) into key
        template <int NDIM>
        inline unsigned long key(Level n, const Displacement<NDIM>& d) const {
            MADNESS_ASSERT((6+NDIM*4) <= sizeof(unsigned long)*8);
            unsigned long k = n<<2;
            for (int i=0; i<NDIM; i++) k = (k<<4) | (d[i]+8);
            return k;
        }
        
    public:
        SimpleCache() : cache(), end(const_cast<const mapT*>(&cache)->end()) {};
        
        SimpleCache(const SimpleCache& c) : cache(c.cache), end(const_cast<const mapT*>(&cache)->end()) {};
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
            if (test == end) return 0;
            return &((*test).second);
        }
        

        /// Set value associated with (n,index)
        template <typename indexT>
        inline void set(Level n, const indexT& index, const Q& val) {
            cache.insert(pairT(key(n,index),val));
        }
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
        mutable SimpleCache< Tensor<Q> > rnlij_cache;
        mutable SimpleCache< Tensor<Q> > ns_cache;
        mutable SimpleCache< Tensor<Q> > ns_T_cache;
        
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
        Tensor<Q> rnlp(long n, long l) const {
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
            if (nbox < 1) nbox = 1;       // If the exponent is complex ??
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
            R(Slice(0,twok-1)) = rnlp(n,lx-1);
            R(Slice(twok,2*twok-1)) = rnlp(n,lx);
            R.scale(pow(0.5,0.5*n));        
            R = inner(c,R);
            // Enforce symmetry because it seems important ... is it?
            // What about a complex exponents?
            if (lx == 0) 
                for (int i=0; i<k; i++) 
                    for (int j=0; j<i; j++) 
                        R(i,j) = R(j,i) = ((i+j)&1) ? 0.0 : 0.5*(R(i,j)+R(j,i));
    
            rnlij_cache.set(n,lx,R);
            return *rnlij_cache.getptr(n,lx);
        };

        /// Returns a reference to the nonstandard form of the operator
        const Tensor<Q>& nonstandard(long n, long lx) const {
            const Tensor<Q> *p = ns_cache.getptr(n,lx);
            if (p) return *p;
    
            long lx2 = lx*2;
            Tensor<Q> R(2*k,2*k);
            Slice s0(0,k-1), s1(k,2*k-1);
            
            R(s0,s0) = R(s1,s1) = rnlij(n+1,lx2);
            R(s1,s0) = rnlij(n+1,lx2+1);
            R(s0,s1) = rnlij(n+1,lx2-1);
            
            R = transform(R,hgT);
            // Enforce symmetry because it seems important ... what about complex?????????
            if (lx == 0) 
                for (int i=0; i<2*k; i++) 
                    for (int j=0; j<i; j++) 
                        R(i,j) = R(j,i) = ((i+j)&1) ? 0.0 : 0.5*(R(i,j)+R(j,i));
            
            ns_cache.set(n,lx,transpose(R));
            return *(ns_cache.getptr(n,lx));
        };


        /// Returns a reference to the T block of the nonstandard form
        const Tensor<Q>& nonstandard_T(long n, long lx) const {
            const Tensor<Q> *p = ns_T_cache.getptr(n,lx);
            if (p) return *p;
            Slice s0(0,k-1);
            ns_T_cache.set(n,lx,copy(nonstandard(n,lx)(s0,s0)));
            return *(ns_T_cache.getptr(n,lx));
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
            
            return (beta*ll*ll > 49.0);      // 49 -> 5e-22
        };
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
        std::vector<Slice> s0;
        mutable SimpleCache<double> norms;


        /// Apply the specified block of the mu'th term to a vector accumulating into the result
        template <typename T>
        void muopxv(const long mu, 
                    Level n,
                    const Displacement<NDIM>& shift,
                    const Tensor<T>& f, Tensor<TENSOR_RESULT_TYPE(T,Q)>& result,
                    const double tol) const {

            Tensor<TENSOR_RESULT_TYPE(T,Q)> tmp = inner(f,ops[mu].nonstandard(n,shift[0]),0,0);
            for (int d=1; d<NDIM; d++) {
                tmp = inner(tmp,ops[mu].nonstandard(n,shift[d]),0,0);
            }
            result.gaxpy(1.0,tmp,signs[mu]);
            
            if (n > 0) {
                tmp = inner(copy(f(s0)),ops[mu].nonstandard_T(n,shift[0]),0,0);
                for (int d=1; d<NDIM; d++) {
                    tmp = inner(tmp,ops[mu].nonstandard_T(n,shift[d]),0,0);
                }
                result(s0).gaxpy(1.0,tmp,-signs[mu]);
            }
        }


        /// Apply the transpose of the specified block of the mu'th term to a vector accumulating into the result
        template <typename T>
        void muopxvT(const long mu, 
                     Level n,
                     const Displacement<NDIM>& shift,
                     const Tensor<T>& f, Tensor<TENSOR_RESULT_TYPE(T,Q)>& result,
                     const double tol) const {

            Tensor<TENSOR_RESULT_TYPE(T,Q)> tmp = inner(f,ops[mu].nonstandard(n,shift[0]),0,1);
            for (int d=1; d<NDIM; d++) {
                tmp = inner(tmp,ops[mu].nonstandard(n,shift[d]),0,1);
            }
            result.gaxpy(1.0,tmp,signs[mu]);
            
            if (n > 0) {
                tmp = inner(copy(f(s0)),ops[mu].nonstandard_T(n,shift[0]),0,1);
                for (int d=1; d<NDIM; d++) {
                    tmp = inner(tmp,ops[mu].nonstandard_T(n,shift[d]),0,1);
                }
                result(s0).gaxpy(1.0,tmp,-signs[mu]);
            }
        }

        /// Returns the norm of the specified block of the NS operator of the mu'th Gaussian
        double munorm(long mu, Level n, const Displacement<NDIM>& shift) const {
            std::vector<long> v2k(NDIM);
            for (int i=0; i<NDIM; i++) v2k[i] = 2*k;            
            Tensor<Q> f(v2k), ff(v2k);
            
            double tol = 1e-20;
            
            f.fillrandom();
            f.scale(1.0/f.normf());
            double evalp = 1.0, eval, ratio=99.0;
            for (int iter=0; iter<100; iter++) {
                ff.fill(0.0); 
                muopxv(mu,n,shift,f,ff,tol);
                f.fill(0.0);  
                muopxvT(mu,n,shift,ff,f,tol);
                
                eval = f.normf();
                if (eval == 0.0) break;
                f.scale(1.0/eval);
                eval = sqrt(eval);
                ratio = eval/evalp;
                std::printf("munorm: %d %10.2e %10.2e %10.2e \n", iter, eval, evalp, ratio);
                if (iter>0 && ratio<1.2 && ratio>=1.0) break; // 1.2 was 1.02
                if (iter>10 && eval<1e-30) break;
                evalp = eval;
                if (iter == 99) throw "munorm failed";
            }
            return eval*ratio;
        }




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
            , s0(NDIM)
        {
            for (int i=0; i<rank; i++) {
                Q coeff = coeffs(i);
                signs[i] = munge_sign(coeff); 
                coeff = std::pow(coeff,1.0/NDIM);
                ops[i] = GaussianConvolution<Q>(k, coeff, expnts(i));
            }

            for (int i=0; i<NDIM; i++) s0[i] = Slice(0,k-1);
            this->process_pending();
        }

        
        double norm(const Key<NDIM>& key, const Displacement<NDIM>& d) const {
            const double *n = norms.getptr(key.level(),d);
            if (n) return *n;
            double thenorm = munorm(0, key.level(), d);
            norms.set(key.level(), d, thenorm);
            return thenorm;
        }

        template <typename T>
        Tensor<TENSOR_RESULT_TYPE(T,Q)> apply(const Key<NDIM>& source,
                                              const Displacement<NDIM>& shift,
                                              const Tensor<T>& coeff) const {
            typedef TENSOR_RESULT_TYPE(T,Q) resultT;
            print("sepop",source,shift);
            std::vector<long> v2k(NDIM);
            for (int i=0; i<NDIM; i++) v2k[i] = 2*k;            
            
            Tensor<resultT> r = Tensor<resultT>(v2k);
            if (coeff.dim[0] != k) muopxv(0L, source.level(), shift, coeff, r, tol);
            return r;
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
