#ifndef MAD_OPERATOR_H
#define MAD_OPERATOR_H

#include <mra/mra.h>
#include <tensor/mtxmq.h>

namespace madness {

    extern void bsh_fit(double mu, double lo, double hi, double eps, 
                        Tensor<double> *pcoeff, Tensor<double> *pexpnt, bool prnt=false);


    /// Simplified interface around hash_map to cache stuff

    /// Since insertions into STL containers have the nasty habit of
    /// invalidating iterators we actually store shared pointers
    template <typename Q>
    class SimpleCache {
    private:
        typedef HASH_MAP_NAMESPACE::hash_map< unsigned long, SharedPtr<Q> > mapT;
        typedef std::pair< unsigned long, SharedPtr<Q> > pairT;
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
        Tensor<Q> R, T;  // R=ns, T=T part of ns

        ConvolutionData1D(const Tensor<Q>& R, const Tensor<Q>& T) : R(R), T(T) {}
    };

    /// Stores info about 1D Gaussian convolution
    template <typename Q>
    class GaussianConvolution1D {
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
        mutable SimpleCache< ConvolutionData1D<Q> > ns_cache;
        
        GaussianConvolution1D() : k(-1), npt(-1), coeff(0.0), expnt(0.0) {}; 
        GaussianConvolution1D(int k, Q coeff, double expnt)
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
        Tensor<Q> rnlp(long n, long lx) const {
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

        /// Returns a pointer to the cached nonstandard form of the operator
        const ConvolutionData1D<Q>* nonstandard(long n, long lx) const {
            const ConvolutionData1D<Q>* p = ns_cache.getptr(n,lx);
            if (p) return p;
    
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

            R = transpose(R);
            Tensor<Q> T = copy(R(s0,s0));
            
            ns_cache.set(n,lx,ConvolutionData1D<Q>(R,T));
            return ns_cache.getptr(n,lx);
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


    /// Convolutions in separated form (presently Gaussian)
    template <typename Q, int NDIM>
    class SeparatedConvolution : public WorldObject< SeparatedConvolution<Q,NDIM> > {
    public:
        typedef Q opT;  ///< The apply function uses this to infer resultT=opT*inputT
    private:
        const int k;
        const int rank;
        const std::vector<long> vk;
        const std::vector<long> v2k;
        const std::vector<Slice> s0;
        std::vector<Q> signs;
        mutable std::vector< GaussianConvolution1D<Q> > ops;
        mutable SimpleCache< SeparatedConvolutionData<Q,NDIM> > data;
        mutable double nflop[64];

        /// Apply one of the separated terms, accumulating into the result
        template <typename T>
        void muopxv_fast(Level n,
                         const ConvolutionData1D<Q>* const ops[NDIM],
                         const Tensor<T>& f, const Tensor<T>& f0,
                         Tensor<TENSOR_RESULT_TYPE(T,Q)>& result,
                         Tensor<TENSOR_RESULT_TYPE(T,Q)>& result0,
                         const double tol,
                         const double musign, 
                         Tensor<TENSOR_RESULT_TYPE(T,Q)>& work1,
                         Tensor<TENSOR_RESULT_TYPE(T,Q)>& work2,
                         Tensor<TENSOR_RESULT_TYPE(T,Q)>& work3,
                         Tensor<TENSOR_RESULT_TYPE(T,Q)>& work4) const {

            long dimj = 2*k;
            long dimk = 2*k;
            long dimi = 1;
            for (int i=1; i<NDIM; i++) dimi *= 2*k;
            
            TENSOR_RESULT_TYPE(T,Q) *w1=work1.ptr(), *w2=work2.ptr();
            mTxmq(dimi, dimj, dimk, w1, f.ptr(), ops[0]->R.ptr());
            for (int d=1; d<NDIM; d++) {
                mTxmq(dimi, dimj, dimk, w2, w1, ops[d]->R.ptr());
                std::swap(w1,w2);
            }

            nflop[n] += NDIM*dimi*dimj*dimk*2.0;

            
            if (NDIM&1)
                result.gaxpy(1.0,work1,musign);
            else
                result.gaxpy(1.0,work2,musign);
            
            if (n > 0) {
                if (k&1) {
                    Tensor<TENSOR_RESULT_TYPE(T,Q)> tmp = inner(f0,ops[0]->T,0,0);
                    for (int d=1; d<NDIM; d++) {
                        tmp = inner(tmp,ops[d]->T,0,0);
                    }
                    result0.gaxpy(1.0,tmp,-musign);
                }
                else {
                    long dimj = k;
                    long dimk = k;
                    long dimi = 1;
                    for (int i=1; i<NDIM; i++) dimi *= k;

                    TENSOR_RESULT_TYPE(T,Q) *w3=work3.ptr(), *w4=work4.ptr();
                    mTxmq(dimi, dimj, dimk, w3, f0.ptr(), ops[0]->T.ptr());
                    for (int d=1; d<NDIM; d++) {
                        mTxmq(dimi, dimj, dimk, w4, w3, ops[d]->T.ptr());
                        std::swap(w3,w4);
                    }
                    
                    if (NDIM&1)
                        result0.gaxpy(1.0,work3,-musign);
                    else
                        result0.gaxpy(1.0,work4,-musign);
                    
                    nflop[n] += NDIM*dimi*dimj*dimk*2.0;
                }
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


        /// Computes the norm of one of the separated terms
        double munorm(Level n, const ConvolutionData1D<Q>* ops[]) const {
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
                if (iter>0 && ratio<1.2 && ratio>=1.0) break; // 1.2 was 1.02
                if (iter>10 && eval<tol) break;
                evalp = eval;
                if (iter == 99) throw "munorm failed";
            }
            return eval*ratio;
        }

        double munorm2(Level n, const ConvolutionData1D<Q>* ops[]) const {
            double prod=1.0, sum=0.0;
            for (int d=0; d<NDIM; d++) {
                Tensor<Q> q = copy(ops[d]->R);
                q(s0) = 0.0;
                double a = q.normf();
                double b = ops[d]->T.normf();
                double aa = std::min(a,b);
                double bb = std::max(a,b);
                prod *= bb;
                if (bb > 0.0) sum +=(aa/bb);
            }
            if (n) prod *= sum;

            return prod; 
        }


        const SeparatedConvolutionInternal<Q,NDIM> getmuop(int mu, Level n, const Displacement<NDIM>& disp) const {
            SeparatedConvolutionInternal<Q,NDIM> op;
            for (int d=0; d<NDIM; d++) {
                op.ops[d] = ops[mu].nonstandard(n, disp[d]);
            }
            double newnorm = munorm2(n, op.ops);
//             double oldnorm = munorm(n, op.ops);
//             if (newnorm < oldnorm) print("munorm", mu, newnorm, oldnorm);
            op.norm = newnorm;
            return op;
        }

        
        /// Returns pointer to cached operator
        const SeparatedConvolutionData<Q,NDIM>* getop(Level n, const Displacement<NDIM>& d) const {
            const SeparatedConvolutionData<Q,NDIM>* p = data.getptr(n,d);
            if (p) return p;

            SeparatedConvolutionData<Q,NDIM> op(rank);
            for (int mu=0; mu<rank; mu++) {
                op.muops[mu] = getmuop(mu, n, d);   // ???? SEGV HERE ????
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


    public:

        SeparatedConvolution(World& world,
                             long k, 
                             const Tensor<Q>& coeffs, 
                             const Tensor<double>& expnts)    // restriction to double must be relaxed
            : WorldObject< SeparatedConvolution<Q,NDIM> >(world)
            , k(k)
            , rank(coeffs.dim[0])
            , vk(NDIM,k)
            , v2k(NDIM,2*k)
            , s0(NDIM,Slice(0,k-1))
            , signs(rank) 
            , ops(rank)
        {
            // !!! NB ... cell volume obtained from global defaults
            Tensor<double> cell_width(FunctionDefaults<NDIM>::cell(_,1)-FunctionDefaults<NDIM>::cell(_,0));
            // Check that the cell is cubic since currently is assumed
            for (long d=1; d<NDIM; d++) {
                MADNESS_ASSERT(fabs(cell_width(d)-cell_width(0L)) < 1e-14*cell_width(0L));
            }

            for (int i=0; i<rank; i++) {
                Q coeff = coeffs(i);
                signs[i] = munge_sign(coeff); 
                coeff = std::pow(coeff,1.0/NDIM);
                // !!! Note that this assumes cubic box
                ops[i] = GaussianConvolution1D<Q>(k, 
                                                  coeff*cell_width(0L), 
                                                  expnts(i)*cell_width(0L)*cell_width(0L));
            }

            for (int i=0; i<64; i++) nflop[i] = 0.0;

            this->process_pending();
        }

        const double* get_nflop() const {
            return nflop;
        }

        double norm(Level n, const Displacement<NDIM>& d) const {
            return getop(n, d)->norm;
        }

        template <typename T>
        Tensor<TENSOR_RESULT_TYPE(T,Q)> apply(const Key<NDIM>& source,
                                              const Displacement<NDIM>& shift,
                                              const Tensor<T>& coeff,
                                              double tol) const {
            tol = tol/rank;
            typedef TENSOR_RESULT_TYPE(T,Q) resultT;
            const SeparatedConvolutionData<Q,NDIM>* op = getop(source.level(), shift);
            //print("sepop",source,shift,op->norm,tol);
            Tensor<resultT> r(v2k), r0(vk);
            Tensor<resultT> work1(v2k), work2(v2k), work3(vk), work4(vk);

            const Tensor<T>* input = &coeff;
            Tensor<T> dummy;

            if (coeff.dim[0] == k) {
                dummy = Tensor<T>(v2k);
                dummy(s0) = coeff;
                input = &dummy;
            }
            const Tensor<T> f0 = copy(coeff(s0));
            for (int mu=0; mu<rank; mu++) {
                const SeparatedConvolutionInternal<Q,NDIM>& muop =  op->muops[mu];
                if (muop.norm > tol) { 
                    muopxv_fast(source.level(), muop.ops, *input, f0, r, r0, tol, signs[mu], 
                                work1, work2, work3, work4);
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
        Tensor<double> cell_width(FunctionDefaults<NDIM>::cell(_,1)-FunctionDefaults<NDIM>::cell(_,0));
        double hi = sqrt(double(NDIM))*cell_width.normf(); // Diagonal width of cell
        const double pi = 3.14159265358979323846264338328;

        // bsh_fit generates representation for 1/4Pir but we want 1/r
        // so have to scale eps by 1/4Pi
        Tensor<double> coeff, expnt;
        bsh_fit(0.0, lo, hi, eps/(4.0*pi), &coeff, &expnt, false);
        coeff.scale(4.0*pi);
        return SeparatedConvolution<Q,NDIM>(world, k, coeff, expnt);
    }

    /// Factory function generating separated kernel for convolution with exp(-mu*r)/(4*pi*r)
    template <typename Q, int NDIM>
    SeparatedConvolution<Q,NDIM> BSHOperator(World& world,
                                             double mu,
                                             long k,
                                             double lo,
                                             double eps) {
        Tensor<double> cell_width(FunctionDefaults<NDIM>::cell(_,1)-FunctionDefaults<NDIM>::cell(_,0));
        double hi = sqrt(double(NDIM))*cell_width.normf(); // Diagonal width of cell
        Tensor<double> coeff, expnt;
        bsh_fit(mu, lo, hi, eps, &coeff, &expnt, false);
        return SeparatedConvolution<Q,NDIM>(world, k, coeff, expnt);
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
