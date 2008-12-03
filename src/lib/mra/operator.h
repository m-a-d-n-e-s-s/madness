#ifndef MAD_OPERATOR_H
#define MAD_OPERATOR_H

#include <mra/mra.h>
#include <limits.h>
#include <mra/adquad.h>
#include <tensor/mtxmq.h>
#include <tensor/aligned.h>
#include <linalg/tensor_lapack.h>

#include <mra/simplecache.h>
#include <mra/convolution1d.h>
#include <mra/displacements.h>

namespace madness {


    extern void bsh_fit(double mu, double lo, double hi, double eps,
                        Tensor<double> *pcoeff, Tensor<double> *pexpnt, bool prnt=false);

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
        bool isperiodicsum;///< If true the operator 1D kernels have been summed over lattice translations and may be non-zero at both ends of the unit cell
    private:
        const int k;
        const int rank;
        const std::vector<long> vk;
        const std::vector<long> v2k;
        const std::vector<Slice> s0;
        mutable std::vector< SharedPtr< Convolution1D<Q> > > ops;

        mutable SimpleCache< SeparatedConvolutionData<Q,NDIM>, NDIM > data;

        struct Transformation {
            long r;             // Effective rank of transformation
            const Q* U;         // Ptr to matrix
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

            R* restrict w1=work1.ptr();
            R* restrict w2=work2.ptr();
            Q* restrict w3=work3.ptr();

            const Q* U;

            U = (trans[0].r == dimk) ? trans[0].U : shrink(dimk,dimk,trans[0].r,trans[0].U,w3);
            mTxmq(dimi, trans[0].r, dimk, w1, f.ptr(), U);
            size = trans[0].r * size / dimk;
            dimi = size/dimk;
            for (int d=1; d<NDIM; d++) {
                U = (trans[d].r == dimk) ? trans[d].U : shrink(dimk,dimk,trans[d].r,trans[d].U,w3);
                mTxmq(dimi, trans[d].r, dimk, w2, w1, U);
                size = trans[d].r * size / dimk;
                dimi = size/dimk;
                std::swap(w1,w2);
            }

            // If all blocks are full rank we can skip the transposes
            bool doit = false;
            for (int d=0; d<NDIM; d++) doit = doit || trans[d].VT;

            if (doit) {
                for (int d=0; d<NDIM; d++) {
                    if (trans[d].VT) {
                        dimi = size/trans[d].r;
                        mTxmq(dimi, dimk, trans[d].r, w2, w1, trans[d].VT);
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
            R* restrict p = result.ptr();
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
                         Tensor<TENSOR_RESULT_TYPE(T,Q)>& work2,
                         Tensor<Q>& work5) const {

            PROFILE_MEMBER_FUNC(SeparatedConvolution);
            Transformation trans[NDIM];

            double Rnorm = 1.0;
            for (int d=0; d<NDIM; d++) Rnorm *= ops[d]->Rnorm;
            if (Rnorm == 0.0) return;

            tol = tol/(Rnorm*NDIM);  // Errors are relative within here

            // Determine rank of SVD to use or if to use the full matrix
            const long twok = 2*k;
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

            double Rnorm = 1.0;
            for (int d=0; d<NDIM; d++) Rnorm *= ops[d]->Rnorm;
            if (Rnorm == 0.0) return;

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

            double Rnorm = 1.0;
            for (int d=0; d<NDIM; d++) Rnorm *= ops[d]->Rnorm;
            if (Rnorm == 0.0) return;

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
                             std::vector< SharedPtr< Convolution1D<Q> > >& ops,
                             bool doleaves = false,
                             bool isperiodicsum = false)
            : WorldObject< SeparatedConvolution<Q,NDIM> >(world)
            , doleaves(doleaves)
            , isperiodicsum(isperiodicsum)
            , k(k)
            , rank(ops.size())
            , vk(NDIM,k)
            , v2k(NDIM,2*k)
            , s0(std::max(2,NDIM),Slice(0,k-1))
            , ops(ops)
        {

            check_cubic();

            this->process_pending();
        }

        /// Constructor for Gaussian Convolutions (mostly for backward compatability)
        SeparatedConvolution(World& world,
                             int k,
                             const Tensor<Q>& coeff, const Tensor<double>& expnt,
                             bool doleaves = false)
            : WorldObject< SeparatedConvolution<Q,NDIM> >(world)
            , doleaves(doleaves)
            , isperiodicsum(false)
            , k(k)
            , rank(coeff.dim[0])
            , vk(NDIM,k)
            , v2k(NDIM,2*k)
            , s0(std::max(2,NDIM),Slice(0,k-1))
            , ops(coeff.dim[0])
        {
            check_cubic();
            double width = FunctionDefaults<NDIM>::get_cell_width()(0L);

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

        const std::vector< Key<NDIM> > get_disp(Level n) const {
            return Displacements<NDIM>().get_disp(n, isperiodicsum);
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
            Tensor<Q> work5(2*k,2*k);

            const Tensor<T> f0 = copy(coeff(s0));
            for (int mu=0; mu<rank; mu++) {
                const SeparatedConvolutionInternal<Q,NDIM>& muop =  op->muops[mu];
                //print(source, shift, mu, muop.norm);
                if (muop.norm > tol) {
                    muopxv_fast(source.level(), muop.ops, *input, f0, r, r0, tol, ops[mu]->sign,
                                work1, work2, work5);
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
