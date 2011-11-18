/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id$
*/
#ifndef MADNESS_MRA_OPERATOR_H__INCLUDED
#define MADNESS_MRA_OPERATOR_H__INCLUDED

/// \file mra/operator.h
/// \brief Implements most functionality of separated operators

/// \ingroup function

//extern "C" void daxpy_(const long*, const double*, const double*, const long*, double*, const long*);

#include <mra/mra.h>
#include <limits.h>
#include <mra/adquad.h>
#include <tensor/mtxmq.h>
#include <tensor/aligned.h>
#include <linalg/tensor_lapack.h>
#include <constants.h>

#include <mra/simplecache.h>
#include <mra/convolution1d.h>
#include <mra/displacements.h>

namespace madness {


    extern void truncate_periodic_expansion(Tensor<double>& c, Tensor<double>& e,
      double L, bool discardG0);

    extern void bsh_fit(double mu, double lo, double hi, double eps,
                            Tensor<double> *pcoeff, Tensor<double> *pexpnt, bool prnt=false);

    extern void bsh_fit_ndim(int ndim, double mu, double lo, double hi, double eps,
                                 Tensor<double> *pcoeff, Tensor<double> *pexpnt, bool prnt=false);

    template <typename Q, std::size_t NDIM>
    struct SeparatedConvolutionInternal {
        double norm;
        const ConvolutionData1D<Q>* ops[NDIM];
    };

    template <typename Q, std::size_t NDIM>
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
    template <typename Q, std::size_t NDIM>
    class SeparatedConvolution : public WorldObject< SeparatedConvolution<Q,NDIM> > {
    public:
        typedef Q opT;  ///< The apply function uses this to infer resultT=opT*inputT
        bool doleaves;  ///< If should be applied to leaf coefficients ... false by default
        bool isperiodicsum;///< If true the operator 1D kernels have been summed over lattice translations and may be non-zero at both ends of the unit cell
    private:
        mutable std::vector< ConvolutionND<Q,NDIM> > ops;
        const BoundaryConditions<NDIM> bc;
        const int k;
        const int rank;
        const std::vector<long> vk;
        const std::vector<long> v2k;
        const std::vector<Slice> s0;

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
                                  const Q mufac,
                                  Tensor<R>& result) const {

            PROFILE_MEMBER_FUNC(SeparatedConvolution);
            long size = 1;
            for (std::size_t i=0; i<NDIM; ++i) size *= dimk;
            long dimi = size/dimk;

            R* restrict w1=work1.ptr();
            R* restrict w2=work2.ptr();
            Q* restrict w3=work3.ptr();

            const Q* U;

            U = (trans[0].r == dimk) ? trans[0].U : shrink(dimk,dimk,trans[0].r,trans[0].U,w3);
            mTxmq(dimi, trans[0].r, dimk, w1, f.ptr(), U);
            size = trans[0].r * size / dimk;
            dimi = size/dimk;
            for (std::size_t d=1; d<NDIM; ++d) {
                U = (trans[d].r == dimk) ? trans[d].U : shrink(dimk,dimk,trans[d].r,trans[d].U,w3);
                mTxmq(dimi, trans[d].r, dimk, w2, w1, U);
                size = trans[d].r * size / dimk;
                dimi = size/dimk;
                std::swap(w1,w2);
            }

            // If all blocks are full rank we can skip the transposes
            bool doit = false;
            for (std::size_t d=0; d<NDIM; ++d) doit = doit || trans[d].VT;

            if (doit) {
                for (std::size_t d=0; d<NDIM; ++d) {
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
            // Assuming here that result is contiguous and aligned
            aligned_axpy(size, result.ptr(), w1, mufac);
            //    long one = 1;
            //daxpy_(&size, &mufac, w1, &one, result.ptr(), &one);
        }

        /// Apply one of the separated terms, accumulating into the result
        template <typename T>
        void muopxv_fast(Level n,
                         const ConvolutionData1D<Q>* const ops[NDIM],
                         const Tensor<T>& f, const Tensor<T>& f0,
                         Tensor<TENSOR_RESULT_TYPE(T,Q)>& result,
                         Tensor<TENSOR_RESULT_TYPE(T,Q)>& result0,
                         double tol,
                         const Q mufac,
                         Tensor<TENSOR_RESULT_TYPE(T,Q)>& work1,
                         Tensor<TENSOR_RESULT_TYPE(T,Q)>& work2,
                         Tensor<Q>& work5) const {

            PROFILE_MEMBER_FUNC(SeparatedConvolution);
            Transformation trans[NDIM];

            double Rnorm = 1.0;
            for (std::size_t d=0; d<NDIM; ++d) Rnorm *= ops[d]->Rnorm;
            if (Rnorm == 0.0) return;

            tol = tol/(Rnorm*NDIM);  // Errors are relative within here

            // Determine rank of SVD to use or if to use the full matrix
            const long twok = 2*k;
            long break_even;
            if (NDIM==1) break_even = long(0.5*twok);
            else if (NDIM==2) break_even = long(0.6*twok);
            else if (NDIM==3) break_even=long(0.65*twok);
            else break_even=long(0.7*twok);
            for (std::size_t d=0; d<NDIM; ++d) {
                long r;
                for (r=0; r<twok; ++r) {
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
            apply_transformation(n, twok, trans, f, work1, work2, work5, mufac, result);

            if (n > 0) {
                if (NDIM==1) break_even = long(0.5*k);
                else if (NDIM==2) break_even = long(0.6*k);
                else if (NDIM==3) break_even=long(0.65*k);
                else break_even=long(0.7*k);
                for (std::size_t d=0; d<NDIM; ++d) {
                    long r;
                    for (r=0; r<k; ++r) {
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
                apply_transformation(n, k, trans, f0, work1, work2, work5, -mufac, result0);
            }
        }


        /// Computes the Frobenius norm of one of the separated terms ... WITHOUT FACTOR INCLUDED
        double munorm2(Level n, const ConvolutionData1D<Q>* ops[]) const {
            PROFILE_MEMBER_FUNC(SeparatedConvolution);
            double prodR=1.0, prodT=1.0;
            for (std::size_t d=0; d<NDIM; ++d) {
                prodR *= ops[d]->Rnormf;
                prodT *= ops[d]->Tnormf;

            }
            if (n) prodR = sqrt(std::max(prodR*prodR - prodT*prodT,0.0));

            if (prodR < 1e-8*prodT) {
                double prod=1.0, sum=0.0;
                for (std::size_t d=0; d<NDIM; ++d) {
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
            for (std::size_t d=0; d<NDIM; ++d) {
                op.ops[d] = ops[mu].getop(d)->nonstandard(n, disp.translation()[d]);
            }
            op.norm = munorm2(n, op.ops)*std::abs(ops[mu].getfac());

//             double newnorm = munorm2(n, op.ops);
//             // This rescaling empirically based upon BSH separated expansion
//             // ... needs more testing.  OK also for TDSE.
//             // All is good except for some 000 blocks which are up to sqrt(k^d) off.
//             for (int d=0; d<NDIM; ++d)  {
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
            for (int mu=0; mu<rank; ++mu) {
                op.muops[mu] = getmuop(mu, n, d);
            }

            double norm = 0.0;
            for (int mu=0; mu<rank; ++mu) {
                const double munorm = op.muops[mu].norm;
                norm += munorm*munorm;
            }
	    //print("getop", n, d, norm);
            op.norm = sqrt(norm);
            data.set(n, d, op);
            return data.getptr(n,d);
        }

        void check_cubic() {
            // !!! NB ... cell volume obtained from global defaults
            const Tensor<double>& cell_width = FunctionDefaults<NDIM>::get_cell_width();
            // Check that the cell is cubic since currently is assumed
            for (std::size_t d=1; d<NDIM; ++d) {
                MADNESS_ASSERT(fabs(cell_width(d)-cell_width(0L)) < 1e-14*cell_width(0L));
            }
        }

    public:
        // Default constructor for invoking compute member functions on aggregate arguments
        SeparatedConvolution(World * w) :WorldObject<SeparatedConvolution<Q,NDIM> >(*w), k(0), rank(0) {
        }

        // For separated convolutions with same operator in each direction (isotropic)
        SeparatedConvolution(World& world,
                             std::vector< std::shared_ptr< Convolution1D<Q> > >& argops,
                             const BoundaryConditions<NDIM>& bc = FunctionDefaults<NDIM>::get_bc(),
                             long k = FunctionDefaults<NDIM>::get_k(),
                             bool doleaves = false)
                : WorldObject< SeparatedConvolution<Q,NDIM> >(world)
                , doleaves(doleaves)
                , isperiodicsum(bc(0,0)==BC_PERIODIC)
                , bc(bc)
                , k(k)
                , rank(argops.size())
                , vk(NDIM,k)
                , v2k(NDIM,2*k)
                , s0(std::max<std::size_t>(2,NDIM),Slice(0,k-1))
        {
            // Presently we must have periodic or non-periodic in all dimensions.
            for (std::size_t d=1; d<NDIM; ++d) {
                MADNESS_ASSERT(bc(d,0)==bc(0,0));
            }
            check_cubic();

            for (unsigned int mu=0; mu < argops.size(); ++mu) {
              this->ops.push_back(ConvolutionND<Q,NDIM>(argops[mu]));
            }

            this->process_pending();
        }

        //copy constructor
        /*SeparatedConvolution(const SeparatedConvolution& sc)
                : WorldObject< SeparatedConvolution<Q,NDIM> >(sc.world)
                , doleaves(sc.doleaves)
                , isperiodicsum(sc.isperiodicsum)
                , ops(sc.ops)
                , bc(sc.bc)
                , k(sc.k)
                , rank(sc.rank)
                , vk(sc.vk)
                , v2k(sc.v2k)
                , s0(sc.s0)
        {
            this->process_pending();
        }*/
        // For general convolutions
        SeparatedConvolution(World& world,
                             std::vector< ConvolutionND<Q,NDIM> >& argops,
                             const BoundaryConditions<NDIM>& bc = FunctionDefaults<NDIM>::get_bc(),
                             long k = FunctionDefaults<NDIM>::get_k(),
                             bool doleaves = false)
                : WorldObject< SeparatedConvolution<Q,NDIM> >(world)
                , doleaves(doleaves)
                , isperiodicsum(bc(0,0)==BC_PERIODIC)
                , ops(argops)
                , bc(bc)
                , k(k)
                , rank(argops.size())
                , vk(NDIM,k)
                , v2k(NDIM,2*k)
                , s0(std::max<std::size_t>(2,NDIM),Slice(0,k-1))
        {
            // Presently we must have periodic or non-periodic in all dimensions.
            for (std::size_t d=1; d<NDIM; ++d) {
                MADNESS_ASSERT(bc(d,0)==bc(0,0));
            }
            this->process_pending();
        }

        /// Constructor for Gaussian Convolutions (mostly for backward compatability)
        SeparatedConvolution(World& world,
                             const Tensor<Q>& coeff, const Tensor<double>& expnt,
                             const BoundaryConditions<NDIM>& bc = FunctionDefaults<NDIM>::get_bc(),
                             int k=FunctionDefaults<NDIM>::get_k(),
                             bool doleaves = false)
                : WorldObject< SeparatedConvolution<Q,NDIM> >(world)
                , doleaves(doleaves)
                , isperiodicsum(bc(0,0)==BC_PERIODIC)
                , ops(coeff.dim(0))
                , bc(bc)
                , k(k)
                , rank(coeff.dim(0))
                , vk(NDIM,k)
                , v2k(NDIM,2*k)
                , s0(std::max<std::size_t>(2,NDIM),Slice(0,k-1))
        {
            // Presently we must have periodic or non-periodic in all dimensions.
            for (std::size_t d=1; d<NDIM; ++d) {
                MADNESS_ASSERT(bc(d,0)==bc(0,0));
            }

            const Tensor<double>& width = FunctionDefaults<NDIM>::get_cell_width();
            const double pi = constants::pi;

            for (int mu=0; mu<rank; ++mu) {
                Q c = std::pow(sqrt(expnt(mu)/pi),static_cast<int>(NDIM)); // Normalization coeff

                // We cache the normalized operator so the factor is the value we must multiply
                // by to recover the coeff we want.
                ops[mu].setfac(coeff(mu)/c);

                for (std::size_t d=0; d<NDIM; ++d) {
                  ops[mu].setop(d,GaussianConvolution1DCache<Q>::get(k, expnt(mu)*width[d]*width[d], 0, isperiodicsum));
                }
            }
        }

        /// WSTHORNTON Constructor for Gaussian Convolutions (mostly for backward compatability)
        SeparatedConvolution(World& world,
                             Vector<double,NDIM> args,
                             const Tensor<Q>& coeff, const Tensor<double>& expnt,
                             const BoundaryConditions<NDIM>& bc = FunctionDefaults<NDIM>::get_bc(),
                             int k=FunctionDefaults<NDIM>::get_k(),
                             bool doleaves=false)
                : WorldObject< SeparatedConvolution<Q,NDIM> >(world)
                , doleaves(doleaves)
                , isperiodicsum(bc(0,0)==BC_PERIODIC)
                , ops(coeff.dim(0))
                , bc(bc)
                , k(k)
                , rank(coeff.dim(0))
                , vk(NDIM,k)
                , v2k(NDIM,2*k)
                , s0(std::max<std::size_t>(2,NDIM),Slice(0,k-1))
        {
            // Presently we must have periodic or non-periodic in all dimensions.
            for (std::size_t d=1; d<NDIM; ++d) {
                MADNESS_ASSERT(bc(d,0)==bc(0,0));
            }

            const Tensor<double>& width = FunctionDefaults<NDIM>::get_cell_width();

            for (int mu=0; mu<rank; ++mu) {
                for (std::size_t d=0; d<NDIM; ++d) {
                  Q c = pow(coeff[mu],1.0/NDIM);
                  std::shared_ptr<GaussianConvolution1D<Q> >
                      gcptr(new GaussianConvolution1D<Q>(k, c*width[d], expnt(mu)*width[d]*width[d],
                              0, isperiodicsum, args[d]));
                  ops[mu].setop(d,gcptr);
                }
            }
        }

        virtual ~SeparatedConvolution() { }

        const BoundaryConditions<NDIM>& get_bc() const {return bc;}

        const std::vector< Key<NDIM> > get_disp(Level n) const {
            return Displacements<NDIM>().get_disp(n, isperiodicsum);
        }

        double norm(Level n, const Key<NDIM>& d) const {
            return getop(n, d)->norm;
        }

        template <typename T>
        Function<TENSOR_RESULT_TYPE(T,Q),NDIM> operator()(const Function<T,NDIM>& f) const {
            return madness::apply(*this, f);
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

            if (coeff.dim(0) == k) {
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
                MADNESS_ASSERT(coeff.dim(0)==2*k);
            }

            tol = tol/rank; // Error is per separated term

            const SeparatedConvolutionData<Q,NDIM>* op = getop(source.level(), shift);

            //print("sepop",source,shift,op->norm,tol);

            Tensor<resultT> r(v2k), r0(vk);
            Tensor<resultT> work1(v2k,false), work2(v2k,false);
            Tensor<Q> work5(2*k,2*k);

            const Tensor<T> f0 = copy(coeff(s0));
            for (int mu=0; mu<rank; ++mu) {
                const SeparatedConvolutionInternal<Q,NDIM>& muop =  op->muops[mu];
                //print("muop",source, shift, mu, muop.norm);
                if (muop.norm > tol) {
                    Q fac = ops[mu].getfac();
                    muopxv_fast(source.level(), muop.ops, *input, f0, r, r0, tol/std::abs(fac), fac,
                                work1, work2, work5);
                }
            }
            r(s0).gaxpy(1.0,r0,1.0);
            return r;
        }

        typedef Key<NDIM> keyT;


        template <typename T>
        std::tr1::tuple< Tensor<TENSOR_RESULT_TYPE(T,Q)> *, Tensor<TENSOR_RESULT_TYPE(T,Q)> *,
                         WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> >, keyT, double, double> 
        apply_computept(std::tr1::tuple<keyT, keyT, keyT, 
                                      double, double, double, 
                                      Tensor<TENSOR_RESULT_TYPE(T,Q)>, 
                                      WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > >& t1) const{
            
          typedef TENSOR_RESULT_TYPE(T,Q) resultT;
          typedef resultT R;
          typedef  WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > dcT;

          /*const*/ keyT/*&*/ argskey = std::tr1::get<0>(t1);
          /*const*/ keyT/*&*/ argsd = std::tr1::get<1>(t1);
          keyT/*&*/ argsdest = std::tr1::get<2>(t1);
          double/*&*/ argstol = std::tr1::get<3>(t1);
          double/*&*/ argsfac = std::tr1::get<4>(t1);
          double/*&*/ argscnorm = std::tr1::get<5>(t1);
          /*const*/ Tensor<R>/*&*/ coeff = std::tr1::get<6>(t1);
          dcT/*&*/ coeffs = std::tr1::get<7>(t1);

          /*const*/ Key<NDIM>/*&*/ source = argskey;
          /*const*/ Key<NDIM>/*&*/ shift = argsd;
          double tol = argstol/argsfac/argscnorm;

            PROFILE_MEMBER_FUNC(SeparatedConvolution);
            //print("inlined_apply \n");
            const Tensor<T>* input = &coeff;
            Tensor<T> dummy;

            if (coeff.dim(0) == k) {
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
                MADNESS_ASSERT(coeff.dim(0)==2*k);
            }

            tol = tol/rank; // Error is per separated term

            const SeparatedConvolutionData<Q,NDIM>* op = getop(source.level(), shift);

            //print("sepop",source,shift,op->norm,tol);

            Tensor<resultT> r(v2k);
            Tensor<resultT> r0(vk);
            Tensor<resultT> work1(v2k,false), work2(v2k,false);
            Tensor<Q> work5(2*k,2*k);


            //std::vector<Slice> s1(std::max<std::size_t>(2,NDIM),Slice(0,k-1))

            const Tensor<T> f0 = copy(coeff(s0));
            
            Level n = source.level();
            const Tensor<T>& f = *input;
            Transformation trans[NDIM];
	    const long twok = 2*k;
	    long break_even;
	    
            if (NDIM==1) break_even = long(0.5*twok);
	    else if (NDIM==2) break_even = long(0.6*twok);
	    else if (NDIM==3) break_even=long(0.65*twok);
	    else break_even=long(0.7*twok);
	    
            long break_even2;
            if (NDIM==1) break_even2 = long(0.5*k);
	    else if (NDIM==2) break_even2 = long(0.6*k);
	    else if (NDIM==3) break_even2=long(0.65*k);
	    else break_even2=long(0.7*k);

            R* restrict w1=work1.ptr();
            R* restrict w2=work2.ptr();
	    Q* restrict w3=work5.ptr();

            for (int mu=0; mu<rank; ++mu) {
                const SeparatedConvolutionInternal<Q,NDIM>& muop =  op->muops[mu];
                //print("muop",source, shift, mu, muop.norm);
                if (muop.norm > tol) {
                    Q fac = ops[mu].getfac();
                    //muopxv_fast(source.level(), muop.ops, *input, f0, r, r0, tol/std::abs(fac), fac,
                    //            work1, work2, work5);

                    //glue
                    const ConvolutionData1D<Q>* const* ops/*[NDIM]*/ = muop.ops;
                    //const Tensor<T>& f0
                    Tensor<TENSOR_RESULT_TYPE(T,Q)>& result = r;
                    Tensor<TENSOR_RESULT_TYPE(T,Q)>& result0 = r0;
                    double tol1 = tol/std::abs(fac);
                    const Q mufac = fac;
                    //Tensor<TENSOR_RESULT_TYPE(T,Q)>& work1
                    //Tensor<TENSOR_RESULT_TYPE(T,Q)>& work2
                    //Tensor<Q>& work5
                     
		    double Rnorm = 1.0;
		    for (std::size_t d=0; d<NDIM; ++d) Rnorm *= ops[d]->Rnorm;
		    if (Rnorm == 0.0) continue;

		    tol1 = tol1/(Rnorm*NDIM);  // Errors are relative within here

		    // Determine rank of SVD to use or if to use the full matrix
		    for (std::size_t d=0; d<NDIM; ++d) {
			long r1;
			for (r1=0; r1<twok; ++r1) {
			    if (ops[d]->Rs[r1] < tol1) break;
			}
			if (r1 >= break_even) {
			    trans[d].r = twok;
			    trans[d].U = ops[d]->R.ptr();
			    trans[d].VT = 0;
			}
			else {
			    r1 += (r1&1L);
			    trans[d].r = std::max(2L,r1);
			    trans[d].U = ops[d]->RU.ptr();
			    trans[d].VT = ops[d]->RVT.ptr();
			}
		    }
		    ////apply_transformation(n, twok, trans, f, work1, work2, work5, mufac, result);

		    long dimk = twok;
		   
                    long size = 1;
		    for (std::size_t i=0; i<NDIM; ++i) size *= dimk;
		    long dimi = size/dimk;

		    const Q* U;

		    U = (trans[0].r == dimk) ? trans[0].U : shrink(dimk,dimk,trans[0].r,trans[0].U,w3);
                    ////GPU
		    mTxmq(dimi, trans[0].r, dimk, w1, f.ptr(), U);
		    size = trans[0].r * size / dimk;
		    dimi = size/dimk;
		    for (std::size_t d=1; d<NDIM; ++d) {
			U = (trans[d].r == dimk) ? trans[d].U : shrink(dimk,dimk,trans[d].r,trans[d].U,w3);
                        ////GPU
			mTxmq(dimi, trans[d].r, dimk, w2, w1, U);
			size = trans[d].r * size / dimk;
			dimi = size/dimk;
                        ////GPU
			std::swap(w1,w2);
		    }

		    // If all blocks are full rank we can skip the transposes
		    bool doit = false;
		    for (std::size_t d=0; d<NDIM; ++d) doit = doit || trans[d].VT;

		    if (doit) {
			for (std::size_t d=0; d<NDIM; ++d) {
			    if (trans[d].VT) {
				dimi = size/trans[d].r;
                                ////GPU
				mTxmq(dimi, dimk, trans[d].r, w2, w1, trans[d].VT);
				size = dimk*size/trans[d].r;
			    }
			    else {
                                ////GPU
				fast_transpose(dimk, dimi, w1, w2);
			    }
                            ////GPU
			    std::swap(w1,w2);
			}
		    }
		    // Assuming here that result is contiguous and aligned
                    ////GPU
		    aligned_axpy(size, result.ptr(), w1, mufac);
		    //    long one = 1;
		    //daxpy_(&size, &mufac, w1, &one, result.ptr(), &one);

		    if (n > 0) {
			for (std::size_t d=0; d<NDIM; ++d) {
			    long r1;
			    for (r1=0; r1<k; ++r1) {
				if (ops[d]->Ts[r1] < tol1) break;
			    }
			    if (r1 >= break_even2) {
				trans[d].r = k;
				trans[d].U = ops[d]->T.ptr();
				trans[d].VT = 0;
			    }
			    else {
				r1 += (r1&1L);
				trans[d].r = std::max(2L,r1);
				trans[d].U = ops[d]->TU.ptr();
				trans[d].VT = ops[d]->TVT.ptr();
			    }
			}
			////apply_transformation(n, k, trans, f0, work1, work2, work5, -mufac, result0);
		        dimk = k;
                        const Tensor<T>& f1 = f0;
                        const Q mufac1 = -mufac;
                        Tensor<R>& result1 = result0;

			size = 1;
			for (std::size_t i=0; i<NDIM; ++i) size *= dimk;
			dimi = size/dimk;

			const Q* U1;

			U1 = (trans[0].r == dimk) ? trans[0].U : shrink(dimk,dimk,trans[0].r,trans[0].U,w3);
                        ////GPU
			mTxmq(dimi, trans[0].r, dimk, w1, f1.ptr(), U1);
			size = trans[0].r * size / dimk;
			dimi = size/dimk;
			for (std::size_t d=1; d<NDIM; ++d) {
	                    U1 = (trans[d].r == dimk) ? trans[d].U : shrink(dimk,dimk,trans[d].r,trans[d].U,w3);
                            ////GPU
			    mTxmq(dimi, trans[d].r, dimk, w2, w1, U1);
			    size = trans[d].r * size / dimk;
			    dimi = size/dimk;
                            ////GPU
		            std::swap(w1,w2);
			}

			// If all blocks are full rank we can skip the transposes
			bool doit = false;
			for (std::size_t d=0; d<NDIM; ++d) doit = doit || trans[d].VT;

			if (doit) {
			    for (std::size_t d=0; d<NDIM; ++d) {
				if (trans[d].VT) {
			            dimi = size/trans[d].r;
                                    ////GPU
				    mTxmq(dimi, dimk, trans[d].r, w2, w1, trans[d].VT);
				    size = dimk*size/trans[d].r;
				}
				else {
                                    ////GPU
			            fast_transpose(dimk, dimi, w1, w2);
				}
                                ////GPU
				std::swap(w1,w2);
		            }
			 }
			 // Assuming here that result is contiguous and aligned
                         ////GPU
			 aligned_axpy(size, result1.ptr(), w1, mufac1);
			 //    long one = 1;
			 //daxpy_(&size, &mufac, w1, &one, result.ptr(), &one);
                    }
                }
            }
            Tensor<R> * r1 = new Tensor<R>(r); 
            Tensor<R> * r01 = new Tensor<R>(r0); 
            std::tr1::tuple<Tensor<R>*, Tensor<R>*, dcT, keyT, double, double> t2(r1, r01, coeffs, argsdest, argstol, argsfac);
            return t2;  
        }

        //op->opt_inlined_apply(args.key, args.d, c, args.tol/args.fac/args.cnorm);
        template <typename T>
        std::tr1::tuple< Tensor<TENSOR_RESULT_TYPE(T,Q)> *, Tensor<TENSOR_RESULT_TYPE(T,Q)> *,
                         WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> >&, keyT&, double&, double&> 
        apply_compute(std::tr1::tuple<keyT&, keyT&, keyT&, 
                                      double&, double&, double&, 
                                      Tensor<TENSOR_RESULT_TYPE(T,Q)>&, 
                                      WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> >& >& t1) const{
            
          typedef TENSOR_RESULT_TYPE(T,Q) resultT;
          typedef resultT R;
          typedef  WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > dcT;

          const keyT& argskey = std::tr1::get<0>(t1);
          const keyT& argsd = std::tr1::get<1>(t1);
          keyT& argsdest = std::tr1::get<2>(t1);
          double& argstol = std::tr1::get<3>(t1);
          double& argsfac = std::tr1::get<4>(t1);
          double& argscnorm = std::tr1::get<5>(t1);
          const Tensor<R>& coeff = std::tr1::get<6>(t1);
          dcT& coeffs = std::tr1::get<7>(t1);

          const Key<NDIM>& source = argskey;
          const Key<NDIM>& shift = argsd;
          double tol = argstol/argsfac/argscnorm;

            PROFILE_MEMBER_FUNC(SeparatedConvolution);
            //print("inlined_apply \n");
            const Tensor<T>* input = &coeff;
            Tensor<T> dummy;

            if (coeff.dim(0) == k) {
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
                MADNESS_ASSERT(coeff.dim(0)==2*k);
            }

            tol = tol/rank; // Error is per separated term

            const SeparatedConvolutionData<Q,NDIM>* op = getop(source.level(), shift);

            //print("sepop",source,shift,op->norm,tol);

            Tensor<resultT> r(v2k);
            Tensor<resultT> r0(vk);
            Tensor<resultT> work1(v2k,false), work2(v2k,false);
            Tensor<Q> work5(2*k,2*k);


            const Tensor<T> f0 = copy(coeff(s0));
            
            Level n = source.level();
            const Tensor<T>& f = *input;
            Transformation trans[NDIM];
	    const long twok = 2*k;
	    long break_even;
	    
            if (NDIM==1) break_even = long(0.5*twok);
	    else if (NDIM==2) break_even = long(0.6*twok);
	    else if (NDIM==3) break_even=long(0.65*twok);
	    else break_even=long(0.7*twok);
	    
            long break_even2;
            if (NDIM==1) break_even2 = long(0.5*k);
	    else if (NDIM==2) break_even2 = long(0.6*k);
	    else if (NDIM==3) break_even2=long(0.65*k);
	    else break_even2=long(0.7*k);

            R* restrict w1=work1.ptr();
            R* restrict w2=work2.ptr();
	    Q* restrict w3=work5.ptr();

            for (int mu=0; mu<rank; ++mu) {
                const SeparatedConvolutionInternal<Q,NDIM>& muop =  op->muops[mu];
                //print("muop",source, shift, mu, muop.norm);
                if (muop.norm > tol) {
                    Q fac = ops[mu].getfac();
                    //muopxv_fast(source.level(), muop.ops, *input, f0, r, r0, tol/std::abs(fac), fac,
                    //            work1, work2, work5);

                    //glue
                    const ConvolutionData1D<Q>* const* ops/*[NDIM]*/ = muop.ops;
                    //const Tensor<T>& f0
                    Tensor<TENSOR_RESULT_TYPE(T,Q)>& result = r;
                    Tensor<TENSOR_RESULT_TYPE(T,Q)>& result0 = r0;
                    double tol1 = tol/std::abs(fac);
                    const Q mufac = fac;
                    //Tensor<TENSOR_RESULT_TYPE(T,Q)>& work1
                    //Tensor<TENSOR_RESULT_TYPE(T,Q)>& work2
                    //Tensor<Q>& work5
                     
		    double Rnorm = 1.0;
		    for (std::size_t d=0; d<NDIM; ++d) Rnorm *= ops[d]->Rnorm;
		    if (Rnorm == 0.0) continue;

		    tol1 = tol1/(Rnorm*NDIM);  // Errors are relative within here

		    // Determine rank of SVD to use or if to use the full matrix
		    for (std::size_t d=0; d<NDIM; ++d) {
			long r1;
			for (r1=0; r1<twok; ++r1) {
			    if (ops[d]->Rs[r1] < tol1) break;
			}
			if (r1 >= break_even) {
			    trans[d].r = twok;
			    trans[d].U = ops[d]->R.ptr();
			    trans[d].VT = 0;
			}
			else {
			    r1 += (r1&1L);
			    trans[d].r = std::max(2L,r1);
			    trans[d].U = ops[d]->RU.ptr();
			    trans[d].VT = ops[d]->RVT.ptr();
			}
		    }
		    ////apply_transformation(n, twok, trans, f, work1, work2, work5, mufac, result);

		    long dimk = twok;
		   
                    long size = 1;
		    for (std::size_t i=0; i<NDIM; ++i) size *= dimk;
		    long dimi = size/dimk;

		    const Q* U;

		    U = (trans[0].r == dimk) ? trans[0].U : shrink(dimk,dimk,trans[0].r,trans[0].U,w3);
                    ////GPU
		    mTxmq(dimi, trans[0].r, dimk, w1, f.ptr(), U);
		    size = trans[0].r * size / dimk;
		    dimi = size/dimk;
		    for (std::size_t d=1; d<NDIM; ++d) {
			U = (trans[d].r == dimk) ? trans[d].U : shrink(dimk,dimk,trans[d].r,trans[d].U,w3);
                        ////GPU
			mTxmq(dimi, trans[d].r, dimk, w2, w1, U);
			size = trans[d].r * size / dimk;
			dimi = size/dimk;
                        ////GPU
			std::swap(w1,w2);
		    }

		    // If all blocks are full rank we can skip the transposes
		    bool doit = false;
		    for (std::size_t d=0; d<NDIM; ++d) doit = doit || trans[d].VT;

		    if (doit) {
			for (std::size_t d=0; d<NDIM; ++d) {
			    if (trans[d].VT) {
				dimi = size/trans[d].r;
                                ////GPU
				mTxmq(dimi, dimk, trans[d].r, w2, w1, trans[d].VT);
				size = dimk*size/trans[d].r;
			    }
			    else {
                                ////GPU
				fast_transpose(dimk, dimi, w1, w2);
			    }
                            ////GPU
			    std::swap(w1,w2);
			}
		    }
		    // Assuming here that result is contiguous and aligned
                    ////GPU
		    aligned_axpy(size, result.ptr(), w1, mufac);
		    //    long one = 1;
		    //daxpy_(&size, &mufac, w1, &one, result.ptr(), &one);

		    if (n > 0) {
			for (std::size_t d=0; d<NDIM; ++d) {
			    long r1;
			    for (r1=0; r1<k; ++r1) {
				if (ops[d]->Ts[r1] < tol1) break;
			    }
			    if (r1 >= break_even2) {
				trans[d].r = k;
				trans[d].U = ops[d]->T.ptr();
				trans[d].VT = 0;
			    }
			    else {
				r1 += (r1&1L);
				trans[d].r = std::max(2L,r1);
				trans[d].U = ops[d]->TU.ptr();
				trans[d].VT = ops[d]->TVT.ptr();
			    }
			}
			////apply_transformation(n, k, trans, f0, work1, work2, work5, -mufac, result0);
		        dimk = k;
                        const Tensor<T>& f1 = f0;
                        const Q mufac1 = -mufac;
                        Tensor<R>& result1 = result0;

			size = 1;
			for (std::size_t i=0; i<NDIM; ++i) size *= dimk;
			dimi = size/dimk;

			const Q* U1;

			U1 = (trans[0].r == dimk) ? trans[0].U : shrink(dimk,dimk,trans[0].r,trans[0].U,w3);
                        ////GPU
			mTxmq(dimi, trans[0].r, dimk, w1, f1.ptr(), U1);
			size = trans[0].r * size / dimk;
			dimi = size/dimk;
			for (std::size_t d=1; d<NDIM; ++d) {
	                    U1 = (trans[d].r == dimk) ? trans[d].U : shrink(dimk,dimk,trans[d].r,trans[d].U,w3);
                            ////GPU
			    mTxmq(dimi, trans[d].r, dimk, w2, w1, U1);
			    size = trans[d].r * size / dimk;
			    dimi = size/dimk;
                            ////GPU
		            std::swap(w1,w2);
			}

			// If all blocks are full rank we can skip the transposes
			bool doit = false;
			for (std::size_t d=0; d<NDIM; ++d) doit = doit || trans[d].VT;

			if (doit) {
			    for (std::size_t d=0; d<NDIM; ++d) {
				if (trans[d].VT) {
			            dimi = size/trans[d].r;
                                    ////GPU
				    mTxmq(dimi, dimk, trans[d].r, w2, w1, trans[d].VT);
				    size = dimk*size/trans[d].r;
				}
				else {
                                    ////GPU
			            fast_transpose(dimk, dimi, w1, w2);
				}
                                ////GPU
				std::swap(w1,w2);
		            }
			 }
			 // Assuming here that result is contiguous and aligned
                         ////GPU
			 aligned_axpy(size, result1.ptr(), w1, mufac1);
			 //    long one = 1;
			 //daxpy_(&size, &mufac, w1, &one, result.ptr(), &one);
                    }
                }
            }
            Tensor<R> * r1 = new Tensor<R>(r); 
            Tensor<R> * r01 = new Tensor<R>(r0); 
            std::tr1::tuple<Tensor<R>*, Tensor<R>*, dcT&, keyT&, double&, double&> t2(r1, r01, coeffs, argsdest, argstol, argsfac);
            return t2;  
        }

/*
typedef std::tr1::tuple<const keyT&, const keyT&, const keyT&, const double&, const double&, const double&, const Tensor<R>&, dcT&> tuple1T;
            
tuple1T t1(args.key, args.d, args.dest, args.tol, args.fac, args.cnorm, c, coeffs);
            
typedef std::tr1::tuple< Tensor<R> *, Tensor<R> *,dcT&, const keyT&, const double&, const double&> tuple2T;

            typedef std::vector<tuple2T> (opT::*memfun2T)(std::vector< tuple1T& >, std::vector<opT*> );
            typedef Void(opT::*memfun3T)(tuple2T );

*/
	template<typename T, typename R>
	int foo() const{
		return 1;
	}	
	
        template<typename T, typename R>
        int foof1(std::vector<std::tr1::tuple<keyT&, keyT&, keyT&, double&, double&, double&, Tensor<R>&, WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> >& > > t1) const{
                return 1;
        }

        template<typename T, typename R>
        int foof2(std::vector<std::tr1::tuple<keyT&, keyT&, keyT&, double&, double&, double&, Tensor<R>&, WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> >& > > t1, int i) const{
                return 1;
        }

        template <typename T, typename opT>
        std::vector< std::tr1::tuple< Tensor<TENSOR_RESULT_TYPE(T,Q)> *, Tensor<TENSOR_RESULT_TYPE(T,Q)> *,
                         WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> >, keyT, double, double> >
        apply_allComputepreGPU(std::vector<std::tr1::tuple<keyT, keyT, keyT, 
                                      double, double, double, 
                                      Tensor<TENSOR_RESULT_TYPE(T,Q)>, 
                                      WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > > > inArgs, 
                      std::vector< SeparatedConvolution<Q,NDIM>* > inObj) const {

            print(inArgs.size());
            
            std::vector< std::tr1::tuple< Tensor<TENSOR_RESULT_TYPE(T,Q)> *, Tensor<TENSOR_RESULT_TYPE(T,Q)> *,
                         WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> >, keyT, double, double> > outArg/*(inArgs.size(), inObj.at(0)->apply_compute(inArgs.at(0)))*/;
           typedef TENSOR_RESULT_TYPE(T,Q) resultT;
	   typedef resultT R;
	   typedef  WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > dcT;
           

            for (unsigned i = 0; i < inArgs.size(); i++){
               // std::tr1::tuple<Tensor<TENSOR_RESULT_TYPE(T,Q)> *, Tensor<TENSOR_RESULT_TYPE(T,Q)> *,
               //          WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> >, keyT, double, double> temp = inObj.at(i)->apply_computept(inArgs.at(i));

                  std::tr1::tuple<keyT, keyT, keyT,
                                  double, double, double,
                                  Tensor<TENSOR_RESULT_TYPE(T,Q)>,
                                  WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > > arg1 = inArgs.at(i);

		  keyT& source = std::tr1::get<0>(arg1);
		  keyT& shift = std::tr1::get<1>(arg1);
		  keyT& argsdest = std::tr1::get<2>(arg1);
		  double argstol = std::tr1::get<3>(arg1);
		  double argsfac = std::tr1::get<4>(arg1);
		  double argscnorm = std::tr1::get<5>(arg1);
		  Tensor<R>& coeff = std::tr1::get<6>(arg1);
		  dcT& coeffs = std::tr1::get<7>(arg1);
                  int kref = inObj.at(i)->k;
                  int rankref = inObj.at(i)->rank;
                  const std::vector<long>& vkref = inObj.at(i)->vk;
                  const std::vector<long>& v2kref = inObj.at(i)->v2k;
                  const std::vector<Slice>& s0ref = inObj.at(i)->s0;
                  //const SimpleCache< SeparatedConvolutionData<Q,NDIM>, NDIM >& data = inObj.at(i)->data;

		  //Key<NDIM>& source = argskey;
		  //Key<NDIM>& shift = argsd;
		  double tol = argstol/argsfac/argscnorm;

		  const SeparatedConvolutionData<Q,NDIM>* op = inObj.at(i)->getop(source.level(), shift);

		  Tensor<resultT> r(v2kref);
		  Tensor<resultT> r0(vkref);
		  Tensor<resultT> work1(v2kref,false), work2(v2kref,false);
		  Tensor<Q> work5(2*kref,2*kref);

		  const Tensor<T>* input = &coeff;
		  Tensor<T> dummy;

		  if (coeff.dim(0) == kref) {
			// This processes leaf nodes with only scaling
			// coefficients ... FuncImpl::apply by default does not
			// apply the operator to these since for smoothing operators
			// it is not necessary.  It is necessary for operators such
			// as differentiation and time evolution and will also occur
			// if the application of the operator widens the tree.
			dummy = Tensor<T>(v2kref);
			dummy(s0ref) = coeff;
			input = &dummy;
		  }
		  else {
	              MADNESS_ASSERT(coeff.dim(0)==2*kref);
		  }

		  tol = tol/rankref; // Error is per separated term


		  //print("sepop",source,shift,op->norm,tol);


		  const Tensor<T> f0 = copy(coeff(s0ref));
		    
		  Level n = source.level();
		  const Tensor<T>& f = *input;
                  //trans can be made local
		  Transformation trans[NDIM];
		  const long twok = 2*kref;
		  long break_even;
		    
		  if (NDIM==1) break_even = long(0.5*twok);
		  else if (NDIM==2) break_even = long(0.6*twok);
		  else if (NDIM==3) break_even=long(0.65*twok);
		  else break_even=long(0.7*twok);
		    
		  long break_even2;
		  if (NDIM==1) break_even2 = long(0.5*kref);
		  else if (NDIM==2) break_even2 = long(0.6*kref);
		  else if (NDIM==3) break_even2=long(0.65*kref);
		  else break_even2=long(0.7*kref);

		  R* restrict w1=work1.ptr();
		  R* restrict w2=work2.ptr();
		  Q* restrict w3=work5.ptr();
			    
                  Tensor<TENSOR_RESULT_TYPE(T,Q)>& result = r;
	          Tensor<TENSOR_RESULT_TYPE(T,Q)>& result0 = r0;

                  R* resultptr = result.ptr();
                  R* result0ptr = result0.ptr();

		    for (int mu=0; mu<rankref; ++mu) {
			const SeparatedConvolutionInternal<Q,NDIM>& muop =  op->muops[mu];
			//print("muop",source, shift, mu, muop.norm);
			if (muop.norm > tol) {
			    Q fac = inObj.at(i)->ops[mu].getfac();
			    //muopxv_fast(source.level(), muop.ops, *input, f0, r, r0, tol/std::abs(fac), fac,
			    //            work1, work2, work5);

			    //glue
			    const ConvolutionData1D<Q>* const* ops/*[NDIM]*/ = muop.ops;
			    //const Tensor<T>& f0
			    //Tensor<TENSOR_RESULT_TYPE(T,Q)>& result = r;
			    //Tensor<TENSOR_RESULT_TYPE(T,Q)>& result0 = r0;
			    double tol1 = tol/std::abs(fac);
			    const Q mufac = fac;
			    //Tensor<TENSOR_RESULT_TYPE(T,Q)>& work1
			    //Tensor<TENSOR_RESULT_TYPE(T,Q)>& work2
			    //Tensor<Q>& work5
			     
			    double Rnorm = 1.0;
			    for (std::size_t d=0; d<NDIM; ++d) Rnorm *= ops[d]->Rnorm;
			    if (Rnorm == 0.0) continue;

			    tol1 = tol1/(Rnorm*NDIM);  // Errors are relative within here

			    // Determine rank of SVD to use or if to use the full matrix
			    for (std::size_t d=0; d<NDIM; ++d) {
				long r1;
				for (r1=0; r1<twok; ++r1) {
				    if (ops[d]->Rs[r1] < tol1) break;
				}
				if (r1 >= break_even) {
				    trans[d].r = twok;
				    trans[d].U = ops[d]->R.ptr();
				    trans[d].VT = 0;
				}
				else {
				    r1 += (r1&1L);
				    trans[d].r = std::max(2L,r1);
				    trans[d].U = ops[d]->RU.ptr();
				    trans[d].VT = ops[d]->RVT.ptr();
				}
			    }
			    ////apply_transformation(n, twok, trans, f, work1, work2, work5, mufac, result);

			    long dimk = twok;
			   
			    long size = 1;
			    for (std::size_t ii=0; ii<NDIM; ++ii) size *= dimk;
			    long dimi = size/dimk;

			    const Q* U;

			    U = (trans[0].r == dimk) ? trans[0].U : shrink(dimk,dimk,trans[0].r,trans[0].U,w3);
			    ////GPU
			    mTxmq(dimi, trans[0].r, dimk, w1, f.ptr(), U);
			    size = trans[0].r * size / dimk;
			    dimi = size/dimk;
			    for (std::size_t d=1; d<NDIM; ++d) {
				U = (trans[d].r == dimk) ? trans[d].U : shrink(dimk,dimk,trans[d].r,trans[d].U,w3);
				////GPU
				mTxmq(dimi, trans[d].r, dimk, w2, w1, U);
				size = trans[d].r * size / dimk;
				dimi = size/dimk;
				////GPU
			        std::swap(w1,w2);
			    }

			    // If all blocks are full rank we can skip the transposes
			    bool doit = false;
			    for (std::size_t d=0; d<NDIM; ++d) doit = doit || trans[d].VT;

			    if (doit) {
				for (std::size_t d=0; d<NDIM; ++d) {
				    if (trans[d].VT) {
					dimi = size/trans[d].r;
					////GPU
					mTxmq(dimi, dimk, trans[d].r, w2, w1, trans[d].VT);
					size = dimk*size/trans[d].r;
				    }
				    else {
					////GPU
					fast_transpose(dimk, dimi, w1, w2);
				    }
				    ////GPU
				    std::swap(w1,w2);
				}
			    }
			    // Assuming here that result is contiguous and aligned
			    ////GPU
			    aligned_axpy(size, /*result.ptr()*/resultptr, w1, mufac);
			    //    long one = 1;
			    //daxpy_(&size, &mufac, w1, &one, result.ptr(), &one);

			    if (n > 0) {
				for (std::size_t d=0; d<NDIM; ++d) {
				    long r1;
				    for (r1=0; r1<kref; ++r1) {
					if (ops[d]->Ts[r1] < tol1) break;
				    }
				    if (r1 >= break_even2) {
					trans[d].r = kref;
					trans[d].U = ops[d]->T.ptr();
					trans[d].VT = 0;
				    }
				    else {
					r1 += (r1&1L);
					trans[d].r = std::max(2L,r1);
					trans[d].U = ops[d]->TU.ptr();
					trans[d].VT = ops[d]->TVT.ptr();
				    }
				}
				////apply_transformation(n, k, trans, f0, work1, work2, work5, -mufac, result0);
				dimk = kref;
				//const Tensor<T>& f1 = f0;
				const Q mufac1 = -mufac;
				//Tensor<R>& result1 = result0;

				size = 1;
				for (std::size_t i=0; i<NDIM; ++i) size *= dimk;
				dimi = size/dimk;

				const Q* U1;

				U1 = (trans[0].r == dimk) ? trans[0].U : shrink(dimk,dimk,trans[0].r,trans[0].U,w3);
				////GPU
				mTxmq(dimi, trans[0].r, dimk, w1, f0.ptr(), U1);
				size = trans[0].r * size / dimk;
				dimi = size/dimk;
				for (std::size_t d=1; d<NDIM; ++d) {
				    U1 = (trans[d].r == dimk) ? trans[d].U : shrink(dimk,dimk,trans[d].r,trans[d].U,w3);
				    ////GPU
				    mTxmq(dimi, trans[d].r, dimk, w2, w1, U1);
				    size = trans[d].r * size / dimk;
				    dimi = size/dimk;
				    ////GPU
				    std::swap(w1,w2);
				}

				// If all blocks are full rank we can skip the transposes
				bool doit = false;
				for (std::size_t d=0; d<NDIM; ++d) doit = doit || trans[d].VT;

				if (doit) {
				    for (std::size_t d=0; d<NDIM; ++d) {
					if (trans[d].VT) {
					    dimi = size/trans[d].r;
					    ////GPU
					    mTxmq(dimi, dimk, trans[d].r, w2, w1, trans[d].VT);
					    size = dimk*size/trans[d].r;
					}
					else {
					    ////GPU
					    fast_transpose(dimk, dimi, w1, w2);
					}
					////GPU
					std::swap(w1,w2);
				    }
				 }
				 // Assuming here that result is contiguous and aligned
				 ////GPU
				 aligned_axpy(size, /*result0.ptr()*/result0ptr, w1, mufac1);
				 //    long one = 1;
				 //daxpy_(&size, &mufac, w1, &one, result.ptr(), &one);
			    }
			}
		    }
		    Tensor<R> * r1 = new Tensor<R>(r); 
		    Tensor<R> * r01 = new Tensor<R>(r0); 
		    std::tr1::tuple<Tensor<R>*, Tensor<R>*, dcT, keyT, double, double> t2(r1, r01, coeffs, argsdest, argstol, argsfac);
                    outArg.push_back(t2);
            }
 
            return outArg;
        }

        template <typename T, typename opT>
        std::vector< std::tr1::tuple< Tensor<TENSOR_RESULT_TYPE(T,Q)> *, Tensor<TENSOR_RESULT_TYPE(T,Q)> *,
                         WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> >, keyT, double, double> >
        apply_allComputeGPU(std::vector<std::tr1::tuple<keyT, keyT, keyT, 
                                      double, double, double, 
                                      Tensor<TENSOR_RESULT_TYPE(T,Q)>, 
                                      WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > > > inArgs, 
                      std::vector< SeparatedConvolution<Q,NDIM>* > inObj) const {

            print("      apply_allComputeGPU              ",inArgs.size());
            
            std::vector< std::tr1::tuple< Tensor<TENSOR_RESULT_TYPE(T,Q)> *, Tensor<TENSOR_RESULT_TYPE(T,Q)> *,
                         WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> >, keyT, double, double> > outArg/*(inArgs.size(), inObj.at(0)->apply_compute(inArgs.at(0)))*/;
           typedef TENSOR_RESULT_TYPE(T,Q) resultT;
	   typedef resultT R;
	   typedef  WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > dcT;
            R* w1_array;
            R* w2_array;
            Q* w5_array;

            R* rptr_array;
            R* r0ptr_array;
            T* f0ptr_array;
            T* fptr_array;

            unsigned int* w1_offarray = new unsigned int[inArgs.size()];
            unsigned int* w2_offarray = new unsigned int[inArgs.size()];
            unsigned int* w5_offarray = new unsigned int[inArgs.size()];
            
            unsigned int* rptr_offarray = new unsigned int[inArgs.size()];
            unsigned int* r0ptr_offarray = new unsigned int[inArgs.size()];
            unsigned int* f0ptr_offarray = new unsigned int[inArgs.size()];
            unsigned int* fptr_offarray = new unsigned int[inArgs.size()];


            unsigned int w1_off = 0;
            unsigned int w2_off = 0;
            unsigned int w5_off = 0;
            
            unsigned int rptr_off = 0;
            unsigned int r0ptr_off = 0;
            unsigned int f0ptr_off = 0;
            unsigned int fptr_off = 0;


            std::tr1::tuple<keyT, keyT, keyT,
                  double, double, double,
                  Tensor<TENSOR_RESULT_TYPE(T,Q)>,
                  WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > > ** args = new std::tr1::tuple<keyT, keyT, keyT,
                                                                                           double, double, double,
                                                                                           Tensor<TENSOR_RESULT_TYPE(T,Q)>,
                                                                                           WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > > *[inArgs.size()];

            //on CPU
            Tensor<R> * r_array = new Tensor<R>[inArgs.size()];
            Tensor<R> * r0_array = new Tensor<R>[inArgs.size()];
            double * argstol_array = new double[inArgs.size()]; 
            double * argsfac_array = new double[inArgs.size()]; 
            keyT * argsdest_array = new keyT[inArgs.size()];
            dcT ** coeffs_array = new dcT*[inArgs.size()]; 
            int * kref_array = new int[inArgs.size()]; 
            int * rankref_array = new int[inArgs.size()]; 
            double * tol_array = new double[inArgs.size()]; 
            Tensor<T> * f_array = new Tensor<T>[inArgs.size()];
            Tensor<T> * f0_array = new Tensor<T>[inArgs.size()];
            
            unsigned int i;
            for (i = 0; i < inArgs.size(); i++){
               // std::tr1::tuple<Tensor<TENSOR_RESULT_TYPE(T,Q)> *, Tensor<TENSOR_RESULT_TYPE(T,Q)> *,
               //          WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> >, keyT, double, double> temp = inObj.at(i)->apply_computept(inArgs.at(i));

                  /*
                  std::tr1::tuple<keyT, keyT, keyT,
                                  double, double, double,
                                  Tensor<TENSOR_RESULT_TYPE(T,Q)>,
                                  WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > > arg1 = inArgs.at(i);*/
                  args[i] = &(inArgs.at(i));

		  //keyT& source = std::tr1::get<0>(arg1);
		  //keyT& shift = std::tr1::get<1>(arg1);
		  //keyT& argsdest = std::tr1::get<2>(arg1);
                  argsdest_array[i] = std::tr1::get<2>(*args[i]);
		  //double argstol = std::tr1::get<3>(/*arg1*/*args[i]);
		  argstol_array[i] = std::tr1::get<3>(/*arg1*/*args[i]);
		  //double argsfac = std::tr1::get<4>(/*arg1*/*args[i]);
		  argsfac_array[i] = std::tr1::get<4>(/*arg1*/*args[i]);
		  double argscnorm = std::tr1::get<5>(/*arg1*/*args[i]);
		  Tensor<R>& coeff = std::tr1::get<6>(/*arg1*/*args[i]);
		  //dcT& coeffs = std::tr1::get<7>(arg1);
                  coeffs_array[i] = &(std::tr1::get<7>(*args[i]));
                  //int kref = inObj.at(i)->k;
                  kref_array[i] = inObj.at(i)->k;
                  //int rankref = inObj.at(i)->rank;
                  rankref_array[i] = inObj.at(i)->rank;
                  const std::vector<long>& vkref = inObj.at(i)->vk;
                  const std::vector<long>& v2kref = inObj.at(i)->v2k;
                  const std::vector<Slice>& s0ref = inObj.at(i)->s0;
                  //const SimpleCache< SeparatedConvolutionData<Q,NDIM>, NDIM >& data = inObj.at(i)->data;

		  //Key<NDIM>& source = argskey;
		  //Key<NDIM>& shift = argsd;
		  //double tol = argstol/argsfac/argscnorm;
		  tol_array[i] = argstol_array[i]/argsfac_array[i]/argscnorm;

		  //const SeparatedConvolutionData<Q,NDIM>* op = inObj.at(i)->getop(source.level(), shift);

		  Tensor<resultT> r(v2kref);
		  Tensor<resultT> r0(vkref);
                  r_array[i] = r;
                  r0_array[i] = r0;
		  Tensor<resultT> work1(v2kref,false), work2(v2kref,false);
		  Tensor<Q> work5(2*kref_array[i] /*kref*/,2*kref_array[i] /*kref*/);

		  const Tensor<T>* input = &coeff;
		  Tensor<T> dummy;

		  if (coeff.dim(0) == kref_array[i]/*kref*/) {
			// This processes leaf nodes with only scaling
			// coefficients ... FuncImpl::apply by default does not
			// apply the operator to these since for smoothing operators
			// it is not necessary.  It is necessary for operators such
			// as differentiation and time evolution and will also occur
			// if the application of the operator widens the tree.
			dummy = Tensor<T>(v2kref);
			dummy(s0ref) = coeff;
			input = &dummy;
		  }
		  else {
	              MADNESS_ASSERT(coeff.dim(0)==2*kref_array[i] /*kref*/);
		  }

		  //tol = tol/rankref; // Error is per separated term


		  //print("sepop",source,shift,op->norm,tol);


		  const Tensor<T> f0 = copy(coeff(s0ref));
                  f0_array[i] = f0;
		    
		  //Level n = source.level();
		  const Tensor<T>& f = *input;
                  f_array[i] = f;
                  //trans can be made local
		  //Transformation trans[NDIM];
		  /*const long twok = 2*kref;
		  long break_even;
		    
		  if (NDIM==1) break_even = long(0.5*twok);
		  else if (NDIM==2) break_even = long(0.6*twok);
		  else if (NDIM==3) break_even=long(0.65*twok);
		  else break_even=long(0.7*twok);
		    
		  long break_even2;
		  if (NDIM==1) break_even2 = long(0.5*kref);
		  else if (NDIM==2) break_even2 = long(0.6*kref);
		  else if (NDIM==3) break_even2=long(0.65*kref);
		  else break_even2=long(0.7*kref);*/

		  //R* restrict w1=work1.ptr();
		  //R* restrict w2=work2.ptr();
		  //Q* restrict w3=work5.ptr();
			    
                  Tensor<TENSOR_RESULT_TYPE(T,Q)>& result = /*r*/r_array[i];
	          Tensor<TENSOR_RESULT_TYPE(T,Q)>& result0 = /*r0*/r0_array[i];

                  //R* resultptr = result.ptr();
                  //R* result0ptr = result0.ptr();
                  
                  w1_offarray[i] = w1_off;
                  w2_offarray[i] = w2_off;
                  w5_offarray[i] = w5_off;
                  
                  rptr_offarray[i] = rptr_off;
                  r0ptr_offarray[i] = r0ptr_off;
                  f0ptr_offarray[i] = f0ptr_off;
                  fptr_offarray[i] = fptr_off;

                  w1_off += work1.size();
                  w2_off += work2.size();
                  w5_off += work5.size();
                  
                  rptr_off += result.size();
                  r0ptr_off += result0.size();
                  f0ptr_off += f0.size();
                  fptr_off += f.size();
            }

            w1_array = new R[w1_off]; 
            w2_array = new R[w2_off]; 
            w5_array = new Q[w5_off]; 

            rptr_array = new R[rptr_off];
            r0ptr_array = new R[r0ptr_off];
            f0ptr_array = new T[f0ptr_off];
            fptr_array = new T[fptr_off];


            w1_off = 0;
            w2_off = 0;
            w5_off = 0;
            
            rptr_off = 0;
            r0ptr_off = 0;
            f0ptr_off = 0;
            fptr_off = 0;

            for (i = 0; i < inArgs.size(); i++){
               // std::tr1::tuple<Tensor<TENSOR_RESULT_TYPE(T,Q)> *, Tensor<TENSOR_RESULT_TYPE(T,Q)> *,
               //          WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> >, keyT, double, double> temp = inObj.at(i)->apply_computept(inArgs.at(i));

                  //std::tr1::tuple<keyT, keyT, keyT,
                  //                double, double, double,
                  //                Tensor<TENSOR_RESULT_TYPE(T,Q)>,
                  //                WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > > arg1 = inArgs.at(i);

		  //keyT& source = std::tr1::get<0>(arg1);
		  //keyT& shift = std::tr1::get<1>(arg1);
		  //keyT& argsdest = std::tr1::get<2>(arg1);
		  //double argstol = std::tr1::get<3>(arg1);
		  //double argsfac = std::tr1::get<4>(arg1);
		  //double argscnorm = std::tr1::get<5>(arg1);
		  //Tensor<R>& coeff = std::tr1::get<6>(arg1);
		  //dcT& coeffs = std::tr1::get<7>(arg1);
                  //int kref = inObj.at(i)->k;
                  //int rankref = inObj.at(i)->rank;
                  //const std::vector<long>& vkref = inObj.at(i)->vk;
                  const std::vector<long>& v2kref = inObj.at(i)->v2k;
                  //const std::vector<Slice>& s0ref = inObj.at(i)->s0;
                  //const SimpleCache< SeparatedConvolutionData<Q,NDIM>, NDIM >& data = inObj.at(i)->data;

		  //Key<NDIM>& source = argskey;
		  //Key<NDIM>& shift = argsd;
		  //double tol = argstol/argsfac/argscnorm;

		  //const SeparatedConvolutionData<Q,NDIM>* op = inObj.at(i)->getop(source.level(), shift);

		  ////Tensor<resultT> r(v2kref);
		  ////Tensor<resultT> r0(vkref);
		  Tensor<resultT> work1(v2kref,false), work2(v2kref,false);
		  Tensor<Q> work5(2*kref_array[i] /*kref*/,2*kref_array[i] /*kref*/);

		  /*const Tensor<T>* input = &coeff;
		  Tensor<T> dummy;

		  if (coeff.dim(0) == kref) {
			// This processes leaf nodes with only scaling
			// coefficients ... FuncImpl::apply by default does not
			// apply the operator to these since for smoothing operators
			// it is not necessary.  It is necessary for operators such
			// as differentiation and time evolution and will also occur
			// if the application of the operator widens the tree.
			dummy = Tensor<T>(v2kref);
			dummy(s0ref) = coeff;
			input = &dummy;
		  }
		  else {
	              MADNESS_ASSERT(coeff.dim(0)==2*kref);
		  }*/

		  //tol = tol/rankref; // Error is per separated term


		  //print("sepop",source,shift,op->norm,tol);


		  //const Tensor<T> f0 = copy(coeff(s0ref));
		    
		  //Level n = source.level();
		  //const Tensor<T>& f = *input;
                  //trans can be made local
		  //Transformation trans[NDIM];
		  /*const long twok = 2*kref;
		  long break_even;
		    
		  if (NDIM==1) break_even = long(0.5*twok);
		  else if (NDIM==2) break_even = long(0.6*twok);
		  else if (NDIM==3) break_even=long(0.65*twok);
		  else break_even=long(0.7*twok);
		    
		  long break_even2;
		  if (NDIM==1) break_even2 = long(0.5*kref);
		  else if (NDIM==2) break_even2 = long(0.6*kref);
		  else if (NDIM==3) break_even2=long(0.65*kref);
		  else break_even2=long(0.7*kref);*/

		  R* restrict w1=work1.ptr();
		  R* restrict w2=work2.ptr();
		  Q* restrict w3=work5.ptr();
			    
                  Tensor<TENSOR_RESULT_TYPE(T,Q)>& result = /*r*/r_array[i];
	          Tensor<TENSOR_RESULT_TYPE(T,Q)>& result0 = /*r0*/r0_array[i];

                  R* resultptr = result.ptr();
                  R* result0ptr = result0.ptr();
                  T* f0ptr = const_cast<T*>(f0_array[i].ptr());
                  T* fptr = const_cast<T*>(f_array[i].ptr());
                  
                  memcpy(w1_array + w1_offarray[i], w1, work1.size()*sizeof(R));
                  memcpy(w2_array + w2_offarray[i], w2, work2.size()*sizeof(R));
                  memcpy(w5_array + w5_offarray[i], w3, work5.size()*sizeof(Q));
                  
                  memcpy(rptr_array + rptr_offarray[i], resultptr, result.size()*sizeof(R));
                  memcpy(r0ptr_array + r0ptr_offarray[i], result0ptr, result0.size()*sizeof(R));
                  memcpy(f0ptr_array + f0ptr_offarray[i], f0ptr, /*f0*/f0_array[i].size()*sizeof(T));
                  memcpy(fptr_array + fptr_offarray[i], fptr, /*f*/f_array[i].size()*sizeof(T));

                  //memcpy(w1_array + w1_off, w1, w1.size()*sizeof(R));
                  //memcpy(w2_array + w2_off, w2, w2.size()*sizeof(R));
                  //memcpy(w5_array + w5_off, w5, w5.size()*sizeof(Q));
                  
                  //memcpy(rptr_array + rptr_off, resultptr, resultptr.size()*sizeof(R));
                  //memcpy(r0ptr_array + r0ptr_off, result0ptr, result0ptr.size()*sizeof(R));

                  //w1_off += work1.size();
                  //w2_off += work2.size();
                  //w5_off += work5.size();
                  
                  //rptr_off += result.size();
                  //r0ptr_off += result0.size();
                  //f0ptr_off += /*f0*/f0_array[i].size();
                  //fptr_off += /*f*/f_array[i].size();
            }

            for (i = 0; i < inArgs.size(); i++){
               // std::tr1::tuple<Tensor<TENSOR_RESULT_TYPE(T,Q)> *, Tensor<TENSOR_RESULT_TYPE(T,Q)> *,
               //          WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> >, keyT, double, double> temp = inObj.at(i)->apply_computept(inArgs.at(i));

                  //std::tr1::tuple<keyT, keyT, keyT,
                  //                double, double, double,
                  //                Tensor<TENSOR_RESULT_TYPE(T,Q)>,
                  //                WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > > arg1 = inArgs.at(i);

		  keyT& source = std::tr1::get<0>(/*arg1*/*args[i]);
		  keyT& shift = std::tr1::get<1>(/*arg1*/*args[i]);
		  //keyT& argsdest = std::tr1::get<2>(arg1);
		  //double argstol = std::tr1::get<3>(arg1);
		  //double argsfac = std::tr1::get<4>(arg1);
		  //double argscnorm = std::tr1::get<5>(arg1);
		  //Tensor<R>& coeff = std::tr1::get<6>(arg1);
		  //dcT& coeffs = std::tr1::get<7>(arg1);
                  //int kref = inObj.at(i)->k;
                  //int rankref = inObj.at(i)->rank;
                  //const std::vector<long>& vkref = inObj.at(i)->vk;
                  //const std::vector<long>& v2kref = inObj.at(i)->v2k;
                  //const std::vector<Slice>& s0ref = inObj.at(i)->s0;
                  //const SimpleCache< SeparatedConvolutionData<Q,NDIM>, NDIM >& data = inObj.at(i)->data;

		  //Key<NDIM>& source = argskey;
		  //Key<NDIM>& shift = argsd;
		  //double tol = argstol/argsfac/argscnorm;

		  const SeparatedConvolutionData<Q,NDIM>* op = inObj.at(i)->getop(source.level(), shift);

		  ////Tensor<resultT> r(v2kref);
		  ////Tensor<resultT> r0(vkref);
		  ////Tensor<resultT> work1(v2kref,false), work2(v2kref,false);
		  //Tensor<Q> work5(2*kref,2*kref);

                  /*
		  const Tensor<T>* input = &coeff;
		  Tensor<T> dummy;

		  if (coeff.dim(0) == kref) {
			// This processes leaf nodes with only scaling
			// coefficients ... FuncImpl::apply by default does not
			// apply the operator to these since for smoothing operators
			// it is not necessary.  It is necessary for operators such
			// as differentiation and time evolution and will also occur
			// if the application of the operator widens the tree.
			dummy = Tensor<T>(v2kref);
			dummy(s0ref) = coeff;
			input = &dummy;
		  }
		  else {
	              MADNESS_ASSERT(coeff.dim(0)==2*kref);
		  }
                  */

		  //tol = tol/rankref; // Error is per separated term


		  //print("sepop",source,shift,op->norm,tol);


		  //const Tensor<T> f0 = copy(coeff(s0ref));
		    
		  Level n = source.level();
		  ////const Tensor<T>& f = *input;
                  //trans can be made local
		  Transformation trans[NDIM];
		  const long twok = 2*kref_array[i] /*kref*/;
		  long break_even;
		    
		  if (NDIM==1) break_even = long(0.5*twok);
		  else if (NDIM==2) break_even = long(0.6*twok);
		  else if (NDIM==3) break_even=long(0.65*twok);
		  else break_even=long(0.7*twok);
		    
		  long break_even2;
		  if (NDIM==1) break_even2 = long(0.5*kref_array[i] /*kref*/);
		  else if (NDIM==2) break_even2 = long(0.6*kref_array[i] /*kref*/);
		  else if (NDIM==3) break_even2=long(0.65*kref_array[i] /*kref*/);
		  else break_even2=long(0.7*kref_array[i] /*kref*/);

		  //R* restrict w1=work1.ptr();
		  //R* restrict w2=work2.ptr();
		  //Q* restrict w3=work5.ptr();
			    
                  Tensor<TENSOR_RESULT_TYPE(T,Q)>& result = /*r*/r_array[i];
	          Tensor<TENSOR_RESULT_TYPE(T,Q)>& result0 = /*r0*/r0_array[i];

                  R* resultptr;// = result.ptr();
                  R* result0ptr;// = result0.ptr();

                  R* w1ptr;
                  R* w2ptr;
                  Q* w3ptr /*= w3*/  ;
                  T* f0ptr;
                  T* fptr;
                  
                  //R* rptr;
                  //R* r0ptr;

		    for (int mu=0; mu<rankref_array[i] /*rankref*/; ++mu) {
			const SeparatedConvolutionInternal<Q,NDIM>& muop =  op->muops[mu];
			//print("muop",source, shift, mu, muop.norm);
			if (muop.norm > tol_array[i]/*tol*/) {
			    Q fac = inObj.at(i)->ops[mu].getfac();
			    //muopxv_fast(source.level(), muop.ops, *input, f0, r, r0, tol/std::abs(fac), fac,
			    //            work1, work2, work5);

			    //glue
			    const ConvolutionData1D<Q>* const* ops/*[NDIM]*/ = muop.ops;
			    //const Tensor<T>& f0
			    //Tensor<TENSOR_RESULT_TYPE(T,Q)>& result = r;
			    //Tensor<TENSOR_RESULT_TYPE(T,Q)>& result0 = r0;
			    double tol1 = /*tol*/tol_array[i]/std::abs(fac);
			    const Q mufac = fac;
			    //Tensor<TENSOR_RESULT_TYPE(T,Q)>& work1
			    //Tensor<TENSOR_RESULT_TYPE(T,Q)>& work2
			    //Tensor<Q>& work5
			     
			    double Rnorm = 1.0;
			    for (std::size_t d=0; d<NDIM; ++d) Rnorm *= ops[d]->Rnorm;
			    if (Rnorm == 0.0) continue;

			    tol1 = tol1/(Rnorm*NDIM);  // Errors are relative within here

			    // Determine rank of SVD to use or if to use the full matrix
			    for (std::size_t d=0; d<NDIM; ++d) {
				long r1;
				for (r1=0; r1<twok; ++r1) {
				    if (ops[d]->Rs[r1] < tol1) break;
				}
				if (r1 >= break_even) {
				    trans[d].r = twok;
				    trans[d].U = ops[d]->R.ptr();
				    trans[d].VT = 0;
				}
				else {
				    r1 += (r1&1L);
				    trans[d].r = std::max(2L,r1);
				    trans[d].U = ops[d]->RU.ptr();
				    trans[d].VT = ops[d]->RVT.ptr();
				}
			    }
			    ////apply_transformation(n, twok, trans, f, work1, work2, work5, mufac, result);

			    long dimk = twok;
			   
			    long size = 1;
			    for (std::size_t ii=0; ii<NDIM; ++ii) size *= dimk;
			    long dimi = size/dimk;

			    const Q* U;

                            w3ptr = &w5_array[w5_offarray[i]];
			    U = (trans[0].r == dimk) ? trans[0].U : shrink(dimk,dimk,trans[0].r,trans[0].U,/*w3*/w3ptr);
                            w1ptr = &w1_array[w1_offarray[i]];
                            w2ptr = &w2_array[w2_offarray[i]];
                            f0ptr = &f0ptr_array[f0ptr_offarray[i]];
                            fptr = &fptr_array[fptr_offarray[i]];
                            resultptr = &rptr_array[rptr_offarray[i]];
                            result0ptr = &r0ptr_array[r0ptr_offarray[i]];
			    ////GPU
			    mTxmq(dimi, trans[0].r, dimk, /*w1*/w1ptr, /*f.ptr()*/fptr, U);
			    size = trans[0].r * size / dimk;
			    dimi = size/dimk;
			    for (std::size_t d=1; d<NDIM; ++d) {
				U = (trans[d].r == dimk) ? trans[d].U : shrink(dimk,dimk,trans[d].r,trans[d].U,/*w3*/w3ptr);
				////GPU
				mTxmq(dimi, trans[d].r, dimk, /*w2*/w2ptr, /*w1*/w1ptr, U);
				size = trans[d].r * size / dimk;
				dimi = size/dimk;
				////GPU
                                std::swap(w1ptr,w2ptr);
				////std::swap(w1,w2);
			    }

			    // If all blocks are full rank we can skip the transposes
			    bool doit = false;
			    for (std::size_t d=0; d<NDIM; ++d) doit = doit || trans[d].VT;

			    if (doit) {
				for (std::size_t d=0; d<NDIM; ++d) {
				    if (trans[d].VT) {
					dimi = size/trans[d].r;
					////GPU
					mTxmq(dimi, dimk, trans[d].r, /*w2*/w2ptr, /*w1*/w1ptr, trans[d].VT);
					size = dimk*size/trans[d].r;
				    }
				    else {
					////GPU
					fast_transpose(dimk, dimi, /*w1*/w1ptr, /*w2*/w2ptr);
				    }
				    ////GPU
                                    std::swap(w1ptr,w2ptr);
				    ////std::swap(w1,w2);
				}
			    }
			    // Assuming here that result is contiguous and aligned
			    ////GPU
			    aligned_axpy(size, /*result.ptr()*/resultptr, /*w1*/w1ptr, mufac);
			    //    long one = 1;
			    //daxpy_(&size, &mufac, w1, &one, result.ptr(), &one);

			    if (n > 0) {
				for (std::size_t d=0; d<NDIM; ++d) {
				    long r1;
				    for (r1=0; r1< kref_array[i] /*kref*/; ++r1) {
					if (ops[d]->Ts[r1] < tol1) break;
				    }
				    if (r1 >= break_even2) {
					trans[d].r = kref_array[i] /*kref*/;
					trans[d].U = ops[d]->T.ptr();
					trans[d].VT = 0;
				    }
				    else {
					r1 += (r1&1L);
					trans[d].r = std::max(2L,r1);
					trans[d].U = ops[d]->TU.ptr();
					trans[d].VT = ops[d]->TVT.ptr();
				    }
				}
				////apply_transformation(n, k, trans, f0, work1, work2, work5, -mufac, result0);
				dimk = kref_array[i] /*kref*/;
				//const Tensor<T>& f1 = f0;
				const Q mufac1 = -mufac;
				//Tensor<R>& result1 = result0;

				size = 1;
				for (std::size_t i=0; i<NDIM; ++i) size *= dimk;
				dimi = size/dimk;

				const Q* U1;

				U1 = (trans[0].r == dimk) ? trans[0].U : shrink(dimk,dimk,trans[0].r,trans[0].U,/*w3*/w3ptr);
				////GPU
				mTxmq(dimi, trans[0].r, dimk, /*w1*/w1ptr, /*f0.ptr()*/f0ptr, U1);
				size = trans[0].r * size / dimk;
				dimi = size/dimk;
				for (std::size_t d=1; d<NDIM; ++d) {
				    U1 = (trans[d].r == dimk) ? trans[d].U : shrink(dimk,dimk,trans[d].r,trans[d].U,/*w3*/w3ptr);
				    ////GPU
				    mTxmq(dimi, trans[d].r, dimk, /*w2*/w2ptr, /*w1*/w1ptr, U1);
				    size = trans[d].r * size / dimk;
				    dimi = size/dimk;
				    ////GPU
                                    std::swap(w1ptr,w2ptr);
				    ////std::swap(w1,w2);
				}

				// If all blocks are full rank we can skip the transposes
				bool doit = false;
				for (std::size_t d=0; d<NDIM; ++d) doit = doit || trans[d].VT;

				if (doit) {
				    for (std::size_t d=0; d<NDIM; ++d) {
					if (trans[d].VT) {
					    dimi = size/trans[d].r;
					    ////GPU
					    mTxmq(dimi, dimk, trans[d].r, /*w2*/w2ptr, /*w1*/w1ptr, trans[d].VT);
					    size = dimk*size/trans[d].r;
					}
					else {
					    ////GPU
					    fast_transpose(dimk, dimi, /*w1*/w1ptr, /*w2*/w2ptr);
					}
					////GPU
                                        std::swap(w1ptr,w2ptr);
					////std::swap(w1,w2);
				    }
				 }
				 // Assuming here that result is contiguous and aligned
				 ////GPU
				 aligned_axpy(size, /*result0.ptr()*/result0ptr, w1ptr, mufac1);
				 //    long one = 1;
				 //daxpy_(&size, &mufac, w1, &one, result.ptr(), &one);
			    }
			}
		    }
                    memcpy(result.ptr(), resultptr, result.size()*sizeof(R));
                    memcpy(result0.ptr(), result0ptr, result0.size()*sizeof(R));
		    Tensor<R> * r1 = new Tensor<R>(/*r*/r_array[i]); 
		    Tensor<R> * r01 = new Tensor<R>(/*r0*/r0_array[i]); 
		    std::tr1::tuple<Tensor<R>*, Tensor<R>*, dcT, keyT, double, double> t2(r1, r01, /*coeffs*/*coeffs_array[i], /*argsdest*/argsdest_array[i], /*argstol*/argstol_array[i], /*argsfac*/argsfac_array[i]);
                    outArg.push_back(t2);
            }

            delete[] w1_array; //GPU 
            delete[] w2_array; //GPU
            delete[] w5_array; 

            delete[] rptr_array; //GPU
            delete[] r0ptr_array; //GPU
            delete[] f0ptr_array; //GPU
            delete[] fptr_array;  //GPU

            delete[] args; 
            
            delete[] r_array;
            delete[] r0_array;
            delete[] argstol_array;
            delete[] argsfac_array;
            delete[] argsdest_array;
            delete[] coeffs_array;
            delete[] kref_array;
            delete[] rankref_array;
            delete[] tol_array;
            delete[] f_array;
            delete[] f0_array;
            

            return outArg;
        }
        
        template <typename T, typename opT>
        std::vector< std::tr1::tuple< Tensor<TENSOR_RESULT_TYPE(T,Q)> *, Tensor<TENSOR_RESULT_TYPE(T,Q)> *,
                         WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> >, keyT, double, double> >
        apply_allCompute(std::vector<std::tr1::tuple<keyT, keyT, keyT, 
                                      double, double, double, 
                                      Tensor<TENSOR_RESULT_TYPE(T,Q)>, 
                                      WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > > > inArgs, 
                      std::vector< SeparatedConvolution<Q,NDIM>* > inObj) const {

            print(inArgs.size());
            
            std::vector< std::tr1::tuple< Tensor<TENSOR_RESULT_TYPE(T,Q)> *, Tensor<TENSOR_RESULT_TYPE(T,Q)> *,
                         WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> >, keyT, double, double> > outArg/*(inArgs.size(), inObj.at(0)->apply_compute(inArgs.at(0)))*/;
            for (unsigned int i = 0; i < inArgs.size(); i++){
                std::tr1::tuple<Tensor<TENSOR_RESULT_TYPE(T,Q)> *, Tensor<TENSOR_RESULT_TYPE(T,Q)> *,
                         WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> >, keyT, double, double> temp = inObj.at(i)->apply_computept(inArgs.at(i));
                outArg.push_back(temp);
            }
 
            return outArg;
        }
        
        template <typename T>
        Void apply_postprocesspt(std::tr1::tuple< Tensor<TENSOR_RESULT_TYPE(T,Q)> *, Tensor<TENSOR_RESULT_TYPE(T,Q)> *, 
                               WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> >, keyT, double, double> t1) const{
            typedef Tensor<TENSOR_RESULT_TYPE(T,Q)> resultT;
            typedef  WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > dcT;
            resultT* r = std::tr1::get<0>(t1);
            resultT* r0 = std::tr1::get<1>(t1);
            dcT& coeffs = std::tr1::get<2>(t1);
            keyT argsdest = std::tr1::get<3>(t1);
            double argstol = std::tr1::get<4>(t1);
            double argsfac = std::tr1::get<5>(t1);

            (*r)(s0).gaxpy(1.0,*r0,1.0);
            if (r->normf()> 0.3*argstol/argsfac) {
                // OPTIMIZATION NEEDED HERE ... CHANGING THIS TO TASK NOT SEND REMOVED
                // BUILTIN OPTIMIZATION TO SHORTCIRCUIT MSG IF DATA IS LOCAL
                coeffs.task(argsdest, &FunctionNode<T,NDIM>::accumulate, *r, coeffs, argsdest, TaskAttributes::hipri());
            }
            delete r0;
            delete r;
            return None;
        } 


        template <typename T>
        Void apply_postprocess(std::tr1::tuple< Tensor<TENSOR_RESULT_TYPE(T,Q)> *, Tensor<TENSOR_RESULT_TYPE(T,Q)> *, 
                               WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> >&, keyT&, double&, double&> t1) const{
            typedef Tensor<TENSOR_RESULT_TYPE(T,Q)> resultT;
            typedef  WorldContainer<Key<NDIM> , FunctionNode<T, NDIM> > dcT;
            resultT* r = std::tr1::get<0>(t1);
            resultT* r0 = std::tr1::get<1>(t1);
            dcT& coeffs = std::tr1::get<2>(t1);
            const keyT& argsdest = std::tr1::get<3>(t1);
            const double& argstol = std::tr1::get<4>(t1);
            const double& argsfac = std::tr1::get<5>(t1);

            (*r)(s0).gaxpy(1.0,*r0,1.0);
            if (r->normf()> 0.3*argstol/argsfac) {
                // OPTIMIZATION NEEDED HERE ... CHANGING THIS TO TASK NOT SEND REMOVED
                // BUILTIN OPTIMIZATION TO SHORTCIRCUIT MSG IF DATA IS LOCAL
                coeffs.task(argsdest, &FunctionNode<T,NDIM>::accumulate, *r, coeffs, argsdest, TaskAttributes::hipri());
            }
            delete r0;
            delete r;
            return None;
        } 

        template <typename T>
        Tensor<TENSOR_RESULT_TYPE(T,Q)> opt_inlined_apply(const Key<NDIM>& source,
                                              const Key<NDIM>& shift,
                                              const Tensor<T>& coeff,
                                              double tol) const {
            PROFILE_MEMBER_FUNC(SeparatedConvolution);
            print("inlined_apply \n");
            typedef TENSOR_RESULT_TYPE(T,Q) resultT;
            typedef resultT R;
            const Tensor<T>* input = &coeff;
            Tensor<T> dummy;

            if (coeff.dim(0) == k) {
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
                MADNESS_ASSERT(coeff.dim(0)==2*k);
            }

            tol = tol/rank; // Error is per separated term

            const SeparatedConvolutionData<Q,NDIM>* op = getop(source.level(), shift);

            //print("sepop",source,shift,op->norm,tol);

            Tensor<resultT> r(v2k), r0(vk);
            Tensor<resultT> work1(v2k,false), work2(v2k,false);
            Tensor<Q> work5(2*k,2*k);


            const Tensor<T> f0 = copy(coeff(s0));
            
            Level n = source.level();
            const Tensor<T>& f = *input;
            Transformation trans[NDIM];
	    const long twok = 2*k;
	    long break_even;
	    
            if (NDIM==1) break_even = long(0.5*twok);
	    else if (NDIM==2) break_even = long(0.6*twok);
	    else if (NDIM==3) break_even=long(0.65*twok);
	    else break_even=long(0.7*twok);
	    
            long break_even2;
            if (NDIM==1) break_even2 = long(0.5*k);
	    else if (NDIM==2) break_even2 = long(0.6*k);
	    else if (NDIM==3) break_even2=long(0.65*k);
	    else break_even2=long(0.7*k);

            R* restrict w1=work1.ptr();
            R* restrict w2=work2.ptr();
	    Q* restrict w3=work5.ptr();

            for (int mu=0; mu<rank; ++mu) {
                const SeparatedConvolutionInternal<Q,NDIM>& muop =  op->muops[mu];
                //print("muop",source, shift, mu, muop.norm);
                if (muop.norm > tol) {
                    Q fac = ops[mu].getfac();
                    //muopxv_fast(source.level(), muop.ops, *input, f0, r, r0, tol/std::abs(fac), fac,
                    //            work1, work2, work5);

                    //glue
                    const ConvolutionData1D<Q>* const* ops/*[NDIM]*/ = muop.ops;
                    //const Tensor<T>& f0
                    Tensor<TENSOR_RESULT_TYPE(T,Q)>& result = r;
                    Tensor<TENSOR_RESULT_TYPE(T,Q)>& result0 = r0;
                    double tol1 = tol/std::abs(fac);
                    const Q mufac = fac;
                    //Tensor<TENSOR_RESULT_TYPE(T,Q)>& work1
                    //Tensor<TENSOR_RESULT_TYPE(T,Q)>& work2
                    //Tensor<Q>& work5
                     
		    double Rnorm = 1.0;
		    for (std::size_t d=0; d<NDIM; ++d) Rnorm *= ops[d]->Rnorm;
		    if (Rnorm == 0.0) continue;

		    tol1 = tol1/(Rnorm*NDIM);  // Errors are relative within here

		    // Determine rank of SVD to use or if to use the full matrix
		    for (std::size_t d=0; d<NDIM; ++d) {
			long r1;
			for (r1=0; r1<twok; ++r1) {
			    if (ops[d]->Rs[r1] < tol1) break;
			}
			if (r1 >= break_even) {
			    trans[d].r = twok;
			    trans[d].U = ops[d]->R.ptr();
			    trans[d].VT = 0;
			}
			else {
			    r1 += (r1&1L);
			    trans[d].r = std::max(2L,r1);
			    trans[d].U = ops[d]->RU.ptr();
			    trans[d].VT = ops[d]->RVT.ptr();
			}
		    }
		    ////apply_transformation(n, twok, trans, f, work1, work2, work5, mufac, result);

		    long dimk = twok;
		   
                    long size = 1;
		    for (std::size_t i=0; i<NDIM; ++i) size *= dimk;
		    long dimi = size/dimk;

		    const Q* U;

		    U = (trans[0].r == dimk) ? trans[0].U : shrink(dimk,dimk,trans[0].r,trans[0].U,w3);
                    ////GPU
		    mTxmq(dimi, trans[0].r, dimk, w1, f.ptr(), U);
		    size = trans[0].r * size / dimk;
		    dimi = size/dimk;
		    for (std::size_t d=1; d<NDIM; ++d) {
			U = (trans[d].r == dimk) ? trans[d].U : shrink(dimk,dimk,trans[d].r,trans[d].U,w3);
                        ////GPU
			mTxmq(dimi, trans[d].r, dimk, w2, w1, U);
			size = trans[d].r * size / dimk;
			dimi = size/dimk;
                        ////GPU
			std::swap(w1,w2);
		    }

		    // If all blocks are full rank we can skip the transposes
		    bool doit = false;
		    for (std::size_t d=0; d<NDIM; ++d) doit = doit || trans[d].VT;

		    if (doit) {
			for (std::size_t d=0; d<NDIM; ++d) {
			    if (trans[d].VT) {
				dimi = size/trans[d].r;
                                ////GPU
				mTxmq(dimi, dimk, trans[d].r, w2, w1, trans[d].VT);
				size = dimk*size/trans[d].r;
			    }
			    else {
                                ////GPU
				fast_transpose(dimk, dimi, w1, w2);
			    }
                            ////GPU
			    std::swap(w1,w2);
			}
		    }
		    // Assuming here that result is contiguous and aligned
                    ////GPU
		    aligned_axpy(size, result.ptr(), w1, mufac);
		    //    long one = 1;
		    //daxpy_(&size, &mufac, w1, &one, result.ptr(), &one);

		    if (n > 0) {
			for (std::size_t d=0; d<NDIM; ++d) {
			    long r1;
			    for (r1=0; r1<k; ++r1) {
				if (ops[d]->Ts[r1] < tol1) break;
			    }
			    if (r1 >= break_even2) {
				trans[d].r = k;
				trans[d].U = ops[d]->T.ptr();
				trans[d].VT = 0;
			    }
			    else {
				r1 += (r1&1L);
				trans[d].r = std::max(2L,r1);
				trans[d].U = ops[d]->TU.ptr();
				trans[d].VT = ops[d]->TVT.ptr();
			    }
			}
			////apply_transformation(n, k, trans, f0, work1, work2, work5, -mufac, result0);
		        dimk = k;
                        const Tensor<T>& f1 = f0;
                        const Q mufac1 = -mufac;
                        Tensor<R>& result1 = result0;

			size = 1;
			for (std::size_t i=0; i<NDIM; ++i) size *= dimk;
			dimi = size/dimk;

			const Q* U1;

			U1 = (trans[0].r == dimk) ? trans[0].U : shrink(dimk,dimk,trans[0].r,trans[0].U,w3);
                        ////GPU
			mTxmq(dimi, trans[0].r, dimk, w1, f1.ptr(), U1);
			size = trans[0].r * size / dimk;
			dimi = size/dimk;
			for (std::size_t d=1; d<NDIM; ++d) {
	                    U1 = (trans[d].r == dimk) ? trans[d].U : shrink(dimk,dimk,trans[d].r,trans[d].U,w3);
                            ////GPU
			    mTxmq(dimi, trans[d].r, dimk, w2, w1, U1);
			    size = trans[d].r * size / dimk;
			    dimi = size/dimk;
                            ////GPU
		            std::swap(w1,w2);
			}

			// If all blocks are full rank we can skip the transposes
			bool doit = false;
			for (std::size_t d=0; d<NDIM; ++d) doit = doit || trans[d].VT;

			if (doit) {
			    for (std::size_t d=0; d<NDIM; ++d) {
				if (trans[d].VT) {
			            dimi = size/trans[d].r;
                                    ////GPU
				    mTxmq(dimi, dimk, trans[d].r, w2, w1, trans[d].VT);
				    size = dimk*size/trans[d].r;
				}
				else {
                                    ////GPU
			            fast_transpose(dimk, dimi, w1, w2);
				}
                                ////GPU
				std::swap(w1,w2);
		            }
			 }
			 // Assuming here that result is contiguous and aligned
                         ////GPU
			 aligned_axpy(size, result1.ptr(), w1, mufac1);
			 //    long one = 1;
			 //daxpy_(&size, &mufac, w1, &one, result.ptr(), &one);
                    }
                }
            }
            r(s0).gaxpy(1.0,r0,1.0);
            return r;
        }

        template <typename T>
        Tensor<TENSOR_RESULT_TYPE(T,Q)> inlined_apply(const Key<NDIM>& source,
                                              const Key<NDIM>& shift,
                                              const Tensor<T>& coeff,
                                              double tol) const {
            PROFILE_MEMBER_FUNC(SeparatedConvolution);
            print("inlined_apply \n");
            typedef TENSOR_RESULT_TYPE(T,Q) resultT;
            typedef resultT R;
            const Tensor<T>* input = &coeff;
            Tensor<T> dummy;

            if (coeff.dim(0) == k) {
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
                MADNESS_ASSERT(coeff.dim(0)==2*k);
            }

            tol = tol/rank; // Error is per separated term

            const SeparatedConvolutionData<Q,NDIM>* op = getop(source.level(), shift);

            //print("sepop",source,shift,op->norm,tol);

            Tensor<resultT> r(v2k), r0(vk);
            Tensor<resultT> work1(v2k,false), work2(v2k,false);
            Tensor<Q> work5(2*k,2*k);


            const Tensor<T> f0 = copy(coeff(s0));
            for (int mu=0; mu<rank; ++mu) {
                const SeparatedConvolutionInternal<Q,NDIM>& muop =  op->muops[mu];
                //print("muop",source, shift, mu, muop.norm);
                if (muop.norm > tol) {
                    Q fac = ops[mu].getfac();
                    //muopxv_fast(source.level(), muop.ops, *input, f0, r, r0, tol/std::abs(fac), fac,
                    //            work1, work2, work5);

                    //glue
                    Level n = source.level();
                    const ConvolutionData1D<Q>* const* ops/*[NDIM]*/ = muop.ops;
                    const Tensor<T>& f = *input;
                    //const Tensor<T>& f0
                    Tensor<TENSOR_RESULT_TYPE(T,Q)>& result = r;
                    Tensor<TENSOR_RESULT_TYPE(T,Q)>& result0 = r0;
                    double tol1 = tol/std::abs(fac);
                    const Q mufac = fac;
                    //Tensor<TENSOR_RESULT_TYPE(T,Q)>& work1
                    //Tensor<TENSOR_RESULT_TYPE(T,Q)>& work2
                    //Tensor<Q>& work5

                    Transformation trans[NDIM];
                     
		    double Rnorm = 1.0;
		    for (std::size_t d=0; d<NDIM; ++d) Rnorm *= ops[d]->Rnorm;
		    if (Rnorm == 0.0) continue;

		    tol1 = tol1/(Rnorm*NDIM);  // Errors are relative within here

		    // Determine rank of SVD to use or if to use the full matrix
		    const long twok = 2*k;
		    long break_even;
		    if (NDIM==1) break_even = long(0.5*twok);
		    else if (NDIM==2) break_even = long(0.6*twok);
		    else if (NDIM==3) break_even=long(0.65*twok);
		    else break_even=long(0.7*twok);
		    for (std::size_t d=0; d<NDIM; ++d) {
			long r1;
			for (r1=0; r1<twok; ++r1) {
			    if (ops[d]->Rs[r1] < tol1) break;
			}
			if (r1 >= break_even) {
			    trans[d].r = twok;
			    trans[d].U = ops[d]->R.ptr();
			    trans[d].VT = 0;
			}
			else {
			    r1 += (r1&1L);
			    trans[d].r = std::max(2L,r1);
			    trans[d].U = ops[d]->RU.ptr();
			    trans[d].VT = ops[d]->RVT.ptr();
			}
		    }
		    ////apply_transformation(n, twok, trans, f, work1, work2, work5, mufac, result);

		    long dimk = twok;
		    Tensor<Q>& work3 = work5;
		    
                    long size = 1;
		    for (std::size_t i=0; i<NDIM; ++i) size *= dimk;
		    long dimi = size/dimk;

		    R* restrict w1=work1.ptr();
		    R* restrict w2=work2.ptr();
		    Q* restrict w3=work3.ptr();

		    const Q* U;

		    U = (trans[0].r == dimk) ? trans[0].U : shrink(dimk,dimk,trans[0].r,trans[0].U,w3);
		    mTxmq(dimi, trans[0].r, dimk, w1, f.ptr(), U);
		    size = trans[0].r * size / dimk;
		    dimi = size/dimk;
		    for (std::size_t d=1; d<NDIM; ++d) {
			U = (trans[d].r == dimk) ? trans[d].U : shrink(dimk,dimk,trans[d].r,trans[d].U,w3);
			mTxmq(dimi, trans[d].r, dimk, w2, w1, U);
			size = trans[d].r * size / dimk;
			dimi = size/dimk;
			std::swap(w1,w2);
		    }

		    // If all blocks are full rank we can skip the transposes
		    bool doit = false;
		    for (std::size_t d=0; d<NDIM; ++d) doit = doit || trans[d].VT;

		    if (doit) {
			for (std::size_t d=0; d<NDIM; ++d) {
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
		    // Assuming here that result is contiguous and aligned
		    aligned_axpy(size, result.ptr(), w1, mufac);
		    //    long one = 1;
		    //daxpy_(&size, &mufac, w1, &one, result.ptr(), &one);

		    if (n > 0) {
			if (NDIM==1) break_even = long(0.5*k);
			else if (NDIM==2) break_even = long(0.6*k);
			else if (NDIM==3) break_even=long(0.65*k);
			else break_even=long(0.7*k);
			for (std::size_t d=0; d<NDIM; ++d) {
			    long r1;
			    for (r1=0; r1<k; ++r1) {
				if (ops[d]->Ts[r1] < tol1) break;
			    }
			    if (r1 >= break_even) {
				trans[d].r = k;
				trans[d].U = ops[d]->T.ptr();
				trans[d].VT = 0;
			    }
			    else {
				r1 += (r1&1L);
				trans[d].r = std::max(2L,r1);
				trans[d].U = ops[d]->TU.ptr();
				trans[d].VT = ops[d]->TVT.ptr();
			    }
			}
			////apply_transformation(n, k, trans, f0, work1, work2, work5, -mufac, result0);
		        dimk = k;
                        const Tensor<T>& f1 = f0;
			Tensor<Q>& work31 = work5;
                        const Q mufac1 = -mufac;
                        Tensor<R>& result1 = result0;

			size = 1;
			for (std::size_t i=0; i<NDIM; ++i) size *= dimk;
			dimi = size/dimk;

			w1=work1.ptr();
			w2=work2.ptr();
			w3=work31.ptr();

			const Q* U1;

			U1 = (trans[0].r == dimk) ? trans[0].U : shrink(dimk,dimk,trans[0].r,trans[0].U,w3);
			mTxmq(dimi, trans[0].r, dimk, w1, f1.ptr(), U1);
			size = trans[0].r * size / dimk;
			dimi = size/dimk;
			for (std::size_t d=1; d<NDIM; ++d) {
	                    U1 = (trans[d].r == dimk) ? trans[d].U : shrink(dimk,dimk,trans[d].r,trans[d].U,w3);
			    mTxmq(dimi, trans[d].r, dimk, w2, w1, U1);
			    size = trans[d].r * size / dimk;
			    dimi = size/dimk;
		            std::swap(w1,w2);
			}

			// If all blocks are full rank we can skip the transposes
			bool doit = false;
			for (std::size_t d=0; d<NDIM; ++d) doit = doit || trans[d].VT;

			if (doit) {
			    for (std::size_t d=0; d<NDIM; ++d) {
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
			 // Assuming here that result is contiguous and aligned
			 aligned_axpy(size, result1.ptr(), w1, mufac1);
			 //    long one = 1;
			 //daxpy_(&size, &mufac, w1, &one, result.ptr(), &one);
                    }
                }
            }
            r(s0).gaxpy(1.0,r0,1.0);
            return r;
        }

    }; 

    /// Factory function generating separated kernel for convolution with 1/r in 3D.
    static
    inline
    SeparatedConvolution<double_complex,3> PeriodicHFExchangeOperator(World& world,
                                                   Vector<double,3> args,
                                                   double lo,
                                                   double eps,
                                                   const BoundaryConditions<3>& bc=FunctionDefaults<3>::get_bc(),
                                                   int k=FunctionDefaults<3>::get_k())
    {
        const Tensor<double>& cell_width = FunctionDefaults<3>::get_cell_width();
        double hi = cell_width.normf(); // Diagonal width of cell
        if (bc(0,0) == BC_PERIODIC) hi *= 100; // Extend range for periodic summation
        const double pi = constants::pi;

        // bsh_fit generates representation for 1/4Pir but we want 1/r
        // so have to scale eps by 1/4Pi

        Tensor<double> coeff, expnt;
        bsh_fit(0.0, lo, hi, eps/(4.0*pi), &coeff, &expnt, false);

        if (bc(0,0) == BC_PERIODIC) {
            truncate_periodic_expansion(coeff, expnt, cell_width.max(), true);
        }
        coeff.scale(4.0*pi);
        return SeparatedConvolution<double_complex,3>(world, args, coeff, expnt, bc, k, false);
//        return SeparatedConvolution<double_complex,3>(world, coeff, expnt, bc, k);

    }

    /// Factory function generating separated kernel for convolution with 1/r in 3D.
    static
    inline
    SeparatedConvolution<double,3> CoulombOperator(World& world,
                                                   double lo,
                                                   double eps,
                                                   const BoundaryConditions<3>& bc=FunctionDefaults<3>::get_bc(),
                                                   int k=FunctionDefaults<3>::get_k())
    {
        const Tensor<double>& cell_width = FunctionDefaults<3>::get_cell_width();
        double hi = cell_width.normf(); // Diagonal width of cell
        if (bc(0,0) == BC_PERIODIC) hi *= 100; // Extend range for periodic summation
        const double pi = constants::pi;

        // bsh_fit generates representation for 1/4Pir but we want 1/r
        // so have to scale eps by 1/4Pi

        Tensor<double> coeff, expnt;
        bsh_fit(0.0, lo, hi, eps/(4.0*pi), &coeff, &expnt, false);

        if (bc(0,0) == BC_PERIODIC) {
            truncate_periodic_expansion(coeff, expnt, cell_width.max(), true);
        }
        coeff.scale(4.0*pi);
        return SeparatedConvolution<double,3>(world, coeff, expnt, bc, k);
    }


    /// Factory function generating separated kernel for convolution with 1/r in 3D.
    static
    inline
    SeparatedConvolution<double,3>* CoulombOperatorPtr(World& world,
                                                       double lo,
                                                       double eps,
                                                       const BoundaryConditions<3>& bc=FunctionDefaults<3>::get_bc(),
                                                       int k=FunctionDefaults<3>::get_k())
    {
        const Tensor<double>& cell_width = FunctionDefaults<3>::get_cell_width();
        double hi = cell_width.normf(); // Diagonal width of cell
        if (bc(0,0) == BC_PERIODIC) hi *= 100; // Extend range for periodic summation
        const double pi = constants::pi;

        // bsh_fit generates representation for 1/4Pir but we want 1/r
        // so have to scale eps by 1/4Pi
        Tensor<double> coeff, expnt;
        bsh_fit(0.0, lo, hi, eps/(4.0*pi), &coeff, &expnt, false);
        if (bc(0,0) == BC_PERIODIC) {
            truncate_periodic_expansion(coeff, expnt, cell_width.max(), true);
        }
        coeff.scale(4.0*pi);
        return new SeparatedConvolution<double,3>(world, coeff, expnt, bc, k);
    }


    /// Factory function generating separated kernel for convolution with BSH kernel in general NDIM
    template <std::size_t NDIM>
    static
    inline
    SeparatedConvolution<double,NDIM> BSHOperator(World& world,
                                                  double mu,
                                                  double lo,
                                                  double eps,
                                                  const BoundaryConditions<NDIM>& bc=FunctionDefaults<NDIM>::get_bc(),
                                                  int k=FunctionDefaults<NDIM>::get_k())
    {
        const Tensor<double>& cell_width = FunctionDefaults<NDIM>::get_cell_width();
        double hi = cell_width.normf(); // Diagonal width of cell
        if (bc(0,0) == BC_PERIODIC) hi *= 100; // Extend range for periodic summation
        Tensor<double> coeff, expnt;
        bsh_fit_ndim(NDIM, mu, lo, hi, eps, &coeff, &expnt, false);
        if (bc(0,0) == BC_PERIODIC) {
            truncate_periodic_expansion(coeff, expnt, cell_width.max(), false);
        }

        return SeparatedConvolution<double,NDIM>(world, coeff, expnt, bc, k);
    }


    /// Factory function generating separated kernel for convolution with exp(-mu*r)/(4*pi*r) in 3D
    static
    inline
    SeparatedConvolution<double,3> BSHOperator3D(World& world,
                                                 double mu,
                                                 double lo,
                                                 double eps,
                                                 const BoundaryConditions<3>& bc=FunctionDefaults<3>::get_bc(),
                                                 int k=FunctionDefaults<3>::get_k())

    {
        const Tensor<double>& cell_width = FunctionDefaults<3>::get_cell_width();
        double hi = cell_width.normf(); // Diagonal width of cell
        if (bc(0,0) == BC_PERIODIC) hi *= 100; // Extend range for periodic summation
        Tensor<double> coeff, expnt;
        bsh_fit(mu, lo, hi, eps, &coeff, &expnt, false);
        if (bc(0,0) == BC_PERIODIC) {
            truncate_periodic_expansion(coeff, expnt, cell_width.max(), false);
        }
        return SeparatedConvolution<double,3>(world, coeff, expnt, bc, k);
    }

    /// Factory function generating separated kernel for convolution with exp(-mu*r)/(4*pi*r) in 3D
    static
    inline
    SeparatedConvolution<double,3>* BSHOperatorPtr3D(World& world,
                                                     double mu,
                                                     double lo,
                                                     double eps,
                                                     const BoundaryConditions<3>& bc=FunctionDefaults<3>::get_bc(),
                                                     int k=FunctionDefaults<3>::get_k())
    {
        const Tensor<double>& cell_width = FunctionDefaults<3>::get_cell_width();
        double hi = cell_width.normf(); // Diagonal width of cell
        if (bc(0,0) == BC_PERIODIC) hi *= 100; // Extend range for periodic summation
        Tensor<double> coeff, expnt;
        bsh_fit(mu, lo, hi, eps, &coeff, &expnt, false);
        if (bc(0,0) == BC_PERIODIC) {
            truncate_periodic_expansion(coeff, expnt, cell_width.max(), false);
        }
        return new SeparatedConvolution<double,3>(world, coeff, expnt, bc, k);
    }


    /// Factory function generating operator for convolution with grad(1/r) in 3D
    
    /// Returns a 3-vector containing the convolution operator for the
    /// x, y, and z components of grad(1/r)
    static
    inline
    std::vector< std::shared_ptr< SeparatedConvolution<double,3> > >
    GradCoulombOperator(World& world,
                        double lo,
                        double eps,
                        const BoundaryConditions<3>& bc=FunctionDefaults<3>::get_bc(),
                        int k=FunctionDefaults<3>::get_k())
    {
        typedef SeparatedConvolution<double,3> real_convolution_3d;
        typedef std::shared_ptr<real_convolution_3d> real_convolution_3d_ptr;
        const double pi = constants::pi;
        const Tensor<double> width = FunctionDefaults<3>::get_cell_width();
        double hi = width.normf(); // Diagonal width of cell
        const bool isperiodicsum = (bc(0,0)==BC_PERIODIC);
        if (isperiodicsum) hi *= 100; // Extend range for periodic summation
        
        // bsh_fit generates representation for 1/4Pir but we want 1/r
        // so have to scale eps by 1/4Pi
        Tensor<double> coeff, expnt;
        bsh_fit(0.0, lo, hi, eps/(4.0*pi), &coeff, &expnt, false);
        
        if (bc(0,0) == BC_PERIODIC) {
            truncate_periodic_expansion(coeff, expnt, width.max(), true);
        }
        
        coeff.scale(4.0*pi);
        int rank = coeff.dim(0);
        
        std::vector<real_convolution_3d_ptr> gradG(3);
        
        for (int dir=0; dir<3; dir++) {
            std::vector< ConvolutionND<double,3> > ops(rank);
            for (int mu=0; mu<rank; mu++) {
                // We cache the normalized operator so the factor is the value we must multiply
                // by to recover the coeff we want.
                double c = std::pow(sqrt(expnt(mu)/pi),3); // Normalization coeff
                ops[mu].setfac(coeff(mu)/c/width[dir]);
                
                for (int d=0; d<3; d++) {
                    if (d != dir)
                        ops[mu].setop(d,GaussianConvolution1DCache<double>::get(k, expnt(mu)*width[d]*width[d], 0, isperiodicsum));
                }
                ops[mu].setop(dir,GaussianConvolution1DCache<double>::get(k, expnt(mu)*width[dir]*width[dir], 1, isperiodicsum));
            }
            gradG[dir] = real_convolution_3d_ptr(new SeparatedConvolution<double,3>(world, ops));
        }
        
        return gradG;
    }

    namespace archive {
        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveLoadImpl<Archive,const SeparatedConvolution<T,NDIM>*> {
            static inline void load(const Archive& ar, const SeparatedConvolution<T,NDIM>*& ptr) {
                WorldObject< SeparatedConvolution<T,NDIM> >* p;
                ar & p;
                ptr = static_cast< const SeparatedConvolution<T,NDIM>* >(p);
            }
        };

        template <class Archive, class T, std::size_t NDIM>
        struct ArchiveStoreImpl<Archive,const SeparatedConvolution<T,NDIM>*> {
            static inline void store(const Archive& ar, const SeparatedConvolution<T,NDIM>*const& ptr) {
                ar & static_cast< const WorldObject< SeparatedConvolution<T,NDIM> >* >(ptr);
            }
        };
    }

}




#endif // MADNESS_MRA_OPERATOR_H__INCLUDED
