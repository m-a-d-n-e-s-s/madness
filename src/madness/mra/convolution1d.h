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
*/

#ifndef MADNESS_MRA_CONVOLUTION1D_H__INCLUDED
#define MADNESS_MRA_CONVOLUTION1D_H__INCLUDED

#include <madness/world/vector.h>
#include <madness/constants.h>
#include <limits.h>
#include <madness/tensor/tensor.h>
#include <madness/mra/simplecache.h>
#include <madness/mra/adquad.h>
#include <madness/mra/twoscale.h>
#include <madness/tensor/aligned.h>
#include <madness/tensor/tensor_lapack.h>
#include <madness/misc/kahan_accumulator.h>
#include <algorithm>

/// \file mra/convolution1d.h
/// \brief Computes most matrix elements over 1D operators (including Gaussians)

/// \ingroup function

namespace madness {

    void aligned_add(long n, double* MADNESS_RESTRICT a, const double* MADNESS_RESTRICT b);
    void aligned_sub(long n, double* MADNESS_RESTRICT a, const double* MADNESS_RESTRICT b);
    void aligned_add(long n, double_complex* MADNESS_RESTRICT a, const double_complex* MADNESS_RESTRICT b);
    void aligned_sub(long n, double_complex* MADNESS_RESTRICT a, const double_complex* MADNESS_RESTRICT b);

    template <typename T>
    static void copy_2d_patch(T* MADNESS_RESTRICT out, long ldout, const T* MADNESS_RESTRICT in, long ldin, long n, long m) {
        for (long i=0; i<n; ++i, out+=ldout, in+=ldin) {
            for (long j=0; j<m; ++j) {
                out[j] = in[j];
            }
        }
    }

    /// a(n,m) --> b(m,n) ... optimized for smallish matrices
    template <typename T>
    inline void fast_transpose(long n, long m, const T* a, T* MADNESS_RESTRICT b) {
        // n will always be k or 2k (k=wavelet order) and m will be anywhere
        // from 2^(NDIM-1) to (2k)^(NDIM-1).

//                  for (long i=0; i<n; ++i)
//                      for (long j=0; j<m; ++j)
//                          b[j*n+i] = a[i*m+j];
//                  return;

        if (n==1 || m==1) {
            long nm=n*m;
            for (long i=0; i<nm; ++i) b[i] = a[i];
            return;
        }

        long n4 = (n>>2)<<2;
        long m4 = m<<2;
        const T* a0 = a;
        for (long i=0; i<n4; i+=4, a0+=m4) {
            const T* a1 = a0+m;
            const T* a2 = a1+m;
            const T* a3 = a2+m;
            T* MADNESS_RESTRICT bi = b+i;
            for (long j=0; j<m; ++j, bi+=n) {
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

        for (long i=n4; i<n; ++i)
            for (long j=0; j<m; ++j)
                b[j*n+i] = a[i*m+j];

    }

    // /// a(i,j) --> b(i,j) for i=0..n-1 and j=0..r-1 noting dimensions are a(n,m) and b(n,r).

    // /// returns b
    // template <typename T>
    // inline T* shrink(long n, long m, long r, const T* a, T* MADNESS_RESTRICT b) {
    //     T* result = b;
    //     if (r == 0) {
    //         ;
    //     }
    //     else if (r == 1) {
    //         for (long i=0; i<n; ++i) {
    //             b[i] = a[i];
    //         }
    //     }
    //     else if (r == 2) {
    //         for (long i=0; i<n; ++i, a+=m, b+=2) {
    //             b[0] = a[0];
    //             b[1] = a[1];
    //         }
    //     }
    //     else if (r == 3) {
    //         for (long i=0; i<n; ++i, a+=m, b+=3) {
    //             b[0] = a[0];
    //             b[1] = a[1];
    //             b[2] = a[2];
    //         }
    //     }
    //     else if (r == 4) {
    //         for (long i=0; i<n; ++i, a+=m, b+=4) {
    //             b[0] = a[0];
    //             b[1] = a[1];
    //             b[2] = a[2];
    //             b[3] = a[3];
    //         }
    //     }
    //     else {
    //         for (long i=0; i<n; ++i, a+=m, b+=r) {
    //             for (long j=0; j<r; j++) {
    //                 b[j] = a[j];
    //             }
    //         }
    //     }
    //     return result;
    // }

    /// actual data for 1 dimension and for 1 term and for 1 displacement for a convolution operator
    /// here we keep the transformation matrices

    /// !!! Note that if Rnormf is zero then ***ALL*** of the tensors are empty
    template <typename Q>
    struct ConvolutionData1D {


        Tensor<Q> R, T;                 ///< if NS: R=ns, T=T part of ns; if modified NS: T=\uparrow r^(n-1)
        Tensor<Q> RU, RVT, TU, TVT;     ///< SVD approximations to R and T
        Tensor<typename Tensor<Q>::scalar_type> Rs, Ts;     ///< hold relative errors, NOT the singular values..

        // norms for NS form
        double Rnorm, Tnorm, Rnormf, Tnormf, NSnormf;

        // norms for modified NS form
        double N_up, N_diff, N_F;               ///< the norms according to Beylkin 2008, Eq. (21) ff


        /// ctor for NS form
        /// make the operator matrices r^n and \uparrow r^(n-1)
        /// @param[in]  R   operator matrix of the requested level;     NS: unfilter(r^(n+1)); modified NS: r^n
        /// @param[in]  T   upsampled operator matrix from level n-1;   NS: r^n; modified NS: filter( r^(n-1) )
        ConvolutionData1D(const Tensor<Q>& R, const Tensor<Q>& T) : R(R), T(T) {
            Rnormf = R.normf();
            // Making the approximations is expensive ... only do it for
            // significant components
            if (Rnormf > 1e-20) {
                Tnormf = T.normf();
                make_approx(T, TU, Ts, TVT, Tnorm);
                make_approx(R, RU, Rs, RVT, Rnorm);
                int k = T.dim(0);

                Tensor<Q> NS = copy(R);
                for (int i=0; i<k; ++i)
                    for (int j=0; j<k; ++j)
                        NS(i,j) = 0.0;
                NSnormf = NS.normf();

            }
            else {
                Rnorm = Tnorm = Rnormf = Tnormf = NSnormf = 0.0;
                N_F = N_up = N_diff = 0.0;
            }
        }

        /// ctor for modified NS form
        /// make the operator matrices r^n and \uparrow r^(n-1)
        /// @param[in]  R   operator matrix of the requested level;     NS: unfilter(r^(n+1)); modified NS: r^n
        /// @param[in]  T   upsampled operator matrix from level n-1;   NS: r^n; modified NS: filter( r^(n-1) )
        /// @param[in]  modified    use (un) modified NS form
        ConvolutionData1D(const Tensor<Q>& R, const Tensor<Q>& T, const bool modified) : R(R), T(T) {

            // note that R can be small, but T still be large

            Rnormf = R.normf();
            Tnormf = T.normf();
            // Making the approximations is expensive ... only do it for
            // significant components
            if (Rnormf > 1e-20) make_approx(R, RU, Rs, RVT, Rnorm);
            if (Tnormf > 1e-20) make_approx(T, TU, Ts, TVT, Tnorm);

            // norms for modified NS form: follow Beylkin, 2008, Eq. (21) ff
            N_F=Rnormf;
            N_up=Tnormf;
            N_diff=(R-T).normf();
        }



        /// approximate the operator matrices using SVD, and abuse Rs to hold the error instead of
        /// the singular values (seriously, who named this??)
        void make_approx(const Tensor<Q>& R,
                         Tensor<Q>& RU, Tensor<typename Tensor<Q>::scalar_type>& Rs, Tensor<Q>& RVT, double& norm) {
            int n = R.dim(0);
            svd(R, RU, Rs, RVT);
            for (int i=0; i<n; ++i) {
                for (int j=0; j<n; ++j) {
                    RVT(i,j) *= Rs[i];
                }
            }
            for (int i=n-1; i>1; --i) { // Form cumulative sum of norms
                Rs[i-1] += Rs[i];
            }

            norm = Rs[0];
            if (Rs[0]>0.0) { // Turn into relative errors
                double rnorm = 1.0/norm;
                for (int i=0; i<n; ++i) {
                    Rs[i] *= rnorm;
                }
            }
        }
    };

    /// Provides the common functionality/interface of all 1D convolutions

    /// interface for 1 term and for 1 dimension;
    /// the actual data are kept in ConvolutionData1D
    /// Derived classes must implement rnlp, issmall, natural_level
    template <typename Q>
    class Convolution1D {
    public:
        typedef Q opT;  ///< The apply function uses this to infer resultT=opT*inputT
        int k;          ///< Wavelet order
        int npt;        ///< Number of quadrature points (is this used?)
        int maxR;       ///< Number of lattice translations for sum
        double bloch_k;  ///< k in exp(i k R) Bloch phase factor folded into lattice sum
        std::optional<unsigned int> D;  ///< if D is nonnull, kernel range limited to [-D/2,D/2] (in simulation cell units), useful for finite-range convolutions with periodic functions; for infinite-range use lattice summation (maxR > 0)
        Tensor<double> quad_x;
        Tensor<double> quad_w;
        Tensor<double> c;
        Tensor<double> hgT, hg;
        Tensor<double> hgT2k;

        mutable SimpleCache<Tensor<Q>, 1> rnlp_cache;
        mutable SimpleCache<Tensor<Q>, 1> rnlij_cache;
        mutable SimpleCache<ConvolutionData1D<Q>, 1> ns_cache;
        mutable SimpleCache<ConvolutionData1D<Q>, 2> mod_ns_cache;

        bool lattice_summed() const { return maxR != 0; }
        bool range_restricted() const { return D.has_value(); }

        virtual ~Convolution1D() {};

        Convolution1D(int k, int npt, int maxR,
                      double bloch_k = 0.0,
                      std::optional<unsigned int> D = {})
                : k(k)
                , npt(npt)
                , maxR(maxR)
                , quad_x(npt)
                , quad_w(npt)
                , bloch_k(bloch_k)
                , D(D)
        {
            if (range_restricted()) MADNESS_CHECK(!lattice_summed());

            auto success = autoc(k,&c);
            MADNESS_CHECK(success);

            gauss_legendre(npt,0.0,1.0,quad_x.ptr(),quad_w.ptr());
            success = two_scale_hg(k,&hg);
            MADNESS_ASSERT(success);
            hgT = transpose(hg);
            success = two_scale_hg(2*k,&hgT2k);
            MADNESS_ASSERT(success);
            hgT2k = transpose(hgT2k);

            // Cannot construct the coefficients here since the
            // derived class is not yet constructed so cannot call
            // (even indirectly) a virtual method
        }

        /// Compute the projection of the operator onto the double order polynomials
        virtual Tensor<Q> rnlp(Level n, Translation lx) const = 0;

        /// @return true if the block of [r^n_l]_ij is expected to be small
        virtual bool issmall(Level n, Translation lx) const = 0;

        /// @return true if the block of [r^n_l]_ij is expected to be small
        /// @note unlike issmall(), this handles periodicity and range restriction
        bool get_issmall(Level n, Translation lx) const {
          if (lattice_summed()) {
            Translation twon = Translation(1) << n;
            for (int R = -maxR; R <= maxR; ++R) {
              if (!issmall(n, R * twon + lx))
                return false;
            }
            return true;
          } else { // !lattice_summed
            if (!range_restricted())
              return issmall(n, lx);
            else {
              // [r^n_l]_ij = superposition of [r^n_l]_p and [r^n_l-1]_p
              // so rnlij is out of the range only if both rnlp contributions are
              const auto closest_lx = lx<=0 ? lx : lx-1;
              return rnlp_is_zero(n, closest_lx) || issmall(n, lx);
            }
          }
        }

        /// @return true if `[r^n_l]` is zero due to range restriction \p D
        bool rnlp_is_zero(Level n, Translation l) const {
          bool result = false;
          if (range_restricted()) {
            if (n == 0) {
              result = l > 0 || l < -1;
            } else { // n > 0
              if (l >= 0)
                result = (1 << (n - 1)) * Translation(*D) <= l;
              else
                result = (-(1 << (n - 1)) * Translation(*D)) > l;
            }
          }
          return result;
        }

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
        /// \note if `this->range_restricted()==true`, `θ(D/2 - |x-y|) K(x-y)` is used as the kernel
        const Tensor<Q>& rnlij(Level n, Translation lx, bool do_transpose=false) const {
            const Tensor<Q>* p=rnlij_cache.getptr(n,lx);
            if (p) return *p;

            // PROFILE_MEMBER_FUNC(Convolution1D); // Too fine grain for routine profiling

            long twok = 2*k;
            Tensor<Q> R(2*twok);
            R(Slice(0,twok-1)) = get_rnlp(n,lx-1);
            R(Slice(twok,2*twok-1)) = get_rnlp(n,lx);

            R.scale(pow(0.5,0.5*n));
            R = inner(c,R);
            if (do_transpose) R = transpose(R);
            rnlij_cache.set(n,lx,R);
            return *rnlij_cache.getptr(n,lx);
        };


        /// Returns a pointer to the cached modified make_nonstandard form of the operator

        /// @param[in]  op_key  holds the scale and the source and target translations
        /// @return     a pointer to the cached modified make_nonstandard form of the operator
        const ConvolutionData1D<Q>* mod_nonstandard(const Key<2>& op_key) const {

            const Level& n=op_key.level();
            const Translation& sx=op_key.translation()[0];      // source translation
            const Translation& tx=op_key.translation()[1];      // target translation
            const Translation  lx=tx-sx;                        // displacement
            const Translation  s_off=sx%2;
            const Translation  t_off=tx%2;

            // we cache translation and source offset
            const Key<2> cache_key(n, Vector<Translation,2>{lx, s_off} );
            const ConvolutionData1D<Q>* p = mod_ns_cache.getptr(cache_key);
            if (p) return p;

            // for paranoid me
            MADNESS_ASSERT(sx>=0 and tx>=0);


            Tensor<Q> R, T, Rm;
//            if (!get_issmall(n, lx)) {
//                print("no issmall", lx, source, n);

                const Translation lx_half = tx/2 - sx/2;
                const Slice s0(0,k-1), s1(k,2*k-1);
//                print("sx, tx",lx,lx_half,sx, tx,"off",s_off,t_off);

                // this is the operator matrix in its actual level
                R = rnlij(n,lx);

                // this is the upsampled operator matrix
                Rm = Tensor<Q>(2*k,2*k);
                if (n>0) Rm(s0,s0)=rnlij(n-1,lx_half);
                {
                    // PROFILE_BLOCK(Convolution1D_nstran); // Too fine grain for routine profiling
                    Rm = transform(Rm,hg);
                }

                {
                    // PROFILE_BLOCK(Convolution1D_nscopy); // Too fine grain for routine profiling
                    T=Tensor<Q>(k,k);
                    if (t_off==0 and s_off==0) T=copy(Rm(s0,s0));
                    if (t_off==0 and s_off==1) T=copy(Rm(s0,s1));
                    if (t_off==1 and s_off==0) T=copy(Rm(s1,s0));
                    if (t_off==1 and s_off==1) T=copy(Rm(s1,s1));
//                    if (t_off==0 and s_off==0) T=copy(Rm(s0,s0));
//                    if (t_off==1 and s_off==0) T=copy(Rm(s0,s1));
//                    if (t_off==0 and s_off==1) T=copy(Rm(s1,s0));
//                    if (t_off==1 and s_off==1) T=copy(Rm(s1,s1));
                }

                {
                    // PROFILE_BLOCK(Convolution1D_trans); // Too fine grain for routine profiling

                    Tensor<Q> RT(k,k), TT(k,k);
                    fast_transpose(k,k,R.ptr(), RT.ptr());
                    fast_transpose(k,k,T.ptr(), TT.ptr());
                    R = RT;
                    T = TT;
                }

//            }

            mod_ns_cache.set(cache_key,ConvolutionData1D<Q>(R,T,true));
            return mod_ns_cache.getptr(cache_key);
        }

        /// Returns a pointer to the cached make_nonstandard form of the operator
        const ConvolutionData1D<Q>* nonstandard(Level n, Translation lx) const {
            const ConvolutionData1D<Q>* p = ns_cache.getptr(n,lx);
            if (p) return p;

            // PROFILE_MEMBER_FUNC(Convolution1D); // Too fine grain for routine profiling

            Tensor<Q> R, T;
            if (!get_issmall(n, lx)) {
                Translation lx2 = lx*2;
#if 0 // UNUSED VARIABLES
                Slice s0(0,k-1), s1(k,2*k-1);
#endif
                const Tensor<Q> r0 = rnlij(n+1,lx2);
                const Tensor<Q> rp = rnlij(n+1,lx2+1);
                const Tensor<Q> rm = rnlij(n+1,lx2-1);

                R = Tensor<Q>(2*k,2*k);

                // this does not match Eq 18 of the MRAQC appendix .. parity off??
                // either rnlij dox swap bra/ket or hgT include extra parity phases
//                 R(s0,s0) = r0;
//                 R(s1,s1) = r0;
//                 R(s1,s0) = rp;
//                 R(s0,s1) = rm;

                {
                    // PROFILE_BLOCK(Convolution1D_nscopy); // Too fine grain for routine profiling
                    copy_2d_patch(R.ptr(),           2*k, r0.ptr(), k, k, k);
                    copy_2d_patch(R.ptr()+2*k*k + k, 2*k, r0.ptr(), k, k, k);
                    copy_2d_patch(R.ptr()+2*k*k,     2*k, rp.ptr(), k, k, k);
                    copy_2d_patch(R.ptr()       + k, 2*k, rm.ptr(), k, k, k);
                }

                //print("R ", n, lx, R.normf(), r0.normf(), rp.normf(), rm.normf());


                {
                    // PROFILE_BLOCK(Convolution1D_nstran); // Too fine grain for routine profiling
                    R = transform(R,hgT);
                }

                //print("RX", n, lx, R.normf(), r0.normf(), rp.normf(), rm.normf());

                {
                    // PROFILE_BLOCK(Convolution1D_trans); // Too fine grain for routine profiling

                    Tensor<Q> RT(2*k,2*k);
                    fast_transpose(2*k, 2*k, R.ptr(), RT.ptr());
                    R = RT;

                    //print("RT", n, lx, R.normf(), r0.normf(), rp.normf(), rm.normf());

                    //T = copy(R(s0,s0));
                    T = Tensor<Q>(k,k);
                    copy_2d_patch(T.ptr(), k, R.ptr(), 2*k, k, k);
                }

                //print("NS", n, lx, R.normf(), T.normf());
            }

            ns_cache.set(n,lx,ConvolutionData1D<Q>(R,T));

            return ns_cache.getptr(n,lx);
        };

        Q phase(double R) const {
          if constexpr (std::is_arithmetic_v<Q>)
            return 1;
          else
            return exp(Q(0.0,bloch_k*R));
        }


        const Tensor<Q>& get_rnlp(Level n, Translation lx) const {
            const Tensor<Q>* p=rnlp_cache.getptr(n,lx);
            if (p) return *p;

            // PROFILE_MEMBER_FUNC(Convolution1D); // Too fine grain for routine profiling

            long twok = 2*k;
            Tensor<Q> r;

            if (get_issmall(n, lx)) {
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
                // PROFILE_BLOCK(Convolution1Drnlp); // Too fine grain for routine profiling

                if (lattice_summed()) {
                    Translation twon = Translation(1)<<n;
                    r = Tensor<Q>(2*k);
                    for (int R=-maxR; R<=maxR; ++R) {
                        r.gaxpy(1.0, rnlp(n,R*twon+lx), phase(R));
                    }
                }
                else {
                    r = rnlp(n, lx);
                }
            }

            rnlp_cache.set(n, lx, r);
            //print("   SET rnlp", n, lx, r);
            return *rnlp_cache.getptr(n,lx);
        }
    };

    /// Array of 1D convolutions (one / dimension)

    /// data for 1 term and all dimensions
    template <typename Q, int NDIM>
    class ConvolutionND {
        std::array<std::shared_ptr<Convolution1D<Q> >, NDIM> ops;
        Q fac;

    public:
        ConvolutionND() : fac(1.0) {}

        ConvolutionND(const ConvolutionND& other) : fac(other.fac)
        {
          ops = other.ops;
        }

        ConvolutionND(std::shared_ptr<Convolution1D<Q> > op, Q fac=1.0) : fac(fac)
        {
            std::fill(ops.begin(), ops.end(), op);
        }

        void setop(int dim, const std::shared_ptr<Convolution1D<Q> >& op)  {
            ops[dim] = op;
        }

        std::shared_ptr<Convolution1D<Q> > getop(int dim) const  {
            return ops[dim];
        }

        void setfac(Q value) {
            fac = value;
        }

        Q getfac() const {
            return fac;
        }
    };

    // To test generic convolution by comparing with GaussianConvolution1D
    template <typename Q>
    class GaussianGenericFunctor {
    private:
        Q coeff;
        double exponent;
        int m;
        Level natlev;

    public:
        // coeff * exp(-exponent*x^2) * x^m
        GaussianGenericFunctor(Q coeff, double exponent, int m=0)
            : coeff(coeff), exponent(exponent), m(m),
              natlev(Level(0.5*log(exponent)/log(2.0)+1)) {}

        Q operator()(double x) const {
            Q ee = coeff*exp(-exponent*x*x);
            for (int mm=0; mm<m; ++mm) ee *= x;
            return ee;
        }
        Level natural_level() const {return natlev;}
    };


    /// Generic 1D convolution using brute force (i.e., slow) adaptive quadrature for rnlp

    /// Calls op(x) with x in *simulation coordinates* to evaluate the function.
    template <typename Q, typename opT>
    class GenericConvolution1D : public Convolution1D<Q> {
    private:
        opT op;
        long maxl;    ///< At natural level is l beyond which operator is zero
    public:

        GenericConvolution1D() {}

        GenericConvolution1D(int k, const opT& op, int maxR, double bloch_k = 0.0)
            : Convolution1D<Q>(k, 20, maxR, bloch_k), op(op), maxl(LONG_MAX-1) {
            // PROFILE_MEMBER_FUNC(GenericConvolution1D); // Too fine grain for routine profiling

            // For efficiency carefully compute outwards at the "natural" level
            // until several successive boxes are determined to be zero.  This
            // then defines the future range of the operator and also serves
            // to precompute the values used in the rnlp cache.

            Level natl = natural_level();
            int nzero = 0;
            for (Translation lx=0; lx<(1L<<natl); ++lx) {
                const Tensor<Q>& rp = this->get_rnlp(natl, lx);
                const Tensor<Q>& rm = this->get_rnlp(natl,-lx);
                if (rp.normf()<1e-12 && rm.normf()<1e-12) ++nzero;
                if (nzero == 3) {
                    maxl = lx-2;
                    break;
                }
            }
        }

        virtual Level natural_level() const final {return op.natural_level();}

        struct Shmoo {
            typedef Tensor<Q> returnT;
            Level n;
            Translation lx;
            const GenericConvolution1D<Q,opT>& q;

            Shmoo(Level n, Translation lx, const GenericConvolution1D<Q,opT>* q)
                    : n(n), lx(lx), q(*q) {}

            returnT operator()(double x) const {
                int twok = q.k*2;
                double fac = std::pow(0.5,n);
                double phix[twok];
                legendre_scaling_functions(x-lx,twok,phix);
                Q f = q.op(fac*x)*sqrt(fac);
                returnT v(twok);
                for (long p=0; p<twok; ++p) v(p) += f*phix[p];
                return v;
            }
        };

        Tensor<Q> rnlp(Level n, Translation lx) const final {
            return adq1(lx, lx+1, Shmoo(n, lx, this), 1e-12,
                        this->npt, this->quad_x.ptr(), this->quad_w.ptr(), 0);
        }

        bool issmall(Level n, Translation lx) const final {
            if (lx < 0) lx = 1 - lx;
            // Always compute contributions to nearest neighbor coupling
            // ... we are two levels below so 0,1 --> 0,1,2,3 --> 0,...,7
            if (lx <= 7) return false;

            n = natural_level()-n;
            if (n >= 0) lx = lx << n;
            else lx = lx >> n;

            return lx >= maxl;
        }
    };


    /// 1D convolution with (derivative) Gaussian; coeff and expnt given in *simulation* coordinates [0,1]

    /// Note that the derivative is computed in *simulation* coordinates so
    /// you must be careful to scale the coefficients correctly.
    template <typename Q>
    class GaussianConvolution1D : public Convolution1D<Q> {
        // Returns range of Gaussian for periodic lattice sum in simulation coords
        static int maxR(bool periodic, double expnt) {
            if (periodic) {
                return std::max(1,int(sqrt(16.0*2.3/expnt)+1));
            }
            else {
                return 0;
            }
        }
    public:
        const Q coeff;          ///< Coefficient
        const double expnt;     ///< Exponent
        const Level natlev;     ///< Level to evaluate
        const int m;            ///< Order of derivative (0, 1, or 2 only)

        explicit GaussianConvolution1D(int k, Q coeff, double expnt,
        		int m, bool periodic, double bloch_k = 0.0,
                        std::optional<unsigned int> D = {})
            : Convolution1D<Q>(k,k+11,maxR(periodic,expnt),bloch_k, D)
            , coeff(coeff)
            , expnt(expnt)
            , natlev(Level(0.5*log(expnt)/log(2.0)+1))
            , m(m)
        {
            MADNESS_ASSERT(m>=0 && m<=2);
            // std::cout << "GC expnt=" << expnt << " coeff="  << coeff << " natlev=" << natlev << " maxR=" << maxR(periodic,expnt) << std::endl;
            // for (Level n=0; n<5; n++) {
            //     for (Translation l=0; l<(1<<n); l++) {
            //         std::cout << "RNLP " << n << " " << l << " " << this->get_rnlp(n,l).normf() << std::endl;
            //     }
            //     std::cout << std::endl;
            // }
        }

        virtual ~GaussianConvolution1D() {}

        virtual Level natural_level() const final {
            return natlev;
        }

        /// Compute the projection of the operator onto the double order polynomials

        /// The returned reference is to a cached tensor ... if you want to
        /// modify it, take a copy first.
        ///
        /// Return in \c v[p] \c p=0..2*k-1
        /// \code
        /// r(n,l,p) = 2^(-n) * int(K(2^(-n)*(z+l)) * phi(p,z), z=0..1)
        /// \endcode
        /// The kernel is coeff*exp(-expnt*z^2)*z^m (with m>0).  This is equivalent to
        /// \code
        /// r(n,l,p) = 2^(-n*(m+1))*coeff * int( ((d/dz)^m exp(-beta*z^2)) * phi(p,z-l), z=l..l+1)
        /// \endcode
        /// where
        /// \code
        /// beta = alpha * 2^(-2*n)
        /// \endcode
        Tensor<Q> rnlp(Level n, const Translation lx) const final {
            int twok = 2*this->k;
            Tensor<Q> v(twok);       // Can optimize this away by passing in
            KahanAccumulator<Q> v_accumulator[twok];
            constexpr bool use_kahan = false;  // change to true to use Kahan accumulator

            // if outside the range, early return, else update the integration limits
            std::pair<double, double> integration_limits{0,1};
            if (this->range_restricted()) {
              const auto two_to_nm1 = (1ul << n) * 0.5;
              if (lx < 0) {
                integration_limits = std::make_pair(
                    std::min(std::max(-two_to_nm1 * this->D.value() - lx, 0.), 1.), 1.);
              } else {
                integration_limits = std::make_pair(
                    0., std::max(std::min(two_to_nm1 * this->D.value() - lx, 1.), 0.));
              }
              // early return if empty integration range (this indicates that
              // the range restriction makes the kernel zero everywhere in the box)
              if (integration_limits.first == integration_limits.second) {
                MADNESS_ASSERT(this->rnlp_is_zero(n, lx));
                return v;
              }
              else {
                MADNESS_ASSERT(!this->rnlp_is_zero(n, lx));
              }
            }
            // integration range lower bound, upper bound, length
            const auto x0 = integration_limits.first;
            const auto x1 = integration_limits.second;
            const auto L = x1 - x0;

            /* Apply high-order Gauss Legendre onto subintervals

               coeff*int(exp(-beta(x+l)**2) * z^m * phi[p](x),x=0..1);

               The translations internally considered are all +ve, so
               significant pieces will be on the left.  Finish after things
               become insignificant.

               The resulting coefficients are accurate to about 1e-20.
            */

            // Rescale expnt & coeff onto level n so integration range
            // is [l,l+1]
            Q scaledcoeff = coeff*pow(0.5,0.5*n*(2*m+1));

            // Subdivide interval into nbox boxes of length h
            // ... estimate appropriate size from the exponent.  A
            // Gaussian with real-part of the exponent beta falls in
            // magnitude by a factor of 1/e at x=1/sqrt(beta), and by
            // a factor of e^-49 ~ 5e-22 at x=7/sqrt(beta) (and with
            // polyn of z^2 it is 1e-20).  So, if we use a box of size
            // 1/sqrt(beta) we will need at most 7 boxes.  Incorporate
            // the coefficient into the screening since it may be
            // large.  We can represent exp(-x^2) over [l,l+1] with a
            // polynomial of order 21 to a relative
            // precision of better than machine precision for
            // l=0,1,2,3 and for l>3 the absolute error is less than
            // 1e-23.  We want to compute matrix elements with
            // polynomials of order 2*k-1+m, so the total order is
            // 2*k+20+m, which can be integrated with a quadrature rule
            // of npt=k+11+(m+1)/2.  npt is set in the constructor.

            double fourn = std::pow(4.0,double(n));
            double beta = expnt * pow(0.25,double(n));
            double h = 1.0/sqrt(beta);  // 2.0*sqrt(0.5/beta);
            long nbox = long(1.0/h);
            if (nbox < 1) nbox = 1;
            h = L/nbox;

            // Find argmax such that h*scaledcoeff*exp(-argmax)=1e-22 ... if
            // beta*xlo*xlo is already greater than argmax we can neglect this
            // and subsequent boxes.

            // The derivatives add a factor of expnt^m to the size of
            // the function at long range.
            double sch = std::abs(scaledcoeff*h);
            if (m == 1) sch *= expnt;
            else if (m == 2) sch *= expnt*expnt;
            double argmax = std::abs(log(1e-22/sch)); // perhaps should be -log(1e-22/sch) ?

            // to screen need to iterate over boxes in the order of decreasing kernel values
            const bool left_to_right = lx >= 0;
            // if going left-to-right, start at left, else at right
            const double xstartedge = left_to_right ? x0+lx : lx + 1;

            // with oscillatory integrands the heuristic for reducing roundoff
            // is to sum from large to small, i.e. proceed in same direction as the order of boxes
            // WARNING: the grid points in quad_{x,w} are in order of decreasing x!
            // hence decrement grid point indices for left_to_right, increment otherwise
            const long first_pt = left_to_right ? this->npt-1: 0;
            const long sentinel_pt = left_to_right ? -1 : this->npt;
            const auto next_pt = [lx, left_to_right](auto i) { return left_to_right ? i-1 : i+1; };

            double xlo = left_to_right ? xstartedge : xstartedge-h;
            double xhi;
            for (long box=0; box!=nbox; ++box, xlo = (left_to_right ? xhi : xlo-h)) {

                // can ignore this and rest of boxes if the Gaussian has decayed enough at the side of the box closest to the origin
                xhi=xlo+h;
                const auto xabs_min = std::min(std::abs(xhi),std::abs(xlo));
                if (beta*xabs_min*xabs_min > argmax) break;

                for (long i=first_pt; i!=sentinel_pt; i=next_pt(i)) {
#ifdef IBMXLC
                    double phix[80];
#else
                    double phix[twok];
#endif
                    double xx = xlo + h*this->quad_x(i);
                    Q ee = scaledcoeff*exp(-beta*xx*xx)*this->quad_w(i)*h;

                    // Differentiate as necessary
                    if (m == 1) {
                        ee *= -2.0*expnt*xx;
                    }
                    else if (m == 2) {
                        ee *= (4.0*xx*xx*expnt*expnt - 2.0*expnt*fourn);
                    }

                    legendre_scaling_functions(xx-lx,twok,phix);
                    for (long p=0; p<twok; ++p) {
                      if constexpr (use_kahan)
                        v_accumulator[p] += ee * phix[p];
                      else
                        v(p) += ee * phix[p];
                    }
                }
            }

            if constexpr (use_kahan) {
              for (long p = 0; p < twok; ++p)
                v(p) = static_cast<Q>(v_accumulator[p]);
            }

            return v;
        }

        /// @return true if the block of [r^n_l]_ij is expected to be small
        bool issmall(Level n, Translation lx) const final {
            // [r^n_l]_ij = superposition of [r^n_l]_p and [r^n_l-1]_ij
            // lx>0? the nearest box is lx-1 -> the edge closest to the origin is lx - 1
            // lx<0? the nearest box is lx -> the edge closest to the origin is lx + 1
            // lx==0? interactions within same box  are never small
            if (lx == 0)
              return false;
            else {
              const double beta = expnt * pow(0.25,double(n));
              const double overly_large_beta_r2 = 49.0;      // 49 -> 5e-22     69 -> 1e-30
              if (lx > 0) {
                const auto ll = lx - 1;
                return beta * ll * ll > overly_large_beta_r2;
              } else {
                const auto ll = lx + 1;
                return beta * ll * ll > overly_large_beta_r2;
              }
            }
        };
    };


    template <typename Q>
    struct GaussianConvolution1DCache {
        static ConcurrentHashMap<hashT, std::shared_ptr< GaussianConvolution1D<Q> > > map;
        typedef typename ConcurrentHashMap<hashT, std::shared_ptr< GaussianConvolution1D<Q> > >::iterator iterator;
        typedef typename ConcurrentHashMap<hashT, std::shared_ptr< GaussianConvolution1D<Q> > >::datumT datumT;

        static std::shared_ptr< GaussianConvolution1D<Q> > get(int k, double expnt, int m, bool periodic,
                                                               double bloch_k = 0.0,
                                                               std::optional<unsigned int> D = {}) {
            hashT key = hash_value(expnt);
            hash_combine(key, k);
            hash_combine(key, m);
            hash_combine(key, int(periodic));
            hash_combine(key, bloch_k);
            if (D) hash_combine(key, *D);

            MADNESS_PRAGMA_CLANG(diagnostic push)
            MADNESS_PRAGMA_CLANG(diagnostic ignored "-Wundefined-var-template")

            iterator it = map.find(key);
            if (it == map.end()) {
                [[maybe_unused]] auto&& [tmpit, inserted] = map.insert(datumT(key, std::make_shared< GaussianConvolution1D<Q> >(k,
                                                                                    Q(sqrt(expnt/constants::pi)),
                                                                                    expnt,
                                                                                    m,
                                                                                    periodic,
                                                                                    bloch_k,
                                                                                    D
                                                                                    )));
                MADNESS_ASSERT(inserted);
                it = map.find(key);
                //printf("conv1d: making  %d %.8e\n",k,expnt);
            }
            else {
                //printf("conv1d: reusing %d %.8e\n",k,expnt);
            }
            auto& result = it->second;
            MADNESS_ASSERT(result->expnt == expnt &&
                           result->k == k &&
                           result->m == m &&
                           result->lattice_summed() == periodic &&
                           result->D == D &&
                           result->bloch_k == bloch_k);
            return result;

            MADNESS_PRAGMA_CLANG(diagnostic pop)

        }
    };

    // instantiated in mra1.cc
    template <>
    ConcurrentHashMap< hashT, std::shared_ptr< GaussianConvolution1D<double> > >
        GaussianConvolution1DCache<double>::map;

    // instantiated in mra1.cc
    template <>
    ConcurrentHashMap< hashT, std::shared_ptr< GaussianConvolution1D<double_complex> > >
        GaussianConvolution1DCache<double_complex>::map;

}

#endif // MADNESS_MRA_CONVOLUTION1D_H__INCLUDED
