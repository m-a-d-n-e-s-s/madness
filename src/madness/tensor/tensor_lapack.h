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

  
#ifndef MADNESS_LINALG_TENSOR_LAPACK_H__INCLUDED
#define MADNESS_LINALG_TENSOR_LAPACK_H__INCLUDED

#include <madness/tensor/tensor.h>
#include <madness/fortran_ctypes.h>

/*!
  \file tensor_lapack.h
  \brief Prototypes for a partial interface from Tensor to LAPACK
  \ingroup linalg
@{
*/

namespace madness {

    /// Computes singular value decomposition of matrix
    
    /// \ingroup linalg
    template <typename T>
    void svd(const Tensor<T>& a, Tensor<T>& U,
             Tensor< typename Tensor<T>::scalar_type >& s, Tensor<T>& VT);

    /// \ingroup linalg
    template <typename T>
    void svd_result(Tensor<T>& a, Tensor<T>& U,
             Tensor< typename Tensor<T>::scalar_type >& s, Tensor<T>& VT, Tensor<T>& work);

    /// SVD - MATLAB syntax

    /// call as
    /// auto [U,s,VT] = svd(A);
    /// with A=U*S*VT
    template <typename T>
    std::tuple<Tensor<T>, Tensor< typename Tensor<T>::scalar_type >, Tensor<T>>
    svd(const Tensor<T>& A) {
        Tensor<T> U,VT;
        Tensor< typename Tensor<T>::scalar_type > s;
        svd(A,U,s,VT);
        return std::make_tuple(U,s,VT);
    }

    /// Solves linear equations
    
    /// \ingroup linalg
    template <typename T>
    void gesv(const Tensor<T>& a, const Tensor<T>& b, Tensor<T>& x);

    /// Solves linear equations using least squares
    
    /// \ingroup linalg
    template <typename T>
    void gelss(const Tensor<T>& a, const Tensor<T>& b, double rcond,
               Tensor<T>& x, Tensor< typename Tensor<T>::scalar_type >& s,
               long &rank, Tensor<typename Tensor<T>::scalar_type>& sumsq);

    /// Solves symmetric or Hermitian eigenvalue problem
    
    /// \ingroup linalg
    template <typename T>
    void syev(const Tensor<T>& A,
              Tensor<T>& V, Tensor< typename Tensor<T>::scalar_type >& e);

    /// Solves symmetric or Hermitian eigenvalue problem - MATLAB syntax

    /// call as
    /// auto [eval, evec] = syev(A);
    template <typename T>
    std::tuple<Tensor< typename Tensor<T>::scalar_type >, Tensor<T>>
    syev(const Tensor<T>& A) {
    	Tensor<T> V;
    	Tensor< typename Tensor<T>::scalar_type > e;
    	syev(A,V,e);
    	return std::make_tuple(e,V);
    }
    // START BRYAN ADDITION
    /// Solves non-symmetric or non-Hermitian eigenvalue problem

    template <typename T>
    void geev(const Tensor<T>& A, Tensor<T>& VR, Tensor<std::complex<T>>& e);

    /// Solves non-symmetric or non-Hermitian generalized eigenvalue problem

    template <typename T>
    struct complex_type {
        typedef std::complex<T> type;
    };

    template <typename T>
    struct complex_type<std::complex<T>> {
        typedef std::complex<T> type;
    };
    
    template <typename T>
    struct real_type {
        typedef T type;
    };

    template <typename T>
    struct real_type<std::complex<T>> {
        typedef T type;
    };
    
    template <typename T>
    void ggev(const Tensor<T>& A, const Tensor<T>& B, Tensor<typename complex_type<T>::type>& VR,
              Tensor<typename complex_type<T>::type>& e);
    // template <typename T>
    // void ggev(const Tensor<T>& A, Tensor<T>& B, Tensor<T>& VR,
    //           Tensor<std::complex<T>>& e);
    // END BRYAN ADDITIONS
    /// Solves symmetric or Hermitian eigenvalue problem - MATLAB syntax

    /// Solves linear equations
    
    /// \ingroup linalg
    template <typename T>
    void gesv(const Tensor<T>& a, const Tensor<T>& b, Tensor<T>& x);


    /// Solves symmetric or Hermitian generalized eigenvalue problem
    
    /// \ingroup linalg
    template <typename T>
    void sygv(const Tensor<T>& A, const Tensor<T>& B, int itype,
              Tensor<T>& V, Tensor< typename Tensor<T>::scalar_type >& e);

    class World; // UGH!
    /// Solves symmetric or Hermitian generalized eigenvalue problem
    

    // !!!!!!!!!! sygvp and gesvp are now in the ELEMENTAL inteface
    // /// \ingroup linalg
    // template <typename T>
    // void sygvp(World& world, const Tensor<T>& A, const Tensor<T>& B, int itype,
    //           Tensor<T>& V, Tensor< typename Tensor<T>::scalar_type >& e);

    // /// Solves linear equations
    
    // /// \ingroup linalg
    // template <typename T>
    // void gesvp(World& world, const Tensor<T>& a, const Tensor<T>& b, Tensor<T>& x);

    /// Cholesky factorization
    
    /// \ingroup linalg
    template <typename T>
    void cholesky(Tensor<T>& A);

    /// rank-revealing Cholesky factorization

    /// \ingroup linalg
    template <typename T>
    void rr_cholesky(Tensor<T>& A, typename Tensor<T>::scalar_type tol, Tensor<integer>& piv, int& rank);

    /// \ingroup linalg
    template <typename T>
    Tensor<T> inverse(const Tensor<T>& A);

    /// QR decomposition
    template<typename T>
    void qr(Tensor<T>& A, Tensor<T>& R);

    /// LQ decomposition
    template<typename T>
    void lq(Tensor<T>& A, Tensor<T>& L);
    /// LQ decomposition
    template<typename T>
    void lq_result(Tensor<T>& A, Tensor<T>& R, Tensor<T>& tau, Tensor<T>& work,bool do_qr);

    template <typename T>
    void geqp3(Tensor<T>& A, Tensor<T>& tau, Tensor<integer>& jpvt);

    /// orgqr generates an M-by-N complex matrix Q with orthonormal columns

    /// which is defined as the first N columns of a product of K elementary
    /// reflectors of order M
    ///       Q  =  H(1) H(2) . . . H(k)
    /// as returned by ZGEQRF.
    template <typename T>
    void orgqr(Tensor<T>& A, const Tensor<T>& tau);


    /// Dunno
    
//     /// \ingroup linalg
//     template <typename T>
//     void triangular_solve(const Tensor<T>& L, Tensor<T>& B, 
//                           const char* side, const char* transa);

    /// Runs the tensor test code, returns true on success
    
    /// \ingroup linalg
    bool test_tensor_lapack();

    /// World/MRA initialization calls this before going multithreaded due to static data in \c dlamch
    
    /// \ingroup linalg
    void init_tensor_lapack();
}

#include <madness/tensor/elem.h>

#endif // MADNESS_LINALG_TENSOR_LAPACK_H__INCLUDED
