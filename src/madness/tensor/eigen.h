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


#ifdef MADNESS_HAS_EIGEN3

#include <madness/tensor/tensor.h>
using madness::Tensor;

#include <iostream>
using std::cout;
using std::endl;

#include <algorithm>
using std::min;
using std::max;

#include <Eigen/Dense>
using namespace Eigen;


#include <madness/world/print.h>

#ifdef MADNESS_FORINT
typedef MADNESS_FORINT integer;
#else
typedef long integer;
#endif //MADNESS_FORINT


/*
 The bigest problem with the Eigen3 and Madness interface 
 is in the transposed copy of the matrixes.
*/

namespace madness {
    /** \brief   Translate between MADNESS Tensor to EIGEN3 Matrix classes and viceversa

    Given N elements, this function copy the values from a Tensor in
    MADNESS to a Matrix in Eigen, and viceversa. The coefficients of 
    both systems must be of the same type.

    */

//vama   template <typename T>
//vama   void copy_mad2eigen2(long n, const T * restrict p1, T * restrict  p0){
//vama    for (long j=0; j<n; ++j,++p1,++p0) {
//vama      *p0 = *p1;
//vama    }
//vama   }

   template <typename T>
   void copy_mad2eigen2(integer n, const T * restrict p1, T * restrict  p0){
#if 0
    //vama #ifdef HAVE_MEMMOVE
    /* Jeff assumes that HAVE_MEMMOVE is a good approximation to HAVE_MEMCPY */
    memcpy(p0,p1,n); /* memcpy is often significantly faster than what the compiler comes up with */
#else
    for (integer j=0; j<n; ++j,++p1,++p0) { *p0 = *p1; } /*JEFF: Is there a reason to not use memcpy? */
#endif
   }

    /** \brief   Compute the singluar value decomposition of an n-by-m matrix using JacobiSVD

    Returns via arguments U, s, VT where

    A = U * diag(s) * VT    for A real
    A = U * diag(s) * VH    for A complex

    or

    UT * A * V = diag(s)   for A real
    UH * A * V = diag(s)   for A complex

    If A is [m,n] and r=min(m,n) then we have U[m,r], s[r], and VT[r,n]
  
    ColPivHouseholderQRPreconditioner is the default. In practice it's very safe.
    It uses column-pivoting QR.

    */


    template <typename T>
    void svd(const Tensor<T>& a, Tensor<T>& U,
             Tensor< typename Tensor<T>::scalar_type >& s, Tensor<T>& VT) {
        TENSOR_ASSERT(a.ndim() == 2, "svd requires matrix",a.ndim(),&a);
        integer m = a.dim(0), n = a.dim(1), rmax = min<int>(m,n);

        s = Tensor< typename Tensor<T>::scalar_type >(rmax);
        U = Tensor<T>(m,rmax);
        VT = Tensor<T>(rmax,n);

        Eigen::Matrix<T, Dynamic, Dynamic> g(n,m);

        //a = transpose(a);
        copy_mad2eigen2( a.size(), a.ptr(), g.data());

        //g.transposeInPlace();
        JacobiSVD< Eigen::Matrix<T, Dynamic, Dynamic>,HouseholderQRPreconditioner> svdm(g, ComputeThinU | ComputeThinV);

        Eigen::Matrix<T, Dynamic, Dynamic> dU = svdm.matrixV();
        Eigen::Matrix<T, Dynamic, Dynamic> dV = svdm.matrixU();

        copy_mad2eigen2( s.size(), svdm.singularValues().data(), s.ptr());
        //g.transposeInPlace();
        dU.transposeInPlace();
        copy_mad2eigen2( VT.size(), dV.data(), VT.ptr());
        copy_mad2eigen2( U.size(), dU.data(), U.ptr());

        //U = transpose(U);
        //U = conj_transpose(U);  //fool

    }




    /** \brief  Solve Ax = b for general A using the EIGEN3 JacobiSVD routine.

    A should be a matrix (float, double, float_complex,
    double_complex) and b should be either a vector, or a matrix with
    each vector stored in a column (i.e., b[n,nrhs]).

    The input A and b are unchanged.  
    It will solve Ax=b as written.

    RANK (output) REAL The number of singular values that are not exactly 0.
    S    (output) REAL Singular values of A.

    The returned solution is guaranteed to minimize the Euclidian norm
    ||Ax-b||

    */
    template <typename T>
    void gelss(const Tensor<T>& a, const Tensor<T>& b, double rcond,
               Tensor<T>& x, Tensor< typename Tensor<T>::scalar_type >& s,
               long& rank, Tensor<typename Tensor<T>::scalar_type>& sumsq) {
        TENSOR_ASSERT(a.ndim() == 2, "gelss requires matrix",a.ndim(),&a);
        integer n = a.dim(1), nrhs = b.dim(1);
        TENSOR_ASSERT(b.ndim() <= 2, "gelss require a vector or matrix for the RHS",b.ndim(),&b);
        TENSOR_ASSERT(a.dim(0) == b.dim(0), "gelss matrix and RHS must conform",b.ndim(),&b);
//
        s = Tensor< typename Tensor<T>::scalar_type >(n);

        Tensor<T> AT = transpose(a);
        Eigen::Matrix<T, Dynamic, Dynamic> g(n,n);
        copy_mad2eigen2( AT.size(), AT.ptr(), g.data());

       Tensor<T> bT;

       if (nrhs == 1) {
            x = Tensor<T>(n);
            bT = Tensor<T>(n);
            bT = copy(b);
            x = bT; //creating the correct size
       }
       else {
            x = Tensor<T>(n,nrhs);
            bT = Tensor<T>(n,nrhs);
            bT =  transpose(b);
       }

        Eigen::Matrix<T, Dynamic, Dynamic> h(n,nrhs);
        copy_mad2eigen2( b.size(), bT.ptr(), h.data());

        rank = 0;
        JacobiSVD< Eigen::Matrix<T, Dynamic, Dynamic> > svdm(g, ComputeThinU | ComputeThinV);
        Eigen::Matrix<T, Dynamic, Dynamic> sol = svdm.solve(h);
        sol.transposeInPlace();

        copy_mad2eigen2( b.size(), sol.data(), x.ptr());

        rank =  svdm.singularValues().nonZeros(); //    in lapack this is different
                                                  // S(i) <= rcond*S(1) are treated
                                                  //                        as zero 
        copy_mad2eigen2( s.size(), svdm.singularValues().data(), s.ptr());
     }

    /** \brief   Real-symmetric or complex-Hermitian eigenproblem.

    A is a real symmetric or complex Hermitian matrix.  Return V and e
    where V is a matrix whose columns are the eigenvectors and e is a
    vector containing the corresponding eigenvalues.  
    The eigenvalues are sorted into ascending
    order. 
    The EIGEN3 SelfAdjointEigenSolver routine is used to address this problem.

    The reults will satisfy A*V(_,i) = V(_,i)*e(i).
    */
    template <typename T>
    void syev(const Tensor<T>& a,
              Tensor<T>& V, Tensor< typename Tensor<T>::scalar_type >& e) {
        TENSOR_ASSERT(a.ndim() == 2, "syev requires a matrix",a.ndim(),&a);
        TENSOR_ASSERT(a.dim(0) == a.dim(1), "syev requires square matrix",0,&a);
        integer n = a.dim(0);

        V = transpose(a);               // For Hermitian case
////
        Eigen::Matrix<T, Dynamic, Dynamic> g(n,n);
        copy_mad2eigen2( a.size(), a.ptr(), g.data());
        Eigen::Matrix<T, Dynamic, Dynamic> gT = g.transpose();
////
        SelfAdjointEigenSolver< Eigen::Matrix<T, Dynamic, Dynamic> > sol(gT);
////
        Eigen::Matrix<T, Dynamic, Dynamic> ev = sol.eigenvectors();
//        cout << " ahy" << endl;
//        cout << ev.transpose() << endl;
//        cout << endl <<sol.eigenvectors() << endl;
////
        V = Tensor<T>(n,n); //fool
        e = Tensor< typename Tensor<T>::scalar_type  >(n);
////

        copy_mad2eigen2( e.size(), sol.eigenvalues().data(), e.ptr());
        //copy_mad2eigen2( e.size(), sol.eigenvalues().data(), e.ptr());

        //ev.transposeInPlace();
        copy_mad2eigen2( V.size(), ev.data(), V.ptr());
//
        //V = transpose(V);
        V = conj_transpose(V);
    }




    /** \brief   Cholesky factorization

    Compute the Cholesky factorization of the symmetric positive definite matrix A
    This function uses the EIGEN3 LLT routine.
    For memory efficiency A is modified inplace.  Its upper
    triangle will hold the result and the lower trianlge will be
    zeroed such that input = inner(transpose(output),output).

    */
    template <typename T>
    void cholesky(Tensor<T>& a) {
        integer n = a.dim(0);

        Eigen::Matrix<T, Dynamic, Dynamic> g(n,n);
        copy_mad2eigen2( a.size(), a.ptr(), g.data());
        g.transposeInPlace();

        LLT< Eigen::Matrix<T, Dynamic, Dynamic> > ColDec(g);
        Eigen::Matrix<T, Dynamic, Dynamic> L = ColDec.matrixU(); // remember the transpose effect U->L ;^

        copy_mad2eigen2( a.size(), L.data(), a.ptr());
        a = transpose(a); 

        for (int i=0; i<n; ++i)
            for (int j=0; j<i; ++j)
                a(i,j) = 0.0;
    }

    /** \brief  Solve Ax = b for general A using the EIGEN3 PartialPivLU routine.

    A should be a square matrix (float, double, float_complex,
    double_complex) and b should be either a vector, or a matrix with
    each vector stored in a column (i.e., b[n,nrhs]).

    It will solve Ax=b as written.

    Use LU decomposition of a matrix with partial pivoting.
    A must be invertible. A is decomposed as A=PLU
    b can be a vector or a matrix, the only restriction is that satisfies b.rows()==A.rows()

    */
    template <typename T>
    void gesv(const Tensor<T>& a, const Tensor<T>& b, Tensor<T>& x) {
        TENSOR_ASSERT(a.ndim() == 2, "gesve requires matrix",a.ndim(),&a);
        integer n = a.dim(0), m = a.dim(1), nrhs = b.dim(1);
        TENSOR_ASSERT(m == n, "gesve requires square matrix",0,&a);
        TENSOR_ASSERT(b.ndim() <= 2, "gesve require a vector or matrix for the RHS",b.ndim(),&b);
        TENSOR_ASSERT(a.dim(0) == b.dim(0), "gesve matrix and RHS must conform",b.ndim(),&b);
/*
   A must be invertible, if you dont know better use FullPivLu
   (check it before)
*/

        typedef typename TensorTypeData<T>::scalar_type scalar_type;
        Tensor<T> AT = transpose(a);

        Eigen::Matrix<T, Dynamic, Dynamic> g(n,n);
        copy_mad2eigen2( AT.size(), AT.ptr(), g.data());

        Tensor<T> bT;
        if (nrhs == 1) {
             x = Tensor<T>(n);
             bT = Tensor<T>(n);
             bT = copy(b);
             x = bT; //creating the correct size
        } 
        else {
             x = Tensor<T>(n,nrhs);
             bT =  transpose(b);
        }

        Eigen::Matrix<T, Dynamic, Dynamic> h(n,nrhs);
        copy_mad2eigen2( b.size(), bT.ptr(), h.data());

        Eigen::Matrix<T, Dynamic, Dynamic> sol = g.householderQr().solve(h);

//vama        Matrix<T, Dynamic, Dynamic> sol = g.lu().solve(h);
        sol.transposeInPlace();
        copy_mad2eigen2( b.size(), sol.data(), x.ptr());
    }

    /** \brief  Generalized real-symmetric or complex-Hermitian eigenproblem.

    This function uses the EIGEN3 GeneralizedSelfAdjointEigenSolver routine.

    A should be selfadjoint and B positive definide.

    \verbatim
    Specifies the problem type to be solved:
    = 1:  A*x = (lambda)*B*x
    = 2:  A*B*x = (lambda)*x
    = 3:  B*A*x = (lambda)*x
    \endverbatim

    */
    template <typename T>
    void sygv(const Tensor<T>& a, const Tensor<T>& B, int itype,
              Tensor<T>& V, Tensor< typename Tensor<T>::scalar_type >& e) {
        TENSOR_ASSERT(a.ndim() == 2, "sygv requires a matrix",a.ndim(),&a);
        TENSOR_ASSERT(a.dim(0) == a.dim(1), "sygv requires square matrix",0,&a);
        TENSOR_ASSERT(B.ndim() == 2, "sygv requires a matrix",B.ndim(),&a);
        TENSOR_ASSERT(B.dim(0) == B.dim(1), "sygv requires square matrix",0,&a);
        integer n = a.dim(0);

        Tensor<T> b = transpose(B);     // For Hermitian case
        V = transpose(a);               // For Hermitian case
        e = Tensor<typename Tensor<T>::scalar_type>(n);

        Eigen::Matrix<T, Dynamic, Dynamic> g(n,n);
        Tensor<T> aT=transpose(a);

        copy_mad2eigen2( a.size(), aT.ptr(), g.data());
        Eigen::Matrix<T, Dynamic, Dynamic> gT = g.transpose();

        Eigen::Matrix<T, Dynamic, Dynamic> h(n,n);
        Tensor<T> bT=transpose(b);

        copy_mad2eigen2( b.size(), bT.ptr(), h.data());
        Eigen::Matrix<T, Dynamic, Dynamic> hT = h.transpose();

        GeneralizedSelfAdjointEigenSolver<Eigen::Matrix<T, Dynamic, Dynamic>> sol(g, h);

        Eigen::Matrix<T, Dynamic, Dynamic> ev = sol.eigenvectors();

        V = Tensor<T>(n,n); //fool


        ev.transposeInPlace();
        copy_mad2eigen2( V.size(), ev.data(), V.ptr());
        copy_mad2eigen2( e.size(), sol.eigenvalues().data(), e.ptr());


    }
}
#endif //MADNESS_HAS_EIGEN3
