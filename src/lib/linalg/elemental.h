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


  $Id: eigen.h 2615 2011-10-23 13:24:06Z jeff.science@gmail.com $
*/


#ifdef MADNESS_HAS_ELEMENTAL



#include <mra/mra.h>

#include <tensor/tensor.h>
using madness::Tensor;


#include <mra/mra.h>
#include <iostream>
using std::cout;
using std::endl;

#include <algorithm>
using std::min;
using std::max;

#include "elemental.hpp"
using namespace elem;


#ifdef MADNESS_FORINT
typedef MADNESS_FORINT integer;
#else
typedef long integer;
#endif //MADNESS_FORINT

namespace madness {
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
        const integer n = a.dim(1);

        e = Tensor<typename Tensor<T>::scalar_type>(n);
        V = Tensor<T>(n,n);
//madness
        madness::World world(MPI::COMM_WORLD);
//elemental
        mpi::Comm comm = mpi::COMM_WORLD;
        try {
            const Grid GG( comm );
            elem::DistMatrix<T> gd( n, n, GG );
    
            const int colShift = gd.ColShift(); // first row we own
            const int rowShift = gd.RowShift(); // first col we own
            const int colStride =gd.ColStride();
            const int rowStride = gd.RowStride();
            const int localHeight = gd.LocalHeight();
            const int localWidth = gd.LocalWidth();
            {
               const T * buffer = a.ptr();
               for( int jLocal=0; jLocal<localWidth; ++jLocal )
               {
                   for( int iLocal=0; iLocal<localHeight; ++iLocal )
                   {
                         const int i = colShift + iLocal*colStride;
                         const int j = rowShift + jLocal*rowStride;
                         gd.SetLocal( iLocal, jLocal, buffer[i+j*n] );
                   }
               }
            }
            //gd.Print("gs");
            elem::DistMatrix<T> hd( n, n, GG );
            {
               const T * buffer = B.ptr();
               for( int jLocal=0; jLocal<localWidth; ++jLocal )
               {
                   for( int iLocal=0; iLocal<localHeight; ++iLocal )
                   {
                         const int i = colShift + iLocal*colStride;
                         const int j = rowShift + jLocal*rowStride;
                         hd.SetLocal( iLocal, jLocal, buffer[i+j*(n)] );
                   }
               }
            }
     
            mpi::Barrier( GG.Comm() );
     
            HermitianGenDefiniteEigType eigType = AXBX;
            char* uu="U";
            UpperOrLower uplo = CharToUpperOrLower( *uu);
            elem::DistMatrix<T> Xd( n, n, GG );
            elem::DistMatrix<T,VR,STAR> wd( n, n);
            HermitianGenDefiniteEig( eigType, uplo, gd, hd, wd, Xd );
     
            mpi::Barrier( GG.Comm() );
            //Xd.Print("Xs");
     
     //retrive eigenvalues
            {
               const int  localHeight1 = wd.LocalHeight();
               const int colShift1 = wd.ColShift(); // first row we own
               const int colStride1 =wd.ColStride();
               T * buffer = e.ptr();
               for( int iLocal=0; iLocal<localHeight1; ++iLocal )
               {
                   const int jLocal=0;
                   const int i = colShift1 + iLocal*colStride1;
                   buffer[i]= wd.GetLocal( iLocal, jLocal);
               }
            }
            world.gop.sum(e.ptr(),n);
     //retrive eigenvectors
            {
               T * buffer = V.ptr();
               for( int jLocal=0; jLocal<localWidth; ++jLocal )
               {
                  for( int iLocal=0; iLocal<localHeight; ++iLocal )
                  {
                     const int i = colShift + iLocal*colStride;
                     const int j = rowShift + jLocal*rowStride;
                     buffer[i+j*n]= Xd.GetLocal( iLocal, jLocal);
                  }
               }
            }
            world.gop.sum(V.ptr(), n*n);
            V=transpose(V);
        }
        catch (TensorException S) {
           cout << "Caught a tensor exception in sygv elemental\n";
           cout << S;
        }
    }
}
#endif //MADNESS_HAS_ELEMENTAL
