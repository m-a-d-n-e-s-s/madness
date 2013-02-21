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


#include <iostream>


#include <tensor/tensor.h>
#include <world/print.h>
#include <iostream>
#include <linalg/tensor_lapack.h>


#include <mra/mra.h>
#include <mra/lbdeux.h>


#include <Eigen/Dense>


#include <ctime>
#include "elemental.hpp"
using namespace elem;


//1#include <linalg/clapack.h>

using namespace madness;
using namespace std;
using namespace Eigen;
using namespace SafeMPI;


#ifndef SCALAR
// #define SCALAR std::complex<float>
#define SCALAR double
#endif
#ifndef REALX
// #define SCALAR std::complex<float>
#define REALX real8
#endif

#ifndef STATIC
#define STATIC static
#endif



//vama11#include <iostream>
//vama11using std::cout;
//vama11using std::endl;
//vama11
//vama11
//vama11#include <tensor/tensor.h>
//vama11#include <world/print.h>
//vama11#include <algorithm>
//vama11#include <linalg/tensor_lapack.h>
//vama11
//vama11#include <mra/mra.h>
//vama11#include <mra/lbdeux.h>
//vama11
//vama11#ifdef MADNESS_HAS_EIGEN3
//vama11# include <Eigen/Dense>
//vama11  using namespace Eigen;
//vama11#endif
//vama11
//vama11#include <ctime>
//vama11#include "elemental.hpp"
//vama11using namespace elem;
//vama11
//vama11//vama#include <mpi.h>
//vama11//vamausing namespace SafeMPI;
//vama11
//vama11using namespace madness;
//vama11using namespace std;
//vama11using namespace SafeMPI;
//vama11
//vama11
//vama11
//vama11
//vama11#ifdef MADNESS_FORINT
//vama11typedef MADNESS_FORINT integer;
//vama11#else
//vama11typedef long integer;
//vama11#endif //MADNESS_FORINT




namespace madness {
    /** \brief  Generalized real-symmetric or complex-Hermitian eigenproblem.

    This function uses the Elemental HermitianGenDefiniteEig routine.

    A should be selfadjoint and B positive definide.

    \verbatim
    Specifies the problem type to be solved:
    = 1:  A*x = (lambda)*B*x
    = 2:  A*B*x = (lambda)*x (TODO)
    = 3:  B*A*x = (lambda)*x (TODO)
    \endverbatim

    */
    template <typename T>
    void sygvp(const Tensor<T>& a, const Tensor<T>& B, int itype,
              Tensor<T>& V, Tensor< typename Tensor<T>::scalar_type >& e) {
        TENSOR_ASSERT(a.ndim() == 2, "sygv requires a matrix",a.ndim(),&a);
        TENSOR_ASSERT(a.dim(0) == a.dim(1), "sygv requires square matrix",0,&a);
        TENSOR_ASSERT(B.ndim() == 2, "sygv requires a matrix",B.ndim(),&a);
        TENSOR_ASSERT(B.dim(0) == B.dim(1), "sygv requires square matrix",0,&a);
        const int n = a.dim(1);

        e = Tensor<typename Tensor<T>::scalar_type>(n);
        V = Tensor<T>(n,n);
//elemental
        //mpi::Comm comm = mpi::COMM_WORLD;
        mpi::Comm comm = COMM_WORLD.comm;
        //mpi::Comm comm = mpi::COMM_WORLD;
//madness
        World world(SafeMPI::COMM_WORLD);

        //const int myrank = mpi::CommRank( COMM_WORLD.comm );

        int blocksize = 128;
        elem::InitBlocksize(blocksize, comm);

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
                         //gd.Set( iLocal, jLocal, buffer[i+j*n] );
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
                         //hd.Set( iLocal, jLocal, buffer[i+j*(n)] );
                         hd.SetLocal( iLocal, jLocal, buffer[i+j*(n)] );
                   }
               }
            }
           // hd.Print("hd");

            mpi::Barrier( GG.Comm() );
     
            elem::HermitianGenDefiniteEigType eigType = AXBX;
            char* uu="U";
            elem::UpperOrLower uplo = CharToUpperOrLower( *uu);
            elem::DistMatrix<T> Xd( n, n, GG );
            elem::DistMatrix<T,VR,STAR> wd( n, n, GG);
            elem::HermitianGenDefiniteEig( eigType, uplo, gd, hd, wd, Xd );
            elem::SortEig( wd, Xd );
     
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
                   //buffer[i]= wd.Get( iLocal, jLocal);
                   buffer[i]= wd.GetLocal( iLocal, jLocal);
               }
            }
            //world.gop.broadcast(e.ptr(),e.size(), 0);
            world.gop.sum(e.ptr(),n);
            //if(myrank ==0) cout<< e << endl;
     //retrive eigenvectors
            {
               T * buffer = V.ptr();
               for( int jLocal=0; jLocal<localWidth; ++jLocal )
               {
                  for( int iLocal=0; iLocal<localHeight; ++iLocal )
                  {
                     const int i = colShift + iLocal*colStride;
                     const int j = rowShift + jLocal*rowStride;
                     //buffer[i+j*n]= Xd.Get( iLocal, jLocal);
                     buffer[i+j*n]= Xd.GetLocal( iLocal, jLocal);
                  }
               }
            }
            world.gop.sum(V.ptr(), n*n);
            V=madness::transpose(V);
            //world.gop.broadcast(V.ptr(),V.size(), 0);
            //if(myrank ==0)cout<< V << endl;
        }
        catch (TensorException S) {
           cout << "Caught a tensor exception in sygv elemental\n";
           cout << S;
        }
        elem::FinBlocksize();
    }

    /** \brief  Solve Ax = b for general A using the Elemental. 
    The solution is computed through (partially pivoted) Gaussian elimination.

    A should be a square matrix (float, double, float_complex,
    double_complex) and b should be either a vector, or a matrix with
    each vector stored in a column (i.e., b[n,nrhs]).

    It will solve Ax=b as written.

    b can be a vector or a matrix, the only restriction is that satisfies b.rows()==A.rows()

    */
    template <typename T>
    void gesvp(const Tensor<T>& a, const Tensor<T>& b, Tensor<T>& x) {

        TENSOR_ASSERT(a.ndim() == 2, "gesv requires matrix",a.ndim(),&a);

        int n = a.dim(0), m = a.dim(1), nrhs = b.dim(1);

        TENSOR_ASSERT(m == n, "gesv requires square matrix",0,&a);
        TENSOR_ASSERT(b.ndim() <= 2, "gesv require a vector or matrix for the RHS",b.ndim(),&b);
        TENSOR_ASSERT(a.dim(0) == b.dim(0), "gesv matrix and RHS must conform",b.ndim(),&b);

        Tensor<T> AT = transpose(a);

//elemental
        //mpi::Comm comm = mpi::COMM_WORLD;
        mpi::Comm comm = COMM_WORLD.comm;
//madness
        World world(SafeMPI::COMM_WORLD);
        //World world(mpi::COMM_WORLD);

        //const int myrank = mpi::CommRank( COMM_WORLD.comm );
        int blocksize = 128;
        InitBlocksize(blocksize, comm);

        try {
            const elem::Grid GG( comm );
            elem::DistMatrix<T> gd( n, n, GG );

    
            {
                 const int colShift = gd.ColShift(); // 1st row local
                 const int rowShift = gd.RowShift(); // 1st col local
                 const int colStride =gd.ColStride();
                 const int rowStride = gd.RowStride();
                 const int localHeight = gd.LocalHeight();
                 const int localWidth = gd.LocalWidth();
                 {
                    const T * buffer = AT.ptr();
                    for( int jLocal=0; jLocal<localWidth; ++jLocal )
                    {
                        for( int iLocal=0; iLocal<localHeight; ++iLocal )
                        {
                              const int i = colShift + iLocal*colStride;
                              const int j = rowShift + jLocal*rowStride;
                              //gd.Set( iLocal, jLocal, buffer[i+j*n] );
                              gd.SetLocal( iLocal, jLocal, buffer[i+j*n] );
                        }
                    }
                }
            }
            Tensor<T> bT;
            if (nrhs == 1) {
                 x = Tensor<T>(n);
                 bT = Tensor<T>(n);
                 bT = copy(b);
                 x = copy(b); //creating the correct size
            }
            else {
                 x = Tensor<T>(n,nrhs);
                 bT =  transpose(b);
            }

           for(int i=0; i< x.size(); ++i) {
            T * buffer = x.ptr() ;
            buffer[i] = 0.0;
            }

           //cout << "Caught a tensor exception \n";
           //cout << bT <<endl;
            elem::DistMatrix<T> hd( n, nrhs, GG );
            {
                 const int colShift = hd.ColShift(); // 1st row local
                 const int rowShift = hd.RowShift(); // 1st col local
                 const int colStride =hd.ColStride();
                 const int rowStride = hd.RowStride();
                 const int localHeight = hd.LocalHeight();
                 const int localWidth = hd.LocalWidth();
                 {
                     const T * buffer = bT.ptr();
                     for( int jLocal=0; jLocal<localWidth; ++jLocal )
                     {
                         for( int iLocal=0; iLocal<localHeight; ++iLocal )
                         {
                               const int i = colShift + iLocal*colStride;
                               const int j = rowShift + jLocal*rowStride;
                               //hd.Set( iLocal, jLocal, buffer[i+j*(n)]);
                               hd.SetLocal( iLocal, jLocal, buffer[i+j*(n)]);
        
                         }
                     }
                }
           }
           mpi::Barrier( GG.Comm() );
    
           elem::GaussianElimination(gd, hd);
    
           mpi::Barrier( GG.Comm() );
           {
                const int colShift = hd.ColShift(); // 1st row local
                const int rowShift = hd.RowShift(); // 1st col local
                const int colStride =hd.ColStride();
                const int rowStride = hd.RowStride();
                const int localHeight = hd.LocalHeight();
                const int localWidth = hd.LocalWidth();
      
                T * buffer = x.ptr();
                for( int jLocal=0; jLocal<localWidth; ++jLocal )
                {
                    for( int iLocal=0; iLocal<localHeight; ++iLocal )
                    {
                        const int i = colShift + iLocal*colStride;
                        const int j = rowShift + jLocal*rowStride;
                        //buffer[j+i*nrhs]= hd.Get( iLocal, jLocal);
                        buffer[j+i*nrhs]= hd.GetLocal( iLocal, jLocal);
                    }
                }
           }
           world.gop.sum(x.ptr(), n*nrhs);
            //world.gop.broadcast(x.ptr(),x.size(), 0);
            //if(myrank ==0) cout<< x << endl;
           
       }
       catch (TensorException S) {
           cout << "Caught a tensor exception in gesv elemental\n";
           cout << S;
       }
       elem::FinBlocksize();
    }
}
#endif //MADNESS_HAS_ELEMENTAL
