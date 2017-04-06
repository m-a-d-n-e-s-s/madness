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

#ifndef MADNESS_TENSOR_ELEM_H__INCLUDED
#define MADNESS_TENSOR_ELEM_H__INCLUDED

#include <madness/madness_config.h>
#include <madness/world/MADworld.h>

#ifdef MADNESS_HAS_ELEMENTAL_EMBEDDED

#include <madness/tensor/tensor.h>
#include <madness/tensor/tensor_lapack.h>
#include <madness/tensor/distributed_matrix.h>
#include <madness/world/print.h>
#include <madness/world/binsorter.h>

#include <ctime>
#if defined(HAVE_EL_H)
# include <El.hpp>
namespace elem = El;
#elif defined(HAVE_ELEMENTAL_H)
# include <elemental.hpp>
#else
# error "MADNESS_HAS_ELEMENTAL set but neither HAVE_EL_H nor HAVE_ELEMENTAL_H set: file an issue at " MADNESS_PACKAGE_URL
#endif

namespace madness {

    namespace detail {

        template <typename T> 
        struct Value {
            int i;
            int j;
            T t;
            Value(int i, int j, T t) : i(i), j(j), t(t) {}
            Value(){} // required for use in STL container
            template <typename Archive>
            void serialize(const Archive& ar) {
                ar & i & j & t;
            }
        };

        template <typename T>  
        class MadToElemDistCopy {
            elem::DistMatrix<T>& d;
        public:
            MadToElemDistCopy(elem::DistMatrix<T>& d) : d(d) {}

            void operator()(const Value<T>& v) {
                d.Set(v.i, v.j, v.t);
            }
        };

        template <typename T>  
        class ElemToMadDistCopy {
            DistributedMatrix<T>& d;
        public:
            ElemToMadDistCopy(DistributedMatrix<T>& d) : d(d) {}

            void operator()(const Value<T>& v) {
                d.set(v.i, v.j, v.t);
            }
        };

    }


    /// Backport of more recent Elemental DistMatrix API
    template <typename T>
    ProcessID Owner(const elem::DistMatrix<T>& d, int i, int j) {
        int RowOwner = (i+d.ColAlign()) % d.ColStride(); // is Col/Row Align in later versions ... no ment
        int ColOwner = (j+d.RowAlign()) % d.RowStride(); 
        return RowOwner+ColOwner*d.ColStride();        
    }


    /// Copy a MADNESS distributed matrix into an Elemental distributed matrix

    /// Should work for any distribution of either
    template <typename T>
    void copy_to_elemental(const DistributedMatrix<T>& din, elem::DistMatrix<T>& dout) {
        BinSorter< detail::Value<T> , detail::MadToElemDistCopy<T> > s(din.get_world(), detail::MadToElemDistCopy<T>(dout));

        int64_t ilo, ihi, jlo, jhi;
        din.local_colrange(ilo,ihi);
        din.local_rowrange(jlo,jhi);
        const Tensor<T>& t = din.data();
        for (int64_t i=ilo; i<=ihi; i++) {
            for (int64_t j=jlo; j<=jhi; j++) {
                ProcessID owner = Owner(dout, i, j);
                s.insert(owner, detail::Value<T>(i,j,t(i-ilo, j-jlo)));
            }
        }

        s.finish();
    }

    /// Copy a MADNESS distributed matrix from an Elemental distributed matrix

    /// Should work for any distribution of either
    template <typename T>
    void copy_from_elemental(const elem::DistMatrix<T>& din, DistributedMatrix<T>& dout) {
        BinSorter< detail::Value<T> , detail::ElemToMadDistCopy<T> > s(dout.get_world(), detail::ElemToMadDistCopy<T>(dout));
        
        const int64_t colShift =    din.ColShift(); // first row we own
        const int64_t rowShift =    din.RowShift(); // first col we own
        const int64_t colStride =   din.ColStride();
        const int64_t rowStride =   din.RowStride();
        const int64_t localHeight = din.LocalHeight();
        const int64_t localWidth =  din.LocalWidth();
        
        for( int64_t jLocal=0; jLocal<localWidth; ++jLocal ) {
            for( int64_t iLocal=0; iLocal<localHeight; ++iLocal ) {
                const int64_t i = colShift + iLocal*colStride;
                const int64_t j = rowShift + jLocal*rowStride;

                const ProcessID owner = dout.owner(i,j);
                s.insert(owner, detail::Value<T>(i,j,din.GetLocal(iLocal,jLocal)));
            }
        }

        s.finish();
    }

    /** \brief  Generalized real-symmetric or complex-Hermitian eigenproblem.

    This function uses the Elemental HermitianGenDefiniteEig routine.

    A should be selfadjoint and B positive definite.

    \verbatim
    Specifies the problem type to be solved:
    = 1:  A*x = (lambda)*B*x
    = 2:  A*B*x = (lambda)*x (TODO)
    = 3:  B*A*x = (lambda)*x (TODO)
    \endverbatim

    */
    template <typename T>
    void sygv(const DistributedMatrix<T>& A, 
               const DistributedMatrix<T>& B, 
               int itype,
               DistributedMatrix<T>& X, 
               Tensor< typename Tensor<T>::scalar_type >& e) 
    {
        const int64_t n = A.coldim();
        const int64_t m = A.rowdim();
        MADNESS_ASSERT(n==m);
        MADNESS_ASSERT(n==B.coldim() && m==B.rowdim());
        
        const int blocksize = 128;
        const elem::Grid GG(A.get_world().mpi.comm().Get_mpi_comm() );
        elem::SetBlocksize(blocksize);
        
        elem::DistMatrix<T> EA(n,n,GG);
        elem::DistMatrix<T> EB(n,n,GG);
        
        copy_to_elemental(A,EA);
        copy_to_elemental(B,EB);
        
        elem::HermitianGenDefiniteEigType eigType = elem::AXBX;
        elem::UpperOrLower uplo = elem::CharToUpperOrLower('U');
        elem::DistMatrix<T> Xd(n,n,GG);
        elem::DistMatrix<T,elem::VR,elem::STAR> wd( n, n, GG);
        
        // 0.83+ ???
        elem::HermitianGenDefiniteEig(eigType, uplo, EA, EB, wd, Xd, elem::SortType::ASCENDING);
        
        // 0.79-0.82 ?
        //elem::HermitianGenDefiniteEig(eigType, uplo, EA, EB, wd, Xd);
        //elem::hermitian_eig::Sort(wd, Xd);
        
        A.get_world().mpi.Barrier();
        
        X = DistributedMatrix<T>(A.distribution());
        e = Tensor<typename Tensor<T>::scalar_type>(n);
        
        copy_from_elemental(Xd, X);
        
        const int  localHeight1 = wd.LocalHeight();
        const int colShift1 = wd.ColShift(); // first row we own
        const int colStride1= wd.ColStride();
        T * buffer = e.ptr();
        for( int iLocal=0; iLocal<localHeight1; ++iLocal ) {
            const int jLocal=0;
            const int i = colShift1 + iLocal*colStride1;
            //buffer[i]= wd.Get( iLocal, jLocal);
            buffer[i]= wd.GetLocal( iLocal, jLocal);
        }
        
        A.get_world().gop.sum(e.ptr(),n);
    }

    /** \brief  Generalized real-symmetric or complex-Hermitian eigenproblem.

    This function uses the Elemental HermitianGenDefiniteEig routine.

    A should be selfadjoint and B positive definite.

    \verbatim
    Specifies the problem type to be solved:
    = 1:  A*x = (lambda)*B*x
    = 2:  A*B*x = (lambda)*x (TODO)
    = 3:  B*A*x = (lambda)*x (TODO)
    \endverbatim

    */
    template <typename T>
    void sygvp(World& world,
               const Tensor<T>& a, const Tensor<T>& B, int itype,
              Tensor<T>& V, Tensor< typename Tensor<T>::scalar_type >& e) {
        TENSOR_ASSERT(a.ndim() == 2, "sygv requires a matrix",a.ndim(),&a);
        TENSOR_ASSERT(a.dim(0) == a.dim(1), "sygv requires square matrix",0,&a);
        TENSOR_ASSERT(B.ndim() == 2, "sygv requires a matrix",B.ndim(),&a);
        TENSOR_ASSERT(B.dim(0) == B.dim(1), "sygv requires square matrix",0,&a);
        TENSOR_ASSERT(a.iscontiguous(),"sygvp requires a contiguous matrix (a)",0,&a);
        TENSOR_ASSERT(B.iscontiguous(),"sygvp requires a contiguous matrix (B)",0,&B);

        world.gop.broadcast(a.ptr(), a.size(), 0);
        world.gop.broadcast(B.ptr(), B.size(), 0);

        const int n = a.dim(1);

        e = Tensor<typename Tensor<T>::scalar_type>(n);
        V = Tensor<T>(n,n);

        if (a.dim(0) <= 4) {
            // Work around bug in elemental/pmrrr/mkl/something for n=2,3
            sygv(a, B, itype, V, e);
            return;
        }

        world.gop.fence(); //<<<<<< Essential to quiesce MADNESS threads/comms

        const int blocksize = 128;
        try {
            const elem::Grid GG( world.mpi.comm().Get_mpi_comm() );
            elem::SetBlocksize(blocksize);

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

            world.mpi.Barrier();
     
            elem::HermitianGenDefiniteEigType eigType = elem::AXBX;
            //const char* uu="U";
            elem::UpperOrLower uplo = elem::CharToUpperOrLower('U');
            elem::DistMatrix<T> Xd( n, n, GG );
            elem::DistMatrix<T,elem::VR,elem::STAR> wd( n, n, GG);


            // 0.83+ ???
            elem::HermitianGenDefiniteEig( eigType, uplo, gd, hd, wd, Xd,elem::SortType::ASCENDING);

            // 0.79-0.82 ?
            //elem::HermitianGenDefiniteEig( eigType, uplo, gd, hd, wd, Xd);
            //elem::hermitian_eig::Sort( wd, Xd );

            world.mpi.Barrier();
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
            std::cerr << S << std::endl;
        }

        world.gop.fence(); //<<<<<< Essential to quiesce MADNESS threads/comms
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
    void gesvp(World& world, const Tensor<T>& a, const Tensor<T>& b, Tensor<T>& x) {

        TENSOR_ASSERT(a.ndim() == 2, "gesv requires matrix",a.ndim(),&a);

        int n = a.dim(0), m = a.dim(1), nrhs = b.dim(1);

        TENSOR_ASSERT(m == n, "gesv requires square matrix",0,&a);
        TENSOR_ASSERT(b.ndim() <= 2, "gesv require a vector or matrix for the RHS",b.ndim(),&b);
        TENSOR_ASSERT(a.dim(0) == b.dim(0), "gesv matrix and RHS must conform",b.ndim(),&b);
        TENSOR_ASSERT(a.iscontiguous(),"gesvp requires a contiguous matrix (a)",0,&a);
        TENSOR_ASSERT(b.iscontiguous(),"gesvp requires a contiguous matrix (b)",0,&b);
        world.gop.broadcast(a.ptr(), a.size(), 0);
        world.gop.broadcast(b.ptr(), b.size(), 0);

        Tensor<T> AT = transpose(a);

        world.gop.fence(); //<<<<<< Essential to quiesce MADNESS threads/comms

        int blocksize = 128;
        try {
            const elem::Grid GG( world.mpi.comm().Get_mpi_comm() );
            elem::SetBlocksize(blocksize);
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
            world.mpi.Barrier();
    
           elem::GaussianElimination(gd, hd);
    
            world.mpi.Barrier();
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
           std::cerr << S << std::endl;
       }
        world.gop.fence(); //<<<<<< Essential to quiesce MADNESS threads/comms
    }
}

#else

namespace madness {
    // sequential fall back code
    template <typename T>
    void sygvp(World& world,
               const Tensor<T>& a, const Tensor<T>& B, int itype,
               Tensor<T>& V, Tensor< typename Tensor<T>::scalar_type >& e) {
        sygv(a, B, itype, V, e);
	world.gop.broadcast_serializable(V,0);
	world.gop.broadcast_serializable(e,0);
    }
    
    // sequential fall back code
    template <typename T>
    void gesvp(World& world, const Tensor<T>& a, const Tensor<T>& b, Tensor<T>& x) {
        gesv(a, b, x);
    }
}

#endif //MADNESS_HAS_ELEMENTAL_EMBEDDED

#endif // MADNESS_TENSOR_ELEM_H__INCLUDED
