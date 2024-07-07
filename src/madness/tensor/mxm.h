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

#ifndef MADNESS_TENSOR_MXM_H__INCLUDED
#define MADNESS_TENSOR_MXM_H__INCLUDED

#include <madness/madness_config.h>

// This just to check if config is actually working
//#ifndef HAVE_MTXMQ
//#error "MTXMQ missing"
//#endif

#define HAVE_FAST_BLAS
#ifdef  HAVE_FAST_BLAS
//#ifdef HAVE_INTEL_MKL
#include <madness/tensor/cblas.h>
#endif

/// \file tensor/mxm.h
/// \brief Internal use only

// This file is ONLY included into tensor.cc ... separated here just
// to shrink file size.  Don't try to include anywhere else

// Due to both flakey compilers and performance concerns,
// we use a simple reference implementation of the mxm
// routines for all except T=double.


namespace madness {

    // Start with reference implementations.  Then provide optimized implementations, falling back to reference if not available on specific platforms

    /// Matrix \c += Matrix * matrix reference implementation (slow but correct)
    template <typename T, typename Q, typename S>
    static inline void mxm_reference(long dimi, long dimj, long dimk,
                                     T* MADNESS_RESTRICT c, const Q* MADNESS_RESTRICT a,
                                     const S* MADNESS_RESTRICT b) {
        /*
          c(i,j) = c(i,j) + sum(k) a(i,k)*b(k,j)
          
          where it is assumed that the last index in each array is has unit
          stride and the dimensions are as provided.
        */
        
        for (long i=0; i<dimi; ++i) {
            for (long k=0; k<dimk; ++k) {
                for (long j=0; j<dimj; ++j) {
                    c[i*dimj+j] += a[i*dimk+k]*b[k*dimj+j];
                }
            }
        }
    }
    

    /// Matrix \c += Matrix transpose * matrix ... reference implementation (slow but correct)
    template <typename T, typename Q, typename S>
    static inline
    void mTxm_reference(long dimi, long dimj, long dimk,
                        T* MADNESS_RESTRICT c, const Q* MADNESS_RESTRICT a,
                        const S* MADNESS_RESTRICT b) {
        /*
          c(i,j) = c(i,j) + sum(k) a(k,i)*b(k,j)
          
          where it is assumed that the last index in each array is has unit
          stride and the dimensions are as provided.
          
          i loop might be long in anticipated application
        */
        
        for (long k=0; k<dimk; ++k) {
            for (long j=0; j<dimj; ++j) {
                for (long i=0; i<dimi; ++i) {
                    c[i*dimj+j] += a[k*dimi+i]*b[k*dimj+j];
                }
            }
        }
    }

    /// Matrix \c += Matrix * matrix transpose ... reference implementation (slow but correct)
    template <typename T, typename Q, typename S>
    static inline void mxmT_reference (long dimi, long dimj, long dimk,
                                       T* MADNESS_RESTRICT c, const Q* MADNESS_RESTRICT a,
                                       const S* MADNESS_RESTRICT b) {
        /*
          c(i,j) = c(i,j) + sum(k) a(i,k)*b(j,k)
          
          where it is assumed that the last index in each array is has unit
          stride and the dimensions are as provided.
          
          i loop might be long in anticipated application
        */
        
        for (long i=0; i<dimi; ++i) {
            for (long j=0; j<dimj; ++j) {
                T sum = 0;
                for (long k=0; k<dimk; ++k) {
                    sum += a[i*dimk+k]*b[j*dimk+k];
                }
                c[i*dimj+j] += sum;
            }
        }
    }
    
    /// Matrix \c += Matrix transpose * matrix transpose reference implementation (slow but correct)
    template <typename T, typename Q, typename S>
    static inline void mTxmT_reference(long dimi, long dimj, long dimk,
                                       T* MADNESS_RESTRICT c, const Q* MADNESS_RESTRICT a,
                                       const S* MADNESS_RESTRICT b) {
        /*
          c(i,j) = c(i,j) + sum(k) a(k,i)*b(j,k)
          
          where it is assumed that the last index in each array is has unit
          stride and the dimensions are as provided.
        */
        
        for (long i=0; i<dimi; ++i) {
            for (long j=0; j<dimj; ++j) {
                for (long k=0; k<dimk; ++k) {
                    c[i*dimj+j] += a[k*dimi+i]*b[j*dimk+k];
                }
            }
        }
    }

    /// Matrix = Matrix transpose * matrix ... slow reference implementation
    
    /// This routine does \c C=AT*B whereas mTxm does C=C+AT*B.
    /// \code
    ///    c(i,j) = sum(k) a(k,i)*b(k,j)  <------ does not accumulate into C
    /// \endcode
    ///
    /// \c ldb is the last dimension of b in C storage (the leading dimension
    /// in fortran storage).  It is here to accomodate multiplying by a matrix
    /// stored with \c ldb>dimj which happens in madness when transforming with
    /// low rank matrices.  A matrix in dense storage has \c ldb=dimj which is
    /// the default for backward compatibility.
    template <typename aT, typename bT, typename cT>
    void mTxmq_reference(long dimi, long dimj, long dimk,
                         cT* MADNESS_RESTRICT c, const aT* a, const bT* b, long ldb=-1) {
        if (ldb == -1) ldb=dimj;
        MADNESS_ASSERT(ldb>=dimj);
        //std::cout << "IN GENERIC mTxmq " << tensor_type_names[TensorTypeData<aT>::id] << " " << tensor_type_names[TensorTypeData<bT>::id] << " " << tensor_type_names[TensorTypeData<cT>::id] << "\n";
        for (long i=0; i<dimi; ++i,c+=dimj,++a) {
            for (long j=0; j<dimj; ++j) c[j] = 0.0;
            const aT *aik_ptr = a;
            for (long k=0; k<dimk; ++k,aik_ptr+=dimi) {
                aT aki = *aik_ptr;
                for (long j=0; j<dimj; ++j) {
                    c[j] += aki*b[k*ldb+j];
                }
            }
        }
    }
    

#if defined(HAVE_FAST_BLAS) && !defined(HAVE_INTEL_MKL)
    // MKL provides support for mixed real/complex operations but most other libraries do not
    
    /// Matrix += Matrix * matrix ... BLAS/MKL interface version
    
    /// Does \c C=C+A*B 
    /// \code
    ///    c(i,j) = c(i,j) + sum(k) a(i,k)*b(k,j)
    /// \endcode
    template <typename T>
    void mxm(long dimi, long dimj, long dimk,
              T* MADNESS_RESTRICT c, const T* a, const T* b) {
        const T one = 1.0;  // alpha in *gemm
        cblas::gemm(cblas::NoTrans,cblas::NoTrans,dimj,dimi,dimk,one,b,dimj,a,dimk,one,c,dimj);
    }

    /// Matrix += Matrix transpose * matrix ... MKL interface version
    
    /// Does \c C=C+AT*B 
    /// \code
    ///    c(i,j) = c(i,j) + sum(k) a(k,i)*b(k,j)
    /// \endcode
    template <typename T>
    void mTxm(long dimi, long dimj, long dimk,
              T* MADNESS_RESTRICT c, const T* a, const T* b) {
        const T one = 1.0;  // alpha in *gemm
        cblas::gemm(cblas::NoTrans,cblas::Trans,dimj,dimi,dimk,one,b,dimj,a,dimi,one,c,dimj);
    }

    /// Matrix += Matrix * matrix transpose ... MKL interface version
    
    /// Does \c C=C+A*BT
    /// \code
    ///    c(i,j) = c(i,j) + sum(k) a(i,k)*b(j,k)
    /// \endcode
    template <typename T>
    void mxmT(long dimi, long dimj, long dimk,
              T* MADNESS_RESTRICT c, const T* a, const T* b) {
        const T one = 1.0;  // alpha in *gemm
        cblas::gemm(cblas::Trans,cblas::NoTrans,dimj,dimi,dimk,one,b,dimk,a,dimk,one,c,dimj);
    }
    
    /// Matrix += Matrix transpose * matrix transpose ... MKL interface version
    
    /// Does \c C=C+AT*BT
    /// \code
    ///    c(i,j) = c(i,j) + sum(k) a(k,i)*b(j,k)
    /// \endcode
    template <typename T>
    void mTxmT(long dimi, long dimj, long dimk,
               T* MADNESS_RESTRICT c, const T* a, const T* b) {
        const T one = 1.0;  // alpha in *gemm
        cblas::gemm(cblas::Trans,cblas::Trans,dimj,dimi,dimk,one,b,dimk,a,dimi,one,c,dimj);
    }

    /// Matrix = Matrix transpose * matrix ... MKL interface version
    
    /// Does \c C=AT*B whereas mTxm does C=C+AT*B.  
    /// \code
    ///    c(i,j) = sum(k) a(k,i)*b(k,j)  <------ does not accumulate into C
    /// \endcode
    ///
    /// \c ldb is the last dimension of b in C storage (the leading dimension
    /// in fortran storage).  It is here to accomodate multiplying by a matrix
    /// stored with \c ldb>dimj which happens in madness when transforming with
    /// low rank matrices.  A matrix in dense storage has \c ldb=dimj which is
    /// the default for backward compatibility.
    template <typename T>
    void mTxmq(long dimi, long dimj, long dimk,
               T* MADNESS_RESTRICT c, const T* a, const T* b, long ldb=-1) {
        if (ldb == -1) ldb=dimj;
        MADNESS_ASSERT(ldb>=dimj);

        if (dimi==0 || dimj==0) return; // nothing to do and *GEMM will complain
        if (dimk==0) {
            for (long i=0; i<dimi*dimj; i++) c[i] = 0.0;
        }
        
        const T one = 1.0;  // alpha in *gemm
        const T zero = 0.0; // beta  in *gemm
        cblas::gemm(cblas::NoTrans,cblas::Trans,dimj,dimi,dimk,one,b,ldb,a,dimi,zero,c,dimj);
    }  

#ifdef HAVE_MTXMQ
    template <>
    void mTxmq(long dimi, long dimj, long dimk, double* MADNESS_RESTRICT c, const double* a, const double* b, long ldb);

    // Bootstrap complex*real from real*real
    template <typename T>
    void mTxmq(long dimi, long dimj, long dimk, std::complex<T>* MADNESS_RESTRICT c, const std::complex<T>* a, const T* b, long ldb) {
      T* Rc = new T[dimi*dimj];
      T* Ic = new T[dimi*dimj];
      T* Ra = new T[dimi*dimk];
      T* Ia = new T[dimi*dimk];

      for (long i=0; i<dimi*dimk; i++) {
	Ra[i] = a[i].real();
	Ia[i] = a[i].imag();
      }
      mTxmq(dimi,dimj,dimk,Rc,Ra,b,ldb);
      mTxmq(dimi,dimj,dimk,Ic,Ia,b,ldb);
      for (long i=0; i<dimi*dimj; i++) c[i] = std::complex<T>(Rc[i],Ic[i]);
      
      delete[] Rc;
      delete[] Ic;
      delete[] Ra;
      delete[] Ia;
    }  
  
#endif

#endif
    
#ifdef HAVE_INTEL_MKL
    /// Matrix += Matrix * matrix ... MKL interface version
    
    /// Does \c C=C+A*B 
    /// \code
    ///    c(i,j) = c(i,j) + sum(k) a(i,k)*b(k,j)
    /// \endcode
    template <typename aT, typename bT, typename cT>
    void mxm(long dimi, long dimj, long dimk,
              cT* MADNESS_RESTRICT c, const aT* a, const bT* b) {
        const cT one = 1.0;  // alpha in *gemm
        cblas::gemm(cblas::NoTrans,cblas::NoTrans,dimj,dimi,dimk,one,b,dimj,a,dimk,one,c,dimj);
    }

    /// Matrix += Matrix transpose * matrix ... MKL interface version
    
    /// Does \c C=C+AT*B 
    /// \code
    ///    c(i,j) = c(i,j) + sum(k) a(k,i)*b(k,j)
    /// \endcode
    template <typename aT, typename bT, typename cT>
    void mTxm(long dimi, long dimj, long dimk,
              cT* MADNESS_RESTRICT c, const aT* a, const bT* b) {
        const cT one = 1.0;  // alpha in *gemm
        cblas::gemm(cblas::NoTrans,cblas::Trans,dimj,dimi,dimk,one,b,dimj,a,dimi,one,c,dimj);
    }

    /// Matrix += Matrix * matrix transpose ... MKL interface version
    
    /// Does \c C=C+A*BT
    /// \code
    ///    c(i,j) = c(i,j) + sum(k) a(i,k)*b(j,k)
    /// \endcode
    template <typename aT, typename bT, typename cT>
    void mxmT(long dimi, long dimj, long dimk,
              cT* MADNESS_RESTRICT c, const aT* a, const bT* b) {
        const cT one = 1.0;  // alpha in *gemm
        cblas::gemm(cblas::Trans,cblas::NoTrans,dimj,dimi,dimk,one,b,dimk,a,dimk,one,c,dimj);
    }
    
    /// Matrix += Matrix transpose * matrix transpose ... MKL interface version
    
    /// Does \c C=C+AT*BT
    /// \code
    ///    c(i,j) = c(i,j) + sum(k) a(k,i)*b(j,k)
    /// \endcode
    template <typename aT, typename bT, typename cT>
    void mTxmT(long dimi, long dimj, long dimk,
               cT* MADNESS_RESTRICT c, const aT* a, const bT* b) {
        const cT one = 1.0;  // alpha in *gemm
        cblas::gemm(cblas::Trans,cblas::Trans,dimj,dimi,dimk,one,b,dimk,a,dimi,one,c,dimj);
    }

    /// Matrix = Matrix transpose * matrix ... MKL interface version
    
    /// Does \c C=AT*B whereas mTxm does C=C+AT*B.  
    /// \code
    ///    c(i,j) = sum(k) a(k,i)*b(k,j)  <------ does not accumulate into C
    /// \endcode
    ///
    /// \c ldb is the last dimension of b in C storage (the leading dimension
    /// in fortran storage).  It is here to accomodate multiplying by a matrix
    /// stored with \c ldb>dimj which happens in madness when transforming with
    /// low rank matrices.  A matrix in dense storage has \c ldb=dimj which is
    /// the default for backward compatibility.
    template <typename aT, typename bT, typename cT>
    void mTxmq(long dimi, long dimj, long dimk,
               cT* MADNESS_RESTRICT c, const aT* a, const bT* b, long ldb=-1) {
        if (ldb == -1) ldb=dimj;
        MADNESS_ASSERT(ldb>=dimj);

        if (dimi==0 || dimj==0) return; // nothing to do and *GEMM will complain
        if (dimk==0) {
            for (long i=0; i<dimi*dimj; i++) c[i] = 0.0;
        }
        
        const cT one = 1.0;  // alpha in *gemm
        const cT zero = 0.0; // beta  in *gemm
        cblas::gemm(cblas::NoTrans,cblas::Trans,dimj,dimi,dimk,one,b,ldb,a,dimi,zero,c,dimj);
    }

#ifdef HAVE_MTXMQ
template <>
void mTxmq(long dimi, long dimj, long dimk, double* MADNESS_RESTRICT c, const double* a, const double* b, long ldb);
#endif
    
#else

    // Fall back to reference implementations

    template <typename T, typename Q, typename S>
    static inline void mxm(long dimi, long dimj, long dimk,
                           T* MADNESS_RESTRICT c, const Q* MADNESS_RESTRICT a,
                           const S* MADNESS_RESTRICT b) {
        mxm_reference(dimi, dimj, dimk, c, a, b);
    }
    
    template <typename T, typename Q, typename S>
    static inline
    void mTxm(long dimi, long dimj, long dimk,
              T* MADNESS_RESTRICT c, const Q* MADNESS_RESTRICT a,
              const S* MADNESS_RESTRICT b) {
        mTxm_reference(dimi, dimj, dimk, c, a, b);
    }

    template <typename T, typename Q, typename S>
    static inline void mxmT(long dimi, long dimj, long dimk,
                            T* MADNESS_RESTRICT c, const Q* MADNESS_RESTRICT a,
                            const S* MADNESS_RESTRICT b) {
        mxmT_reference(dimi, dimj, dimk, c, a, b);
    }
    
    template <typename T, typename Q, typename S>
    static inline void mTxmT(long dimi, long dimj, long dimk,
                             T* MADNESS_RESTRICT c, const Q* MADNESS_RESTRICT a,
                             const S* MADNESS_RESTRICT b) {
        mTxmT_reference(dimi, dimj, dimk, c, a, b);
    }

    template <typename aT, typename bT, typename cT>
    void mTxmq(long dimi, long dimj, long dimk,
               cT* MADNESS_RESTRICT c, const aT* a, const bT* b, long ldb=-1) {
        mTxmq_reference(dimi, dimj, dimk, c, a, b, ldb);
    }

    // The following are restricted to double only
    
    /// Matrix transpose * matrix (hand unrolled version)
    
    template <>
    inline void mTxm(long dimi, long dimj, long dimk,
                     double* MADNESS_RESTRICT c, const double* MADNESS_RESTRICT a,
                     const double* MADNESS_RESTRICT b) {
        /*
          c(i,j) = c(i,j) + sum(k) a(k,i)*b(k,j)  <--- NOTE ACCUMULATION INTO C
          
          where it is assumed that the last index in each array is has unit
          stride and the dimensions are as provided.
          
          i loop might be long in anticipated application
          
          4-way unrolled k loop ... empirically fastest on PIII
          compared to 2/3 way unrolling (though not by much).
        */
        
        long dimk4 = (dimk/4)*4;
        for (long i=0; i<dimi; ++i,c+=dimj) {
            const double* ai = a+i;
            const double* p = b;
            for (long k=0; k<dimk4; k+=4,ai+=4*dimi,p+=4*dimj) {
                double ak0i = ai[0   ];
                double ak1i = ai[dimi];
                double ak2i = ai[dimi+dimi];
                double ak3i = ai[dimi+dimi+dimi];
                const double* bk0 = p;
                const double* bk1 = p+dimj;
                const double* bk2 = p+dimj+dimj;
                const double* bk3 = p+dimj+dimj+dimj;
                for (long j=0; j<dimj; ++j) {
                    c[j] += ak0i*bk0[j] + ak1i*bk1[j] + ak2i*bk2[j] + ak3i*bk3[j];
                }
            }
            for (long k=dimk4; k<dimk; ++k) {
                double aki = a[k*dimi+i];
                const double* bk = b+k*dimj;
                for (long j=0; j<dimj; ++j) {
                    c[j] += aki*bk[j];
                }
            }
        }
    }
    
    
    
    /// Matrix * matrix transpose (hand unrolled version)
    
    template <>
    inline void mxmT(long dimi, long dimj, long dimk,
                     double* MADNESS_RESTRICT c,
                     const double* MADNESS_RESTRICT a, const double* MADNESS_RESTRICT b) {
        /*
          c(i,j) = c(i,j) + sum(k) a(i,k)*b(j,k)
          
          where it is assumed that the last index in each array is has unit
          stride and the dimensions are as provided.
          
          j loop might be long in anticipated application
          
          Unrolled i loop.  Empirically fastest on PIII compared
          to unrolling j, or both i&j.
        */
        
        long dimi2 = (dimi/2)*2;
        for (long i=0; i<dimi2; i+=2) {
            const double* ai0 = a+i*dimk;
            const double* ai1 = a+i*dimk+dimk;
            double* MADNESS_RESTRICT ci0 = c+i*dimj;
            double* MADNESS_RESTRICT ci1 = c+i*dimj+dimj;
            for (long j=0; j<dimj; ++j) {
                double sum0 = 0;
                double sum1 = 0;
                const double* bj = b + j*dimk;
                for (long k=0; k<dimk; ++k) {
                    sum0 += ai0[k]*bj[k];
                    sum1 += ai1[k]*bj[k];
                }
                ci0[j] += sum0;
                ci1[j] += sum1;
            }
        }
        for (long i=dimi2; i<dimi; ++i) {
            const double* ai = a+i*dimk;
            double* MADNESS_RESTRICT ci = c+i*dimj;
            for (long j=0; j<dimj; ++j) {
                double sum = 0;
                const double* bj = b+j*dimk;
                for (long k=0; k<dimk; ++k) {
                    sum += ai[k]*bj[k];
                }
                ci[j] += sum;
            }
        }
    }
    
    /// Matrix * matrix (hand unrolled version)
    template <>
    inline void mxm(long dimi, long dimj, long dimk,
                    double* MADNESS_RESTRICT c, const double* MADNESS_RESTRICT a, const double* MADNESS_RESTRICT b) {
        /*
          c(i,j) = c(i,j) + sum(k) a(i,k)*b(k,j)
          
          where it is assumed that the last index in each array is has unit
          stride and the dimensions are as provided.
          
          4-way unrolled k loop ... empirically fastest on PIII
          compared to 2/3 way unrolling (though not by much).
        */
        
        long dimk4 = (dimk/4)*4;
        for (long i=0; i<dimi; ++i, c+=dimj,a+=dimk) {
            const double* p = b;
            for (long k=0; k<dimk4; k+=4,p+=4*dimj) {
                double aik0 = a[k  ];
                double aik1 = a[k+1];
                double aik2 = a[k+2];
                double aik3 = a[k+3];
                const double* bk0 = p;
                const double* bk1 = bk0+dimj;
                const double* bk2 = bk1+dimj;
                const double* bk3 = bk2+dimj;
                for (long j=0; j<dimj; ++j) {
                    c[j] += aik0*bk0[j] + aik1*bk1[j] + aik2*bk2[j] + aik3*bk3[j];
                }
            }
            for (long k=dimk4; k<dimk; ++k) {
                double aik = a[k];
                for (long j=0; j<dimj; ++j) {
                    c[j] += aik*b[k*dimj+j];
                }
            }
        }
    }
    
    /// Matrix transpose * matrix transpose (hand tiled and unrolled)
    template <>
    inline void mTxmT(long dimi, long dimj, long dimk,
                      double* MADNESS_RESTRICT csave, const double* MADNESS_RESTRICT asave, const double* MADNESS_RESTRICT b) {
        /*
          c(i,j) = c(i,j) + sum(k) a(k,i)*b(j,k)
          
          where it is assumed that the last index in each array is has unit
          stride and the dimensions are as provided.
          
          Tiled k, copy row of a into temporary, and unroll j once.
        */
        
        const int ktile=32;
        double ai[ktile];
        long dimj2 = (dimj/2)*2;
        
        for (long klo=0; klo<dimk; klo+=ktile, asave+=ktile*dimi, b+=ktile) {
            long khi = klo+ktile;
            if (khi > dimk) khi = dimk;
            long nk = khi-klo;
            
            const double * MADNESS_RESTRICT a = asave;
            double * MADNESS_RESTRICT c = csave;
            for (long i=0; i<dimi; ++i,c+=dimj,++a) {
                const double* q = a;
                for (long k=0; k<nk; ++k,q+=dimi) ai[k] = *q;
                
                const double* bj0 = b;
                for (long j=0; j<dimj2; j+=2,bj0+=2*dimk) {
                    const double* bj1 = bj0+dimk;
                    double sum0 = 0;
                    double sum1 = 0;
                    for (long k=0; k<nk; ++k) {
                        sum0 += ai[k]*bj0[k];
                        sum1 += ai[k]*bj1[k];
                    }
                    c[j  ] += sum0;
                    c[j+1] += sum1;
                }
                
                for (long j=dimj2; j<dimj; ++j,bj0+=dimk) {
                    double sum = 0;
                    for (long k=0; k<nk; ++k) {
                        sum += ai[k]*bj0[k];
                    }
                    c[j] += sum;
                }
            }
        }
    }

    /*
     * mtxm, but with padded buffers.
     *
     * ext_b is the extent of the b array, so shrink() isn't needed.
     */
    template <typename aT, typename bT, typename cT>
    void mTxmq_padding(long dimi, long dimj, long dimk, long ext_b,
               cT* c, const aT* a, const bT* b) {
        const int alignment = 4;
        bool free_b = false;
        long effj = dimj;

        /* Setup a buffer for c if needed */
        cT* c_buf = c;
        if (dimj%alignment) {
            effj = (dimj | 3) + 1;
            c_buf = (cT*)malloc(sizeof(cT)*dimi*effj);
        }

        /* Copy b into a buffer if needed */
        if (ext_b%alignment) {
            free_b = true;
            bT* b_buf = (bT*)malloc(sizeof(bT)*dimk*effj);

            bT* bp = b_buf;
            for (long k=0; k<dimk; k++, bp += effj, b += ext_b)
                memcpy(bp, b, sizeof(bT)*dimj);

            b = b_buf;
            ext_b = effj;
        }

        cT* c_work = c_buf;
        /* mTxm */
        for (long i=0; i<dimi; ++i,c_work+=effj,++a) {
            for (long j=0; j<dimj; ++j) c_work[j] = 0.0;
            const aT *aik_ptr = a;
            for (long k=0; k<dimk; ++k,aik_ptr+=dimi) {
                aT aki = *aik_ptr;
                for (long j=0; j<dimj; ++j) {
                    c_work[j] += aki*b[k*ext_b+j];
                }
            }
        }

        /* Copy c out if needed */
        if (dimj%alignment) {
            cT* ct = c_buf;
            for (long i=0; i<dimi; i++, ct += effj, c += dimj)
                memcpy(c, ct, sizeof(cT)*dimj);

            free(c_buf);
        }

        /* Free the buffer for b */
        if (free_b) free((bT*)b);
    }

#ifdef HAVE_IBMBGQ
    extern void bgq_mtxmq_padded(long ni, long nj, long nk, long ej, 
            double* c, const double* a, const double* b);
    extern void bgq_mtxmq_padded(long ni, long nj, long nk, long ej, 
            __complex__ double* c, const __complex__ double* a, const __complex__ double* b);
    extern void bgq_mtxmq_padded(long ni, long nj, long nk, long ej, 
            __complex__ double* c, const double* a, const __complex__ double* b);
    extern void bgq_mtxmq_padded(long ni, long nj, long nk, long ej, 
            __complex__ double* c, const __complex__ double* a, const double* b);

    template <>
        inline void mTxmq_padding(long ni, long nj, long nk, long ej, 
                double* c, const double* a, const double* b) {
            bgq_mtxmq_padded(ni, nj, nk, ej, c, a, b);
        }

    template <>
        inline void mTxmq_padding(long ni, long nj, long nk, long ej, 
                __complex__ double* c, const __complex__ double* a, const __complex__ double* b) {
            bgq_mtxmq_padded(ni, nj, nk, ej, c, a, b);
        }

    template <>
        inline void mTxmq_padding(long ni, long nj, long nk, long ej, 
                __complex__ double* c, const double* a, const __complex__ double* b) {
            bgq_mtxmq_padded(ni, nj, nk, ej, c, a, b);
        }

    template <>
        inline void mTxmq_padding(long ni, long nj, long nk, long ej, 
                __complex__ double* c, const __complex__ double* a, const double* b) {
            bgq_mtxmq_padded(ni, nj, nk, ej, c, a, b);
        }
#elif defined(HAVE_IBMBGP)
    extern void bgpmTxmq(long ni, long nj, long nk, double* MADNESS_RESTRICT c, 
                         const double* a, const double* b);
    extern void bgpmTxmq(long ni, long nj, long nk, double_complex* MADNESS_RESTRICT c, 
                         const double_complex* a, const double_complex* b);
 
    template <>
    inline void mTxmq(long ni, long nj, long nk, double* MADNESS_RESTRICT c, const double* a, const double* b) {
        bgpmTxmq(ni, nj, nk, c, a, b);
    }

    template <>
    inline void mTxmq(long ni, long nj, long nk, double_complex* MADNESS_RESTRICT c, const double_complex* a, const double_complex* b) {
        bgpmTxmq(ni, nj, nk, c, a, b);
    }

// #elif defined(X86_64) && !defined(DISABLE_SSE3)
//     template <>
//     void mTxmq(long dimi, long dimj, long dimk,
//                double* MADNESS_RESTRICT c, const double* a, const double* b);

//     template <>
//     void mTxmq(long dimi, long dimj, long dimk,
//                double_complex* MADNESS_RESTRICT c, const double_complex* a, const double_complex* b);

// #ifndef __INTEL_COMPILER
//     template <>
//     void mTxmq(long dimi, long dimj, long dimk,
//                double_complex* MADNESS_RESTRICT c, const double_complex* a, const double* b);
// #endif

// #elif defined(X86_32)
//     template <>
//     void mTxmq(long dimi, long dimj, long dimk,
//                double* MADNESS_RESTRICT c, const double* a, const double* b);
#endif // HAVE_IBMBGQ

#endif // HAVE_INTEL_MKL
    
}    
#endif // MADNESS_TENSOR_MXM_H__INCLUDED
