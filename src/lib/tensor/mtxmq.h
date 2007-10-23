#ifndef MTXMQ_H
#define MTXMQ_H

#include <madness_config.h>

namespace madness {
    
    /// Matrix = Matrix transpose * matrix ... reference implementation
    
    /// Does \c C=AT*B whereas mTxm does C=C+AT*B.  It also supposed
    /// to be fast which it achieves thru restrictions
    ///   * All dimensions even
    ///   * All pointers aligned
    /// \code
    ///    c(i,j) = sum(k) a(k,i)*b(k,j)  <------ does not accumulate into C
    /// \endcode
    template <typename aT, typename bT, typename cT>
    void mTxmq(long dimi, long dimj, long dimk,
               cT* RESTRICT c, const aT* a, const bT* b) {
        for (long i=0; i<dimi; i++,c+=dimj,a++) {
            for (long j=0; j<dimj; j++) c[j] = 0.0;
            const aT *aik_ptr = a;
            for (long k=0; k<dimk; k++,aik_ptr+=dimi) {
                aT aki = *aik_ptr;
                for (long j=0; j<dimj; j++) {
                    c[j] += aki*b[k*dimj+j];
                }
            }
        }
      
    }

#if defined(X86_32) || defined(x86_64)
    template <>
    void mTxmq(long dimi, long dimj, long dimk,
               double* RESTRICT c, const double* a, const double* b);
#endif
}

#endif
