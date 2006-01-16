#ifndef MAD_MPIAR_H
#define MAD_MPIAR_H

#include <serialize/archive.h>
#include <misc/communicator.h>

namespace madness {
    namespace archive {
        
        class MPIOutputArchive : public BaseOutputArchive {
            mutable Communicator& comm;
            ProcessID dest;
        public:
            MPIOutputArchive(Communicator& comm, const ProcessID& dest) : comm(comm), dest(dest) {};
            
            template <class T>
            inline 
            typename madness::enable_if< madness::is_fundamental<T>, void >::type
            store(const T* t, long n) const {
                comm.Send(t, n, dest, 5447468);
            };
        };
        
        
        class MPIInputArchive : public BaseInputArchive {
            mutable Communicator& comm;
            ProcessID src;
        public:
            MPIInputArchive(Communicator& comm, const ProcessID& src) : comm(comm), src(src) {};

            template <class T>
            inline 
            typename madness::enable_if< madness::is_fundamental<T>, void >::type
            load(T* t, long n) const {
                comm.Recv(t, n, src, 5447468);
            }
        };

        // No type checking over MPI stream for efficiency
        template <class T> 
        struct ArchivePrePostImpl<MPIOutputArchive,T> {
            static void preamble_store(const MPIOutputArchive& ar) {};
            static inline void postamble_store(const MPIOutputArchive& ar) {};
        };
        
        // No type checking over MPI stream for efficiency
        template <class T> 
        struct ArchivePrePostImpl<MPIInputArchive,T> {
            static inline void preamble_load(const MPIInputArchive& ar) {};
            static inline void postamble_load(const MPIInputArchive& ar) {};
        };
        
    }  
}
#endif
