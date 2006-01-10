#ifndef MAD_MPIAR_H
#define MAD_MPIAR_H

#include <serialize/archive.h>
#include <misc/communicator.h>

namespace madness {
    namespace archive {
        
        class MPIOutputArchive : public BaseOutputArchive {
            Communicator& comm;
            ProcessID dest;
        public:
            MPIOutputArchive(Communicator& comm, const ProcessID& dest) : comm(comm), dest(dest) {};
            
            template <class T>
            inline 
            typename madness::enable_if< madness::is_fundamental<T>, void >::type
            store(const T* t, long n) {
                cout << "sending to " << dest << endl;
                comm.send((const char *) t, n*sizeof(T), dest, 0);
            };
        };
        
        
        class MPIInputArchive : public BaseInputArchive {
            Communicator& comm;
            ProcessID src;
        public:
            MPIInputArchive(Communicator& comm, const ProcessID& src) : comm(comm), src(src) {};

            template <class T>
            inline 
            typename madness::enable_if< madness::is_fundamental<T>, void >::type
            load(T* t, long n) {
                cout << "receiving from " << src << endl;
                comm.recv((char *) t, n*sizeof(T), src, 0);
            }
        };

        // No type checking over MPI stream for efficiency
        template <class T> 
        struct ArchivePrePostImpl<MPIOutputArchive,T> {
            static void preamble_store(MPIOutputArchive& ar) {};
            static inline void postamble_store(MPIOutputArchive& ar) {};
        };
        
        // No type checking over MPI stream for efficiency
        template <class T> 
        struct ArchivePrePostImpl<MPIInputArchive,T> {
            static inline void preamble_load(MPIInputArchive& ar) {};
            static inline void postamble_load(MPIInputArchive& ar) {};
        };
        
    }  
}
#endif
