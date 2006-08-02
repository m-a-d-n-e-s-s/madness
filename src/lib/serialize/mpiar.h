#ifndef MAD_MPIAR_H
#define MAD_MPIAR_H

#include <serialize/archive.h>
#include <serialize/vecar.h>
#include <misc/communicator.h>

namespace madness {
    namespace archive {

        class MPIRawOutputArchive : public BaseOutputArchive {
            mutable Communicator& comm;
            ProcessID dest;
            int tag;
        public:
            MPIRawOutputArchive(Communicator& comm, const ProcessID& dest, int tag=5447468)
                    : comm(comm), dest(dest), tag(tag) {};

            template <class T>
            inline
            typename madness::enable_if< madness::is_fundamental<T>, void >::type
            store(const T* t, long n) const {
                comm.Send(t, n, dest, tag);
            };
        };

        class MPIRawInputArchive : public BaseInputArchive {
            mutable Communicator& comm;
            ProcessID src;
            int tag;
        public:
            MPIRawInputArchive(Communicator& comm, const ProcessID& src, int tag=5447468)
                    : comm(comm), src(src), tag(tag) {};

            template <class T>
            inline
            typename madness::enable_if< madness::is_fundamental<T>, void >::type
            load(T* t, long n) const {
                comm.Recv(t, n, src, tag);
            }
        };


        class MPIOutputArchive : public BaseOutputArchive {
            mutable Communicator& comm;
            ProcessID dest;
            int tag;
	    const std::size_t bufsize;
	    mutable std::vector<unsigned char> v;
	    madness::archive::VectorOutputArchive var;
        public:
            MPIOutputArchive(Communicator& comm, const ProcessID& dest, int tag=5447468)
                    : comm(comm), dest(dest), tag(tag), bufsize(1024*1024), v(), var(v) 
		{v.reserve(2*bufsize);};

            template <class T>
            inline
            typename madness::enable_if< madness::is_fundamental<T>, void >::type
            store(const T* t, long n) const {
		var.store(t,n);
		if (v.size() > bufsize) flush();
            };

	    void flush() const {
		if (v.size()) {
		    comm.Send(v.size(), dest, tag);
		    comm.Send(&v[0], v.size(), dest, tag);
		    v.clear();
		    if (v.capacity() == 0) v.reserve(2*bufsize);
		}
	    };

	    void close() {flush();};

	    ~MPIOutputArchive() {close();};
        };

        class MPIInputArchive : public BaseInputArchive {
            mutable Communicator& comm;
            ProcessID src;
            int tag;
	    mutable std::vector<unsigned char> v;
	    madness::archive::VectorInputArchive var;
        public:
            MPIInputArchive(Communicator& comm, const ProcessID& src, int tag=5447468)
                    : comm(comm), src(src), tag(tag), v(), var(v) {};

            template <class T>
            inline
            typename madness::enable_if< madness::is_fundamental<T>, void >::type
            load(T* t, long n) const {
		if (var.nbyte_avail()) {
		    var.load(t,n);
		}
		else {
		    std::size_t m;
		    comm.Recv(m, src, tag);
		    var.rewind();
		    v.clear();
		    if (v.capacity() < m) v.reserve(m);
		    comm.Recv(&v[0], m, src, tag);
		}		    
            }
        };

        // No type checking over MPI stream for efficiency
        template <class T>
        struct ArchivePrePostImpl<MPIRawOutputArchive,T> {
            static void preamble_store(const MPIRawOutputArchive& ar) {};
            static inline void postamble_store(const MPIRawOutputArchive& ar) {};
        };

        // No type checking over MPI stream for efficiency
        template <class T>
        struct ArchivePrePostImpl<MPIRawInputArchive,T> {
            static inline void preamble_load(const MPIRawInputArchive& ar) {};
            static inline void postamble_load(const MPIRawInputArchive& ar) {};
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
