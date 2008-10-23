#ifndef WORLD_SAFE_MPI_H
#define WORLD_SAFE_MPI_H

/// \file safempi.h
/// \brief Serializes calls to MPI in case it does not support THREAD_MULTIPLE

#include <mpi.h>
#include <world/typestuff.h>
#include <world/worldthread.h>

namespace SafeMPI {

#define SERIALIZE_MPI

#ifdef SERIALIZE_MPI
    extern madness::Mutex charon;      // Inside safempi.cc
#define GLOBAL_MUTEX madness::ScopedMutex obolus(charon)
#else
#define GLOBAL_MUTEX
#endif

    class Request : private MPI::Request{
    public:
        Request() : MPI::Request() {}

        Request(const MPI::Request& request) : MPI::Request(request) {}

        bool Test() {
            GLOBAL_MUTEX;
            return MPI::Request::Test();
        }

        static int Testsome(int n, SafeMPI::Request* request, int* ind, MPI::Status* status) {
            GLOBAL_MUTEX;
            return MPI::Request::Testsome(n, static_cast<MPI::Request*>(request), ind, status);
        }

        static int Testsome(int n, MPI::Request* request, int* ind, MPI::Status* status) {
            GLOBAL_MUTEX;
            return MPI::Request::Testsome(n, request, ind, status);
        }
    };

    class Intracomm {
        MPI::Intracomm& comm;
        int me;
        int numproc;
        int mpi_tag_ub;

    public:
        Intracomm(MPI::Intracomm& comm) : comm(comm)
        {
            GLOBAL_MUTEX;
            me = comm.Get_rank();
            numproc = comm.Get_size();

            long value; // Must be 64-bit on 64-bit machines
            comm.Get_attr(MPI::TAG_UB, &value);
            mpi_tag_ub = value;
        }
        
        int Get_rank() const {
            return me;
        }

        int Get_size() const {
            return numproc;
        }

        MPI::Intracomm& Get_comm() const {
            return comm;
        }

        ::SafeMPI::Request Isend(const void* buf, size_t count, const MPI::Datatype& datatype, int dest, int tag) const {
            GLOBAL_MUTEX;
            return comm.Isend(buf,count,datatype,dest,tag);
        }

        ::SafeMPI::Request Irecv(void* buf, size_t count, const MPI::Datatype& datatype, int src, int tag) const {
            GLOBAL_MUTEX;
            return comm.Irecv(buf, count, datatype, src, tag);
        }

        void Send(const void* buf, size_t count, const MPI::Datatype& datatype, int dest, int tag) const {
            GLOBAL_MUTEX;
            comm.Send(buf,count,datatype,dest,tag);
        }

        void Recv(void* buf, int count, const MPI::Datatype& datatype, int source, int tag, MPI::Status& status) const {
            GLOBAL_MUTEX;
            comm.Recv(buf,count,datatype,source,tag,status);
        }

        void Recv(void* buf, int count, const MPI::Datatype& datatype, int source, int tag) const {
            GLOBAL_MUTEX;
            comm.Recv(buf,count,datatype,source,tag);
        }

        void Bcast(void* buf, size_t count, const MPI::Datatype& datatype, int root) const {
            GLOBAL_MUTEX;
            return comm.Bcast(buf, count, datatype, root);
        }

        void Reduce(void* sendbuf, void* recvbuf, int count, const MPI::Datatype& datatype, const MPI::Op& op, int root) const {
            GLOBAL_MUTEX;
            comm.Reduce(sendbuf, recvbuf, count, datatype, op, root);
        };

        void Allreduce(void* sendbuf, void* recvbuf, int count, const MPI::Datatype& datatype, const MPI::Op& op) const {
            GLOBAL_MUTEX;
            comm.Allreduce(sendbuf, recvbuf, count, datatype, op);
        };
        void Get_attr(int key, void* value) const {
            GLOBAL_MUTEX;
            comm.Get_attr(key, value);
        }

        void Abort(int code=1) const {
            comm.Abort(code);
        }


        /// Returns a unique tag for general use (tag>1023)

        /// Unique is slightly optimistic.  The method simply
        /// increments/wraps a counter and returns the next legal
        /// value. For most MPIs with maximum tag values equal to
        /// 2^31-1 this is adequate.  For MPIs such as LAM which only
        /// provide the bare minimum of 32768 tags you have a greater
        /// chance of collision.
        int unique_tag() {
            GLOBAL_MUTEX;
            static volatile int tag = 1024;
            int result = tag++;
            if (tag > mpi_tag_ub) tag = 1024;
            return result;
        }

        /// Returns a unique tag reserved for system use (tag<1023)

        /// Note that the MPI RMI server thread reserves tag 1023 so 
        /// user tags begin above that number.
        int unique_reserved_tag() {
            GLOBAL_MUTEX;
            static volatile int tag = 1;
            int result = tag++;
            if (tag > 1022) tag = 1;
            return result;
        }

        // Below here are extensions to MPI but have to be here since overload resolution
        // is not applied across different class scopes ... hence even having Isend with
        // a different signature in a derived class (WorldMpiInterface) hides this interface.
        // Thus, they must all be here.
        //
        // !! All of the routines below call the protected interfaces provided above.
        // !! Please ensure any additional routines follow this convention.

        /// Isend one element ... disabled for pointers to reduce accidental misuse.
        template <class T>
        typename madness::enable_if_c< !madness::is_pointer<T>::value, SafeMPI::Request>::type
        Isend(const T& datum, int dest, int tag=1) const {
            return Isend(&datum, sizeof(T), MPI::BYTE, dest, tag);
        }

        /// Async receive data of up to lenbuf elements from process dest
        template <class T>
        SafeMPI::Request
        Irecv(T* buf, int count, int source, int tag) const {
            return Irecv(buf, count*sizeof(T), MPI::BYTE, source, tag);
        }
        
        
        /// Async receive datum from process dest with default tag=1
        template <class T>
        typename madness::enable_if_c< !madness::is_pointer<T>::value, SafeMPI::Request>::type
        Irecv(T& buf, int source, int tag=1) const {
            return Irecv(&buf, sizeof(T), MPI::BYTE, source, tag);
        }
        

        /// Send array of lenbuf elements to process dest 
        template <class T>
        void Send(const T* buf, long lenbuf, int dest, int tag) const {
            Send((void* )buf, lenbuf*sizeof(T), MPI::BYTE, dest, tag);
        }
     

        /// Send element to process dest with default tag=1
        
        /// Disabled for pointers to reduce accidental misuse.
        template <class T>
        typename madness::enable_if_c< !madness::is_pointer<T>::value, void>::type
        Send(const T& datum, int dest, int tag=1) const {
            Send((void* )&datum, sizeof(T), MPI::BYTE, dest, tag);
        }
     

        /// Receive data of up to lenbuf elements from process dest
        template <class T>
            void
            Recv(T* buf, long lenbuf, int src, int tag) const {
            Recv(buf, lenbuf*sizeof(T), MPI::BYTE, src, tag);
        }
     
        /// Receive data of up to lenbuf elements from process dest with status
        template <class T>
            void
            Recv(T* buf, long lenbuf, int src, int tag, MPI::Status& status) const {
            Recv(buf, lenbuf*sizeof(T), MPI::BYTE, src, tag, status);
        }
     
     
        /// Receive datum from process src with default tag=1
        template <class T>
        typename madness::enable_if_c< !madness::is_pointer<T>::value, void>::type
        Recv(T& buf, int src, int tag=1) const {
            Recv(&buf, sizeof(T), MPI::BYTE, src, tag);
        }


        /// MPI broadcast an array of count elements
        
        /// NB.  Read documentation about interaction of MPI collectives and AM/task handling.
        template <class T>
        void Bcast(T* buffer, int count, int root) const {
            Bcast(buffer,count*sizeof(T),MPI::BYTE,root);
        }
        
        
        /// MPI broadcast a datum
        
        /// NB.  Read documentation about interaction of MPI collectives and AM/task handling.
        template <class T>
        typename madness::enable_if_c< !madness::is_pointer<T>::value, void>::type
        Bcast(T& buffer, int root) const {
            Bcast(&buffer, sizeof(T), MPI::BYTE,root);
        }

        int rank() const {return Get_rank();}
        
        int nproc() const {return Get_size();}
        
        int size() const {return Get_size();}


        /// Construct info about a binary tree with given root
        
        /// Constructs a binary tree spanning the communicator with
        /// process root as the root of the tree.  Returns the logical
        /// parent and children in the tree of the calling process.  If
        /// there is no parent/child the value -1 will be set.
        void binary_tree_info(int root, int& parent, int& child0, int& child1) {
            int np = nproc();
            int me = (rank()+np-root)%np;   // Renumber processes so root has me=0
            parent = (((me-1)>>1)+root)%np;        // Parent in binary tree
            child0 = (me<<1)+1+root;        // Left child
            child1 = (me<<1)+2+root;        // Right child
            if (child0 >= np && child0<(np+root)) child0 -= np;
            if (child1 >= np && child1<(np+root)) child1 -= np;
            
            if (me == 0) parent = -1;
            if (child0 >= np) child0 = -1;
            if (child1 >= np) child1 = -1;
        }

    };
#undef GLOBAL_MUTEX
}

#endif
