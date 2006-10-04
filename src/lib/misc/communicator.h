#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

/// \file communicator.h
/// \brief Defines Communicator (interprocess communication + topology)

#include <iostream>
#include <unistd.h>
#include <vector>
#include <complex>
#include <cstdlib>
#include <unistd.h>

#include <madness_config.h>
#include <mad_types.h>
#include <typestuff.h>
#include <misc/print.h>
#include <misc/mpitags.h>
#include <misc/madexcept.h>


#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
#include <mpi.h>

namespace madness {
    
    // The brain dead std::mem_fun won't take a function with reference arguments,
    // hence this crappy little binder
    template <class T, typename argT, typename resultT> 
    class BindMemFun {
    private:
        T* t;
        resultT (T::*op)(argT);
    public:
        BindMemFun(T* t, resultT (T::*op)(argT)) : t(t), op(op) {};
        resultT operator()(argT arg) {return (t->*op)(arg);};
    };
    
    // Partial specialization of crappy binder for void return
    template <class T, typename argT> 
    class BindMemFun<T,argT,void> {
    private:
        T* t;
        void (T::*op)(argT);
    public:
        BindMemFun(T* t, void (T::*op)(argT)) : t(t), op(op) {};
        void operator()(argT arg) {(t->*op)(arg);};
    };
    
    // Corresponding crappy factory function for above crappy binder
    template <class T, typename argT, typename resultT>
    inline BindMemFun<T,argT,resultT> bind_mem_fun(T* t, resultT (T::*op)(argT)) {
        return BindMemFun<T,argT,resultT>(t,op);
    };
    
    static inline void noop(){};
    
    class Communicator;

    extern Communicator* comm_default;

    void xterm_debug(const Communicator& comm, const char* path=0, const char* display=0);

    /// Holds arguments sent via active messages ... deliberately small.
    class AMArg {
    private:
        typedef unsigned long ulong;
        static const ulong bad=0xffffffff;
        friend class Communicator;
        mutable ulong function;   ///< Handle to AM function
        
    public:
        ulong arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11;
        
        AMArg() : function(bad) {};
        AMArg(ulong arg0) : function(bad), arg0(arg0) {};
        AMArg(ulong arg0, ulong arg1) : function(bad), arg0(arg0), arg1(arg1) {};
        AMArg(ulong arg0, ulong arg1, ulong arg2) : function(bad), arg0(arg0), arg1(arg1), arg2(arg2) {};
        AMArg(ulong arg0, ulong arg1, ulong arg2, ulong arg3) : function(bad), arg0(arg0), arg1(arg1), arg2(arg2), arg3(arg3) {};
        AMArg(ulong arg0, ulong arg1, ulong arg2, ulong arg3, ulong arg4) : function(bad), arg0(arg0), arg1(arg1), arg2(arg2), arg3(arg3), arg4(arg4) {};
        AMArg(ulong arg0, ulong arg1, ulong arg2, ulong arg3, ulong arg4, ulong arg5) : function(bad), arg0(arg0), arg1(arg1), arg2(arg2), arg3(arg3), arg4(arg4), arg5(arg5) {};
        AMArg(ulong arg0, ulong arg1, ulong arg2, ulong arg3, ulong arg4, ulong arg5, ulong arg6) : function(bad), arg0(arg0), arg1(arg1), arg2(arg2), arg3(arg3), arg4(arg4), arg5(arg5), arg6(arg6) {};
        AMArg(const void* p, int nbyte) {
            copyin(p,nbyte);
        };
        void copyin(const void* p, unsigned int nbyte) {
            MADNESS_ASSERT(nbyte <= sizeof(ulong)*12);
            memcpy(&arg0,p,nbyte); 
        };
        void copyout(void* p, unsigned int nbyte) const {
            MADNESS_ASSERT(nbyte <= sizeof(ulong)*12);
            memcpy(p,&arg0,nbyte); 
        };        
    };

    template <typename T> MPI::Datatype MPITypeFromType();
    template<> inline MPI::Datatype MPITypeFromType<int>() {return MPI::INT;};
    template<> inline MPI::Datatype MPITypeFromType<unsigned int>() {return MPI::UNSIGNED;};
    template<> inline MPI::Datatype MPITypeFromType<long>() {return MPI::LONG;};
    template<> inline MPI::Datatype MPITypeFromType< long long >() {return MPI::LONG;};
    template<> inline MPI::Datatype MPITypeFromType<unsigned long>() {return MPI::UNSIGNED_LONG;};
    template<> inline MPI::Datatype MPITypeFromType<double>() {return MPI::DOUBLE;};
    template<> inline MPI::Datatype MPITypeFromType< std::complex<double> >() {return MPI_COMPLEX;};

    typedef void (*am_handlerT)(Communicator&, ProcessID, const AMArg&);
    void am_ndiff_handler(Communicator& comm, ProcessID proc, const AMArg& arg);
    
    
    /// Manages mapping between general stuff (usually pointers) and handles
    
    /// Hash table with linear probling.  Lightweight structure intended for 
    /// infrequent writes and fast reads of a small amount of data.
    template <typename ptrT>
    class HandleManager {
    private:
        typedef std::pair<long,ptrT> pairT;
        const int size; //  must be a power of 2
        const unsigned long mask ;
        std::vector<pairT> v;  ///< For quick hashing of pointer to handle
        std::vector<ptrT> h;   ///< Map from handle to pointer

        
        inline int hash_ptr(ptrT p) const {
            // Returns lowest masked bits of 8-byte aligned pointer
            unsigned long i = (unsigned long) p;
            return (int) ((i>>3) & mask);            
        };
        
    public:
        HandleManager() : size(32768), mask(size-1), v(size), h() {
            h.reserve(size);
            for (int i=0; i<size; i++) v[i] = pairT(-1,0);  // -1 handle indicates empty
        };
        
        /// Registers the pointer and returns the handle
        long insert(ptrT pointer) {
            long handle = h.size();
            MADNESS_ASSERT(handle < (size/8));  // In part efficiency, also awareness
            h.push_back(pointer);
            int i = hash_ptr(pointer);
            while(v[i].first != -1) i = (i+1)&mask;
            v[i]= pairT(handle,pointer);
            return handle;        
        };
        
        /// Map handle to pointer, throws MadnessException on failure
        ptrT get_pointer(long handle) const {
            if (handle<0 || handle>=(int) h.size() || h[handle]<0)
                MADNESS_EXCEPTION("HandleManager: invalid handle",handle);
            return h[handle];
        };
        
        /// Map pointer to handle, throws MadnessException on failure
        long get_handle(ptrT pointer) const {
            int i = hash_ptr(pointer);
            while(v[i].second != pointer) {
                if (v[i].first < 0) 
                    MADNESS_EXCEPTION("HandleManager: pointer is not registered",0);
                i = (i+1)&mask;
            }
            return v[i].first;
        };
        
        void dump() const {
            madness::print("Registered Handlers");
            for (long i=0; i<(long) h.size(); i++) {
                madness::print("    ", i,
                               "--->", (void *) get_pointer(i),
                               "--->", get_handle(get_pointer(i)));
            }
        };
    };
        
    // Forward decl
    class TaskQueue;
    void task_add_am(am_handlerT op, ProcessID src, const AMArg& arg);
    
    
    /// Wraps an MPI communicator and provides topology, AM routines, etc.

    /// This class wraps an MPI communicator and provides info about
    /// process rank, number of processes along with basic send/recv
    /// capabilities.
    ///
    /// This class now exactly replicates a subset of the MPI C++
    /// communicator and provides additional functionality.
    ///
    /// Also managed by this class is the logical layout of processes as a
    /// 3D mesh (currently P=2^n only).  This might be better separated
    /// from the communicator, but usually if you want one you want the
    /// other.
    ///
    /// Now includes a crude "active" message interface.  In the
    /// absence of preemptive scheduling, interrupts, or threads, we
    /// have to use polling.  This is portable, but there are big
    /// implications w.r.t. performance (poll often, but not too
    /// often) and correctness.
    ///
    /// First, place calls to \c comm.am_poll() before and after any
    /// expensive (milli+ second) operation.   But, ...
    ///
    /// How to write a correct program using the AM interface?
    ///
    /// The following only applies to code segments that are using
    /// the AM interface.  However, since MADNESS is ultimately
    /// a library we may have to assume that all parts of the code
    /// must adhere to the following advice.
    ///
    /// The active messages are in the same data stream as all of
    /// other messages.   An example illustrates the consequences
    /// of this.  In the folllowing, another piece of code registered
    /// the handler \c action().
    /// \code
    /// void action(Communicator& comm, const AMArg& arg, ProcessID from) {
    ///    int reply = arg.arg0+1;
    ///    comm.send(reply,from);
    /// }
    /// 
    /// !!! THIS IS NOW OUT OF DATE !!!!
    ///
    /// int test(Communicator& comm, int handle, ProcessID p) {
    ///    int reply;
    ///    comm.am_send(AMArg(handle,0),p);
    ///    comm.recv(reply,p);
    ///    return reply;
    /// }
    /// \endcode
    /// The intent of the above code is to provide and use a service
    /// to which you can send an integer and receive a reply
    /// containing the incremented value.  Sadly, it is \em incorrect
    /// and will occasionally \em hang.
    ///
    /// The problem above is that if two processes send each other an
    /// AM request and then immediately enter the \c recv(), the AM
    /// request will never be processed.  Shoving an \c am_poll()
    /// between the \c am_send() and \c recv() might appear to work
    /// sometimes, but it is not sufficient, and the problem is
    /// deeper.
    ///
    /// While in a block of code that requires AM requests continue to
    /// be processed, a process must \em not block on any operation
    /// that depends upon another process for completion.  Otherwise,
    /// there is no easy guarantee that deadlock cannot occur.  To be
    /// completely safe, and for best performance, all recieves
    /// (whether associated with an active message, or otherwise)
    /// should be posted in advance of the expected message arrival.
    /// This will avoid buffer full problems and ensure active messages
    /// and their replies will be visible.
    ///
    /// A correct version of the above example is (the AM handler need
    /// not change).
    /// \code
    /// int test(Communicator& comm, int handle, ProcessID p) {
    ///    int reply;
    ///    MPI::Request req = comm.irecv(reply,p);
    ///    comm.am_send(AMArg(handle,0),p);
    ///    comm.am_wait(req);
    ///    return reply;
    /// }
    /// \endcode
    /// This is rather cumbersome and since request/response is
    /// a common motif you can also write this.
    /// \code
    /// int test(Communicator& comm) {
    ///    int reply;
    ///    comm.am_send_recv(AMArg(handle,0),p,&reply,sizeof(reply),p,tag);
    ///    return reply;
    /// }
    /// \endcode
    class Communicator {
        friend class TaskQueue;
    private:
        long _nproc;                ///< total no. of processes
        long _npx, _npy, _npz;      ///< #processes in each dimension
        long _mx, _my, _mz;         ///< mx=log2(npx)
        long _px, _py, _pz;         ///< coords of this process in mesh
        ProcessID _rank;            ///< rank of this process
        bool debug;                 ///< if true print send/receive traces
        int _unique_tag_next;       ///< used for dynamic unique tag generation


        // On the SGI Altix the MPI copy constructor for communicators
        // does not seem to work (even using Dup) so we have to store
        // a reference.  This is not satisfactory and we'll no doubt
        // return to this topic.
        mutable MPI::Intracomm& _comm;

        AMArg _am_arg;              // Used for async recv of AM
        AMArg _am_send_arg;         // Used for async send of AM
        MPI::Request _am_req;       // Request for async recv in server
        MPI::Request _am_send_req;  // Request for async send to remote server
        static const int _am_tag = AM_TAG;
        HandleManager<am_handlerT> _am_handle_manager;
        long _am_processing;    ///< If non-zero then am_messages are blocked
        volatile long _am_nsent;         ///< No. of AM messages sent
        volatile long _am_nrecv;         ///< No. of AM messages received
        bool _am_send_active;   ///< True if _am_arg is still in use by _am_send_req
        
        /// Given p=2^n processes make as close as possible to a cubic grid.

        /// Each dimension will be a power of 2 with npx >= npy >= npz
        void setup() {
//            _nproc = _comm.Get_size();
//            _rank = _comm.Get_rank();
            _nproc = MPI::COMM_WORLD.Get_size();
            _rank = MPI::COMM_WORLD.Get_rank();
            debug = false;
            _unique_tag_next = 2048;

            // Register am_ndiff then post AM receive buffer
            _am_processing = _am_nsent = _am_nrecv = 0;
            am_register(am_ndiff_handler);
            _am_req = Irecv(_am_arg, MPI::ANY_SOURCE, _am_tag);
            _am_send_active = false;

            _npz = 1;
            _mz = 0;
            while (8*_npz*_npz*_npz <= _nproc) {
                _npz *= 2;
                _mz++;
            }
            _npy = _npz;
            _my = _mz;
            while (4*_npy*_npy*_npz <= _nproc) {
                _npy *= 2;
                _my++;
            }
            _npx = _npy;
            _mx = _my;
            while (2*_npx*_npy*_npz <= _nproc) {
                _npx *= 2;
                _mx++;
            }

            if (_npx*_npy*_npz != _nproc) throw "nproc not a power of 2";

            // Determine coords of calling process in process mesh
            // Relies on npx >= npy >= npz.
            long pyz = _rank/_npx;
            _px = _rank - pyz*_npx;
            _py = pyz - (pyz/_npy)*_npy;
            _pz = (pyz - _py)/_npy;
            if (_px<0 || _py<0 || _pz<0 || _px>=_npx || _py>=_npy || _pz>=_npz ||
                    rank_from_coords(_px,_py,_pz)!=_rank)
                throw "p3_coords failed";
        };
        
        // packing is such that zero bits above function handle give immediate AM
        // call without broadcast
        inline long pack_handle(long handle, ProcessID root, bool immediate) const {
            handle &= 1023L; // To get rid of any junk after 10 bits
            if (!immediate) handle |= 1024L; // 11th bit
            return handle | ((root+1) << 11); // remainder is root+1 in broadcast, root=-1 means not broadcast
        };
        
        inline long unpack_function(long handle) const {
            return handle & 1023L;
        };
        
        inline ProcessID unpack_root(long handle) const {
            return (handle >> 11) - 1;
        };
        
        inline bool unpack_immediate(long handle) const {
            return (handle&1024) == 0;
        };
        
        
#ifdef USE_REAL_USLEEP
        inline void backoff(unsigned long& count) {
            count++;
            if (count < 3) return;
            else if (count < 1000) usleep(1000);  // 1s at 1ms intervals
            else if (count < 2000) usleep(10000); // 10s at 10ms intervals
            else if (count < 3000) usleep(100000); // 100s at 100ms intervals
            else throw "Deadlocked inside backoff"; // FOR DEBUGGING AM ...
        };
#else
        inline void backoff(unsigned long& count) {};

        // Jaguar has really bad problems if the real usleep is called
        inline void usleep(int i) {};
#endif


        Communicator(const Communicator&); // Forbidden
        Communicator& operator=(const Communicator&); // Forbidden

    public:
    
        /// Yield control to other threads
        
        /// Only really useful when oversubscribing processors. 
        inline void yield() {usleep(10);};

        /// Use given communicator and setup topology as 3D mesh (P=2^n)
        Communicator(MPI::Intracomm comm) : _comm(comm) {
            setup();
        };

        /// Use MPI_COMM_WORLD and setup topology as 3D mesh (P=2^n)
        Communicator() : _comm(MPI::COMM_WORLD) {
            setup();
        };
        
        /// Set debug flag to new value and return old value
        bool set_debug(bool value) {
            bool status = debug;
            debug = value;
            return status;
        };

        /// Return the process rank within communnicator
        inline ProcessID rank() const {
            return _rank;
        };

        /// Return the number of processes within communicator
        inline long nproc() const {
            return _nproc;
        };

        /// Return the number of processes within communicator
        inline long size() const {
            return _nproc;
        };
        
        /// Return the MPI communicator
        inline MPI::Intracomm& mpi_comm() {
            return _comm;
        };

        /// Return coords of this process in the process mesh
        inline void coords(long& px, long& py, long& pz) const {
            px=_px;
            py=_py;
            pz=_pz;
        };

        /// Return the dimensions of the process mesh
        inline void mesh(long& npx, long& npy, long& npz) const {
            npx=_npx;
            npy=_npy;
            npz=_npz;
        };

        /// Return rank of process given coords in process mesh
        inline ProcessID rank_from_coords(long px, long py, long pz) const {
            return px + _npx*(py + _npy*pz);
        };

        void print() const {
            std::cout << "process topology nproc=" << _nproc
            << " rank=" << _rank
            << " mesh=(" << _npx << "," << _npy << "," << _npz << ")"
            << " coords=(" << _px << "," << _py << "," << _pz << ")" << std::endl;
            std::cout.flush();
        };
        
        void print_handlers() const {
            _am_handle_manager.dump();
        };
        
        inline int unique_tag() {
            int result = _unique_tag_next;
            _unique_tag_next++;
            if (_unique_tag_next > MPI::TAG_UB) _unique_tag_next=2048;
            return result;
        };

        /// Same as MPI::Intracomm::Send
        inline void Send(const void* buf, int count, const MPI::Datatype& datatype,
                         ProcessID dest, int tag) const {
            if (debug) madness::print("Comm: sending",count,"bytes to",dest,"with tag",tag);
//            _comm.Send(buf,count,datatype,dest,tag);
            MPI::COMM_WORLD.Send(buf,count,datatype,dest,tag);
            if (debug) madness::print("Comm: sent");
        };

        /// Same as MPI::Intracomm::Isend
        inline MPI::Request Isend(const void* buf, int count, const MPI::Datatype& datatype,
                         ProcessID dest, int tag) const {
            if (debug) madness::print("Comm: Isending",count,"bytes to",dest,"with tag",tag);
            return _comm.Isend(buf,count,datatype,dest,tag);
        };


        /// Send array of lenbuf elements to process dest 
        template <class T>
        inline void Send(const T* buf, long lenbuf, ProcessID dest, int tag) const {
            Send((void* )buf, lenbuf*sizeof(T), MPI::BYTE, dest, tag);
        }


        /// Send datum to process dest 
        template <class T>
        inline
        typename madness::enable_if_c< !madness::is_pointer<T>::value, void>::type
        Send(const T& datum, ProcessID dest, int tag) const {
            Send((void* )&datum, sizeof(T), MPI::BYTE, dest, tag);
        }


        /// Same as MPI::Intracomm::Recv
        inline void Recv(void* buf, int count, const MPI::Datatype& datatype,
                         ProcessID source, int tag, MPI::Status& status) const {
            if (debug) madness::print("Comm: receiving",count,"bytes from",source,"with tag",tag);
//            _comm.Recv(buf,count,datatype,source,tag,status);
            MPI::COMM_WORLD.Recv(buf,count,datatype,source,tag,status);
            if (debug) madness::print("Comm: received");
        };


        /// Same as MPI::Intracomm::Recv
        inline void Recv(void* buf, int count, const MPI::Datatype& datatype,
                         ProcessID source, int tag) const {
            if (debug) madness::print("Comm: receiving",count,"bytes from",source,"with tag",tag);
//            _comm.Recv(buf,count,datatype,source,tag);
            MPI::COMM_WORLD.Recv(buf,count,datatype,source,tag);
            if (debug) madness::print("Comm: received");
        };


        /// Receive data of up to lenbuf elements from process dest with default tag=665
        template <class T>
        inline void
        Recv(T* buf, long lenbuf, ProcessID src, int tag=665) const {
            Recv(buf, lenbuf*sizeof(T), MPI::BYTE, src, tag);
        }


        /// Receive data of up to lenbuf elements from process dest with status
        template <class T>
        inline void
        Recv(T* buf, long lenbuf, ProcessID src, int tag, MPI::Status& status) const {
            Recv(buf, lenbuf*sizeof(T), MPI::BYTE, src, tag, status);
        }


        /// Receive datum from process src with default tag=666
        template <class T>
        inline
        typename madness::enable_if_c< !madness::is_pointer<T>::value, void>::type
        Recv(T& buf, ProcessID src, int tag=666) const {
            Recv(&buf, sizeof(T), MPI::BYTE, src, tag);
        }


        /// Receive datum from process src with status
        template <class T>
        inline
        typename madness::enable_if_c< !madness::is_pointer<T>::value, void>::type
        Recv(T& buf, ProcessID src, int tag, MPI::Status& status) const {
            Recv(&buf, sizeof(T), MPI::BYTE, src, tag, status);
        }


        /// Same as MPI::Intracomm::Irecv
        inline MPI::Request Irecv(void* buf, int count, const MPI::Datatype& datatype,
                                  ProcessID source, int tag) const {
            if (debug) madness::print("Comm: posting async receive",count,"elements from",source,"with tag",tag);
//            return _comm.Irecv(buf, count, datatype, source, tag);
            return MPI::COMM_WORLD.Irecv(buf, count, datatype, source, tag);
        };


        /// Async receive data of up to lenbuf elements from process dest
        template <class T>
        inline MPI::Request
        Irecv(T* buf, int count, ProcessID source, int tag) const {
            if (debug) madness::print("Comm: posting async receive",count,"bytes from",source,"with tag",tag);
//            return _comm.Irecv(buf, count*sizeof(T), MPI::BYTE, source, tag);
            return MPI::COMM_WORLD.Irecv(buf, count*sizeof(T), MPI::BYTE, source, tag);
        }


        /// Async receive datum from process dest
        template <class T>
        inline
        typename madness::enable_if_c< !madness::is_pointer<T>::value, MPI::Request>::type
        Irecv(T& buf, ProcessID source, int tag) const {
//            return _comm.Irecv(&buf, sizeof(T), MPI::BYTE, source, tag);
            return MPI::COMM_WORLD.Irecv(&buf, sizeof(T), MPI::BYTE, source, tag);
        }


	/// Same as MPI::Intracomm::Allreduce
	inline void Allreduce(void* sendbuf, void* recvbuf, int count, const MPI::Datatype& datatype,
			 const MPI::Op& op) const {
//	    _comm.Allreduce(sendbuf, recvbuf, count, datatype, op);
	    MPI::COMM_WORLD.Allreduce(sendbuf, recvbuf, count, datatype, op);
	};


	/// Same as MPI::Intracomm::Reduce
	inline void Reduce(void* sendbuf, void* recvbuf, int count, const MPI::Datatype& datatype,
			  const MPI::Op& op, int root) const {
//	    _comm.Reduce(sendbuf, recvbuf, count, datatype, op, root);
	    MPI::COMM_WORLD.Reduce(sendbuf, recvbuf, count, datatype, op, root);
	};

        /// Same as MPI::Intracomm::Bcast
        inline void Bcast(void* buffer, int count, const MPI::Datatype& datatype,
                          int root) const {
//            _comm.Bcast(buffer,count,datatype,root);
            MPI::COMM_WORLD.Bcast(buffer,count,datatype,root);
        };


        /// Broadcast an array of count elements
        template <class T>
        inline void Bcast(T* buffer, int count, int root) const {
//            _comm.Bcast(buffer,count*sizeof(T),MPI::BYTE,root);
            MPI::COMM_WORLD.Bcast(buffer,count*sizeof(T),MPI::BYTE,root);
        }


        /// Broadcast a datum
        template <class T>
        inline void Bcast(T& buffer, int root) const {
//            _comm.Bcast(&buffer, sizeof(T), MPI::BYTE,root);
            MPI::COMM_WORLD.Bcast(&buffer, sizeof(T), MPI::BYTE,root);
        }



        /// Same as MPI::Intracomm::Iprobe
        inline bool Iprobe(ProcessID source, int tag, MPI::Status& status) const {
//            return _comm.Iprobe(source, tag, status);
            return MPI::COMM_WORLD.Iprobe(source, tag, status);
        };


        /// Same as MPI::Intracomm::Iprobe
        inline bool Iprobe(ProcessID source, int tag) const {
//            return _comm.Iprobe(source,tag);
            return MPI::COMM_WORLD.Iprobe(source,tag);
        };


        /// Error abort with integer code
        void Abort(int code=1) const {
//            _comm.Abort(code);
            MPI::COMM_WORLD.Abort(code);
        };
        
        /// Group synchronization
        void Barrier() {
            if (debug) madness::print(rank(),"Comm: entering barrier");
            _comm.Barrier();
            if (debug) madness::print(rank(),"Comm: leaving barrier");
        };


        /// Global reduction of a scalar via MPI::Allreduce
        template <typename T>
        T global_sum(const T t, const MPI::Op& op = MPI::SUM) {
            T result;
            if (debug) madness::print(rank(),"Comm: global sum of scalar");
            _comm.Allreduce(&t, &result, 1, MPITypeFromType<T>(), op);
            if (debug) madness::print(rank(),"Comm: global sum done");
            return result;
        }
        
        
        /// Inplace global reduction of an array via MPI::Allreduce
        template <typename T>
        void global_sum(T* t, int count, const MPI::Op& op = MPI::SUM) {
            T* result = new T[count];
            if (debug) madness::print("Comm: global sum of vector",count);
            _comm.Allreduce(t, result, count, MPITypeFromType<T>(), op);
            for (int i=0; i<count; i++) t[i] = result[i];
            delete [] result;
            if (debug) madness::print("Comm: global sum done");
        }


        /// Register an "active message" handler
        int am_register(void (*handler)(Communicator&, ProcessID, const AMArg&)) {
            long handle = _am_handle_manager.insert(handler); 
            if (debug) madness::print("Comm:  registering handler", (void *) handler,"as",handle);
            return handle;
        };

        /// Suspend processing of active messages for critical sections
        inline void am_suspend() {
            // So that calls to suspend/resume nest correctly increment 
            // a counter rather than toggle true/false.  
            // Inside am_poll, if (_am_processing) then don't do anything.
            _am_processing++;
        };
        
        /// Resume processing of active messages 
        inline void am_resume() {_am_processing--;}; 

        /// Poll for and execute "active message" sent to this process
        inline void am_poll() {
            if (_am_processing) return;  // To short circuit recursive calls
            MPI::Status status;
            while (_am_req.Test(status)) {
                am_suspend();
                _am_nrecv++;
                // To make this routine reentrant: copy the amarg, repost the irecv
                // before calling the handler and remove am_suspend/resume.
                ProcessID src = status.Get_source();
                long handle = unpack_function(_am_arg.function); 
                ProcessID root = unpack_root(_am_arg.function);
                bool immediate = unpack_immediate(_am_arg.function);
                am_handlerT op = _am_handle_manager.get_pointer(handle);
                if (root >= 0) {
                    if (debug) madness::print(rank(),"am_poll handling broadcast of function",handle,src,root,immediate);
                    am_broadcast(handle, _am_arg, immediate, root);
                    src = root;
                }
                if (debug) madness::print(rank(),"am_poll calling function",handle);
                
                if (immediate) op(*this, src, _am_arg);
                else task_add_am(op, src, _am_arg);
                
                _am_req = Irecv(_am_arg, MPI::ANY_SOURCE, _am_tag);
                am_resume();
            }
        };


       private:
        /// Send an "active message" to process dest using a handle to specify handler
        
        /// If immediate is true, handler is invoked immediately upon receipt (usual AM behaviour).  Otherwise,
        /// it is inserted as a task in the task q of node dest.
        /// Root is the root node of a broadcast tree, or -1 if it is not a broadcast.
        inline void am_send(ProcessID dest, long handle, const AMArg& arg, bool immediate, ProcessID root) {
            if (debug) std::cout << rank() << " am_send for function " << handle << " to " << dest << std::endl;
            if (_am_send_active) am_wait(_am_send_req);
            // Since we are now using a non-blocking send we must take a copy of arg
            _am_send_arg = arg;
            _am_send_arg.function = pack_handle(handle,root,immediate);
            _am_send_req = _comm.Isend(&_am_send_arg, sizeof(_am_send_arg), MPI::BYTE, dest, _am_tag);
            _am_send_active = true;
            _am_nsent++;
        };
        
        /// Send an "active message" to all other nodes. 
        void am_broadcast(long handle, const AMArg& arg, bool immediate, ProcessID root) {
            ProcessID me=(rank()+nproc()-root)%nproc();
            ProcessID child0=(me<<1)+1+root, child1=(me<<1)+2+root;
            if (child0 >= nproc() && child0<(nproc()+root)) child0 -= nproc();
            if (child1 >= nproc() && child1<(nproc()+root)) child1 -= nproc();
            //madness::print(rank(),"amb",root,me,child0,child1,nproc());
            // Store the root of the tree into the top 20 bits of the function handle
            handle = pack_handle(handle,root,immediate);
            if (child0<nproc()) am_send(child0,handle,arg,immediate,root);
            if (child1<nproc()) am_send(child1,handle,arg,immediate,root);
        };
                    
        
        
        public:
        
        /// Send an "active message" to process dest using function pointer to specify handler
        inline void am_send(ProcessID dest, am_handlerT handler, const AMArg& arg) {
            am_send(dest, _am_handle_manager.get_handle(handler), arg, true, -1);
        };
        
 
        /// Send an "active message" to p and wait for a reply from q
        inline void am_send_recv(ProcessID p, am_handlerT handler, const AMArg& arg, 
                                 void* buf, long count, ProcessID q, int tag) {
            MPI::Request req = Irecv(buf, count, MPI::BYTE, q, tag);
            am_send(p,handler,arg);
            am_wait(req); 
        };

        
        /// Send an "active message" to all other nodes.
        void am_broadcast(am_handlerT handler, const AMArg& arg) {
            am_broadcast(_am_handle_manager.get_handle(handler),arg,true,rank());
        };
        
        /// Use this to wait for an MPI request to complete while processing AMs
        inline void am_wait(MPI::Request& req) {
            unsigned long count = 0;
            while (!req.Test()) {
                backoff(count);
                am_poll();
            }
        };

            
        /// AM-safe broadcast in SPMD mode of value from node zero to all other nodes
        
        /// This routine is also called by the Task class which must keep
        /// processing tasks as well as AM hence need wait() as an argument.
        template <typename T, typename waitT>
        void am_broadcast_value_spmd_doit(T& value, waitT wait) {
            // In binary tree with 0 as root, determine parent & children
            ProcessID me=rank(),parent=(me-1)>>1,child0=(me<<1)+1,child1=(me<<1)+2;
            if (me > 0) {
                MPI::Request req = Irecv(&value, sizeof(value), MPI::BYTE, parent, AM_BVALUE_TAG);
                wait(req);
            }
            MPI::Request req0, req1;
            if (child0<nproc()) req0 = Isend(&value,sizeof(value),MPI::BYTE, child0, AM_BVALUE_TAG);
            if (child1<nproc()) req1 = Isend(&value,sizeof(value),MPI::BYTE, child1, AM_BVALUE_TAG);
            if (child0<nproc()) wait(req0);
            if (child1<nproc()) wait(req1);
        }
        
        template <typename T>
        inline void am_broadcast_value_spmd(T& value) {
            am_broadcast_value_spmd_doit(value, bind_mem_fun(this,&Communicator::am_wait));
        }
        
        /// Computes on node 0 the global sum of #AM sent - #AM recvieved
        
        /// If do_am_send is true, active messages are sent to invoke this on
        /// other nodes.  Otherwise, we are in SPMD mode and all processes will
        /// eventually invoke this routine.
        ///
        /// This routine is also called by the Task class which must keep
        /// processing tasks as well as AM hence need wait()/fence() as an argument.
        template <typename waitT, typename fenceT>
        long am_ndiff_single_threaded(bool do_am_send, waitT wait, fenceT fence) {
            // In binary tree with 0 as root, determine parent & children
            ProcessID me=rank(),parent=(me-1)>>1,child0=(me<<1)+1,child1=(me<<1)+2;
            AMArg arg;
            long tmp0=0, tmp1=0;
            MPI::Request req0, req1;
            if (child0<nproc()) {
                req0 = Irecv(&tmp0, sizeof(tmp0), MPI::BYTE, child0, AM_RIGHT_TAG);
                if (do_am_send) am_send(child0,am_ndiff_handler,arg);
            }
            if (child1<nproc()) {
                req1 = Irecv(&tmp1, sizeof(tmp1), MPI::BYTE, child1, AM_LEFT_TAG);
                if (do_am_send) am_send(child1,am_ndiff_handler,arg);
            }
            if (child0<nproc()) wait(req0);
            if (child1<nproc()) wait(req1);
            
            fence();
            long sum = tmp0 + tmp1 + _am_nsent - _am_nrecv;
            if (me > 0) {
                MPI::Request req = Isend(&sum, sizeof(sum), MPI::BYTE, parent, AM_LEFT_TAG+(me&1)); // Note ugly assumption
                wait(req);
            }
            //madness::print("tmp",tmp0,tmp1,"sent/recv",_am_nsent,_am_nrecv,"sum",sum,"processing",_am_processing);                
            return sum; 
        }
        
        /// Computes on all nodes in SPMD mode the global sum of #AM sent - #AM recvieved
        
        /// This routine is also called by the Task class which must keep
        /// processing tasks as well as AM hence need wait() as an argument.
        template <typename waitT, typename fenceT>
        long am_ndiff_spmd(waitT wait, fenceT fence) {
            //madness::print("ENTERING ANST");
            //std::cout.flush();
            long sum = am_ndiff_single_threaded(false,wait,fence);
            //madness::print("ENTERING BCAST");
            //std::cout.flush();
            am_broadcast_value_spmd_doit(sum,wait);
            //madness::print("LEAVING BCAST");
            //std::cout.flush();
            return sum;
        }
        
        /// In SPMD mode, all nodes wait until global sum (#AM sent - #AM recveived) = 0
        void am_global_fence_spmd() {
            unsigned long count = 0;
            while (am_ndiff_spmd(bind_mem_fun(this,&Communicator::am_wait),noop)) {
                backoff(count);
                am_poll();
            }
        };    

        /// In single-threaded mode, node 0 waits until global sum (#AM sent - #AM recveived) = 0
        void am_global_fence_single_thread() {
            MADNESS_ASSERT(rank() == 0);
            unsigned long count = 0;
            while (am_ndiff_single_threaded(true, bind_mem_fun(this,&Communicator::am_wait),noop)) {
                backoff(count);
                am_poll();
            }
        };    
        

        void close() {
            _am_req.Cancel();
        };
    };
}

#endif
