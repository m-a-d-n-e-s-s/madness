#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

/// \file communicator.h
/// \brief Defines Communicator (interprocess communication + topology)

#include <iostream>
#include <unistd.h>
#include <vector>
#include <madness_config.h>
#include <mad_types.h>
#include <typestuff.h>

#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
#include <mpi.h>

namespace madness {
    class Communicator;

    extern Communicator* comm_default;

    void xterm_debug(const Communicator& comm, const char* path, const char* display);

    /// Holds arguments sent via active messages ... deliberately small.
    class AMArg {
    private:
        friend class Communicator;
        long function;
        AMArg() : function(-1) {};
    public:
        long arg0, arg1, arg2, arg3;

        AMArg(long function)
                : function(function) {};
        AMArg(long function, long arg0)
                : function(function), arg0(arg0) {};
        AMArg(long function, long arg0, long arg1)
                : function(function), arg0(arg0), arg1(arg1) {};
        AMArg(long function, long arg0, long arg1, long arg2)
                : function(function), arg0(arg0), arg1(arg1), arg2(arg2) {};
        AMArg(long function, long arg0, long arg1, long arg2, long arg3)
                : function(function), arg0(arg0), arg1(arg1), arg2(arg2),
        arg3(arg3) {};
    };


    void am_barrier_handler(Communicator& comm, ProcessID proc, const AMArg& arg);
    long am_barrier_nchild_registered();
    void am_barrier_zero_nchild_registered();

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
    /// While in a block of code that requres AM requests continue to
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
    ///    comm.recv(reply,p);
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
    private:
        long _nproc;                ///< total no. of processes
        long _npx, _npy, _npz;      ///< #processes in each dimension
        long _mx, _my, _mz;         ///< mx=log2(npx)
        long _px, _py, _pz;         ///< coords of this process in mesh
        ProcessID _rank;            ///< rank of this process


        // On the SGI Altix the MPI copy constructor for communicators
        // does not seem to work (even using Dup) so we have to store
        // a reference.  This is not satisfactory and we'll no doubt
        // return to this topic.
        mutable MPI::Intracomm& _comm;

        AMArg _am_arg;              // Used for async recv of AM
        MPI::Request _am_req;
        static const int _am_tag = 37919;
        int _am_barrier_handle;
        std::vector<void (*)(Communicator&, ProcessID, const AMArg&)> _am_handlers;


        /// Given p=2^n processes make as close as possible to a cubic grid.

        /// Each dimension will be a power of 2 with npx >= npy >= npz
        void setup() {
            _nproc = _comm.Get_size();
            _rank = _comm.Get_rank();

            // Register am_barrier then post AM receive buffer
            _am_barrier_handle = am_register(am_barrier_handler);
            _am_req = Irecv(_am_arg, MPI::ANY_SOURCE, _am_tag);

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

    void usleep(int i) {};
        inline void backoff(unsigned long& count) {
            count++;
            if (count < 100) return;
            else if (count < 10000) usleep(1000);  // 10s at 1ms intervals
            else if (count < 20000) usleep(10000); // 100s at 10ms intervals
            else throw "Deadlocked inside backoff"; // FOR DEBUGGING AM ...
        };

        Communicator(const Communicator&); // Forbidden
        Communicator& operator=(const Communicator&); // Forbidden

    public:

        /// Use given communicator and setup topology as 3D mesh (P=2^n)
        Communicator(MPI::Intracomm comm) : _comm(comm) {
            setup();
        };

        /// Use MPI_COMM_WORLD and setup topology as 3D mesh (P=2^n)
        Communicator() : _comm(MPI::COMM_WORLD) {
            setup();
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

        /// Same as MPI::Intracomm::Send
        inline void Send(const void* buf, int count, const MPI::Datatype& datatype,
                         ProcessID dest, int tag) const {
            _comm.Send(buf,count,datatype,dest,tag);
        };


        /// Send array of lenbuf elements to process dest with default tag=665
        template <class T>
        inline void Send(const T* buf, long lenbuf, ProcessID dest, int tag=665) const {
            Send((void* )buf, lenbuf*sizeof(T), MPI::BYTE, dest, tag);
        }


        /// Send datum to process dest with default tag=666
        template <class T>
        inline
        typename madness::enable_if_c< !madness::is_pointer<T>::value, void>::type
        Send(const T& datum, ProcessID dest, int tag=666) const {
            Send((void* )&datum, sizeof(T), MPI::BYTE, dest, tag);
        }


        /// Same as MPI::Intracomm::Recv
        inline void Recv(void* buf, int count, const MPI::Datatype& datatype,
                         ProcessID source, int tag, MPI::Status& status) const {
            _comm.Recv(buf,count,datatype,source,tag,status);
        };


        /// Same as MPI::Intracomm::Recv
        inline void Recv(void* buf, int count, const MPI::Datatype& datatype,
                         ProcessID source, int tag) const {
            _comm.Recv(buf,count,datatype,source,tag);
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
            return _comm.Irecv(buf, count, datatype, source, tag);
        };


        /// Async receive data of up to lenbuf elements from process dest
        template <class T>
        inline MPI::Request
        Irecv(T* buf, int count, ProcessID source, int tag) const {
            return _comm.Irecv(buf, count*sizeof(T), MPI::BYTE, source, tag);
        }


        /// Async receive datum from process dest
        template <class T>
        inline
        typename madness::enable_if_c< !madness::is_pointer<T>::value, MPI::Request>::type
        Irecv(T& buf, ProcessID source, int tag) const {
            return _comm.Irecv(&buf, sizeof(T), MPI::BYTE, source, tag);
        }


        /// Same as MPI::Intracomm::Bcast
        inline void Bcast(void* buffer, int count, const MPI::Datatype& datatype,
                          int root) const {
            _comm.Bcast(buffer,count,datatype,root);
        };


        /// Broadcast an array of count elements
        template <class T>
        inline void Bcast(T* buffer, int count, int root) const {
            _comm.Bcast(buffer,count*sizeof(T),MPI::BYTE,root);
        }


        /// Broadcast a datum
        template <class T>
        inline void Bcast(T& buffer, int root) const {
            _comm.Bcast(&buffer, sizeof(T), MPI::BYTE,root);
        }



        /// Same as MPI::Intracomm::Iprobe
        inline bool Iprobe(ProcessID source, int tag, MPI::Status& status) const {
            return _comm.Iprobe(source, tag, status);
        };


        /// Same as MPI::Intracomm::Iprobe
        inline bool Iprobe(ProcessID source, int tag) const {
            return _comm.Iprobe(source,tag);
        };


        /// Error abort with integer code
        void Abort(int code=1) const {
            _comm.Abort(code);
        };


        /// Typed global sum ... NOT YET ACTUALLY IMPLEMENTED!
        template <typename T>
        inline void global_sum(T* t, long n) const {}


        /// Register an "active" message handler
        int am_register(void (*handler)(Communicator&, ProcessID, const AMArg&)) {
            //std::cout << rank() << " registering " << (void *) handler << std::endl;
            _am_handlers.push_back(handler);
            return _am_handlers.size()-1;
        };


        /// Poll for and execute "active" message sent to this process
        inline void am_poll() {
            MPI::Status status;
            while (_am_req.Test(status)) {
                ProcessID src = status.Get_source();
                if (_am_arg.function<0 || _am_arg.function>=(long)_am_handlers.size())
                    throw "AM_POLL: invalid function index received";
                (*_am_handlers[_am_arg.function])(*this, src, _am_arg);
                _am_req = Irecv(_am_arg, MPI::ANY_SOURCE, _am_tag);
            }
        };


        /// Send an active message to someone
        inline void am_send(const AMArg& arg, ProcessID dest) {
            //std::cout << rank() << " am_send for function " << arg.function << " to " << dest << std::endl;
            _comm.Send(&arg, sizeof(arg), MPI::BYTE, dest, _am_tag);
        };


        /// Send an active message to p and wait for a reply from q
        inline void am_send_recv(const AMArg& arg, ProcessID p,
                                 void* buf, long count, ProcessID q, int tag) {
            MPI::Request req = Irecv(buf, count, MPI::BYTE, q, tag);
            am_send(arg,p);
            am_wait(req);
        };


        /// Use this to wait for an MPI request to complete while processing AMs
        inline void am_wait(MPI::Request& req) {
            unsigned long count = 0;
            while (!req.Test()) {
                backoff(count);
                am_poll();
            }
        }


        /// Use this barrier to keep processing AM messages while blocking

        /// Currently busy waits but better might be waiting inside MPI_Probe
        void am_barrier() {
            // In binary tree, determine parent & children
            ProcessID me=rank(),parent=(me-1)>>1,child0=(me<<1)+1,child1=(me<<1)+2,nchild=0;
            if (child0 < nproc()) nchild++;
            if (child1 < nproc()) nchild++;

            // Wait for children to check in
            unsigned long count = 0;
            while (am_barrier_nchild_registered() < nchild) {
                backoff(count);
                am_poll();
            }
            am_barrier_zero_nchild_registered();

            if (me != 0) {
                // Notify parent and wait for its reply
                am_send(AMArg(_am_barrier_handle),parent);
                count = 0;
                while (am_barrier_nchild_registered() < 1) {
                    backoff(count);
                    am_poll();
                }
                am_barrier_zero_nchild_registered();
            }

            // Notify children
            if (child0<nproc()) am_send(AMArg(_am_barrier_handle),child0);
            if (child1<nproc()) am_send(AMArg(_am_barrier_handle),child1);
        };

        void close() {
            _am_req.Cancel();
        };
    };
}

#endif
