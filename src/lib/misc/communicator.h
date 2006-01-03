#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <iostream>

#include <madness_config.h>
#include <mad_types.h>

/// \file communicator.h
/// \brief Defines Communicator (interprocess communication + topology)

#ifdef HAVE_MPI
    // The C++ interface to MPI2 will be used
#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
#include <mpi.h>
#define PBEGIN_ ::MPI::Init
#define PEND_ ::MPI::Finalize
    
#elif defined(USE_TCGMSG)
    // The old TCGMSG library will be used
#include "sndrcv.h"
    
#else
    // Run sequentially by redefining the TCGMSG interface
    inline long NNODES_() {return 1;}
    inline long NODEID_() {return 0;}
    inline void SND_(long* type, char* buf, long* lenbuf, long* dest, long* sync){};
    inline void RCV_(long* type, char* buf, long* lenbuf, long* lenmes, long* src, 
                     long* nodefrom, long* sync) {};
    inline void BRDCST_(long* type, void* buf, long* lenbuf, long* root) {};
    inline void PBEGIN_(int argc, char* argv[]) {};
    inline void PEND_(){};
    inline void Error(const char* msg, long code) {
      std::cout<<"Fatal error " << msg << " (" << code << ")" << std::endl;
      std::exit(1);
    }
#endif
    
    
namespace madness {
    
    /// Holds info about process topology and provides message passing
    
    /// This class wraps an MPI or TCGMSG communicator and provides
    /// info about process rank, number of processes along with basic
    /// send/recv capabilities.
    /// 
    /// Also managed by this class is the logical layout of processes as a
    /// 3D mesh (currently P=2^n only).  This might be better separated
    /// from the communicator, but usually if you want one you want the
    /// other.
    /// 
    /// Expect this class to change as the overall design is completed.
    class Communicator {
    private:
        long _nproc;                ///< total no. of processes
        long _npx, _npy, _npz;      ///< #processes in each dimension
        long _mx, _my, _mz;         ///< mx=log2(npx)
        long _px, _py, _pz;         ///< coords of this process in mesh
        ProcessID _rank;            ///< rank of this process
#ifdef HAVE_MPI
        MPI::Intracomm _comm;
#endif
        
        /// Given p=2^n processes make as close as possible to a cubic grid.
        
        /// Each dimension will be a power of 2 with npx >= npy >= npz
        void setup() {
#ifdef HAVE_MPI
            _nproc = _comm.Get_size();
            _rank = _comm.Get_rank();
#else
            _nproc = NNODES_();
            _rank = NODEID_();
#endif
            _npz = 1; _mz = 0;
            while (8*_npz*_npz*_npz <= _nproc) {_npz *= 2; _mz++;}
            _npy = _npz; _my = _mz;
            while (4*_npy*_npy*_npz <= _nproc) {_npy *= 2; _my++;}
            _npx = _npy; _mx = _my;
            while (2*_npx*_npy*_npz <= _nproc) {_npx *= 2; _mx++;}
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
        
    public:
        
#ifdef HAVE_MPI
        /// Use given communicator and setup topology as 3D mesh (P=2^n)
        Communicator(const MPI::Intracomm comm) : _comm(comm) {setup();};
        
        /// Use MPI_COMM_WORLD and setup topology as 3D mesh (P=2^n)
        Communicator() {_comm=MPI::COMM_WORLD; setup();};
#else
        /// TCGMSG or sequential and setup topology as 3D mesh (P=2^n)
        Communicator() {setup();};
#endif
        /// Return the process rank (within this communicator)
        inline ProcessID rank() const {return _rank;};
        
        /// Return the total number of processes (within this communicator)
        inline long nproc() const {return _nproc;};
        
        /// Return the total number of processes (within this communicator)
        inline long size() const {return _nproc;};
        
        /// Return coords of this process in the process mesh (within this communicator)
        inline void coords(long& px, long& py, long& pz) const {
            px=_px; py=_py; pz=_pz;
        };
        
        /// Return the dimensions of the process mesh (within this communicator)
        inline void mesh(long& npx, long& npy, long& npz) const {
            npx=_npx; npy=_npy; npz=_npz;
        };
        
        /// Return rank of process given coords in process mesh (within this communicator)
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
        
        /// Send typed data of lenbuf elements to process dest
        
        /// Generic version assumes simple types that can be sent as a
        /// byte stream.  Complex data types will need to define
        /// specialization of this routine.  Currently this is a stupid
        /// version that maps to untyped data of length lenbuf*sizeof(T).
        template <class T> 
        inline void send(const T* buf, long lenbuf, ProcessID dest, long tag) const {
            send((void* )buf, lenbuf*sizeof(T), dest, tag);
        };
        
        /// Send typed datum to process dest
        
        /// Generic version assumes simple types that can be sent as a
        /// byte stream. 
        template <class T> 
        inline void send(const T& datum, ProcessID dest, long tag) const {
            send(&datum, 1, dest, tag);
        };
        
        /// Receive typed data of up to lenbuf elements from process dest
        
        /// If count is non-null it will be assigned the number of elements
        /// actually received. If from is non-null it will be assigned the
        /// process the data actually came from which may be of interest
        /// if src is a wild card (src=ANY_SOURCE).
        /// 
        /// This generic version will need to be specialized for complex types
        /// and currently maps directly to an untyped byte stream.
        template <class T> 
        inline void recv(T* buf, long lenbuf, ProcessID src, long tag, 
                         long *count=0, long *from=0) const {
            recv((void *) buf, lenbuf*sizeof(T), src, tag, count, from);
            if (count) *count /= sizeof(T);
        }
        
        /// Receive typed datum from process dest
        
        /// If from is non-null it will be assigned the process the data
        /// actually came from which may be of interest if src is a wild
        /// card (src=ANY_SOURCE).  
        template <class T> 
        inline void recv(T& buf, ProcessID src, long tag, long *from=0) const {
            recv(&buf, 1, src, tag, 0, from);
        }
        
        /// Broadcast typed array of lenbuf elements from rooot
        template <class T> 
        inline void bcast(T *buf, long lenbuf, ProcessID root) {
            bcast((void *) buf, lenbuf*sizeof(T), root);
        };
        
        /// Broadcast typed datum from rooot
        template <class T> 
        inline void bcast(T& buf, ProcessID root) {
            bcast(&buf, 1, root);
        };
        
        /// Broadcast untyped array of bytes from root
        inline void bcast(void *buf, long lenbuf, ProcessID root) {
#ifdef HAVE_MPI
            _comm.Bcast(buf, lenbuf, MPI::BYTE, root);
#else
            long tag = 99;
            BRDCST_(&tag, &buf, &lenbuf, &root);
#endif
        };
        
        /// Error abort with given integer code
        void abort(int code=1) {
#ifdef HAVE_MPI
            _comm.Abort(code);
#else
            Error("madness error termination",code);
#endif
        }; 
        
        /// Typed global sum ... NOT YET ACTUALLY IMPLEMENTED!
        template <typename T>
        inline void global_sum(T* t, long n) {};
               
    };
    
    /// Send untyped data of lenbuf bytes to process dest
    template <>
    inline void Communicator::send(const void* buf, long lenbuf, ProcessID dest, 
                                   long tag) const {
#ifdef HAVE_MPI
        _comm.Send(buf, lenbuf, MPI::BYTE, dest, tag);
#else
        long sync = 1;
        SND_(&tag, (char *) buf, &lenbuf, &dest, &sync);
#endif
    }

    /// Receive untyped data of up to lenbuf bytes from process dest
    
    /// If count is non-null it will be assigned the number of bytes
    /// actually received.  If from is non-null it will be assigned the
    /// process the data actually came from which may be of interest
    /// if src is a wild card (src=ANY_SOURCE).
    template <>
    inline void Communicator::recv(void* buf, long lenbuf, ProcessID src, long tag, 
                                   long *count, long *from) const {
#ifdef HAVE_MPI
        MPI::Status status;
        _comm.Recv(buf, lenbuf, MPI::BYTE, src, tag, status);
        if (count) *count = status.Get_count(MPI::BYTE);
        if (from) *from = status.Get_source();
#else
        long sync=1, nodefrom, lenmes;
        RCV_(&tag, (char *) buf, &lenbuf, &lenmes, &src, &nodefrom, &sync);
        if (count) *count = lenmes;
        if (from) *from = nodefrom;
#endif
    }
    
}

#endif
