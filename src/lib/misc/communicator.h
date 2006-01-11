#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <iostream>

#include <madness_config.h>
#include <mad_types.h>
#include <typestuff.h>

/// \file communicator.h
/// \brief Defines Communicator (interprocess communication + topology)

#define MAD_HAVE_MPI_CXX

#ifdef MAD_HAVE_MPI_CXX
// Directly use the C++ MPI interface

#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
#include <mpi.h>
    typedef MPI::Intracomm MADMPIIntracomm;
    typedef MPI::Status MADMPIStatus;
    typedef MPI::Datatype MADMPIDatatype ;
    static const MADMPIDatatype MADMPIByte = MPI::BYTE;
#define MADMPIInit(a,b) MPI::Init(a,b)
#define MADMPIFinalize() MPI::Finalize()

#else

// Use the C MPI interface thru a crude wrapper This is mostly used
// for development under mingw using MPICH.
#include <misc/madmpiwrap.h>
#define MADMPIInit(a,b) madmpi_init(&(a),&(b))
#define MADMPIFinalize() madmpi_finalize()
    typedef int MADMPIDatatype;
    static const MADMPIDatatype MADMPIByte = 5551212;
    
class MADMPIStatus { /* Must be >= size used by MPI */
        int data[64];
    public:
        inline int Get_count(MADMPIDatatype type) const {
            if (type != MADMPIByte) throw "MADMPIStatus:stupid";
            return madmpi_status_count(data);
        };
        inline int Get_source() const {return madmpi_status_source(data);};
        inline int Get_tag() const {return madmpi_status_tag(data);};
        inline int Get_err() const {return madmpi_status_err(data);};
    };
    
    /// Wraps a C MPI communicator for the machines (mingw) without MPI C++
    class MADMPIIntracomm {
    private:
        int comm;               /* The MPICH communicator */
    public:
        MADMPIIntracomm(int comm) : comm(comm) {};
        inline int Get_size() const {return madmpi_comm_size(comm);};
        inline int Get_rank() const {return madmpi_comm_rank(comm);};
        inline void Bcast(void *buf, long lenbuf, MADMPIDatatype type, 
                          ProcessID root) const {
            if (type != MADMPIByte) throw "MADMPIIntracomm: stupid";
            madmpi_comm_bcast(buf, lenbuf, root, comm);
        };
        inline void Abort(int) const {throw "MPI:abort!";};
        inline void Send(const void* buf, long lenbuf, MADMPIDatatype type, 
                         ProcessID dest, int tag) const {
            if (type != MADMPIByte) throw "MADMPIIntracomm: stupid";
            madmpi_comm_send(buf, lenbuf, dest, tag, comm);
        };
        inline void Recv(void* buf, long lenbuf, MADMPIDatatype type , 
                         ProcessID src, int tag, MADMPIStatus& status) const {
            if (type != MADMPIByte) throw "MADMPIIntracomm: stupid";
            madmpi_comm_recv(buf, lenbuf, src, tag, comm, &status);
        };
    };
#endif

namespace madness {
    /// Holds info about process topology and provides message passing
    
    /// This class wraps an MPI communicator and provides info about
    /// process rank, number of processes along with basic send/recv
    /// capabilities.
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
        
        MADMPIIntracomm _comm;
        
        /// Given p=2^n processes make as close as possible to a cubic grid.
        
        /// Each dimension will be a power of 2 with npx >= npy >= npz
        void setup() {
            _nproc = _comm.Get_size();
            _rank = _comm.Get_rank();

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
        
        /// Use given communicator and setup topology as 3D mesh (P=2^n)
        Communicator(const MADMPIIntracomm comm) : _comm(comm) {setup();};
        
#ifdef MAD_HAVE_MPI_CXX
        /// Use MPI_COMM_WORLD and setup topology as 3D mesh (P=2^n)
        Communicator() : _comm(MPI::COMM_WORLD) {setup();};
#else
        /// Use MPI_COMM_WORLD and setup topology as 3D mesh (P=2^n)
        Communicator() : _comm(madmpi_comm_world()) {setup();};
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
        /// specialization of this routine.  Currently this is a
        /// stupid version that maps to untyped data of length
        /// lenbuf*sizeof(T).
        template <class T> 
        inline 
        void 
        send(const T* buf, long lenbuf, ProcessID dest, long tag=665) const {
            send((void* )buf, lenbuf*sizeof(T), dest, tag);
        };
        
        /// Send typed datum to process dest
        template <class T> 
        inline
        typename madness::enable_if_c< !madness::is_pointer<T>::value, void>::type
        send(const T& datum, ProcessID dest, long tag=666) const {
            send(&datum, 1, dest, tag);
        };
        
        /// Receive data of up to lenbuf elements from process dest
        
        /// If count is non-null it will be assigned the number of elements
        /// actually received. If from is non-null it will be assigned the
        /// process the data actually came from which may be of interest
        /// if src is a wild card (src=ANY_SOURCE).
        /// 
        /// This generic version will need to be specialized for complex types
        /// and currently maps directly to an untyped byte stream.
        template <class T> 
        inline 
        void 
        recv(T* buf, long lenbuf, ProcessID src, long tag=665, 
             long *count=0, long *from=0) const {
            recv((void *) buf, lenbuf*sizeof(T), src, tag, count, from);
            if (count) *count /= sizeof(T);
        }
        

        /// Receive datum from process dest
        
        /// If from is non-null it will be assigned the process the data
        /// actually came from which may be of interest if src is a wild
        /// card (src=ANY_SOURCE).  
        template <class T> 
        inline
        typename madness::enable_if_c< !madness::is_pointer<T>::value, void>::type
        recv(T& buf, ProcessID src, long tag=666, long *from=0) const {
            recv(&buf, 1, src, tag, 0, from);
        }
        
        /// Broadcast array of lenbuf elements from root
        template <class T> 
        inline 
        void
        bcast(T *buf, long lenbuf, ProcessID root) {
            bcast((void *) buf, lenbuf*sizeof(T), root);
        };
        
        /// Broadcast datum from root
        template <class T> 
        inline 
        typename madness::enable_if_c< !madness::is_pointer<T>::value, void>::type
        bcast(T& buf, ProcessID root) {
            bcast(&buf, 1, root);
        };
        
        /// Error abort with given integer code
        void abort(int code=1) {
            _comm.Abort(code);
        }; 
        
        /// Typed global sum ... NOT YET ACTUALLY IMPLEMENTED!
        template <typename T>
        inline void global_sum(T* t, long n) {};
               
    };
    
    /// Broadcast untyped array of bytes from root
    template <>
    inline void Communicator::bcast(void *buf, long lenbuf, ProcessID root) {
        _comm.Bcast(buf, lenbuf, MADMPIByte, root);
    };
        

    /// Send untyped data of lenbuf bytes to process dest
    template <>
    inline void Communicator::send(const void* buf, long lenbuf, ProcessID dest, 
                                   long tag) const {
        _comm.Send(buf, lenbuf, MADMPIByte, dest, tag);
    }

    /// Receive untyped data of up to lenbuf bytes from process dest
    
    /// If count is non-null it will be assigned the number of bytes
    /// actually received.  If from is non-null it will be assigned the
    /// process the data actually came from which may be of interest
    /// if src is a wild card (src=ANY_SOURCE).
    template <>
    inline void Communicator::recv(void* buf, long lenbuf, ProcessID src, long tag, 
                                   long *count, long *from) const {
        MADMPIStatus status;
        _comm.Recv(buf, lenbuf, MADMPIByte, src, tag, status);
        if (count) *count = status.Get_count(MADMPIByte);
        if (from) *from = status.Get_source();
    }
    
}

#endif
