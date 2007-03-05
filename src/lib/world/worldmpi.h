#ifndef WORLDMPI_H
#define WORLDMPI_H

/// \file worldmpi.h
/// \brief Implements WorldMpiInterface

/*
// If include mpi.h BEFORE stdio/iostream should not need undefs
#ifdef SEEK_CUR
#undef SEEK_CUR
#endif
#ifdef SEEK_SET
#undef SEEK_SET
#endif
#ifdef SEEK_END
#undef SEEK_END
#endif
*/

typedef int ProcessID; //< Used to clearly identify process number/rank
typedef int Tag;       //< Used to clearly identify message tag/type

namespace madness {

    /**
      
    MPI tag management:

    0-1023            ... available for general use (not used internally)
    1024-MPI::TAG_UB  ... dynamically assigned by unique_tag() 

    TAG_UB is guaranteed by the MPI standard to be at least 32767,
    which is what LAM provides.  Most other MPIs provide the full
    range of an int.

    (odd dynamically generated tags are reserved for internal use and
    so unique_tag() only returns even tags to the user).
    
    **/
    static const Tag DYNAMIC_TAG_BASE = 1024;
        
    class WorldAmInterface;
    class WorldGopInterface;

    /// This class wraps/extends the MPI interface for World
    class WorldMpiInterface {
        friend class WorldAmInterface;
        friend class WorldGopInterface;

    private:
        // These are initialized in world.cc
        static Tag dynamic_tag_reserved; //< Used to deliver "unique" reserved tags
        static Tag dynamic_tag_general;  //< Used to deliver "unique" general tags
        
        MPI::Intracomm& _comm;  //< Associated communicator
        const ProcessID _rank;  //< MPI rank of current process
        const int _nproc;       //< No. of processes in communicator
        bool debug;             //< If true, print debug information
        unsigned int myrand_next;

    public:
        WorldMpiInterface(MPI::Intracomm& comm) 
            : _comm(comm)
            , _rank(comm.Get_rank())
            , _nproc(comm.Get_size())
            , debug(false)
            , myrand_next(comm.Get_rank())
        {
            for (int i=0; i<1000; i++) rand(); // Warmup
        };
            
        /// Set debug flag to new value and return old value
        bool set_debug(bool value) {
            bool status = debug;
            debug = value;
            return status;
        };
        
        ///////////////////////////////////////////////////////////////////
        /*      Routines that replicate and extend the MPI interface     */
        ///////////////////////////////////////////////////////////////////
        
        /// Returns the associated MPI communicator
        
        /// Note that most MPI functions are replicated and extended in
        /// class World both for convenience and to intercept these calls
        /// for instrumentation, etc.. The MPI versions should usually not
        /// be called directly.
        inline MPI::Intracomm& comm() const {return _comm;};
        
        /// Returns MPI rank of this process
        inline ProcessID rank() const {return _rank;};
        
        /// Returns size of the MPI communicator
        inline ProcessID nproc() const {return _nproc;};
        
        /// Returns size of the MPI communicator
        inline ProcessID size() const {return _nproc;};
        
        /// Returns MPI rank of this process
        inline ProcessID Get_rank() const {return _rank;};
        
        /// Returns size of the MPI communicator
        inline ProcessID Get_size() const {return _nproc;};
        

        /// Same as MPI::Intracomm::Sendrecv
        inline void Sendrecv(const void* sendbuf, int sendcount,
                             const MPI::Datatype& sendtype, ProcessID dest,
                             Tag sendtag, void* recvbuf, int recvcount,
                             const MPI::Datatype& recvtype, ProcessID source,
                             Tag recvtag) {
            _comm.Sendrecv(sendbuf, sendcount,
                           sendtype, dest,
                           sendtag, recvbuf, recvcount,
                           recvtype, source,
                           recvtag);
        };
        

        /// Same as MPI::Intracomm::Send
        inline void Send(const void* buf, int count, const MPI::Datatype& datatype,
                         ProcessID dest, Tag tag) const {
            if (debug) madness::print("World:",rank(),"sending",count,"bytes to",dest,"with tag",tag);
            _comm.Send(buf,count,datatype,dest,tag);
            if (debug) madness::print("World: sent");
        };
        
        /// Same as MPI::Intracomm::Isend
        inline MPI::Request Isend(const void* buf, int count, const MPI::Datatype& datatype,
                                  ProcessID dest, Tag tag) const {
            if (debug) madness::print("World:",rank(),"Isending",count,"bytes to",dest,"with tag",tag);
            return _comm.Isend(buf,count,datatype,dest,tag);
        };
        
        /// Send array of lenbuf elements to process dest 
        template <class T>
            inline void Send(const T* buf, long lenbuf, ProcessID dest, Tag tag) const {
            Send((void* )buf, lenbuf*sizeof(T), MPI::BYTE, dest, tag);
        }
        
        /// Send datum to process dest with default tag=1
        
        /// Disabled for pointers to reduce accidental misuse.
        template <class T>
            inline
            typename enable_if_c< !is_pointer<T>::value, void>::type
            Send(const T& datum, ProcessID dest, Tag tag=1) const {
            Send((void* )&datum, sizeof(T), MPI::BYTE, dest, tag);
        }
        
        /// Same as MPI::Intracomm::Recv with status
        inline void Recv(void* buf, int count, const MPI::Datatype& datatype,
                         ProcessID source, Tag tag, MPI::Status& status) const {
            if (debug) madness::print("World:",rank(),"receiving",count,"bytes from",source,"with tag",tag);
            _comm.Recv(buf,count,datatype,source,tag,status);
            if (debug) madness::print("World:",rank(),"received");
        };
        
        
        /// Same as MPI::Intracomm::Recv
        inline void Recv(void* buf, int count, const MPI::Datatype& datatype,
                         ProcessID source, Tag tag) const {
            if (debug) madness::print("World:",rank(),"receiving",count,"bytes from",source,"with tag",tag);
            _comm.Recv(buf,count,datatype,source,tag);
            if (debug) madness::print("World:",rank(),"received");
        };
        
        
        /// Receive data of up to lenbuf elements from process dest
        template <class T>
            inline void
            Recv(T* buf, long lenbuf, ProcessID src, Tag tag) const {
            Recv(buf, lenbuf*sizeof(T), MPI::BYTE, src, tag);
        }
        
        /// Receive data of up to lenbuf elements from process dest with status
        template <class T>
            inline void
            Recv(T* buf, long lenbuf, ProcessID src, Tag tag, MPI::Status& status) const {
            Recv(buf, lenbuf*sizeof(T), MPI::BYTE, src, tag, status);
        }
        
        
        /// Receive datum from process src with default tag=1
        template <class T>
            inline
            typename enable_if_c< !is_pointer<T>::value, void>::type
            Recv(T& buf, ProcessID src, Tag tag=1) const {
            Recv(&buf, sizeof(T), MPI::BYTE, src, tag);
        }
        
        
        /// Same as MPI::Intracomm::Irecv
        inline MPI::Request Irecv(void* buf, int count, const MPI::Datatype& datatype,
                                  ProcessID source, Tag tag) const {
            if (debug) madness::print("World:",rank(),"posting async receive",count,"elements from",source,"with tag",tag);
            return _comm.Irecv(buf, count, datatype, source, tag);
        };
        
        
        /// Async receive data of up to lenbuf elements from process dest
        template <class T>
            inline MPI::Request
            Irecv(T* buf, int count, ProcessID source, Tag tag) const {
            if (debug) madness::print("World:",rank(),"posting async receive",count,"bytes from",source,"with tag",tag);
            return _comm.Irecv(buf, count*sizeof(T), MPI::BYTE, source, tag);
        }
        
        
        /// Async receive datum from process dest with default tag=1
        template <class T>
            inline
            typename enable_if_c< !is_pointer<T>::value, MPI::Request>::type
            Irecv(T& buf, ProcessID source, Tag tag=1) const {
            return _comm.Irecv(&buf, sizeof(T), MPI::BYTE, source, tag);
        }
        
        
        /// Same as MPI::Intracomm::Allreduce
        
        /// NB.  Read documentation about interaction of MPI collectives and AM/task handling.
        inline void Allreduce(void* sendbuf, void* recvbuf, int count, const MPI::Datatype& datatype,
                              const MPI::Op& op) const {
            _comm.Allreduce(sendbuf, recvbuf, count, datatype, op);
        };
        
        
        /// Same as MPI::Intracomm::Reduce
        
        /// NB.  Read documentation about interaction of MPI collectives and AM/task handling.
        inline void Reduce(void* sendbuf, void* recvbuf, int count, const MPI::Datatype& datatype,
                           const MPI::Op& op, int root) const {
            _comm.Reduce(sendbuf, recvbuf, count, datatype, op, root);
        };
        
        /// Same as MPI::Intracomm::Bcast
        
        /// NB.  Read documentation about interaction of MPI collectives and AM/task handling.
        inline void Bcast(void* buffer, int count, const MPI::Datatype& datatype,
                          int root) const {
            _comm.Bcast(buffer,count,datatype,root);
        };
        
        
        /// MPI broadcast an array of count elements
        
        /// NB.  Read documentation about interaction of MPI collectives and AM/task handling.
        template <class T>
            inline void Bcast(T* buffer, int count, int root) const {
            _comm.Bcast(buffer,count*sizeof(T),MPI::BYTE,root);
        }
        
        
        /// MPI broadcast a datum
        
        /// NB.  Read documentation about interaction of MPI collectives and AM/task handling.
        template <class T>
            inline
            typename enable_if_c< !is_pointer<T>::value, void>::type
            Bcast(T& buffer, int root) const {
            _comm.Bcast(&buffer, sizeof(T), MPI::BYTE,root);
        }
        
        /// Same as MPI::Intracomm::Iprobe
        inline bool Iprobe(ProcessID source, Tag tag, MPI::Status& status) const {
            return _comm.Iprobe(source, tag, status);
        };
        
        
        /// Same as MPI::Intracomm::Iprobe
        inline bool Iprobe(ProcessID source, Tag tag) const {
            return _comm.Iprobe(source,tag);
        };
        
        /// Same as MPI::Intracomm::Abort
        void Abort(int code=1) const {
            _comm.Abort(code);
        };


        /// Returns a unique tag for general use (tag is even and > 1023)

        /// Unique is slightly optimistic.  The method simply
        /// increments/wraps a counter and returns the next legal
        /// value. For most MPIs with maximum tag values equal to
        /// 2^31-1 this is adequate.  For MPIs such as LAM which only
        /// provide the bare minimum of 32768 tags (of which half are
        /// reserved for internal use) you have a greater chance of
        /// collision.
        ///
        /// Ideally, we would have a separate count per World / MPI
        /// communicator but currently we have one count shared by all
        /// communicators.
        static Tag unique_tag() {
            Tag result = dynamic_tag_general;
            dynamic_tag_general += 2;
            if (dynamic_tag_general > MPI::TAG_UB) dynamic_tag_general = DYNAMIC_TAG_BASE;
            return result;
        };

    private:
        /// Private: Returns a unique tag for internal use (tag is odd and > 1023)
        static Tag unique_reserved_tag() {
            Tag result = dynamic_tag_reserved;
            dynamic_tag_reserved += 2;
            if (dynamic_tag_reserved > MPI::TAG_UB) dynamic_tag_reserved = DYNAMIC_TAG_BASE+1;
            return result;
        };
    public:

        
        /// Returns a CRUDE, LOW-QUALITY, random number uniformly distributed in [0,2**24).

        /// Each process has a distinct seed for the generator.
        int rand() {
            myrand_next = myrand_next * 1103515245UL + 12345UL;
            return int((myrand_next>>8) & 0xfffffful);
        };


        /// Returns a CRUDE, LOW-QUALITY, random number uniformly distributed in [0,1).
        double drand() {
            return rand()/16777216.0;
        };


        /// Returns a random process number [0,world.size())
        ProcessID random_proc() {
            return rand()%size();
        };


        /// Returns a random process number [0,world.size()) != current process

        /// Makes no sense to call this with just one process, but just in case you 
        /// do it returns -1 in the hope that you won't actually use the result.
        ProcessID random_proc_not_me() {
            if (size() == 1) return -1;
            ProcessID p;
            do {p = rand()%size();} while (p == rank());
            return p;
        };


        /// Construct info about a binary tree with given root
        
        /// Constructs a binary tree spanning the communicator with
        /// process root as the root of the tree.  Returns the logical
        /// parent and children in the tree of the calling process.  If
        /// there is no parent/child the value -1 will be set.
        void binary_tree_info(ProcessID root, ProcessID& parent, ProcessID& child0, ProcessID& child1) {
            int np = nproc();
            ProcessID me = (rank()+np-root)%np;   // Renumber processes so root has me=0
            parent = (((me-1)>>1)+root)%np;        // Parent in binary tree
            child0 = (me<<1)+1+root;        // Left child
            child1 = (me<<1)+2+root;        // Right child
            if (child0 >= np && child0<(np+root)) child0 -= np;
            if (child1 >= np && child1<(np+root)) child1 -= np;
            
            if (me == 0) parent = -1;
            if (child0 >= np) child0 = -1;
            if (child1 >= np) child1 = -1;
        };
    };
    
}

#endif
