#ifndef SAV_H_
#define SAV_H_

/// \file sav.h
/// \brief Basic functionality for single assignment variables

#include <misc/shared_ptr.h>
#include <misc/communicator.h>
#include <tensor/tensor.h>

namespace madness {

    // NB ... we need to rework this to reuse the serialization
    // interface by adding an asynchronous archive.  This will require
    // adding a new method to the ArchiveLoadStoreImpl class and
    // maintaining additional state.  However, it should eliminate
    // any need to specialize SAVImpl.

    /*
     * SAV provide the following
     * - constructor
     * - virtual destructor
     * - probe()
     * - void set(const T&)
     * - T& get()
     * - forward(??)
     * - rank and tag used by SAV<T>
     * 
     * any_tag and any_source are often implemented as -1, so we cannot
     * use those as invalid values. sigh. 
     */
    
    template <class T> class SAV;
    
    /// Single assignment variable with support for MPI asynchronous messages
    template <typename T>
    class SAVImpl {
    protected:
        friend class SAV<T>;
        
        mutable bool assigned; ///< True if assigned
        bool posted; ///< True if an async receive has been posted
        int rank; ///< Rank of related MPI process, or INVALID_RANK if none
        int tag; ///< Tag of related MPI message, or INVALID_TAG if none
        mutable T t;  ///< The actual value (simple data type)
        mutable MPI::Request handle;

    public:
        /// Makes a local, unassigned variable
        SAVImpl() : assigned(false), posted(false), rank(INVALID_RANK), tag(INVALID_TAG), t(), handle() {};

        /// Makes a local, assigned variable
        SAVImpl(const T& t) : assigned(true), posted(false), rank(INVALID_RANK), tag(INVALID_TAG), t(t), handle() {};

        /// Makes a remote, unassigned variable.  If (postrecv) posts an async receive.

        /// If (postrecv) this variable will eventually be assigned
        /// by a message from process rank with the given tag.  Otherwise,
        /// this is an unassigned variable, which when locally set/assigned to
        /// will send a message of the given tag to process rank.
        SAVImpl(int rank, int tag, bool postrecv)
                : assigned(false), posted(postrecv), rank(rank), tag(tag), t(), handle() {
            if (postrecv) handle = madness::comm_default->Irecv(t, rank, tag);
        };

        /// Probe returns true if assigned and does not block
        inline bool probe() const {    // Logically const ...
            if (!assigned && posted) assigned = handle.Test();
            return assigned;
        };

        /// Returns const reference to the value. Throws exception if not yet assigned.
        const T& get() const {
                if (!probe()) throw "SAV: trying to get unassigned SAV";
                return t;
            };

        /// Sets the value. Throws exception if already assigned
        void set(const T& value) {
            if (probe()) throw "SAV: trying to set assigned SAV";
            if (posted) throw "SAV: trying to set an SAV with a posted recv";
            t = value;
            if (rank >= 0) madness::comm_default->Send(t, rank, tag);
            assigned = true;
        };
        
        virtual ~SAVImpl() {
            if (posted && !assigned) throw "SAV: destructor for unset remote SAV";
        };
    };

    /// Specialization for vectors of fixed & unknown length
    template <typename T>
    class SAVImpl< std::vector<T> > {
    protected:
        friend class SAV< std::vector<T> >;
        mutable bool assigned;
        bool posted;
        int rank;
        int tag;
        mutable std::vector<T> t;
        mutable int state;
        mutable MPI::Request handle;
        mutable std::size_t n;

    public:
        SAVImpl() : assigned(false), posted(false), rank(INVALID_RANK), tag(INVALID_TAG), t(), state(0), handle() {};

        SAVImpl(const std::vector<T>& t) :
        assigned(true), posted(false), rank(INVALID_RANK), tag(INVALID_TAG), t(t), state(0), handle() {};

        SAVImpl(int rank, int tag, bool postrecv)
                : assigned(false), posted(postrecv), rank(rank), tag(tag), t(), state(0), handle() {
            if (postrecv) handle = madness::comm_default->Irecv(n, rank, tag);
        };

        template <class argT>
        SAVImpl(int rank, int tag, bool postrecv, const argT& args)
                : assigned(false), posted(postrecv), rank(rank), tag(tag), t(args), state(1), handle() {
            if (postrecv) handle = madness::comm_default->Irecv(&t[0], n, rank, tag);
        };

        inline bool probe() const {
            if (!assigned && posted) {
                bool status = handle.Test();
                if (status) {
                    state++;
                    if (state == 1) {
                        t.resize(n);
                        handle = madness::comm_default->Irecv(&t[0], n, rank, tag);
                    } else {
                        assigned = true;
                    }
                }
            }
            return assigned;
        };

        const std::vector<T>& get() const {
                if (!probe()) throw "SAV: trying to get unassigned SAV";
                return t;
            };

        void set(const std::vector<T>& value) {
            if (probe()) throw "SAV: trying to set assigned SAV";
            t = value;
            if (rank >= 0) {
                n = t.size();
                madness::comm_default->Send(n, rank, tag);
                madness::comm_default->Send(&t[0], n, rank, tag);
            }
            assigned = true;
        };
        
        virtual ~SAVImpl() {
            if (posted && !assigned) throw "SAV: destructor for unset remote SAV";
        };
        
    };

    /// Specialization for Tensors of both fixed unknown size and dimension
    template <typename T>
    class SAVImpl< Tensor<T> > {
    protected:
        friend class SAV< Tensor<T> >;
        mutable bool assigned;
        bool posted;
        int rank;
        int tag;
        mutable Tensor<T> t;
        mutable int state;
        mutable MPI::Request handle;
        mutable long info[3+TENSOR_MAXDIM];

    public:
        SAVImpl() : assigned(false), posted(false), rank(INVALID_RANK), tag(INVALID_TAG), t(), state(0), handle() {};

        SAVImpl(const Tensor<T>& t) :
        assigned(true), posted(false), rank(INVALID_RANK), tag(INVALID_TAG), t(t), state(0), handle() {};

        SAVImpl(int rank, int tag, bool postrecv)
                : assigned(false), posted(postrecv), rank(rank), tag(tag), t(), state(0), handle() {
            if (postrecv) handle = madness::comm_default->Irecv(info, 3+TENSOR_MAXDIM, rank, tag);
        };

        template <class argT>
        SAVImpl(int rank, int tag, bool postrecv, const argT& args)
                : assigned(false), posted(postrecv), rank(rank), tag(tag), t(args), state(1), handle() {
            if (postrecv) {
                //t = Tensor<T>(args);
                handle = madness::comm_default->Irecv(t.ptr(), t.size, rank, tag);
            }
        };

        inline bool probe() const {
            if (!assigned && posted) {
                MPI::Status status;
                bool gotit = handle.Test(status);
                if (gotit) {
                    state++;
                    if (state == 1) {
                        if (t.id != info[1]) throw "SAV: receiving tensor of wrong type?";
                        std::cout<<info[2]<<std::endl;
                        t = Tensor<T>(info[2],info+3,false);
                        if (t.size != info[0]) throw "SAV: diasgree on size of tensor";
                        handle = madness::comm_default->Irecv(t.ptr(), t.size, rank, tag);
                    } else {
                        assigned = true;
                        if (status.Get_count(MPI::BYTE) == 0) t = Tensor<T>();
                    }
                }
            }
            return assigned;
        };

        const Tensor<T>& get() const {
                if (!probe()) throw "SAV: trying to get unassigned SAV";
                return t;
            };

        void set(const Tensor<T>& value) {
            if (probe()) throw "SAV: trying to set assigned SAV";
            t = value;
            if (rank >= 0) {
                if (state == 0) madness::comm_default->Send((long *) &t.size,
                            3+TENSOR_MAXDIM, rank, tag);
                if (t.iscontiguous())
                    madness::comm_default->Send(t.ptr(), t.size, rank, tag);
                else {
                    Tensor<T> c = copy(t);
                    madness::comm_default->Send(c.ptr(), t.size, rank, tag);
                }
            }
            assigned = true;
        };

        virtual ~SAVImpl() {
            if (posted && !assigned) throw "SAV: destructor for unset remote SAV";
        };
    };


    /// Public wrapper for single assignment variables

    /// In addition to easy memory management and privacy,
    /// the wrapping of the implementation in a shared pointer
    /// makes this quick to copy.  It is also safe since behind the
    /// scenes there is only one instance of the variable
    /// which ensures the correct semantics.
    template <class T>
    class SAV {
    public:
        typedef madness::SharedPtr< SAVImpl<T> > ptrT;

    protected:
        mutable ptrT p;  ///< If null, the variable is not yet initialized

        inline void init() const {
            if (!p) p = ptrT(new SAVImpl<T>());
        };

    public:
        /// Makes a local, unassigned variable

        /// Note that actual construction is deferred until the variable is used.
        /// This makes the default constructor very inexpensive, and permits
        /// a container to be made and its elements subsequently initialized. E.g.,
        /// \code
        /// SAV<double> args[5];
        /// for (int i=0; i<5; i++) args[i] = SAV<double>(whatever you want);
        /// \endcode
    SAV() : p() {};

        /// Copy constructor is shallow
        SAV(const SAV<T>& v) {
            v.init();
            p = v.p;
        };

        /// Makes a local, assigned variable
        SAV(const T& t) : p(new SAVImpl<T>(t)) {};

        /// Makes a remote, unassigned variable.  If (postrecv) posts an async receive.

        /// If (postrecv) this variable will eventually be assigned
        /// by a message from process rank with the given tag.  Otherwise,
        /// this is an unassigned variable, which when set/assigned to (locally)
        /// will send a message of the given tag to process rank.
        SAV(int rank, int tag, bool postrecv) : p(new SAVImpl<T>(rank,tag,postrecv)) {};

        /// Makes a remote, unassigned variable.  If (postrecv) posts an async receive.

        /// Differs from SAV(rank,tag,postrecv) only in that args is passed to the
        /// constructor of the value (i.e., t = T(args)).
        /// This cumbersome method may be revised once we have experience with its use.
        /// (an alternative might be to pass a functor/function_factory that
        /// will be called to make the object, but this implies an extra copy).
        /// If (postrecv) this variable will eventually be assigned
        /// by a message from process rank with the given tag.  Otherwise,
        /// this is an unassigned variable, which when set/assigned to (locally)
        /// will send a message of the given tag to process rank.
        template <class argT>
        SAV(int rank, int tag, bool postrecv, const argT& args)
                : p(new SAVImpl<T>(rank,tag,postrecv,args)) {};


        /// Assignment makes a shallow copy
        SAV<T>& operator=(const SAV<T>& v) {
            if (&v != this) {
                v.init();
                p = v.p;
            }
            return *this;
        };

        /// Sets or assigns a value to the underlying variable
        inline void set(const T& value) {
            init();
            p->set(value);
        };

        /// Returns a const reference to the underlying variable
        inline const T& get() const {
            init();
            return p->get();
        };

        inline bool probe() const {
            init();
            return p->probe();
        };
        
        inline void msginfo(ProcessID& rank, int& tag) const {
            init();
            rank = p->rank;
            tag = p->tag;
        };
        
        inline bool islocal() const {
            init();
            return (p->tag==INVALID_TAG) && (p->rank==INVALID_RANK);
        };

        virtual ~SAV() {};
    };
}
#endif /*SAV_H_*/
