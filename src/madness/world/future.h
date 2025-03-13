/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

/**
 \file future.h
 \brief Implements \c Future and related items.
 \ingroup futures
*/

#ifndef MADNESS_WORLD_FUTURE_H__INCLUDED
#define MADNESS_WORLD_FUTURE_H__INCLUDED

#include <atomic>
#include <vector>
#include <stack>
#include <new>
#include <madness/world/nodefaults.h>
#include <madness/world/dependency_interface.h>
#include <madness/world/stack.h>
#include <madness/world/worldref.h>
#include <madness/world/world.h>

/// \addtogroup futures
/// @{
namespace madness {

    //extern SharedCounter future_count; // For tracking memory leak

    // forward decl
    template <typename T> class Future;

    /// Human readable printing of a \c Future to a stream.

    /// \tparam T The type of future.
    /// \param[in,out] out The output stream.
    /// \param[in] f The future.
    /// \return The output stream.
    template <typename T>
    std::ostream& operator<<(std::ostream& out, const Future<T>& f);


    /// Implements the functionality of futures.

    /// \tparam T The type of future.
    template <typename T>
    class FutureImpl : private Spinlock {
        friend class Future<T>;
        friend std::ostream& operator<< <T>(std::ostream& out, const Future<T>& f);

    private:
        /// The static stack size used for callbacks and future assignment used
        /// for small stack size optimizations.
        static const int MAXCALLBACKS = 4;

        typedef Stack<CallbackInterface*, MAXCALLBACKS> callbackT; ///< Callback type.
        typedef Stack<std::shared_ptr<FutureImpl<T> >,MAXCALLBACKS> assignmentT; ///< Assignment type.

        /// A stack that stores callbacks that are invoked once the future has
        /// been assigned.
        volatile callbackT callbacks;

        /// A stack that stores future objects that are set to the same value
        /// as this future, once it has been set.
        volatile mutable assignmentT assignments;

        /// A flag indicating if the future has been set.
        std::atomic<bool> assigned;  // Use of atomic for necessary memory barriers/fences

        /// Reference to a remote future pimpl.
        RemoteReference< FutureImpl<T> > remote_ref;

        volatile T t; ///< The future data.

        /// AM handler for remote set operations.

        /// \todo Description needed.
        /// \param[in] arg Description needed.
        static void set_handler(const AmArg& arg) {
            RemoteReference< FutureImpl<T> > ref;
            archive::BufferInputArchive input_arch = arg & ref;
            // The remote reference holds a copy of the shared_ptr, so no need
            // to take another.
            {
                FutureImpl<T>* pimpl = ref.get();

                ScopedMutex<Spinlock> fred(pimpl);
                if(pimpl->remote_ref) {
                    // Unarchive the value to a temporary since it is going to
                    // be forwarded to another node.
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wuninitialized-const-reference"
#endif
  		    T value;
                    input_arch & value;
#ifdef __clang__
#pragma clang diagnostic pop
#endif
	

                    // Copy world and owner from remote_ref since sending remote_ref
                    // will invalidate it.
                    World& world = pimpl->remote_ref.get_world();
                    const ProcessID owner = pimpl->remote_ref.owner();
                    world.am.send(owner, FutureImpl<T>::set_handler,
                            new_am_arg(pimpl->remote_ref, value));

                    pimpl->set_assigned(value);
                } else {
                    // Unarchive the value of the future
                    input_arch & const_cast<T&>(pimpl->t);

                    pimpl->set_assigned(const_cast<const T&>(pimpl->t));
                }
            }
            ref.reset();
        }


        /// \todo Brief description needed.

        /// Invoked locally by set routine after assignment.
        /// \todo Description needed.
        /// \param[in] value Description needed.
        inline void set_assigned(const T& value) {
            // Assume that whoever is invoking this routine is holding
            // a copy of our shared pointer on its *stack* so that
            // if this future is destroyed as a result of a callback
            // the destructor of this object is not invoked until
            // we return.
            //
            // Also assume that the caller either has the lock
            // or is sure that we are single threaded.
            MADNESS_ASSERT(!assigned);
            assigned = true;

            assignmentT& as = const_cast<assignmentT&>(assignments);
            callbackT& cb = const_cast<callbackT&>(callbacks);

            while (!as.empty()) {
                MADNESS_ASSERT(as.top());
                as.top()->set(value);
                as.pop();
            }

            while (!cb.empty()) {
                MADNESS_ASSERT(cb.top());
                cb.top()->notify();
                cb.pop();
            }

            as.reset();
            cb.reset();
        }

        /// Pass by value with implied copy to manage lifetime of \c f.

        /// \todo Description needed.
        /// \param[in] f Description needed.
        inline void add_to_assignments(const std::shared_ptr< FutureImpl<T> > f) {
            // ASSUME lock is already acquired
            if (assigned) {
                f->set(const_cast<T&>(t));
            }
            else {
                assignmentT* as = const_cast<assignmentT*>(&assignments);
                as->push(f);
            }
        }


    public:

        /// Constructor that uses a local unassigned value.
        FutureImpl()
                : callbacks()
                , assignments()
                , assigned(false)
                , remote_ref()
                , t()
        {
        }


        /// Constructor that uses a wrapper for a remote future.

        /// \todo Description needed.
        /// \param[in] remote_ref Description needed.
        FutureImpl(const RemoteReference< FutureImpl<T> >& remote_ref)
                : callbacks()
                , assignments()
                , assigned(false)
                , remote_ref(remote_ref)
                , t()
        {
        }


        /// Checks if the value has been assigned.

        /// \return True if the value has been assigned; false otherwise.
        inline bool probe() const {
            return assigned;
        }


        /// Registers a function to be invoked when future is assigned.

        /// Callbacks are invoked in the order registered. If the
        /// future is already assigned, the callback is immediately
        /// invoked.
        /// \todo Description needed.
        /// \param callback Description needed.
        inline void register_callback(CallbackInterface* callback) {
            ScopedMutex<Spinlock> fred(this);
            if (assigned) callback->notify();
            else const_cast<callbackT&>(callbacks).push(callback);
        }


        /// Sets the value of the future (assignment).

        /// \todo Descriptions needed.
        /// \tparam U Description needed.
        /// \param[in] value Description needed.
        template <typename U>
        void set(U&& value) {
            ScopedMutex<Spinlock> fred(this);
            if(remote_ref) {
                // Copy world and owner from remote_ref since sending remote_ref
                // will invalidate it.
                World& world = remote_ref.get_world();
                const ProcessID owner = remote_ref.owner();
                world.am.send(owner, FutureImpl<T>::set_handler,
                        new_am_arg(remote_ref, value));
                set_assigned(std::forward<U>(value));
            } else {
                set_assigned((const_cast<T&>(t) = std::forward<U>(value)));
            }
        }


        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \param[in] input_arch Description needed.
        void set(const archive::BufferInputArchive& input_arch) {
            ScopedMutex<Spinlock> fred(this);
            MADNESS_ASSERT(! remote_ref);
            input_arch & const_cast<T&>(t);
            set_assigned(const_cast<T&>(t));
        }


        /// Gets/forces the value, waiting if necessary.

        /// \param dowork If true (default) and the future is not ready
        /// this thread will execute other tasks while waiting
        /// \warning PaRSEC backend does not support `dowork=true`
        /// from a callstack containing a running MADNESS task already,
        /// hence must use `dowork=false` when need to call from within a task
        /// \return reference to the contained value
        /// \pre `this->is_local()`
        T& get(bool dowork = true) {
            MADNESS_ASSERT(! remote_ref);  // Only for local futures
            World::await([this] () -> bool { return this->probe(); }, dowork);
            return *const_cast<T*>(&t);
        }


        /// Gets/forces the value, waiting if necessary.

        /// \param dowork If true (default) and the future is not ready
        /// this thread will execute other tasks while waiting
        /// \warning PaRSEC backend does not support `dowork=true`
        /// from a callstack containing a running MADNESS task already,
        /// hence must use `dowork=false` when need to call from within a task
        /// \return reference to the contained value
        /// \pre `this->is_local()`
        const T& get(bool dowork = true) const {
            MADNESS_ASSERT(! remote_ref);  // Only for local futures
            World::await([this] () -> bool { return this->probe(); }, dowork);
            return *const_cast<const T*>(&t);
        }

        /// \todo Brief description needed.

        /// \todo Description needed.
        /// \return Description needed.
        bool is_local() const {
            return ! remote_ref;
        }

        /// \todo Brief description needed.

        /// \todo Is this function needed?
        /// \todo Details needed.
        /// \param f Description needed.
        /// \return Description needed.
        bool replace_with(FutureImpl<T>* f) {
            MADNESS_EXCEPTION("IS THIS WORKING? maybe now we have the mutex", 0);
//            ScopedMutex<Spinlock> fred(this);
//             MADNESS_ASSERT(!world); // was return false;
//             MADNESS_ASSERT(!assigned || f->assigned);
//             if (f->world) {
//                 world = f->world;
//                 remote_ref = f->remote_ref;
//                 f->world = 0;
//             }
//             while(f->callbacks.size()) callbacks.push(f->callbacks.pop());
//             while(f->assignments.size()) assignments.push(f->assignments.pop());
            return true;
        }

        /// Destructor.

        /// \todo Perhaps a comment about its behavior.
        virtual ~FutureImpl() {
            if (const_cast<callbackT&>(callbacks).size()) {
                print("Future: uninvoked callbacks being destroyed?", assigned);
                abort();
            }
            if (const_cast<assignmentT&>(assignments).size()) {
                print("Future: uninvoked assignment being destroyed?", assigned);
                abort();
            }
        }
    }; // class FutureImpl


    /// A future is a possibly yet unevaluated value.

    /// Uses delegation to \c FutureImpl to provide desired copy/assignment
    /// semantics, as well as safe reference counting for remote futures.
    ///
    /// Since we are using futures a lot to store local values coming
    /// from containers and inside task wrappers for messages, we
    /// included in this class a value. If a future is assigned
    /// before a copy/remote-reference is taken, the shared pointer is
    /// never made. The point of this is to eliminate the two `malloc`s
    /// that must be performed for every new \c shared_ptr.
    /// \tparam T The type of future.
    /// \todo Can this detailed description be made clearer?
    template <typename T>
    class Future {

        friend std::ostream& operator<< <T>(std::ostream& out, const Future<T>& f);

    private:

        /// Pointer to the implementation object.
        std::shared_ptr< FutureImpl<T> > f;
        char buffer[sizeof(T)]; ///< Buffer to hold a single \c T object.
        T* const value; ///< Pointer to buffer when it holds a \c T object.

        /// \todo Has something to do with the "Gotchas" section in \ref futures. More detail needed.

        /// \todo Perhaps more detail here, too... At the very least, can we give it a better name?
        class dddd {};


        /// \todo Constructor for ...

        /// \todo Description needed.
        /// \param[in] blah Description needed.
        explicit Future(const dddd& blah) : f(), value(nullptr) { }

    public:
        /// \todo Brief description needed.
        typedef RemoteReference< FutureImpl<T> > remote_refT;

        /// Makes an unassigned future.
        Future() :
            f(new FutureImpl<T>()), value(nullptr)
        {
        }

        /// Makes an assigned future.

        /// \todo Description needed.
        /// \param[in] t Description needed.
        explicit Future(const T& t) :
            f(), value(new(static_cast<void*>(buffer)) T(t))
        {
        }


        /// Makes a future wrapping a remote reference.

        /// \param[in] remote_ref The remote reference.
        explicit Future(const remote_refT& remote_ref) :
                f(remote_ref.is_local() ?
                        remote_ref.get_shared() :
                        std::make_shared<FutureImpl<T> >(remote_ref)),
                //                        std::shared_ptr<FutureImpl<T> >(new FutureImpl<T>(remote_ref))),
                value(nullptr)
        {
        }


        /// Makes an assigned future from an input archive.

        /// \param[in] input_arch The input archive.
        explicit Future(const archive::BufferInputArchive& input_arch) :
            f(), value(new(static_cast<void*>(buffer)) T())
        {
            input_arch & (*value);
        }


        /// Shallow copy constructor.

        /// \param[in] other The future to copy.
        Future(const Future<T>& other) :
            f(other.f),
            value(other.value ?
                new(static_cast<void*>(buffer)) T(* other.value) :
                nullptr)
        {
            if(other.is_default_initialized())
                f.reset(new FutureImpl<T>()); // Other was default constructed so make a new f
        }

        /// Destructor.
        ~Future() {
            if(value)
                value->~T();
        }


        /// \todo Informative description needed.

        /// See "Gotchas" on \ref futures about why this exists and how to use it.
        static const Future<T> default_initializer() {
            return Future<T>(dddd());
        }

        /// Check if the future is default initialized.

        /// \return True if this future was constructed with
        ///     \c default_initializer(); false otherwise.
        bool is_default_initialized() const {
            return ! (f || value);
        }


        /// Shallow assignment operator.

        /// \param[in] other The future to copy.
        /// \return This.
        Future<T>& operator=(const Future<T>& other) {
            if(this != &other) {
                MADNESS_ASSERT(!probe());
                if(f && other.value)
                    set(other);
                else
                    f = other.f;
            }
            return *this;
        }


        /// \brief `A.set(B)`, where `A` and `B` are futures ensures `A`
        ///     has/will have the same value as `B`.

        /// An exception is thrown if `A` is already assigned since a
        /// \c Future is a single assignment variable. We don't yet
        /// track multiple assignments from unassigned futures.
        ///
        /// If `B` is already assigned, this is the same as `A.set(B.get())`,
        /// which sets `A` to the value of `B`.
        ///
        /// If `B` has not yet been assigned, the behavior is to ensure
        /// that, when `B` is assigned, both `A` and `B` will be assigned
        /// and have the same value (though they may/may not refer to
        /// the same underlying copy of the data and indeed may even
        /// be in different processes).
        /// \todo Verification needed in the param statement.
        /// \param[in] other The future `B` described above. `*this` is `A`.
        void set(const Future<T>& other) {
            MADNESS_ASSERT(f);
            if(f != other.f) {
                MADNESS_ASSERT(! f->probe());
                if (other.probe()) {
                    set(other.get());     // The easy case
                } else {
                    // Assignment is supposed to happen just once so
                    // safe to assume that this is not being messed
                    // with ... also other might invoke the assignment
                    // callback since it could have been assigned
                    // between the test above and now (and this does
                    // happen)
                    std::shared_ptr< FutureImpl<T> > ff = f; // manage lifetime of me
                    std::shared_ptr< FutureImpl<T> > of = other.f; // manage lifetime of other

                    { // BEGIN CRITICAL SECTION
                        ScopedMutex<Spinlock> fred(of.get());
                        of->add_to_assignments(ff); // Recheck of assigned is performed in here
                    } // END CRITICAL SECTION
                }
            }
        }


        /// Assigns the value.

        /// The value can only be set \em once.
        /// \param[in] value The value to be assigned.
        inline void set(const T& value) {
            MADNESS_CHECK(f);
            std::shared_ptr< FutureImpl<T> > ff = f; // manage lifetime of f
            ff->set(value);
        }

        /// Assigns the value.

        /// The value can only be set \em once.
        /// \param[in] value The value to be assigned.
        inline void set(T&& value) {
            MADNESS_CHECK(f);
            std::shared_ptr< FutureImpl<T> > ff = f; // manage lifetime of f
            ff->set(std::move(value));
        }

        /// Assigns the value.

        /// The value can only be set \em once.
        /// \todo Description needed.
        /// \param[in] input_arch Description needed.
        inline void set(const archive::BufferInputArchive& input_arch) {
            MADNESS_CHECK(f);
            std::shared_ptr< FutureImpl<T> > ff = f; // manage lifetime of f
            ff->set(input_arch);
        }


        /// Gets the value, waiting if necessary.

        /// \param dowork If true (default) and the future is not ready
        /// this thread will execute other tasks while waiting
        /// \warning PaRSEC backend does not support `dowork=true`
        /// from a callstack containing a running MADNESS task already,
        /// hence must use `dowork=false` when need to call from within a task
        /// \return reference to the contained value
        /// \pre `this->is_local()`
        inline T& get(bool dowork = true) & {
            MADNESS_CHECK(this->is_local());
            return (f ? f->get(dowork) : *value);
        }

        /// Gets the value, waiting if necessary.

        /// \param dowork If true (default) and the future is not ready
        /// this thread will execute other tasks while waiting
        /// \warning PaRSEC backend does not support `dowork=true`
        /// from a callstack containing a running MADNESS task already,
        /// hence must use `dowork=false` when need to call from within a task
        /// \return reference to the contained value
        /// \pre `this->is_local()`
        inline const T& get(bool dowork = true) const & {
            MADNESS_CHECK(this->is_local());
            return (f ? f->get(dowork) : *value);
        }

        /// Gets the value, waiting if necessary.

        /// \param dowork If true (default) and the future is not ready
        /// this thread will execute other tasks while waiting
        /// \warning PaRSEC backend does not support `dowork=true`
        /// from a callstack containing a running MADNESS task already,
        /// hence must use `dowork=false` when need to call from within a task
        /// \return the contained value
        /// \pre `this->is_local()`
        inline T get(bool dowork = true) && {
            // According to explanation here
            // https://stackoverflow.com/questions/4986673/c11-rvalues-and-move-semantics-confusion-return-statement
            // should not return && since that will still be a ref to
            // local value.  Also, does not seem as if reference
            // lifetime extension applies (see https://abseil.io/tips/107)
            MADNESS_CHECK(this->is_local());
            if (f) {
                // N.B. it is only safe to move the data out if f is the owner (i.e. it is the sole ref)
                // it's not clear how to safely detect ownership beyond the obvious cases (see https://stackoverflow.com/questions/41142315/why-is-stdshared-ptrunique-deprecated)
                // e.g. in the common case f is set by a task, hence its refcount is 2 until the task has completed *and* been destroyed
                // so will tread lightly here ... don't move unless safe
                auto &value_ref = f->get(dowork);
                // if the task setting f is gone *and* this is the only ref, move out
                // to prevent a race with another thread that will make a copy of this Future while we are moving the data
                // atomically swap f with nullptr, then check use count of f's copy
                auto fcopy = std::atomic_exchange(&f, {});  // N.B. deprecated in C++20! need to convert f to std::atomic<std::shared_ptr<FutureImpl>>
                // f is now null and no new ref to the value can be created
                if (fcopy.use_count() == 1)
                    return std::move(value_ref);
                else
                    return value_ref;
            } else
                return std::move(*value);
        }

        /// Check whether this future has been assigned.

        /// \return True if the future has been assigned; false otherwise.
        inline bool probe() const {
            return (f ? f->probe() : bool(value));
        }

        /// Same as \c get()&.

        /// \return An lvalue reference to the value.
        inline operator T&() & { // Must return & since taskfn uses it in argument holder
            return static_cast<Future&>(*this).get();
        }

        /// Same as `get() const&`.

        /// \return A const lvalue reference to the value.
        inline operator const T&() const& { // Must return & since taskfn uses it in argument holder
            return static_cast<const Future&>(*this).get();
        }

        /// An rvalue analog of \c get().

        /// \return An rvalue reference to the value. 
        inline operator T() && {
            return static_cast<Future&&>(*this).get();
        }

        /// Returns a structure used to pass references to another process.

        /// This is used for passing pointers/references to another
        /// process. To make remote references completely safe, the
        /// \c RemoteReference increments the internal reference count of
        /// the \c Future. The counter is decremented by either
        /// assigning to the remote \c Future or its destructor if it is
        /// never assigned. The remote \c Future is \em only useful for
        /// setting the future. It will \em not be notified if the value
        /// is set elsewhere.
        ///
        /// If this is already a reference to a remote future, the
        /// actual remote reference is returned; that is, \em not a
        /// a reference to the local future. Therefore, the local
        /// future will not be notified when the result is set
        /// (i.e., the communication is short circuited).
        /// \param[in,out] world The communication world.
        /// \todo Verify the return comment.
        /// \return The remote reference.
        inline remote_refT remote_ref(World& world) const {
            MADNESS_ASSERT(!probe());
            if (f->remote_ref)
                return f->remote_ref;
            else
                return RemoteReference< FutureImpl<T> >(world, f);
        }


        /// Checks the future remoteness trait

        /// \return true if the future refers to an already available or
        /// to-be-locally-produced result
        inline bool is_local() const {
            return (f && f->is_local()) || value;
        }


        /// Checks the future remoteness trait

        /// \return true if the future refers to another future on a another rank, i.e. it will be set when the
        /// referred-to future is set
        inline bool is_remote() const {
            return !is_local();
        }


        /// Registers an object to be called when future is assigned.

        /// Callbacks are invoked in the order registered. If the
        /// future is already assigned, the callback is immediately
        /// invoked.
        /// \param[in] callback The callback to be invoked.
        inline void register_callback(CallbackInterface* callback) {
            if(probe()) {
                callback->notify();
            } else {
                MADNESS_ASSERT(f);
                f->register_callback(callback);
            }
        }
    }; // class Future


    /// A future of a future is forbidden (by deleted constructor).

    /// \tparam T The type of future.
    template <typename T>
    class Future< Future<T> > {
        Future() = delete;
    };


    /// \brief Specialization of \c FutureImpl<void> for internal convenience.
    ///     This does nothing useful!
    template <>
    class FutureImpl<void> {};

    /// \brief Specialization of \c Future<void> for internal convenience.
    ///     This does nothing useful!
    template <> class
    Future<void> {
    public:
        /// \todo Brief description needed.
        typedef RemoteReference< FutureImpl<void> > remote_refT;

        /// \todo Brief description needed.
        static const Future<void> value;


        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \param[in,out] world Description needed.
        /// \return Description needed.
        static remote_refT remote_ref(World& world) {
            return remote_refT();
        }

        Future() {}


        /// \todo Brief description needed.

        /// \todo Description needed.
        /// \param[in] remote_ref Description needed.
        Future(const RemoteReference< FutureImpl<void> >& remote_ref) {}


        /// Construct from an input archive.

        /// \param[in] input_arch The input archive.
        Future(const archive::BufferInputArchive& input_arch) {
            input_arch & *this;
        }


        /// Assignment operator.

        /// \param[in] other The future to copy.
        /// \return This.
        inline Future<void>& operator=(const Future<void>& other) {
            return *this;
        }

        /// Set the future from another \c void future.

        /// In this specialization, do nothing.
        /// \param[in] f The other future.
        static void set(const Future<void>& f) { }


        /// Set the future.

        /// In this specialization, do nothing.
        static void set() { }


        /// Check if this future has been assigned.

        /// \return True (in this specialization).
        static bool probe() {
            return true;
        }

    }; // class Future<void>


    /// Specialization of \c Future for a vector of `Future`s.

    /// Enables passing a vector of futures into a task and having the
    /// dependencies correctly tracked. Does not directly support most
    /// operations that other futures do; these are the responsibility of the
    /// individual futures in the vector.
    /// \tparam T The type of future.
    template <typename T>
    class Future< std::vector< Future<T> > > : public DependencyInterface, private NO_DEFAULTS {
    private:
        /// Alias for a vector of futures.
        typedef typename std::vector< Future<T> > vectorT;

        /// The vector of futures.
        vectorT v;

    public:
        Future() : v() { }

        /// \todo Brief description needed.

        /// \todo Description needed.
        /// \param[in] v Vector of something...
        Future(const vectorT& v) : DependencyInterface(v.size()), v(v) {
            for (int i=0; i<(int)v.size(); ++i) {
                this->v[i].register_callback(this);
            }
        }

        /// \todo Brief description needed.

        /// \todo Description needed.
        /// \param[in] input_arch Description needed.
        ///
        /// \todo Not implemented. If this is deliberate, specify why and change the tag to \\attention.
        explicit Future(const archive::BufferInputArchive& input_arch) {
            input_arch & v;
        }


        /// Access the vector of futures.

        /// \return The vector of futures.
        vectorT& get() {
            return v;
        }


        /// Access the const vector of futures.

        /// \return The vector of futures.
        const vectorT& get() const {
            return v;
        }


        /// Access the vector of futures.

        /// \return The vector of futures.
        operator vectorT& () {
            return get();
        }


        /// Access the const vector of futures.

        /// \return The vector of futures.
        operator const vectorT& () const {
            return get();
        }


        /// Check if all of the futures in the vector have been assigned.

        /// \return True if all futures have been assigned; false otherwise.
        bool probe() const {
            return DependencyInterface::probe();
        }

    }; // class Future< std::vector< Future<T> > >


    /// Factory for a vectors of futures.

    /// Rationale for this function can be found in \ref futures.
    /// \tparam T The type of future in the vector.
    /// \param[in] n The size of the vector to create.
    /// \return A vector of futures, as described in \ref futures.
    template <typename T>
    std::vector< Future<T> > future_vector_factory(std::size_t n) {
        return std::vector< Future<T> >(n, Future<T>::default_initializer());
    }


    namespace archive {

        /// Serialize an assigned future.

        /// \tparam Archive Archive type.
        /// \tparam T Future type.
        template <class Archive, typename T>
        struct ArchiveStoreImpl< Archive, Future<T>, std::enable_if_t<!std::is_void_v<T>> > {

            /// Store the assigned future in an archive.

            /// \param[in,out] ar The archive.
            /// \param[in] f The future.
            template <typename U = T, typename A = Archive>
            static inline
                std::enable_if_t<is_serializable_v<A, U>,void>
                store(const Archive& ar, const Future<T>& f) {
                MAD_ARCHIVE_DEBUG(std::cout << "serializing future" << std::endl);
                MADNESS_ASSERT(f.probe());
                ar & f.get();
            }
        };


        /// Deserialize a future into an unassigned future.

        /// \tparam Archive Archive type.
        /// \tparam T Future type.
        template <class Archive, typename T>
        struct ArchiveLoadImpl< Archive, Future<T>, std::enable_if_t<!std::is_void_v<T>> > {

            /// Read into an unassigned future.

            /// \param[in,out] ar The archive.
            /// \param[out] f The future.
            template <typename U = T, typename A = Archive>
            static inline
                std::enable_if_t<is_serializable_v<A, U>,void> load(const Archive& ar, Future<T>& f) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserializing future" << std::endl);
                MADNESS_ASSERT(!f.probe());
                T value;
                ar & value;
                f.set(value);
            }
        };


        /// Serialize an assigned future (\c void specialization).

        /// \tparam Archive Archive type.
        template <class Archive>
        struct ArchiveStoreImpl< Archive, Future<void> > {

            /// Store the assigned \c void future in the archive (do nothing).

            /// \param[in,out] ar The archive.
            /// \param[in] f The \c void future.
            static inline void store(const Archive& ar, const Future<void>& f)
            { }
        };


        /// Deserialize a future into an unassigned future (\c void specialization).

        /// \tparam Archive Archive type.
        template <class Archive>
        struct ArchiveLoadImpl< Archive, Future<void> > {

            /// Read into an unassigned \c void future.

            /// \param[in,out] ar The archive.
            /// \param[out] f The \c void future.
            static inline void load(const Archive& ar, const Future<void>& f)
            { }
        };

        /// Serialize a vector of assigned futures.

        /// \tparam Archive Archive type.
        /// \tparam T Future type.
        template <class Archive, typename T>
        struct ArchiveStoreImpl< Archive, std::vector<Future<T> > > {

            /// Store the vector of assigned futures in the archive.

            /// \param[in,out] ar The archive.
            /// \param[in] v The vector of futures.
            static inline void store(const Archive& ar, const std::vector<Future<T> >& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "serializing vector of futures" << std::endl);
                ar & v.size();
                for(typename std::vector<Future<T> >::const_iterator it = v.begin(); it != v.end(); ++it) {
                    MADNESS_ASSERT(it->probe());
                    ar & it->get();
                }
            }
        };


        /// Deserialize a vector of futures into a vector of unassigned futures.

        /// \tparam Archive Archive type.
        /// \tparam T Future type.
        template <class Archive, typename T>
        struct ArchiveLoadImpl< Archive, std::vector<Future<T> > > {

            /// Read into a vector of unassigned futures.

            /// \param[in,out] ar The archive.
            /// \param[out] v The vector of futures.
            static inline void load(const Archive& ar, std::vector<Future<T> >& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserializing vector of futures" << std::endl);
                std::size_t n = 0;
                ar & n;
                if(v.size() < n)
                    v.reserve(n);
                if(v.size() > n)
                    v.resize(n);
                for(typename std::vector<Future<T> >::iterator it = v.begin(); it < v.end(); ++it, --n) {
                    MADNESS_ASSERT(! it->probe());
                    it->set(ar);
                }
                for(; n != 0; --n)
                    v.push_back(Future<T>(ar));
            }
        };
    } // namespace archive


    // Friendly I/O to streams for futures

    /// Stream output operator for a future.

    /// \tparam T The type of future.
    /// \param[in,out] out The output stream.
    /// \param[in] f The future.
    /// \return The output stream.
    template <typename T>
    std::ostream& operator<<(std::ostream& out, const Future<T>& f);


    /// Stream output operator for a \c void future.

    /// \param[in,out] out The output stream.
    /// \param[in] f The future.
    /// \return The output stream.
    template <>
    std::ostream& operator<<(std::ostream& out, const Future<void>& f);


    /// \todo Brief description needed.

    /// \todo Descriptions needed.
    /// \tparam T Description needed.
    /// \param[in,out] The output stream.
    /// \param[in] f The future to be output.
    /// \return The output stream.
    template <typename T>
    inline std::ostream& operator<<(std::ostream& out, const Future<T>& f) {
        if (f.probe()) {
          if constexpr (is_ostreammable_v<T>) {
            out << f.get();
          }
        }
        else if (f.is_remote()) out << f.f->remote_ref;
        else if (f.f) out << "<unassigned refcnt=" << f.f.use_count() << ">";
        else out << "<unassigned>";
        return out;
    }

} // namespace madness

/// @}

#endif // MADNESS_WORLD_FUTURE_H__INCLUDED
