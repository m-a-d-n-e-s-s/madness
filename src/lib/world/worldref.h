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


  $Id$
*/


#ifndef MADNESS_WORLD_WORLDREF_H__INCLUDED
#define MADNESS_WORLD_WORLDREF_H__INCLUDED

/// \file worldref.h
/// \brief Implements RemoteReference which is for internal use

#include <world/world.h>        // for World
#include <world/atomicint.h>    // for AtomicInt
#include <world/sharedptr.h>    // for shared_ptr
#include <world/worldtypes.h>   // for ProcessID
#include <world/print.h>        // for print
#include <world/archive.h>      // for wrap_opaque
#include <world/worldam.h>      // for new_am_arg
#include <map>                  // for map

namespace madness {

    class World;
    class AmArg;
    class Mutex;
    template <typename> class ScopedMutex;
    template <typename T> class RemoteReference;

    template <typename T>
    std::ostream& operator<<(std::ostream& s, const RemoteReference<T>& ref);

    namespace detail {

        template<typename U>
        struct ptr_traits {
            typedef U & reference;
        };

        template<>
        struct ptr_traits<void> {
            typedef void reference;
        };

        /// A global pointer, valid anywhere in the world.

        /// Stores a globally addressable pointer. It can safely be sent to any
        /// process in the world.
        /// \tparam T The pointer type
        template <typename T>
        class WorldPointer {
        private:
            World* world_;      ///< A pointer to the local world
            ProcessID rank_;    ///< The rank of the node that the pointer belongs to
            T* pointer_;        ///< The pointer being referenced

        public:

            typedef T* pointer;
            typedef typename ptr_traits<T>::reference reference;

            /// Default constructor

            /// Creates a NULL pointer
            /// \throw nothing
            WorldPointer() :
                world_(NULL), rank_(-1), pointer_(NULL)
            { }

            /// World pointer constructor

            /// Construct a world pointer form a local pointer.
            /// \param w A reference to the local world.
            /// \param p The local pointer.
            /// \throw nothing
            WorldPointer(World& w, T* p) :
                world_(&w), rank_(w.mpi.rank()), pointer_(p)
            { }

            /// Copy constructor

            /// \param other The world pointer to be copied
            /// \throw nothing
            WorldPointer(const WorldPointer<T>& other) :
                world_(other.world_), rank_(other.rank_), pointer_(other.pointer_)
            { }

            /// Copy conversion constructor

            /// Copy and convert pointer from \c U* to \c T* type.
            /// \tparam U The pointer type of the \c other pointer
            /// \param other The world pointer to be copied
            /// \throw nothing
            /// \note \c U* must be implicitly convertible to T* type.
            template <typename U>
            WorldPointer(const WorldPointer<U>& other) :
                world_(other.world_), rank_(other.rank_), pointer_(other.pointer_)
            { }

            /// Copy assignment operator

            /// \param other The world pointer to be copied
            /// \return A reference to this object
            /// \throw nothing
            WorldPointer<T>& operator=(const WorldPointer<T>& other) {
                world_ = other.world_;
                rank_ = other.rank_;
                pointer_ = other.pointer_;
            }

            /// Copy conversion assignment operator

            /// Copy and convert pointer from \c U* to \c T* type.
            /// \tparam U The pointer type of the \c other pointer
            /// \param other The world pointer to be copied
            /// \return A reference to this object
            /// \throw nothing
            /// \note \c U* must be implicitly convertible to T* type.
            template <typename U>
            WorldPointer<T>& operator=(const WorldPointer<U>& other) {
                world_ = other.world_;
                rank_ = other.rank_;
                pointer_ = other.pointer_;
            }

            /// Check that the world pointer references a local pointer.

            /// \return \c true when the pointer points to a local address, and
            /// \c false when it points to a remote address or is NULL.
            /// \throw nothing
            bool is_local() const { return world_->mpi.rank() == rank_; }

            /// Pointer accessor

            /// Get the pointer referenced by the
            /// \return The local pointer.
            /// \throw MadnessException When the pointer references a remote
            /// address.
            pointer get() const {
                // It is not safe to access this pointer remotely.
                MADNESS_ASSERT(is_local());
                return (is_local() ? pointer_ : NULL);
            }

            /// Dereference operator

            /// Dereference the local pointer.
            /// \return A reference to the local pointer
            /// \throw MadnessException When the pointer references a remote
            /// address, or when the pointer is NULL.
            reference operator*() const {
                // It is not safe to access a NULL pointer with this operator.
                MADNESS_ASSERT(pointer_ != NULL);
                // It is not safe to access this pointer remotely.
                MADNESS_ASSERT(is_local());
                return *pointer_;
            }

            /// Pointer arrow operator

            /// Access members of the pointer
            /// \return The local pointer
            /// \throw MadnessException When the pointer references a remote
            /// address, or when the pointer is NULL.
            pointer operator->() const {
                // It is not safe to access a NULL pointer with this operator.
                MADNESS_ASSERT(pointer_ != NULL);
                // It is not safe to access this pointer remotely.
                MADNESS_ASSERT(is_local());
                return pointer_;
            }

            /// Boolean conversion operator

            /// \return \c true when the pointer is non-null and \c false otherwise
            /// \throw nothing
            operator bool () const { return pointer_; }

            /// Boolean conversion operator

            /// \return \c true when the pointer is null and \c false otherwise
            /// \throw nothing
            bool operator ! () const { return !pointer_; }

            // The comparison operators do not compare worlds because, at the
            // moment, multiple worlds is not supported. If that changes, we
            // need to add the world comparisons.

            /// Equality comparison operator

            /// \tparam U Another pointer type
            /// \return \c true when the pointers refer to the same address from
            /// the same node in the same world.
            /// \throw nothing
            template <typename U>
            bool operator==(const WorldPointer<U>& other) const {
                return (pointer_ == other.pointer_) && (rank_ == other.rank_);
                    //&& (world_   == other.world_);
            }

            /// Inequality comparison operator

            /// \tparam U Another pointer type
            /// \return \c true when the pointers refer to different addresses or
            /// different nodes or different worlds.
            /// \throw nothing
            template <typename U>
            bool operator!=(const WorldPointer<U>& other) const {
                return (pointer_ != other.pointer_) || (rank_ != other.rank_);
                    //|| (world_ != other.world_);
            }

            /// Less-than comparison operator

            /// This operator does a lexicographical comparison of world ID,
            /// rank, and pointer (in that order).
            /// \tparam U Another pointer type
            /// \return \c true when the lexicographical comparison of world ID,
            /// rank, and pointer is true, and false otherwise.
            /// \throw nothing
            template <typename U>
            bool operator<(const WorldPointer<U>& other) const {
                return (rank_ < other.rank_) || ((rank_ == other.rank_)
                    && (pointer_ < other.pointer_));
//                return (world_->id() < other.world_->id()) ||
//                        ((world_ == other.world_) && ((rank_ < other.rank_) ||
//                        ((rank_ == other.rank_) && (pointer_ < other.pointer_))));
            }

            /// World accessor

            /// \return A pointer to the world (may be NULL)
            /// \throw nothing
            World& get_world() const { return *world_; }

            /// Rank accessor

            /// \return The rank of the process that owns the pointer
            /// \throw nothing
            ProcessID owner() const { return rank_; }

            void print(const char* msg = "") const {
                madness::print("RemoteReference:",msg,": ptr rank id :", pointer_, rank_, world_->id());
            }

            /// Swap the content of \c this with \c other

            /// \tparam U The other world pointer type
            /// \throw nothing
            /// \note \c U* must be implicitly convertible to T* type.
            template <typename U>
            void swap(WorldPointer<U>& other) {
                std::swap(world_, other.world_);
                std::swap(rank_, other.rank_);
                std::swap(pointer_, other.pointer_);
            }

            /// Serialize/deserialize the world pointer.

            /// Serialize/deserialize the world pointer for remote communication
            /// or write to disk.
            /// \tparam Archive The archive object type.
            template <class Archive>
            inline void serialize(const Archive& ar) {
                ar & world_ & rank_ & archive::wrap_opaque(pointer_);
            }
        }; // class WorldPointer

        /// Swap the content of \c l with \c r

        /// \tparam U The other world pointer type
        /// \throw nothing
        /// \note \c U* must be implicitly convertible to T* type.
        template <typename T, typename U>
        void swap(WorldPointer<T>& l, WorldPointer<U>& r) {
            l.swap(r);
        }

        class RemoteCounterBase {
        private:
            madness::AtomicInt count_;          ///< reference count

            RemoteCounterBase(const RemoteCounterBase&);
            RemoteCounterBase& operator=(const RemoteCounterBase&);

            virtual void* get_ptr() const = 0;

        public:

            RemoteCounterBase() {
                count_ = 1;
            }

            virtual ~RemoteCounterBase() { }

            void* get() const { return this->get_ptr(); }
            long use_count() const { return count_; }
            void add_ref() { count_++; }
            bool relase() {
                return count_.dec_and_test();
            }
        }; // class RemoteCounterBase

        template <typename T>
        class RemoteCounterImpl : public RemoteCounterBase {
        private:
            typedef std::allocator<RemoteCounterImpl<T> > A;

            // Keep a copy of the shared pointer to make sure it stays in memory
            // while we have outstanding remote references to it.
            std::shared_ptr<T> pointer_; ///< pointer that is remotely referenced

            virtual void* get_ptr() const {
                return static_cast<void*>(pointer_.get());
            }

        public:
            explicit RemoteCounterImpl(const std::shared_ptr<T>& p) :
                RemoteCounterBase(), pointer_(p)
            { }

            virtual ~RemoteCounterImpl() { }

            void* operator new(std::size_t) {
                return A().allocate(1);
            }

            void operator delete(void * p) {
                A().deallocate(static_cast<RemoteCounterImpl<T> *>(p), 1);
            }
        }; // clast class RemoteCounterImpl

        class RemoteCounter {
        private:
            typedef RemoteCounterBase* pimplT;
            typedef std::map<void*, pimplT> pimpl_mapT;

            static Mutex global_mutex_;
            static pimpl_mapT pimpl_map_;

            /// Pointer to the shared counter implementation object
            pimplT pimpl_;

            /// Clean-up the implementation object

            /// Here we check that the pimpl has been initialized, and if so, we
            /// release the current reference. If the count drops to zero, then
            /// this is the last reference to the pimpl and it should be deleted.
            void destroy();

        public:

            RemoteCounter();

            RemoteCounter(const RemoteCounter& other);

            template <typename T>
            RemoteCounter(const std::shared_ptr<T>& p) :
                pimpl_(NULL)
            {
                if(p.get() != NULL) {
                    ScopedMutex<Mutex> buckleup(&global_mutex_);
                    pimpl_mapT::const_iterator it = pimpl_map_.find(static_cast<void*>(p.get()));

                    if(it == pimpl_map_.end()) {
                        // The pointer is not registered so we need to make a
                        // new counter.
                        pimpl_ = new RemoteCounterImpl<T>(p);
                        try {
                            pimpl_map_.insert(std::make_pair(static_cast<void* const>(p.get()), pimpl_));
                        } catch(...) {
                            delete pimpl_;
                            throw;
                        }
                    } else {
                        // The pointer is already registered, so we just need
                        // increment the counter.
                        pimpl_ = it->second;
                        pimpl_->add_ref();
                    }
                }
            }

            template <typename T>
            RemoteCounter(const WorldPointer<T>& wp) :
                pimpl_(NULL)
            {
                if(wp.is_local()) {
                    if(wp.get() != NULL) {
                        ScopedMutex<Mutex> buckleup(&global_mutex_);
                        pimpl_mapT::const_iterator it = pimpl_map_.find(static_cast<void*>(wp.get()));

                        MADNESS_ASSERT(it != pimpl_map_.end());
                        pimpl_ = it->second;
                    }
                }
            }

            ~RemoteCounter();

            RemoteCounter& operator=(const RemoteCounter& other);

            long use_count() const;
            bool unique() const;
            bool empty() const;
            void swap(RemoteCounter& other);

            template <typename Archive, typename T>
            void load(const Archive&, const detail::WorldPointer<T>& wp) {
                if(wp.is_local())
                    RemoteCounter(wp).swap(*this);
            }

            template <typename Archive, typename T>
            void store(const Archive&, const detail::WorldPointer<T>& wp) const {
                if(wp.is_local()) {
                    // This should not be serialized if there is no local counter
                    MADNESS_ASSERT(pimpl_ != NULL);
                    pimpl_->add_ref();
                }
            }
        }; // class RemoteCounter

        void swap(RemoteCounter& l, RemoteCounter& r);

    } // namespace detail


    /// Simple structure used to manage references/pointers to remote instances

    /// This class was intended only for internal use and is still rather
    /// poorly thought through, however, it seems to fill a wider need.
    ///
    /// Can be copied/sent as a simple, contiguous block of memory.
    ///
    /// !!! Unlike SharedPtr, copying or assigning does NOT cause any
    /// reference count to be incremented.  A good policy is that
    /// ownership is transferred to whomever you give a copy.
    ///
    /// !!! Every time you send a reference to another process use a
    /// new RemoteReference so that any SharedPtr reference count is
    /// correctly incremented.
    ///
    /// !!! If a shared pointer is passed into the constructor, its
    /// reference count is incremented.  It is YOUR RESPONSIBILITY to
    /// eventually invoke (locally or remotely) RemoteReference.dec()
    /// to decrement the reference count.  This is so that you can
    /// eliminate the implied communcation by invoking it locally as a
    /// side-effect of other communication.  If dec() is not called,
    /// you will have a memory leak.
    template <typename T>
    class RemoteReference {
        friend std::ostream& operator<< <T> (std::ostream&, const RemoteReference<T>&);

    private:
        detail::WorldPointer<T> wpointer_;  ///< World pointer
        detail::RemoteCounter counter_;     ///< Remote reference counter

        static void reset_handler(const AmArg& arg) {
            RemoteReference<T> r;
            arg & r;
            // Reset when the function ends.
        }

        /// Destroy the remote reference

        /// Calling this function resets the remote reference to its default,
        /// initialized state. This will prevent the remote reference from
        /// initiating any communication on destruction.
        /// \throw nothing
        void destroy() {
            detail::WorldPointer<T>().swap(wpointer_);
            detail::RemoteCounter().swap(counter_);
        }

    public:
        /// Makes a non-shared (no reference count) null pointer
        RemoteReference() :
            wpointer_(), counter_() {};

        /// Makes a shared reference and increments ptr reference count
        RemoteReference(World& w, std::shared_ptr<T>& p) :
            wpointer_(w, p.get()), counter_(p)
        { }

        template <typename U>
        RemoteReference(const detail::WorldPointer<U>& wp) :
            wpointer_(wp), counter_(wp)
        { }

        RemoteReference(const RemoteReference<T>& other) :
            wpointer_(other.wpointer_), counter_(other.counter_)
        { }

        template <typename U>
        RemoteReference(const RemoteReference<U>& other) :
            wpointer_(other.wpointer_), counter_(other.counter_)
        { }

        ~RemoteReference() {
            if(wpointer_ && (! is_local()) && (owner() != -1)) {
                get_world().am.send(owner(), RemoteReference<T>::reset_handler, new_am_arg(*this));
            }
        }

        RemoteReference<T>& operator=(const RemoteReference<T>& other) {
            RemoteReference<T>(other).swap(*this);
            return *this;
        }

        template <typename U>
        RemoteReference<T>& operator=(const RemoteReference<U>& other) {
            RemoteReference<T>(other).swap(*this);
            return *this;
        }


        /// Call this when you logically release the remote reference

        /// Can be called locally or remotely.
        ///
        /// When invoked, the internal information is destroyed by
        /// default since the reference is no longer valid.
        void reset() {
            // invalidate the reference
            RemoteReference<T>().swap(*this);
        }

        template <typename U>
        void reset(const std::shared_ptr<U>& p) {
            RemoteReference<T>(p).swap(*this);
        }

        template <typename U>
        void reset(const detail::WorldPointer<U>& wp) {
            RemoteReference<T>(wp).swap(*this);
        }


        /// Returns true if holding a non-zero pointer
        inline operator bool() const {
            return wpointer_;
        }

        /// Returns possibly remote pointer which will be 0 if not initialized
        inline T* get() const {
            MADNESS_ASSERT(is_local());
            return wpointer_.get();
        }

        void print(const char* msg = "") const {
            wpointer_.print(msg);
        }

        template <typename U>
        void swap(RemoteReference<U>& other) {
            madness::detail::swap(wpointer_, other.wpointer_);
            madness::detail::swap(counter_, other.counter_);
        }

        inline bool is_local() const { return /*wpointer_.is_local()*/ true; }

        /// Returns rank of owning process, or -1 if not initialized
        inline ProcessID owner() const {
            return wpointer_.owner();
        }

        World& get_world() const { return wpointer_.get_world(); }

        template <typename Archive>
        void load(const Archive& ar) {
            ar & wpointer_;
            counter_.load(ar, wpointer_);
        }

        template <typename Archive>
        void store(const Archive& ar) const {
            ar & wpointer_;
            counter_.store(ar, wpointer_);
        }

    }; // class RemoteReference

    template <typename T, typename U>
    void swap(RemoteReference<T>& l, RemoteReference<U>& r) {
        l.swap(r);
    }

#ifdef WORLD_INSTANTIATE_STATIC_TEMPLATES
    template <typename T>
    std::ostream& operator<<(std::ostream& s, const detail::WorldPointer<T>& p) {
        s << "<world pointer: ptr=" << static_cast<void*>(p.get()) << ", rank=" << p.owner() << ", id=" << p.get_world().id() << ">";
        return s;
    }

    template <typename T>
    std::ostream& operator<<(std::ostream& s, const RemoteReference<T>& ref) {
        s << "<remote: ptr=" << ref.wpointer_ << ">";
        return s;
    }
#endif

    namespace archive {
        template <typename, typename>
        struct ArchiveLoadImpl;

        template <typename, typename>
        struct ArchiveStoreImpl;

        template <typename Archive, typename T>
        struct ArchiveLoadImpl<Archive, RemoteReference<T> > {
            static inline void load(const Archive& ar, RemoteReference<T>& r) {
                r.load(ar);
            }
        };

        template <typename Archive, typename T>
        struct ArchiveStoreImpl<Archive, RemoteReference<T> > {
            static inline void store(const Archive& ar, const RemoteReference<T>& r) {
                r.store(ar);
            }
        };

    } // namespace archive
} // namespace madness

#endif // MADNESS_WORLD_WORLDREF_H__INCLUDED
