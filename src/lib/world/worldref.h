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
#include <world/worldptr.h>     // for WorldPtr
#include <map>                  // for map

namespace madness {

    template <typename T> class RemoteReference;

    template <typename T>
    std::ostream& operator<<(std::ostream& s, const RemoteReference<T>& ref);

    namespace detail {

        class RemoteCounterBase {
        private:
            madness::AtomicInt count_;          ///< reference count

            // Copy not allowed
            RemoteCounterBase(const RemoteCounterBase&);
            RemoteCounterBase& operator=(const RemoteCounterBase&);

            virtual void* get_ptr() const = 0;

        public:

            RemoteCounterBase() { count_ = 1; }
            virtual ~RemoteCounterBase() { }

            void* get() const { return this->get_ptr(); }
            long use_count() const { return count_; }
            void add_ref() { count_++; }
            bool relase() { return count_.dec_and_test(); }
        }; // class RemoteCounterBase

        template <typename T>
        class RemoteCounterImpl : public RemoteCounterBase {
        private:
//            typedef std::allocator<RemoteCounterImpl<T> > A;

            // Keep a copy of the shared pointer to make sure it stays in memory
            // while we have outstanding remote references to it.
            std::shared_ptr<T> pointer_; ///< pointer that is remotely referenced

            virtual void* get_ptr() const { return static_cast<void*>(pointer_.get()); }

        public:
            explicit RemoteCounterImpl(const std::shared_ptr<T>& p) :
                RemoteCounterBase(), pointer_(p)
            { }

            virtual ~RemoteCounterImpl() { }

//            void* operator new(std::size_t) {
//                return A().allocate(1);
//            }
//
//            void operator delete(void * p) {
//                A().deallocate(static_cast<RemoteCounterImpl<T> *>(p), 1);
//            }
        }; // clast class RemoteCounterImpl

        class RemoteCounter {
        private:
            typedef RemoteCounterBase implT;
            typedef std::map<void*, WorldPtr<implT> > pimpl_mapT;

            static Mutex mutex_;
            static pimpl_mapT pimpl_map_;

            /// Pointer to the shared counter implementation object
            mutable WorldPtr<implT> pimpl_;

            /// Clean-up the implementation object

            /// Here we check that the pimpl has been initialized, and if so, we
            /// release the current reference. If the count drops to zero, then
            /// this is the last reference to the pimpl and it should be deleted.
            void destroy();

            static void unregister_ptr_(void* k);

            /// Register a local shared pointer

            /// This function will first search the local pointer register for
            /// the shared pointer \c p. If found the pimpl for that pointer
            /// will be returned. Otherwise a new pimpl will be created and
            /// returned.
            /// \tparam T The shared pointer type to register
            /// \param w The world where the pointer lives
            /// \param p The shared pointer to register
            /// \return A world pointer to the pimpl
            /// \throw std::bad_alloc If pimpl allocation fails.
            /// \throw madness::MadnessException If pointer cannot be inserted
            /// into the pointer registration map.
            template <typename T>
            static WorldPtr<implT> register_ptr_(World& w, const std::shared_ptr<T>& p) {
                WorldPtr<implT> result;
                if(p.get() != NULL) {
                    // Pointer is local and non-null
                    pimpl_mapT::const_iterator it =
                        pimpl_map_.find(static_cast<void*>(p.get()));

                    if(it == pimpl_map_.end()) {
                        // The pointer is not registered so we need to make a new one.
                        implT* pimpl = new RemoteCounterImpl<T>(p);
                        try {
                            MADNESS_ASSERT(pimpl != NULL);
                            result = WorldPtr<implT>(w, pimpl);
                            std::pair<pimpl_mapT::const_iterator, bool> insert_result
                                = pimpl_map_.insert(std::make_pair(static_cast<void*>(p.get()), result));
                            MADNESS_ASSERT(insert_result.second);
                        } catch(...) {
                            delete pimpl;
                            throw;
                        }
                    } else {
                        // The pointer is already registered, so we just need
                        // increment the counter.
                        result = it->second;
                        result->add_ref();
                    }
                }

                return result;
            }

        private:
            RemoteCounter(const WorldPtr<implT>& p);

        public:

            RemoteCounter();

            RemoteCounter(const RemoteCounter& other);

            template <typename T>
            explicit RemoteCounter(World& w, const std::shared_ptr<T>& p) :
                pimpl_(register_ptr_(w, p))
            { }

            ~RemoteCounter();

            RemoteCounter& operator=(const RemoteCounter& other);

            long use_count() const;
            bool unique() const;
            bool empty() const;
            void swap(RemoteCounter& other);

            template <typename Archive>
            void load_internal_(const Archive& ar) {
                WorldPtr<implT> p;
                ar & p;

                RemoteCounter(p).swap(*this);
            }

            template <typename Archive>
            void store_internal_(const Archive& ar) const {
                if(pimpl_.is_local())
                    pimpl_->add_ref();
                else
                    pimpl_ = WorldPtr<implT>();

                ar & pimpl_;
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
    private:
        friend std::ostream& operator<< <T> (std::ostream&, const RemoteReference<T>&);

        mutable detail::WorldPtr<T> wpointer_;  ///< World pointer
        detail::RemoteCounter counter_;     ///< Remote reference counter

        static void reset_handler(const AmArg& arg) {
            RemoteReference<T> r;
            arg & r;
            // r resets on scope exit.
        }

        /// Destroy the remote reference

        /// Calling this function resets the remote reference to its default,
        /// initialized state. This will prevent the remote reference from
        /// initiating any communication on destruction.
        /// \throw nothing
        void destroy() {
            detail::WorldPtr<T>().swap(wpointer_);
            detail::RemoteCounter().swap(counter_);
        }

    public:
        /// Makes a non-shared (no reference count) null pointer
        RemoteReference() :
            wpointer_(), counter_() {};

        /// Makes a shared reference and increments ptr reference count
        RemoteReference(World& w, std::shared_ptr<T>& p) :
            wpointer_(w, p.get()), counter_(w, p)
        { }

        RemoteReference(const RemoteReference<T>& other) :
            wpointer_(other.wpointer_), counter_(other.counter_)
        { }

        template <typename U>
        RemoteReference(const RemoteReference<U>& other) :
            wpointer_(other.wpointer_), counter_(other.counter_)
        { }

        ~RemoteReference() {
            if((owner() != -1) && (! is_local()))
                get_world().am.send(owner(), RemoteReference<T>::reset_handler, new_am_arg(*this));
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
        void reset(const detail::WorldPtr<U>& wp) {
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
            std::cout << msg << wpointer_;
        }

        template <typename U>
        void swap(RemoteReference<U>& other) {
            madness::detail::swap(wpointer_, other.wpointer_);
            madness::detail::swap(counter_, other.counter_);
        }

        inline bool is_local() const { return wpointer_.is_local(); }

        /// Returns rank of owning process, or -1 if not initialized
        inline ProcessID owner() const {
            return wpointer_.owner();
        }

        World& get_world() const { return wpointer_.get_world(); }

        template <typename Archive>
        void load_internal_(const Archive& ar) {
            ar & wpointer_ & counter_;
        }

        template <typename Archive>
        void store_internal_(const Archive& ar) const {
            ar & wpointer_ & counter_;
            if(! wpointer_.is_local())
                wpointer_ = detail::WorldPtr<T>();
        }

    }; // class RemoteReference

    template <typename T, typename U>
    void swap(RemoteReference<T>& l, RemoteReference<U>& r) {
        l.swap(r);
    }

    template <typename T>
    std::ostream& operator<<(std::ostream& s, const RemoteReference<T>& ref) {
        s << "<remote: ptr=" << ref.wpointer_ << ">";
        return s;
    }


    namespace archive {
        template <typename, typename>
        struct ArchiveLoadImpl;

        template <typename, typename>
        struct ArchiveStoreImpl;

        template <typename Archive, typename T>
        struct ArchiveLoadImpl<Archive, RemoteReference<T> > {
            static inline void load(const Archive& ar, RemoteReference<T>& r) {
                r.load_internal_(ar);
            }
        };

        template <typename Archive, typename T>
        struct ArchiveStoreImpl<Archive, RemoteReference<T> > {
            static inline void store(const Archive& ar, const RemoteReference<T>& r) {
                r.store_internal_(ar);
            }
        };

        template <typename Archive>
        struct ArchiveLoadImpl<Archive, detail::RemoteCounter > {
            static inline void load(const Archive& ar, detail::RemoteCounter& c) {
                c.load_internal_(ar);
            }
        };

        template <typename Archive>
        struct ArchiveStoreImpl<Archive, detail::RemoteCounter > {
            static inline void store(const Archive& ar, const detail::RemoteCounter& c) {
                c.store_internal_(ar);
            }
        };

    } // namespace archive
} // namespace madness

#endif // MADNESS_WORLD_WORLDREF_H__INCLUDED
