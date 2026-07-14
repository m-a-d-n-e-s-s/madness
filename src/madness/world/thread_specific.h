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

#ifndef MADNESS_WORLD_THREAD_SPECIFIC_H__INCLUDED
#define MADNESS_WORLD_THREAD_SPECIFIC_H__INCLUDED

/**
 \file thread_specific.h
 \brief Reclaimable thread-specific storage (a `thread_local` you can free).
 \ingroup threads
*/

#include <madness/world/madness_exception.h>

#include <pthread.h>

#include <atomic>
#include <cerrno>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <thread>
#include <utility>

namespace madness {
    namespace detail {

        /// Thread-specific storage of one \c Item per touching thread that,
        /// unlike a \c thread_local, can be reclaimed on demand.

        /// The items are owned by this object — a \c std::map keyed on
        /// \c std::thread::id, guarded by a mutex, with a per-thread
        /// \c pthread_key caching the item pointer so the steady-state
        /// \c local() is lock-free (one \c pthread_getspecific plus a
        /// generation check).  Because ownership lives in the object rather
        /// than in thread-static storage, destroying — or \c clear()ing — the
        /// pool frees every item at a defined point.
        ///
        /// This matters for large scratch buffers: a plain \c thread_local
        /// pins its storage for the life of the OS thread (≈ the process), so
        /// under the TBB / PaRSEC task backends — whose arenas can touch more
        /// distinct threads than \c MAD_NUM_THREADS — it leaks one full item
        /// per touching thread with no way to reclaim it.  Modeled on
        /// \c tbb::enumerable_thread_specific, trimmed to \c local()/\c clear().
        ///
        /// \c Item must be copy-constructible; the seed passed to the
        /// constructor is copied to initialize each thread's instance.
        /// \c local() is safe to call concurrently from different threads.
        /// \c clear() is \em not safe to call concurrently with \c local() (or
        /// while any reference previously returned by \c local() is still in
        /// use); call it only at a quiescent point, e.g. after a fence.
        template <typename Item>
        class thread_specific {
        public:
            explicit thread_specific(Item init = Item())
                    : init_(std::move(init)), gen_(0),
                      key_ptr_(new pthread_key_t, &pthread_key_deleter) {
                const int rc = pthread_key_create(key_ptr_.get(), &slot_deleter);
                if (rc == EAGAIN) {
                    MADNESS_EXCEPTION("thread_specific: PTHREAD_KEYS_MAX reached", rc);
                } else if (rc == ENOMEM) {
                    MADNESS_EXCEPTION("thread_specific: out of memory creating pthread key", rc);
                } else if (rc != 0) {
                    MADNESS_EXCEPTION("thread_specific: pthread_key_create failed", rc);
                }
            }

            // Non-copyable and non-movable: holds a pthread key and an atomic,
            // and is meant to live where it is constructed (e.g. a static).
            thread_specific(const thread_specific&) = delete;
            thread_specific& operator=(const thread_specific&) = delete;
            thread_specific(thread_specific&&) = delete;
            thread_specific& operator=(thread_specific&&) = delete;

            /// \return reference to the calling thread's \c Item, creating it on
            /// first use.  Lock-free once the thread has been seen at the
            /// current generation.
            Item& local() {
                Slot* s = static_cast<Slot*>(pthread_getspecific(key()));
                if (s != nullptr && s->gen == gen_.load(std::memory_order_acquire))
                    return *s->item;
                return local_slow(s);
            }

            /// Free every thread's \c Item.  A thread that later calls
            /// \c local() gets a freshly seeded one.  Not safe to call
            /// concurrently with \c local(); see the class note.
            void clear() {
                std::lock_guard<std::mutex> lock(mtx_);
                items_.clear();
                gen_.fetch_add(1, std::memory_order_acq_rel);
            }

            /// \return number of live per-thread items (diagnostic; takes the lock).
            std::size_t size() const {
                std::lock_guard<std::mutex> lock(mtx_);
                return items_.size();
            }

        private:
            /// Per-thread fast-path cache stored in the pthread key.  \c gen
            /// stamps which \c clear() generation \c item is valid for.
            struct Slot {
                std::uint64_t gen;
                Item* item;
            };

            pthread_key_t& key() { return *key_ptr_; }

            Item& local_slow(Slot* s) {
                std::lock_guard<std::mutex> lock(mtx_);
                const auto tid = std::this_thread::get_id();
                auto it = items_.find(tid);
                if (it == items_.end())
                    it = items_.emplace(tid, std::make_unique<Item>(init_)).first;
                if (s == nullptr) {
                    s = new Slot;
                    pthread_setspecific(key(), s);
                }
                s->item = it->second.get();
                s->gen = gen_.load(std::memory_order_relaxed);
                return *s->item;
            }

            static void slot_deleter(void* p) { delete static_cast<Slot*>(p); }
            static void pthread_key_deleter(const pthread_key_t* key_ptr) {
                if (key_ptr) {
                    pthread_key_delete(*key_ptr);
                    delete key_ptr;
                }
            }

            Item init_;                                     ///< seed for new items
            std::atomic<std::uint64_t> gen_;                ///< bumped by clear()
            mutable std::mutex mtx_;                        ///< guards items_
            std::map<std::thread::id, std::unique_ptr<Item>> items_;
            std::unique_ptr<pthread_key_t, decltype(&pthread_key_deleter)> key_ptr_;
        };

    }  // namespace detail
}  // namespace madness

#endif  // MADNESS_WORLD_THREAD_SPECIFIC_H__INCLUDED
