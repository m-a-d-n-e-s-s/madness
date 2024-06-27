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

#ifndef MADNESS_WORLD_DEFERRED_CLEANUP_H__INCLUDED
#define MADNESS_WORLD_DEFERRED_CLEANUP_H__INCLUDED

#include <madness/world/worldmutex.h>  // For Mutex
#include <list>
#include <memory>

namespace madness {

    // Forward declarations
    template <typename, typename>
    class DeferredDeleter;
    class World;

    namespace detail {

        template <typename objT>
        inline void deferred_cleanup(World& world, const std::shared_ptr<objT>& p, bool assume_p_is_unique = false);

        /// Deferred cleanup of shared_ptr's

        /// Holds dynamically allocated pointers until it is ready for cleanup.
        /// \note This object is considered an implementation detail and should
        /// not be used directly, instead use \c DeferredDeleter with a
        /// \c std::shared_ptr.
        class DeferredCleanup {
        public:
            typedef std::shared_ptr<void> void_ptr; ///< input pointer type

        private:
            typedef std::list<void_ptr> void_ptr_list;

            RecursiveMutex mutex_;      ///< Worldwide mutex
            void_ptr_list deferred_;    ///< List of pointers to cleanup
            bool destroy_;              ///< Object destroy mode
                                        ///< true = destroy immediate
                                        ///< false = destroy deferred (default)

            // not allowed
            DeferredCleanup(const DeferredCleanup&);
            DeferredCleanup& operator=(const DeferredCleanup&);

            template <typename objT>
            friend void deferred_cleanup(World&, const std::shared_ptr<objT>&, bool);

            /// Access deferred cleanup object of world

            /// \param w The \c World object that holds the \c DeferredCleanup object
            /// \return A shared pointer to the deferred cleanup object of world \c w.
            static std::shared_ptr<DeferredCleanup> get_deferred_cleanup(const World& w);

        public:
            /// Construct a deferred deleter object.
            DeferredCleanup() : mutex_(), deferred_(), destroy_(false) { }

            /// Set the destruction mode

            /// \param mode true for immediate destruction, false for deferred
            /// destruction.
            void destroy(bool mode);

            /// Get the current destruction mode mode

            /// \return true for immediate destruction and false for deferred
            /// destruction.
            bool destroy() const;

            /// Adds \c item to cleanup list

            /// If destroy mode is true then the pointer is destroyed immediately.
            /// Otherwise it is stored until \c do_cleanup() is called.
            /// \param obj Object that is ready for destruction
            void add(const void_ptr& obj);

            /// Adds \c item to cleanup list

            /// If destroy mode is true then the pointer is destroyed immediately.
            /// Otherwise it is stored until \c do_cleanup() is called.
            /// \tparam objT The object pointer type
            /// \param obj Object that is ready for destruction
            template <typename objT>
            void add(const std::shared_ptr<objT>& obj) {
                add(std::static_pointer_cast<void>(obj));
            }



            /// Deletes/frees any pointers that are in the list
            void do_cleanup();
        }; // class DeferredCleanup


        /// Defer the cleanup of a shared pointer to the end of the next fence

        /// Call this function before destroying a shared pointer. If the shared
        /// pointer is the last reference to the object,
        /// or \p assume_p_is_unique is true,
        /// it is placed in the deferred deletion list.
        /// Otherwise, nothing is done with the pointer.
        template <typename objT>
        inline void deferred_cleanup(World& world, const std::shared_ptr<objT>& p,
                                     bool assume_p_is_unique) {
            const auto p_is_unique = p.use_count() == 1;
            if(p_is_unique || assume_p_is_unique) {
                // This is the last local pointer so we will place it in the
                // deferred deleter list for later cleanup.
                DeferredCleanup::get_deferred_cleanup(world)->add(p);
            }
        }

    }  // namespace detail
}  // namespace madness

#endif // MADNESS_WORLD_DEFERRED_CLEANUP_H__INCLUDED
