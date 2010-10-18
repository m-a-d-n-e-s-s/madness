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

#include <world/sharedptr.h>
#include <world/worldmutex.h>  // For Mutex
#include <list>

namespace madness {
    namespace detail {

        class DeferredCleanup {
        private:
            typedef std::shared_ptr<void> void_ptr;
            typedef std::list<void_ptr> void_ptr_list;

            Mutex mutex;            ///< Worldwide mutex
            void_ptr_list deferred; ///< List of pointers to cleanup


        public:
            DeferredCleanup() : mutex(), deferred() { }

            ~DeferredCleanup() { do_cleanup(); }

            /// Adds item to list of stuff to be deleted at next global_fence()

            /// The item must be derived from DeferredCleanupInterface so that the
            /// pointer type T* can be statically cast to DeferredCleanupInterface*
            void add(const void_ptr& item);

            /// Does any deferred cleanup and returns true if cleaning was necessary
            bool do_cleanup();
        };

    }  // namespace detail
}  // namespace madness

#endif // MADNESS_WORLD_DEFERRED_CLEANUP_H__INCLUDED
