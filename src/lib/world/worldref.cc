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

#include <world/world.h>
#include <world/worldref.h>
#include <world/worldmutex.h>

namespace madness {
    namespace detail {
        madness::Mutex RemoteCounter::global_mutex_;
        RemoteCounter::pimpl_mapT RemoteCounter::pimpl_map_;

        /// Clean-up the implementation object

        /// Here we check that the pimpl has been initialized, and if so, we
        /// release the current reference. If the count drops to zero, then
        /// this is the last reference to the pimpl and it should be deleted.
        void RemoteCounter::destroy() {
            if(pimpl_ != NULL) {
                if(pimpl_->relase()) {
                    {
                        ScopedMutex<Mutex> buckleup(&global_mutex_);
                        std::size_t ereased = pimpl_map_.erase(pimpl_->get());
                        MADNESS_ASSERT(ereased == 1);
                    }

                    // No one else is referencing this pointer.
                    // We can safely delete it.
                    delete pimpl_;
                }
            }

            pimpl_ = NULL;
        }


        RemoteCounter::RemoteCounter() : pimpl_(NULL) { }

        RemoteCounter::RemoteCounter(const RemoteCounter& other) :
            pimpl_(other.pimpl_)
        {
            if(pimpl_ != NULL)
                pimpl_->add_ref();
        }

        RemoteCounter::~RemoteCounter() {
            destroy();
        }

        RemoteCounter& RemoteCounter::operator=(const RemoteCounter& other) {
            RemoteCounterBase* temp = other.pimpl_;

            if(pimpl_ != temp) {
                if(temp != NULL)
                    temp->add_ref();
                destroy();
                pimpl_ = temp;
            }

            return *this;
        }

        long RemoteCounter::use_count() const { return (pimpl_ != NULL ? pimpl_->use_count() : 0); }
        bool RemoteCounter::unique() const { return use_count() == 1; }
        bool RemoteCounter::empty() const { return pimpl_ == NULL; }
        void RemoteCounter::swap(RemoteCounter& other) {
            std::swap(pimpl_, other.pimpl_);
        }

        void swap(RemoteCounter& l, RemoteCounter& r) {
            l.swap(r);
        }
    } // namespace detail
} // namespace madness

