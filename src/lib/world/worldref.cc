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
#include <world/worldmutex.h>

namespace madness {
    namespace detail {
        madness::Mutex RemoteCounter::mutex_;
        RemoteCounter::pimpl_mapT RemoteCounter::pimpl_map_;

        void RemoteCounter::unregister_ptr_(void* k) {
            ScopedMutex<Mutex> buckleup(&mutex_);
            std::size_t ereased = pimpl_map_.erase(k);
            MADNESS_ASSERT(ereased > 0);
        }

        /// Clean-up the implementation object

        /// Here we check that the pimpl has been initialized, and if so, we
        /// release the current reference. If the count drops to zero, then
        /// this is the last reference to the pimpl and it should be deleted.
        void RemoteCounter::destroy() {
            if(pimpl_ && pimpl_.is_local()) {
                if(pimpl_->release()) {
                    // No one else is referencing this pointer.
                    // We can safely dispose of it.
                    std::cout << ">>> RemoteCounter::unregister_ptr_: key= " << pimpl_->key() << ", value= " << pimpl_ << std::endl;
                    unregister_ptr_(pimpl_->key());
                    delete pimpl_.get();
                }
            }

            pimpl_ = WorldPtr<implT>();
        }

        RemoteCounter::RemoteCounter(const WorldPtr<implT>& p) :
            pimpl_(p)
        {
#ifndef NDEBUG
            // Check to make sure the pimpl still exists.
            if(p.is_local()) {
                ScopedMutex<Mutex> buckleup(&mutex_);
                pimpl_mapT::const_iterator it;
                for(it = pimpl_map_.begin(); it != pimpl_map_.end(); ++it)
                    if(it->second == p)
                        break;

                MADNESS_ASSERT(it != pimpl_map_.end());
            }
#endif
        }

        RemoteCounter::RemoteCounter() : pimpl_() { }

        RemoteCounter::RemoteCounter(const RemoteCounter& other) :
            pimpl_(other.pimpl_)
        {
            if(pimpl_ && pimpl_.is_local())
                pimpl_->add_ref();
        }

        RemoteCounter::~RemoteCounter() {
            destroy();
        }

        RemoteCounter& RemoteCounter::operator=(const RemoteCounter& other) {
            WorldPtr<implT> temp = other.pimpl_;

            if(pimpl_ != temp) {
                if(temp)
                    temp->add_ref();
                destroy();
                pimpl_ = temp;
            }

            return *this;
        }

        long RemoteCounter::use_count() const { return (pimpl_ && pimpl_.is_local() ? pimpl_->use_count() : 0); }
        bool RemoteCounter::unique() const { return use_count() == 1; }
        bool RemoteCounter::empty() const { return pimpl_; }
        bool RemoteCounter::is_local() const { return pimpl_.is_local(); }
        bool RemoteCounter::has_owner() const { return pimpl_.has_owner(); }
        ProcessID RemoteCounter::owner() const { return pimpl_.owner(); }
        WorldPtr<RemoteCounter::implT>::worldidT
        RemoteCounter::get_worldid() const { return pimpl_.get_worldid(); }
        World& RemoteCounter::get_world() const { return pimpl_.get_world(); }
        void RemoteCounter::swap(RemoteCounter& other) {
            madness::detail::swap(pimpl_, other.pimpl_);
        }

        void swap(RemoteCounter& l, RemoteCounter& r) {
            l.swap(r);
        }

        std::ostream& operator<<(std::ostream& out, const RemoteCounter& counter) {
            out << "RemoteCounter( owner=" << counter.owner() << " worldid=" <<
                    counter.get_worldid() << " use_count=" << counter.use_count() << ")";
            return out;
        }
    } // namespace detail
} // namespace madness

