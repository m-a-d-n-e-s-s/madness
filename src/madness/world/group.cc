/*
  This file is part of MADNESS.

  Copyright (C) 2013  Virginia Tech

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


  $Id: group.cc 3266 2013-05-17 12:53:42Z justus.c79@gmail.com $
 */

#include <iostream>
#include <madness/world/group.h>
#include <madness/world/worldgop.h>

namespace madness {

    namespace {

        /// Group registry container type
        typedef ConcurrentHashMap<DistributedID, Future<Group> >
                group_registry_container;

        /// Group registry container
        group_registry_container group_registry;

    } // namespace

    /// Register a group

    /// Register a group so that it can be used in active messages and tasks
    /// spawned on remote nodes.
    /// \param group The group to be registered
    /// \throw TiledArray::Exception When the group is empty
    /// \throw TiledArray::Exception When the group is already in the registry
    void Group::register_group() const {
        // Get/insert the group into the registry
        group_registry_container::accessor acc;
        if(group_registry.insert(acc, group_registry_container::datumT(pimpl_->id(),
                Future<Group>::default_initializer()))) {
            acc->second = Future<Group>();
        }

        // Initialize the group object that will be used by remote tasks in
        // global operations.
        acc->second.set(Group(pimpl_.get()));
    }

    /// Remove the given group from the registry

    /// Groups are removed via a lazy sync operation, which will only remove the
    /// group from the registry once unregistered has been called on all processes
    /// in the group.
    /// \param group The group to be removed from the registry
    void Group::unregister_group(const DistributedID& did) {
        group_registry_container::accessor acc;
        [[maybe_unused]] auto found = group_registry.find(acc, did);
        MADNESS_ASSERT(found);
        group_registry.erase(acc);
    }

    /// Get a registered group

    /// This function is used to acquire the group in an active message handler.
    /// \param did The id associated with the group
    /// \return A future to the group
    Future<Group> Group::get_group(const DistributedID& did) {
        group_registry_container::accessor acc;
        if(group_registry.insert(acc, group_registry_container::datumT(did,
                Future<Group>::default_initializer())))
        {
            acc->second = Future<Group>();
        }

        return acc->second;
    }

} // namespace madness
