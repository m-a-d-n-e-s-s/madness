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


  $Id$
 */

#ifndef MADNESS_WORLD_GROUP_H__INCLUDED
#define MADNESS_WORLD_GROUP_H__INCLUDED

#include <world/dist_keys.h>
#include <world/worldexc.h>
#include <world/worldfwd.h>

namespace madness {

    /// A collection of processes

    /// \c Group is a light weight object that can be used to specify a set of
    /// processes that will participate in Gop collective operations. The
    /// advantage of Group over MPI (or SafeMPI) groups is that it eliminates
    /// the need to construct new communicator and the associated barrier.
    /// \note This is \b NOT an MPI or SafeMPI group.
    class Group {
    private:

      class Impl {
      private:
        World& world_; ///< Parent world for this group
        DistributedID did_; ///< Group id
        std::vector<ProcessID> group_to_world_map_; ///< List of nodes in the group
        ProcessID group_rank_; ///< The group rank of this process
        volatile bool registered_; ///< Flag that is true when the group is in the register

        /// Array begin iterator accessor

        /// \tparam T The array type
        /// \tparam N The size of the array
        /// \param a A c-style array
        /// \return A pointer to the first element of \c a
        template <typename T, std::size_t N>
        static T* begin(T (&a)[N]) { return a; }

        /// Array end iterator accessor

        /// \tparam T The array type
        /// \tparam N The size of the array
        /// \param a A c-style array
        /// \return A pointer to one past the last element of \c a
        template <typename T, std::size_t N>
        static T* end(T (&a)[N]) { return (a + N); }

        /// Array size accessor

        /// \tparam T The array type
        /// \tparam N The size of the array
        /// \param a A c-style array
        /// \return The size of array \c a
        template <typename T, std::size_t N>
        static std::size_t size(T (&)[N]) { return N; }

        /// Array const begin iterator accessor

        /// \tparam T The array type
        /// \param a An array object
        /// \return The begin const_iterator of \c a
        template <typename vectorT>
        static typename vectorT::const_iterator begin(const vectorT &v) {
            return v.begin();
        }

        /// Array const end iterator accessor

        /// \tparam T The array type
        /// \param a An array object
        /// \return The end cosnt_iterator of \c a
        template <typename vectorT>
        static typename vectorT::const_iterator end(const vectorT &v) {
            return v.end();
        }

        /// Array size accessor

        /// \tparam T The array type
        /// \param a An array object
        /// \return The size of array \c a
        template <typename vectorT>
        static typename disable_if<std::is_array<vectorT>, std::size_t>::type
        size(const vectorT &v) { return v.size(); }

      public:
        /// Constructor

        /// \tparam A An std compliant array (e.g. \c std::array or <tt>std::vector</tt>)
        /// \param world The world that is the basis for this group
        /// \param did The distributed id associated with this group
        /// \param group An array of Processes in world
        template <typename A>
        Impl(World& world, const A& group, const DistributedID& did) :
          world_(world), did_(did), group_to_world_map_(begin(group), end(group)),
          group_rank_(-1), registered_(false)
        {
          // Check that there is at least one process in group
          MADNESS_ASSERT(size(group) > 0ul);

          // Sort and remove duplicates from group
          std::sort(group_to_world_map_.begin(), group_to_world_map_.end());
          group_to_world_map_.erase(std::unique(group_to_world_map_.begin(),
              group_to_world_map_.end()), group_to_world_map_.end());

          // Check that all processes in the group map are contained by world
          MADNESS_ASSERT(group_to_world_map_.front() >= 0);
          MADNESS_ASSERT(group_to_world_map_.back() < world_.size());

          // Get the group rank for this process
          group_rank_ = rank(world_.rank());

          // Check that this process is in the group
          MADNESS_ASSERT(group_rank_ != group_to_world_map_.size());
        }

        /// Parent world accessor

        /// \return A reference to the parent world of this group
        World& get_world() const { return world_; }

        /// Group id accessor

        /// \return A const reference to the group id
        const DistributedID& id() const { return did_; }

        /// Set the registered status

        /// \param status The new registered status
        void set_register_status(const bool status) { registered_ = status; }

        /// Registration status query

        /// \return \c true if the group has been registered, otherwise \c false.
        bool is_registered() const { return registered_; }

        /// Group rank accessor

        /// \return The rank of this process in the group
        ProcessID rank() const { return group_rank_; }

        /// Map world rank to group rank

        /// \param world_rank The world rank to be mapped
        /// \return The group rank of \c world_rank when it is a member of this
        /// group, otherwise \c -1.
        ProcessID rank(const ProcessID world_rank) const {
          ProcessID result = std::distance(group_to_world_map_.begin(),
              std::find(group_to_world_map_.begin(), group_to_world_map_.end(),
              world_rank));
          if(static_cast<std::size_t>(result) == group_to_world_map_.size())
            result = -1;
          return result;
        }

        /// Group size accessor

        /// \return The number of processes in the group
        ProcessID size() const { return group_to_world_map_.size(); }

        /// Map group rank to world rank

        /// \return The rank of this process in the world
        ProcessID world_rank(const ProcessID group_rank) const {
          MADNESS_ASSERT(group_rank >= 0);
          MADNESS_ASSERT(group_rank < ProcessID(group_to_world_map_.size()));
          return group_to_world_map_[group_rank];
        }

       /// Compute the binary tree parent and children

       /// \param[out] parent The parent node of the binary tree
       /// \param[out] child1 The left child node of the binary tree
       /// \param[out] child2 The right child node of the binary tree
       /// \param[in] group_root The head node of the binary tree
       void make_tree(const ProcessID group_root, ProcessID& parent,
           ProcessID& child0, ProcessID& child1) const
       {
         const ProcessID group_size = group_to_world_map_.size();

         // Check that root is in the range of the group
         MADNESS_ASSERT(group_root >= 0);
         MADNESS_ASSERT(group_root < group_size);

         // Renumber processes so root has me == 0
         const ProcessID me = (group_rank_ + group_size - group_root) % group_size;

         // Compute the group parent
         parent = (me == 0 ? -1 : group_to_world_map_[(((me - 1) >> 1) + group_root) % group_size]);

         // Compute children
         child0 = (me << 1) + 1 + group_root;
         child1 = child0 + 1;

         const ProcessID end = group_size + group_root;
         if(child0 < end)
           child0 = group_to_world_map_[child0 % group_size];
         else
           child0 = -1;
         if(child1 < end)
           child1 = group_to_world_map_[child1 % group_size];
         else
           child1 = -1;
       }

      }; // struct Impl

      class UnregisterGroup {
      private:
        DistributedID did_;

      public:
        UnregisterGroup() : did_(madness::uniqueidT(), 0ul) { }

        UnregisterGroup(const DistributedID& did) : did_(did) { }

        UnregisterGroup(const UnregisterGroup& other) : did_(other.did_) { }

        UnregisterGroup& operator=(const UnregisterGroup& other) {
          did_ = other.did_;
          return *this;
        }

        void operator()() const;

      }; // class UnregisterGroup


      std::shared_ptr<Impl> pimpl_;

    public:

      /// Default constructor

      /// Create an empty group
      Group() : pimpl_() { }

      /// Copy constructor

      /// \param other The group to be copied
      /// \note Copy is shallow.
      Group(const Group& other) : pimpl_(other.pimpl_) { }

      /// Create a new group

      /// \tparam A An array type
      /// \param world The parent world for this group
      /// \param group An array with a list of process to be included in the
      /// \param did The distributed id associated with this group
      /// \note All processes in the \c group list must be included in the
      /// parent world.
      template <typename A>
      Group(World& world, const A& group, const DistributedID& did) :
        pimpl_(new Impl(world, group, did))
      { }

      /// Create a new group

      /// \tparam A An array type
      /// \param world The parent world for this group
      /// \param group An array with a list of process to be included in the
      /// \param tag The tag associated with this group
      /// group.
      /// \note All processes in the \c group list must be included in the
      /// parent world.
      template <typename A>
      Group(World& world, const A& group, const std::size_t tag) :
        pimpl_(new Impl(world, group, DistributedID(uniqueidT(), tag)))
      { }

      /// Create a new group

      /// \tparam A An array type
      /// \param world The parent world for this group
      /// \param group An array with a list of process to be included in the
      /// \param uid The unique id (used by \c WorldObject ) associated with this group
      /// \param tag The tag associated with this group
      /// group.
      /// \note All processes in the \c group list must be included in the
      /// parent world.
      template <typename A>
      Group(World& world, const A& group, const uniqueidT& uid, const std::size_t tag) :
        pimpl_(new Impl(world, group, DistributedID(uid, tag)))
      { }

      /// Copy assignment operator

      /// \param other The group to be copied
      /// \note Copy is shallow.
      Group& operator=(const Group& other) {
        pimpl_ = other.pimpl_;
        return *this;
      }

      /// Register a group

      /// Register a group so that it can be used in active messages and tasks
      /// spawned on remote nodes.
      /// \param group The group to be registered
      /// \throw TiledArray::Exception When the group is empty
      /// \throw TiledArray::Exception When the group is already in the registry
      void register_group() const;

      /// Remove the given group from the registry

      /// Groups are removed via a lazy sync operation, which will only remove the
      /// group from the registry once \c unregister_group() has been called on
      /// all processes in the group.
      void unregister_group() const;

      /// Get a registered group

      /// This function is used to acquire the group in an active message handler.
      /// \param did The id associated with the group
      /// \return A future to the group
      static madness::Future<Group> get_group(const DistributedID& did);


      /// Registration status query

      /// \return \c true if the group has been registered, otherwise \c false.
      bool is_registered() const { return (pimpl_ && pimpl_->is_registered()); }

      /// Quary empty group

      /// \return \c true when this group is empty
      bool empty() const { return !pimpl_; }

      /// Group id accessor

      /// \return A const reference to the group id
      const DistributedID& id() const {
        MADNESS_ASSERT(pimpl_);
        return pimpl_->id();
      }

      /// Parent world accessor

      /// \return A reference to the parent world of this group
      World& get_world() const {
        MADNESS_ASSERT(pimpl_);
        return pimpl_->get_world();
      }

      /// Group rank accessor

      /// \return The rank of this process in the group
      ProcessID rank() const {
        MADNESS_ASSERT(pimpl_);
        return pimpl_->rank();
      }

      /// Map world rank to group rank

      /// \param world_rank The world rank to be mapped
      /// \return The rank of \c world_rank process in the group
      ProcessID rank(const ProcessID world_rank) const {
        MADNESS_ASSERT(pimpl_);
        return pimpl_->rank(world_rank);
      }

      /// Group size accessor

      /// \return The number of processes in the group
      ProcessID size() const {
        return (pimpl_ ? pimpl_->size() : 0);
      }

      /// Map group rank to world rank

      /// \param group_rank The group rank to be mapped to a world rank
      /// \return The parent world rank of group_rank.
      ProcessID world_rank(const ProcessID group_rank) const {
        MADNESS_ASSERT(pimpl_);
        return pimpl_->world_rank(group_rank);
      }

      /// Compute the binary tree parents and children

      /// \param[in] root The head node of the binary tree in the group
      /// \param[out] parent The parent node of the binary tree
      /// \param[out] child1 The left child node of the binary tree
      /// \param[out] child2 The right child node of the binary tree
      /// \note Output ranks are in the parent world.
      void make_tree(const ProcessID group_root, ProcessID& parent,
          ProcessID& child1, ProcessID& child2) const
      {
        MADNESS_ASSERT(pimpl_);
        pimpl_->make_tree(group_root, parent, child1, child2);
      }

      template <typename Archive>
      void serialize(const Archive&) {
        MADNESS_ASSERT(false); // not allowed
      }
    }; // class Group

} // namespace madness

#endif // MADNESS_WORLD_GROUP_H__INCLUDED
