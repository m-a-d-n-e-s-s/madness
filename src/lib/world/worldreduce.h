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

#ifndef MADNESS_WOLRD_REDUCTION_H__INCLUDED
#define MADNESS_WOLRD_REDUCTION_H__INCLUDED

#include <world/worldtypes.h>
#include <world/worldobj.h>
#include <world/worldtask.h>
#include <algorithm>
#include <iterator>

namespace madness {

    namespace detail {

        /// Interface class for GroupReduction objects
        class ReductionInterface {
        private:

            unsigned int count_;    ///< The number of reductions that will be done on this node
            ProcessID parent_;      ///< This node's parent process for remote reduction
            ProcessID child0_;      ///< This node's children process for remote reduction
            ProcessID child1_;      ///< This node's children process for remote reduction

            // Not allowed
            ReductionInterface(const ReductionInterface&);
            ReductionInterface& operator=(const ReductionInterface&);

        public:

            /// Construct the reduction object base
            ReductionInterface() :
                count_(0), parent_(-1), child0_(-1), child1_(-1)
            { }

            virtual ~ReductionInterface() { }

            /// Define the group members

            /// \tparam initerT Input iterator type that dereferences to ProcessID
            /// \param world_rank The rank of this process in the world.
            /// \param first An input iterator that points to a list of processes
            /// that are included in the group.
            /// \param last An input
            template <typename initerT>
            void set_group(const ProcessID world_rank, initerT first, initerT last, ProcessID root) {
                // This function should only be called once so count will be 0.
                MADNESS_ASSERT(count_ == 0);

                // Count for this process
                count_ = 1;

                // Sort and remove duplicate elements from process map
                std::vector<ProcessID> group(first, last);

                std::sort(group.begin(), group.end());
                std::unique(group.begin(), group.end());

                // Get the size of the group
                const int size = group.size();

                // Get the group rank of this process
                const ProcessID rank = std::distance(group.begin(),
                        std::find_if(group.begin(), group.end(),
                        std::bind1st(std::equal_to<ProcessID>(), world_rank)));
                MADNESS_ASSERT(rank != size);

                // Get the root of the binary tree
                if(root != -1)
                    root = std::distance(group.begin(), std::find_if(group.begin(),
                        group.end(), std::bind1st(std::equal_to<ProcessID>(), root)));
                else
                    root = group.front();
                MADNESS_ASSERT(root != size);

                // Renumber processes so root has me=0
                int me = (rank + size-root)%size;

                // Parent in binary tree
                if(me != 0) {
                    parent_ = (((me - 1) >> 1) + root) % size;
                    parent_ = group[parent_];
                } else
                    parent_ = -1;

                // Left child
                child0_ = (me << 1) + 1 + root;
                if (child0_ >= size && child0_<(size+root))
                    child0_ -= size;
                if (child0_ < size) {
                    child0_ = group[child0_];
                    ++count_;
                } else
                    child0_ = -1;

                // Right child
                child1_ = (me << 1) + 2 + root;
                if (child1_ >= size && child1_<(size+root))
                    child1_ -= size;
                if (child1_ < size) {
                    child1_ = group[child1_];
                    ++count_;
                } else
                    child1_ = -1;
            }

            ProcessID parent() const { return parent_; }
            ProcessID child0() const { return child0_; }
            ProcessID child1() const { return child1_; }
            bool is_root() const { return parent_ == -1; }
            unsigned int count() const { return count_; }

            template <typename T>
            void reduce(const T& value) {
                MADNESS_ASSERT(this->type() == typeid(T));
                this->reduce_value(reinterpret_cast<const void*>(&value));
            }

            template <typename T>
            void reduce(const Future<T>& fut) {
                MADNESS_ASSERT(this->type() == typeid(T));
                this->reduce_future(reinterpret_cast<const void*>(&fut));
            }

        private:

            virtual void reduce_value(const void* p) = 0;
            virtual void reduce_future(const void* p) = 0;
            virtual const std::type_info& type() const = 0;

        }; // class ReductionInterface

        template <typename T>
        class GroupReduction : public ReductionInterface {
        private:
            Future<T> r0_;
            Future<T> r1_;
            Future<T> r2_;
            Future<T> result_;

        public:

            GroupReduction() { }

            virtual ~GroupReduction() { }

            template <typename opT>
            Future<T> set_op(World& w, const opT& op) {
                MADNESS_ASSERT(ReductionInterface::count() != 0);

                TaskAttributes attr;
                attr.set_highpriority(true);
                switch(ReductionInterface::count()) {
                case 3:
                    {
                        Future<T> temp = w.taskq.add(op, r0_, r1_, attr);
                        result_ = w.taskq.add(op, temp, r2_, attr);
                    }
                    break;

                case 2:
                    result_ = w.taskq.add(op, r0_, r1_, attr);
                    break;

                case 1:
                    result_.set(r0_);
                    break;
                }

                return result_;
            }

        private:

            virtual void reduce_value(const void* p) {
                const T* v = reinterpret_cast<const T*>(p);
                set_future(*v);
            }

            virtual void reduce_future(const void * p) {
                const Future<T>* f = reinterpret_cast<const Future<T>* >(p);
                set_future(*f);
            }

            template <typename U>
            void set_future(const U& u) {
                if(! r0_.probe())
                    r0_.set(u);
                else if(! r1_.probe())
                    r1_.set(u);
                else if(! r2_.probe())
                    r2_.set(u);
                else
                    MADNESS_EXCEPTION("All reduction values have already been set.", false);
            }

            virtual const std::type_info& type() const { return typeid(T); }

        }; // class GroupReduction

    }  // namespace detail

    /// WorldObject with group reduction functionality.

    /// \tparam derivedT The derived class type.
    /// \tparam keyT Reduction group key type.
    /// \note This is still a \c WorldObject , which means all the rules for
    /// \c WorldObjects apply to this object. Most notably, the order of
    /// object construction must be the same on all nodes and the derived
    /// class must call \c WorldObject::process_pending() .
    template <typename derivedT, typename keyT = std::size_t>
    class WorldReduce : public WorldObject<derivedT> {
    public:
        typedef keyT key;                           ///< The key type used to track reductions

    private:
        typedef WorldObject<derivedT> WorldObject_; ///< The base class type
        typedef WorldReduce<derivedT> WorldReduce_; ///< This class type
        typedef ConcurrentHashMap<keyT, std::shared_ptr<detail::ReductionInterface> >
                reduce_container;                   ///< The container type that holds reductions

        reduce_container reductions_;               ///< Stores reduction objects


        template <typename valueT>
        Void reduce_value(const key& k, const valueT& value) {
            typename reduce_container::accessor acc;
            insert<valueT>(acc, k);
            // The accessor ensures the data is write locked.
            acc->second->reduce(value);
            return None;
        }

        template <typename valueT>
        Void send_to_parent(ProcessID parent, const key& k, const valueT& value) {
            if(parent != -1)
                WorldObject_::send(parent, & WorldReduce_::template reduce_value<valueT>, k, value);

            // We do not need the local reduction anymore
            reductions_.erase(k);
            return None;
        }

        template <typename valueT>
        void insert(typename reduce_container::accessor& acc, const key& k) {
            std::shared_ptr<detail::GroupReduction<valueT> > reduction;
            if(reductions_.insert(acc, k)) {
                reduction.reset(new detail::GroupReduction<valueT>());
                acc->second = reduction;
            }
        }

    public:

        /// Construct a world reduction object

        /// \param world The world where this reducer lives.
        /// \note Derived constructors must call \c WorldObject::process_pending() .
        WorldReduce(World& world) :
            WorldObject_(world)
        { }

        virtual ~WorldReduce() { }

        /// Create a group reduction.

        /// The reduction groups are differentiated by key value. There fore it
        /// is the user's responsibility to ensure that the key is unique to the
        /// group. All reductions are guaranteed to be complete after a global
        /// fence. You may safely reuse group keys after a fence.
        /// \tparam valueT The type that will be reduced
        /// \tparam opT The reduction operation type
        /// \tparam initerT Input iterator type for list of process in group
        /// \param k The key that will be used for the group
        /// \param value The local value to be reduced
        /// \param first An iterator to the list of nodes that will be included
        /// in the group.
        /// \param last An iterator to the end of the list of group nodes
        /// \param root The root process [default = the node with the lowest
        /// rank in the node list]
        /// \return A future to the reduction results on the local node
        /// \note The return future will contain an incomplete reduction on all
        /// nodes other than the root.
        template <typename valueT, typename opT, typename initerT>
        Future<valueT> reduce(std::size_t k, const valueT& value, opT op, initerT first, initerT last, ProcessID root = -1) {

            // Make sure nodes in group are >= 0 and < w.size()
            MADNESS_ASSERT(std::find_if(first, last,
                    std::bind2nd(std::less<ProcessID>(), 0)) == last);
            MADNESS_ASSERT(std::find_if(first, last,
                    std::bind2nd(std::greater_equal<ProcessID>(),
                    WorldObject_::get_world().size())) == last);

            // Create/find the reduction object
            typename reduce_container::accessor acc;
            insert<valueT>(acc, k);

            // Set the reduction group
            acc->second->set_group(WorldObject_::get_world().rank(), first, last, root);

            // Set the reduction operation
            std::shared_ptr<detail::GroupReduction<valueT> > reduction =
                std::static_pointer_cast<detail::GroupReduction<valueT> >(acc->second);
            Future<valueT> result = reduction->set_op(WorldObject_::get_world(), op);

            // Forward the results of the local reductions to the parent
            if(! acc->second->is_root())
                WorldObject_::get_world().taskq.add(*this,
                    & WorldReduce_::template send_to_parent<valueT>,
                    acc->second->parent(), k, result);

            // Reduce the local value
            acc->second->reduce(value);

            return result;
        }
    }; // class WorldReduce

}  // namespace madness

#endif // MADNESS_WOLRD_REDUCTION_H__INCLUDED
