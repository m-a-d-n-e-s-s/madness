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

#ifndef MADNESS_WORLD_WORLDGOP_H__INCLUDED
#define MADNESS_WORLD_WORLDGOP_H__INCLUDED

/// \file worldgop.h
/// \brief Implements global operations

/// If you can recall the Intel hypercubes, their comm lib used GOP as
/// the abbreviation.

#include <functional>
#include <type_traits>
#include <madness/world/worldtypes.h>
#include <madness/world/buffer_archive.h>
#include <madness/world/world.h>
#include <madness/world/deferred_cleanup.h>
#include <madness/world/world_task_queue.h>
#include <madness/world/group.h>
#include <madness/world/dist_cache.h>
#include <madness/world/units.h>

namespace madness {

    // Forward declarations
    class World;
    class WorldAmInterface;
    class WorldTaskQueue;
    namespace detail {

        class DeferredCleanup;

    }  // namespace detail

    template <typename T>
    struct WorldSumOp {
        inline T operator()(const T& a, const T& b) const {
            return a+b;
        }
    };

    template <typename T>
    struct WorldMultOp {
        inline T operator()(const T& a, const T& b) const {
            return a*b;
        }
    };

    template <typename T>
    struct WorldMaxOp {
        inline T operator()(const T& a, const T& b) const {
            return a>b? a : b;
        }
    };

    template <typename T>
    struct WorldAbsMaxOp {
        inline T operator()(const T& a, const T& b) const {
            return abs(a)>abs(b)? abs(a) : abs(b);
        }
    };

    template <typename T>
    struct WorldMinOp {
        inline T operator()(const T& a, const T& b) const {
            return a<b? a : b;
        }
    };

    template <typename T>
    struct WorldAbsMinOp {
        inline T operator()(const T& a, const T& b) const {
            return abs(a)<abs(b)? abs(a) : abs(b);
        }
    };

    template <typename T>
    struct WorldBitAndOp {
        inline T operator()(const T& a, const T& b) const {
            return a & b;
        }
    };

    template <typename T>
    struct WorldBitOrOp {
        inline T operator()(const T& a, const T& b) const {
            return a | b;
        }
    };

    template <typename T>
    struct WorldBitXorOp {
        inline T operator()(const T& a, const T& b) const {
            return a ^ b;
        }
    };

    template <typename T>
    struct WorldLogicAndOp {
        inline T operator()(const T& a, const T& b) const {
            return a && b;
        }
    };

    template <typename T>
    struct WorldLogicOrOp {
        inline T operator()(const T& a, const T& b) const {
            return a || b;
        }
    };


    /// Provides collectives that interoperate with the AM and task interfaces

    /// If native AM interoperates with MPI we probably should map these to MPI.
    class WorldGopInterface {
    private:
        World& world_; ///< World object that this is a part of
        std::shared_ptr<detail::DeferredCleanup> deferred_; ///< Deferred cleanup object.
        bool debug_; ///< Debug mode
        bool forbid_fence_=false; ///< forbid calling fence() in case of several active worlds
        int max_reducebcast_msg_size_ = std::numeric_limits<int>::max();  ///< maximum size of messages (in bytes) sent by reduce and broadcast

        friend class detail::DeferredCleanup;

        // Message tags
        struct PointToPointTag { };
        struct LazySyncTag { };
        struct GroupLazySyncTag { };
        struct BcastTag { };
        struct GroupBcastTag { };
        struct ReduceTag { };
        struct GroupReduceTag { };
        struct AllReduceTag { };
        struct GroupAllReduceTag { };


        /// Delayed send callback object

        /// This callback object is used to send local data to a remove process
        /// once it has been set.
        /// \tparam keyT The data key
        /// \tparam valueT The type of data to be sent
        template <typename keyT, typename valueT>
        class DelayedSend : public CallbackInterface {
        private:
            World& world_; ///< The communication world
            const ProcessID dest_; ///< The destination process id
            const keyT key_; ///< The distributed id associated with \c value_
            Future<valueT> value_; ///< The data to be sent

            // Not allowed
            DelayedSend(const DelayedSend<keyT, valueT>&);
            DelayedSend<keyT, valueT>& operator=(const DelayedSend<keyT, valueT>&);

        public:

            /// Constructor
            DelayedSend(World& world, const ProcessID dest,
                    const keyT& key, const Future<valueT>& value) :
                        world_(world), dest_(dest), key_(key), value_(value)
            { }

            virtual ~DelayedSend() { }

            /// Notify this object that the future has been set.

            /// This will set the value of the future on the remote node and delete
            /// this callback object.
            virtual void notify() {
                MADNESS_ASSERT(value_.probe());
                world_.gop.send_internal(dest_, key_, value_.get());
                delete this;
            }
        }; // class DelayedSend

        /// Receive data from remote node

        /// \tparam valueT The data type stored in cache
        /// \param key The distributed ID
        /// \return A future to the data
        template <typename valueT, typename keyT>
        static Future<valueT> recv_internal(const keyT& key) {
            return detail::DistCache<keyT>::template get_cache_value<valueT>(key);
        }

        /// Send \c value to \c dest

        /// Send non-future data to \c dest.
        /// \tparam keyT The key type
        /// \tparam valueT The value type
        /// \param dest The node where the data will be sent
        /// \param key The key that is associated with the data
        /// \param value The data to be sent to \c dest
        template <typename keyT, typename valueT>
        typename std::enable_if<!is_future<valueT>::value >::type
        send_internal(const ProcessID dest, const keyT& key, const valueT& value) const {
            typedef detail::DistCache<keyT> dist_cache;

            if(world_.rank() == dest) {
                // When dest is this process, skip the task and set the future immediately.
                dist_cache::set_cache_value(key, value);
            } else {
                // Spawn a remote task to set the value
                world_.taskq.add(dest, dist_cache::template set_cache_value<valueT>, key,
                        value, TaskAttributes::hipri());
            }
        }

        /// Send \c value to \c dest

        /// Send data that is stored in a future to \c dest. The data in
        /// \c value is only sent to the remote process once it has been set.
        /// \tparam keyT The key type
        /// \tparam valueT The value type
        /// \param dest The node where the data will be sent
        /// \param key The key that is associated with the data
        /// \param value The data to be sent to \c dest
        template <typename keyT, typename valueT>
        void send_internal(ProcessID dest, const keyT& key, const Future<valueT>& value) const {
            typedef detail::DistCache<keyT> dist_cache;

            if(world_.rank() == dest) {
                dist_cache::set_cache_value(key, value);
            } else {
                // The destination is not this node, so send it to the destination.
                if(value.probe()) {
                    // Spawn a remote task to set the value
                    world_.taskq.add(dest, dist_cache::template set_cache_value<valueT>, key,
                            value.get(), TaskAttributes::hipri());
                } else {
                    // The future is not ready, so create a callback object that will
                    // send value to the destination node when it is ready.
                    DelayedSend<keyT, valueT>* delayed_send_callback =
                            new DelayedSend<keyT, valueT>(world_, dest, key, value);
                    const_cast<Future<valueT>&>(value).register_callback(delayed_send_callback);

                }
            }
        }

        /// Lazy sync parent task

        /// Send signal to the parent process in the binary tree for a lazy sync
        /// operation.
        /// \tparam keyT The key type
        /// \param parent The parent process of this process in the binary tree
        /// \param key The lazy sync key
        template <typename keyT>
        void lazy_sync_parent(const ProcessID parent, const keyT& key,
                const ProcessID, const ProcessID) const
        {
            send_internal(parent, key, key.proc());
        }

        /// Lazy sync parent task

        /// Send signal to the child processes in the binary tree for a lazy
        /// sync operation. After the signal has been sent to the children, the
        /// sync operation, \c op, will be run.
        /// \tparam keyT The key type
        /// \tparam opT The sync operation type
        /// \param child0 The first child process of this process in the binary tree
        /// \param child1 The second child process of this process in the binary tree
        /// \param key The key associated with the sync operation
        /// \param op The sync operation that will be run
        template <typename keyT, typename opT>
        void lazy_sync_children(const ProcessID child0, const ProcessID child1,
                const keyT& key, opT& op, const ProcessID) const
        {
            // Signal children to execute the operation.
            if(child0 != -1)
                send_internal(child0, key, 1);
            if(child1 != -1)
                send_internal(child1, key, 1);

            // Execute the operation on this process.
            op();
        }

        /// Start a distributed lazy sync operation

        /// \param key The sync key
        /// \param op The sync operation to be executed on this process
        template <typename tagT, typename keyT, typename opT>
        void lazy_sync_internal(const ProcessID parent, const ProcessID child0,
                const ProcessID child1, const keyT& key, const opT& op) const {
            typedef ProcessKey<keyT, tagT> key_type;

            // Get signals from parent and children.
            madness::Future<ProcessID> child0_signal = (child0 != -1 ?
                    recv_internal<ProcessID>(key_type(key, child0)) :
                    madness::Future<ProcessID>(-1));
            madness::Future<ProcessID> child1_signal = (child1 != -1 ?
                    recv_internal<ProcessID>(key_type(key, child1)) :
                    madness::Future<ProcessID>(-1));
            madness::Future<ProcessID> parent_signal = (parent != -1 ?
                    recv_internal<ProcessID>(key_type(key, parent)) :
                    madness::Future<ProcessID>(-1));

            // Construct the task that notifies children to run the operation
            key_type my_key(key, world_.rank());
            auto lazy_sync_children_fn = & WorldGopInterface::template lazy_sync_children<key_type, opT>;
            world_.taskq.add(*this, lazy_sync_children_fn,
                    child0_signal, child1_signal, my_key, op, parent_signal,
                    TaskAttributes::hipri());

            // Send signal to parent
            if(parent != -1) {
                if(child0_signal.probe() && child1_signal.probe())
                    send_internal(parent, my_key, world_.rank());
                else {
                    auto lazy_sync_parent_fn = & WorldGopInterface::template lazy_sync_parent<key_type>;
                    world_.taskq.add(*this, lazy_sync_parent_fn,
                            parent, my_key, child0_signal, child1_signal,
                            TaskAttributes::hipri());
                }
            }
        }


        template <typename keyT, typename valueT, typename taskfnT>
        static void bcast_handler(const AmArg& arg) {
            // Deserialize message arguments
            taskfnT taskfn;
            keyT key;
            valueT value;
            ProcessID root;

            arg & taskfn & key & value & root;

            // Add task to queue
            arg.get_world()->taskq.add(arg.get_world()->gop, taskfn, key,
                    value, root, TaskAttributes::hipri());
        }

        template <typename keyT, typename valueT, typename taskfnT>
        static void group_bcast_handler(const AmArg& arg) {
            // Deserialize message arguments
            taskfnT taskfn;
            keyT key;
            valueT value;
            ProcessID group_root;
            DistributedID group_key;

            arg & taskfn & key & value & group_root & group_key;

            // Get the local group
            const Future<Group> group = Group::get_group(group_key);

            // Add task to queue
            arg.get_world()->taskq.add(arg.get_world()->gop, taskfn, key, value,
                    group_root, group, TaskAttributes::hipri());
        }


        /// Broadcast task

        /// This task will set the local cache with the broadcast data and send
        /// it to child processes in the binary tree.
        template <typename keyT, typename valueT>
        void bcast_task(const keyT& key, const valueT& value, const ProcessID root) const {
            typedef void (WorldGopInterface::*taskfnT)(const keyT&, const valueT&,
                    const ProcessID) const;

            // Compute binary tree data
            ProcessID parent = -1, child0 = -1, child1 = -1;
            world_.mpi.binary_tree_info(root, parent, child0, child1);

            // Set the local data, except on the root process
            if(parent != -1)
                detail::DistCache<keyT>::set_cache_value(key, value);

            if(child0 != -1) { // Check that this process has children in the binary tree

                // Get handler function and arguments
                void (*handler)(const AmArg&) =
                        & WorldGopInterface::template bcast_handler<keyT, valueT, taskfnT>;
                AmArg* const args0 = new_am_arg(
                        & WorldGopInterface::template bcast_task<keyT, valueT>,
                        key, value, root);

                // Send active message to children
                if(child1 != -1) {
                    AmArg* const args1 = copy_am_arg(*args0);
                    world_.am.send(child1, handler, args1, RMI::ATTR_UNORDERED);
                }
                world_.am.send(child0, handler, args0, RMI::ATTR_UNORDERED);
            }
        }

        template <typename keyT, typename valueT>
        void group_bcast_task(const keyT& key, const valueT& value,
                const ProcessID group_root, const Group& group) const
        {
            typedef void (WorldGopInterface::*taskfnT)(const keyT&, const valueT&,
                    const ProcessID, const Group&) const;

            // Set the local data, except on the root process
            ProcessID parent = -1, child0 = -1, child1 = -1;
            group.make_tree(group_root, parent, child0, child1);

            // Set the local data
            if(parent != -1) {
                detail::DistCache<keyT>::set_cache_value(key, value);
                group.remote_update();
            }

            if(child0 != -1) { // Check that this process has children in the binary tree

                // Get handler function and arguments
                void (*handler)(const AmArg&) =
                        & WorldGopInterface::template group_bcast_handler<keyT, valueT, taskfnT>;
                AmArg* const args0 = new_am_arg(
                        & WorldGopInterface::template group_bcast_task<keyT, valueT>,
                        key, value, group_root, group.id());

                // Send active message to children
                if(child1 != -1) {
                    AmArg* const args1 = copy_am_arg(*args0);
                    world_.am.send(child1, handler, args1, RMI::ATTR_UNORDERED);
                }
                world_.am.send(child0, handler, args0, RMI::ATTR_UNORDERED);
            }
        }

        /// Broadcast

        /// Broadcast data from the \c root process to all processes in \c world.
        /// The input/output data is held by \c value.
        /// \tparam tagT The tag type that is attached to \c keyT
        /// \tparam keyT The base key type
        /// \tparam valueT The value type that will be broadcast
        /// \param[in] key The key associated with this broadcast
        /// \param[in,out] value On the \c root process, this is used as the input
        /// data that will be broadcast to all other processes in the group.
        /// On other processes it is used as the output to the broadcast
        /// \param root The process that owns the data to be broadcast
        /// \throw madness::Exception When \c value has been set, except on the
        /// \c root process.
        template <typename tagT, typename keyT, typename valueT>
        void bcast_internal(const keyT& key, Future<valueT>& value, const ProcessID root) const {
            MADNESS_ASSERT((root >= 0) && (root < world_.size()));
            MADNESS_ASSERT((world_.rank() == root) || (! value.probe()));

            // Add operation tag to key
            typedef TaggedKey<keyT, tagT> key_type;
            const key_type tagged_key(key);

            if(world_.size() > 1) { // Do nothing for the trivial case
                if(world_.rank() == root) {
                    // This process owns the data to be broadcast.

                    // Spawn remote tasks that will set the local cache for this
                    // broadcast on other nodes.
                    if(value.probe())
                        // The value is ready so send it now
                        bcast_task(tagged_key, value.get(), root);
                    else {
                        // The value is not ready so spawn a task to send the
                        // data when it is ready.
                        auto bcast_task_fn = & WorldGopInterface::template bcast_task<key_type, valueT>;
                        world_.taskq.add(*this, bcast_task_fn,
                                tagged_key, value, root, TaskAttributes::hipri());
                    }
                } else {
                    MADNESS_ASSERT(! value.probe());

                    // Get the broadcast value from local cache
                    detail::DistCache<key_type>::get_cache_value(tagged_key, value);
                }
            }
        }

        /// Group broadcast

        /// Broadcast data from the \c group_root process to all processes in
        /// \c group. The input/output data is held by \c value.
        /// \tparam tagT The tag type that is attached to \c keyT
        /// \tparam keyT The base key type
        /// \tparam valueT The value type that will be broadcast
        /// \param[in] key The key associated with this broadcast
        /// \param[in,out] value On the \c group_root process, this is used as the
        /// input data that will be broadcast to all other processes in the group.
        /// On other processes it is used as the output to the broadcast
        /// \param group_root The process in \c group that owns the data to be
        /// broadcast
        /// \param group The process group where value will be broadcast
        /// \throw madness::Exception When \c value has been set, except on the
        /// \c group_root process.
        template <typename tagT, typename keyT, typename valueT>
        void bcast_internal(const keyT& key, Future<valueT>& value,
                const ProcessID group_root, const Group& group) const
        {
            // Construct the internal broadcast key
            typedef TaggedKey<keyT, tagT> key_type;
            const key_type tagged_key(key);

            if(group.rank() == group_root) {
                // This process owns the data to be broadcast.
                if(value.probe())
                    group_bcast_task(tagged_key, value.get(), group_root, group);
                else {
                    auto group_bcast_task_fn = & WorldGopInterface::template group_bcast_task<key_type, valueT>;
                    world_.taskq.add(this, group_bcast_task_fn,
                            tagged_key, value, group_root, group,
                            TaskAttributes::hipri());
                }
            } else {
                MADNESS_ASSERT(! value.probe());

                // This is not the root process, so retrieve the broadcast data
                detail::DistCache<key_type>::get_cache_value(tagged_key, value);

                // Increment local use counter for group
                group.local_update();
            }
        }

        template <typename valueT, typename opT>
        static typename detail::result_of<opT>::type
        reduce_task(const valueT& value, const opT& op) {
            typename detail::result_of<opT>::type result = op();
            op(result, value);
            return result;
        }

        template <typename opT>
        static typename detail::result_of<opT>::type
        reduce_result_task(const std::vector<Future<typename detail::result_of<opT>::type> >& results,
                const opT& op)
        {
            MADNESS_ASSERT(results.size() != 0ul);
            Future<typename detail::result_of<opT>::type> result = results.front();
            for(std::size_t i = 1ul; i < results.size(); ++i)
                op(result.get(), results[i].get());
            return result.get();
        }

        /// Distributed reduce

        /// \tparam tagT The tag type to be added to the key type
        /// \tparam keyT The key type
        /// \tparam valueT The data type to be reduced
        /// \tparam opT The reduction operation type
        /// \param key The key associated with this reduction
        /// \param value The local value to be reduced
        /// \param op The reduction operation to be applied to local and remote data
        /// \param root The process that will receive the result of the reduction
        /// \return A future to the reduce value on the root process, otherwise an
        /// uninitialized future that may be ignored.
        template <typename tagT, typename keyT, typename valueT, typename opT>
        Future<typename detail::result_of<opT>::type>
        reduce_internal(const ProcessID parent, const ProcessID child0,
                const ProcessID child1, const ProcessID root, const keyT& key,
                const valueT& value, const opT& op)
        {
            // Create tagged key
            typedef ProcessKey<keyT, tagT> key_type;
            typedef typename detail::result_of<opT>::type result_type;
            typedef typename remove_future<valueT>::type value_type;
            std::vector<Future<result_type> > results;
            results.reserve(3);

            // Add local data to vector of values to reduce
            results.push_back(world_.taskq.add(WorldGopInterface::template reduce_task<value_type, opT>,
                    value, op, TaskAttributes::hipri()));

            // Reduce child data
            if(child0 != -1)
                results.push_back(recv_internal<result_type>(key_type(key, child0)));
            if(child1 != -1)
                results.push_back(recv_internal<result_type>(key_type(key, child1)));

            // Submit the local reduction task
            Future<result_type> local_result =
                    world_.taskq.add(WorldGopInterface::template reduce_result_task<opT>,
                            results, op, TaskAttributes::hipri());

            // Send reduced value to parent or, if this is the root process, set the
            // result future.
            if(parent == -1)
                return local_result;
            else
                send_internal(parent, key_type(key, world_.rank()), local_result);

            return Future<result_type>::default_initializer();
        }

        /// Implementation of fence

        /// \param[in] epilogue the action to execute (by the calling thread) immediately after the fence
        /// \param[in] pause_during_epilogue whether to suspend work while executing epilogue
        /// \param[in] debug set to true to print progress statistics using madness::print(); the default is false.
        /// \warning currently only \c pause_during_epilogue=false is supported
        void fence_impl(std::function<void()> epilogue = []{},
                        bool pause_during_epilogue = false,
                        bool debug = false);

        int initial_max_reducebcast_msg_size() {
          int result = std::numeric_limits<int>::max();
          const auto* initial_max_reducebcast_msg_size_cstr = std::getenv("MAD_MAX_REDUCEBCAST_MSG_SIZE");
          if (initial_max_reducebcast_msg_size_cstr) {
            auto result_u64 = cstr_to_memory_size(initial_max_reducebcast_msg_size_cstr);
            const auto do_print = SafeMPI::COMM_WORLD.Get_rank() == 0 && !madness::quiet();
            if (result_u64>std::numeric_limits<int>::max()) {
              if (do_print)
                std::cout
                    << "!!MADNESS WARNING: Invalid value for environment variable MAD_MAX_REDUCEBCAST_MSG_SIZE.\n"
                    << "!!MADNESS WARNING: MAD_MAX_REDUCEBCAST_MSG_SIZE = "
                    << result_u64 << "\n";
              result = std::numeric_limits<int>::max();
            }
            result = static_cast<int>(result_u64);
            if(do_print) {
              std::cout
                  << "MADNESS max msg size for GOP reduce/broadcast set to "
                  << result << " bytes.\n";
            }
          }
          return result;
        }

    public:

        // In the World constructor can ONLY rely on MPI and MPI being initialized
        WorldGopInterface(World& world) :
            world_(world), deferred_(new detail::DeferredCleanup()), debug_(false), max_reducebcast_msg_size_(initial_max_reducebcast_msg_size())
        { }

        ~WorldGopInterface() {
            deferred_->destroy(true);
            deferred_->do_cleanup();
        }


        /// Set debug flag to new value and return old value
        bool set_debug(bool value) {
            bool status = debug_;
            debug_ = value;
            return status;
        }

        /// Set forbid_fence flag to new value and return old value
        bool set_forbid_fence(bool value) {
            bool status = forbid_fence_;
            forbid_fence_ = value;
            return status;
        }

        /// Set the maximum size of messages (in bytes) sent by reduce and broadcast

        /// \param sz the maximum size of messages (in bytes) sent by reduce and broadcast
        /// \return the previous maximum size of messages (in bytes) sent by reduce and broadcast
        /// \pre `sz>0`
        int set_max_reducebcast_msg_size(int sz) {
          MADNESS_ASSERT(sz>0);
          std::swap(max_reducebcast_msg_size_,sz);
          return max_reducebcast_msg_size_;
        }


        /// Returns the maximum size of messages (in bytes) sent by reduce and broadcast

        /// \return the maximum size of messages (in bytes) sent by reduce and broadcast
        int max_reducebcast_msg_size() const {
          return max_reducebcast_msg_size_;
        }

        /// Synchronizes all processes in communicator ... does NOT fence pending AM or tasks
        void barrier() {
            long i = world_.rank();
            sum(i);
            if (i != world_.size()*(world_.size()-1)/2) error("bad value after sum in barrier");
        }


        /// Synchronizes all processes in communicator AND globally ensures no pending AM or tasks

        /// \internal Runs Dykstra-like termination algorithm on binary tree by
        /// locally ensuring ntask=0 and all am sent and processed,
        /// and then participating in a global sum of nsent and nrecv.
        /// Then globally checks that nsent=nrecv and that both are
        /// constant over two traversals.  We are then sure
        /// that all tasks and AM are processed and there no AM in
        /// flight.
        /// \param[in] debug set to true to print progress statistics using madness::print(); the default is false.
        void fence(bool debug = false);

        /// Executes an action on single (this) thread after ensuring all other work is done

        /// \param[in] action the action to execute (by the calling thread)
        void serial_invoke(std::function<void()> action);

        /// Broadcasts bytes from process root while still processing AM & tasks

        /// Optimizations can be added for long messages
        void broadcast(void* buf, size_t nbyte, ProcessID root, bool dowork = true, Tag bcast_tag = -1);


        /// Broadcasts typed contiguous data from process root while still processing AM & tasks

        /// Optimizations can be added for long messages
        template <typename T, typename = std::enable_if_t<madness::is_trivially_copyable_v<T>>>
        inline void broadcast(T* buf, size_t nelem, ProcessID root) {
            broadcast((void *) buf, nelem*sizeof(T), root);
        }

        /// Broadcast of a scalar from node 0 to all other nodes
        template <typename T, typename = std::enable_if_t<madness::is_trivially_copyable_v<T>>>
        void broadcast(T& t) {
            broadcast(&t, 1, 0);
        }

        /// Broadcast of a scalar from node root to all other nodes
        template <typename T, typename = std::enable_if_t<madness::is_trivially_copyable_v<T>>>
        void broadcast(T& t, ProcessID root) {
            broadcast(&t, 1, root);
        }

        /// Broadcast a serializable object
        template <typename objT,
                  typename = std::void_t<decltype(std::declval<archive::BufferInputArchive&>()&std::declval<objT&>())>,
                  typename = std::void_t<decltype(std::declval<archive::BufferOutputArchive&>()&std::declval<const objT&>())>>
        void broadcast_serializable(objT& obj, ProcessID root) {
            MADNESS_ASSERT(root < world_.size());
            if (world_.size() == 1) return;

            size_t BUFLEN;
            if (world_.rank() == root) {
                archive::BufferOutputArchive count;
                count & obj;
                BUFLEN = count.size();
            }
            broadcast(BUFLEN, root);

            unsigned char* buf = new unsigned char[BUFLEN];
            if (world_.rank() == root) {
                archive::BufferOutputArchive ar(buf,BUFLEN);
                ar & obj;
            }
            broadcast(buf, BUFLEN, root);
            if (world_.rank() != root) {
                archive::BufferInputArchive ar(buf,BUFLEN);
                ar & obj;
            }
            delete [] buf;
        }

        /// Inplace global reduction (like MPI all_reduce) while still processing AM & tasks

        /// Optimizations can be added for long messages and to reduce the memory footprint
        template <typename T, class opT>
            void reduce(T* buf, std::size_t nelem, opT op) {
          static_assert(madness::is_trivially_copyable_v<T>, "T must be trivially copyable");

          ProcessID parent, child0, child1;
          world_.mpi.binary_tree_info(0, parent, child0, child1);
          const std::size_t nelem_per_maxmsg =
              max_reducebcast_msg_size() / sizeof(T);

          const auto buf_size = ((sizeof(T) * std::min(nelem_per_maxmsg, nelem) +
                                 std::alignment_of_v<T> - 1) /
                                std::alignment_of_v<T>) * std::alignment_of_v<T>;
          struct free_dtor {
            void operator()(T *ptr) {
              if (ptr != nullptr)
                std::free(ptr);
            };
          };
          using sptr_t = std::unique_ptr<T[], free_dtor>;

          auto aligned_buf_alloc = [&]() -> T* {
            // posix_memalign requires alignment to be an integer multiple of sizeof(void*)!! so ensure that
            const std::size_t alignment =
                ((std::alignment_of_v<T> + sizeof(void *) - 1) /
                 sizeof(void *)) *
                sizeof(void *);
#ifdef HAVE_POSIX_MEMALIGN
            void *ptr;
            if (posix_memalign(&ptr, alignment, buf_size) != 0) {
              throw std::bad_alloc();
            }
            return static_cast<T *>(ptr);
#else
            return static_cast<T *>(std::aligned_alloc(alignment, buf_size));
#endif
          };

          sptr_t buf0;
          if (child0 != -1)
            buf0 = sptr_t(aligned_buf_alloc(),
                          free_dtor{});
          sptr_t buf1(nullptr);
          if (child1 != -1)
            buf1 = sptr_t(aligned_buf_alloc(),
                          free_dtor{});

          auto reduce_impl = [&,this](T* buf, size_t nelem) {
            MADNESS_ASSERT(nelem <= nelem_per_maxmsg);
            SafeMPI::Request req0, req1;
            Tag gsum_tag = world_.mpi.unique_tag();

            if (child0 != -1)
              req0 = world_.mpi.Irecv(buf0.get(), nelem * sizeof(T), MPI_BYTE,
                                      child0, gsum_tag);
            if (child1 != -1)
              req1 = world_.mpi.Irecv(buf1.get(), nelem * sizeof(T), MPI_BYTE,
                                      child1, gsum_tag);

            if (child0 != -1) {
              World::await(req0);
              for (long i = 0; i < (long)nelem; ++i)
                buf[i] = op(buf[i], buf0[i]);
            }
            if (child1 != -1) {
              World::await(req1);
              for (long i = 0; i < (long)nelem; ++i)
                buf[i] = op(buf[i], buf1[i]);
            }

            if (parent != -1) {
              req0 = world_.mpi.Isend(buf, nelem * sizeof(T), MPI_BYTE, parent,
                                      gsum_tag);
              World::await(req0);
            }

            broadcast(buf, nelem, 0);
          };

          while (nelem) {
            const int n = std::min(nelem_per_maxmsg, nelem);
            reduce_impl(buf, n);
            nelem -= n;
            buf += n;
          }
        }

        /// Inplace global sum while still processing AM & tasks
        template <typename T>
        inline void sum(T* buf, size_t nelem) {
            reduce< T, WorldSumOp<T> >(buf, nelem, WorldSumOp<T>());
        }

        /// Inplace global min while still processing AM & tasks
        template <typename T>
        inline void min(T* buf, size_t nelem) {
            reduce< T, WorldMinOp<T> >(buf, nelem, WorldMinOp<T>());
        }

        /// Inplace global max while still processing AM & tasks
        template <typename T>
        inline void max(T* buf, size_t nelem) {
            reduce< T, WorldMaxOp<T> >(buf, nelem, WorldMaxOp<T>());
        }

        /// Inplace global absmin while still processing AM & tasks
        template <typename T>
        inline void absmin(T* buf, size_t nelem) {
            reduce< T, WorldAbsMinOp<T> >(buf, nelem, WorldAbsMinOp<T>());
        }

        /// Inplace global absmax while still processing AM & tasks
        template <typename T>
        inline void absmax(T* buf, size_t nelem) {
            reduce< T, WorldAbsMaxOp<T> >(buf, nelem, WorldAbsMaxOp<T>());
        }

        /// Inplace global product while still processing AM & tasks
        template <typename T>
        inline void product(T* buf, size_t nelem) {
            reduce< T, WorldMultOp<T> >(buf, nelem, WorldMultOp<T>());
        }

        template <typename T>
        inline void bit_and(T* buf, size_t nelem) {
            reduce< T, WorldBitAndOp<T> >(buf, nelem, WorldBitAndOp<T>());
        }

        template <typename T>
        inline void bit_or(T* buf, size_t nelem) {
            reduce< T, WorldBitOrOp<T> >(buf, nelem, WorldBitOrOp<T>());
        }

        template <typename T>
        inline void bit_xor(T* buf, size_t nelem) {
            reduce< T, WorldBitXorOp<T> >(buf, nelem, WorldBitXorOp<T>());
        }

        template <typename T>
        inline void logic_and(T* buf, size_t nelem) {
            reduce< T, WorldLogicAndOp<T> >(buf, nelem, WorldLogicAndOp<T>());
        }

        template <typename T>
        inline void logic_or(T* buf, size_t nelem) {
            reduce< T, WorldLogicOrOp<T> >(buf, nelem, WorldLogicOrOp<T>());
        }

        /// Global sum of a scalar while still processing AM & tasks
        template <typename T>
        void sum(T& a) {
            sum(&a, 1);
        }

        /// Global max of a scalar while still processing AM & tasks
        template <typename T>
        void max(T& a) {
            max(&a, 1);
        }

        /// Global min of a scalar while still processing AM & tasks
        template <typename T>
        void min(T& a) {
            min(&a, 1);
        }

        /// Concatenate an STL vector of serializable stuff onto node 0

        /// \param[in] v input vector
        /// \param[in] bufsz the max number of bytes in the result; must be less than std::numeric_limits<int>::max()
        /// \return on rank 0 returns the concatenated vector, elsewhere returns an empty vector
        template <typename T>
        std::vector<T> concat0(const std::vector<T>& v, size_t bufsz=1024*1024) {
            MADNESS_ASSERT(bufsz <= std::numeric_limits<int>::max());
            // bufsz must be multiple of alignment!!! so ensure that
            bufsz = ((bufsz + sizeof(void*) - 1) / sizeof(void*)) * sizeof(void*);

            ProcessID parent, child0, child1;
            world_.mpi.binary_tree_info(0, parent, child0, child1);
            int child0_nbatch = 0, child1_nbatch = 0;

            struct free_dtor {
              void operator()(std::byte *ptr) {
                if (ptr != nullptr)
                  std::free(ptr);
              };
            };
            using sptr_t = std::unique_ptr<std::byte[], free_dtor>;

            auto aligned_buf_alloc = [&]() -> std::byte* {
#ifdef HAVE_POSIX_MEMALIGN
                void *ptr;
                if (posix_memalign(&ptr, sizeof(void *), bufsz) != 0) {
                    throw std::bad_alloc();
                }
                return static_cast<std::byte *>(ptr);
#else
                return static_cast<std::byte *>(
                    std::aligned_alloc(sizeof(void *), bufsz));
#endif
            };

            auto buf0 = sptr_t(aligned_buf_alloc(),
                               free_dtor{});
            auto buf1 = sptr_t(aligned_buf_alloc(),
                               free_dtor{});

            // transfer data in chunks at most this large
            const int batch_size = static_cast<int>(
                std::min(static_cast<size_t>(max_reducebcast_msg_size()), bufsz));

            // precompute max # of tags any node ... will need, and allocate them on every node to avoid tag counter divergence
            const int max_nbatch = bufsz / batch_size;
            // one tag is reserved for sending the number of messages to expect and the size of the last message
            const int max_ntags = max_nbatch + 1;
            MADNESS_ASSERT(max_nbatch < world_.mpi.unique_tag_period());
            std::vector<Tag> tags;  // stores tags used to send each batch
            tags.reserve(max_nbatch);
            for(int t=0; t<max_ntags; ++t) tags.push_back(world_.mpi.unique_tag());

            if (child0 != -1 || child1 != -1) {
              // receive # of batches

              auto receive_nbatch = [&,this]() {
                if (child0 != -1) {
                  world_.mpi.Recv(&child0_nbatch, 1, MPI_INT, child0,
                                          tags[0]);
                }
                if (child1 != -1) {
                  world_.mpi.Recv(&child1_nbatch, 1, MPI_INT, child1,
                                          tags[0]);
                }
              };

              receive_nbatch();

              // receive data in batches

              auto receive_batch = [&,this](const int batch, const size_t buf_offset) {
                SafeMPI::Request req0, req1;
                if (child0 != -1 && batch < child0_nbatch) {
                  int msg_size = batch_size;
                  // if last batch, receive # of bytes to expect
                  if (batch + 1 == child0_nbatch) {
                    auto req = world_.mpi.Irecv(
                        &msg_size, 1, MPI_INT, child0, tags[0]);
                    World::await(req);
                  }

                  req0 = world_.mpi.Irecv(buf0.get() + buf_offset,
                                          msg_size, MPI_BYTE, child0,
                                          tags[batch + 1]);
                }
                if (child1 != -1 && batch < child1_nbatch) {
                  int msg_size = batch_size;
                  // if last batch, receive # of bytes to expect
                  if (batch + 1 == child1_nbatch) {
                    auto req = world_.mpi.Irecv(
                        &msg_size, 1, MPI_INT, child1, tags[0]);
                    World::await(req);
                  }
                  req1 = world_.mpi.Irecv(buf1.get() + buf_offset,
                                          msg_size, MPI_BYTE, child1,
                                          tags[batch + 1]);
                }

                if (child0 != -1 && batch < child0_nbatch) {
                  World::await(req0);
                }
                if (child1 != -1 && batch < child1_nbatch) {
                  World::await(req1);
                }
              };

              size_t buf_offset = 0;
              int batch = 0;
              while (buf_offset < bufsz) {
                receive_batch(batch, buf_offset);
                buf_offset += batch_size;
                buf_offset = std::min(buf_offset, bufsz);
                ++batch;
              }
            }


            std::vector<T> left, right;
            if (child0 != -1) {
              archive::BufferInputArchive ar(buf0.get(), bufsz);
              ar & left;
            }
            if (child1 != -1) {
              archive::BufferInputArchive ar(buf1.get(), bufsz);
              ar & right;
              for (unsigned int i = 0; i < right.size(); ++i)
                left.push_back(right[i]);
            }
            for (unsigned int i=0; i<v.size(); ++i) left.push_back(v[i]);

            // send data in batches
            if (parent != -1) {
              archive::BufferOutputArchive ar(buf0.get(), bufsz);
              ar & left;
              const auto total_nbytes_to_send = ar.size();

              // send nbatches to expect
              const int nbatch = (total_nbytes_to_send + batch_size - 1) / batch_size;
              world_.mpi.Send(&nbatch, 1, MPI_INT, parent,
                              tags[0]);

              size_t buf_offset = 0;
              int batch = 0;
              while (buf_offset < bufsz) {

                // send data in batches
                auto send_batch = [&,this](const int batch, const size_t buf_offset) {
                  const int nbytes_to_send = static_cast<int>(
                      std::min(static_cast<size_t>(batch_size),
                               total_nbytes_to_send - buf_offset));
                  // if last batch, send # of bytes to expect
                  if (batch + 1 == nbatch) {
                    auto req = world_.mpi.Isend(
                        &nbytes_to_send, 1, MPI_INT, parent, tags[0]);
                    World::await(req);
                  }
                  auto req0 =
                      world_.mpi.Isend(buf0.get() + buf_offset, nbytes_to_send,
                                       MPI_BYTE, parent, tags[batch + 1]);
                  World::await(req0);
                };

                send_batch(batch, buf_offset);
                buf_offset += batch_size;
                buf_offset = std::min(buf_offset, bufsz);
                ++batch;
              }
            }

            if (parent == -1) return left;
            else return std::vector<T>();
        }

        /// Receive data from \c source

        /// \tparam valueT The data type stored in cache
        /// \tparam keyT The key type
        /// \param source The process that is sending the data to this process
        /// \param key The key associated with the received data
        /// \return A future that will be set with the received data
        /// \note It is the user's responsibility to ensure that \c key does not
        /// conflict with other calls to \c recv. Keys may be reused after the
        /// associated operation has finished.
        template <typename valueT, typename keyT>
        static Future<valueT> recv(const ProcessID source, const keyT& key) {
            return recv_internal<valueT>(ProcessKey<keyT, PointToPointTag>(key, source));
        }

        /// Send value to \c dest

        /// \tparam keyT The key type
        /// \tparam valueT The value type (this may be a \c Future type)
        /// \param dest The process where the data will be sent
        /// \param key The key that is associated with the data
        /// \param value The data to be sent to \c dest
        /// \note It is the user's responsibility to ensure that \c key does not
        /// conflict with other calls to \c send. Keys may be reused after the
        /// associated operation has finished.
        template <typename keyT, typename valueT>
        void send(const ProcessID dest, const keyT& key, const valueT& value) const {
            send_internal(dest, ProcessKey<keyT, PointToPointTag>(key, world_.rank()), value);
        }

        /// Lazy sync

        /// Lazy sync functions are asynchronous barriers with a nullary functor
        /// that is called after all processes have called it with the same
        /// key. You can think of lazy_sync as an asynchronous barrier. The
        /// lazy_sync functor must have the following signature:
        /// \code
        /// class SyncFunc {
        /// public:
        ///     // typedefs
        ///     typedef void result_type;
        ///
        ///     // Constructors
        ///     SyncFunc(const SyncFunc&);
        ///
        ///     // The function that performs the sync operation
        ///     void operator()();
        ///
        /// }; // class SyncFunc
        /// \endcode
        /// \tparam keyT The key type
        /// \tparam opT The operation type
        /// \param key The sync key
        /// \param op The sync operation to be executed on this process
        /// \note It is the user's responsibility to ensure that \c key does not
        /// conflict with other calls to \c lazy_sync. Keys may be reused after
        /// the associated operation has finished.
        template <typename keyT, typename opT>
        void lazy_sync(const keyT& key, const opT& op) const {
            if(world_.size() > 1) { // Do nothing for the trivial case
                // Get the binary tree data
                Hash<keyT> hasher;
                const ProcessID root = hasher(key) % world_.size();
                ProcessID parent = -1, child0 = -1, child1 = -1;
                world_.mpi.binary_tree_info(root, parent, child0, child1);

                lazy_sync_internal<LazySyncTag>(parent, child0, child1, key, op);
            } else {
                auto lazy_sync_children_fn = & WorldGopInterface::template lazy_sync_children<keyT, opT>;
                // There is only one process, so run the sync operation now.
                world_.taskq.add(*this, lazy_sync_children_fn,
                        -1, -1, key, op, -1, TaskAttributes::hipri());
            }
        }


        /// Group lazy sync

        /// Lazy sync functions are asynchronous barriers with a nullary functor
        /// that is called after all processes in the group have called it with
        /// the same key. You can think of lazy_sync as an asynchronous barrier.
        /// The \c op functor must have the following signature:
        /// \code
        /// class SyncFunc {
        /// public:
        ///     // typedefs
        ///     typedef void result_type;
        ///
        ///     // Constructors
        ///     SyncFunc(const SyncFunc&);
        ///
        ///     // The function that performs the sync operation
        ///     void operator()();
        ///
        /// }; // class SyncFunc
        /// \endcode
        /// \tparam keyT The key type
        /// \tparam opT The operation type
        /// \param key The sync key
        /// \param op The sync operation to be executed on this process
        /// \note It is the user's responsibility to ensure that \c key does not
        /// conflict with other calls to \c lazy_sync. Keys may be reused after
        /// the associated operation has finished.
        template <typename keyT, typename opT>
        void lazy_sync(const keyT& key, const opT& op, const Group& group) const {
            MADNESS_ASSERT(! group.empty());
            MADNESS_ASSERT(group.get_world().id() == world_.id());

            if(group.size() > 1) { // Do nothing for the trivial case
                // Get the binary tree data
                Hash<keyT> hasher;
                const ProcessID group_root = hasher(key) % group.size();
                ProcessID parent = -1, child0 = -1, child1 = -1;
                group.make_tree(group_root, parent, child0, child1);

                lazy_sync_internal<GroupLazySyncTag>(parent, child0, child1, key, op);
            } else {
                auto lazy_sync_children_fn = & WorldGopInterface::template lazy_sync_children<keyT, opT>;
                world_.taskq.add(*this, lazy_sync_children_fn,
                        -1, -1, key, op, -1, TaskAttributes::hipri());
            }
        }

        /// Broadcast

        /// Broadcast data from the \c root process to all processes. The input/
        /// output data is held by \c value.
        /// \param[in] key The key associated with this broadcast
        /// \param[in,out] value On the \c root process, this is used as the input
        /// data that will be broadcast to all other processes. On other
        /// processes it is used as the output to the broadcast.
        /// \param root The process that owns the data to be broadcast
        /// \throw madness::Exception When \c root is less than 0 or greater
        /// than or equal to the world size.
        /// \throw madness::Exception When \c value has been set, except on the
        /// \c root process.
        /// \note It is the user's responsibility to ensure that \c key does not
        /// conflict with other calls to \c bcast. Keys may be reused after
        /// the associated operation has finished.
        template <typename keyT, typename valueT>
        void bcast(const keyT& key, Future<valueT>& value, const ProcessID root) const {
            MADNESS_ASSERT((root >= 0) && (root < world_.size()));
            MADNESS_ASSERT((world_.rank() == root) || (! value.probe()));

            if(world_.size() > 1) // Do nothing for the trivial case
                bcast_internal<BcastTag>(key, value, root);
        }

        /// Group broadcast

        /// Broadcast data from the \c group_root process to all processes in
        /// \c group. The input/output data is held by \c value.
        /// \param[in] key The key associated with this broadcast
        /// \param[in,out] value On the \c group_root process, this is used as the
        /// input data that will be broadcast to all other processes in the group.
        /// On other processes it is used as the output to the broadcast
        /// \param group_root The process in \c group that owns the data to be
        /// broadcast
        /// \param group The process group where value will be broadcast
        /// \throw madness::Exception When \c group is empty
        /// \throw madness::Exception When \c group is not registered
        /// \throw madness::Exception When the world id of \c group is not
        /// equal to that of the world used to construct this object
        /// \throw madness::Exception When this process is not in the group
        /// \throw madness::Exception When \c group_root is less than 0 or
        /// greater than or equal to \c group size
        /// \throw madness::Exception When \c data has been set except on the
        /// \c root process
        /// \note It is the user's responsibility to ensure that \c key does not
        /// conflict with other calls to \c bcast. Keys may be reused after
        /// the associated operation has finished.
        template <typename keyT, typename valueT>
        void bcast(const keyT& key, Future<valueT>& value,
                const ProcessID group_root, const Group& group) const
        {
            MADNESS_ASSERT(! group.empty());
            MADNESS_ASSERT(group.get_world().id() == world_.id());
            MADNESS_ASSERT((group_root >= 0) && (group_root < group.size()));
            MADNESS_ASSERT((group.rank() == group_root) || (! value.probe()));

            if(group.size() > 1) // Do nothing for the trivial case
                bcast_internal<GroupBcastTag>(key, value, group_root, group);
        }

        /// Distributed reduce

        /// The reduce functor must have the following signature:
        /// \code
        /// class ReduceFunc {
        /// public:
        ///     // Typedefs
        ///     typedef ... result_type;
        ///     typedef ... argument_type;
        ///
        ///     // Constructors
        ///     ReduceFunc(const ReduceFunc&);
        ///
        ///     // Initialization operation, which returns a default result object
        ///     result_type operator()() const;
        ///
        ///     // Reduce two result objects
        ///     void operator()(result_type&, const result_type&) const;
        ///
        ///     // Reduce a result object and an argument object
        ///     void operator()(result_type&, const argument_type&) const;
        /// }; // class ReduceFunc
        /// \endcode
        /// \tparam keyT The key type
        /// \tparam valueT The data type to be reduced
        /// \tparam opT The reduction operation type
        /// \param key The key associated with this reduction
        /// \param value The local value to be reduced
        /// \param op The reduction operation to be applied to local and remote data
        /// \param root The process that will receive the result of the reduction
        /// \return A future to the reduce value on the root process, otherwise an
        /// uninitialized future that may be ignored.
        /// \note It is the user's responsibility to ensure that \c key does not
        /// conflict with other calls to \c reduce. Keys may be reused after
        /// the associated operation has finished.
        template <typename keyT, typename valueT, typename opT>
        Future<typename detail::result_of<opT>::type>
        reduce(const keyT& key, const valueT& value, const opT& op, const ProcessID root) {
            MADNESS_ASSERT((root >= 0) && (root < world_.size()));

            // Get the binary tree data
            ProcessID parent = -1, child0 = -1, child1 = -1;
            world_.mpi.binary_tree_info(root, parent, child0, child1);

            return reduce_internal<ReduceTag>(parent, child0, child1, root, key,
                    value, op);
        }

        /// Distributed group reduce

        /// The reduce functor must have the following signature:
        /// \code
        /// class ReduceFunc {
        /// public:
        ///     // Typedefs
        ///     typedef ... result_type;
        ///     typedef ... argument_type;
        ///
        ///     // Constructors
        ///     ReduceFunc(const ReduceFunc&);
        ///
        ///     // Initialization operation, which returns a default result object
        ///     result_type operator()() const;
        ///
        ///     // Reduce two result objects
        ///     void operator()(result_type&, const result_type&) const;
        ///
        ///     // Reduce a result object and an argument object
        ///     void operator()(result_type&, const argument_type&) const;
        /// }; // class ReduceFunc
        /// \endcode
        /// \tparam keyT The key type
        /// \tparam valueT The data type to be reduced
        /// \tparam opT The reduction operation type
        /// \param key The key associated with this reduction
        /// \param value The local value to be reduced
        /// \param op The reduction operation to be applied to local and remote data
        /// \param group_root The group process that will receive the result of the reduction
        /// \param group The group that will preform the reduction
        /// \return A future to the reduce value on the root process, otherwise an
        /// uninitialized future that may be ignored.
        /// \throw madness::Exception When \c group is empty
        /// \throw madness::Exception When \c group is not registered
        /// \throw madness::Exception When the world id of \c group is not
        /// equal to that of the world used to construct this object
        /// \throw madness::Exception When this process is not in the group
        /// \throw madness::Exception When \c group_root is less than zero or
        /// greater than or equal to \c group size.
        /// \note It is the user's responsibility to ensure that \c key does not
        /// conflict with other calls to \c reduce. Keys may be reused after
        /// the associated operation has finished.
        template <typename keyT, typename valueT, typename opT>
        Future<typename detail::result_of<opT>::type>
        reduce(const keyT& key, const valueT& value, const opT& op,
                const ProcessID group_root, const Group& group)
        {
            MADNESS_ASSERT(! group.empty());
            MADNESS_ASSERT(group.get_world().id() == world_.id());
            MADNESS_ASSERT((group_root >= 0) && (group_root < group.size()));

            // Get the binary tree data
            ProcessID parent = -1, child0 = -1, child1 = -1;
            group.make_tree(group_root, parent, child0, child1);

            return reduce_internal<ReduceTag>(parent, child0, child1, group_root,
                    key, value, op);
        }

        /// Distributed all reduce

        /// The reduce functor must have the following signature:
        /// \code
        /// class ReduceFunc {
        /// public:
        ///     // Typedefs
        ///     typedef ... result_type;
        ///     typedef ... argument_type;
        ///
        ///     // Constructors
        ///     ReduceFunc(const ReduceFunc&);
        ///
        ///     // Initialization operation, which returns a default result object
        ///     result_type operator()() const;
        ///
        ///     // Reduce two result objects
        ///     void operator()(result_type&, const result_type&) const;
        ///
        ///     // Reduce a result object and an argument object
        ///     void operator()(result_type&, const argument_type&) const;
        /// }; // class ReduceFunc
        /// \endcode
        /// \tparam keyT The key type
        /// \tparam valueT The data type to be reduced
        /// \tparam opT The reduction operation type
        /// \param key The key associated with this reduction
        /// \param value The local value to be reduced
        /// \param op The reduction operation to be applied to local and remote data
        /// \return A future to the reduce value on the root process, otherwise an
        /// uninitialized future that may be ignored.
        /// \note It is the user's responsibility to ensure that \c key does not
        /// conflict with other calls to \c all_reduce. Keys may be reused after
        /// the associated operation has finished.
        template <typename keyT, typename valueT, typename opT>
        Future<typename detail::result_of<opT>::type>
        all_reduce(const keyT& key, const valueT& value, const opT& op) {
            // Compute the parent and child processes of this process in a binary tree.
            Hash<keyT> hasher;
            const ProcessID root = hasher(key) % world_.size();
            ProcessID parent = -1, child0 = -1, child1 = -1;
            world_.mpi.binary_tree_info(root, parent, child0, child1);

            // Reduce the data
            Future<typename detail::result_of<opT>::type> reduce_result =
                    reduce_internal<AllReduceTag>(parent, child0, child1, root,
                            key, value, op);

            if(world_.rank() != root)
                reduce_result = Future<typename detail::result_of<opT>::type>();

            // Broadcast the result of the reduction to all processes
            bcast_internal<AllReduceTag>(key, reduce_result, root);

            return reduce_result;
        }

        /// Distributed, group all reduce

        /// The reduce functor must have the following signature:
        /// \code
        /// class ReduceFunc {
        /// public:
        ///     // Typedefs
        ///     typedef ... result_type;
        ///     typedef ... argument_type;
        ///
        ///     // Constructors
        ///     ReduceFunc(const ReduceFunc&);
        ///
        ///     // Initialization operation, which returns a default result object
        ///     result_type operator()() const;
        ///
        ///     // Reduce two result objects
        ///     void operator()(result_type&, const result_type&) const;
        ///
        ///     // Reduce a result object and an argument object
        ///     void operator()(result_type&, const argument_type&) const;
        /// }; // class ReduceFunc
        /// \endcode
        /// \tparam keyT The key type
        /// \tparam valueT The data type to be reduced
        /// \tparam opT The reduction operation type
        /// \param key The key associated with this reduction
        /// \param value The local value to be reduced
        /// \param op The reduction operation to be applied to local and remote data
        /// \param group The group that will preform the reduction
        /// \return A future to the reduce value on the root process, otherwise an
        /// uninitialized future that may be ignored
        /// \throw madness::Exception When \c group is empty
        /// \throw madness::Exception When \c group is not registered
        /// \throw madness::Exception When the world id of \c group is not
        /// equal to that of the world used to construct this object
        /// \throw madness::Exception When this process is not in the group
        /// \note It is the user's responsibility to ensure that \c key does not
        /// conflict with other calls to \c reduce. Keys may be reused after
        /// the associated operation has finished.
        template <typename keyT, typename valueT, typename opT>
        Future<typename detail::result_of<opT>::type>
        all_reduce(const keyT& key, const valueT& value, const opT& op, const Group& group) {
            MADNESS_ASSERT(! group.empty());
            MADNESS_ASSERT(group.get_world().id() == world_.id());

            // Compute the parent and child processes of this process in a binary tree.
            Hash<keyT> hasher;
            const ProcessID group_root = hasher(key) % group.size();
            ProcessID parent = -1, child0 = -1, child1 = -1;
            group.make_tree(group_root, parent, child0, child1);

            // Reduce the data
            Future<typename detail::result_of<opT>::type> reduce_result =
                    reduce_internal<GroupAllReduceTag>(parent, child0, child1,
                            group_root, key, value, op);


            if(group.rank() != group_root)
                reduce_result = Future<typename detail::result_of<opT>::type>();

            // Broadcast the result of the reduction to all processes in the group
            bcast_internal<GroupAllReduceTag>(key, reduce_result, 0, group);

            return reduce_result;
        }
    }; // class WorldGopInterface

} // namespace madness

#endif // MADNESS_WORLD_WORLDGOP_H__INCLUDED
