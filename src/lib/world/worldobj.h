/*
  This file is part of MADNESS.

  Copyright (C) <2007> <Oak Ridge National Laboratory>

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


#ifndef WORLD_OBJ_H
#define WORLD_OBJ_H

#include <world/worldthread.h>

namespace madness {

    namespace detail {

        // Common base class for pending messages to ensure in-order processing

        // To eliminate synchronization when a distributed object is first
        // constructed, we buffer pending messages for containers that
        // don't have their id yet registered.
        struct PendingMsg {
            uniqueidT id;
            am_handlerT handler;
            AmArg* arg;

            PendingMsg(uniqueidT id, am_handlerT handler, const AmArg& arg)
                    : id(id), handler(handler), arg(copy_am_arg(arg)) {}

            void invokehandler() {
                handler(*arg);
                free_am_arg(arg);
            };
        };

        // It is annoying that we must replicate the task forwarding stuff here but we must
        // so that pending messages creating tasks are correctly handled.  Cannot merge
        // easily with send handlers since task layer is more restrictive on the
        // copy capability of arguments.

        /// Info stored for AM method forwarding
        template <typename memfunT>
        struct info {
            typedef Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > futureT;
            typedef RemoteReference< FutureImpl< REMFUTURE(MEMFUN_RETURNT(memfunT)) > > refT;
            uniqueidT id; // Must be at front ... see peek.
            ProcessID requestor;
            memfunT memfun;
            refT ref;
            TaskAttributes attr;

            info() {};

            info(const AmArg& arg) {
                arg & *this;
            }

            info(const uniqueidT& id, ProcessID requestor,  memfunT memfun, const refT& ref,
                 const TaskAttributes& attr=TaskAttributes())
                    : id(id)
                    , requestor(requestor)
                    , memfun(memfun)
                    , ref(ref)
                    , attr(attr) {};

            template <typename Archive>
            void serialize(const Archive& ar) {
                ar & archive::wrap_opaque(*this); // Must be opaque ... see peek.
            }
        };

        /// Extract the unique object ID from an incoming active message header

        /// We deserialize the header and all arguments at the same
        /// time to simplify the code.  However, it is common that
        /// when sending a message to an item in a container to
        /// include a pointer to the container itself.  But this
        /// breaks if the container is not initialized since the
        /// deserialization throws if the object is not initialized
        /// (which seems preferable to hidden race condition).  Hence,
        /// we use this routine to extract the unique ID from the very
        /// front of the info structure.  For efficiency we here rely
        /// upon the serialization of info being opaque and the
        /// id being at the front of info.
        static inline const uniqueidT& peek(const AmArg& arg) {
            return *((uniqueidT*)(arg.buf()));
        }
    }


    /// Implements most parts of a globally addressable object (via unique ID)

    /// 1) Derived class has WorldObject<Derived> as a public base class
    /// 2) Derived constructor
    ///    a) invokes WorldObject<Derived>(world) constructor
    ///    b) invokes process_pending()
    /// 3) Derived destructor must either be deferred or preceeded by gop.fence()
    ///
    /// This class is deliberately not default constructible and does
    /// not support assignment or copying.  This ensures that each instance
    /// is unique.  Have a look at the WorldContainer for an example
    /// of wrapping this using the PIMPL idiom and a shared pointer.
    ///
    /// Note that world is exposed for convenience as a public data member.
    template <class Derived>
    class WorldObject : public DeferredCleanupInterface {
        typedef WorldObject<Derived> objT;
        typedef std::list<detail::PendingMsg> pendingT;
    public:
        World& world;                              ///< Think globally act locally
    private:
        // The order here matters in a multi-threaded world
        volatile bool ready;                       ///< True if ready to rock 'n roll
        ProcessID me;                              ///< Rank of self
        static Spinlock pending_mutex;
        static volatile pendingT pending;          ///< Buffers pending messages
        uniqueidT objid;                           ///< Sense of self

        // Handler for incoming AM with 0 arguments
        template <typename memfunT>
        static void handler(const AmArg& arg) {
            const uniqueidT& id = detail::peek(arg);
            Derived* obj = arg.get_world()-> template ptr_from_id<Derived>(id);
            am_handlerT ptr = handler<memfunT>;
            if (is_ready(id,obj,arg,ptr)) {
                detail::info<memfunT> info;
                arg & info;
                typename detail::info<memfunT>::futureT result(info.ref);
                result.set((obj->*info.memfun)());
            }
        }

        // Handler for incoming AM with 1 argument
        template <typename memfunT, typename arg1T>
        static void handler(const AmArg& arg) {
            const uniqueidT& id = detail::peek(arg);
            Derived* obj = arg.get_world()-> template ptr_from_id<Derived>(id);
            am_handlerT ptr = handler<memfunT,arg1T>;
            if (is_ready(id,obj,arg,ptr)) {
                detail::info<memfunT> info;
                arg1T arg1;
                arg & info & arg1;
                typename detail::info<memfunT>::futureT result(info.ref);
                result.set((obj->*info.memfun)(arg1));
            }
        }

        // Handler for incoming AM with 2 arguments
        template <typename memfunT, typename arg1T, typename arg2T>
        static void handler(const AmArg& arg) {
            const uniqueidT& id = detail::peek(arg);
            Derived* obj = arg.get_world()-> template ptr_from_id<Derived>(id);
            am_handlerT ptr = handler<memfunT,arg1T,arg2T>;
            if (is_ready(id,obj,arg,ptr)) {
                detail::info<memfunT> info;
                arg1T arg1;
                arg2T arg2;
                arg & info & arg1 & arg2;
                typename detail::info<memfunT>::futureT result(info.ref);
                result.set((obj->*info.memfun)(arg1,arg2));
            }
        }

        // Handler for incoming AM with 3 arguments
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        static void handler(const AmArg& arg) {
            const uniqueidT& id = detail::peek(arg);
            Derived* obj = arg.get_world()-> template ptr_from_id<Derived>(id);
            am_handlerT ptr = handler<memfunT,arg1T,arg2T,arg3T>;
            if (is_ready(id,obj,arg,ptr)) {
                detail::info<memfunT> info;
                arg1T arg1;
                arg2T arg2;
                arg3T arg3;
                arg & info & arg1 & arg2 & arg3;
                typename detail::info<memfunT>::futureT result(info.ref);
                result.set((obj->*info.memfun)(arg1,arg2,arg3));
            }
        }

        // Handler for incoming AM with 4 arguments
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        static void handler(const AmArg& arg) {
            const uniqueidT& id = detail::peek(arg);
            Derived* obj = arg.get_world()-> template ptr_from_id<Derived>(id);
            am_handlerT ptr = handler<memfunT,arg1T,arg2T,arg3T,arg4T>;
            if (is_ready(id,obj,arg,ptr)) {
                detail::info<memfunT> info;
                arg1T arg1;
                arg2T arg2;
                arg3T arg3;
                arg4T arg4;
                arg & info & arg1 & arg2 & arg3 & arg4;
                typename detail::info<memfunT>::futureT result(info.ref);
                result.set((obj->*info.memfun)(arg1,arg2,arg3,arg4));
            }
        }

        // Handler for incoming AM with 5 arguments
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T,
        typename arg5T>
        static void handler(const AmArg& arg) {
            const uniqueidT& id = detail::peek(arg);
            Derived* obj = arg.get_world()-> template ptr_from_id<Derived>(id);
            am_handlerT ptr = handler<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T>;
            if (is_ready(id,obj,arg,ptr)) {
                detail::info<memfunT> info;
                arg1T arg1;
                arg2T arg2;
                arg3T arg3;
                arg4T arg4;
                arg5T arg5;
                arg & info & arg1 & arg2 & arg3 & arg4 & arg5;
                typename detail::info<memfunT>::futureT result(info.ref);
                result.set((obj->*info.memfun)(arg1,arg2,arg3,arg4,arg5));
            }
        }

        // Handler for incoming AM with 6 arguments
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T,
        typename arg5T, typename arg6T>
        static void handler(const AmArg& arg) {
            const uniqueidT& id = detail::peek(arg);
            Derived* obj = arg.get_world()-> template ptr_from_id<Derived>(id);
            am_handlerT ptr = handler<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T,arg6T>;
            if (is_ready(id,obj,arg,ptr)) {
                detail::info<memfunT> info;
                arg1T arg1;
                arg2T arg2;
                arg3T arg3;
                arg4T arg4;
                arg5T arg5;
                arg6T arg6;
                arg & info & arg1 & arg2 & arg3 & arg4 & arg5 & arg6;
                typename detail::info<memfunT>::futureT result(info.ref);
                result.set((obj->*info.memfun)(arg1,arg2,arg3,arg4,arg5,arg6));
            }
        }

        // Handler for incoming AM with 7 arguments
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T,
        typename arg5T, typename arg6T, typename arg7T>
        static void handler(const AmArg& arg) {
            const uniqueidT& id = detail::peek(arg);
            Derived* obj = arg.get_world()-> template ptr_from_id<Derived>(id);
            am_handlerT ptr = handler<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T>;
            if (is_ready(id,obj,arg,ptr)) {
                detail::info<memfunT> info;
                arg1T arg1;
                arg2T arg2;
                arg3T arg3;
                arg4T arg4;
                arg5T arg5;
                arg6T arg6;
                arg7T arg7;
                arg & info & arg1 & arg2 & arg3 & arg4 & arg5 & arg6 & arg7;
                typename detail::info<memfunT>::futureT result(info.ref);
                result.set((obj->*info.memfun)(arg1,arg2,arg3,arg4,arg5,arg6,arg7));
            }
        }

        // Handler for incoming task with 0 arguments
        template <typename memfunT>
        static void handler_task(const AmArg& arg) {
            const uniqueidT& id = detail::peek(arg);
            Derived* obj = arg.get_world()-> template ptr_from_id<Derived>(id);
            am_handlerT ptr = handler_task<memfunT>;
            if (is_ready(id,obj,arg,ptr)) {
                detail::info<memfunT> info;
                arg & info;
                typename detail::info<memfunT>::futureT result(info.ref);
                arg.get_world()->taskq.add(new TaskMemfun<memfunT>(result,*obj,info.memfun,info.attr));
            }
        }

        // Handler for incoming task with 1 argument
        template <typename memfunT, typename arg1T>
        static void handler_task(const AmArg& arg) {
            const uniqueidT& id = detail::peek(arg);
            Derived* obj = arg.get_world()-> template ptr_from_id<Derived>(id);
            am_handlerT ptr = handler_task<memfunT,arg1T>;
            if (is_ready(id,obj,arg,ptr)) {
                detail::info<memfunT> info;
                arg1T arg1;
                arg & info & arg1;
                typename detail::info<memfunT>::futureT result(info.ref);
                arg.get_world()->taskq.add(new TaskMemfun<memfunT>(result,*obj,info.memfun,arg1,info.attr));
            }
        }

        // Handler for incoming task with 2 arguments
        template <typename memfunT, typename arg1T, typename arg2T>
        static void handler_task(const AmArg& arg) {
            const uniqueidT& id = detail::peek(arg);
            Derived* obj = arg.get_world()-> template ptr_from_id<Derived>(id);
            am_handlerT ptr = handler_task<memfunT,arg1T,arg2T>;
            if (is_ready(id,obj,arg,ptr)) {
                detail::info<memfunT> info;
                arg1T arg1;
                arg2T arg2;
                arg & info & arg1 & arg2;
                typename detail::info<memfunT>::futureT result(info.ref);
                arg.get_world()->taskq.add(new TaskMemfun<memfunT>(result,*obj,info.memfun,arg1,arg2,info.attr));
            }
        }

        // Handler for incoming task with 3 arguments
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        static void handler_task(const AmArg& arg) {
            const uniqueidT& id = detail::peek(arg);
            Derived* obj = arg.get_world()-> template ptr_from_id<Derived>(id);
            am_handlerT ptr = handler_task<memfunT,arg1T,arg2T,arg3T>;
            if (is_ready(id,obj,arg,ptr)) {
                detail::info<memfunT> info;
                arg1T arg1;
                arg2T arg2;
                arg3T arg3;
                arg & info & arg1 & arg2 & arg3;
                typename detail::info<memfunT>::futureT result(info.ref);
                arg.get_world()->taskq.add(new TaskMemfun<memfunT>(result,*obj,info.memfun,arg1,arg2,arg3,info.attr));
            }
        }

        // Handler for incoming task with 4 arguments
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        static void handler_task(const AmArg& arg) {
            const uniqueidT& id = detail::peek(arg);
            Derived* obj = arg.get_world()-> template ptr_from_id<Derived>(id);
            am_handlerT ptr = handler_task<memfunT,arg1T,arg2T,arg3T,arg4T>;
            if (is_ready(id,obj,arg,ptr)) {
                detail::info<memfunT> info;
                arg1T arg1;
                arg2T arg2;
                arg3T arg3;
                arg4T arg4;
                arg & info & arg1 & arg2 & arg3 & arg4;
                typename detail::info<memfunT>::futureT result(info.ref);
                arg.get_world()->taskq.add(new TaskMemfun<memfunT>(result,*obj,info.memfun,arg1,arg2,arg3,arg4,info.attr));
            }
        }

        // Handler for incoming task with 5 arguments
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T,
        typename arg5T>
        static void handler_task(const AmArg& arg) {
            const uniqueidT& id = detail::peek(arg);
            Derived* obj = arg.get_world()-> template ptr_from_id<Derived>(id);
            am_handlerT ptr = handler_task<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T>;
            if (is_ready(id,obj,arg,ptr)) {
                detail::info<memfunT> info;
                arg1T arg1;
                arg2T arg2;
                arg3T arg3;
                arg4T arg4;
                arg5T arg5;
                arg & info & arg1 & arg2 & arg3 & arg4 & arg5;
                typename detail::info<memfunT>::futureT result(info.ref);
                arg.get_world()->taskq.add(new TaskMemfun<memfunT>(result,*obj,info.memfun,arg1,arg2,arg3,arg4,arg5,info.attr));
            }
        }

        // Handler for incoming task with 6 arguments
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T,
        typename arg5T, typename arg6T>
        static void handler_task(const AmArg& arg) {
            const uniqueidT& id = detail::peek(arg);
            Derived* obj = arg.get_world()-> template ptr_from_id<Derived>(id);
            am_handlerT ptr = handler_task<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T,arg6T>;
            if (is_ready(id,obj,arg,ptr)) {
                detail::info<memfunT> info;
                arg1T arg1;
                arg2T arg2;
                arg3T arg3;
                arg4T arg4;
                arg5T arg5;
                arg6T arg6;
                arg & info & arg1 & arg2 & arg3 & arg4 & arg5 & arg6;
                typename detail::info<memfunT>::futureT result(info.ref);
                arg.get_world()->taskq.add(new TaskMemfun<memfunT>(result,*obj,info.memfun,arg1,arg2,arg3,arg4,arg5,arg6,info.attr));
            }
        }

        // Handler for incoming task with 6 arguments
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T,
        typename arg5T, typename arg6T, typename arg7T>
        static void handler_task(const AmArg& arg) {
            const uniqueidT& id = detail::peek(arg);
            Derived* obj = arg.get_world()-> template ptr_from_id<Derived>(id);
            am_handlerT ptr = handler_task<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T>;
            if (is_ready(id,obj,arg,ptr)) {
                detail::info<memfunT> info;
                arg1T arg1;
                arg2T arg2;
                arg3T arg3;
                arg4T arg4;
                arg5T arg5;
                arg6T arg6;
                arg7T arg7;
                arg & info & arg1 & arg2 & arg3 & arg4 & arg5 & arg6 & arg7;
                typename detail::info<memfunT>::futureT result(info.ref);
                arg.get_world()->taskq.add(new TaskMemfun<memfunT>(result,*obj,info.memfun,arg1,arg2,arg3,arg4,arg5,arg6,arg7,info.attr));
            }
        }

        // This slightly convoluted logic is to ensure ordering when
        // processing pending messages.  If a new message arrives
        // while processing incoming messages it must be queued.
        //
        // If the object does not exist ---> not ready
        // If the object exists and is ready ---> ready
        // If the object exists and is not ready then
        //    if we are doing a queued/pending message --> ready
        //    else this is a new message --> not ready
        static bool is_ready(const uniqueidT& id, Derived* obj, const AmArg& arg, am_handlerT ptr) {
            WorldObject* p=0;
            if (obj) {
                p = static_cast<WorldObject*>(obj);
                if (p->ready || arg.is_pending()) return true;
            }

            pending_mutex.lock(); // BEGIN CRITICAL SECTION
            if (obj && p->ready) {
                pending_mutex.unlock();  // ... LEAVE CRITICAL SECTION
                return true;
            }
            arg.set_pending();
            MADNESS_ASSERT(arg.is_pending());
            const_cast<pendingT&>(pending).push_back(detail::PendingMsg(id, ptr, arg));
            MADNESS_ASSERT(const_cast<pendingT&>(pending).back().arg->is_pending());
            pending_mutex.unlock(); // END CRITICAL SECTION

            return false;
        };

    protected:

        /// To be called from \em derived constructor to process pending messages

        /// Cannot call this from the WorldObject constructor since the
        /// derived class would not yet be fully constructed.
        ///
        /// !! No incoming messages are processed until this routine is invoked
        /// so the Derived class may rely upon well defined state until this routine
        /// is invoked.
        void process_pending() {
            // Messages may be arriving while we are processing the
            // pending queue.  To maximize concurrency copy messages
            // out of queue before processing outside critical section.
            while (!ready) {
                pendingT tmp;

                pending_mutex.lock(); // BEGIN CRITICAL SECTION
                pendingT& nv = const_cast<pendingT&>(pending);
                for (pendingT::iterator it=nv.begin(); it!=nv.end();) {
                    detail::PendingMsg& p = *it;
                    if (p.id == objid) {
                        tmp.push_back(p);
                        it = nv.erase(it);
                    }
                    else {
                        ++it;
                    }
                }
                if (tmp.size() == 0) ready=true;
                pending_mutex.unlock(); // END CRITICAL SECTION

                while (tmp.size()) {
                    tmp.front().invokehandler();
                    tmp.pop_front();
                }
            }
        };


    public:
        /// Associates object with globally unique ID

        /// !! The derived class MUST call process_pending both to
        /// process any messages that arrived prior to construction
        /// and to enable processing of future messages.
        WorldObject(World& world)
                : world(world)
                , ready(false)
                , me(world.rank())
                , objid(world.register_ptr(static_cast<Derived*>(this))) {};


        /// Returns the globally unique object ID
        const uniqueidT& id() const {
            return objid;
        }

        /// Returns a reference to the world
        World& get_world() const {
            return const_cast<WorldObject<Derived>*>(this)->world;
        }


        /// Sends active message to derived class method "returnT (this->*memfun)()"
        template <typename memfunT>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(ProcessID dest, memfunT memfun) {
            Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
            Derived* obj = static_cast<Derived*>(this);
            if (dest == me) {
                result.set((obj->*memfun)());
            }
            else {
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world));
                world.am.send(dest,handler<memfunT>,
                              new_am_arg(info));
            }
            return result;
        }

        /// Sends active message to derived class method "returnT (this->*memfun)(arg1)"
        template <typename memfunT, typename arg1T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(ProcessID dest, memfunT memfun, const arg1T& arg1) {
            Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
            Derived* obj = static_cast<Derived*>(this);
            if (dest == me) {
                result.set((obj->*memfun)(arg1));
            }
            else {
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world));
                world.am.send(dest,handler<memfunT,arg1T>,
                              new_am_arg(info,arg1));
            }
            return result;
        }


        /// Sends active message to derived class method "returnT (this->*memfun)(arg1,arg2)"
        template <typename memfunT, typename arg1T, typename arg2T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2) {
            Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
            Derived* obj = static_cast<Derived*>(this);
            if (dest == me) {
                result.set((obj->*memfun)(arg1,arg2));
            }
            else {
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world));
                world.am.send(dest, handler<memfunT,arg1T,arg2T>,
                              new_am_arg(info,arg1,arg2));
            }
            return result;
        }


        /// Sends active message to derived class method "returnT (this->*memfun)(arg1,arg2,arg3)"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3) {
            Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
            Derived* obj = static_cast<Derived*>(this);
            if (dest == me) {
                result.set((obj->*memfun)(arg1,arg2,arg3));
            }
            else {
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world));
                world.am.send(dest, handler<memfunT,arg1T,arg2T,arg3T>,
                              new_am_arg(info,arg1,arg2,arg3));
            }
            return result;
        }


        /// Sends active message to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4)"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4) {
            Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
            Derived* obj = static_cast<Derived*>(this);
            if (dest == me) {
                result.set((obj->*memfun)(arg1,arg2,arg3,arg4));
            }
            else {
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world));
                world.am.send(dest, handler<memfunT,arg1T,arg2T,arg3T,arg4T>,
                              new_am_arg(info,arg1,arg2,arg3,arg4));
            }
            return result;
        }

        /// Sends active message to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4,arg5)"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4,
             const arg5T& arg5) {
            Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
            Derived* obj = static_cast<Derived*>(this);
            if (dest == me) {
                result.set((obj->*memfun)(arg1,arg2,arg3,arg4,arg5));
            }
            else {
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world));
                world.am.send(dest, handler<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T>,
                              new_am_arg(info,arg1,arg2,arg3,arg4,arg5));
            }
            return result;
        }


        /// Sends active message to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4,arg5,arg6)"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T,
        typename arg5T, typename arg6T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4,
             const arg5T& arg5, const arg6T& arg6) {
            Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
            Derived* obj = static_cast<Derived*>(this);
            if (dest == me) {
                result.set((obj->*memfun)(arg1,arg2,arg3,arg4,arg5,arg6));
            }
            else {
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world));
                world.am.send(dest, handler<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T,arg6T>,
                              new_am_arg(info,arg1,arg2,arg3,arg4,arg5,arg6));
            }
            return result;
        }


        /// Sends active message to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4,arg5,arg6,arg7)"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T,
        typename arg5T, typename arg6T, typename arg7T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4,
             const arg5T& arg5, const arg6T& arg6, const arg7T& arg7) {
            Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
            Derived* obj = static_cast<Derived*>(this);
            if (dest == me) {
                result.set((obj->*memfun)(arg1,arg2,arg3,arg4,arg5,arg6,arg7));
            }
            else {
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world));
                world.am.send(dest, handler<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T>,
                              new_am_arg(info,arg1,arg2,arg3,arg4,arg5,arg6,arg7));
            }
            return result;
        }


        /// Sends active message to derived class method "returnT (this->*memfun)() const"
        template <typename memfunT, typename arg1T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(ProcessID dest, memfunT memfun) const {
            return const_cast<objT*>(this)->send(dest,memfun);
        }

        /// Sends active message to derived class method "returnT (this->*memfun)(arg1) const"
        template <typename memfunT, typename arg1T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(ProcessID dest, memfunT memfun, const arg1T& arg1) const {
            return const_cast<objT*>(this)->send(dest,memfun,arg1);
        }

        /// Sends active message to derived class method "returnT (this->*memfun)(arg1,arg2) const"
        template <typename memfunT, typename arg1T, typename arg2T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2) const {
            return const_cast<objT*>(this)->send(dest,memfun,arg1,arg2);
        }

        /// Sends active message to derived class method "returnT (this->*memfun)(arg1,arg2,arg3) const"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3) const {
            return const_cast<objT*>(this)->send(dest,memfun,arg1,arg2,arg3);
        }

        /// Sends active message to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4) const"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4) const {
            return const_cast<objT*>(this)->send(dest,memfun,arg1,arg2,arg3,arg4);
        }

        /// Sends active message to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4,arg5) const"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4,
             const arg5T& arg5) const {
            return const_cast<objT*>(this)->send(dest,memfun,arg1,arg2,arg3,arg4,arg5);
        }


        /// Sends active message to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4,arg5, arg6) const"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4,
             const arg5T& arg5, const arg6T& arg6) const {
            return const_cast<objT*>(this)->send(dest,memfun,arg1,arg2,arg3,arg4,arg5,arg6);
        }


        /// Sends active message to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4,arg5, arg6, arg7) const"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T,
        typename arg5T, typename arg6T, typename arg7T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        send(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4,
             const arg5T& arg5, const arg6T& arg6, const arg7T& arg7) const {
            return const_cast<objT*>(this)->send(dest,memfun,arg1,arg2,arg3,arg4,arg5,arg6,arg7);
        }


        /// Sends task to derived class method "returnT (this->*memfun)()"
        template <typename memfunT>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(ProcessID dest, memfunT memfun, const TaskAttributes& attr = TaskAttributes()) {
            if (dest == me) {
                return world.taskq.add(*static_cast<Derived*>(this), memfun, attr);
            }
            else {
                Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world), attr);
                world.am.send(dest,handler_task<memfunT>, new_am_arg(info));
                return result;
            }
        }

        /// Sends task to derived class method "returnT (this->*memfun)(arg1)"
        template <typename memfunT, typename arg1T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const TaskAttributes& attr = TaskAttributes()) {
            if (dest == me) {
                return world.taskq.add(*static_cast<Derived*>(this), memfun, arg1, attr);
            }
            else {
                Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world), attr);
                world.am.send(dest,handler_task<memfunT,arg1T>,
                              new_am_arg(info,arg1));
                return result;
            }
        }

        /// Sends task to derived class method "returnT (this->*memfun)(arg1,arg2)"
        template <typename memfunT, typename arg1T, typename arg2T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2,
             const TaskAttributes& attr = TaskAttributes()) {
            if (dest == me) {
                return world.taskq.add(*static_cast<Derived*>(this), memfun, arg1, arg2, attr);
            }
            else {
                Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world), attr);
                world.am.send(dest,handler_task<memfunT,arg1T,arg2T>,
                              new_am_arg(info,arg1, arg2));
                return result;
            }
        }

        /// Sends task to derived class method "returnT (this->*memfun)(arg1,arg2,arg3)"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3,
             const TaskAttributes& attr = TaskAttributes()) {
            if (dest == me) {
                return world.taskq.add(*static_cast<Derived*>(this), memfun, arg1, arg2, arg3, attr);
            }
            else {
                Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world), attr);
                world.am.send(dest,handler_task<memfunT,arg1T,arg2T,arg3T>,
                              new_am_arg(info,arg1, arg2, arg3));
                return result;
            }
        }

        /// Sends task to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4)"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3,
             const arg4T& arg4, const TaskAttributes& attr = TaskAttributes()) {
            if (dest == me) {
                return world.taskq.add(*static_cast<Derived*>(this), memfun, arg1, arg2, arg3, arg4, attr);
            }
            else {
                Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world), attr);
                world.am.send(dest,handler_task<memfunT,arg1T,arg2T,arg3T,arg4T>,
                              new_am_arg(info,arg1, arg2, arg3, arg4));
                return result;
            }
        }

        /// Sends task to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4,arg5)"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3,
             const arg4T& arg4, const arg5T& arg5, const TaskAttributes& attr = TaskAttributes()) {
            if (dest == me) {
                return world.taskq.add(*static_cast<Derived*>(this), memfun, arg1, arg2, arg3, arg4, arg5, attr);
            }
            else {
                Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world), attr);
                world.am.send(dest,handler_task<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T>,
                              new_am_arg(info,arg1, arg2, arg3, arg4, arg5));
                return result;
            }
        }

        /// Sends task to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4,arg5,arg6)"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3,
             const arg4T& arg4, const arg5T& arg5, const arg6T& arg6, const TaskAttributes& attr = TaskAttributes()) {
            if (dest == me) {
                return world.taskq.add(*static_cast<Derived*>(this), memfun, arg1, arg2, arg3, arg4, arg5, arg6, attr);
            }
            else {
                Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world), attr);
                world.am.send(dest,handler_task<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T,arg6T>,
                              new_am_arg(info,arg1, arg2, arg3, arg4, arg5, arg6));
                return result;
            }
        }

        /// Sends task to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4,arg5,arg6,arg7)"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3,
             const arg4T& arg4, const arg5T& arg5, const arg6T& arg6, const arg7T& arg7, const TaskAttributes& attr = TaskAttributes()) {
            if (dest == me) {
                return world.taskq.add(*static_cast<Derived*>(this), memfun, arg1, arg2, arg3, arg4, arg5, arg6, arg7, attr);
            }
            else {
                Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world), attr);
                world.am.send(dest,handler_task<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T>,
                              new_am_arg(info,arg1, arg2, arg3, arg4, arg5, arg6, arg7));
                return result;
            }
        }


        /// Sends task to derived class method "returnT (this->*memfun)() const"
        template <typename memfunT>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(ProcessID dest, memfunT memfun, const TaskAttributes& attr = TaskAttributes()) const {
            return const_cast<objT*>(this)->task(dest,memfun,attr);
        }

        /// Sends task to derived class method "returnT (this->*memfun)(arg1) const"
        template <typename memfunT, typename arg1T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const TaskAttributes& attr = TaskAttributes()) const {
            return const_cast<objT*>(this)->task(dest,memfun,arg1,attr);
        }

        /// Sends task to derived class method "returnT (this->*memfun)(arg1,arg2) const"
        template <typename memfunT, typename arg1T, typename arg2T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2,
             const TaskAttributes& attr = TaskAttributes()) const {
            return const_cast<objT*>(this)->task(dest,memfun,arg1,arg2,attr);
        }

        /// Sends task to derived class method "returnT (this->*memfun)(arg1,arg2,arg3) const"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3,
             const TaskAttributes& attr = TaskAttributes()) const {
            return const_cast<objT*>(this)->task(dest,memfun,arg1,arg2,arg3,attr);
        }

        /// Sends task to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4) const"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3,
             const arg4T& arg4, const TaskAttributes& attr = TaskAttributes()) const {
            return const_cast<objT*>(this)->task(dest,memfun,arg1,arg2,arg3,arg4,attr);
        }

        /// Sends task to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4,arg5) const"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3,
             const arg4T& arg4, const arg5T& arg5, const TaskAttributes& attr = TaskAttributes()) const {
            return const_cast<objT*>(this)->task(dest,memfun,arg1,arg2,arg3,arg4,arg5,attr);
        }

        /// Sends task to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4,arg5,arg6) const"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) >
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3,
             const arg4T& arg4, const arg5T& arg5, const arg6T& arg6, const TaskAttributes& attr = TaskAttributes()) const {
            return const_cast<objT*>(this)->task(dest,memfun,arg1,arg2,arg3,arg4,arg5,arg6,attr);
        }

        virtual ~WorldObject() {
            world.unregister_ptr(static_cast<Derived*>(this));
        }
    };

    namespace archive {
        template <class Derived>
        struct ArchiveLoadImpl<BufferInputArchive,WorldObject<Derived>*> {
            static inline void load(const BufferInputArchive& ar, WorldObject<Derived>*& ptr) {
                uniqueidT id;
                ar & id;
                World* world = World::world_from_id(id.get_world_id());
                MADNESS_ASSERT(world);
                ptr = world->ptr_from_id< WorldObject<Derived> >(id);
                if (!ptr) MADNESS_EXCEPTION("WorldObj: remote operation attempting to use a locally uninitialized object",0);
            }
        };

        template <class Derived>
        struct ArchiveStoreImpl<BufferOutputArchive,WorldObject<Derived>*> {
            static inline void store(const BufferOutputArchive& ar, const WorldObject<Derived>*const& ptr) {
                ar & ptr->id();
            }
        };

        template <class Derived>
        struct ArchiveLoadImpl<BufferInputArchive,const WorldObject<Derived>*> {
            static inline void load(const BufferInputArchive& ar, const WorldObject<Derived>*& ptr) {
                uniqueidT id;
                ar & id;
                World* world = World::world_from_id(id.get_world_id());
                MADNESS_ASSERT(world);
                ptr = world->ptr_from_id< WorldObject<Derived> >(id);
                if (!ptr) MADNESS_EXCEPTION("WorldObj: remote operation attempting to use a locally uninitialized object",0);
            }
        };

        template <class Derived>
        struct ArchiveStoreImpl<BufferOutputArchive,const WorldObject<Derived>*> {
            static inline void store(const BufferOutputArchive& ar, const WorldObject<Derived>*const& ptr) {
                ar & ptr->id();
            }
        };
    }
}

#ifdef WORLD_INSTANTIATE_STATIC_TEMPLATES
template <typename Derived>
volatile std::list<madness::detail::PendingMsg> madness::WorldObject<Derived>::pending;

template <typename Derived>  madness::Spinlock madness::WorldObject<Derived>::pending_mutex;

#endif

#endif
