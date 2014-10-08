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

/// \file world/worldobj.h
/// \brief Defines and implements WorldObject
/// \ingroup worldobj

#ifndef MADNESS_WORLD_WORLDOBJ_H__INCLUDED
#define MADNESS_WORLD_WORLDOBJ_H__INCLUDED

#include <madness/world/worldthread.h>
#include <madness/world/worldtask.h>

namespace madness {

    template <typename> class WorldObject;

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
            }
        };

        // It is annoying that we must replicate the task forwarding stuff here but we must
        // so that pending messages creating tasks are correctly handled.  Cannot merge
        // easily with send handlers since task layer is more restrictive on the
        // copy capability of arguments.

        // It is also annoying that info needs to be broken into two parts so
        // that it id and ref are properly serialized. We need to have id
        // correctly aligned via opaque_wrap, but ref cannot be serialized that
        // way. Thus we break the class into two parts.

        // Info stored for AM method forwarding
        template <typename memfunT>
        struct info_base {
            uniqueidT id; // Must be at front ... see peek.
            ProcessID requestor;
            memfunT memfun;
            TaskAttributes attr;

        protected:

            info_base() {}

            info_base(const uniqueidT& id, ProcessID requestor,  memfunT memfun,
                 const TaskAttributes& attr=TaskAttributes())
                    : id(id)
                    , requestor(requestor)
                    , memfun(memfun)
                    , attr(attr) {}


            template <typename Archive>
            void serialize(const Archive& ar) {
                ar & archive::wrap_opaque(*this); // Must be opaque ... see peek.
            }
        }; // struct info_base

        template <typename memfunT>
        struct info : public info_base<memfunT> {
            typedef Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > futureT;
            typedef RemoteReference< FutureImpl< REMFUTURE(MEMFUN_RETURNT(memfunT)) > > refT;
            refT ref;

            info() : info_base<memfunT>() {}

            info(const AmArg& arg) :
                info_base<memfunT>()
            {
                arg & *this;
            }

            info(const uniqueidT& id, ProcessID requestor,  memfunT memfun, const refT& ref,
                 const TaskAttributes& attr=TaskAttributes()) :
                     info_base<memfunT>(id, requestor, memfun, attr), ref(ref)
            {}

            template <typename Archive>
            void serialize(const Archive& ar) {
                info_base<memfunT>::serialize(ar);
                ar & ref;
            }
        }; // struct info

        // Extract the unique object ID from an incoming active message header

        // We deserialize the header and all arguments at the same
        // time to simplify the code.  However, it is common that
        // when sending a message to an item in a container to
        // include a pointer to the container itself.  But this
        // breaks if the container is not initialized since the
        // deserialization throws if the object is not initialized
        // (which seems preferable to hidden race condition).  Hence,
        // we use this routine to extract the unique ID from the very
        // front of the info structure.  For efficiency we here rely
        // upon the serialization of info being opaque and the
        // id being at the front of info.
        static inline const uniqueidT& peek(const AmArg& arg) {
            return *((uniqueidT*)(arg.buf()));
        }



        template <typename objT, typename memfnT, typename Enabler = void>
        struct WorldObjectTaskHelper {
            typedef typename if_c<memfunc_traits<memfnT>::constness,
                    const objT*, objT*>::type ptrT;
            typedef MemFuncWrapper<ptrT, memfnT, typename result_of<memfnT>::type> wrapperT;


            static wrapperT make_task_fn(const objT* obj, memfnT memfn) {
                return wrapperT(const_cast<ptrT>(obj), memfn);
            }

            static wrapperT make_task_fn(objT* obj, memfnT memfn) {
                return wrapperT(const_cast<ptrT>(obj), memfn);
            }

            static wrapperT make_task_fn(const WorldObject<objT>* obj, memfnT memfn) {
                return make_task_fn(static_cast<const objT*>(obj), memfn);
            }

            static wrapperT make_task_fn(WorldObject<objT>* obj, memfnT memfn) {
                return make_task_fn(static_cast<objT*>(obj), memfn);
            }
        }; // struct WorldObjectTaskHelpe


#ifndef MADNESS_DISABLE_SHARED_FROM_THIS
        // Disable the use of std::enable_shared_from_this if we are using MADNESS's
        // implementation since weak_ptr is not fully implemented.

        template <typename objT, typename memfnT>
        struct WorldObjectTaskHelper<objT, memfnT,
                typename std::is_base_of<std::enable_shared_from_this<objT>, objT>::type>
        {
            typedef typename if_c<memfunc_traits<memfnT>::constness,
                    const objT*, objT*>::type ptrT;
            typedef typename if_c<memfunc_traits<memfnT>::constness,
                    std::shared_ptr<const objT>, std::shared_ptr<objT> >::type shared_ptrT;
            typedef MemFuncWrapper<shared_ptrT, memfnT, typename result_of<memfnT>::type> wrapperT;

            static wrapperT make_task_fn(const objT* obj, memfnT memfn) {
                return wrapperT(const_cast<ptrT>(obj)->shared_from_this(), memfn);
            }

            static wrapperT make_task_fn(objT* obj, memfnT memfn) {
                return wrapperT(const_cast<ptrT>(obj)->shared_from_this(), memfn);
            }

            static wrapperT make_task_fn(const WorldObject<objT>* obj, memfnT memfn) {
                return make_task_fn(static_cast<const objT*>(obj), memfn);
            }

            static wrapperT make_task_fn(WorldObject<objT>* obj, memfnT memfn) {
                return make_task_fn(static_cast<objT*>(obj), memfn);
            }
        }; // struct WorldObjectTaskHelpe

#endif // MADNESS_DISABLE_SHARED_FROM_THIS

    } // namespace detail


    /// Implements most parts of a globally addressable object (via unique ID)

    /// \ingroup worldobj
    /// 1) Derived class has WorldObject<Derived> as a public base class
    /// 2) Derived constructor
    ///    a) invokes WorldObject<Derived>(world) constructor
    ///    b) invokes process_pending()
    /// 3) Derived destructor must either be deferred or preceeded by gop.fence()
    /// 4) Derived class must have at least one virtual function for serialization
    ///    of derived class pointers to be cast to the appropriate type.
    ///
    /// This class is deliberately not default constructible and does
    /// not support assignment or copying.  This ensures that each instance
    /// is unique.  Have a look at the WorldContainer for an example
    /// of wrapping this using the PIMPL idiom and a shared pointer.
    ///
    /// Note that world is exposed for convenience as a public data member.
    template <class Derived>
    class WorldObject {
    public:
        // Typedefs
        typedef WorldObject<Derived> objT;

    private:
        typedef std::list<detail::PendingMsg> pendingT;
        typedef detail::voidT voidT;

        World& world;                              ///< Think globally act locally
        // The order here matters in a multi-threaded world
        volatile bool ready;                       ///< True if ready to rock 'n roll
        ProcessID me;                              ///< Rank of self
        uniqueidT objid;                           ///< Sense of self


        static Spinlock pending_mutex;
        static volatile pendingT pending;          ///< Buffers pending messages

        // This slightly convoluted logic is to ensure ordering when
        // processing pending messages.  If a new message arrives
        // while processing incoming messages it must be queued.
        //
        // If the object does not exist ---> not ready
        // If the object exists and is ready ---> ready
        // If the object exists and is not ready then
        //    if we are doing a queued/pending message --> ready
        //    else this is a new message --> not ready
        static bool is_ready(const uniqueidT& id, objT*& obj, const AmArg& arg, am_handlerT ptr) {
            obj = static_cast<objT*>(arg.get_world()->template ptr_from_id<Derived>(id));

            if (obj) {
                if (obj->ready || arg.is_pending()) return true;
            }

            ScopedMutex<Spinlock> lock(pending_mutex); // BEGIN CRITICAL SECTION
            if (!obj) obj = static_cast<objT*>(arg.get_world()->template ptr_from_id<Derived>(id));

            if (obj) {
                if (obj->ready || arg.is_pending())
                    return true; // END CRITICAL SECTION
            }
            const_cast<AmArg&>(arg).set_pending();
            const_cast<pendingT&>(pending).push_back(detail::PendingMsg(id, ptr, arg));

            return false; // END CRITICAL SECTION
        }


        // Handler for incoming AM
        template <typename memfnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T,
        typename arg5T, typename arg6T, typename arg7T, typename arg8T, typename arg9T>
        static void handler(const AmArg& arg) {
            typedef detail::WorldObjectTaskHelper<Derived, memfnT> task_helper;

            const uniqueidT& id = detail::peek(arg);
            am_handlerT ptr = handler<memfnT,arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T,arg9T>;
            objT* obj;
            if (is_ready(id,obj,arg,ptr)) {
                detail::info<memfnT> info;
                typename detail::task_arg<arg1T>::type arg1;
                typename detail::task_arg<arg2T>::type arg2;
                typename detail::task_arg<arg3T>::type arg3;
                typename detail::task_arg<arg4T>::type arg4;
                typename detail::task_arg<arg5T>::type arg5;
                typename detail::task_arg<arg6T>::type arg6;
                typename detail::task_arg<arg7T>::type arg7;
                typename detail::task_arg<arg8T>::type arg8;
                typename detail::task_arg<arg9T>::type arg9;
                arg & info & arg1 & arg2 & arg3 & arg4 & arg5 & arg6 & arg7 & arg8 & arg9;
                typename detail::info<memfnT>::futureT result(info.ref);
                detail::run_function(result, task_helper::make_task_fn(obj, info.memfun),
                        arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
            }
        }

        // Handler for remote arguments
        template <typename taskT>
        static void spawn_remote_task_handler(const AmArg& arg) {
            typedef detail::WorldObjectTaskHelper<Derived,
                    typename taskT::functionT::memfn_type> task_helper;

            MADNESS_ASSERT(taskT::arity <= 9u);

            const uniqueidT& id = detail::peek(arg);
            am_handlerT ptr = & objT::template spawn_remote_task_handler<taskT>;
            objT* obj;
            if (is_ready(id,obj,arg,ptr)) {
                detail::info<typename taskT::functionT::memfn_type> info;

                archive::BufferInputArchive input_arch = arg & info;

                // Construct task
                taskT* task = new taskT(typename taskT::futureT(info.ref),
                        task_helper::make_task_fn(obj, info.memfun), info.attr, input_arch);

                // Add task to queue
                arg.get_world()->taskq.add(task);
            }
        }

        template <typename T>
        static inline const T& am_arg(const Future<T>& f) {
            MADNESS_ASSERT(f.probe()); // Cannot serialize unassigned futures
            return f.get();
        }

        template <typename T> static inline const T& am_arg(const T& t) { return t; }

        template <typename memfnT, typename a1T, typename a2T, typename a3T,
                typename a4T, typename a5T, typename a6T, typename a7T,
                typename a8T, typename a9T>
        typename detail::task_result_type<memfnT>::futureT
        send_am(ProcessID dest, memfnT memfn, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6,
                const a7T& a7, const a8T& a8, const a9T& a9) const
        {
            typedef detail::WorldObjectTaskHelper<Derived, memfnT> task_helper;
            typename detail::task_result_type<memfnT>::futureT result;
            if (dest == me)
                detail::run_function(result, task_helper::make_task_fn(this, memfn),
                        a1, a2, a3, a4, a5, a6, a7, a8, a9);
            else {
                detail::info<memfnT> info(objid, me, memfn, result.remote_ref(world));
                world.am.send(dest, & objT::template handler<memfnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T>,
                        new_am_arg(info, a1, a2, a3, a4, a5, a6, a7, a8, a9));
            }

            return result;
        }

        template <typename taskT, typename memfnT, typename a1T, typename a2T, typename a3T,
                typename a4T, typename a5T, typename a6T, typename a7T,
                typename a8T, typename a9T>
        typename taskT::futureT
        send_task(ProcessID dest, memfnT memfn, const a1T& a1,
                const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5,
                const a6T& a6, const a7T& a7, const a8T& a8, const a9T& a9,
                const TaskAttributes& attr) const
        {
            typename taskT::futureT result;
            detail::info<memfnT> info(objid, me, memfn, result.remote_ref(world), attr);
            world.am.send(dest, & objT::template spawn_remote_task_handler<taskT>,
                    new_am_arg(info, a1, a2, a3, a4, a5, a6, a7, a8, a9));

            return result;
        }

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
            //int ndone = 0;
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
                    //++ndone;
                }
            }
            //if (ndone) std::cout << world.rank() << ":pending:" << ndone << std::endl;
        }


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

        template <typename memfnT>
        typename detail::task_result_type<memfnT>::futureT
        send(ProcessID dest, memfnT memfn) const {
            return send_am(dest, memfn, voidT::value, voidT::value, voidT::value,
                    voidT::value, voidT::value, voidT::value, voidT::value,
                    voidT::value, voidT::value);
        }

        template <typename memfnT, typename a1T>
        typename detail::task_result_type<memfnT>::futureT
        send(ProcessID dest, memfnT memfn, const a1T& a1) const {
            return send_am(dest, memfn, am_arg(a1), voidT::value, voidT::value,
                    voidT::value, voidT::value, voidT::value, voidT::value,
                    voidT::value, voidT::value);
        }

        template <typename memfnT, typename a1T, typename a2T>
        typename detail::task_result_type<memfnT>::futureT
        send(ProcessID dest, memfnT memfn, const a1T& a1, const a2T& a2) const {
            return send_am(dest, memfn, am_arg(a1), am_arg(a2), voidT::value,
                    voidT::value, voidT::value, voidT::value, voidT::value,
                    voidT::value, voidT::value);
        }

        template <typename memfnT, typename a1T, typename a2T, typename a3T>
        typename detail::task_result_type<memfnT>::futureT
        send(ProcessID dest, memfnT memfn, const a1T& a1, const a2T& a2,
                const a3T& a3) const
        {
            return send_am(dest, memfn, am_arg(a1), am_arg(a2), am_arg(a3),
                    voidT::value, voidT::value, voidT::value, voidT::value,
                    voidT::value, voidT::value);
        }

        template <typename memfnT, typename a1T, typename a2T, typename a3T,
                typename a4T>
        typename detail::task_result_type<memfnT>::futureT
        send(ProcessID dest, memfnT memfn, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4) const
        {
            return send_am(dest, memfn, am_arg(a1), am_arg(a2), am_arg(a3),
                    am_arg(a4), voidT::value, voidT::value, voidT::value,
                    voidT::value, voidT::value);
        }

        template <typename memfnT, typename a1T, typename a2T, typename a3T,
                typename a4T, typename a5T>
        typename detail::task_result_type<memfnT>::futureT
        send(ProcessID dest, memfnT memfn, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5) const
        {
            return send_am(dest, memfn, am_arg(a1), am_arg(a2), am_arg(a3),
                    am_arg(a4), am_arg(a5), voidT::value, voidT::value,
                    voidT::value, voidT::value);
        }

        template <typename memfnT, typename a1T, typename a2T, typename a3T,
                typename a4T, typename a5T, typename a6T>
        typename detail::task_result_type<memfnT>::futureT
        send(ProcessID dest, memfnT memfn, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6) const
        {
            return send_am(dest, memfn, am_arg(a1), am_arg(a2), am_arg(a3),
                    am_arg(a4), am_arg(a5), am_arg(a6), voidT::value, voidT::value,
                    voidT::value);
        }

        template <typename memfnT, typename a1T, typename a2T, typename a3T,
                typename a4T, typename a5T, typename a6T, typename a7T>
        typename detail::task_result_type<memfnT>::futureT
        send(ProcessID dest, memfnT memfn, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6,
                const a7T& a7) const
        {
            return send_am(dest, memfn, am_arg(a1), am_arg(a2), am_arg(a3),
                    am_arg(a4), am_arg(a5), am_arg(a6), am_arg(a7), voidT::value,
                    voidT::value);
        }

        template <typename memfnT, typename a1T, typename a2T, typename a3T,
                typename a4T, typename a5T, typename a6T, typename a7T, typename a8T>
        typename detail::task_result_type<memfnT>::futureT
        send(ProcessID dest, memfnT memfn, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6,
                const a7T& a7, const a8T& a8) const
        {
            return send_am(dest, memfn, am_arg(a1), am_arg(a2), am_arg(a3),
                    am_arg(a4), am_arg(a5), am_arg(a6), am_arg(a7), am_arg(a8),
                    voidT::value);
        }

        template <typename memfnT, typename a1T, typename a2T, typename a3T,
                typename a4T, typename a5T, typename a6T, typename a7T, typename a8T,
                typename a9T>
        typename detail::task_result_type<memfnT>::futureT
        send(ProcessID dest, memfnT memfn, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6,
                const a7T& a7, const a8T& a8, const a9T& a9) const
        {
            return send_am(dest, memfn, am_arg(a1), am_arg(a2), am_arg(a3),
                    am_arg(a4), am_arg(a5), am_arg(a6), am_arg(a7), am_arg(a8),
                    am_arg(a9));
        }

        /// Sends task to derived class method "returnT (this->*memfn)(a1,a2,a3,a4,a5,a6,a7,a8,a9)"
        template <typename memfnT>
        typename detail::task_result_type<memfnT>::futureT
        task(ProcessID dest, memfnT memfn, const TaskAttributes& attr = TaskAttributes()) const {
            typedef detail::WorldObjectTaskHelper<Derived, memfnT> task_helper;
            typedef TaskFn<typename task_helper::wrapperT> taskT;
            if (dest == me)
                return world.taskq.add(task_helper::make_task_fn(this, memfn), attr);
            else
                return send_task<taskT>(dest, memfn, voidT::value, voidT::value,
                        voidT::value, voidT::value, voidT::value, voidT::value,
                        voidT::value, voidT::value, voidT::value, attr);
        }

        /// Sends task to derived class method "returnT (this->*memfn)(a1)"
        template <typename memfnT, typename a1T>
        typename detail::task_result_type<memfnT>::futureT
        task(ProcessID dest, memfnT memfn, const a1T& a1,
                const TaskAttributes& attr = TaskAttributes()) const
        {
            typedef detail::WorldObjectTaskHelper<Derived, memfnT> task_helper;
            typedef TaskFn<typename task_helper::wrapperT,
                    typename detail::task_arg<a1T>::type> taskT;
            if (dest == me)
                return world.taskq.add(task_helper::make_task_fn(this, memfn),
                        a1, attr);
            else
                return send_task<taskT>(dest, memfn, am_arg(a1), voidT::value,
                        voidT::value, voidT::value, voidT::value, voidT::value,
                        voidT::value, voidT::value, voidT::value, attr);
        }

        /// Sends task to derived class method "returnT (this->*memfn)(a1,a2)"
        template <typename memfnT, typename a1T, typename a2T>
        typename detail::task_result_type<memfnT>::futureT
        task(ProcessID dest, memfnT memfn, const a1T& a1, const a2T& a2,
                const TaskAttributes& attr = TaskAttributes()) const
        {
            typedef detail::WorldObjectTaskHelper<Derived, memfnT> task_helper;
            typedef TaskFn<typename task_helper::wrapperT,
                    typename detail::task_arg<a1T>::type,
                    typename detail::task_arg<a2T>::type> taskT;
            if (dest == me)
                return world.taskq.add(task_helper::make_task_fn(this, memfn),
                        a1, a2, attr);
            else
                return send_task<taskT>(dest, memfn, am_arg(a1), am_arg(a2),
                        voidT::value, voidT::value, voidT::value, voidT::value,
                        voidT::value, voidT::value, voidT::value, attr);
        }

        /// Sends task to derived class method "returnT (this->*memfn)(a1,a2,a3)"
        template <typename memfnT, typename a1T, typename a2T, typename a3T>
        typename detail::task_result_type<memfnT>::futureT
        task(ProcessID dest, memfnT memfn, const a1T& a1, const a2T& a2,
                const a3T& a3, const TaskAttributes& attr = TaskAttributes()) const
        {
            typedef detail::WorldObjectTaskHelper<Derived, memfnT> task_helper;
            typedef TaskFn<typename task_helper::wrapperT,
                    typename detail::task_arg<a1T>::type,
                    typename detail::task_arg<a2T>::type,
                    typename detail::task_arg<a3T>::type> taskT;
            if (dest == me)
                return world.taskq.add(task_helper::make_task_fn(this, memfn),
                        a1, a2, a3, attr);
            else
                return send_task<taskT>(dest, memfn, am_arg(a1), am_arg(a2),
                        am_arg(a3), voidT::value, voidT::value, voidT::value,
                        voidT::value, voidT::value, voidT::value, attr);
        }

        /// Sends task to derived class method "returnT (this->*memfn)(a1,a2,a3,a4)"
        template <typename memfnT, typename a1T, typename a2T, typename a3T,
                typename a4T>
        typename detail::task_result_type<memfnT>::futureT
        task(ProcessID dest, memfnT memfn, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const TaskAttributes& attr = TaskAttributes()) const
        {
            typedef detail::WorldObjectTaskHelper<Derived, memfnT> task_helper;
            typedef TaskFn<typename task_helper::wrapperT,
                    typename detail::task_arg<a1T>::type,
                    typename detail::task_arg<a2T>::type,
                    typename detail::task_arg<a3T>::type,
                    typename detail::task_arg<a4T>::type> taskT;
            if (dest == me)
                return world.taskq.add(task_helper::make_task_fn(this, memfn),
                        a1, a2, a3, a4, attr);
            else
                return send_task<taskT>(dest, memfn, am_arg(a1), am_arg(a2),
                        am_arg(a3), am_arg(a4), voidT::value, voidT::value,
                        voidT::value, voidT::value, voidT::value, attr);
        }

        /// Sends task to derived class method "returnT (this->*memfn)(a1,a2,a3,a4,a5)"
        template <typename memfnT, typename a1T, typename a2T, typename a3T,
                typename a4T, typename a5T>
        typename detail::task_result_type<memfnT>::futureT
        task(ProcessID dest, memfnT memfn, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5,
                const TaskAttributes& attr = TaskAttributes()) const
        {
            typedef detail::WorldObjectTaskHelper<Derived, memfnT> task_helper;
            typedef TaskFn<typename task_helper::wrapperT,
                    typename detail::task_arg<a1T>::type,
                    typename detail::task_arg<a2T>::type,
                    typename detail::task_arg<a3T>::type,
                    typename detail::task_arg<a4T>::type,
                    typename detail::task_arg<a5T>::type> taskT;
            if (dest == me)
                return world.taskq.add(task_helper::make_task_fn(this, memfn),
                        a1, a2, a3, a4, a5, attr);
            else
                return send_task<taskT>(dest, memfn, am_arg(a1), am_arg(a2),
                        am_arg(a3), am_arg(a4), am_arg(a5), voidT::value,
                        voidT::value, voidT::value, voidT::value, attr);
        }

        /// Sends task to derived class method "returnT (this->*memfn)(a1,a2,a3,a4,a5,a6)"
        template <typename memfnT, typename a1T, typename a2T, typename a3T,
                typename a4T, typename a5T, typename a6T>
        typename detail::task_result_type<memfnT>::futureT
        task(ProcessID dest, memfnT memfn, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6,
                const TaskAttributes& attr = TaskAttributes()) const
        {
            typedef detail::WorldObjectTaskHelper<Derived, memfnT> task_helper;
            typedef TaskFn<typename task_helper::wrapperT,
                    typename detail::task_arg<a1T>::type,
                    typename detail::task_arg<a2T>::type,
                    typename detail::task_arg<a3T>::type,
                    typename detail::task_arg<a4T>::type,
                    typename detail::task_arg<a5T>::type,
                    typename detail::task_arg<a6T>::type> taskT;
            if (dest == me)
                return world.taskq.add(task_helper::make_task_fn(this, memfn),
                        a1, a2, a3, a4, a5, a6, attr);
            else {
                return send_task<taskT>(dest, memfn, am_arg(a1), am_arg(a2),
                        am_arg(a3), am_arg(a4), am_arg(a5), am_arg(a6),
                        voidT::value, voidT::value, voidT::value, attr);
            }
        }

        /// Sends task to derived class method "returnT (this->*memfn)(a1,a2,a3,a4,a5,a6,a7)"
        template <typename memfnT, typename a1T, typename a2T, typename a3T,
                typename a4T, typename a5T, typename a6T, typename a7T>
        typename detail::task_result_type<memfnT>::futureT
        task(ProcessID dest, memfnT memfn, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6,
                const a7T& a7, const TaskAttributes& attr = TaskAttributes()) const
        {
            typedef detail::WorldObjectTaskHelper<Derived, memfnT> task_helper;
            typedef TaskFn<typename task_helper::wrapperT,
                    typename detail::task_arg<a1T>::type,
                    typename detail::task_arg<a2T>::type,
                    typename detail::task_arg<a3T>::type,
                    typename detail::task_arg<a4T>::type,
                    typename detail::task_arg<a5T>::type,
                    typename detail::task_arg<a6T>::type,
                    typename detail::task_arg<a7T>::type> taskT;
            if (dest == me)
                return world.taskq.add(task_helper::make_task_fn(this, memfn),
                        a1, a2, a3, a4, a5, a6, a7, attr);
            else
                return send_task<taskT>(dest, memfn, am_arg(a1), am_arg(a2),
                        am_arg(a3), am_arg(a4), am_arg(a5), am_arg(a6),
                        am_arg(a7), voidT::value, voidT::value, attr);
        }

        /// Sends task to derived class method "returnT (this->*memfn)(a1,a2,a3,a4,a5,a6,a7,a8)"
        template <typename memfnT, typename a1T, typename a2T, typename a3T,
                typename a4T, typename a5T, typename a6T, typename a7T, typename a8T>
        typename detail::task_result_type<memfnT>::futureT
        task(ProcessID dest, memfnT memfn, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6,
                const a7T& a7, const a8T& a8,
                const TaskAttributes& attr = TaskAttributes()) const
        {
            typedef detail::WorldObjectTaskHelper<Derived, memfnT> task_helper;
            typedef TaskFn<typename task_helper::wrapperT,
                    typename detail::task_arg<a1T>::type,
                    typename detail::task_arg<a2T>::type,
                    typename detail::task_arg<a3T>::type,
                    typename detail::task_arg<a4T>::type,
                    typename detail::task_arg<a5T>::type,
                    typename detail::task_arg<a6T>::type,
                    typename detail::task_arg<a7T>::type,
                    typename detail::task_arg<a8T>::type> taskT;
            if (dest == me)
                return world.taskq.add(task_helper::make_task_fn(this, memfn),
                        a1, a2, a3, a4, a5, a6, a7, a8, attr);
            else {
                return send_task<taskT>(dest, memfn, am_arg(a1), am_arg(a2),
                        am_arg(a3), am_arg(a4), am_arg(a5), am_arg(a6),
                        am_arg(a7), am_arg(a8), voidT::value, attr);
            }
        }

        /// Sends task to derived class method "returnT (this->*memfn)(a1,a2,a3,a4,a5,a6,a7,a8,a9)"
        template <typename memfnT, typename a1T, typename a2T, typename a3T,
                typename a4T, typename a5T, typename a6T, typename a7T, typename a8T,
                typename a9T>
        typename detail::task_result_type<memfnT>::futureT
        task(ProcessID dest, memfnT memfn, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6,
                const a7T& a7, const a8T& a8, const a9T& a9,
                const TaskAttributes& attr = TaskAttributes()) const
        {
            typedef detail::WorldObjectTaskHelper<Derived, memfnT> task_helper;
            typedef TaskFn<typename task_helper::wrapperT,
                    typename detail::task_arg<a1T>::type,
                    typename detail::task_arg<a2T>::type,
                    typename detail::task_arg<a3T>::type,
                    typename detail::task_arg<a4T>::type,
                    typename detail::task_arg<a5T>::type,
                    typename detail::task_arg<a6T>::type,
                    typename detail::task_arg<a7T>::type,
                    typename detail::task_arg<a8T>::type,
                    typename detail::task_arg<a9T>::type> taskT;
            if (dest == me)
                return world.taskq.add(task_helper::make_task_fn(this, memfn),
                        a1, a2, a3, a4, a5, a6, a7, a8, a9, attr);
            else
                return send_task<taskT>(dest, memfn, am_arg(a1), am_arg(a2),
                        am_arg(a3), am_arg(a4), am_arg(a5), am_arg(a6),
                        am_arg(a7), am_arg(a8), am_arg(a9), attr);
        }

        virtual ~WorldObject() {
            if(initialized())
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
            static inline void store(const BufferOutputArchive& ar, WorldObject<Derived>* const&  ptr) {
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
            static inline void store(const BufferOutputArchive& ar, const WorldObject<Derived>* const& ptr) {
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

#endif // MADNESS_WORLD_WORLDOBJ_H__INCLUDED
