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

namespace madness {

    namespace detail {

        // Common base class for pending messages to ensure in-order processing
        
        // To eliminate synchronization when a distributed object is first
        // constructed, we buffer pending messages for containers that
        // don't have their id yet registered.
        struct PendingMsg {
            uniqueidT id;
            
            PendingMsg(uniqueidT id) : id(id) {};
            
            virtual void invokehandler(World& world) = 0;
            
            virtual ~PendingMsg(){};
        };
        

        /// Pending short AM
        struct PendingShortMsg : public PendingMsg {
            ProcessID src;
            am_handlerT handler;
            AmArg arg;

            PendingShortMsg(uniqueidT id, am_handlerT handler, ProcessID src, const AmArg& arg)
                : PendingMsg(id), src(src), handler(handler), arg(arg)
            {};

            void invokehandler(World& world) {
                handler(world, src, arg);
            };
        };

        /// Pending long AM
        struct PendingLongMsg : public PendingMsg {
            ProcessID src;
            am_long_handlerT handler;
            size_t nbyte;
            unsigned char* buf;

            PendingLongMsg(uniqueidT id, am_long_handlerT handler, ProcessID src, void *arg, size_t nbyte)
                : PendingMsg(id), src(src), handler(handler), nbyte(nbyte)
            {
                this->buf = new unsigned char[nbyte];
                std::memcpy(buf,arg,nbyte);
            };

            void invokehandler(World& world) {
                handler(world, src, buf, nbyte);
                delete [] buf;
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
            uniqueidT id;
            ProcessID requestor;
            memfunT memfun;
            refT ref; 
            TaskAttributes attr;

            info() {};

            info(const uniqueidT& id, ProcessID requestor,  memfunT memfun, const refT& ref, 
                 const TaskAttributes& attr=TaskAttributes())
                : id(id)
                , requestor(requestor)
                , memfun(memfun)
                , ref(ref)
                , attr(attr)
            {};

            template <typename Archive>
            void serialize(const Archive& ar) {
                ar & archive::wrap_opaque(*this);
            }
        };
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
    template <class Derived>
    class WorldObject : public DeferredCleanupInterface {
    private:
        typedef WorldObject<Derived> objT;
        static std::list<detail::PendingMsg*> pending;   //< Holds pending short/long messages
        World& world;                              //< Think globally act locally
        uniqueidT objid;                           //< Sense of self
        ProcessID me;                              //< Rank of self
        bool ready;                                //< True when process_pending has been called
        bool doing_pending;                        //< Temporary hack to aid in processing pending msg.
    
        // Handler for incoming AM with 0 arguments
        template <typename memfunT>
        static void handler(World& world, ProcessID src, const AmArg& arg) {
            detail::info<memfunT> info = arg;
            Derived* obj = world.ptr_from_id<Derived>(info.id);
            if (is_ready(obj)) {
                typename detail::info<memfunT>::futureT result(info.ref);
                result.set((obj->*info.memfun)());
            }
            else {
                pending.push_back(new detail::PendingShortMsg(info.id, handler<memfunT>, src, arg));
            }
        }

    
        // Handler for incoming AM with 1 argument
        template <typename memfunT, typename arg1T>
        static void handler(World& world, ProcessID src, void* buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg *) buf;
            detail::info<memfunT> info;
            arg1T arg1;
            arg->unstuff(nbyte, info, arg1);
            Derived* obj = world.ptr_from_id<Derived>(info.id);
            if (is_ready(obj)) {
                typename detail::info<memfunT>::futureT result(info.ref);
                result.set((obj->*info.memfun)(arg1));
            }
            else {
                pending.push_back(new detail::PendingLongMsg(info.id, handler<memfunT,arg1T>, src, buf, nbyte));
            }
        }

        // Handler for incoming AM with 2 arguments
        template <typename memfunT, typename arg1T, typename arg2T>
        static void handler(World& world, ProcessID src, void* buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg *) buf;
            detail::info<memfunT> info;
            arg1T arg1;
            arg2T arg2;
            arg->unstuff(nbyte, info, arg1, arg2);
            Derived* obj = world.ptr_from_id<Derived>(info.id);
            if (is_ready(obj)) {
                typename detail::info<memfunT>::futureT result(info.ref);
                result.set((obj->*info.memfun)(arg1,arg2));
            }
            else {
                pending.push_back(new detail::PendingLongMsg(info.id, handler<memfunT,arg1T,arg2T>, src, buf, nbyte));
            }
        }

        // Handler for incoming AM with 3 arguments
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        static void handler(World& world, ProcessID src, void* buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg *) buf;
            detail::info<memfunT> info;
            arg1T arg1;
            arg2T arg2;
            arg3T arg3;
            arg->unstuff(nbyte, info, arg1, arg2, arg3);
            Derived* obj = world.ptr_from_id<Derived>(info.id);
            if (is_ready(obj)) {
                typename detail::info<memfunT>::futureT result(info.ref);
                result.set((obj->*info.memfun)(arg1,arg2,arg3));
            }
            else {
                pending.push_back(new detail::PendingLongMsg(info.id, handler<memfunT,arg1T,arg2T,arg3T>, src, buf, nbyte));
            }
        }

        // Handler for incoming AM with 4 arguments
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        static void handler(World& world, ProcessID src, void* buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg *) buf;
            detail::info<memfunT> info;
            arg1T arg1;
            arg2T arg2;
            arg3T arg3;
            arg4T arg4;
            arg->unstuff(nbyte, info, arg1, arg2, arg3, arg4);
            Derived* obj = world.ptr_from_id<Derived>(info.id);
            if (is_ready(obj)) {
                typename detail::info<memfunT>::futureT result(info.ref);
                result.set((obj->*info.memfun)(arg1,arg2,arg3,arg4));
            }
            else {
                pending.push_back(new detail::PendingLongMsg(info.id, handler<memfunT,arg1T,arg2T,arg3T,arg4T>, src, buf, nbyte));
            }
        }

        // Handler for incoming AM with 5 arguments
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        static void handler(World& world, ProcessID src, void* buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg *) buf;
            detail::info<memfunT> info;
            arg1T arg1;
            arg2T arg2;
            arg3T arg3;
            arg4T arg4;
            arg5T arg5;
            arg->unstuff(nbyte, info, arg1, arg2, arg3, arg4, arg5);
            Derived* obj = world.ptr_from_id<Derived>(info.id);
            if (is_ready(obj)) {
                typename detail::info<memfunT>::futureT result(info.ref);
                result.set((obj->*info.memfun)(arg1,arg2,arg3,arg4,arg5));
            }
            else {
                pending.push_back(new detail::PendingLongMsg(info.id, handler<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T>, src, buf, nbyte));
            }
        }

        // Handler for incoming task with 0 arguments
        template <typename memfunT>
        static void handler_task(World& world, ProcessID src, const AmArg& arg) {
            detail::info<memfunT> info = arg;
            Derived* obj = world.ptr_from_id<Derived>(info.id);
            if (is_ready(obj)) {
                typename detail::info<memfunT>::futureT result(info.ref);
                world.taskq.add(new TaskMemfun<memfunT>(result,*obj,info.memfun,info.attr));
            }
            else {
                pending.push_back(new detail::PendingShortMsg(info.id, handler_task<memfunT>, src, arg));
            }
        }

    
        // Handler for incoming task with 1 argument
        template <typename memfunT, typename arg1T>
        static void handler_task(World& world, ProcessID src, void* buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg *) buf;
            detail::info<memfunT> info;
            arg1T arg1;
            arg->unstuff(nbyte, info, arg1);
            Derived* obj = world.ptr_from_id<Derived>(info.id);
            if (is_ready(obj)) {
                typename detail::info<memfunT>::futureT result(info.ref);
                world.taskq.add(new TaskMemfun<memfunT>(result,*obj,info.memfun,arg1,info.attr));
            }
            else {
                pending.push_back(new detail::PendingLongMsg(info.id, handler_task<memfunT,arg1T>, src, buf, nbyte));
            }
        }

        // Handler for incoming task with 2 arguments
        template <typename memfunT, typename arg1T, typename arg2T>
        static void handler_task(World& world, ProcessID src, void* buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg *) buf;
            detail::info<memfunT> info;
            arg1T arg1;
            arg2T arg2;
            arg->unstuff(nbyte, info, arg1, arg2);
            Derived* obj = world.ptr_from_id<Derived>(info.id);
            if (is_ready(obj)) {
                typename detail::info<memfunT>::futureT result(info.ref);
                world.taskq.add(new TaskMemfun<memfunT>(result,*obj,info.memfun,arg1,arg2,info.attr));
            }
            else {
                pending.push_back(new detail::PendingLongMsg(info.id, handler_task<memfunT,arg1T,arg2T>, src, buf, nbyte));
            }
        }

        // Handler for incoming task with 3 arguments
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        static void handler_task(World& world, ProcessID src, void* buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg *) buf;
            detail::info<memfunT> info;
            arg1T arg1;
            arg2T arg2;
            arg3T arg3;
            arg->unstuff(nbyte, info, arg1, arg2, arg3);
            Derived* obj = world.ptr_from_id<Derived>(info.id);
            if (is_ready(obj)) {
                typename detail::info<memfunT>::futureT result(info.ref);
                world.taskq.add(new TaskMemfun<memfunT>(result,*obj,info.memfun,arg1,arg2,arg3,info.attr));
            }
            else {
                pending.push_back(new detail::PendingLongMsg(info.id, handler_task<memfunT,arg1T,arg2T,arg3T>, src, buf, nbyte));
            }
        }

        // Handler for incoming task with 4 arguments
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        static void handler_task(World& world, ProcessID src, void* buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg *) buf;
            detail::info<memfunT> info;
            arg1T arg1;
            arg2T arg2;
            arg3T arg3;
            arg4T arg4;
            arg->unstuff(nbyte, info, arg1, arg2, arg3, arg4);
            Derived* obj = world.ptr_from_id<Derived>(info.id);
            if (is_ready(obj)) {
                typename detail::info<memfunT>::futureT result(info.ref);
                world.taskq.add(new TaskMemfun<memfunT>(result,*obj,info.memfun,arg1,arg2,arg3,arg4,info.attr));
            }
            else {
                pending.push_back(new detail::PendingLongMsg(info.id, handler_task<memfunT,arg1T,arg2T,arg3T,arg4T>, src, buf, nbyte));
            }
        }

        // Handler for incoming task with 5 arguments
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        static void handler_task(World& world, ProcessID src, void* buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg *) buf;
            detail::info<memfunT> info;
            arg1T arg1;
            arg2T arg2;
            arg3T arg3;
            arg4T arg4;
            arg5T arg5;
            arg->unstuff(nbyte, info, arg1, arg2, arg3, arg4, arg5);
            Derived* obj = world.ptr_from_id<Derived>(info.id);
            if (is_ready(obj)) {
                typename detail::info<memfunT>::futureT result(info.ref);
                world.taskq.add(new TaskMemfun<memfunT>(result,*obj,info.memfun,arg1,arg2,arg3,arg4,arg5,info.attr));
            }
            else {
                pending.push_back(new detail::PendingLongMsg(info.id, handler_task<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T>, src, buf, nbyte));
            }
        }

        // Returns true if obj is not NULL pointer and is also ready for incoming messages
        static bool is_ready(Derived* obj) {
            if (obj) {
                WorldObject* p = static_cast<WorldObject*>(obj);
                if (p->doing_pending) {
                    p->doing_pending = false;
                    return true;
                }
                else {
                    return p->ready;
                }
            }
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
            // While processing the pendingQ, an incoming AM may still be appended to the
            // Q so we must keep processing until all is done.


            // This loop and use of ndone/ready/doing_pending requires careful reconsideration
            // when going multi-threaded.
            int ndone = 1;
            while (ndone) {
                ndone = 0;
                for (typename std::list<detail::PendingMsg*>::iterator it = pending.begin();
                     it != pending.end();) {
                    detail::PendingMsg* p = *it;
                    if (p->id == objid) {
                        // Ugh!  Override ready=false just for this message.
                        // is_ready() looks at doing_pending and resets it
                        // to false.  This relies on no call to polling
                        // between here and use of this flag, and also
                        // that we are single threaded ... something else
                        // will have to happen for multi-core.
                        doing_pending = true;
                        p->invokehandler(world);
                        delete p;
                        // Having the erase after the invoke ensures any 
                        // incoming AM are also digested in this loop.
                        // But only in polling implementation.
                        it = pending.erase(it);
                        ndone ++;
                    }
                    else {
                        ++it;
                    }
                }
            }

            ready = true;
        };


    public:
        /// Associates object with globally unique ID

        /// !! The derived class MUST call process_pending both to 
        /// process any messages that arrived prior to construction
        /// and to enable processing of future messages.
        WorldObject(World& world) 
            : world(world)
            , objid(world.register_ptr(static_cast<Derived*>(this)))
            , me(world.rank())
            , ready(false)
            , doing_pending(false)
        {};
        
        
        /// Returns the globally unique object ID
        const uniqueidT& id() const {
            return objid;
        };


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
                world.am.send(dest,handler<memfunT>,info);
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
                LongAmArg* arg = new LongAmArg;
                std::size_t nbyte = arg->stuff(info,arg1);
                world.am.send_long_managed(dest,handler<memfunT,arg1T>,arg,nbyte);
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
                LongAmArg* arg = new LongAmArg;
                std::size_t nbyte = arg->stuff(info,arg1,arg2);
                world.am.send_long_managed(dest,handler<memfunT,arg1T,arg2T>,arg,nbyte);
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
                LongAmArg* arg = new LongAmArg;
                std::size_t nbyte = arg->stuff(info,arg1,arg2,arg3);
                world.am.send_long_managed(dest,handler<memfunT,arg1T,arg2T,arg3T>,arg,nbyte);
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
                LongAmArg* arg = new LongAmArg;
                std::size_t nbyte = arg->stuff(info,arg1,arg2,arg3,arg4);
                world.am.send_long_managed(dest,handler<memfunT,arg1T,arg2T,arg3T,arg4T>,arg,nbyte);
            }
            return result;
        }

        /// Sends active message to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4,arg5)"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > 
        send(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5) {
            Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
            Derived* obj = static_cast<Derived*>(this);            
            if (dest == me) {
                result.set((obj->*memfun)(arg1,arg2,arg3,arg4,arg5));
            }
            else {
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world));
                LongAmArg* arg = new LongAmArg;
                std::size_t nbyte = arg->stuff(info,arg1,arg2,arg3,arg4,arg5);
                world.am.send_long_managed(dest,handler<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T>,arg,nbyte);
            }
            return result;
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
                world.am.send(dest,handler_task<memfunT>,info);
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
                LongAmArg* arg = new LongAmArg;
                std::size_t nbyte = arg->stuff(info,arg1);
                world.am.send_long_managed(dest,handler_task<memfunT,arg1T>,arg,nbyte);
                return result;
            }
        }

        /// Sends task to derived class method "returnT (this->*memfun)(arg1,arg2)"
        template <typename memfunT, typename arg1T, typename arg2T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > 
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const TaskAttributes& attr = TaskAttributes()) {
            if (dest == me) {
                return world.taskq.add(*static_cast<Derived*>(this), memfun, arg1, arg2, attr);
            }
            else {
                Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world), attr);
                LongAmArg* arg = new LongAmArg;
                std::size_t nbyte = arg->stuff(info,arg1, arg2);
                world.am.send_long_managed(dest,handler_task<memfunT,arg1T,arg2T>,arg,nbyte);
                return result;
            }
        }

        /// Sends task to derived class method "returnT (this->*memfun)(arg1,arg2,arg3)"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > 
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const TaskAttributes& attr = TaskAttributes()) {
            if (dest == me) {
                return world.taskq.add(*static_cast<Derived*>(this), memfun, arg1, arg2, arg3, attr);
            }
            else {
                Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world), attr);
                LongAmArg* arg = new LongAmArg;
                std::size_t nbyte = arg->stuff(info,arg1, arg2, arg3);
                world.am.send_long_managed(dest,handler_task<memfunT,arg1T,arg2T,arg3T>,arg,nbyte);
                return result;
            }
        }

        /// Sends task to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4)"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > 
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const TaskAttributes& attr = TaskAttributes()) {
            if (dest == me) {
                return world.taskq.add(*static_cast<Derived*>(this), memfun, arg1, arg2, arg3, arg4, attr);
            }
            else {
                Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world), attr);
                LongAmArg* arg = new LongAmArg;
                std::size_t nbyte = arg->stuff(info,arg1, arg2, arg3, arg4);
                world.am.send_long_managed(dest,handler_task<memfunT,arg1T,arg2T,arg3T,arg4T>,arg,nbyte);
                return result;
            }
        }

        /// Sends task to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4,arg5)"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > 
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const TaskAttributes& attr = TaskAttributes()) {
            if (dest == me) {
                return world.taskq.add(*static_cast<Derived*>(this), memfun, arg1, arg2, arg3, arg4, arg5, attr);
            }
            else {
                Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > result;
                detail::info<memfunT> info(objid, me, memfun, result.remote_ref(world), attr);
                LongAmArg* arg = new LongAmArg;
                std::size_t nbyte = arg->stuff(info,arg1, arg2, arg3, arg4, arg5);
                world.am.send_long_managed(dest,handler_task<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T>,arg,nbyte);
                return result;
            }
        }

        virtual ~WorldObject(){
            world.unregister_ptr(static_cast<Derived*>(this));
        };
    };

    namespace archive {
        template <class Archive, class Derived>
        struct ArchiveLoadImpl<Archive,WorldObject<Derived>*> {
            static inline void load(const Archive& ar, WorldObject<Derived>*& ptr) {
                uniqueidT id;
                ar & id;
                World* world = World::world_from_id(id.get_world_id());
                MADNESS_ASSERT(world);
                ptr = world->ptr_from_id< WorldObject<Derived> >(id);
                if (!ptr) MADNESS_EXCEPTION("WorldObj: remote operation attempting to use a locally uninitialized object",0);
            };
        };
        
        template <class Archive, class Derived>
        struct ArchiveStoreImpl<Archive,WorldObject<Derived>*> {
            static inline void store(const Archive& ar, const WorldObject<Derived>*const& ptr) {
                ar & ptr->id();
            };
        };
    }

}

#ifdef WORLD_INSTANTIATE_STATIC_TEMPLATES
template <typename Derived> 
std::list<madness::detail::PendingMsg*>
madness::WorldObject<Derived>::pending;
#endif

#endif
