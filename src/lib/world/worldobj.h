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

        /// Info stored for AM method forwarding
        template <typename memfunT>
        struct info {
            typedef Future< MEMFUN_RETURNT(memfunT) > futureT;
            typedef RemoteReference< FutureImpl< MEMFUN_RETURNT(memfunT) > > refT;
            uniqueidT id;
            ProcessID requestor;
            memfunT memfun;
            bool istask;
            refT ref; 

            info() {};

            info(const uniqueidT& id, ProcessID requestor,  memfunT memfun, bool istask, const refT& ref)
                : id(id)
                , requestor(requestor)
                , memfun(memfun)
                , istask(istask)
                , ref(ref)
            {};

            template <typename Archive>
            void serialize(const Archive& ar) {
                ar & archive::wrap_opaque(*this);
            }
        };
    }


    /// Implements most parts of a globally addressable object (via unique ID)

    /// 1) Derived class has DistributedObject<Derived> as a public base class
    /// 2) Derived constructor 
    ///    a) invokes DistributedObject<Derived>(world) constructor
    ///    b) invokes process_pending()
    /// 3) Derived destructor must either be deferred or preceeded by gop.fence()
    ///
    /// This class is deliberately not default constructible and does
    /// not support assignment or copying.  This ensures that each instance
    /// is unique.  Have a look at the DistributedContainer for an example
    /// of wrapping this using the PIMPL idiom and a shared pointer.
    template <class Derived>
    class DistributedObject : public DeferredCleanupInterface {
    private:
        typedef DistributedObject<Derived> objT;
        static std::list<detail::PendingMsg*> pending;   //< Holds pending short/long messages
        World& world;                              //< Think globally act locally
        uniqueidT objid;                           //< Sense of self
        ProcessID me;                              //< Rank of self
    

        /// Actually runs the task or method
        template <typename memfunT, typename futureT>
        static void doit(World& world, bool istask, futureT& result, Derived* obj, memfunT memfun) {
            if (istask) {
                result.set(world.taskq.add(*obj, memfun));
            }
            else {
                TaskMemfunRun< MEMFUN_RETURNT(memfunT) >::run(result, *obj, memfun);
            }
        }


        /// Actually runs the task or method
        template <typename memfunT, typename futureT, typename arg1T>
        static void doit(World& world, bool istask, futureT& result, Derived* obj, memfunT memfun, arg1T& arg1) {
            if (istask) {
                result.set(world.taskq.add(*obj, memfun, arg1));
            }
            else {
                TaskMemfunRun< MEMFUN_RETURNT(memfunT) >::run(result, *obj, memfun, arg1);
            }
        }


        /// Actually runs the task or method
        template <typename memfunT, typename futureT, typename arg1T, typename arg2T>
        static void doit(World& world, bool istask, futureT& result, Derived* obj, memfunT memfun, arg1T& arg1, arg2T& arg2) {
            if (istask) {
                result.set(world.taskq.add(*obj, memfun, arg1, arg2));
            }
            else {
                TaskMemfunRun< MEMFUN_RETURNT(memfunT) >::run(result, *obj, memfun, arg1, arg2);
            }
        }


        /// Actually runs the task or method
        template <typename memfunT, typename futureT, typename arg1T, typename arg2T, typename arg3T>
        static void doit(World& world, bool istask, futureT& result, Derived* obj, memfunT memfun, arg1T& arg1, arg2T& arg2, arg3T& arg3) {
            if (istask) {
                result.set(world.taskq.add(*obj, memfun, arg1, arg2, arg3));
            }
            else {
                TaskMemfunRun< MEMFUN_RETURNT(memfunT) >::run(result, *obj, memfun, arg1, arg2, arg3);
            }
        }


        /// Actually runs the task or method
        template <typename memfunT, typename futureT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        static void doit(World& world, bool istask, futureT& result, Derived* obj, memfunT memfun, arg1T& arg1, arg2T& arg2, arg3T& arg3, arg4T& arg4) {
            if (istask) {
                result.set(world.taskq.add(*obj, memfun, arg1, arg2, arg3, arg4));
            }
            else {
                TaskMemfunRun< MEMFUN_RETURNT(memfunT) >::run(result, *obj, memfun, arg1, arg2, arg3, arg4);
            }
        }


        /// Actually runs the task or method
        template <typename memfunT, typename futureT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        static void doit(World& world, bool istask, futureT& result, Derived* obj, memfunT memfun, arg1T& arg1, arg2T& arg2, arg3T& arg3, arg4T& arg4, arg5T& arg5) {
            if (istask) {
                result.set(world.taskq.add(*obj, memfun, arg1, arg2, arg3, arg4, arg5));
            }
            else {
                TaskMemfunRun< MEMFUN_RETURNT(memfunT) >::run(result, *obj, memfun, arg1, arg2, arg3, arg4, arg5);
            }
        }


        /// Handler for incoming AM with 0 arguments
        template <typename memfunT>
        static void handler(World& world, ProcessID src, const AmArg& arg) {
            detail::info<memfunT> info = arg;
            Derived* obj = world.ptr_from_id<Derived>(info.id);
            if (obj) {
                typename detail::info<memfunT>::futureT result(info.ref);
                doit(world, info.istask, result, obj, info.memfun);
            }
            else {
                pending.push_back(new detail::PendingShortMsg(info.id, handler<memfunT>, src, arg));
            }
        }

        /// Handler for incoming AM with 1 argument
        template <typename memfunT, typename arg1T>
        static void handler(World& world, ProcessID src, void* buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg *) buf;
            detail::info<memfunT> info;
            arg1T arg1;
            arg->unstuff(nbyte, info, arg1);
            Derived* obj = world.ptr_from_id<Derived>(info.id);
            if (obj) {
                typename detail::info<memfunT>::futureT result(info.ref);
                doit(world, info.istask, result, obj, info.memfun, arg1);
            }
            else {
                pending.push_back(new detail::PendingLongMsg(info.id, handler<memfunT,arg1T>, src, buf, nbyte));
            }
        }

        /// Handler for incoming AM with 2 arguments
        template <typename memfunT, typename arg1T, typename arg2T>
        static void handler(World& world, ProcessID src, void* buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg *) buf;
            detail::info<memfunT> info;
            arg1T arg1;
            arg2T arg2;
            arg->unstuff(nbyte, info, arg1, arg2);
            Derived* obj = world.ptr_from_id<Derived>(info.id);
            if (obj) {
                typename detail::info<memfunT>::futureT result(info.ref);
                doit(world, info.istask, result, obj, info.memfun, arg1, arg2);
            }
            else {
                pending.push_back(new detail::PendingLongMsg(info.id, handler<memfunT,arg1T,arg2T>, src, buf, nbyte));
            }
        }

        /// Handler for incoming AM with 3 arguments
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        static void handler(World& world, ProcessID src, void* buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg *) buf;
            detail::info<memfunT> info;
            arg1T arg1;
            arg2T arg2;
            arg3T arg3;
            arg->unstuff(nbyte, info, arg1, arg2, arg3);
            Derived* obj = world.ptr_from_id<Derived>(info.id);
            if (obj) {
                typename detail::info<memfunT>::futureT result(info.ref);
                doit(world, info.istask, result, obj, info.memfun, arg1, arg2, arg3);
            }
            else {
                pending.push_back(new detail::PendingLongMsg(info.id, handler<memfunT,arg1T,arg2T,arg3T>, src, buf, nbyte));
            }
        }

        /// Handler for incoming AM with 4 arguments
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
            if (obj) {
                typename detail::info<memfunT>::futureT result(info.ref);
                doit(world, info.istask, result, obj, info.memfun, arg1, arg2, arg3, arg4);
            }
            else {
                pending.push_back(new detail::PendingLongMsg(info.id, handler<memfunT,arg1T,arg2T,arg3T,arg4T>, src, buf, nbyte));
            }
        }

        /// Handler for incoming AM with 5 arguments
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
            if (obj) {
                typename detail::info<memfunT>::futureT result(info.ref);
                doit(world, info.istask, result, obj, info.memfun, arg1, arg2, arg3, arg4, arg5);
            }
            else {
                pending.push_back(new detail::PendingLongMsg(info.id, handler<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T>, src, buf, nbyte));
            }
        }

        /// Combines implementation of send/task
        template <typename memfunT>
        Future< MEMFUN_RETURNT(memfunT) > 
        forward(ProcessID dest, bool istask, memfunT memfun) {
            Future< MEMFUN_RETURNT(memfunT) > result;
            Derived* obj = static_cast<Derived*>(this);            
            if (dest == me) {
                doit(world, istask, result, obj, memfun);
            }
            else {
                detail::info<memfunT> info(objid, me, memfun, istask, result.remote_ref(world));
                world.am.send(dest,handler<memfunT>,info);
            }
            return result;
        }

        /// Combines implementation of send/task
        template <typename memfunT, typename arg1T>
        Future< MEMFUN_RETURNT(memfunT) > 
        forward(ProcessID dest, bool istask, memfunT memfun, const arg1T& arg1) {
            Future< MEMFUN_RETURNT(memfunT) > result;
            Derived* obj = static_cast<Derived*>(this);            
            if (dest == me) {
                doit(world, istask, result, obj, memfun, arg1);
            }
            else {
                detail::info<memfunT> info(objid, me, memfun, istask, result.remote_ref(world));
                LongAmArg* arg = new LongAmArg;
                std::size_t nbyte = arg->stuff(info,arg1);
                world.am.send_long_managed(dest,handler<memfunT,arg1T>,arg,nbyte);
            }
            return result;
        }


        /// Combines implementation of send/task
        template <typename memfunT, typename arg1T, typename arg2T>
        Future< MEMFUN_RETURNT(memfunT) > 
        forward(ProcessID dest, bool istask, memfunT memfun, const arg1T& arg1, const arg2T& arg2) {
            Future< MEMFUN_RETURNT(memfunT) > result;
            Derived* obj = static_cast<Derived*>(this);            
            if (dest == me) {
                doit(world, istask, result, obj, memfun, arg1, arg2);
            }
            else {
                detail::info<memfunT> info(objid, me, memfun, istask, result.remote_ref(world));
                LongAmArg* arg = new LongAmArg;
                std::size_t nbyte = arg->stuff(info,arg1,arg2);
                world.am.send_long_managed(dest,handler<memfunT,arg1T,arg2T>,arg,nbyte);
            }
            return result;
        }


        /// Combines implementation of send/task
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        Future< MEMFUN_RETURNT(memfunT) > 
        forward(ProcessID dest, bool istask, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3) {
            Future< MEMFUN_RETURNT(memfunT) > result;
            Derived* obj = static_cast<Derived*>(this);            
            if (dest == me) {
                doit(world, istask, result, obj, memfun, arg1, arg2, arg3);
            }
            else {
                detail::info<memfunT> info(objid, me, memfun, istask, result.remote_ref(world));
                LongAmArg* arg = new LongAmArg;
                std::size_t nbyte = arg->stuff(info,arg1,arg2,arg3);
                world.am.send_long_managed(dest,handler<memfunT,arg1T,arg2T,arg3T>,arg,nbyte);
            }
            return result;
        }


        /// Combines implementation of send/task
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        Future< MEMFUN_RETURNT(memfunT) > 
        forward(ProcessID dest, bool istask, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4) {
            Future< MEMFUN_RETURNT(memfunT) > result;
            Derived* obj = static_cast<Derived*>(this);            
            if (dest == me) {
                doit(world, istask, result, obj, memfun, arg1, arg2, arg3, arg4);
            }
            else {
                detail::info<memfunT> info(objid, me, memfun, istask, result.remote_ref(world));
                LongAmArg* arg = new LongAmArg;
                std::size_t nbyte = arg->stuff(info,arg1,arg2,arg3,arg4);
                world.am.send_long_managed(dest,handler<memfunT,arg1T,arg2T,arg3T,arg4T>,arg,nbyte);
            }
            return result;
        }

        /// Combines implementation of send/task
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        Future< MEMFUN_RETURNT(memfunT) > 
        forward(ProcessID dest, bool istask, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5) {
            Future< MEMFUN_RETURNT(memfunT) > result;
            Derived* obj = static_cast<Derived*>(this);            
            if (dest == me) {
                doit(world, istask, result, obj, memfun, arg1, arg2, arg3, arg4, arg5);
            }
            else {
                detail::info<memfunT> info(objid, me, memfun, istask, result.remote_ref(world));
                LongAmArg* arg = new LongAmArg;
                std::size_t nbyte = arg->stuff(info,arg1,arg2,arg3,arg4,arg5);
                world.am.send_long_managed(dest,handler<memfunT,arg1T,arg2T,arg3T,arg4T,arg5T>,arg,nbyte);
            }
            return result;
        }


    protected:

        /// To be called from DERIVED constructor to process pending messages

        /// Cannot call this from the DistributedObject constructor since the
        /// derived class would not yet be fully constructed.
        void process_pending() {
            for (typename std::list<detail::PendingMsg*>::iterator it = pending.begin();
                 it != pending.end();) {
                detail::PendingMsg* p = *it;
                if (p->id == objid) {
                    p->invokehandler(world);
                    delete p;
                    it = pending.erase(it);
                }
                else {
                    ++it;
                }
            }
            if (!pending.empty()) print("Warning ... remaining pending messages!");
        };


    public:
        /// Associates object with globally unique ID
        DistributedObject(World& world) 
            : world(world)
            , objid(world.register_ptr(static_cast<Derived*>(this)))
            , me(world.rank())
        {};
        
        
        /// Returns the globally unique object ID
        const uniqueidT& id() const {
            return objid;
        };

        /// Sends active message to derived class method "returnT (this->*memfun)()"
        template <typename memfunT>
        Future< MEMFUN_RETURNT(memfunT) > 
        send(ProcessID dest, memfunT memfun) {
            return forward(dest, false, memfun);
        }


        /// Sends active message to derived class method "returnT (this->*memfun)(arg1)"
        template <typename memfunT, typename arg1T>
        Future< MEMFUN_RETURNT(memfunT) > 
        send(ProcessID dest, memfunT memfun, const arg1T& arg1) {
            return forward(dest, false, memfun, arg1);
        }

        /// Sends active message to derived class method "returnT (this->*memfun)(arg1,arg2)"
        template <typename memfunT, typename arg1T, typename arg2T>
        Future< MEMFUN_RETURNT(memfunT) > 
        send(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2) {
            return forward(dest, false, memfun, arg1, arg2);
        }


        /// Sends active message to derived class method "returnT (this->*memfun)(arg1,arg2,arg3)"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        Future< MEMFUN_RETURNT(memfunT) > 
        send(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3) {
            return forward(dest, false, memfun, arg1, arg2, arg3);
        }


        /// Sends active message to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4)"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        Future< MEMFUN_RETURNT(memfunT) > 
        send(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4) {
            return forward(dest, false, memfun, arg1, arg2, arg3, arg4);
        }


        /// Sends active message to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4,arg5)"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        Future< MEMFUN_RETURNT(memfunT) > 
        send(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5) {
            return forward(dest, false, memfun, arg1, arg2, arg3, arg4, arg5);
        }


        /// Sends task to derived class method "returnT (this->*memfun)()"
        template <typename memfunT>
        Future< MEMFUN_RETURNT(memfunT) > 
        task(ProcessID dest, memfunT memfun) {
            return forward(dest, true, memfun);
        }
        
        /// Sends task to derived class method "returnT (this->*memfun)(arg1)"
        template <typename memfunT, typename arg1T>
        Future< MEMFUN_RETURNT(memfunT) > 
        task(ProcessID dest, memfunT memfun, const arg1T& arg1) {
            return forward(dest, true, memfun, arg1);
        }

        /// Sends task to derived class method "returnT (this->*memfun)(arg1,arg2)"
        template <typename memfunT, typename arg1T, typename arg2T>
        Future< MEMFUN_RETURNT(memfunT) > 
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2) {
            return forward(dest, true, memfun, arg1, arg2);
        }

        /// Sends task to derived class method "returnT (this->*memfun)(arg1,arg2,arg3)"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        Future< MEMFUN_RETURNT(memfunT) > 
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3) {
            return forward(dest, true, memfun, arg1, arg2, arg3);
        }


        /// Sends task to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4)"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        Future< MEMFUN_RETURNT(memfunT) > 
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4) {
            return forward(dest, true, memfun, arg1, arg2, arg3, arg4);
        }


        /// Sends task to derived class method "returnT (this->*memfun)(arg1,arg2,arg3,arg4,arg5)"
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        Future< MEMFUN_RETURNT(memfunT) > 
        task(ProcessID dest, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5) {
            return forward(dest, true, memfun, arg1, arg2, arg3, arg4, arg5);
        }




        virtual ~DistributedObject(){};
    };

    namespace archive {
        template <class Archive>
        template <class Derived>
        struct ArchiveLoadImpl<Archive,DistributedObject<Derived>*> {
            static inline void load(const Archive& ar, DistributedObject<Derived>*& ptr) {
                uniqueidT id;
                ar & id;
                World* world = World::world_from_id(id.get_world_id());
                MADNESS_ASSERT(world);
                ptr = world->ptr_from_id< DistributedObject<Derived> >(id);
                MADNESS_ASSERT(ptr);
            };
        };
        
        template <class Archive>
        template <class Derived>
        struct ArchiveStoreImpl<Archive,DistributedObject<Derived>*> {
            static inline void store(const Archive& ar, DistributedObject<Derived>*& ptr) {
                ar & ptr->id();
            };
        };
    }

}

#ifdef WORLD_INSTANTIATE_STATIC_TEMPLATES
template <typename Derived> 
std::list<madness::detail::PendingMsg*>
madness::DistributedObject<Derived>::pending;
#endif

#endif
