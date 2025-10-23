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

/**
 \file world_object.h
 \brief Defines and implements \c WorldObject.
 \ingroup world_object
*/

#ifndef MADNESS_WORLD_WORLD_OBJECT_H__INCLUDED
#define MADNESS_WORLD_WORLD_OBJECT_H__INCLUDED

#include <madness/world/thread.h>
#include <madness/world/world_task_queue.h>

#include <array>
#include <atomic>
#include <cstddef>
#include <type_traits>

/// \addtogroup world_object
/// @{

namespace madness {

    template <typename> class WorldObject;

    namespace detail {

        /// Common base class for pending messages to ensure in-order processing.

        /// To eliminate synchronization when a distributed object is first
        /// constructed, we buffer pending messages for containers that
        /// don't have their ID yet registered.
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

        /// \todo Brief description needed.

        /// We cannot use the normal task forwarding stuff here because pending
        /// messages can creating tasks are must be correctly handled. The
        /// following code does not easily merge with the send handlers since
        /// the task layer is more restrictive on the copy capability of
        /// arguments.
        ///
        /// It is also annoying that info needs to be broken into two parts so
        /// that it \c id and \c ref are properly serialized. We need to have
        /// \c id correctly aligned via \c opaque_wrap, but \c ref cannot be
        /// serialized that way. Thus we break the class into two parts.
        ///
        /// Info stored for AM method forwarding.
        /// \tparam memfunT Description needed.
        /// \todo Verify & complete; what is AM?
        template <typename memfunT>
        struct info_base {
            using memfunT_rel_ptr = decltype(archive::to_rel_memfn_ptr(std::declval<memfunT>()));
            // id must be at front ... see peek.
            uniqueidT id; ///< \todo Description needed. Context with the "see peek" comment above?
            ProcessID requestor; ///< \todo Description needed.
            memfunT_rel_ptr memfun_rel_ptr; ///< \todo Description needed.
            TaskAttributes attr; ///< \todo Description needed.

            /// \return the (absolute) member function pointer
            memfunT memfun() const {
              return archive::to_abs_memfn_ptr<memfunT>(memfun_rel_ptr);
            }

        protected:

            info_base() {}

            /// \todo Constructor that [brief description needed].

            /// \todo Descriptions needed.
            /// \param[in] id Description needed.
            /// \param[in] requestor Description needed.
            /// \param memfun Verify: The member function to be invoked for the task.
            /// \param[in] attr Description needed.
            info_base(const uniqueidT& id, ProcessID requestor, memfunT memfun,
                 const TaskAttributes& attr=TaskAttributes())
                    : id(id)
                    , requestor(requestor)
                    , memfun_rel_ptr(archive::to_rel_memfn_ptr(memfun))
                    , attr(attr) {}

            /// Serializes a \c info_base for I/O.

            /// \tparam Archive The type of I/O archive.
            /// \param[in,out] ar The I/O archive.
            template <typename Archive>
            void serialize(const Archive& ar) {
                ar & archive::wrap_opaque(*this); // Must be opaque ... see peek.
            }
        }; // struct info_base

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam memfunT Verify: Signature of the member function in the derived class to be invoked for the task.
        template <typename memfunT>
        struct info : public info_base<memfunT> {
            /// Future for a return value of the memory function. \todo Verify.
            typedef Future< REMFUTURE(MEMFUN_RETURNT(memfunT)) > futureT;
            /// \todo Description needed.
            typedef RemoteReference< FutureImpl< REMFUTURE(MEMFUN_RETURNT(memfunT)) > > refT;

            refT ref; ///< \todo Description needed.

            info() : info_base<memfunT>() {}

            /// \todo Constructor that [brief description needed].

            /// \todo Descriptions needed.
            /// \param[in] arg Description needed.
            info(const AmArg& arg) :
                info_base<memfunT>()
            {
                arg & *this;
            }

            /// \todo Constructor that [brief description needed].

            /// \todo Descriptions needed.
            /// \param[in] id Description needed.
            /// \param[in] requestor Description needed.
            /// \param memfun Verify: The member function to be invoked for the task.
            /// \param[in] ref Description needed.
            /// \param[in] attr Description needed.
            info(const uniqueidT& id, ProcessID requestor, memfunT memfun,
                 const refT& ref, const TaskAttributes& attr=TaskAttributes())
                : info_base<memfunT>(id, requestor, memfun, attr), ref(ref)
            {}

            /// Serializes a \c info for I/O.

            /// \tparam Archive the type of I/O archive.
            /// \param[in] ar The I/O archive.
            template <typename Archive>
            void serialize(const Archive& ar) {
                info_base<memfunT>::serialize(ar);
                ar & ref;
            }
        }; // struct info

        /// Extract the unique object ID from an incoming active message header.

        /// We deserialize the header and all arguments at the same
        /// time to simplify the code. However, it is common that
        /// when sending a message to an item in a container to
        /// include a pointer to the container itself. But this
        /// breaks if the container is not initialized since the
        /// deserialization throws if the object is not initialized
        /// (which seems preferable to a hidden race condition). Hence,
        /// we use this routine to extract the unique ID from the very
        /// front of the \c info structure. For efficiency we here rely
        /// upon the serialization of \c info being opaque and the
        /// ID being at the front of \c info.
        ///
        /// \todo Verify parameter description.
        /// \param[in] arg The active message header.
        static inline const uniqueidT& peek(const AmArg& arg) {
            return *((uniqueidT*)(arg.buf()));
        }


        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam objT Description needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam Enabler Description needed.
        template <typename objT, typename memfnT, typename Enabler = void>
        struct WorldObjectTaskHelper {
            /// \todo Description needed.
            typedef typename std::conditional<memfunc_traits<memfnT>::constness,
                    const objT*, objT*>::type ptrT;

            /// \todo Description needed.
            typedef MemFuncWrapper<ptrT, memfnT, typename result_of<memfnT>::type> wrapperT;


            /// \todo Brief description needed.

            /// \todo Descriptions needed.
            /// \param[in] obj Description needed.
            /// \param memfn Verify: The member function to be invoked for the task.
            /// \return Description needed.
            static wrapperT make_task_fn(const objT* obj, memfnT memfn) {
                return wrapperT(const_cast<ptrT>(obj), memfn);
            }


            /// \todo Brief description needed.

            /// \todo Descriptions needed.
            /// \param[in] obj Description needed.
            /// \param memfn Verify: The member function to be invoked for the task.
            /// \return Description needed.
            static wrapperT make_task_fn(objT* obj, memfnT memfn) {
                return wrapperT(const_cast<ptrT>(obj), memfn);
            }


            /// \todo Brief description needed.

            /// \todo Descriptions needed.
            /// \param[in] obj Description needed.
            /// \param memfn Verify: The member function to be invoked for the task.
            /// \return Description needed.
            static wrapperT make_task_fn(const WorldObject<objT>* obj, memfnT memfn) {
                return make_task_fn(static_cast<const objT*>(obj), memfn);
            }

            /// \todo Brief description needed.

            /// \todo Descriptions needed.
            /// \param[in] obj Description needed.
            /// \param memfn Verify: The member function to be invoked for the task.
            /// \return Description needed.
            static wrapperT make_task_fn(WorldObject<objT>* obj, memfnT memfn) {
                return make_task_fn(static_cast<objT*>(obj), memfn);
            }
        }; // struct WorldObjectTaskHelper


#ifndef MADNESS_DISABLE_SHARED_FROM_THIS
        // Disable the use of std::enable_shared_from_this if we are using MADNESS's
        // implementation since weak_ptr is not fully implemented.


        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam objT Description needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        template <typename objT, typename memfnT>
        struct WorldObjectTaskHelper<objT, memfnT,
                typename std::enable_if< std::is_base_of<std::enable_shared_from_this<objT>, objT>::value >::type>
        {
            /// \todo Description needed.
            typedef typename std::conditional<memfunc_traits<memfnT>::constness,
                    const objT*, objT*>::type ptrT;

            /// \todo Description needed.
            typedef typename std::conditional<memfunc_traits<memfnT>::constness,
                    std::shared_ptr<const objT>, std::shared_ptr<objT> >::type shared_ptrT;

            /// \todo Description needed.
            typedef MemFuncWrapper<shared_ptrT, memfnT, typename result_of<memfnT>::type> wrapperT;


            /// \todo Brief description needed.

            /// \todo Descriptions needed.
            /// \param[in] obj Description needed.
            /// \param memfn Verify: The member function to be invoked for the task.
            /// \return Description needed.
            static wrapperT make_task_fn(const objT* obj, memfnT memfn) {
                return wrapperT(const_cast<ptrT>(obj)->shared_from_this(), memfn);
            }


            /// \todo Brief description needed.

            /// \todo Descriptions needed.
            /// \param[in] obj Description needed.
            /// \param memfn Verify: The member function to be invoked for the task.
            /// \return Description needed.
            static wrapperT make_task_fn(objT* obj, memfnT memfn) {
                return wrapperT(const_cast<ptrT>(obj)->shared_from_this(), memfn);
            }


            /// \todo Brief description needed.

            /// \todo Descriptions needed.
            /// \param[in] obj Description needed.
            /// \param memfn Verify: The member function to be invoked for the task.
            /// \return Description needed.
            static wrapperT make_task_fn(const WorldObject<objT>* obj, memfnT memfn) {
                return make_task_fn(static_cast<const objT*>(obj), memfn);
            }


            /// \todo Brief description needed.

            /// \todo Descriptions needed.
            /// \param[in] obj Description needed.
            /// \param memfn Verify: The member function to be invoked for the task.
            /// \return Description needed.
            static wrapperT make_task_fn(WorldObject<objT>* obj, memfnT memfn) {
                return make_task_fn(static_cast<objT*>(obj), memfn);
            }
        }; // struct WorldObjectTaskHelper

#endif // MADNESS_DISABLE_SHARED_FROM_THIS

    } // namespace detail

    /// Base class for WorldObject, useful for introspection
    struct WorldObjectBase {
      virtual ~WorldObjectBase() = default;
    public:
        virtual World& get_world() const = 0;
    };

    /// Implements most parts of a globally addressable object (via unique ID).

    /// This class is deliberately not default constructible and does
    /// not support assignment or copying. This ensures that each instance
    /// is unique. Have a look at \c madness::WorldContainer for an example
    /// of wrapping this using the PIMPL idiom and a shared pointer.
    ///
    /// When deriving classes:
    /// -# Derived class has `WorldObject<Derived>` as a public base class.
    /// -# Derived constructor:
    ///    -# invokes `WorldObject<Derived>(world)` constructor.
    ///    -# invokes `process_pending()`.
    /// -# Derived destructor must either be deferred or preceeded by `gop.fence()`.
    /// -# Derived class must have at least one virtual function for serialization
    ///    of derived class pointers to be cast to the appropriate type.
    ///
    /// Note that \c world is exposed for convenience as a public data member.
    /// \tparam Derived The derived class. \c WorldObject is a curiously
    ///     recurring template pattern.
    template <class Derived>
    class WorldObject : public WorldObjectBase {
    public:
        /// \todo Description needed.
        typedef WorldObject<Derived> objT;

        // copy ctor must be enabled to permit RVO; in C++17 will not need this
        WorldObject(const WorldObject& other) : world(other.world) { abort(); }
        // no copy
        WorldObject& operator=(const WorldObject&) = delete;

    private:
        /// \todo Description needed.
        typedef std::list<detail::PendingMsg> pendingT;

        /// \todo Description needed.
        typedef detail::voidT voidT;

        World& world; ///< The \c World this object belongs to. (Think globally, act locally).

        // The order here matters in a multi-threaded world
        volatile bool ready; ///< True if ready to rock 'n roll.
        ProcessID me; ///< Rank of self.
        uniqueidT objid; ///< Sense of self.


        inline static Spinlock pending_mutex; ///< \todo Description needed.
        inline static volatile pendingT pending; ///< Buffer for pending messages.


        /// \todo Complete: Determine if [unknown] is ready (for ...).

        /// The slightly convoluted logic is to ensure ordering when
        /// processing pending messages. If a new message arrives
        /// while processing incoming messages it must be queued.
        ///
        /// - If the object does not exist ---> not ready.
        /// - If the object exists and is ready ---> ready.
        /// - If the object exists and is not ready then
        ///      - if we are doing a queued/pending message --> ready.
        ///      - else this is a new message --> not ready.
        ///
        /// \param[in] id Description needed.
        /// \param[in,out] obj Description needed.
        /// \param[in] arg Description needed.
        /// \param[in,out] ptr Description needed.
        /// \return Description needed.
        /// \todo Parameter/return descriptions needed.
        static bool is_ready(const uniqueidT& id, objT*& obj, const AmArg& arg, am_handlerT ptr) {
          std::optional<Derived *> opt_obj =
              arg.get_world()->template ptr_from_id<Derived>(id);
          if (opt_obj) {
            // if opt_obj == nullptr, then this ID has already been deregistered
            MADNESS_ASSERT(*opt_obj != nullptr);
            obj = static_cast<objT *>(*opt_obj);
          } else
            obj = nullptr;

          if (obj) {
            if (obj->ready || arg.is_pending())
              return true;
          }

          MADNESS_PRAGMA_CLANG(diagnostic push)
          MADNESS_PRAGMA_CLANG(diagnostic ignored "-Wundefined-var-template")

          ScopedMutex<Spinlock> lock(pending_mutex); // BEGIN CRITICAL SECTION

          if (!obj) {
            std::optional<Derived *> opt_obj =
                arg.get_world()->template ptr_from_id<Derived>(id);
            if (opt_obj) {
              // if opt_obj == nullptr, then this ID has already been deregistered
              MADNESS_ASSERT(*opt_obj != nullptr);
              obj = static_cast<objT *>(*opt_obj);
            } else
              obj = nullptr;
          }

          if (obj) {
            if (obj->ready || arg.is_pending())
              return true; // END CRITICAL SECTION
          }
          const_cast<AmArg &>(arg).set_pending();
          const_cast<pendingT &>(pending).push_back(
              detail::PendingMsg(id, ptr, arg));

          MADNESS_PRAGMA_CLANG(diagnostic pop)

          return false; // END CRITICAL SECTION
        }

        /// Handler for an incoming AM.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam arg1T Type of argument 1.
        /// \tparam arg2T Type of argument 2.
        /// \tparam arg3T Type of argument 3.
        /// \tparam arg4T Type of argument 4.
        /// \tparam arg5T Type of argument 5.
        /// \tparam arg6T Type of argument 6.
        /// \tparam arg7T Type of argument 7.
        /// \tparam arg8T Type of argument 8.
        /// \tparam arg9T Type of argument 9.
        /// \param[in] arg Description needed.
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
		MADNESS_PRAGMA_CLANG(diagnostic push)
                MADNESS_PRAGMA_CLANG(diagnostic ignored "-Wuninitialized-const-reference")
                arg & info & arg1 & arg2 & arg3 & arg4 & arg5 & arg6 & arg7 & arg8 & arg9;
		MADNESS_PRAGMA_CLANG(diagnostic pop)
                typename detail::info<memfnT>::futureT result(info.ref);
                detail::run_function(result, task_helper::make_task_fn(obj, info.memfun()),
                        arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);
            }
        }


        /// Handler for remote arguments.

        /// \todo Descriptions needed.
        /// \tparam taskT Description needed.
        /// \param[in] arg Description needed.
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
                        task_helper::make_task_fn(obj, info.memfun()), info.attr, input_arch);

                // Add task to queue
                arg.get_world()->taskq.add(task);
            }
        }


        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam T Description needed.
        /// \param[in] f Description needed.
        /// \return Description needed.
        template <typename T>
        static inline const T& am_arg(const Future<T>& f) {
            MADNESS_ASSERT(f.probe()); // Cannot serialize unassigned futures
            return f.get();
        }


        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam T Description needed.
        /// \param[in] t Description needed.
        /// \return Description needed.
        template <typename T>
        static inline const T& am_arg(const T& t) {
            return t;
        }


        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \tparam a6T Type of argument 6.
        /// \tparam a7T Type of argument 7.
        /// \tparam a8T Type of argument 8.
        /// \tparam a9T Type of argument 9.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param a1 Argument 1.
        /// \param a2 Argument 2.
        /// \param a3 Argument 3.
        /// \param a4 Argument 4.
        /// \param a5 Argument 5.
        /// \param a6 Argument 6.
        /// \param a7 Argument 7.
        /// \param a8 Argument 8.
        /// \param a9 Argument 9.
        /// \return Description needed.
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


        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam taskT Description needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \tparam a6T Type of argument 6.
        /// \tparam a7T Type of argument 7.
        /// \tparam a8T Type of argument 8.
        /// \tparam a9T Type of argument 9.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param a1 Argument 1.
        /// \param a2 Argument 2.
        /// \param a3 Argument 3.
        /// \param a4 Argument 4.
        /// \param a5 Argument 5.
        /// \param a6 Argument 6.
        /// \param a7 Argument 7.
        /// \param a8 Argument 8.
        /// \param a9 Argument 9.
        /// \param attr Description needed.
        /// \return Description needed.
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
                    new_am_arg(info, a1, a2, a3, a4, a5, a6, a7, a8, a9), RMI::ATTR_UNORDERED);

            return result;
        }

    protected:

        /// To be called from \em derived constructor to process pending messages.

        /// Cannot call this from the \c WorldObject constructor since the
        /// derived class would not yet be fully constructed.
        ///
        /// \attention No incoming messages are processed until this routine is
        /// invoked; the derived class may rely upon a well defined state
        /// until this routine is invoked.
        void process_pending() {
            // Messages may be arriving while we are processing the
            // pending queue.  To maximize concurrency copy messages
            // out of queue before processing outside critical section.
            //int ndone = 0;
            while (!ready) {
                pendingT tmp;

                MADNESS_PRAGMA_CLANG(diagnostic push)
                MADNESS_PRAGMA_CLANG(diagnostic ignored "-Wundefined-var-template")

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

                MADNESS_PRAGMA_CLANG(diagnostic pop)

                while (tmp.size()) {
                    tmp.front().invokehandler();
                    tmp.pop_front();
                    //++ndone;
                }
            }
            //if (ndone) std::cout << world.rank() << ":pending:" << ndone << std::endl;
        }


    public:
        /// \brief Constructor that associates an object (via the derived class)
        ///     with a globally unique ID.

        /// \attention The derived class MUST call \c process_pending from
        /// its constructor to both
        /// -# process any messages that arrived prior to construction.
        /// -# to enable processing of future messages.
        /// \param[in,out] world The \c World encapsulating the \"global\" domain.
        WorldObject(World& world)
                : world(world)
                , ready(false)
                , me(world.rank())
                , objid(world.register_ptr(static_cast<Derived*>(this))) {};


        /// Returns the globally unique object ID.
        const uniqueidT& id() const {
            return objid;
        }


        /// Returns a reference to the \c world.
        World& get_world() const {
            return const_cast<WorldObject<Derived>*>(this)->world;
        }


        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \return Description needed.
        template <typename memfnT>
        typename detail::task_result_type<memfnT>::futureT
        send(ProcessID dest, memfnT memfn) const {
            return send_am(dest, memfn, voidT::value, voidT::value, voidT::value,
                    voidT::value, voidT::value, voidT::value, voidT::value,
                    voidT::value, voidT::value);
        }


        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam a1T Type of argument 1.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param a1 Argument 1.
        /// \return Description needed.
        template <typename memfnT, typename a1T>
        typename detail::task_result_type<memfnT>::futureT
        send(ProcessID dest, memfnT memfn, const a1T& a1) const {
            return send_am(dest, memfn, am_arg(a1), voidT::value, voidT::value,
                    voidT::value, voidT::value, voidT::value, voidT::value,
                    voidT::value, voidT::value);
        }


        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param a1 Argument 1.
        /// \param a2 Argument 2.
        /// \return Description needed.
        template <typename memfnT, typename a1T, typename a2T>
        typename detail::task_result_type<memfnT>::futureT
        send(ProcessID dest, memfnT memfn, const a1T& a1, const a2T& a2) const {
            return send_am(dest, memfn, am_arg(a1), am_arg(a2), voidT::value,
                    voidT::value, voidT::value, voidT::value, voidT::value,
                    voidT::value, voidT::value);
        }


        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param a1 Argument 1.
        /// \param a2 Argument 2.
        /// \param a3 Argument 3.
        /// \return Description needed.
        template <typename memfnT, typename a1T, typename a2T, typename a3T>
        typename detail::task_result_type<memfnT>::futureT
        send(ProcessID dest, memfnT memfn, const a1T& a1, const a2T& a2,
                const a3T& a3) const
        {
            return send_am(dest, memfn, am_arg(a1), am_arg(a2), am_arg(a3),
                    voidT::value, voidT::value, voidT::value, voidT::value,
                    voidT::value, voidT::value);
        }


        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param a1 Argument 1.
        /// \param a2 Argument 2.
        /// \param a3 Argument 3.
        /// \param a4 Argument 4.
        /// \return Description needed.
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


        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param a1 Argument 1.
        /// \param a2 Argument 2.
        /// \param a3 Argument 3.
        /// \param a4 Argument 4.
        /// \param a5 Argument 5.
        /// \return Description needed.
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


        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \tparam a6T Type of argument 6.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param a1 Argument 1.
        /// \param a2 Argument 2.
        /// \param a3 Argument 3.
        /// \param a4 Argument 4.
        /// \param a5 Argument 5.
        /// \param a6 Argument 6.
        /// \return Description needed.
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


        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \tparam a6T Type of argument 6.
        /// \tparam a7T Type of argument 7.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param a1 Argument 1.
        /// \param a2 Argument 2.
        /// \param a3 Argument 3.
        /// \param a4 Argument 4.
        /// \param a5 Argument 5.
        /// \param a6 Argument 6.
        /// \param a7 Argument 7.
        /// \return Description needed.
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


        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \tparam a6T Type of argument 6.
        /// \tparam a7T Type of argument 7.
        /// \tparam a8T Type of argument 8.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param a1 Argument 1.
        /// \param a2 Argument 2.
        /// \param a3 Argument 3.
        /// \param a4 Argument 4.
        /// \param a5 Argument 5.
        /// \param a6 Argument 6.
        /// \param a7 Argument 7.
        /// \param a8 Argument 8.
        /// \return Description needed.
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


        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \tparam a6T Type of argument 6.
        /// \tparam a7T Type of argument 7.
        /// \tparam a8T Type of argument 8.
        /// \tparam a9T Type of argument 9.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param a1 Argument 1.
        /// \param a2 Argument 2.
        /// \param a3 Argument 3.
        /// \param a4 Argument 4.
        /// \param a5 Argument 5.
        /// \param a6 Argument 6.
        /// \param a7 Argument 7.
        /// \param a8 Argument 8.
        /// \param a9 Argument 9.
        /// \return Description needed.
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

        /// Sends task to derived class method `returnT (this->*memfn)()`.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param attr Description needed.
        /// \return Description needed.
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

        /// Sends task to derived class method `returnT (this->*memfn)(a1)`.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam a1T Type of argument 1.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param a1 Argument 1.
        /// \param attr Description needed.
        /// \return Description needed.
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

        /// Sends task to derived class method `returnT (this->*memfn)(a1,a2)`.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param a1 Argument 1.
        /// \param a2 Argument 2.
        /// \param attr Description needed.
        /// \return Description needed.
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

        /// Sends task to derived class method `returnT (this->*memfn)(a1,a2,a3)`.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param a1 Argument 1.
        /// \param a2 Argument 2.
        /// \param a3 Argument 3.
        /// \param attr Description needed.
        /// \return Description needed.
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

        /// Sends task to derived class method `returnT (this->*memfn)(a1,a2,a3,a4)`.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param a1 Argument 1.
        /// \param a2 Argument 2.
        /// \param a3 Argument 3.
        /// \param a4 Argument 4.
        /// \param attr Description needed.
        /// \return Description needed.
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

        /// Sends task to derived class method `returnT (this->*memfn)(a1,a2,a3,a4,a5)`.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param a1 Argument 1.
        /// \param a2 Argument 2.
        /// \param a3 Argument 3.
        /// \param a4 Argument 4.
        /// \param a5 Argument 5.
        /// \param attr Description needed.
        /// \return Description needed.
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

        /// Sends task to derived class method `returnT (this->*memfn)(a1,a2,a3,a4,a5,a6)`.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \tparam a6T Type of argument 6.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param a1 Argument 1.
        /// \param a2 Argument 2.
        /// \param a3 Argument 3.
        /// \param a4 Argument 4.
        /// \param a5 Argument 5.
        /// \param a6 Argument 6.
        /// \param attr Description needed.
        /// \return Description needed.
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

        /// Sends task to derived class method `returnT (this->*memfn)(a1,a2,a3,a4,a5,a6,a7)`.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \tparam a6T Type of argument 6.
        /// \tparam a7T Type of argument 7.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param a1 Argument 1.
        /// \param a2 Argument 2.
        /// \param a3 Argument 3.
        /// \param a4 Argument 4.
        /// \param a5 Argument 5.
        /// \param a6 Argument 6.
        /// \param a7 Argument 7.
        /// \param attr Description needed.
        /// \return Description needed.
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

        /// Sends task to derived class method `returnT (this->*memfn)(a1,a2,a3,a4,a5,a6,a7,a8)`.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \tparam a6T Type of argument 6.
        /// \tparam a7T Type of argument 7.
        /// \tparam a8T Type of argument 8.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param a1 Argument 1.
        /// \param a2 Argument 2.
        /// \param a3 Argument 3.
        /// \param a4 Argument 4.
        /// \param a5 Argument 5.
        /// \param a6 Argument 6.
        /// \param a7 Argument 7.
        /// \param a8 Argument 8.
        /// \param attr Description needed.
        /// \return Description needed.
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

        /// Sends task to derived class method `returnT (this->*memfn)(a1,a2,a3,a4,a5,a6,a7,a8,a9)`.

        /// \todo Descriptions needed.
        /// \tparam memfnT Verify: Signature of the member function in the derived class to be invoked for the task.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \tparam a6T Type of argument 6.
        /// \tparam a7T Type of argument 7.
        /// \tparam a8T Type of argument 8.
        /// \tparam a9T Type of argument 9.
        /// \param dest Description needed.
        /// \param memfn Verify: The member function to be invoked for the task.
        /// \param a1 Argument 1.
        /// \param a2 Argument 2.
        /// \param a3 Argument 3.
        /// \param a4 Argument 4.
        /// \param a5 Argument 5.
        /// \param a6 Argument 6.
        /// \param a7 Argument 7.
        /// \param a8 Argument 8.
        /// \param a9 Argument 9.
        /// \param attr Description needed.
        /// \return Description needed.
        ///
        /// \todo Could we use variadic templates to eliminate a lot of this code duplication?
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
            if(initialized()) {
              MADNESS_ASSERT_NOEXCEPT(World::exists(&world) && world.ptr_from_id<Derived>(id()) && *world.ptr_from_id<Derived>(id()) == this);
              world.unregister_ptr(this->id());
            }
        }

#ifdef MADNESS_WORLDOBJECT_FUTURE_TRACE
        /// "traces" future evaluation by counting their assignments in a static table

        /// Counts future assignments in a a statically-sized table to make this as lightweight/lock-free as possible
        /// with minimal effort. Can only trace objects of a single World.
        /// \param[in,out] f the future to be traced; if ready, will be unchanged (but contribute to the trace
        /// statistics of this object), or have a callback registered that will update the tracing statistics on
        /// assignment
        /// \warning this function will trace futures for WorldObjects associated with default world (id=0) only;
        /// use CMake variable `MADNESS_WORLDOBJECT_FUTURE_TRACE_WORLD_ID` to adjust the target World ID.
        /// \warning this function will trace futures for WorldObjects with IDs < 1000000 only;
        /// use CMake variable `MADNESS_WORLDOBJECT_FUTURE_TRACE_MAX_NOBJECTS` to adjust the limit.
        template <typename T>
        std::enable_if_t<!std::is_same_v<T,void>,void> trace(Future<T>& f) const;

        /// \param[in] id a WorldObject ID
        /// \return true if futures associated with \p id are traced
        static bool trace_futures(const uniqueidT &id);

        /// \return true if futures associated with this object are traced
        bool trace_futures() const {
          return trace_futures(this->id());
        }

        /// \param[in] id report tracing stats for this WorldObject
        /// \return number of futures given to trace() of the WorldObject with ID \p id
        static std::size_t trace_status_nfuture_registered(const uniqueidT& id);

        /// \return number of futures given to trace() of this object
        std::size_t trace_status_nfuture_registered() const {
          return trace_status_nfuture_registered(this->id());
        }

        /// \param[in] id report tracing stats for this WorldObject
        /// \return number of assigned futures given to trace() of the WorldObject with ID \p id
        static std::size_t trace_status_nfuture_assigned(const uniqueidT& id);

        /// \return number of assigned futures registered via `this->trace()`
        std::size_t trace_status_nfuture_assigned() const {
          return trace_status_nfuture_assigned(this->id());
        }
#endif
    };

#ifdef MADNESS_WORLDOBJECT_FUTURE_TRACE
    namespace detail {
    template <typename Derived> struct WorldObjectFutureTracer {
      // this value is the world ID to trace
      constexpr static std::size_t world_id =
#ifndef MADNESS_WORLDOBJECT_FUTURE_TRACE_WORLD_ID
          0
#else
          MADNESS_WORLDOBJECT_FUTURE_TRACE_WORLD_ID
#endif
          ;
      // this value is 1 greater than is the highest ID of WorldObjects to trace
      constexpr static std::size_t max_object_id =
#ifndef MADNESS_WORLDOBJECT_FUTURE_TRACE_MAX_NOBJECTS
          1000000
#else
          MADNESS_WORLDOBJECT_FUTURE_TRACE_MAX_NOBJECTS
#endif
          ;

      static constexpr bool do_trace(const uniqueidT& id) {
        return id.get_world_id() == world_id &&
               id.get_obj_id() < max_object_id;
      }

      static std::array<std::atomic<std::size_t>, max_object_id>
          nfuture_registered;
      static std::array<std::atomic<std::size_t>, max_object_id>
          nfuture_assigned;

      struct Initializer {
        Initializer() {
          for (auto &&v : nfuture_registered) {
            v.store(0);
          }
          for (auto &&v : nfuture_assigned) {
            v.store(0);
          }
        }
      };
      static Initializer initializer;

      struct FutureTracer : public CallbackInterface {
        FutureTracer(const uniqueidT &id) : id_(id) {
          if (do_trace(id_)) {
            nfuture_registered[id_.get_obj_id()]++;
          }
        }

        // Not allowed
        FutureTracer(const FutureTracer &) = delete;
        FutureTracer &operator=(const FutureTracer &) = delete;

        virtual ~FutureTracer() {}

        /// Notify this object that the future has been set.

        /// This will set the value of the future on the remote node and delete
        /// this callback object.
        void notify() override {
          if (do_trace(id_)) {
            nfuture_assigned[id_.get_obj_id()]++;
          }
          delete this;
        }

      private:
        uniqueidT id_;
      }; // struct FutureTracer

    }; // struct WorldObjectFutureTracer
    template <typename Derived>
    typename WorldObjectFutureTracer<Derived>::Initializer
        WorldObjectFutureTracer<Derived>::initializer;
    template <typename Derived>
    std::array<std::atomic<std::size_t>, WorldObjectFutureTracer<Derived>::max_object_id>
        WorldObjectFutureTracer<Derived>::nfuture_registered;
    template <typename Derived>
    std::array<std::atomic<std::size_t>, WorldObjectFutureTracer<Derived>::max_object_id>
        WorldObjectFutureTracer<Derived>::nfuture_assigned;
    }  // namespace detail

    template <typename Derived>
    template <typename T>
    std::enable_if_t<!std::is_same_v<T,void>,void> WorldObject<Derived>::trace(Future<T>& f) const {
      f.register_callback(
          new typename detail::WorldObjectFutureTracer<Derived>::FutureTracer(
              this->id()));
    }

    template <typename Derived>
    bool WorldObject<Derived>::trace_futures(const uniqueidT& id) {
      return detail::WorldObjectFutureTracer<Derived>::do_trace(id);
    }

    template <typename Derived>
    std::size_t WorldObject<Derived>::trace_status_nfuture_registered(const uniqueidT& id) {
      if (detail::WorldObjectFutureTracer<
              Derived>::do_trace(id)) {
        return detail::WorldObjectFutureTracer<
            Derived>::nfuture_registered[id.get_obj_id()];
      }
      else return 0;
    }

    template <typename Derived>
    std::size_t WorldObject<Derived>::trace_status_nfuture_assigned(const uniqueidT& id) {
      if (detail::WorldObjectFutureTracer<
              Derived>::do_trace(id)) {
        return detail::WorldObjectFutureTracer<
            Derived>::nfuture_assigned[id.get_obj_id()];
      }
      else return 0;
    }

#endif  // MADNESS_WORLDOBJECT_FUTURE_TRACE

    namespace archive {
        
        /// Specialization of \c ArchiveLoadImpl for globally-addressable objects.

        /// \tparam Derived The derived class of \c WorldObject in a curiously
        ///     repeating template pattern.
        template <class Derived>
        struct ArchiveLoadImpl<BufferInputArchive,WorldObject<Derived>*> {

            /// Read a globally-addressable object from a \c BufferInputArchive.

            /// \param[in,out] ar The archive.
            /// \param[out] ptr The read object.
            static inline void load(const BufferInputArchive& ar, WorldObject<Derived>*& ptr) {
                uniqueidT id;
                ar & id;
                World* world = World::world_from_id(id.get_world_id());
                MADNESS_ASSERT(world);
                auto ptr_opt = world->ptr_from_id< WorldObject<Derived> >(id);
                if (!ptr_opt) MADNESS_EXCEPTION("WorldObj: remote operation attempting to use a locally uninitialized object",0);
                ptr = *ptr_opt;
                if (!ptr) MADNESS_EXCEPTION("WorldObj: remote operation attempting to use a locally deregistered object",0);
            }
        };

        /// Specialization of \c ArchiveStoreImpl for globally-addressable objects.

        /// \tparam Derived The derived class of \c WorldObject in a curiously
        ///     repeating template pattern.
        template <class Derived>
        struct ArchiveStoreImpl<BufferOutputArchive,WorldObject<Derived>*> {

            /// Write a globally-addressable object to a \c BufferOutputArchive.

            /// \param[in,out] ar The archive.
            /// \param[in] ptr The object to store.
            static inline void store(const BufferOutputArchive& ar, WorldObject<Derived>* const& ptr) {
                ar & ptr->id();
            }
        };

        /// Specialization of \c ArchiveLoadImpl for constant, globally-addressable objects.

        /// \tparam Derived The derived class of \c WorldObject in a curiously
        ///     repeating template pattern.
        template <class Derived>
        struct ArchiveLoadImpl<BufferInputArchive, const WorldObject<Derived>*> {

            /// Read a globally-addressable object from a \c BufferInputArchive.

            /// \param[in,out] ar The archive.
            /// \param[out] ptr The read object.
            static inline void load(const BufferInputArchive& ar, const WorldObject<Derived>*& ptr) {
                uniqueidT id;
                ar & id;
                World* world = World::world_from_id(id.get_world_id());
                MADNESS_ASSERT(world);
                auto ptr_opt = world->ptr_from_id< WorldObject<Derived> >(id);
                if (!ptr_opt) MADNESS_EXCEPTION("WorldObj: remote operation attempting to use a locally uninitialized object",0);
                ptr = *ptr_opt;
                if (!ptr) MADNESS_EXCEPTION("WorldObj: remote operation attempting to use a locally deregistered object",0);
            }
        };

        /// Specialization of \c ArchiveStoreImpl for constant, globally-addressable objects.

        /// \tparam Derived The derived class of \c WorldObject in a curiously
        ///     repeating template pattern.
        template <class Derived>
        struct ArchiveStoreImpl<BufferOutputArchive, const WorldObject<Derived>*> {

            /// Write a globally-addressable object to a \c BufferOutputArchive.

            /// \param[in,out] ar The archive.
            /// \param[in] ptr The object to store.
            static inline void store(const BufferOutputArchive& ar, const WorldObject<Derived>* const& ptr) {
                ar & ptr->id();
            }
        };
    }   // namespace archive
}  // namespace madness

#endif // MADNESS_WORLD_WORLD_OBJECT_H__INCLUDED

/// @}
