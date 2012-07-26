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


#ifndef MADNESS_WORLD_WORLDTASK_H__INCLUDED
#define MADNESS_WORLD_WORLDTASK_H__INCLUDED

/// \file worldtask.h
/// \brief Defines TaskInterface and implements WorldTaskQueue and associated stuff.

// To do:
// a) Redo this ontop of serializable tasks which will remote much of the clutter
//    due to multiple length argument lists.
// b) Stealing which pretty much presume a) has been done

#include <iostream>
#include <world/nodefaults.h>
//#include <world/worldtypes.h>
//#include <world/typestuff.h>
//#include <world/worlddep.h>
//#include <world/worldthread.h>
#include <world/worldrange.h>
//#include <world/worldfut.h>
#include <world/worldtime.h>
#include <world/taskfn.h>

#if HAVE_INTEL_TBB
#define MADNESS_ALLOC_CHILD_TASK ( ThreadPool::tbb_parent_task->allocate_child() )
#else
#define MADNESS_ALLOC_CHILD_TASK
#endif

namespace madness {

    // Forward decls
    class World;
    class WorldTaskQueue;
    template <typename functionT> struct TaskFunction;
    template <typename memfunT> struct TaskMemfun;

    namespace detail {

        /// Serialization container for sending tasks to remote nodes

        /// This is for internal use only. You should not use this class directly.
        /// \tparam refT The remote reference type for task result future
        /// \tparam functionT The task function type
        template <typename refT, typename functionT>
        struct TaskHandlerInfo {
            refT ref;               ///< Remote reference for a task result future
            functionT func;         ///< A task function
            TaskAttributes attr;    ///< Task attributes

            /// Construct task info object

            /// \param ref Remote reference to the result future
            /// \param func The task function
            /// \param attr The task attrubutes
            TaskHandlerInfo(const refT& ref, functionT func, const TaskAttributes& attr)
                    : ref(ref), func(func),attr(attr) {}
            TaskHandlerInfo() {}

            /// Serialization of object

            /// \tparam Archive The serialization archive type
            /// \param ar The serialization archive
            template <typename Archive>
            void serialize(const Archive& ar) {
                serialize_internal<functionT>(ar);
            }

        private:

            /// Identify function and member function pointers

            /// \tparam T The function to identify
            template <typename fnT>
            struct is_func_ptr {
                static const bool value =
                    (std::is_function<typename std::remove_pointer<fnT>::type >::value
                    || std::is_member_function_pointer<fnT>::value);
            };

            /// Serialization for function pointers and member function pointers

            /// \tparam fnT The function type
            /// \tparam Archive The serialization archive type
            /// \param ar The serialization archive
            template <typename fnT, typename Archive>
            typename enable_if<is_func_ptr<fnT> >::type
            serialize_internal(const Archive& ar) {
                ar & ref & archive::wrap_opaque(func) & attr;
            }

            /// Serialization for non- function pointers and member function pointers.

            /// \tparam fnT The function type
            /// \tparam Archive The serialization archive type
            /// \param ar The serialization archive
            template <typename fnT, typename Archive>
            typename disable_if<is_func_ptr<fnT> >::type
            serialize_internal(const Archive& ar) {
                ar & ref & func & attr;
            }
        }; // struct TaskHandlerInfo

        template <typename fnT>
        struct function_enabler : public
            lazy_enable_if_c<
                function_traits<fnT>::value || has_result_type<fnT>::value,
                task_result_type<fnT> >
        { };

        template <typename memfnT>
        struct memfunc_enabler : public
            lazy_enable_if<
                std::is_member_function_pointer<memfnT>,
                task_result_type<memfnT> >
        { };

        template <typename T>
        struct DefaultInitPtr {
            static T init() { return NULL; }
        };

        template <typename T>
        struct DefaultInitPtr<std::shared_ptr<T> > {
            static std::shared_ptr<T> init() { return std::shared_ptr<T>(); }
        };

        template <typename ptrT, typename memfnT, typename resT>
        class MemFuncWrapper {
        private:
            ptrT ptr_;
            memfnT memfn_;

        public:
            typedef MemFuncWrapper<ptrT, memfnT, void> MemFuncWrapper_;
            typedef resT result_type;

            MemFuncWrapper() : ptr_(DefaultInitPtr<ptrT>::init()), memfn_() { }

            MemFuncWrapper(const MemFuncWrapper_& other) :
                ptr_(other.ptr_), memfn_(other.memfn_)
            { }

            MemFuncWrapper(ptrT ptr, memfnT memfn) :
                ptr_(ptr), memfn_(memfn)
            { }

            MemFuncWrapper_& operator=(const MemFuncWrapper_& other) {
                ptr_ = other.ptr_;
                memfn_ = other.memfn_;
                return *this;
            }

            resT operator()() const {
                MADNESS_ASSERT(ptr_);
                return (ptr_->*memfn_)();
            }

            template <typename a1T>
            resT operator()(a1T& a1) const {
                MADNESS_ASSERT(ptr_);
                return (ptr_->*memfn_)(a1);
            }

            template <typename a1T, typename a2T>
            resT operator()(a1T& a1, a2T& a2) const {
                MADNESS_ASSERT(ptr_);
                return (ptr_->*memfn_)(a1, a2);
            }

            template <typename a1T, typename a2T, typename a3T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3) const {
                MADNESS_ASSERT(ptr_);
                return (ptr_->*memfn_)(a1, a2, a3);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4) const {
                MADNESS_ASSERT(ptr_);
                return (ptr_->*memfn_)(a1, a2, a3, a4);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5) const {
                MADNESS_ASSERT(ptr_);
                return (ptr_->*memfn_)(a1, a2, a3, a4, a5);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6) const {
                MADNESS_ASSERT(ptr_);
                return (ptr_->*memfn_)(a1, a2, a3, a4, a5, a6);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7) const {
                MADNESS_ASSERT(ptr_);
                return (ptr_->*memfn_)(a1, a2, a3, a4, a5, a6, a7);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T, typename a8T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7, a8T& a8) const {
                MADNESS_ASSERT(ptr_);
                return (ptr_->*memfn_)(a1, a2, a3, a4, a5, a6, a7, a8);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
            resT operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7, a8T& a8, a9T& a9) const {
                MADNESS_ASSERT(ptr_);
                return (ptr_->*memfn_)(a1, a2, a3, a4, a5, a6, a7, a8, a9);
            }

            template <typename Archive>
            void serialize(const Archive& ar) {
                ar & ptr_ & memfn_;
            }
        };

        template <typename ptrT, typename memfnT>
        class MemFuncWrapper<ptrT, memfnT, void> {
        private:
            ptrT ptr_;
            memfnT memfn_;

        public:
            typedef MemFuncWrapper<ptrT, memfnT, void> MemFuncWrapper_;
            typedef void result_type;

            MemFuncWrapper() : ptr_(DefaultInitPtr<ptrT>::init()), memfn_() { }

            MemFuncWrapper(const MemFuncWrapper_& other) :
                ptr_(other.ptr_), memfn_(other.memfn_)
            { }

            MemFuncWrapper(ptrT ptr, memfnT memfn) :
                ptr_(ptr), memfn_(memfn)
            { }

            MemFuncWrapper_& operator=(const MemFuncWrapper_& other) {
                ptr_ = other.ptr_;
                memfn_ = other.memfn_;
                return *this;
            }

            void operator()() const {
                MADNESS_ASSERT(ptr_);
                (ptr_->*memfn_)();
            }

            template <typename a1T>
            void operator()(a1T& a1) const {
                MADNESS_ASSERT(ptr_);
                (ptr_->*memfn_)(a1);
            }

            template <typename a1T, typename a2T>
            void operator()(a1T& a1, a2T& a2) const {
                MADNESS_ASSERT(ptr_);
                (ptr_->*memfn_)(a1, a2);
            }

            template <typename a1T, typename a2T, typename a3T>
            void operator()(a1T& a1, a2T& a2, a3T& a3) const {
                MADNESS_ASSERT(ptr_);
                (ptr_->*memfn_)(a1, a2, a3);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T>
            void operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4) const {
                MADNESS_ASSERT(ptr_);
                (ptr_->*memfn_)(a1, a2, a3, a4);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T>
            void operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5) const {
                MADNESS_ASSERT(ptr_);
                (ptr_->*memfn_)(a1, a2, a3, a4, a5);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T>
            void operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6) const {
                MADNESS_ASSERT(ptr_);
                (ptr_->*memfn_)(a1, a2, a3, a4, a5, a6);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T>
            void operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7) const {
                MADNESS_ASSERT(ptr_);
                (ptr_->*memfn_)(a1, a2, a3, a4, a5, a6, a7);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T, typename a8T>
            void operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7, a8T& a8) const {
                MADNESS_ASSERT(ptr_);
                (ptr_->*memfn_)(a1, a2, a3, a4, a5, a6, a7, a8);
            }

            template <typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
            void operator()(a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7, a8T& a8, a9T& a9) const {
                MADNESS_ASSERT(ptr_);
                (ptr_->*memfn_)(a1, a2, a3, a4, a5, a6, a7, a8, a9);
            }

            template <typename Archive>
            void serialize(const Archive& ar) {
                ar & ptr_ & memfn_;
            }
        };

        template <typename objT, typename memfnT>
        MemFuncWrapper<objT*, memfnT, typename result_of<memfnT>::type>
        wrap_mem_fn(objT& obj, memfnT memfn) {
            return MemFuncWrapper<objT*, memfnT,
                    typename memfunc_traits<memfnT>::result_type>(& obj, memfn);
        }

        template <typename objT, typename memfnT>
        MemFuncWrapper<objT*, memfnT, typename result_of<memfnT>::type>
        wrap_mem_fn(objT* obj, memfnT memfn) {
            return MemFuncWrapper<objT*, memfnT,
                    typename memfunc_traits<memfnT>::result_type>(obj, memfn);
        }

        template <typename objT, typename memfnT>
        MemFuncWrapper<std::shared_ptr<objT>, memfnT, typename result_of<memfnT>::type>
        wrap_mem_fn(const std::shared_ptr<objT>& obj, memfnT memfn) {
            return MemFuncWrapper<std::shared_ptr<objT>, memfnT,
                    typename memfunc_traits<memfnT>::result_type>(obj, memfn);
        }

    }  // namespace detail




    /// Multi-threaded queue to manage and run tasks.
    class WorldTaskQueue : public CallbackInterface, private NO_DEFAULTS {
        friend class TaskInterface;
    private:
        World& world;              ///< The communication context
        const ProcessID me;        ///< This process
        AtomicInt nregistered;     ///< Counts pending tasks

        void notify() { nregistered--; }

        // Used in for_each kernel to check completion
        static bool completion_status(bool left, bool right) {
            return (left && right);
        }


        // Used in reduce kernel
        template <typename resultT, typename opT>
        static resultT sum(const resultT& left, const resultT& right, const opT& op) {
            //std::cout << " REDUCE SUM " << left << " " << right << std::endl;
            return op(left,right);
        }

        template <typename taskT, typename a1T, typename a2T, typename a3T, typename a4T, typename a5T,
                typename a6T, typename a7T, typename a8T, typename a9T>
        static void spawn_remote_task_handler(const AmArg& arg) {
            MADNESS_ASSERT(taskT::arity <= 9u);

            // Get task info and arguments form active message

            detail::TaskHandlerInfo<typename taskT::futureT::remote_refT,
                    typename taskT::functionT> info;
            a1T a1;
            a2T a2;
            a3T a3;
            a4T a4;
            a5T a5;
            a6T a6;
            a7T a7;
            a8T a8;
            a9T a9;

            arg & info & a1 & a2 & a3 & a4 & a5 & a6 & a7 & a8 & a9;

            // Create result future
            typename taskT::futureT result(info.ref);

            // Construct task
            taskT* task = NULL;
            switch(taskT::arity) {
            case 0u:
                task = new MADNESS_ALLOC_CHILD_TASK taskT(result, info.func, info.attr);
                break;
            case 1u:
                task = new MADNESS_ALLOC_CHILD_TASK taskT(result, info.func, a1, info.attr);
                break;
            case 2u:
                task = new MADNESS_ALLOC_CHILD_TASK taskT(result, info.func, a1, a2, info.attr);
                break;
            case 3u:
                task = new MADNESS_ALLOC_CHILD_TASK taskT(result, info.func, a1, a2, a3, info.attr);
                break;
            case 4u:
                task = new MADNESS_ALLOC_CHILD_TASK taskT(result, info.func, a1, a2, a3, a4, info.attr);
                break;
            case 5u:
                task = new MADNESS_ALLOC_CHILD_TASK taskT(result, info.func, a1, a2, a3, a4, a5, info.attr);
                break;
            case 6u:
                task = new MADNESS_ALLOC_CHILD_TASK taskT(result, info.func, a1, a2, a3, a4, a5, a6, info.attr);
                break;
            case 7u:
                task = new MADNESS_ALLOC_CHILD_TASK taskT(result, info.func, a1, a2, a3, a4, a5, a6, a7, info.attr);
                break;
            case 8u:
                task = new MADNESS_ALLOC_CHILD_TASK taskT(result, info.func, a1, a2, a3, a4, a5, a6, a7, a8, info.attr);
                break;
            case 9u:
                task = new MADNESS_ALLOC_CHILD_TASK taskT(result, info.func, a1, a2, a3, a4, a5, a6, a7, a8, a9, info.attr);
                break;
            }

            // Add task to queue
            arg.get_world()->taskq.add(task);
        }

        template <typename T>
        inline const T& am_arg(const Future<T>& f) {
            MADNESS_ASSERT(f.probe());
            return f.get();
        }

        template <typename T> inline const T& am_arg(const T& t) { return t; }


        typedef Future<void> voidT;


        template <typename taskT, typename fnT, typename a1T, typename a2T, typename a3T,
                typename a4T, typename a5T, typename a6T, typename a7T,
                typename a8T, typename a9T>
        inline typename taskT::futureT
        task_sender(ProcessID where, fnT fn, const a1T& a1,
                const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5,
                const a6T& a6, const a7T& a7, const a8T& a8, const a9T& a9,
                const TaskAttributes& attr)
        {
            typename taskT::futureT result;
            typedef detail::TaskHandlerInfo<typename taskT::futureT::remote_refT, typename taskT::functionT> infoT;
            world.am.send(where, & WorldTaskQueue::template spawn_remote_task_handler<taskT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T>,
                    new_am_arg(infoT(result.remote_ref(world), fn, attr),
                    a1, a2, a3, a4, a5, a6, a7, a8, a9));

            return result;
        }

        template <typename fnT>
        typename detail::task_result_type<fnT>::futureT
        send_task(ProcessID dest, fnT fn,
                const TaskAttributes& attr)
        {
            typedef TaskFn<fnT> taskT;
            return task_sender<taskT>(dest, fn, voidT(), voidT(), voidT(), voidT(),
                    voidT(), voidT(), voidT(), voidT(), voidT(), attr);
        }

        template <typename fnT, typename a1T>
        typename detail::task_result_type<fnT>::futureT
        send_task(ProcessID dest, fnT fn, const a1T& a1,
                const TaskAttributes& attr)
        {
            typedef TaskFn<fnT, a1T> taskT;
            return task_sender<taskT>(dest, fn, am_arg(a1), voidT(), voidT(),
                    voidT(), voidT(), voidT(), voidT(), voidT(), voidT(), attr);
        }

        template <typename fnT, typename a1T, typename a2T>
        typename detail::task_result_type<fnT>::futureT
        send_task( ProcessID dest, fnT fn,  const a1T& a1, const a2T& a2,
                const TaskAttributes& attr)
        {
            typedef TaskFn<fnT, a1T, a2T> taskT;
            return task_sender<taskT>(dest, fn, am_arg(a1), am_arg(a2), voidT(),
                    voidT(), voidT(), voidT(), voidT(), voidT(), voidT(), attr);
        }

        template <typename fnT, typename a1T, typename a2T, typename a3T>
        typename detail::task_result_type<fnT>::futureT
        send_task(ProcessID dest, fnT fn,  const a1T& a1, const a2T& a2,
                const a3T& a3, const TaskAttributes& attr)
        {
            typedef TaskFn<fnT, a1T, a2T, a3T> taskT;
            return task_sender<taskT>(dest, fn, am_arg(a1), am_arg(a2), am_arg(a3),
                    voidT(), voidT(), voidT(), voidT(), voidT(), voidT(), attr);
        }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T>
        typename detail::task_result_type<fnT>::futureT
        send_task(ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4,
                const TaskAttributes& attr)
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T> taskT;
            return task_sender<taskT>(dest, fn, am_arg(a1), am_arg(a2), am_arg(a3),
                    am_arg(a4), voidT(), voidT(), voidT(), voidT(), voidT(), attr);
        }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T>
        typename detail::task_result_type<fnT>::futureT
        send_task(ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5,
                const TaskAttributes& attr)
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T> taskT;
            return task_sender<taskT>(dest, fn, am_arg(a1), am_arg(a2), am_arg(a3),
                    am_arg(a4), am_arg(a5), voidT(), voidT(), voidT(), voidT(),
                    attr);
        }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T>
        typename detail::task_result_type<fnT>::futureT
        send_task(ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6,
                const TaskAttributes& attr)
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T> taskT;
            return task_sender<taskT>(dest, fn, am_arg(a1), am_arg(a2), am_arg(a3),
                    am_arg(a4), am_arg(a5), am_arg(a6), voidT(), voidT(), voidT(),
                    attr);
        }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T>
        typename detail::task_result_type<fnT>::futureT
        send_task(ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7,
                const TaskAttributes& attr) {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T> taskT;
            return task_sender<taskT>(dest, fn, am_arg(a1), am_arg(a2), am_arg(a3),
                    am_arg(a4), am_arg(a5), am_arg(a6), am_arg(a7), voidT(), voidT(),
                    attr);
        }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T, typename a8T>
        typename detail::task_result_type<fnT>::futureT
        send_task(ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7,
                const a8T& a8, const TaskAttributes& attr)
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T> taskT;
            return task_sender<taskT>(dest, fn, am_arg(a1), am_arg(a2), am_arg(a3),
                    am_arg(a4), am_arg(a5), am_arg(a6), am_arg(a7), am_arg(a8),
                    voidT(), attr);
        }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
        typename detail::task_result_type<fnT>::futureT
        send_task(ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7,
                const a8T& a8, const a9T& a9, const TaskAttributes& attr)
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T> taskT;
            return task_sender<taskT>(dest, fn, am_arg(a1), am_arg(a2), am_arg(a3),
                    am_arg(a4), am_arg(a5), am_arg(a6), am_arg(a7), am_arg(a8),
                    am_arg(a9), attr);
        }


    public:
        WorldTaskQueue(World& world);

        /// Returns the number of pending tasks
        size_t size() const { return nregistered; }

        class Stealer {
            WorldTaskQueue& q;
            std::vector<TaskInterface*>& v;
            const int nsteal;

        public:
            Stealer(WorldTaskQueue& q, std::vector<TaskInterface*>& v, int nsteal)
                : q(q)
                , v(v)
                , nsteal(nsteal)
            {}

            bool operator()(PoolTaskInterface** pt);
        };


        /// Invoke locally or remotely to send tasks to process P
        std::vector<TaskInterface*> steal(int nsteal);

        /// Add a new local task taking ownership of the pointer

        /// The task pointer (t) is assumed to have been created with
        /// \c new and when the task is eventually run the queue
        /// will call the task's destructor using \c delete.
        ///
        /// Once the task is complete it will execute
        /// task_complete_callback to decrement the number of pending
        /// tasks and be deleted.
        void add(TaskInterface* t)  {
            nregistered++;

            t->set_info(&world, this);       // Stuff info

            if (t->ndep() == 0) {
                ThreadPool::add(t); // If no dependencies directly submit
            } else {
                // With dependencies must use the callback to avoid race condition
                t->register_submit_callback();
                //t->dec();
            }
        }

        template <typename fnT, typename a1T, typename a2T, typename a3T,
            typename a4T, typename a5T, typename a6T, typename a7T, typename a8T,
            typename a9T>
        typename TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T>::futureT
        add(TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T>* t) {
            typename TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T>::futureT
                res(t->result());
            add(static_cast<TaskInterface*>(t));
            return res;
        }

        /// Reduce op(item) for all items in range using op(sum,op(item))

        /// The operation must provide the following interface of
        /// which the \c operator() methods are required by reduce()
        /// and the rest by the task interface.
        /// \code
        /// struct opT {
        ///     opT();
        ///     opT(const &opT);
        ///     resultT operator()(const rangeT::iterator& it) const;
        ///     resultT operator()(const resultT& left, const resultT& right);
        ///     template <typename Archive> void serialize(const Archive& ar);
        /// }
        /// \endcode
        /// Note that the serialize method does not actually have to
        /// work unless you want to have the task be stealable.
        /// Adjust the chunksize in the range to control granularity.
        template <typename resultT, typename rangeT, typename opT>
        Future<resultT> reduce(const rangeT& range, const opT& op) {
            if (range.size() <= range.get_chunksize()) {
                resultT sum = resultT();
                for (typename rangeT::iterator it=range.begin(); it != range.end(); ++it) sum = op(sum,op(it));
                return Future<resultT>(sum);
            } else {
                rangeT left = range;
                rangeT right(left,Split());

                Future<resultT>  leftsum = add(*this, &WorldTaskQueue::reduce<resultT,rangeT,opT>, left,  op);
                Future<resultT> rightsum = add(*this, &WorldTaskQueue::reduce<resultT,rangeT,opT>, right, op);
                return add(&WorldTaskQueue::sum<resultT,opT>, leftsum, rightsum, op);
            }
        }

        /// Apply op(item) for all items in range

        /// The operation must provide the following interface of
        /// which the \c operator() method is required by for_each()
        /// and the rest by the task interface.
        /// \code
        /// struct opT {
        ///     opT();
        ///     opT(const opT&);
        ///     bool operator()(const rangeT::iterator& it) const;
        ///     template <typename Archive> void serialize(const Archive& ar);
        /// };
        /// \endcode
        /// Note that the serialize method does not actually have to
        /// work unless you want to have the task be stealable.
        ///
        /// Adjust the chunksize in the range to control granularity.
        ///
        /// Your operation should return true/false for success failure
        /// and the logical and of all results is returned as the
        /// future result.
        ///
        /// You can ignore the result if you are interested
        /// in neither synchronization nor result status.
        template <typename rangeT, typename opT>
        Future<bool> for_each(const rangeT& range, const opT& op) {
            if (range.size() <= range.get_chunksize()) {
                bool status = true;
                for (typename rangeT::iterator it=range.begin();  it != range.end();  ++it) status &= op(it);
                return Future<bool>(status);
            } else {
                rangeT left = range;
                rangeT right(left,Split());
                Future<bool>  leftsum = add(*this, &WorldTaskQueue::for_each<rangeT,opT>, left,  op, TaskAttributes::hipri());
                Future<bool> rightsum = add(*this, &WorldTaskQueue::for_each<rangeT,opT>, right, op, TaskAttributes::hipri());
                return add(&WorldTaskQueue::completion_status, leftsum, rightsum);
            }
        }


        /// Spawn a local task

        /// Spawns a task on process \c where , which may or may not be this
        /// process.  An argument that is a future may be used to carry
        /// dependencies for local tasks.  An unready future cannot be used as
        /// an argument for a remote tasks -- i.e., remote  tasks must be ready
        /// to execute (you can work around this by making a local task to
        /// submit the remote task once everything is ready).
        /// \tparam functionT A function pointer or functor
        /// \param where The process where the task will be spawned
        /// \param function The function to be called in the task
        /// \param attr The task attributes
        /// \return A future to the task function's results. If the task function
        /// return type is \c void , a \c Future<void> object is return that may
        /// be ignored.
        /// \note Arguments must be (de)serializable and must of course make
        /// sense at the remote destination.  Fundamental types,
        /// simple STL containers, and pointers to World,
        /// WorldContainer, and user-defined types derived from
        /// WorldObject<> are automatically handled.  Anything else is
        /// your problem.
        template <typename fnT>
        typename detail::function_enabler<fnT>::type
        add(fnT fn, const TaskAttributes& attr=TaskAttributes()) {
            typedef TaskFn<fnT> taskT;
            return add(new MADNESS_ALLOC_CHILD_TASK taskT(typename taskT::futureT(),
                    fn, attr));
        }

        template <typename fnT, typename a1T>
        typename detail::function_enabler<fnT>::type
        add(fnT fn, const a1T& a1, const TaskAttributes& attr=TaskAttributes()) {
            typedef TaskFn<fnT, a1T> taskT;
            return add(new MADNESS_ALLOC_CHILD_TASK taskT(typename taskT::futureT(),
                    fn, a1, attr));
        }

        template <typename fnT, typename a1T, typename a2T>
        typename detail::function_enabler<fnT>::type
        add(fnT fn, const a1T& a1, const a2T& a2, const TaskAttributes& attr=TaskAttributes()) {
            typedef TaskFn<fnT, a1T, a2T> taskT;
            return add(new MADNESS_ALLOC_CHILD_TASK taskT(typename taskT::futureT(),
                    fn, a1, a2, attr));
        }

        template <typename fnT, typename a1T, typename a2T, typename a3T>
        typename detail::function_enabler<fnT>::type
        add(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3,
            const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T> taskT;
            return add(new MADNESS_ALLOC_CHILD_TASK taskT(typename taskT::futureT(),
                    fn, a1, a2, a3, attr));
        }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T>
        typename detail::function_enabler<fnT>::type
        add(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
            const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T> taskT;
            return add(new MADNESS_ALLOC_CHILD_TASK taskT(typename taskT::futureT(),
                    fn, a1, a2, a3, a4, attr));
        }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
            typename a5T>
        typename detail::function_enabler<fnT>::type
        add(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
            const a5T& a5, const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T> taskT;
            return add(new MADNESS_ALLOC_CHILD_TASK taskT(typename taskT::futureT(),
                    fn, a1, a2, a3, a4, a5, attr));
        }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
            typename a5T, typename a6T>
        typename detail::function_enabler<fnT>::type
        add(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
            const a5T& a5, const a6T& a6, const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T> taskT;
            return add(new MADNESS_ALLOC_CHILD_TASK taskT(typename taskT::futureT(),
                    fn, a1, a2, a3, a4, a5, a6, attr));
        }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
            typename a5T, typename a6T, typename a7T>
        typename detail::function_enabler<fnT>::type
        add(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
            const a5T& a5, const a6T& a6, const a7T& a7,
            const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T> taskT;
            return add(new MADNESS_ALLOC_CHILD_TASK taskT(typename taskT::futureT(),
                    fn, a1, a2, a3, a4, a5, a6, a7, attr));
        }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
            typename a5T, typename a6T, typename a7T, typename a8T>
        typename detail::function_enabler<fnT>::type
        add(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
            const a5T& a5, const a6T& a6, const a7T& a7, const a8T& a8,
            const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T> taskT;
            return add(new MADNESS_ALLOC_CHILD_TASK taskT(typename taskT::futureT(),
                    fn, a1, a2, a3, a4, a5, a6, a7, a8, attr));
        }


        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
            typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
        typename detail::function_enabler<fnT>::type
        add(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
            const a5T& a5, const a6T& a6, const a7T& a7, const a8T& a8, const a9T& a9,
            const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T> taskT;
            return add(new MADNESS_ALLOC_CHILD_TASK taskT(typename taskT::futureT(),
                    fn, a1, a2, a3, a4, a5, a6, a7, a8, a9, attr));
        }


        /// Spawn a remote task

        /// Spawns a task on process \c where , which may or may not be this
        /// process.
        /// \param where The process where the task will be spawned
        /// \param function The function to be called in the task
        /// \param attr The task attributes
        /// \return A future to the task function's results. If the task function
        /// return type is \c void , a \c Future<void> object is return that may
        /// be ignored.
        template <typename functionT>
        typename detail::function_enabler<functionT>::type
        add(ProcessID where, functionT function, const TaskAttributes attr=TaskAttributes()) {
            if(where == me)
                return add(function, attr);
            else
                return send_task(where, function, attr);
        }

        /// Spawn a remote task

        /// Spawns a task on process \c where , which may or may not be this
        /// process.  An argument that is a future may be used to carry
        /// dependencies for local tasks.  An unready future cannot be used as
        /// an argument for a remote tasks -- i.e., remote  tasks must be ready
        /// to execute (you can work around this by making a local task to
        /// submit the remote task once everything is ready).
        /// \param where The process where the task will be spawned
        /// \param function The function to be called in the task
        /// \param attr The task attributes
        /// \return A future to the task function's results. If the task function
        /// return type is \c void , a \c Future<void> object is return that may
        /// be ignored.
        /// Invoke "resultT (*function)(arg1T)" as a task, local or remote

        /// A future is returned to hold the eventual result of the
        /// task.  Future<void> is an empty class that may be ignored.
        /// \note Arguments must be (de)serializable and must of course make
        /// sense at the remote destination.  Fundamental types,
        /// simple STL containers, and pointers to World,
        /// WorldContainer, and user-defined types derived from
        /// WorldObject<> are automatically handled.  Anything else is
        /// your problem.
        template <typename functionT, typename arg1T>
        typename detail::function_enabler<functionT>::type
        add(ProcessID where, functionT function, const arg1T& arg1,
                const TaskAttributes attr=TaskAttributes()) {
            if(where == me)
                return add(function, arg1, attr);
            else
                return send_task(where, function, arg1, attr);
        }

        /// Invoke "resultT (*function)(arg1T,arg2T)" as a task, local or remote

        /// A future is returned to hold the eventual result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T>
        typename detail::function_enabler<functionT>::type
        add(ProcessID where, functionT function,
            const arg1T& arg1, const arg2T& arg2, const TaskAttributes attr=TaskAttributes()) {
            if (where == me)
                add(function, arg1, arg2, attr);
            else
                return send_task(where, function, arg1, arg2, attr);
        }

        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T)" as a task, local or remote

        /// A future is returned to hold the eventual result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T>
        typename detail::function_enabler<functionT>::type
        add(ProcessID where, functionT function,
            const arg1T& arg1, const arg2T& arg2,
            const arg3T& arg3, const TaskAttributes attr=TaskAttributes()) {
            if (where == me)
                return add(function, arg1, arg2, arg3, attr);
            else
                return send_task(where, function, arg1, arg2, arg3, attr);
        }

        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T,arg4T)" as a task, local or remote

        /// A future is returned to hold the eventual result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        typename detail::function_enabler<functionT>::type
        add(ProcessID where, functionT function,
            const arg1T& arg1, const arg2T& arg2,
            const arg3T& arg3, const arg4T& arg4, const TaskAttributes& attr=TaskAttributes()) {
            if (where == me)
                return add(function, arg1, arg2, arg3, arg4, attr);
            else
                return send_task(where, function, arg1, arg2, arg3, arg4, attr);
        }

        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T,arg4T,arg5T)" as a task, local or remote

        /// A future is returned to hold the eventual result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        typename detail::function_enabler<functionT>::type
        add(ProcessID where, functionT function,
            const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4,
            const arg5T& arg5, const TaskAttributes& attr=TaskAttributes()) {
            if(where == me)
                return add(function, arg1, arg2, arg3, arg4, arg5, attr);
            else
                return send_task(where, function, arg1, arg2, arg3, arg4, arg5, attr);
        }

        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T)" as a task, local or remote

        /// A future is returned to hold the eventual result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
        typename detail::function_enabler<functionT>::type
        add(ProcessID where, functionT function,
            const arg1T& arg1, const arg2T& arg2,
            const arg3T& arg3, const arg4T& arg4,
            const arg5T& arg5, const arg6T& arg6, const TaskAttributes& attr=TaskAttributes()) {
            if(where == me)
                return add(function, arg1, arg2, arg3, arg4, arg5, arg6, attr);
            else
                return send_task(where, function, arg1, arg2, arg3, arg4, arg5, arg6, attr);
        }

        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T)" as a task, local or remote

        /// A future is returned to hold the eventual result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
        typename detail::function_enabler<functionT>::type
        add(ProcessID where, functionT function,
            const arg1T& arg1, const arg2T& arg2,
            const arg3T& arg3, const arg4T& arg4,
            const arg5T& arg5, const arg6T& arg6, const arg7T& arg7, const TaskAttributes& attr=TaskAttributes()) {
            if(where == me)
                return add(function, arg1, arg2, arg3, arg4, arg5, arg6, arg7, attr);
            else
                return send_task(where, function, arg1, arg2, arg3, arg4, arg5, arg6, arg7, attr);
        }


        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T)" as a task, local or remote

        /// A future is returned to hold the eventual result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T, typename arg8T>
        typename detail::function_enabler<functionT>::type
        add(ProcessID where, functionT function,
            const arg1T& arg1, const arg2T& arg2,
            const arg3T& arg3, const arg4T& arg4,
            const arg5T& arg5, const arg6T& arg6,
        const arg7T& arg7, const arg8T& arg8, const TaskAttributes& attr=TaskAttributes()) {
            if(where == me)
                return add(function, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, attr);
            else
                send_task(where, function, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, attr);
        }

        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T,arg9T)" as a task, local or remote

        /// A future is returned to hold the eventual result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T, typename arg8T, typename arg9T>
        typename detail::function_enabler<functionT>::type
        add(ProcessID where, functionT function,
            const arg1T& arg1, const arg2T& arg2,
            const arg3T& arg3, const arg4T& arg4,
            const arg5T& arg5, const arg6T& arg6,
            const arg7T& arg7, const arg8T& arg8,
            const arg9T& arg9, const TaskAttributes& attr=TaskAttributes()) {
            if(where == me)
                return add(function, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, attr);
            else
                return send_task(where, function, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, attr);
        }


        /// Invoke "resultT (obj.*memfun)()" as a local task
        template <typename objT, typename memfunT>
        typename detail::memfunc_enabler<memfunT>::type
        add(objT& obj, memfunT memfun, const TaskAttributes& attr=TaskAttributes()) {
            return add(detail::wrap_mem_fn(obj,memfun),attr);;
        }


        /// Invoke "resultT (obj.*memfun)(arg1T)" as a local task
        template <typename objT, typename memfunT, typename arg1T>
        typename detail::memfunc_enabler<memfunT>::type
        add(objT& obj, memfunT memfun, const arg1T& arg1,
                const TaskAttributes& attr=TaskAttributes()) {
            return add(detail::wrap_mem_fn(obj,memfun),arg1,attr);
        }


        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T)" as a local task
        template <typename objT, typename memfunT, typename arg1T, typename arg2T>
        typename detail::memfunc_enabler<memfunT>::type
        add(objT& obj, memfunT memfun, const arg1T& arg1, const arg2T& arg2,
                const TaskAttributes& attr=TaskAttributes())
        { return add(detail::wrap_mem_fn(obj,memfun),arg1,arg2,attr); }


        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3)" as a local task
        template <typename objT, typename memfunT, typename arg1T, typename arg2T,
            typename arg3T>
        typename detail::memfunc_enabler<memfunT>::type
        add(objT& obj, memfunT memfun, const arg1T& arg1, const arg2T& arg2,
                const arg3T& arg3, const TaskAttributes& attr=TaskAttributes())
        { return add(detail::wrap_mem_fn(obj,memfun),arg1,arg2,arg3,attr); }


        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3,arg4)" as a local task
        template <typename objT, typename memfunT, typename arg1T, typename arg2T,
            typename arg3T, typename arg4T>
        typename detail::memfunc_enabler<memfunT>::type
        add(objT& obj, memfunT memfun, const arg1T& arg1, const arg2T& arg2,
                const arg3T& arg3, const arg4T& arg4, const TaskAttributes& attr=TaskAttributes())
        { return add(detail::wrap_mem_fn(obj,memfun),arg1,arg2,arg3,arg4,attr); }


        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3,arg4,arg5)" as a local task
        template <typename objT, typename memfunT, typename arg1T, typename arg2T,
            typename arg3T, typename arg4T, typename arg5T>
        typename detail::memfunc_enabler<memfunT>::type
        add(objT& obj, memfunT memfun, const arg1T& arg1, const arg2T& arg2,
                const arg3T& arg3, const arg4T& arg4, const arg5T& arg5,
                const TaskAttributes attr=TaskAttributes())
        { return add(detail::wrap_mem_fn(obj,memfun),arg1,arg2,arg3,arg4,arg5,attr); }


        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3,arg4,arg5,arg6)" as a local task
        template <typename objT, typename memfunT, typename arg1T, typename arg2T,
            typename arg3T, typename arg4T, typename arg5T, typename arg6T>
        typename detail::memfunc_enabler<memfunT>::type
        add(objT& obj, memfunT memfun, const arg1T& arg1, const arg2T& arg2,
                const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const arg6T& arg6,
                const TaskAttributes& attr=TaskAttributes())
        { return add(detail::wrap_mem_fn(obj,memfun),arg1,arg2,arg3,arg4,arg5,arg6,attr); }

        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3,arg4,arg5,arg6,arg7)" as a local task
        template <typename objT, typename memfunT, typename arg1T, typename arg2T,
            typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
        typename detail::memfunc_enabler<memfunT>::type
        add(objT& obj, memfunT memfun, const arg1T& arg1, const arg2T& arg2,
                const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const arg6T& arg6,
                const arg7T& arg7, const TaskAttributes& attr=TaskAttributes())
        { return add(detail::wrap_mem_fn(obj,memfun),arg1,arg2,arg3,arg4,arg5,arg6,arg7,attr); }

        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3,arg4,arg5,arg6,arg7,arg8)" as a local task
        template <typename objT, typename memfunT, typename arg1T, typename arg2T,
            typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T,
            typename arg8T>
        typename detail::memfunc_enabler<memfunT>::type
        add(objT& obj, memfunT memfun, const arg1T& arg1, const arg2T& arg2,
                const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const arg6T& arg6,
                const arg7T& arg7, const arg8T& arg8, const TaskAttributes attr=TaskAttributes())
        { return add(detail::wrap_mem_fn(obj,memfun),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,attr); }

        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3,arg4,arg5,arg6,arg7,arg8,arg9)" as a local task
        template <typename objT, typename memfunT, typename arg1T, typename arg2T,
            typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T,
            typename arg8T, typename arg9T>
        typename detail::memfunc_enabler<memfunT>::type
        add(objT& obj, memfunT memfun, const arg1T& arg1, const arg2T& arg2,
                const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const arg6T& arg6,
                const arg7T& arg7, const arg8T& arg8, const arg9T& arg9,
                const TaskAttributes attr=TaskAttributes())
        { return add(detail::wrap_mem_fn(obj,memfun),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,attr); }

        /// Invoke "resultT (obj.*memfun)()" as a local task
        template <typename objT, typename memfunT>
        typename detail::memfunc_enabler<memfunT>::type
        add(objT* obj, memfunT memfun, const TaskAttributes& attr=TaskAttributes()) {
            return add(detail::wrap_mem_fn(obj,memfun),attr);;
        }


        /// Invoke "resultT (obj.*memfun)(arg1T)" as a local task
        template <typename objT, typename memfunT, typename arg1T>
        typename detail::memfunc_enabler<memfunT>::type
        add(objT* obj, memfunT memfun, const arg1T& arg1,
                const TaskAttributes& attr=TaskAttributes()) {
            return add(detail::wrap_mem_fn(obj,memfun),arg1,attr);
        }


        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T)" as a local task
        template <typename objT, typename memfunT, typename arg1T, typename arg2T>
        typename detail::memfunc_enabler<memfunT>::type
        add(objT* obj, memfunT memfun, const arg1T& arg1, const arg2T& arg2,
                const TaskAttributes& attr=TaskAttributes())
        { return add(detail::wrap_mem_fn(obj,memfun),arg1,arg2,attr); }


        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3)" as a local task
        template <typename objT, typename memfunT, typename arg1T, typename arg2T,
            typename arg3T>
        typename detail::memfunc_enabler<memfunT>::type
        add(objT* obj, memfunT memfun, const arg1T& arg1, const arg2T& arg2,
                const arg3T& arg3, const TaskAttributes& attr=TaskAttributes())
        { return add(detail::wrap_mem_fn(obj,memfun),arg1,arg2,arg3,attr); }


        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3,arg4)" as a local task
        template <typename objT, typename memfunT, typename arg1T, typename arg2T,
            typename arg3T, typename arg4T>
        typename detail::memfunc_enabler<memfunT>::type
        add(objT* obj, memfunT memfun, const arg1T& arg1, const arg2T& arg2,
                const arg3T& arg3, const arg4T& arg4, const TaskAttributes& attr=TaskAttributes())
        { return add(detail::wrap_mem_fn(obj,memfun),arg1,arg2,arg3,arg4,attr); }


        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3,arg4,arg5)" as a local task
        template <typename objT, typename memfunT, typename arg1T, typename arg2T,
            typename arg3T, typename arg4T, typename arg5T>
        typename detail::memfunc_enabler<memfunT>::type
        add(objT* obj, memfunT memfun, const arg1T& arg1, const arg2T& arg2,
                const arg3T& arg3, const arg4T& arg4, const arg5T& arg5,
                const TaskAttributes attr=TaskAttributes())
        { return add(detail::wrap_mem_fn(obj,memfun),arg1,arg2,arg3,arg4,arg5,attr); }


        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3,arg4,arg5,arg6)" as a local task
        template <typename objT, typename memfunT, typename arg1T, typename arg2T,
            typename arg3T, typename arg4T, typename arg5T, typename arg6T>
        typename detail::memfunc_enabler<memfunT>::type
        add(objT* obj, memfunT memfun, const arg1T& arg1, const arg2T& arg2,
                const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const arg6T& arg6,
                const TaskAttributes& attr=TaskAttributes())
        { return add(detail::wrap_mem_fn(obj,memfun),arg1,arg2,arg3,arg4,arg5,arg6,attr); }

        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3,arg4,arg5,arg6,arg7)" as a local task
        template <typename objT, typename memfunT, typename arg1T, typename arg2T,
            typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
        typename detail::memfunc_enabler<memfunT>::type
        add(objT* obj, memfunT memfun, const arg1T& arg1, const arg2T& arg2,
                const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const arg6T& arg6,
                const arg7T& arg7, const TaskAttributes& attr=TaskAttributes())
        { return add(detail::wrap_mem_fn(obj,memfun),arg1,arg2,arg3,arg4,arg5,arg6,arg7,attr); }

        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3,arg4,arg5,arg6,arg7,arg8)" as a local task
        template <typename objT, typename memfunT, typename arg1T, typename arg2T,
            typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T,
            typename arg8T>
        typename detail::memfunc_enabler<memfunT>::type
        add(objT* obj, memfunT memfun, const arg1T& arg1, const arg2T& arg2,
                const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const arg6T& arg6,
                const arg7T& arg7, const arg8T& arg8, const TaskAttributes attr=TaskAttributes())
        { return add(detail::wrap_mem_fn(obj,memfun),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,attr); }

        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3,arg4,arg5,arg6,arg7,arg8,arg9)" as a local task
        template <typename objT, typename memfunT, typename arg1T, typename arg2T,
            typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T,
            typename arg8T, typename arg9T>
        typename detail::memfunc_enabler<memfunT>::type
        add(objT* obj, memfunT memfun, const arg1T& arg1, const arg2T& arg2,
                const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const arg6T& arg6,
                const arg7T& arg7, const arg8T& arg8, const arg9T& arg9,
                const TaskAttributes attr=TaskAttributes())
        { return add(detail::wrap_mem_fn(obj,memfun),arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,attr); }

        struct ProbeAllDone {
            WorldTaskQueue* tq;
            double start;
            ProbeAllDone(WorldTaskQueue* tq) : tq(tq),start(cpu_time()) {}
            bool operator()() const {
                if (cpu_time()-start > 1200) {
                    for (int loop = 0; loop<3; ++loop) {
                        std::cout << "HUNG Q? " << tq->size() << " " << ThreadPool::queue_size() << std::endl;
                        std::cout.flush();
                        myusleep(1000000);
                    }
                    MADNESS_ASSERT(cpu_time()-start < 1200);
                }
                return (tq->size() == 0);
            }
        };

        /// Returns after all local tasks have completed

        /// While waiting the calling thread will run tasks.
        void fence()  {
            ProbeAllDone tester(this);
            do {
                world.await(tester);
            }
            while (nregistered);
        }
    };

}


#endif // MADNESS_WORLD_WORLDTASK_H__INCLUDED
