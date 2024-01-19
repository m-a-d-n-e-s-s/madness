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
 \file world_task_queue.h
 \brief Defines \c TaskInterface and implements \c WorldTaskQueue and associated stuff.
 \ingroup taskq
*/

#ifndef MADNESS_WORLD_WORLD_TASK_QUEUE_H__INCLUDED
#define MADNESS_WORLD_WORLD_TASK_QUEUE_H__INCLUDED

#include <type_traits>
#include <iostream>

#include <madness/madness_config.h>

// must be included before world/range.h
#ifdef HAVE_INTEL_TBB
# include <tbb/parallel_reduce.h>
#endif

#include <madness/world/meta.h>
#include <madness/world/nodefaults.h>
#include <madness/world/range.h>
#include <madness/world/timers.h>
#include <madness/world/taskfn.h>
#include <madness/world/mem_func_wrapper.h>

/// \addtogroup taskq
/// @{

namespace madness {

    // Forward decls
    class World;
    class WorldTaskQueue;
    template <typename> struct TaskFunction;
    template <typename> struct TaskMemfun;

    namespace meta {
    template <typename ... argsT>
    struct taskattr_is_last_arg : public std::integral_constant<bool, std::is_same<std::decay_t<typename meta::last_type<argsT...>::type>,
    TaskAttributes>::value> {};
    template <>
    struct taskattr_is_last_arg<> : public std::false_type {};
    }

    namespace detail {

        // a few more forward decls
        template <typename ptrT, typename memfnT, typename resT>
        memfnT get_mem_func_ptr(const MemFuncWrapper<ptrT, memfnT, resT>&);

        template <typename, typename>
        class ForEachRootTask;

        /// Serialization container for sending tasks to remote nodes.

        /// \attention This struct is for internal use only. You should not
        ///     use this class directly.
        /// \tparam refT The remote reference type for task result future.
        /// \tparam functionT The task function type.
        template <typename refT, typename functionT>
        struct TaskHandlerInfo {
            refT ref; ///< Remote reference for a task result future.
            functionT func; ///< A task function.
            TaskAttributes attr; ///< Task attributes.

            TaskHandlerInfo() = default;

            /// Construct task info object.

            /// \param[in] ref Remote reference to the result future.
            /// \param[in] func The task function.
            /// \param[in] attr The task attrubutes.
            TaskHandlerInfo(const refT& ref, functionT func, const TaskAttributes& attr)
                : ref(ref), func(func),attr(attr) {}

            /// Serialization of an object.

            /// \tparam Archive The serialization archive type.
            /// \param[in,out] ar The serialization archive.
            template <typename Archive>
            void serialize(const Archive& ar) {
                serialize_internal<functionT>(ar);
            }

        private:

            /// Serialization for function pointers and member function pointers.

            /// \tparam fnT The function type.
            /// \tparam Archive The serialization archive type.
            /// \param[in,out] ar The serialization archive.
            template <typename fnT, typename Archive>
            typename std::enable_if<is_any_function_pointer_v<fnT>>::type
            serialize_internal(const Archive& ar) {
                ar & ref & func & attr;
            }

            /// Serialization for non- function pointers and member function pointers.

            /// \tparam fnT The function type.
            /// \tparam Archive The serialization archive type.
            /// \param[in,out] ar The serialization archive.
            template <typename fnT, typename Archive>
            typename std::enable_if<!is_any_function_pointer_v<fnT>>::type
            serialize_internal(const Archive& ar) {
                ar & ref & func & attr;
            }
        }; // struct TaskHandlerInfo

        /// Behave like a lazy \c std::enable_if.

        /// Evaluates to \c returnT if \c B is true, otherwise to an invalid type expression
        /// which causes the template expression in which it is used to not be considered for
        /// overload resolution. This "lazy" version is used if \c T is only valid when
        /// B is true. Note: typename T::type is the return type and must be well formed.
        /// \tparam B The bool value.
        /// \tparam returnT The type.
        template <bool B, class returnT>
        struct function_enabler_helper {
          typedef typename returnT::type type;
        };

        /// Specialization that disables \c type when \c B is false.

        /// \tparam returnT The type.
        template <class returnT>
        struct function_enabler_helper<false, returnT> { };

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam fnT Description needed.
        template <typename fnT>
        struct function_enabler : public
            function_enabler_helper<
                function_traits<fnT>::value || is_functor<fnT>::value,
                task_result_type<fnT> >
        { };

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam fnT Description needed.
        template <typename callableT, typename Enabler = void>
        struct callable_enabler;
        template <typename callableT>
        struct callable_enabler<callableT,
        std::enable_if_t<callable_traits<callableT>::value>>
        { using type = typename callable_traits<callableT>::result_type; };

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam objT Description needed.
        /// \tparam memfnT Description needed.
        /// \tparam enableT Description needed.
        template <typename objT, typename memfnT, typename enableT = void>
        struct memfunc_enabler_base
        { };

        /// \todo Brief description needed.

        /// \todo Descriptions needed
        /// \tparam objT Description needed.
        /// \tparam resT Description needed.
        /// \tparam baseT Description needed.
        /// \tparam paramT Description needed.
        template <typename objT, typename resT, typename baseT, typename ... paramT>
        struct memfunc_enabler_base<objT, resT (baseT::*)(paramT...),
            typename std::enable_if<std::is_base_of<baseT, objT>::value>::type >
        {
            /// \todo Brief description needed.
            typedef typename add_future<resT>::type type;
        };

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam objT Description needed.
        /// \tparam resT Description needed.
        /// \tparam baseT Description needed.
        /// \tparam paramT Description needed.
        template <typename objT, typename resT, typename baseT, typename ... paramT>
        struct memfunc_enabler_base<objT, resT (baseT::*)(paramT...) const,
            typename std::enable_if<std::is_base_of<baseT, objT>::value>::type >
        {
            /// \todo Brief description needed.
            typedef typename add_future<resT>::type type;
        };

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam objT Description needed.
        /// \tparam memfnT Description needed.
        template <typename objT, typename memfnT>
        struct memfunc_enabler :
                public memfunc_enabler_base<typename std::decay<objT>::type, memfnT>
        { };

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam objT Description needed.
        /// \tparam memfnT Description needed.
        template <typename objT, typename memfnT>
        struct memfunc_enabler<objT*, memfnT> :
            public memfunc_enabler<objT, memfnT>
        { };

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam objT Description needed.
        /// \tparam memfnT Description needed.
        template <typename objT, typename memfnT>
        struct memfunc_enabler<const objT*, memfnT> :
            public memfunc_enabler<objT, memfnT>
        { };

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam objT Description needed.
        /// \tparam memfnT Description needed.
        template <typename objT, typename memfnT>
        struct memfunc_enabler<objT* const, memfnT> :
            public memfunc_enabler<objT, memfnT>
        { };

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam objT Description needed.
        /// \tparam memfnT Description needed.
        template <typename objT, typename memfnT>
        struct memfunc_enabler<const objT* const, memfnT> :
            public memfunc_enabler<objT, memfnT>
        { };

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam objT Description needed.
        /// \tparam memfnT Description needed.
        template <typename objT, typename memfnT>
        struct memfunc_enabler<std::shared_ptr<objT>&, memfnT> :
            public memfunc_enabler<objT, memfnT>
        { };

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam objT Description needed.
        /// \tparam memfnT Description needed.
        template <typename objT, typename memfnT>
        struct memfunc_enabler<const std::shared_ptr<objT>&, memfnT> :
            public memfunc_enabler<objT, memfnT>
        { };

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam objT Description needed.
        /// \tparam memfnT Description needed.
        template <typename objT, typename memfnT>
        struct memfunc_enabler<std::shared_ptr<objT>, memfnT> :
            public memfunc_enabler<objT, memfnT>
        { };

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam objT Description needed.
        /// \tparam memfnT Description needed.
        template <typename objT, typename memfnT>
        struct memfunc_enabler<const std::shared_ptr<objT>, memfnT> :
            public memfunc_enabler<objT, memfnT>
        { };

    }  // namespace detail


    /// Multi-threaded queue to manage and run tasks.

    /// \todo A concise description of the inner workings...
    class WorldTaskQueue : public CallbackInterface, private NO_DEFAULTS {
        friend class TaskInterface;
    private:
        World& world; ///< The communication context.
        const ProcessID me; ///< This process.
        AtomicInt nregistered; ///< Count of pending tasks.

        /// \todo Brief description needed.
        void notify() {
            nregistered--;
        }

        /// \todo Brief description needed.

        /// This template is used in the reduce kernel.
        /// \todo Template parameter descriptions need verification.
        /// \tparam resultT Return type of the operation.
        /// \tparam opT The operation.
        /// \param[in] left Description needed.
        /// \param[in] right Description needed.
        /// \param[in] op The operation used in the reduce.
        /// \return Reduce of \c left and \c right using \c op.
        template <typename resultT, typename opT>
        static resultT sum(const resultT& left, const resultT& right, const opT& op) {
            //std::cout << " REDUCE SUM " << left << " " << right << std::endl;
            return op(left,right);
        }

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam taskT Description needed.
        /// \param arg Description needed.
        template <typename taskT>
        static void remote_task_handler(const AmArg& arg) {
            MADNESS_ASSERT(taskT::arity <= 9u);

            // Get task info and arguments form active message

            detail::TaskHandlerInfo<typename taskT::futureT::remote_refT,
                    typename taskT::functionT> info;

            archive::BufferInputArchive input_arch = arg & info;

            // Construct task
            taskT* task = new taskT(typename taskT::futureT(info.ref),
                    info.func, info.attr, input_arch);

            // Add task to queue
            arg.get_world()->taskq.add(task);
        }

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam T Description needed.
        /// \param[in] f Description needed.
        template <typename T>
        inline const T& am_arg(const Future<T>& f) {
            MADNESS_ASSERT(f.probe());
            return f.get();
        }

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam T Description needed.
        /// \param[in] t Description needed.
        /// \return Description needed.
        template <typename T>
        inline const T& am_arg(const T& t) {
            return t;
        }

        /// \todo Brief description needed.
        typedef detail::voidT voidT;

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam taskT Description needed.
        /// \tparam fnT Description needed.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \tparam a6T Type of argument 6.
        /// \tparam a7T Type of argument 7.
        /// \tparam a8T Type of argument 8.
        /// \tparam a9T Type of argument 9.
        /// \param where Description needed.
        /// \param fn Description needed.
        /// \param[in] a1 Argument 1.
        /// \param[in] a2 Argument 2.
        /// \param[in] a3 Argument 3.
        /// \param[in] a4 Argument 4.
        /// \param[in] a5 Argument 5.
        /// \param[in] a6 Argument 6.
        /// \param[in] a7 Argument 7.
        /// \param[in] a8 Argument 8.
        /// \param[in] a9 Argument 9.
        /// \param[in] attr Description needed.
        /// \return Description needed.
        template <typename taskT, typename fnT, typename a1T, typename a2T, typename a3T,
                typename a4T, typename a5T, typename a6T, typename a7T,
                typename a8T, typename a9T>
        inline typename taskT::futureT
        send_task(ProcessID where, fnT fn, const a1T& a1,
                const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5,
                const a6T& a6, const a7T& a7, const a8T& a8, const a9T& a9,
                const TaskAttributes& attr)
        {
            typename taskT::futureT result;
            typedef detail::TaskHandlerInfo<typename taskT::futureT::remote_refT, typename taskT::functionT> infoT;
            world.am.send(where, & WorldTaskQueue::template remote_task_handler<taskT>,
                    new_am_arg(infoT(result.remote_ref(world), fn, attr),
                    a1, a2, a3, a4, a5, a6, a7, a8, a9), RMI::ATTR_UNORDERED);

            return result;
        }


    public:
        /// Constructor requiring a communication context (\c World).

        /// \param[in,out] world The communication context.
        WorldTaskQueue(World& world);

        /// Returns the number of pending tasks.

        /// \return The number of pending tasks.
        size_t size() const {
            return nregistered;
        }


        /// Add a new local task, taking ownership of the pointer.

        /// The task pointer (\c t) is assumed to have been created with
        /// \c new and when the task is eventually run the queue
        /// will call the task's destructor using \c delete.
        ///
        /// Once the task is complete it will execute
        /// \c task_complete_callback to decrement the number of pending
        /// tasks and be deleted.
        /// \param[in] t Pointer to the task.
        void add(TaskInterface* t)  {
            nregistered++;

            t->set_info(&world, this);       // Stuff info

            // Always use the callback to avoid race condition
            t->register_submit_callback();
        }

        /// \todo Brief description needed.

        /// \todo Descriptions needed.
        /// \tparam fnT Description needed.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \tparam a6T Type of argument 6.
        /// \tparam a7T Type of argument 7.
        /// \tparam a8T Type of argument 8.
        /// \tparam a9T Type of argument 9.
        /// \param[in] t Description needed.
        /// \return Description needed.
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

        /// Reduce `op(item)` for all items in range using `op(sum,op(item))`.

        /// The operation must provide the following interface, of
        /// which the \c operator() methods are required by \c reduce()
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
        /// \note The serialize method does not actually have to
        /// work unless you want to have the task be stealable.
        ///
        /// Adjust the chunksize in the range to control granularity.
        /// \todo Descriptions needed and/or verified.
        /// \tparam resultT The result type of the operation.
        /// \tparam rangeT Description needed.
        /// \tparam opT Function type of the operation.
        /// \param[in] range The range of items.
        /// \param[in] op The operation.
        /// \return Description needed.
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

        /// Apply `op(item)` on all items in range.

        /// The operation must provide the following interface, of
        /// which the \c operator() method is required by `for_each()`
        /// and the rest by the task interface.
        /// \code
        /// struct opT {
        ///     opT(const opT&);
        ///     bool operator()(const rangeT::iterator& it) const;
        /// };
        /// \endcode
        /// \note The serialize method does not actually have to
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
        /// \todo Descriptions needed and/or verified.
        /// \tparam rangeT Description needed.
        /// \tparam opT Funtion type of the operation. This function should
        ///     have a return type of \c bool.
        /// \param[in] range The range of items.
        /// \param[in] op The operation.
        /// \return Future for a bool which is the logical `and` of all `op(item)` calls.
        template <typename rangeT, typename opT>
        Future<bool> for_each(const rangeT& range, const opT& op) {
            detail::ForEachRootTask<rangeT, opT>* for_each_root =
                    new detail::ForEachRootTask<rangeT, opT>(world, range, op);
            Future<bool> result = for_each_root->result();
            add(for_each_root);
            return result;
        }

#if MADNESS_TASKQ_VARIADICS

        ///////////////////////////////////////////////////////////////////////////////

        /// Create a local task with one argument.

        /// Creates a task in this process. An argument that is a future may be
        /// used to carry dependencies.
        /// \tparam fnT A function pointer or functor.
        /// \tparam a1T Type of argument 1.
        /// \param[in,out] fn The function to be called in the task.
        /// \param[in] a1 Argument 1.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return
        ///     type is \c void, a \c Future<void> object is returned that may
        ///     be ignored.
        /// \note future_to_ref_t is used instead of remove_future_t so that
        ///       argument Future's convert to refs so that X(Y&) can invoked by add(X, Future<T>);
        ///       if remove_future_t were used, X(Y) would be expected instead.
        ///       This is reasonable if we remember the shallow copy semantics of Futures: they act
        ///       as refs.
        template <typename fnT, typename... argsT,
                  typename = std::enable_if_t<
                  meta::taskattr_is_last_arg<argsT...>::value>>
        typename meta::drop_last_arg_and_apply_callable<detail::function_enabler, fnT, future_to_ref_t<argsT>...>::type::type
        add(fnT&& fn, argsT&&... args) {
          using taskT = typename meta::drop_last_arg_and_apply<TaskFn,
                                                               std::decay_t<fnT>,
                                                               std::remove_cv_t<std::remove_reference_t<argsT>>...>::type;
          return add(new taskT(typename taskT::futureT(), std::forward<fnT>(fn),
                               std::forward<argsT>(args)...));
        }

        template <typename fnT, typename... argsT,
                  typename = std::enable_if_t<
                      !meta::taskattr_is_last_arg<argsT...>::value>>
        typename detail::function_enabler<fnT(future_to_ref_t<argsT>...)>::type add(
            fnT&& fn, argsT&&... args) {
          using taskT = TaskFn<std::decay_t<fnT>, std::remove_cv_t<std::remove_reference_t<argsT>>...>;
          return add(new taskT(typename taskT::futureT(), std::forward<fnT>(fn),
                               std::forward<argsT>(args)..., TaskAttributes()));
        }

#else
        /// Create a local task with no arguments.

        /// Creates a task in this process. An argument that is a future may be
        /// used to carry dependencies.
        /// \tparam fnT A function pointer or functor.
        /// \param[in,out] fn The function to be called in the task.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return
        ///     type is \c void, a \c Future<void> object is returned that may
        ///     be ignored.
        template <typename fnT>
        typename detail::function_enabler<fnT()>::type
        add(fnT fn, const TaskAttributes& attr=TaskAttributes()) {
            typedef TaskFn<fnT> taskT;
            return add(new taskT(typename taskT::futureT(),
                    fn, attr));
        }

        /// Create a local task with one argument.

        /// Creates a task in this process. An argument that is a future may be
        /// used to carry dependencies.
        /// \tparam fnT A function pointer or functor.
        /// \tparam a1T Type of argument 1.
        /// \param[in,out] fn The function to be called in the task.
        /// \param[in] a1 Argument 1.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return
        ///     type is \c void, a \c Future<void> object is returned that may
        ///     be ignored.
        template <typename fnT, typename a1T>
        typename detail::function_enabler<fnT(a1T)>::type
        add(fnT fn, const a1T& a1, const TaskAttributes& attr=TaskAttributes()) {
            typedef TaskFn<fnT, a1T> taskT;
            return add(new taskT(typename taskT::futureT(),
                    fn, a1, attr));
        }

        /// Create a local task with two arguments.

        /// Creates a task in this process. An argument that is a future may be
        /// used to carry dependencies.
        /// \tparam fnT A function pointer or functor.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \param[in,out] fn The function to be called in the task.
        /// \param[in] a1 Argument 1.
        /// \param[in] a2 Argument 2.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return
        ///     type is \c void, a \c Future<void> object is returned that may
        ///     be ignored.
        template <typename fnT, typename a1T, typename a2T>
        typename detail::function_enabler<fnT(a1T, a2T)>::type
        add(fnT fn, const a1T& a1, const a2T& a2, const TaskAttributes& attr=TaskAttributes()) {
            typedef TaskFn<fnT, a1T, a2T> taskT;
            return add(new taskT(typename taskT::futureT(),
                    fn, a1, a2, attr));
        }

        /// Create a local task with three arguments.

        /// Creates a task in this process. An argument that is a future may be
        /// used to carry dependencies.
        /// \tparam fnT A function pointer or functor.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \param[in,out] fn The function to be called in the task.
        /// \param[in] a1 Argument 1.
        /// \param[in] a2 Argument 2.
        /// \param[in] a3 Argument 3.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return
        ///     type is \c void, a \c Future<void> object is returned that may
        ///     be ignored.
        template <typename fnT, typename a1T, typename a2T, typename a3T>
        typename detail::function_enabler<fnT(a1T, a2T, a3T)>::type
        add(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3,
            const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T> taskT;
            return add(new taskT(typename taskT::futureT(),
                    fn, a1, a2, a3, attr));
        }

        /// Create a local task with four arguments.

        /// Creates a task in this process. An argument that is a future may be
        /// used to carry dependencies.
        /// \tparam fnT A function pointer or functor.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \param[in,out] fn The function to be called in the task.
        /// \param[in] a1 Argument 1.
        /// \param[in] a2 Argument 2.
        /// \param[in] a3 Argument 3.
        /// \param[in] a4 Argument 4.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return
        ///     type is \c void, a \c Future<void> object is returned that may
        ///     be ignored.
        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T>
        typename detail::function_enabler<fnT(a1T, a2T, a3T, a4T)>::type
        add(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
            const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T> taskT;
            return add(new taskT(typename taskT::futureT(),
                    fn, a1, a2, a3, a4, attr));
        }

        /// Create a local task with five arguments.

        /// Creates a task in this process. An argument that is a future may be
        /// used to carry dependencies.
        /// \tparam fnT A function pointer or functor.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \param[in,out] fn The function to be called in the task.
        /// \param[in] a1 Argument 1.
        /// \param[in] a2 Argument 2.
        /// \param[in] a3 Argument 3.
        /// \param[in] a4 Argument 4.
        /// \param[in] a5 Argument 5.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return
        ///     type is \c void, a \c Future<void> object is returned that may
        ///     be ignored.
        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
            typename a5T>
        typename detail::function_enabler<fnT(a1T, a2T, a3T, a4T, a5T)>::type
        add(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
            const a5T& a5, const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T> taskT;
            return add(new taskT(typename taskT::futureT(),
                    fn, a1, a2, a3, a4, a5, attr));
        }

        /// Create a local task with six arguments.

        /// Creates a task in this process. An argument that is a future may be
        /// used to carry dependencies.
        /// \tparam fnT A function pointer or functor.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \tparam a6T Type of argument 6.
        /// \param[in,out] fn The function to be called in the task.
        /// \param[in] a1 Argument 1.
        /// \param[in] a2 Argument 2.
        /// \param[in] a3 Argument 3.
        /// \param[in] a4 Argument 4.
        /// \param[in] a5 Argument 5.
        /// \param[in] a6 Argument 6.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return
        ///     type is \c void, a \c Future<void> object is returned that may
        ///     be ignored.
        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
            typename a5T, typename a6T>
        typename detail::function_enabler<fnT(a1T, a2T, a3T, a4T, a5T, a6T)>::type
        add(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
            const a5T& a5, const a6T& a6, const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T> taskT;
            return add(new taskT(typename taskT::futureT(),
                    fn, a1, a2, a3, a4, a5, a6, attr));
        }

        /// Create a local task with seven arguments.

        /// Creates a task in this process. An argument that is a future may be
        /// used to carry dependencies.
        /// \tparam fnT A function pointer or functor.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \tparam a6T Type of argument 6.
        /// \tparam a7T Type of argument 7.
        /// \param[in,out] fn The function to be called in the task.
        /// \param[in] a1 Argument 1.
        /// \param[in] a2 Argument 2.
        /// \param[in] a3 Argument 3.
        /// \param[in] a4 Argument 4.
        /// \param[in] a5 Argument 5.
        /// \param[in] a6 Argument 6.
        /// \param[in] a7 Argument 7.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return
        ///     type is \c void, a \c Future<void> object is returned that may
        ///     be ignored.
        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
            typename a5T, typename a6T, typename a7T>
        typename detail::function_enabler<fnT(a1T, a2T, a3T, a4T, a5T, a6T, a7T)>::type
        add(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
            const a5T& a5, const a6T& a6, const a7T& a7,
            const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T> taskT;
            return add(new taskT(typename taskT::futureT(),
                    fn, a1, a2, a3, a4, a5, a6, a7, attr));
        }

        /// Create a local task with eight arguments.

        /// Creates a task in this process. An argument that is a future may be
        /// used to carry dependencies.
        /// \tparam fnT A function pointer or functor.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \tparam a6T Type of argument 6.
        /// \tparam a7T Type of argument 7.
        /// \tparam a8T Type of argument 8.
        /// \param[in,out] fn The function to be called in the task.
        /// \param[in] a1 Argument 1.
        /// \param[in] a2 Argument 2.
        /// \param[in] a3 Argument 3.
        /// \param[in] a4 Argument 4.
        /// \param[in] a5 Argument 5.
        /// \param[in] a6 Argument 6.
        /// \param[in] a7 Argument 7.
        /// \param[in] a8 Argument 8.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return
        ///     type is \c void, a \c Future<void> object is returned that may
        ///     be ignored.
        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
            typename a5T, typename a6T, typename a7T, typename a8T>
        typename detail::function_enabler<fnT(a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T)>::type
        add(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
            const a5T& a5, const a6T& a6, const a7T& a7, const a8T& a8,
            const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T> taskT;
            return add(new taskT(typename taskT::futureT(),
                    fn, a1, a2, a3, a4, a5, a6, a7, a8, attr));
        }

        /// Create a local task with nine arguments.

        /// Creates a task in this process. An argument that is a future may be
        /// used to carry dependencies.
        /// \tparam fnT A function pointer or functor.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \tparam a6T Type of argument 6.
        /// \tparam a7T Type of argument 7.
        /// \tparam a8T Type of argument 8.
        /// \tparam a9T Type of argument 9.
        /// \param[in,out] fn The function to be called in the task.
        /// \param[in] a1 Argument 1.
        /// \param[in] a2 Argument 2.
        /// \param[in] a3 Argument 3.
        /// \param[in] a4 Argument 4.
        /// \param[in] a5 Argument 5.
        /// \param[in] a6 Argument 6.
        /// \param[in] a7 Argument 7.
        /// \param[in] a8 Argument 8.
        /// \param[in] a9 Argument 9.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return
        ///     type is \c void, a \c Future<void> object is returned that may
        ///     be ignored.
        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
            typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
        typename detail::function_enabler<fnT(a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T)>::type
        add(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
            const a5T& a5, const a6T& a6, const a7T& a7, const a8T& a8, const a9T& a9,
            const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T> taskT;
            return add(new taskT(typename taskT::futureT(),
                    fn, a1, a2, a3, a4, a5, a6, a7, a8, a9, attr));
        }
#endif

        /// Create a remote task.

        /// Creates a task in process \c dest, which may or may not be this
        /// process. An argument that is a future may be used to carry
        /// dependencies for local tasks. A future that is not ready cannot be
        /// used as an argument for remote tasks -- i.e., remote tasks must
        /// be ready to execute (you can work around this by making a local task
        /// to submit the remote task once everything is ready).
        /// \tparam fnT A function pointer or functor type.
        /// \param[in] dest The process where the task will be created.
        /// \param[in] fn The function to be called in the task.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return type
        /// is \c void, a \c Future<void> object is return that may
        /// be ignored.
        /// \note Arguments must be (de)serializable and must of course make
        /// sense at the remote destination. Fundamental types, simple STL
        /// containers, and pointers to \c World, \c WorldContainer, and
        /// user-defined types derived from \c WorldObject<> are automatically
        /// handled. Anything else is your problem.
        template <typename fnT>
        typename detail::function_enabler<fnT>::type
        add(ProcessID dest, fnT fn, const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT> taskT;
            if(dest == me)
                return add(fn, attr);
            else
                return send_task<taskT>(dest, fn, voidT::value, voidT::value,
                        voidT::value, voidT::value, voidT::value, voidT::value,
                        voidT::value, voidT::value, voidT::value, attr);
        }

        /// Invoke `resultT (*fn)(a1T)` as a task, local or remote.

        /// A future is returned to hold the eventual result of the task.
        /// \c Future<void> is an empty class that may be ignored.
        /// \tparam fnT A function pointer or functor type.
        /// \tparam a1T Type of argument 1.
        /// \param[in] dest The process where the task will be created.
        /// \param[in] fn The function to be called in the task.
        /// \param[in] a1 Argument 1.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return type
        /// is \c void, a \c Future<void> object is return that may
        /// be ignored.
        /// \note Arguments must be (de)serializable and must of course make
        /// sense at the remote destination. Fundamental types,
        /// simple STL containers, and pointers to \c World,
        /// \c WorldContainer, and user-defined types derived from
        /// \c WorldObject are automatically handled. Anything else is
        /// your problem.
        template <typename fnT, typename a1T>
        typename detail::function_enabler<fnT>::type
        add(ProcessID dest, fnT fn, const a1T& a1, const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T> taskT;
            if(dest == me)
                return add(fn, a1, attr);
            else
                return send_task<taskT>(dest, fn, am_arg(a1), voidT::value,
                        voidT::value, voidT::value, voidT::value, voidT::value,
                        voidT::value, voidT::value, voidT::value, attr);
        }

        /// Invoke `resultT (*fn)(a1T,a2T)` as a task, local or remote.

        /// A future is returned to hold the eventual result of the task.
        /// \c Future<void> is an empty class that may be ignored.
        /// \tparam fnT A function pointer or functor type.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \param[in] dest The process where the task will be created.
        /// \param[in] fn The function to be called in the task.
        /// \param[in] a1 Argument 1.
        /// \param[in] a2 Argument 2.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return type
        /// is \c void, a \c Future<void> object is return that may
        /// be ignored.
        /// \note Arguments must be (de)serializable and must of course make
        /// sense at the remote destination. Fundamental types,
        /// simple STL containers, and pointers to \c World,
        /// \c WorldContainer, and user-defined types derived from
        /// \c WorldObject are automatically handled. Anything else is
        /// your problem.
        template <typename fnT, typename a1T, typename a2T>
        typename detail::function_enabler<fnT>::type
        add(ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
                const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T> taskT;
            if(dest == me)
                return add(fn, a1, a2, attr);
            else
                return send_task<taskT>(dest, fn, am_arg(a1), am_arg(a2),
                        voidT::value, voidT::value, voidT::value, voidT::value,
                        voidT::value, voidT::value, voidT::value, attr);
        }

        /// Invoke `resultT (*fn)(a1T,a2T,a3T)` as a task, local or remote.

        /// A future is returned to hold the eventual result of the task.
        /// \c Future<void> is an empty class that may be ignored.
        /// \tparam fnT A function pointer or functor type.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \param[in] dest The process where the task will be created.
        /// \param[in] fn The function to be called in the task.
        /// \param[in] a1 Argument 1.
        /// \param[in] a2 Argument 2.
        /// \param[in] a3 Argument 3.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return type
        /// is \c void, a \c Future<void> object is return that may
        /// be ignored.
        /// \note Arguments must be (de)serializable and must of course make
        /// sense at the remote destination. Fundamental types,
        /// simple STL containers, and pointers to \c World,
        /// \c WorldContainer, and user-defined types derived from
        /// \c WorldObject are automatically handled. Anything else is
        /// your problem.
        template <typename fnT, typename a1T, typename a2T, typename a3T>
        typename detail::function_enabler<fnT>::type
        add(ProcessID dest, fnT fn, const a1T& a1, const a2T& a2, const a3T& a3,
                const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T> taskT;
            if(dest == me)
                return add(fn, a1, a2, a3, attr);
            else
                return send_task<taskT>(dest, fn, am_arg(a1), am_arg(a2),
                        am_arg(a3), voidT::value, voidT::value, voidT::value,
                        voidT::value, voidT::value, voidT::value, attr);
        }

        /// Invoke `resultT (*fn)(a1T,a2T,a3T,a4T)` as a task, local or remote.

        /// A future is returned to hold the eventual result of the task.
        /// \c Future<void> is an empty class that may be ignored.
        /// \tparam fnT A function pointer or functor type.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \param[in] dest The process where the task will be created.
        /// \param[in] fn The function to be called in the task.
        /// \param[in] a1 Argument 1.
        /// \param[in] a2 Argument 2.
        /// \param[in] a3 Argument 3.
        /// \param[in] a4 Argument 4.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return type
        /// is \c void, a \c Future<void> object is return that may
        /// be ignored.
        /// \note Arguments must be (de)serializable and must of course make
        /// sense at the remote destination. Fundamental types,
        /// simple STL containers, and pointers to \c World,
        /// \c WorldContainer, and user-defined types derived from
        /// \c WorldObject are automatically handled. Anything else is
        /// your problem.
        template <typename fnT, typename a1T, typename a2T, typename a3T,
                typename a4T>
        typename detail::function_enabler<fnT>::type
        add(ProcessID dest, fnT fn, const a1T& a1, const a2T& a2, const a3T& a3,
                const a4T& a4, const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T> taskT;
            if(dest == me)
                return add(fn, a1, a2, a3, a4, attr);
            else
                return send_task<taskT>(dest, fn, am_arg(a1), am_arg(a2),
                        am_arg(a3), am_arg(a4), voidT::value, voidT::value,
                        voidT::value, voidT::value, voidT::value, attr);
        }

        /// Invoke `resultT (*fn)(a1T,a2T,a3T,a4T,a5T)` as a task, local or remote.

        /// A future is returned to hold the eventual result of the task.
        /// \c Future<void> is an empty class that may be ignored.
        /// \tparam fnT A function pointer or functor type.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \param[in] dest The process where the task will be created.
        /// \param[in] fn The function to be called in the task.
        /// \param[in] a1 Argument 1.
        /// \param[in] a2 Argument 2.
        /// \param[in] a3 Argument 3.
        /// \param[in] a4 Argument 4.
        /// \param[in] a5 Argument 5.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return type
        /// is \c void, a \c Future<void> object is return that may
        /// be ignored.
        /// \note Arguments must be (de)serializable and must of course make
        /// sense at the remote destination. Fundamental types,
        /// simple STL containers, and pointers to \c World,
        /// \c WorldContainer, and user-defined types derived from
        /// \c WorldObject are automatically handled. Anything else is
        /// your problem.
        template <typename fnT, typename a1T, typename a2T, typename a3T,
                typename a4T, typename a5T>
        typename detail::function_enabler<fnT>::type
        add(ProcessID dest, fnT fn, const a1T& a1, const a2T& a2, const a3T& a3,
                const a4T& a4, const a5T& a5, const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T> taskT;
            if(dest == me)
                return add(fn, a1, a2, a3, a4, a5, attr);
            else
                return send_task<taskT>(dest, fn, am_arg(a1), am_arg(a2),
                        am_arg(a3), am_arg(a4), am_arg(a5), voidT::value,
                        voidT::value, voidT::value, voidT::value, attr);
        }

        /// Invoke `resultT (*fn)(a1T,a2T,a3T,a4T,a5T,a6T)` as a task, local or remote.

        /// A future is returned to hold the eventual result of the task.
        /// \c Future<void> is an empty class that may be ignored.
        /// \tparam fnT A function pointer or functor type.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \tparam a6T Type of argument 6.
        /// \param[in] dest The process where the task will be created.
        /// \param[in] fn The function to be called in the task.
        /// \param[in] a1 Argument 1.
        /// \param[in] a2 Argument 2.
        /// \param[in] a3 Argument 3.
        /// \param[in] a4 Argument 4.
        /// \param[in] a5 Argument 5.
        /// \param[in] a6 Argument 6.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return type
        /// is \c void, a \c Future<void> object is return that may
        /// be ignored.
        /// \note Arguments must be (de)serializable and must of course make
        /// sense at the remote destination. Fundamental types,
        /// simple STL containers, and pointers to \c World,
        /// \c WorldContainer, and user-defined types derived from
        /// \c WorldObject are automatically handled. Anything else is
        /// your problem.
        template <typename fnT, typename a1T, typename a2T, typename a3T,
                typename a4T, typename a5T, typename a6T>
        typename detail::function_enabler<fnT>::type
        add(ProcessID dest, fnT fn, const a1T& a1, const a2T& a2, const a3T& a3,
                const a4T& a4, const a5T& a5, const a6T& a6,
                const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T> taskT;
            if(dest == me)
                return add(fn, a1, a2, a3, a4, a5, a6, attr);
            else
                return send_task<taskT>(dest, fn, am_arg(a1), am_arg(a2),
                        am_arg(a3), am_arg(a4), am_arg(a5), am_arg(a6),
                        voidT::value, voidT::value, voidT::value, attr);
        }

        /// Invoke `resultT (*fn)(a1T,a2T,a3T,a4T,a5T,a6T,a7T)` as a task, local or remote.

        /// A future is returned to hold the eventual result of the task.
        /// \c Future<void> is an empty class that may be ignored.
        /// \tparam fnT A function pointer or functor type.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \tparam a6T Type of argument 6.
        /// \tparam a7T Type of argument 7.
        /// \param[in] dest The process where the task will be created.
        /// \param[in] fn The function to be called in the task.
        /// \param[in] a1 Argument 1.
        /// \param[in] a2 Argument 2.
        /// \param[in] a3 Argument 3.
        /// \param[in] a4 Argument 4.
        /// \param[in] a5 Argument 5.
        /// \param[in] a6 Argument 6.
        /// \param[in] a7 Argument 7.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return type
        /// is \c void, a \c Future<void> object is return that may
        /// be ignored.
        /// \note Arguments must be (de)serializable and must of course make
        /// sense at the remote destination. Fundamental types,
        /// simple STL containers, and pointers to \c World,
        /// \c WorldContainer, and user-defined types derived from
        /// \c WorldObject are automatically handled. Anything else is
        /// your problem.
        template <typename fnT, typename a1T, typename a2T, typename a3T,
                typename a4T, typename a5T, typename a6T, typename a7T>
        typename detail::function_enabler<fnT>::type
        add(ProcessID dest, fnT fn, const a1T& a1, const a2T& a2, const a3T& a3,
                const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7,
                const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T> taskT;
            if(dest == me)
                return add(fn, a1, a2, a3, a4, a5, a6, a7, attr);
            else
                return send_task<taskT>(dest, fn, am_arg(a1), am_arg(a2),
                        am_arg(a3), am_arg(a4), am_arg(a5), am_arg(a6),
                        am_arg(a7), voidT::value, voidT::value, attr);
        }

        /// Invoke `resultT (*fn)(a1T,a2T,a3T,a4T,a5T,a6T,a7T,a8T)` as a task, local or remote.

        /// A future is returned to hold the eventual result of the task.
        /// \c Future<void> is an empty class that may be ignored.
        /// \tparam fnT A function pointer or functor type.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \tparam a6T Type of argument 6.
        /// \tparam a7T Type of argument 7.
        /// \tparam a8T Type of argument 8.
        /// \param[in] dest The process where the task will be created.
        /// \param[in] fn The function to be called in the task.
        /// \param[in] a1 Argument 1.
        /// \param[in] a2 Argument 2.
        /// \param[in] a3 Argument 3.
        /// \param[in] a4 Argument 4.
        /// \param[in] a5 Argument 5.
        /// \param[in] a6 Argument 6.
        /// \param[in] a7 Argument 7.
        /// \param[in] a8 Argument 8.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return type
        /// is \c void, a \c Future<void> object is return that may
        /// be ignored.
        /// \note Arguments must be (de)serializable and must of course make
        /// sense at the remote destination. Fundamental types,
        /// simple STL containers, and pointers to \c World,
        /// \c WorldContainer, and user-defined types derived from
        /// \c WorldObject are automatically handled. Anything else is
        /// your problem.
        template <typename fnT, typename a1T, typename a2T, typename a3T,
                typename a4T, typename a5T, typename a6T, typename a7T,
                typename a8T>
        typename detail::function_enabler<fnT>::type
        add(ProcessID dest, fnT fn, const a1T& a1, const a2T& a2, const a3T& a3,
                const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7,
                const a8T& a8, const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T> taskT;
            if(dest == me)
                return add(fn, a1, a2, a3, a4, a5, a6, a7, a8, attr);
            else
                return send_task<taskT>(dest, fn, am_arg(a1), am_arg(a2),
                        am_arg(a3), am_arg(a4), am_arg(a5), am_arg(a6),
                        am_arg(a7), am_arg(a8), voidT::value, attr);
        }

        /// Invoke `resultT (*fn)(a1T,a2T,a3T,a4T,a5T,a6T,a7T,a8T,a9T)` as a task, local or remote.

        /// A future is returned to hold the eventual result of the task.
        /// \c Future<void> is an empty class that may be ignored.
        /// \tparam fnT A function pointer or functor type.
        /// \tparam a1T Type of argument 1.
        /// \tparam a2T Type of argument 2.
        /// \tparam a3T Type of argument 3.
        /// \tparam a4T Type of argument 4.
        /// \tparam a5T Type of argument 5.
        /// \tparam a6T Type of argument 6.
        /// \tparam a7T Type of argument 7.
        /// \tparam a8T Type of argument 8.
        /// \tparam a9T Type of argument 9.
        /// \param[in] dest The process where the task will be created.
        /// \param[in] fn The function to be called in the task.
        /// \param[in] a1 Argument 1.
        /// \param[in] a2 Argument 2.
        /// \param[in] a3 Argument 3.
        /// \param[in] a4 Argument 4.
        /// \param[in] a5 Argument 5.
        /// \param[in] a6 Argument 6.
        /// \param[in] a7 Argument 7.
        /// \param[in] a8 Argument 8.
        /// \param[in] a9 Argument 9.
        /// \param[in] attr The task attributes.
        /// \return A future to the result. If the task function return type
        /// is \c void, a \c Future<void> object is return that may
        /// be ignored.
        /// \note Arguments must be (de)serializable and must of course make
        /// sense at the remote destination. Fundamental types,
        /// simple STL containers, and pointers to \c World,
        /// \c WorldContainer, and user-defined types derived from
        /// \c WorldObject are automatically handled. Anything else is
        /// your problem.
        /// \todo Could we use metaprogramming or variadic templates to reduce all of these instances to one template that generates versions for one more (less) parameter?
        template <typename fnT, typename a1T, typename a2T, typename a3T,
                typename a4T, typename a5T, typename a6T, typename a7T,
                typename a8T, typename a9T>
        typename detail::function_enabler<fnT>::type
        add(ProcessID dest, fnT fn, const a1T& a1, const a2T& a2, const a3T& a3,
                const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7,
                const a8T& a8, const a9T& a9, const TaskAttributes& attr=TaskAttributes())
        {
            typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T> taskT;
            if(dest == me)
                return add(fn, a1, a2, a3, a4, a5, a6, a7, a8, a9, attr);
            else
                return send_task<taskT>(dest, fn, am_arg(a1), am_arg(a2),
                        am_arg(a3), am_arg(a4), am_arg(a5), am_arg(a6),
                        am_arg(a7), am_arg(a8), am_arg(a9), attr);
        }

        /// Invoke `resultT (obj.*memfn)(args...)` as a local task.

        /// \todo Verify this documentation.
        /// A future is returned to hold the eventual result of the task.
        /// \c Future<void> is an empty class that may be ignored.
        /// \tparam objT The object type.
        /// \tparam memfnT The member function type.
        /// \tparam argT Variadic template for arguments.
        /// \param[in] obj The associated object for invoking the member function pointer.
        /// \param[in] memfn The member function pointer.
        /// \param[in] args The argument pack.
        /// \return A future to the result. If the task function return type
        /// is \c void, a \c Future<void> object is return that may
        /// be ignored.
        /// \note Arguments must be (de)serializable and must of course make
        /// sense at the remote destination. Fundamental types,
        /// simple STL containers, and pointers to \c World,
        /// \c WorldContainer, and user-defined types derived from
        /// \c WorldObject are automatically handled. Anything else is
        /// your problem.
        template <typename objT, typename memfnT, typename... argT>
        typename detail::memfunc_enabler<objT, memfnT>::type
        add(objT&& obj, memfnT memfn, argT&&... args)
        { return add(detail::wrap_mem_fn(std::forward<objT>(obj),memfn), std::forward<argT>(args)...); }



    private:

        /// \todo Brief description needed.
        struct ProbeAllDone {
            WorldTaskQueue* tq; ///< The task queue.

            /// Construct a \c ProbeAllDone for a given task queue.

            /// \param[in] tq Pointer to the task queue.
            ProbeAllDone(WorldTaskQueue* tq)
                : tq(tq)
            {}

            /// Determine if all tasks in the queue are complete.

            /// \return True if all tasks are complete; false otherwise.
            bool operator()() const {
                return (tq->nregistered == 0);
            }
        };

    public:

        /// Returns after all local tasks have completed.

        /// While waiting, the calling thread will run tasks.
        void fence() {
            try {
                ThreadPool::await(ProbeAllDone(this), true);
            } catch(...) {
                fprintf(stderr, "!!MADNESS ERROR: Exception thrown in WorldTaskQueue::fence() with %i pending task(s)\n", int(nregistered));
                throw;
            }
        }
    };

    namespace detail {

#ifdef HAVE_INTEL_TBB

        /// Apply an operation to a range of iterators.

        /// \tparam rangeT The range of iterators type.
        /// \tparam opT The operation type.
        /// This task creates `for_each` tasks and collects information on the
        /// results of those tasks. Once all tasks have completed it will set
        /// the result future.
        template <typename rangeT, typename opT>
        class ForEachRootTask : public TaskInterface {
        private:
            rangeT range_; ///< The range.
            opT op_; ///< The foreach function.
            Future<bool> completion_status_; ///< The result of this set of tasks.

        public:

            /// Constructor.

            /// \param[in] world The world where the tasks will run.
            /// \param[in] range The range of iterators.
            /// \param[in] op The operation that will be applied to the range of iterators.
            ForEachRootTask(World&, const rangeT range, const opT& op) :
                TaskInterface(0, TaskAttributes::hipri()), range_(range), op_(op)
            { }

            /// Virtual destructor.
            virtual ~ForEachRootTask() = default;

            /// Result accessor.

            /// \return A const reference to the result future.
            const Future<bool>& result() const { return completion_status_; }

            /// Task run work.

            /// Sets the result future based on the status of all iterations.
            virtual void run(World&, const TaskThreadEnv& ) /*override*/ {

                // Note: We use parallel_reduce instead of parallel_foreach to
                // accumulate result flags. Otherwise, this logically functions
                // like parallel_foreach.
                const bool result =
                    tbb::parallel_reduce(range_, true,
                        [this](const rangeT& r, bool init) -> bool {
                            for(typename rangeT::iterator it = r.begin(); it != r.end();  ++it)
                                init = init && op_(it);
                            return init;
                        },
                        [](const bool left, const bool right) -> bool {
                            return left && right;
                        }
                    );
                completion_status_.set(result);
            }

        private:
            /// Get the task ID.

            /// \todo Verify that \c id is an output parameter.
            /// \param[out] id The ID to set for this task.
            virtual void get_id(std::pair<void*,unsigned short>& id) const {
                PoolTaskInterface::make_id(id, *this);
            }
        }; // class ForEachRootTask

#else

        /// Apply an operation to a range of iterators.

        /// This task will progressively split range, creating leaves for each
        /// task, until the range of iterators is smaller than the chunk size.
        /// \tparam rangeT The range of iterators type.
        /// \tparam opT The operation type.
        template <typename rangeT, typename opT>
        class ForEachTask : public TaskInterface {
        private:
            rangeT range_; ///< The range of iterators.
            opT op_; ///< The operation to apply to range.
            ForEachRootTask<rangeT, opT>& root_; ///< The root task that signals completion and status.

            // not allowed
            ForEachTask(const ForEachTask<rangeT, opT>& fet) = delete;
            ForEachTask& operator=(const ForEachTask<rangeT, opT>& fet) = delete;

        public:

            /// Constructor.

            /// \todo Descriptions needed.
            /// \param[in] range The range of tasks.
            /// \param[in] op The operation to perform.
            /// \param root Description needed.
            ForEachTask(const rangeT range, const opT& op, ForEachRootTask<rangeT, opT>& root) :
                TaskInterface(0, TaskAttributes::hipri()), range_(range), op_(op), root_(root)
            {
                // Increment the master task dependency counter for this task
                root_.inc();
            }

            /// Virtual destructor.
            virtual ~ForEachTask() = default;

            /// Run the task.

            /// \todo Descriptions needed.
            /// \param[in] tte Description needed.
            virtual void run(const TaskThreadEnv& tte) {
                // Create leaf tasks and split range until it is less than chuncksize
                while(range_.size() > range_.get_chunksize()) {
                    rangeT right(range_,Split());
                    ForEachTask<rangeT,opT>* leaf = new ForEachTask<rangeT,opT>(right, op_, root_);
                    root_.world().taskq.add(leaf);
                }

                // Iterate over the remaining chunck of range and call op_ for each element
                int status = 0;
                for(typename rangeT::iterator it = range_.begin(); it != range_.end();  ++it)
                    if(op_(it))
                        ++status;

                // Notify the root task that this task is done give the status
                root_.complete(status);
            }

        private:
            /// Get the task ID.

            /// \todo Is \c id an input and/or output parameter?
            /// \param id The ID to set for this task
            virtual void get_id(std::pair<void*,unsigned short>& id) const {
                PoolTaskInterface::make_id(id, *this);
            }
        }; // class ForEachTask


        /// Apply an operation to a range of iterators.

        /// \tparam rangeT The range of iterators type.
        /// \tparam opT The operation type.
        /// This task creates `for_each` tasks and collects information on the
        /// results of those tasks. Once all tasks have completed it will set
        /// the result future.
        template <typename rangeT, typename opT>
        class ForEachRootTask : public TaskInterface {
        private:
            World& world_; ///< The world where this task will run.
            AtomicInt status_; ///< Accumulates the status of each iteration.
            Future<bool> completion_status_; ///< The result of this set of tasks.

        public:

            /// Constructor.

            /// \param[in] world The world where the tasks will run.
            /// \param[in] range The range of iterators.
            /// \param[in] op The operation that will be applied to the range of iterators.
            ForEachRootTask(World& world, const rangeT range, const opT& op) :
                TaskInterface(0, TaskAttributes::hipri()), world_(world)
            {
                status_ = - (range.size());
                // Create the first for each task.
                world_.taskq.add(new ForEachTask<rangeT,opT>(range, op, *this));
            }

            /// Virtual destructor.
            virtual ~ForEachRootTask() = default;

            /// World accessor.

            /// \return A reference to the \c World.
            World& world() const {
                return world_;
            }

            /// Result accessor.

            /// \return A const reference to the result future.
            const Future<bool>& result() const {
                return completion_status_;
            }

            /// Called by child tasks when they are complete.

            /// \param[in] status The number of iterations that returned true.
            void complete(const int status) {
                status_ += status;
                DependencyInterface::notify();
            }

            /// Task run work.

            /// Sets the result future based on the status of all iterations.
            virtual void run(const TaskThreadEnv&) {
                completion_status_.set(status_ == 0);
            }

        private:
            /// Get the task ID.

            /// \todo Verify that \c id is an output parameter.
            /// \param[out] id The ID to set for this task.
            virtual void get_id(std::pair<void*,unsigned short>& id) const {
                PoolTaskInterface::make_id(id, *this);
            }
        }; // class ForEachRootTask

#endif // HAVE_INTEL_TBB

    }  // namespace detail

} // namespace madness

/// @}

#endif // MADNESS_WORLD_WORLD_TASK_QUEUE_H__INCLUDED
