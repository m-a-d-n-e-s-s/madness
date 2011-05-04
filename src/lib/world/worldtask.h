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

#include <world/worldtaskqueue.h>

namespace madness {


    /// Wrap a callable object and its arguments into a task function

    /// The callable object may have up to 10 arguments
    template <typename fnT, typename arg1_type = void,
            typename arg2_type = void, typename arg3_type = void,
            typename arg4_type = void, typename arg5_type = void,
            typename arg6_type = void, typename arg7_type = void,
            typename arg8_type = void, typename arg9_type = void>
    struct TaskFn : public TaskFnBase<fnT> {
    private:

        typedef TaskFnBase<fnT> TaskFnBase_;    ///< Base class type
        typedef TaskFn<fnT, arg1_type, arg2_type, arg3_type, arg4_type,
            arg5_type, arg6_type, arg7_type, arg8_type, arg9_type> TaskFn_;
                                                ///< This class type

    public:
        typedef typename TaskFnBase_::functionT functionT;  ///< The function type
        typedef typename TaskFnBase_::resultT resultT;      ///< The result type (give by the fuction)
        typedef typename TaskFnBase_::futureT futureT;      ///< Result future type
        typedef typename TaskFnBase_::refT refT;            ///< Remote reference to the result future type
        typedef typename detail::remove_fcvr<arg1_type>::type arg1T;///< Argument 1 type
        typedef typename detail::remove_fcvr<arg2_type>::type arg2T;///< Argument 2 type
        typedef typename detail::remove_fcvr<arg3_type>::type arg3T;///< Argument 3 type
        typedef typename detail::remove_fcvr<arg4_type>::type arg4T;///< Argument 4 type
        typedef typename detail::remove_fcvr<arg5_type>::type arg5T;///< Argument 5 type
        typedef typename detail::remove_fcvr<arg6_type>::type arg6T;///< Argument 6 type
        typedef typename detail::remove_fcvr<arg7_type>::type arg7T;///< Argument 7 type
        typedef typename detail::remove_fcvr<arg8_type>::type arg8T;///< Argument 8 type
        typedef typename detail::remove_fcvr<arg9_type>::type arg9T;///< Argument 9 type

        static const unsigned int arity;    ///< The number of arguments given for the function
                                            ///< \note This may not match the arity of the function
                                            ///< if it has default parameter values

    protected:
        // For convinience
        using TaskFnBase_::func_;
        using TaskFnBase_::result_;

    private:
        // If the value of the argument is known at the time the
        // Note: The type argNT for argN, where N  is > arity should be void

        Future<arg1T> arg1_;///< Argument 1 that will be given to the function
        Future<arg2T> arg2_;///< Argument 2 that will be given to the function
        Future<arg3T> arg3_;///< Argument 3 that will be given to the function
        Future<arg4T> arg4_;///< Argument 4 that will be given to the function
        Future<arg5T> arg5_;///< Argument 5 that will be given to the function
        Future<arg6T> arg6_;///< Argument 6 that will be given to the function
        Future<arg7T> arg7_;///< Argument 7 that will be given to the function
        Future<arg8T> arg8_;///< Argument 8 that will be given to the function
        Future<arg9T> arg9_;///< Argument 9 that will be given to the function

        // These functions are here because we have to differentiate the call
        // based on the number of arguments passed to the function and the
        // return type.

        template <unsigned int N, typename R>
        typename enable_if_c<(! std::is_void<R>::value) && (N == 0u)>::type
        runner() { result_.set(func_()); }

        template <unsigned int N, typename R>
        typename enable_if_c<(! std::is_void<R>::value) && (N == 1u)>::type
        runner() { result_.set(func_(arg1_)); }

        template <unsigned int N, typename R>
        typename enable_if_c<(! std::is_void<R>::value) && (N == 2u)>::type
        runner() { result_.set(func_(arg1_,arg2_)); }

        template <unsigned int N, typename R>
        typename enable_if_c<(! std::is_void<R>::value) && (N == 3u)>::type
        runner() { result_.set(func_(arg1_,arg2_,arg3_)); }

        template <unsigned int N, typename R>
        typename enable_if_c<(! std::is_void<R>::value) && (N == 4u)>::type
        runner() { result_.set(func_(arg1_,arg2_,arg3_,arg4_)); }

        template <unsigned int N, typename R>
        typename enable_if_c<(! std::is_void<R>::value) && (N == 5u)>::type
        runner() { result_.set(func_(arg1_,arg2_,arg3_,arg4_,arg5_)); }

        template <unsigned int N, typename R>
        typename enable_if_c<(! std::is_void<R>::value) && (N == 6u)>::type
        runner() { result_.set(func_(arg1_,arg2_,arg3_,arg4_,arg5_,arg6_)); }

        template <unsigned int N, typename R>
        typename enable_if_c<(! std::is_void<R>::value) && (N == 7u)>::type
        runner() { result_.set(func_(arg1_,arg2_,arg3_,arg4_,arg5_,arg6_,arg7_)); }

        template <unsigned int N, typename R>
        typename enable_if_c<(! std::is_void<R>::value) && (N == 8u)>::type
        runner() { result_.set(func_(arg1_,arg2_,arg3_,arg4_,arg5_,arg6_,arg7_,arg8_)); }

        template <unsigned int N, typename R>
        typename enable_if_c<(! std::is_void<R>::value) && (N == 9u)>::type
        runner() { result_.set(func_(arg1_,arg2_,arg3_,arg4_,arg5_,arg6_,arg7_,arg8_,arg9_)); }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 0u)>::type
        runner() { result_.set(func_()); }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 1u)>::type
        runner() { func_(arg1_); result_.set(); }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 2u)>::type
        runner() { func_(arg1_,arg2_); result_.set(); }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 3u)>::type
        runner() { func_(arg1_,arg2_,arg3_); result_.set(); }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 4u)>::type
        runner() { func_(arg1_,arg2_,arg3_,arg4_); result_.set(); }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 5u)>::type
        runner() { func_(arg1_,arg2_,arg3_,arg4_,arg5_); result_.set(); }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 6u)>::type
        runner() { func_(arg1_,arg2_,arg3_,arg4_,arg5_,arg6_); result_.set(); }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 7u)>::type
        runner() { func_(arg1_,arg2_,arg3_,arg4_,arg5_,arg6_,arg7_); result_.set(); }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 8u)>::type
        runner() { func_(arg1_,arg2_,arg3_,arg4_,arg5_,arg6_,arg7_,arg8_); result_.set(); }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 9u)>::type
        runner() { func_(arg1_,arg2_,arg3_,arg4_,arg5_,arg6_,arg7_,arg8_,arg9_); result_.set(); }

        /// Register non-ready future as a dependency

        /// \tparam T The type of the future to check
        /// \param fut The future to check
        template <typename T>
        inline void check_dependency(Future<T>& fut) {
            if (!fut.probe()) {
                DependencyInterface::inc();
                fut.register_callback(this);
            }
        }

        /// Future<void> is always ready => no op
        inline void check_dependency(Future<void>&) { }

        /// Future<void> is always ready => no op
        inline void check_dependency(Future<Void>&) { }

        /// Check dependencies and register callbacks where necessary
        void check_dependencies() {
            // Check only what we need
            switch( arity ) {
                case 9u:
                    check_dependency(arg9_);
                case 8u:
                    check_dependency(arg8_);
                case 7u:
                    check_dependency(arg7_);
                case 6u:
                    check_dependency(arg6_);
                case 5u:
                    check_dependency(arg5_);
                case 4u:
                    check_dependency(arg4_);
                case 3u:
                    check_dependency(arg3_);
                case 2u:
                    check_dependency(arg2_);
                case 1u:
                    check_dependency(arg1_);
                default: ;
            }
        }

        // am_arg and get_arg differentiates the behavior packing and unpacking
        // AM argument lists between void types and non-void types.

        /// Convert a future into the type needed for an AM arg
        template <typename Arg>
        static const Arg& am_arg(const Future<Arg>& arg) {
            MADNESS_ASSERT(arg.probe());
            return arg;
        }

        /// Convert a void future into the type needed for an AM arg
        static const Future<void>& am_arg(const Future<void>& arg) { return arg; }

        /// Convert a Void future into the type needed for an AM arg
        static const Future<Void>& am_arg(const Future<Void>& arg) { return arg; }

        /// Retrieve a void future from an AM argument list
        template <typename Arg>
        static typename enable_if<std::is_void<Arg>, Future<void> >::type
        get_arg(const AmArg& args) {
            Future<void> arg;
            args & arg;
            return arg;
        }

        /// Retrieve an arg from an AM argument list and stuff it in a future

        /// Get the argument from the AM argument list and stuff it in a future
        /// \param args The AM argument list.
        template <typename Arg>
        static typename disable_if<std::is_void<Arg>, Future<Arg> >::type
        get_arg(const AmArg& args) {
            Arg arg;
            args & arg;
            return Future<Arg>(arg);
        }

        /// Construct a task from an AM argument list

        /// The argument list should be constructed by \c send()
        TaskFn(const AmArg& args) :
            TaskFnBase_(TaskFnBase_::get_task_info(args)),
            arg1_(get_arg<arg1T>(args)),
            arg2_(get_arg<arg2T>(args)),
            arg3_(get_arg<arg3T>(args)),
            arg4_(get_arg<arg4T>(args)),
            arg5_(get_arg<arg5T>(args)),
            arg6_(get_arg<arg6T>(args)),
            arg7_(get_arg<arg7T>(args)),
            arg8_(get_arg<arg8T>(args)),
            arg9_(get_arg<arg9T>(args))
        { check_dependencies(); }

        // Not allowed
        TaskFn(const TaskFn_&);
        TaskFn_ operator=(TaskFn_&);

    public:

        TaskFn(const futureT& result, functionT func, const TaskAttributes& attr) :
            TaskFnBase_(result, func, attr),  arg1_(), arg2_(), arg3_(), arg4_(),
            arg5_(), arg6_(), arg7_(), arg8_(), arg9_()
        { check_dependencies(); }

        template <typename a1T>
        TaskFn(const futureT& result, functionT func, const a1T& a1,
                const TaskAttributes& attr) :
            TaskFnBase_(result, func, attr), arg1_(a1), arg2_(), arg3_(), arg4_(),
            arg5_(), arg6_(), arg7_(), arg8_(), arg9_()
        { check_dependencies(); }

        template <typename a1T, typename a2T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const TaskAttributes& attr = TaskAttributes()) :
            TaskFnBase_(result, func, attr), arg1_(a1), arg2_(a2), arg3_(), arg4_(),
            arg5_(), arg6_(), arg7_(), arg8_(), arg9_()
        { check_dependencies(); }

        template <typename a1T, typename a2T, typename a3T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const TaskAttributes& attr) :
            TaskFnBase_(result, func, attr), arg1_(a1), arg2_(a2), arg3_(a3), arg4_(),
            arg5_(), arg6_(), arg7_(), arg8_(), arg9_()
        { check_dependencies(); }

        template <typename a1T, typename a2T, typename a3T, typename a4T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const TaskAttributes& attr) :
            TaskFnBase_(result, func, attr), arg1_(a1), arg2_(a2), arg3_(a3), arg4_(a4),
            arg5_(), arg6_(), arg7_(), arg8_(), arg9_()
        { check_dependencies(); }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const TaskAttributes& attr) :
            TaskFnBase_(result, func, attr), arg1_(a1), arg2_(a2), arg3_(a3), arg4_(a4),
            arg5_(a5), arg6_(), arg7_(), arg8_(), arg9_()
        { check_dependencies(); }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T,
            typename a6T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6,
                const TaskAttributes& attr) :
            TaskFnBase_(result, func, attr), arg1_(a1), arg2_(a2), arg3_(a3), arg4_(a4),
            arg5_(a5), arg6_(a6), arg7_(), arg8_(), arg9_()
        { check_dependencies(); }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T,
            typename a6T, typename a7T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7,
                const TaskAttributes& attr) :
            TaskFnBase_(result, func, attr), arg1_(a1), arg2_(a2), arg3_(a3), arg4_(a4),
            arg5_(a5), arg6_(a6), arg7_(a7), arg8_(), arg9_()
        { check_dependencies(); }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T,
            typename a6T, typename a7T, typename a8T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7,
                const a8T& a8, const TaskAttributes& attr) :
            TaskFnBase_(result, func, attr), arg1_(a1), arg2_(a2), arg3_(a3), arg4_(a4),
            arg5_(a5), arg6_(a6), arg7_(a7), arg8_(a8), arg9_()
        { check_dependencies(); }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T,
            typename a6T, typename a7T, typename a8T, typename a9T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7,
                const a8T& a8, const a9T& a9, const TaskAttributes& attr) :
            TaskFnBase_(result, func, attr), arg1_(a1), arg2_(a2), arg3_(a3), arg4_(a4),
            arg5_(a5), arg6_(a6), arg7_(a7), arg8_(a8), arg9_(a9)
        { check_dependencies(); }

        virtual ~TaskFn() { }

        void run(World&) {
          runner<arity, resultT>();
        }

        /// Send this task to \c dest and submit it to the task queue there.

        /// \param dest The process where the task should be sent to run
        /// \throw MadnessException If the process is not stealable
        /// \note If this task is in the task queue, it must be deleted and set
        /// to NULL or the task will be run twice and the result future will be
        /// set twice.
        virtual void send(ProcessID dest) const {
            MADNESS_ASSERT(TaskAttributes::is_stealable());
            TaskInterface::get_world()->am.send(dest, TaskFn_::handler,
                new_am_arg(TaskFnBase_::task_info(), am_arg(arg1_), am_arg(arg2_),
                am_arg(arg3_), am_arg(arg4_), am_arg(arg5_), am_arg(arg6_),
                am_arg(arg7_), am_arg(arg8_), am_arg(arg9_)));
        }

        /// Construct and add b new task function to the task queue

        /// This function is used for spawning tasks on remote nodes. It should
        /// not be used directly. If you need to spawn a remote task, use
        /// \c remote_task() function.
        /// \param arg The AM arguments
        static void handler(const AmArg& args) {
            args.get_world()->taskq.add(new TaskFn_(args));
        }
    };


    // Task function factory functions

    template <typename fnT>
    TaskFnBase<fnT>* make_task(fnT fn, const TaskAttributes& attr = TaskAttributes()) {
        typedef TaskFn<fnT> task_type;
        typename task_type::futureT r;
        task_type* task = new task_type(r, fn);
        return task;
    }

    template <typename fnT, typename a1T>
    TaskFnBase<fnT>* make_task(fnT fn, const a1T& a1,
            const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T> task_type;
        typename task_type::futureT r;
        task_type* task = new task_type(r, fn, a1);
        return task;
    }

    template <typename fnT, typename a1T, typename a2T>
    TaskFnBase<fnT>* make_task(fnT fn, const a1T& a1, const a2T& a2,
            const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T> task_type;
        typename task_type::futureT r;
        task_type* task = new task_type(r, fn, a1, a2, attr);
        return task;
    }

    template <typename fnT, typename a1T, typename a2T, typename a3T>
    TaskFnBase<fnT>* make_task(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3,
            const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T, a3T> task_type;
        typename task_type::futureT r;
        task_type* task = new task_type(r, fn, a1, a2, a3, attr);
        return task;
    }

    template <typename fnT,typename a1T, typename a2T,
        typename a3T, typename a4T>
    TaskFnBase<fnT>* make_task(fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T> task_type;
        typename task_type::futureT r;
        task_type* task = new task_type(r, fn, a1, a2, a3, a4, attr);
        return task;
    }

    template <typename fnT, typename a1T, typename a2T,
        typename a3T, typename a4T, typename a5T>
    TaskFnBase<fnT>* make_task(fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const a5T& a5, const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T> task_type;
        typename task_type::futureT r;
        task_type* task = new task_type(r, fn, a1, a2, a3, a4, attr);
        return task;
    }

    template <typename fnT, typename a1T, typename a2T, typename a3T,
        typename a4T, typename a5T, typename a6T>
    TaskFnBase<fnT>* make_task(fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6,
            const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<typename result_of<fnT>::type, fnT, a1T, a2T, a3T, a4T, a5T, a6T> task_type;
        typename task_type::futureT r;
        task_type* task = new task_type(r, fn, a1, a2, a3, a4, a5, a6, attr);
        return task;
    }

    template <typename fnT, typename a0T, typename a1T, typename a2T, typename a3T,
        typename a4T, typename a5T, typename a6T, typename a7T>
    TaskFnBase<fnT>* make_task(fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7,
            const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T> task_type;
        typename task_type::futureT r;
        task_type* task = new task_type(r, fn, a1, a2, a3, a4, a5, a6, a7, attr);
        return task;
    }

    template <typename fnT, typename a0T, typename a1T, typename a2T, typename a3T,
    typename a4T, typename a5T, typename a6T, typename a7T, typename a8T>
    TaskFnBase<fnT>* make_task(fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7,
            const a8T& a8, const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T> task_type;
        typename task_type::futureT r;
        task_type* task = new task_type(r, fn, a1, a2, a3, a4, a5, a6, a7, a8, attr);
        return task;
    }

    template <typename fnT, typename a0T, typename a1T, typename a2T, typename a3T,
        typename a4T, typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
    TaskFnBase<fnT>* make_task(fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7,
            const a8T& a8, const a9T& a9, const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T> task_type;
        typename task_type::futureT r;
        task_type* task = new task_type(r, fn, a1, a2, a3, a4, a5, a6, a7, a8, a9, attr);
        return task;
    }

    /// Send the task to \c dest and submit it to the task queue there.

    /// \param world The world that will handle the task
    /// \param dest The process where the task should be sent to run
    /// \param fn The function to be run by the task
    /// \return A future that will be set by the \c dest process when the task
    /// is complete.
    template <typename fnT>
    typename TaskFn<fnT>::functionT
    remote_task(World& world, ProcessID dest, fnT fn, const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT> taskfnT;
        typedef typename taskfnT::futureT futureT;
        typedef typename taskfnT::infoT infotT;
        typedef Future<void> voidT;

        futureT result;
        world.am.send(dest, & taskfnT::handler, new_am_arg(
            infotT(result.remote_ref(), fn, attr), voidT(), voidT(), voidT(),
            voidT(), voidT(), voidT(), voidT(), voidT(), voidT()));

        return result;
    }

    /// Send the task to \c dest and submit it to the task queue there.

    /// \param world The world that will handle the task
    /// \param dest The process where the task should be sent to run
    /// \param fn The function to be run by the task
    /// \param a1 Task argument 1
    /// \return A future that will be set by the \c dest process when the task
    /// is complete.
    template <typename fnT, typename a1T>
    typename TaskFn<fnT, a1T>::functionT
    remote_task(World& world, ProcessID dest, fnT fn, const a1T& a1,
            const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T> taskfnT;
        typedef typename taskfnT::futureT futureT;
        typedef typename taskfnT::infoT infotT;
        typedef Future<void> voidT;

        futureT result;
        world.am.send(dest, & taskfnT::handler, new_am_arg(
            infotT(result.remote_ref(), fn, attr), a1, voidT(), voidT(), voidT(),
            voidT(), voidT(), voidT(), voidT(), voidT()));

        return result;
    }

    /// Send the task to \c dest and submit it to the task queue there.

    /// \param world The world that will handle the task
    /// \param dest The process where the task should be sent to run
    /// \param fn The function to be run by the task
    /// \param a1 Task argument 1
    /// \param a2 Task argument 2
    /// \return A future that will be set by the \c dest process when the task
    /// is complete.
    template <typename fnT, typename a1T, typename a2T>
    typename TaskFn<fnT, a1T, a2T>::functionT
    remote_task(World& world, ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
            const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T> taskfnT;
        typedef typename taskfnT::futureT futureT;
        typedef typename taskfnT::infoT infotT;
        typedef Future<void> voidT;

        futureT result;
        world.am.send(dest, & taskfnT::handler, new_am_arg(
            infotT(result.remote_ref(), fn, attr), a1, a2, voidT(), voidT(),
            voidT(), voidT(), voidT(), voidT(), voidT()));

        return result;
    }

    /// Send the task to \c dest and submit it to the task queue there.

    /// \param world The world that will handle the task
    /// \param dest The process where the task should be sent to run
    /// \param fn The function to be run by the task
    /// \param a1 Task argument 1
    /// \param a2 Task argument 2
    /// \param a3 Task argument 3
    /// \return A future that will be set by the \c dest process when the task
    /// is complete.
    template <typename fnT, typename a1T, typename a2T, typename a3T>
    typename TaskFn<fnT, a1T, a2T, a3T>::functionT
    remote_task(World& world, ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T, a3T> taskfnT;
        typedef typename taskfnT::futureT futureT;
        typedef typename taskfnT::infoT infotT;
        typedef Future<void> voidT;

        futureT result;
        world.am.send(dest, & taskfnT::handler, new_am_arg(
            infotT(result.remote_ref(), fn, attr), a1, a2, a3, voidT(), voidT(),
            voidT(), voidT(), voidT(), voidT()));

        return result;
    }

    /// Send the task to \c dest and submit it to the task queue there.

    /// \param world The world that will handle the task
    /// \param dest The process where the task should be sent to run
    /// \param fn The function to be run by the task
    /// \param a1 Task argument 1
    /// \param a2 Task argument 2
    /// \param a3 Task argument 3
    /// \param a4 Task argument 4
    /// \return A future that will be set by the \c dest process when the task
    /// is complete.
    template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T>
    typename TaskFn<fnT, a1T, a2T, a3T, a4T>::functionT
    remote_task(World& world, ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T> taskfnT;
        typedef typename taskfnT::futureT futureT;
        typedef typename taskfnT::infoT infotT;
        typedef Future<void> voidT;

        typename taskfnT::futureT result;
        world.am.send(dest, & taskfnT::handler, new_am_arg(
            infotT(result.remote_ref(), fn, attr), a1, a2, a3, a4, voidT(),
            voidT(), voidT(), voidT(), voidT()));

        return result;
    }

    /// Send the task to \c dest and submit it to the task queue there.

    /// \param world The world that will handle the task
    /// \param dest The process where the task should be sent to run
    /// \param fn The function to be run by the task
    /// \param a1 Task argument 1
    /// \param a2 Task argument 2
    /// \param a3 Task argument 3
    /// \param a4 Task argument 4
    /// \param a5 Task argument 5
    /// \return A future that will be set by the \c dest process when the task
    /// is complete.
    template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
        typename a5T>
    typename TaskFn<fnT, a1T, a2T, a3T, a4T, a5T>::functionT
    remote_task(World& world, ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const a5T& a5,
            const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T> taskfnT;
        typedef typename taskfnT::futureT futureT;
        typedef typename taskfnT::infoT infotT;
        typedef Future<void> voidT;

        typename taskfnT::futureT result;
        world.am.send(dest, & taskfnT::handler, new_am_arg(
            infotT(result.remote_ref(), fn, attr), a1, a2, a3, a4, a5, voidT(),
            voidT(), voidT(), voidT()));

        return result;
    }

    /// Send the task to \c dest and submit it to the task queue there.

    /// \param world The world that will handle the task
    /// \param dest The process where the task should be sent to run
    /// \param fn The function to be run by the task
    /// \param a1 Task argument 1
    /// \param a2 Task argument 2
    /// \param a3 Task argument 3
    /// \param a4 Task argument 4
    /// \param a5 Task argument 5
    /// \param a6 Task argument 6
    /// \return A future that will be set by the \c dest process when the task
    /// is complete.
    template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
        typename a5T, typename a6T>
    typename TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T>::functionT
    remote_task(World& world, ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6,
            const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T> taskfnT;
        typedef typename taskfnT::futureT futureT;
        typedef typename taskfnT::infoT infotT;
        typedef Future<void> voidT;

        typename taskfnT::futureT result;
        world.am.send(dest, & taskfnT::handler, new_am_arg(
            infotT(result.remote_ref(), fn, attr), a1, a2, a3, a4, a5, a6,
            voidT(), voidT(), voidT()));

        return result;
    }

    /// Send the task to \c dest and submit it to the task queue there.

    /// \param world The world that will handle the task
    /// \param dest The process where the task should be sent to run
    /// \param fn The function to be run by the task
    /// \param a1 Task argument 1
    /// \param a2 Task argument 2
    /// \param a3 Task argument 3
    /// \param a4 Task argument 4
    /// \param a5 Task argument 5
    /// \param a6 Task argument 6
    /// \param a7 Task argument 7
    /// \return A future that will be set by the \c dest process when the task
    /// is complete.
    template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
        typename a5T, typename a6T, typename a7T>
    typename TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T>::functionT
    remote_task(World& world, ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7,
            const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T> taskfnT;
        typedef typename taskfnT::futureT futureT;
        typedef typename taskfnT::infoT infotT;
        typedef Future<void> voidT;

        typename taskfnT::futureT result;
        world.am.send(dest, & taskfnT::handler, new_am_arg(
            infotT(result.remote_ref(), fn, attr), a1, a2, a3, a4, a5, a6, a7,
            voidT(), voidT()));

        return result;
    }

    /// Send the task to \c dest and submit it to the task queue there.

    /// \param world The world that will handle the task
    /// \param dest The process where the task should be sent to run
    /// \param fn The function to be run by the task
    /// \param a1 Task argument 1
    /// \param a2 Task argument 2
    /// \param a3 Task argument 3
    /// \param a4 Task argument 4
    /// \param a5 Task argument 5
    /// \param a6 Task argument 6
    /// \param a7 Task argument 7
    /// \param a8 Task argument 8
    /// \return A future that will be set by the \c dest process when the task
    /// is complete.
    template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
        typename a5T, typename a6T, typename a7T, typename a8T>
    typename TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T>::functionT
    remote_task(World& world, ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7,
            const a8T& a8, const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T> taskfnT;
        typedef typename taskfnT::futureT futureT;
        typedef typename taskfnT::infoT infotT;
        typedef Future<void> voidT;

        typename taskfnT::futureT result;
        world.am.send(dest, & taskfnT::handler, new_am_arg(
            infotT(result.remote_ref(), fn, attr), a1, a2, a3, a4, a5, a6, a7,
            a8, voidT()));

        return result;
    }

    /// Send the task to \c dest and submit it to the task queue there.

    /// \param world The world that will handle the task
    /// \param dest The process where the task should be sent to run
    /// \param fn The function to be run by the task
    /// \param a1 Task argument 1
    /// \param a2 Task argument 2
    /// \param a3 Task argument 3
    /// \param a4 Task argument 4
    /// \param a5 Task argument 5
    /// \param a6 Task argument 6
    /// \param a7 Task argument 7
    /// \param a8 Task argument 8
    /// \param a9 Task argument 9
    /// \return A future that will be set by the \c dest process when the task
    /// is complete.
    template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
        typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
    typename TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T>::functionT
    remote_task(World& world, ProcessID dest, fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7,
            const a8T& a8, const a9T& a9, const TaskAttributes& attr = TaskAttributes())
    {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T> taskfnT;

        typename taskfnT::futureT result;
        world.am.send(dest, & taskfnT::handler, new_am_arg(
            typename taskfnT::infoT(result.remote_ref(), fn, attr), a1, a2, a3,
            a4, a5, a6, a7, a8, a9));

        return result;
    }

    namespace {

        template <typename T>
        struct ArgCountHelper {
            static const unsigned int value = 1;
        };

        template <>
        struct ArgCountHelper<void> {
            static const unsigned int value = 0;
        };

        // Counts the number of arguments that will be given to a task function
        template <typename a1T, typename a2T, typename a3T, typename a4T,
            typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
        struct ArgCount : public std::integral_constant<unsigned int,
                ArgCountHelper<a1T>::value + ArgCountHelper<a2T>::value +
                ArgCountHelper<a3T>::value + ArgCountHelper<a4T>::value +
                ArgCountHelper<a5T>::value + ArgCountHelper<a6T>::value +
                ArgCountHelper<a7T>::value + ArgCountHelper<a8T>::value +
                ArgCountHelper<a9T>::value>
        { };

    }  // namespace

    // static member variables for TaskFn
    template <typename fnT, typename a1T, typename a2T, typename a3T,
            typename a4T, typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
    const unsigned int TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T>::arity =
            ArgCount<a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T>::value;

} // namespace madness


#endif // MADNESS_WORLD_WORLDTASK_H__INCLUDED
