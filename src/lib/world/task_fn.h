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


#ifndef MADNESS_WORLD_TASK_FN_H__INCLUDED
#define MADNESS_WORLD_TASK_FN_H__INCLUDED

/// \file task_fn.h
/// \brief Defines task functions and their interface.

// To do:
// a) Redo this ontop of serializable tasks which will remote much of the clutter
//    due to multiple length argument lists.
// b) Stealing which pretty much presume a) has been done

#include <iostream>
#include <world/nodefaults.h>
#include <world/worlddep.h>
#include <world/worldthread.h>
#include <world/worldfut.h>
#include <world/typestuff.h>

namespace madness {

    // Forward decls
    class World;
    class WorldTaskQueue;
    class TaskInterface;
    template <typename> struct TaskFnBase;

    template <typename fnT>
    TaskFnBase<fnT>* make_task(fnT, const TaskAttributes& = TaskAttributes());

    template <typename fnT, typename a0T>
    TaskFnBase<fnT>* make_task(fnT, const a0T&, const TaskAttributes& = TaskAttributes());

    template <typename fnT, typename a0T, typename a1T>
    TaskFnBase<fnT>* make_task(fnT, const a0T&, const a1T&, const TaskAttributes& = TaskAttributes());

    template <typename fnT, typename a0T, typename a1T, typename a2T>
    TaskFnBase<fnT>* make_task(fnT, const a0T&, const a1T&, const a2T&, const TaskAttributes& = TaskAttributes());

    template <typename fnT, typename a0T, typename a1T, typename a2T,
        typename a3T>
    TaskFnBase<fnT>* make_task(fnT, const a0T&, const a1T&, const a2T&,
            const a3T&, const TaskAttributes& = TaskAttributes());

    template <typename fnT, typename a0T, typename a1T, typename a2T,
        typename a3T, typename a4T>
    TaskFnBase<fnT>* make_task(fnT, const a0T&, const a1T&, const a2T&,
            const a3T&, const a4T&, const TaskAttributes& = TaskAttributes());

    template <typename fnT, typename a0T, typename a1T, typename a2T, typename a3T,
        typename a4T, typename a5T>
    TaskFnBase<fnT>* make_task(fnT, const a0T&, const a1T&, const a2T&,
            const a3T&, const a4T&, const a5T&, const TaskAttributes& = TaskAttributes());

    template <typename fnT, typename a0T, typename a1T, typename a2T, typename a3T,
        typename a4T, typename a5T, typename a6T>
    TaskFnBase<fnT>* make_task(fnT, const a0T&, const a1T&, const a2T&,
            const a3T&, const a4T&, const a5T&, const a6T&, const TaskAttributes& = TaskAttributes());

    template <typename fnT, typename a0T, typename a1T, typename a2T, typename a3T,
    typename a4T, typename a5T, typename a6T, typename a7T>
    TaskFnBase<fnT>* make_task(fnT, const a0T&, const a1T&, const a2T&,
            const a3T&, const a4T&, const a5T&, const a6T&, const a7T&, const TaskAttributes& = TaskAttributes());

    template <typename fnT, typename a0T, typename a1T, typename a2T, typename a3T,
        typename a4T, typename a5T, typename a6T, typename a7T, typename a8T>
    TaskFnBase<fnT>* make_task(fnT, const a0T&, const a1T&, const a2T&,
            const a3T&, const a4T&, const a5T&, const a6T&, const a7T&, const a8T&, const TaskAttributes& = TaskAttributes());

    template <typename fnT, typename a0T, typename a1T, typename a2T, typename a3T,
        typename a4T, typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
    TaskFnBase<fnT>* make_task(fnT, const a0T&, const a1T&, const a2T&,
            const a3T&, const a4T&, const a5T&, const a6T&, const a7T&, const a8T&, const a9T&, const TaskAttributes& = TaskAttributes());

    /// All world tasks must be derived from this public interface

    /// Multiple worlds with independent queues feed tasks into shared task
    /// pool that is mapped to the H/W.
    ///
    /// For simplicity and backward compatibility we maintain two run interfaces
    /// but new code should adopt the multithreaded interface
    ///
    /// \c run(World&) - the user implements this for a single-threaded task
    ///
    /// \c run(World&, \c const \c TaskThreadEnv&) - the user implements this for
    /// a multi-threaded task.
    ///
    class TaskInterface : public DependencyInterface , public PoolTaskInterface {
        friend class WorldTaskQueue;
    private:
        volatile World* world;
        CallbackInterface* completion;

        // Used for submission to underlying queue when all dependencies are satisfied
        struct Submit : public CallbackInterface {
            PoolTaskInterface* p;
            Submit(PoolTaskInterface* p) : p(p) {}
            void notify() {
                ThreadPool::add(p);
            }
        } submit;


        /// Set task info

        /// \param w The world object that contains the task
        /// \param c Call this callback on completion
        void set_info(World* w, CallbackInterface* c) {
            world = w;
            completion = c;
        }

        /// Adds call back to schedule task when outstanding dependencies are satisfied
        void register_submit_callback() { register_callback(&submit); }

    protected:
        virtual void run(const TaskThreadEnv& env);

    public:
        static bool debug;

        /// Create a new task with ndepend dependencies (default 0) and given attributes
        TaskInterface(int ndepend=0, const TaskAttributes attr = TaskAttributes())
                : DependencyInterface(ndepend)
                , PoolTaskInterface(attr)
                , world(0)
                , completion(0)
                , submit(this)
        {}

        /// Create a new task with zero dependencies and given attributes
        explicit TaskInterface(const TaskAttributes& attr)
                : DependencyInterface(0)
                , PoolTaskInterface(attr)
                , world(0)
                , completion(0)
                , submit(this)
        {}

//         void serialize(Buffer& ar) {
//             throw "there is no way this is correct";
//             ar & *static_cast<PoolTaskInterface*>(this) & world;
//         }

        /// Runs a single-threaded task ... derived classes must implement this.

        /// This interface may disappear so new code should use the multi-threaded interface.
        virtual void run(World& /*world*/);

        /// Runs a multi-threaded task
        virtual void run(World& world, const TaskThreadEnv& env) {
            //print("in virtual run(world,env) method", env.nthread(), env.id());
            if (env.nthread() != 1)
                MADNESS_EXCEPTION("World TaskInterface: user did not implement run(world, taskthreadenv) for multithreaded task", 0);
            run(world);
        }

        World* get_world() const { return const_cast<World*>(world); }

        virtual ~TaskInterface() { if (completion) completion->notify(); }

    };



    // Internal: Convenience for serializing
    template <typename refT, typename functionT>
    struct TaskHandlerInfo {
        refT ref;
        functionT func;
        TaskAttributes attr;
        TaskHandlerInfo(const refT& ref, functionT func, const TaskAttributes& attr)
                : ref(ref), func(func),attr(attr) {}
        TaskHandlerInfo() {}

        template <typename Archive>
        void serialize(const Archive& ar) {
            serialize_internal<functionT>(ar);
        }

    private:
        template <typename fnT, typename Archive>
        typename enable_if_c<function_traits<fnT>::value || memfunc_traits<fnT>::value>::type
        serialize_internal(const Archive& ar) {
            ar & ref & archive::wrap_opaque(func) & attr;
        }

        template <typename fnT, typename Archive>
        typename disable_if_c<function_traits<fnT>::value || memfunc_traits<fnT>::value>::type
        serialize_internal(const Archive& ar) {
            ar & ref & func & attr;
        }
    };

    namespace {

        // Remove Future, const, volatile, and reference qualifiers from the type
        template <typename T>
        struct remove_fcvr{
            typedef typename remove_future< typename std::remove_cv<
                    typename std::remove_reference<T>::type >::type >::type type;
        };

    }  // namespace

    template <typename fnT>
    struct TaskFnBase : public TaskInterface {
    public:

        typedef fnT functionT;
        typedef typename remove_future<typename result_of<fnT>::type>::type resultT;
        typedef Future<resultT> futureT;
        typedef typename futureT::remote_refT refT;
        typedef TaskHandlerInfo<refT,functionT> infoT;

    protected:
        futureT result_;
        const functionT func_;


        static infoT get_task_info(const AmArg& args) {
            infoT info;
            args & info;
            return info;
        }

    public:
        TaskFnBase(const futureT& result, functionT func, const TaskAttributes& attr) :
            TaskInterface(attr), result_(result), func_(func)
        { }

        TaskFnBase(const infoT& info) :
            TaskInterface(info.attr), result_(info.ref), func_(info.func)
        { }

        infoT task_info() const {
            return infoT(result_.remote_ref(* TaskInterface::get_world()), func_, *this);
        }

        const futureT& result() const { return result_; }

        virtual ~TaskFnBase() { }
    };

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
        typedef typename remove_fcvr<arg1_type>::type arg1T;///< Argument 1 type
        typedef typename remove_fcvr<arg2_type>::type arg2T;///< Argument 2 type
        typedef typename remove_fcvr<arg3_type>::type arg3T;///< Argument 3 type
        typedef typename remove_fcvr<arg4_type>::type arg4T;///< Argument 4 type
        typedef typename remove_fcvr<arg5_type>::type arg5T;///< Argument 5 type
        typedef typename remove_fcvr<arg6_type>::type arg6T;///< Argument 6 type
        typedef typename remove_fcvr<arg7_type>::type arg7T;///< Argument 7 type
        typedef typename remove_fcvr<arg8_type>::type arg8T;///< Argument 8 type
        typedef typename remove_fcvr<arg9_type>::type arg9T;///< Argument 9 type

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

        Future<arg1T> arg1; ///< Argument 1 that will be given to the function
        Future<arg2T> arg2; ///< Argument 2 that will be given to the function
        Future<arg3T> arg3; ///< Argument 3 that will be given to the function
        Future<arg4T> arg4; ///< Argument 4 that will be given to the function
        Future<arg5T> arg5; ///< Argument 5 that will be given to the function
        Future<arg6T> arg6; ///< Argument 6 that will be given to the function
        Future<arg7T> arg7; ///< Argument 7 that will be given to the function
        Future<arg8T> arg8; ///< Argument 8 that will be given to the function
        Future<arg9T> arg9; ///< Argument 9 that will be given to the function

        // These functions are here because we have to differentiate the call
        // based on the number of arguments passed to the function and the
        // return type.

        template <unsigned int N, typename R>
        typename enable_if_c<(! std::is_void<R>::value) && (N == 0u)>::type
        runner() { result_.set(func_()); }

        template <unsigned int N, typename R>
        typename enable_if_c<(! std::is_void<R>::value) && (N == 1u)>::type
        runner() { result_.set(func_(arg1)); }

        template <unsigned int N, typename R>
        typename enable_if_c<(! std::is_void<R>::value) && (N == 2u)>::type
        runner() { result_.set(func_(arg1,arg2)); }

        template <unsigned int N, typename R>
        typename enable_if_c<(! std::is_void<R>::value) && (N == 3u)>::type
        runner() { result_.set(func_(arg1,arg2,arg3)); }

        template <unsigned int N, typename R>
        typename enable_if_c<(! std::is_void<R>::value) && (N == 4u)>::type
        runner() { result_.set(func_(arg1,arg2,arg3,arg4)); }

        template <unsigned int N, typename R>
        typename enable_if_c<(! std::is_void<R>::value) && (N == 5u)>::type
        runner() { result_.set(func_(arg1,arg2,arg3,arg4,arg5)); }

        template <unsigned int N, typename R>
        typename enable_if_c<(! std::is_void<R>::value) && (N == 6u)>::type
        runner() { result_.set(func_(arg1,arg2,arg3,arg4,arg5,arg6)); }

        template <unsigned int N, typename R>
        typename enable_if_c<(! std::is_void<R>::value) && (N == 7u)>::type
        runner() { result_.set(func_(arg1,arg2,arg3,arg4,arg5,arg6,arg7)); }

        template <unsigned int N, typename R>
        typename enable_if_c<(! std::is_void<R>::value) && (N == 8u)>::type
        runner() { result_.set(func_(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8)); }

        template <unsigned int N, typename R>
        typename enable_if_c<(! std::is_void<R>::value) && (N == 9u)>::type
        runner() { result_.set(func_(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)); }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 0u)>::type
        runner() { result_.set(func_()); }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 1u)>::type
        runner() { func_(arg1); result_.set(); }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 2u)>::type
        runner() { func_(arg1,arg2); result_.set(); }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 3u)>::type
        runner() { func_(arg1,arg2,arg3); result_.set(); }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 4u)>::type
        runner() { func_(arg1,arg2,arg3,arg4); result_.set(); }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 5u)>::type
        runner() { func_(arg1,arg2,arg3,arg4,arg5); result_.set(); }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 6u)>::type
        runner() { func_(arg1,arg2,arg3,arg4,arg5,arg6); result_.set(); }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 7u)>::type
        runner() { func_(arg1,arg2,arg3,arg4,arg5,arg6,arg7); result_.set(); }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 8u)>::type
        runner() { func_(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8); result_.set(); }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 9u)>::type
        runner() { func_(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9); result_.set(); }

        /// Send the task to \c dest and submit it to the task queue there.

        /// \param dest The process where the task should be sent to run
        /// \throw MadnessException If the process is not stealable
        /// \note If this task is in the task queue, it must be deleted and set
        /// to null the or task will be run twice and the result future will be
        /// set twice.
        virtual void send(ProcessID dest) const {
            MADNESS_ASSERT(TaskAttributes::is_stealable());
            TaskInterface::get_world()->am.send(dest, TaskFn_::handler,
                new_am_arg(TaskFnBase_::task_info(), am_arg(arg1), am_arg(arg2),
                am_arg(arg3), am_arg(arg4), am_arg(arg5), am_arg(arg6),
                am_arg(arg7), am_arg(arg8), am_arg(arg9)));
        }

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
            switch( arity ) {
                case 9u:
                    check_dependency(arg9);
                case 8u:
                    check_dependency(arg8);
                case 7u:
                    check_dependency(arg7);
                case 6u:
                    check_dependency(arg6);
                case 5u:
                    check_dependency(arg5);
                case 4u:
                    check_dependency(arg4);
                case 3u:
                    check_dependency(arg3);
                case 2u:
                    check_dependency(arg2);
                case 1u:
                    check_dependency(arg1);
                default: ;
            }
        }

        /// Convert a future into the type needed for an AM arg
        template <typename Arg>
        static const Arg& am_arg(const Future<Arg>& arg) {
            MADNESS_ASSERT(arg.probe());
            return arg;
        }

        /// Convert a void future into the type needed for an AM arg
        static const Future<void>& am_arg(const Future<void>& arg) {
            return arg;
        }

        /// Convert a Void future into the type needed for an AM arg
        static const Future<Void>& am_arg(const Future<Void>& arg) {
            return arg;
        }

        // The function definition is in worldtask.h because it needs the definition
        // of WorldTaskQueue to be inistantiated.
        static void handler(const AmArg& arg);
//        {
//            arg.get_world()->taskq.add(new TaskFn_(arg));
//        }

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
            arg1(get_arg<arg1T>(args)),
            arg2(get_arg<arg2T>(args)),
            arg3(get_arg<arg3T>(args)),
            arg4(get_arg<arg4T>(args)),
            arg5(get_arg<arg5T>(args)),
            arg6(get_arg<arg6T>(args)),
            arg7(get_arg<arg7T>(args)),
            arg8(get_arg<arg8T>(args)),
            arg9(get_arg<arg9T>(args))
        { check_dependencies(); }

        // Not allowed
        TaskFn(const TaskFn_&);
        TaskFn_ operator=(TaskFn_&);

    public:

        TaskFn(const futureT& result, functionT func, const TaskAttributes& attr) :
            TaskFnBase_(result, func, attr),  arg1(), arg2(), arg3(), arg4(), arg5(), arg6(), arg7(), arg8(), arg9()
        { check_dependencies(); }

        template <typename a1T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const TaskAttributes& attr) :
            TaskFnBase_(result, func, attr), arg1(a1), arg2(), arg3(), arg4(), arg5(), arg6(), arg7(), arg8(), arg9()
        { check_dependencies(); }

        template <typename a1T, typename a2T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const TaskAttributes& attr) :
            TaskFnBase_(result, func, attr), arg1(a1), arg2(a2), arg3(), arg4(), arg5(), arg6(), arg7(), arg8(), arg9()
        { check_dependencies(); }

        template <typename a1T, typename a2T, typename a3T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const a3T& a3, const TaskAttributes& attr) :
            TaskFnBase_(result, func, attr), arg1(a1), arg2(a2), arg3(a3), arg4(), arg5(), arg6(), arg7(), arg8(), arg9()
        { check_dependencies(); }

        template <typename a1T, typename a2T, typename a3T, typename a4T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const TaskAttributes& attr) :
            TaskFnBase_(result, func, attr), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(), arg6(), arg7(), arg8(), arg9()
        { check_dependencies(); }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
                const a5T& a5, const TaskAttributes& attr) :
            TaskFnBase_(result, func, attr), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5), arg6(), arg7(), arg8(), arg9()
        { check_dependencies(); }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
                const a5T& a5, const a6T& a6, const TaskAttributes& attr) :
            TaskFnBase_(result, func, attr), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5), arg6(a6), arg7(), arg8(), arg9()
        { check_dependencies(); }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T, typename a7T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
                const a5T& a5, const a6T& a6, const a7T& a7, const TaskAttributes& attr) :
            TaskFnBase_(result, func, attr), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5), arg6(a6), arg7(a7), arg8(), arg9()
        { check_dependencies(); }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T, typename a7T, typename a8T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
                const a5T& a5, const a6T& a6, const a7T& a7, const a8T& a8, const TaskAttributes& attr) :
            TaskFnBase_(result, func, attr), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5), arg6(a6), arg7(a7), arg8(a8), arg9()
        { check_dependencies(); }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4,
                const a5T& a5, const a6T& a6, const a7T& a7, const a8T& a8, const a9T& a9, const TaskAttributes& attr) :
            TaskFnBase_(result, func, attr), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5), arg6(a6), arg7(a7), arg8(a8), arg9(a9)
        { check_dependencies(); }

        virtual ~TaskFn() { }

        void run(World&) {
          runner<arity, resultT>();
        }
    };

    // Task function factory functions

    template <typename fnT>
    TaskFnBase<fnT>* make_task(fnT fn, const TaskAttributes& attr) {
        typedef TaskFn<fnT> task_type;
        typename task_type::futureT r;
        task_type* task = new task_type(r, fn);
        return task;
    }

    template <typename fnT, typename a1T>
    TaskFnBase<fnT>* make_task(fnT fn, const a1T& a1, const TaskAttributes& attr) {
        typedef TaskFn<fnT, a1T> task_type;
        typename task_type::futureT r;
        task_type* task = new task_type(r, fn, a1);
        return task;
    }

    template <typename fnT, typename a1T, typename a2T>
    TaskFnBase<fnT>* make_task(fnT fn, const a1T& a1, const a2T& a2, const TaskAttributes& attr) {
        typedef TaskFn<fnT, a1T, a2T> task_type;
        typename task_type::futureT r;
        task_type* task = new task_type(r, fn, a1, a2, attr);
        return task;
    }

    template <typename fnT, typename a1T, typename a2T, typename a3T>
    TaskFnBase<fnT>* make_task(fnT fn, const a1T& a1, const a2T& a2, const a3T& a3, const TaskAttributes& attr) {
        typedef TaskFn<typename result_of<fnT>::type, fnT, a1T, a2T, a3T> task_type;
        typename task_type::futureT r;
        task_type* task = new task_type(r, fn, a1, a2, a3, attr);
        return task;
    }

    template <typename fnT,typename a1T, typename a2T,
        typename a3T, typename a4T>
    TaskFnBase<fnT>* make_task(fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const TaskAttributes& attr) {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T> task_type;
        typename task_type::futureT r;
        task_type* task = new task_type(r, fn, a1, a2, a3, a4, attr);
        return task;
    }

    template <typename fnT, typename a1T, typename a2T,
        typename a3T, typename a4T, typename a5T>
    TaskFnBase<fnT>* make_task(fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const a5T& a5, const TaskAttributes& attr) {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T> task_type;
        typename task_type::futureT r;
        task_type* task = new task_type(r, fn, a1, a2, a3, a4, attr);
        return task;
    }

    template <typename fnT, typename a1T, typename a2T, typename a3T,
        typename a4T, typename a5T, typename a6T>
    TaskFnBase<fnT>* make_task(fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const TaskAttributes& attr) {
        typedef TaskFn<typename result_of<fnT>::type, fnT, a1T, a2T, a3T, a4T, a5T, a6T> task_type;
        typename task_type::futureT r;
        task_type* task = new task_type(r, fn, a1, a2, a3, a4, a5, a6, attr);
        return task;
    }

    template <typename fnT, typename a0T, typename a1T, typename a2T, typename a3T,
        typename a4T, typename a5T, typename a6T, typename a7T>
    TaskFnBase<fnT>* make_task(fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7, const TaskAttributes& attr) {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T> task_type;
        typename task_type::futureT r;
        task_type* task = new task_type(r, fn, a1, a2, a3, a4, a5, a6, a7, attr);
        return task;
    }

    template <typename fnT, typename a0T, typename a1T, typename a2T, typename a3T,
    typename a4T, typename a5T, typename a6T, typename a7T, typename a8T>
    TaskFnBase<fnT>* make_task(fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7, const a8T& a8, const TaskAttributes& attr) {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T> task_type;
        typename task_type::futureT r;
        task_type* task = new task_type(r, fn, a1, a2, a3, a4, a5, a6, a7, a8, attr);
        return task;
    }

    template <typename fnT, typename a0T, typename a1T, typename a2T, typename a3T,
        typename a4T, typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
    TaskFnBase<fnT>* make_task(fnT fn, const a1T& a1, const a2T& a2,
            const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7, const a8T& a8, const a9T& a9, const TaskAttributes& attr) {
        typedef TaskFn<fnT, a1T, a2T, a3T, a4T, a5T, a6T, a7T, a8T, a9T> task_type;
        typename task_type::futureT r;
        task_type* task = new task_type(r, fn, a1, a2, a3, a4, a5, a6, a7, a8, a9, attr);
        return task;
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


#endif // MADNESS_WORLD_TASK_FN_H__INCLUDED
