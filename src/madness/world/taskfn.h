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

#ifndef MADNESS_WORLD_TASKFN_H__INCLUDED
#define MADNESS_WORLD_TASKFN_H__INCLUDED

#include <type_traits>
#include <madness/world/dependency_interface.h>
#include <madness/world/thread.h>
#include <madness/world/future.h>

namespace madness {


    // Forward decls
    class World;
    class WorldTaskQueue;

    namespace detail {
        template <typename, typename, typename> class MemFuncWrapper;
        template <typename ptrT, typename memfnT, typename resT>
        memfnT get_mem_func_ptr(const MemFuncWrapper<ptrT, memfnT, resT>&);
    }

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
    class TaskInterface : public PoolTaskInterface, public DependencyInterface {
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
                : PoolTaskInterface(attr)
                , DependencyInterface(ndepend)
                , world(0)
                , completion(0)
                , submit(this)
        {}

        /// Create a new task with zero dependencies and given attributes
        explicit TaskInterface(const TaskAttributes& attr)
                : PoolTaskInterface(attr)
                , DependencyInterface(0)
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
        virtual void run(World&) {
            //print("in virtual run(world) method");
            MADNESS_EXCEPTION("World TaskInterface: user did not implement one of run(world) or run(world, taskthreadenv)", 0);
        }

        /// Runs a multi-threaded task
        virtual void run(World& world, const TaskThreadEnv& env) {
            //print("in virtual run(world,env) method", env.nthread(), env.id());
            if (env.nthread() != 1)
                MADNESS_EXCEPTION("World TaskInterface: user did not implement run(world, taskthreadenv) for multithreaded task", 0);
            run(world);
        }

        World* get_world() const { return const_cast<World*>(world); }

        virtual ~TaskInterface() { if (completion) completion->notify(); }

    }; // class TaskInterface

    namespace detail {

        template <typename T>
        struct ArgCountHelper {
            static const unsigned int value = 1u;
        };

        template <>
        struct ArgCountHelper<void> {
            static const unsigned int value = 0u;
        };

        // Counts the number of arguments that will be given to a task function
        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T,
                typename a6T, typename a7T, typename a8T, typename a9T>
        struct ArgCount : public std::integral_constant<unsigned int, ArgCountHelper<a1T>::value
                + ArgCountHelper<a2T>::value + ArgCountHelper<a3T>::value
                + ArgCountHelper<a4T>::value + ArgCountHelper<a5T>::value
                + ArgCountHelper<a6T>::value + ArgCountHelper<a7T>::value
                + ArgCountHelper<a8T>::value + ArgCountHelper<a9T>::value> {
        };


        /// A wrapper object for holding task function objects
        template <typename Arg>
        class ArgHolder : private NO_DEFAULTS{
        private:
            Arg arg_;
        public:

            ArgHolder(const Arg& arg) : arg_(arg) { }

            ArgHolder(const archive::BufferInputArchive& input_arch) :
                arg_()
            {
                input_arch & arg_;
            }

            operator Arg&() { return arg_; }
        }; // class ArgHolder


        template <typename T>
        struct task_arg {
            typedef T type;
            typedef ArgHolder<T> holderT;
        };

        template <typename T>
        struct task_arg<Future<T> > {
            typedef T type;
            typedef Future<T> holderT;
        };

        template <>
        struct task_arg<Future<void> > {
            typedef const Future<void> type;
            typedef const Future<void> holderT;
        };

        template <>
        struct task_arg<void> {
            typedef const Future<void> type;
            typedef const Future<void> holderT;
        };

        template <typename fnT>
        struct task_result_type {
            typedef typename madness::remove_fcvr<typename madness::detail::result_of<fnT>::type>::type resultT;
            typedef Future<resultT> futureT;
            typedef futureT type;
        };

        /// Future<void> is an empty object which is used as a placeholder for
        /// unused arguments.
        typedef Future<void> voidT;


        // These functions are used to differentiate the task function calls
        // based on the number of arguments and return type.

        // Note: voidT arguments must be const or the wrong function will be
        // selected.

        template <typename fnT>
        inline typename std::enable_if<std::is_void<typename detail::result_of<fnT>::type>::value >::type
        run_function(Future<void>& result, fnT fn, const voidT&,
                const voidT&, const voidT&, const voidT&, const voidT&,
                const voidT&, const voidT&, const voidT&, const voidT&)
        { fn(); result.set(); }

        template <typename fnT, typename a1T>
        inline typename std::enable_if<std::is_void<typename detail::result_of<fnT>::type>::value >::type
        run_function(Future<void>& result, fnT fn, a1T& a1,
                const voidT&, const voidT&, const voidT&, const voidT&,
                const voidT&, const voidT&, const voidT&, const voidT&)
        { fn(a1); result.set(); }

        template <typename fnT, typename a1T, typename a2T>
        inline typename std::enable_if<std::is_void<typename detail::result_of<fnT>::type>::value >::type
        run_function(Future<void>& result, fnT fn, a1T& a1, a2T& a2,
                const voidT&, const voidT&, const voidT&, const voidT&,
                const voidT&, const voidT&, const voidT&)
        { fn(a1, a2); result.set(); }

        template <typename fnT, typename a1T, typename a2T, typename a3T>
        inline typename std::enable_if<std::is_void<typename detail::result_of<fnT>::type>::value >::type
        run_function(Future<void>& result, fnT fn, a1T& a1, a2T& a2,
                a3T& a3, const voidT&, const voidT&, const voidT&, const voidT&,
                const voidT&, const voidT&)
        { fn(a1, a2, a3); result.set(); }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T>
        inline typename std::enable_if<std::is_void<typename detail::result_of<fnT>::type>::value >::type
        run_function(Future<void>& result, fnT fn, a1T& a1, a2T& a2,
                a3T& a3, a4T& a4, const voidT&, const voidT&, const voidT&,
                const voidT&, const voidT&)
        { fn(a1, a2, a3, a4); result.set(); }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T>
        inline typename std::enable_if<std::is_void<typename detail::result_of<fnT>::type>::value >::type
        run_function(Future<void>& result, fnT fn, a1T& a1, a2T& a2,
                a3T& a3, a4T& a4, a5T& a5, const voidT&, const voidT&,
                const voidT&, const voidT&)
        { fn(a1, a2, a3, a4, a5); result.set(); }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T>
        inline typename std::enable_if<std::is_void<typename detail::result_of<fnT>::type>::value >::type
        run_function(Future<void>& result, fnT fn, a1T& a1, a2T& a2,
                a3T& a3, a4T& a4, a5T& a5, a6T& a6, const voidT&, const voidT&,
                const voidT&)
        { fn(a1, a2, a3, a4, a5, a6); result.set(); }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T>
        inline typename std::enable_if<std::is_void<typename detail::result_of<fnT>::type>::value >::type
        run_function(Future<void>& result, fnT fn, a1T& a1, a2T& a2,
                a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7, const voidT&,
                const voidT&)
        { fn(a1, a2, a3, a4, a5, a6, a7); result.set(); }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T, typename a8T>
        inline typename std::enable_if<std::is_void<typename detail::result_of<fnT>::type>::value >::type
        run_function(Future<void>& result, fnT fn, a1T& a1, a2T& a2,
                a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7, a8T& a8, const voidT&)
        { fn(a1, a2, a3, a4, a5, a6, a7, a8); result.set(); }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
            typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
        inline typename std::enable_if<std::is_void<typename detail::result_of<fnT>::type>::value >::type
        run_function(Future<void>& result, fnT fn, a1T& a1, a2T& a2,
                a3T& a3, a4T& a4, a5T& a5, a6T& a6, a7T& a7, a8T& a8, a9T& a9)
        { fn(a1, a2, a3, a4, a5, a6, a7, a8, a9); result.set(); }

        template <typename fnT>
        inline typename std::enable_if<!std::is_void<typename detail::result_of<fnT>::type>::value >::type
        run_function(typename task_result_type<fnT>::futureT& result,
                fnT fn, const voidT&, const voidT&, const voidT&, const voidT&,
                const voidT&, const voidT&, const voidT&, const voidT&, const voidT&)
        { result.set(fn()); }

        template <typename fnT, typename a1T>
        inline typename std::enable_if<!std::is_void<typename detail::result_of<fnT>::type>::value >::type
        run_function(typename task_result_type<fnT>::futureT& result,
                fnT fn, a1T& a1, const voidT&, const voidT&, const voidT&,
                const voidT&, const voidT&, const voidT&, const voidT&, const voidT&)
        { result.set(fn(a1)); }

        template <typename fnT, typename a1T, typename a2T>
        inline typename std::enable_if<!std::is_void<typename detail::result_of<fnT>::type>::value >::type
        run_function(typename task_result_type<fnT>::futureT& result,
                fnT fn, a1T& a1, a2T& a2, const voidT&, const voidT&, const voidT&,
                const voidT&, const voidT&, const voidT&, const voidT&)
        { result.set(fn(a1, a2)); }

        template <typename fnT, typename a1T, typename a2T, typename a3T>
        inline typename std::enable_if<!std::is_void<typename detail::result_of<fnT>::type>::value >::type
        run_function(typename task_result_type<fnT>::futureT& result,
                fnT fn, a1T& a1, a2T& a2, a3T& a3, const voidT&, const voidT&,
                const voidT&, const voidT&, const voidT&, const voidT&)
        { result.set(fn(a1, a2, a3)); }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T>
        inline typename std::enable_if<!std::is_void<typename detail::result_of<fnT>::type>::value >::type
        run_function(typename task_result_type<fnT>::futureT& result,
                fnT fn, a1T& a1, a2T& a2, a3T& a3, a4T& a4, const voidT&,
                const voidT&, const voidT&, const voidT&, const voidT&)
        { result.set(fn(a1, a2, a3, a4)); }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T>
        inline typename std::enable_if<!std::is_void<typename detail::result_of<fnT>::type>::value >::type
        run_function(typename task_result_type<fnT>::futureT& result,
                fnT fn, a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, const voidT&,
                const voidT&, const voidT&, const voidT&)
        { result.set(fn(a1, a2, a3, a4, a5)); }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T>
        inline typename std::enable_if<!std::is_void<typename detail::result_of<fnT>::type>::value >::type
        run_function(typename task_result_type<fnT>::futureT& result,
                fnT fn, a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6,
                const voidT&, const voidT&, const voidT&)
        { result.set(fn(a1, a2, a3, a4, a5, a6)); }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T>
        inline typename std::enable_if<!std::is_void<typename detail::result_of<fnT>::type>::value >::type
        run_function(typename task_result_type<fnT>::futureT& result,
                fnT fn, a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6,
                a7T& a7, const voidT&, const voidT&)
        { result.set(fn(a1, a2, a3, a4, a5, a6, a7)); }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T, typename a8T>
        inline typename std::enable_if<!std::is_void<typename detail::result_of<fnT>::type>::value >::type
        run_function(typename task_result_type<fnT>::futureT& result,
                fnT fn, a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6,
                a7T& a7, a8T& a8, const voidT&)
        { result.set(fn(a1, a2, a3, a4, a5, a6, a7, a8)); }

        template <typename fnT, typename a1T, typename a2T, typename a3T, typename a4T,
                typename a5T, typename a6T, typename a7T, typename a8T, typename a9T>
        inline typename std::enable_if<!std::is_void<typename detail::result_of<fnT>::type>::value >::type
        run_function(typename task_result_type<fnT>::futureT& result,
                fnT fn, a1T& a1, a2T& a2, a3T& a3, a4T& a4, a5T& a5, a6T& a6,
                a7T& a7, a8T& a8, a9T& a9)
        { result.set(fn(a1, a2, a3, a4, a5, a6, a7, a8, a9)); }

    } // namespace detail

    /// Wrap a callable object and its arguments into a task function

    /// The callable object may have up to 10 arguments
    template <typename fnT, typename arg1T = void, typename arg2T = void,
            typename arg3T = void, typename arg4T = void, typename arg5T = void,
            typename arg6T = void, typename arg7T = void, typename arg8T = void,
            typename arg9T = void>
    struct TaskFn : public TaskInterface {
    private:
        /// This class type
        typedef TaskFn<fnT, arg1T, arg2T, arg3T, arg4T, arg5T, arg6T, arg7T,
                arg8T, arg9T> TaskFn_;

    public:

        typedef fnT functionT; ///< The task function type
        typedef typename detail::task_result_type<fnT>::resultT resultT;
        ///< The result type of the function
        typedef typename detail::task_result_type<fnT>::futureT futureT;

        // argument value typedefs

        static const unsigned int arity = detail::ArgCount<arg1T, arg2T, arg3T,
                arg4T, arg5T, arg6T, arg7T, arg8T, arg9T>::value;
        ///< The number of arguments given for the function
        ///< \note This may not match the arity of the function
        ///< if it has default parameter values

    private:
        futureT result_; ///< The task Future result
        const functionT func_; ///< The task function

        // If the value of the argument is known at the time the
        // Note: The type argNT for argN, where N  is > arity should be void

        typename detail::task_arg<arg1T>::holderT arg1_;///< Argument 1 that will be given to the function
        typename detail::task_arg<arg2T>::holderT arg2_;///< Argument 2 that will be given to the function
        typename detail::task_arg<arg3T>::holderT arg3_;///< Argument 3 that will be given to the function
        typename detail::task_arg<arg4T>::holderT arg4_;///< Argument 4 that will be given to the function
        typename detail::task_arg<arg5T>::holderT arg5_;///< Argument 5 that will be given to the function
        typename detail::task_arg<arg6T>::holderT arg6_;///< Argument 6 that will be given to the function
        typename detail::task_arg<arg7T>::holderT arg7_;///< Argument 7 that will be given to the function
        typename detail::task_arg<arg8T>::holderT arg8_;///< Argument 8 that will be given to the function
        typename detail::task_arg<arg9T>::holderT arg9_;///< Argument 9 that will be given to the function

        template <typename fT>
        static fT& get_func(fT& f) { return f; }

        template <typename ptrT, typename memfnT, typename resT>
        static memfnT get_func(const detail::MemFuncWrapper<ptrT, memfnT, resT>& wrapper) {
            return detail::get_mem_func_ptr(wrapper);
        }

        virtual void get_id(std::pair<void*,unsigned short>& id) const {
            return make_id(id, get_func(func_));
        }


        /// Register non-ready future as a dependency

        /// \tparam T The type of the future to check
        /// \param fut The future to check
        template <typename T>
        inline void check_dependency(Future<T>& fut) {
            if(!fut.probe()) {
                DependencyInterface::inc();
                fut.register_callback(this);
            }
        }


        /// None future arguments are always ready => no op
        template <typename T>
        inline void check_dependency(detail::ArgHolder<std::vector<Future<T> > >& arg) {
            check_dependency(static_cast<std::vector<Future<T> >&>(arg));
        }

        /// None future arguments are always ready => no op
        template <typename T>
        inline void check_dependency(std::vector<Future<T> >& vec) {
            for(typename std::vector<Future<T> >::iterator it = vec.begin(); it != vec.end(); ++it)
                check_dependency(*it);
        }

        /// Future<void> is always ready => no op
        inline void check_dependency(const std::vector<Future<void> >&) { }

        /// None future arguments are always ready => no op
        template <typename T>
        inline void check_dependency(const detail::ArgHolder<T>&) { }

        /// Future<void> is always ready => no op
        inline void check_dependency(const Future<void>&) { }

        /// Check dependencies and register callbacks where necessary
        void check_dependencies() {
            check_dependency(arg1_);
            check_dependency(arg2_);
            check_dependency(arg3_);
            check_dependency(arg4_);
            check_dependency(arg5_);
            check_dependency(arg6_);
            check_dependency(arg7_);
            check_dependency(arg8_);
            check_dependency(arg9_);
        }

        // Copies are not allowed.
        TaskFn(const TaskFn_&);
        TaskFn_ operator=(TaskFn_&);

    public:

        TaskFn(const futureT& result, functionT func, const TaskAttributes& attr) :
            TaskInterface(attr), result_(result), func_(func), arg1_(), arg2_(),
            arg3_(), arg4_(), arg5_(), arg6_(), arg7_(), arg8_(), arg9_()
        {
            MADNESS_ASSERT(arity == 0u);
            check_dependencies();
        }

        template <typename a1T>
        TaskFn(const futureT& result, functionT func, const a1T& a1,
                const TaskAttributes& attr) :
            TaskInterface(attr), result_(result), func_(func), arg1_(a1), arg2_(),
            arg3_(), arg4_(), arg5_(), arg6_(), arg7_(), arg8_(), arg9_()
        {
            MADNESS_ASSERT(arity == 1u);
            check_dependencies();
        }

        template <typename a1T, typename a2T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const TaskAttributes& attr = TaskAttributes()) :
            TaskInterface(attr), result_(result), func_(func), arg1_(a1), arg2_(a2),
            arg3_(), arg4_(), arg5_(), arg6_(), arg7_(), arg8_(), arg9_()
        {
            MADNESS_ASSERT(arity == 2u);
            check_dependencies();
        }

        template <typename a1T, typename a2T, typename a3T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const TaskAttributes& attr) :
            TaskInterface(attr), result_(result), func_(func), arg1_(a1), arg2_(a2),
            arg3_(a3), arg4_(), arg5_(), arg6_(), arg7_(), arg8_(), arg9_()
        {
            MADNESS_ASSERT(arity == 3u);
            check_dependencies();
        }

        template <typename a1T, typename a2T, typename a3T, typename a4T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const TaskAttributes& attr) :
            TaskInterface(attr), result_(result), func_(func), arg1_(a1), arg2_(a2),
            arg3_(a3), arg4_(a4), arg5_(), arg6_(), arg7_(), arg8_(), arg9_()
        {
            MADNESS_ASSERT(arity == 4u);
            check_dependencies();
        }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const TaskAttributes& attr) :
            TaskInterface(attr), result_(result), func_(func), arg1_(a1), arg2_(a2),
            arg3_(a3), arg4_(a4), arg5_(a5), arg6_(), arg7_(), arg8_(), arg9_()
        {
            MADNESS_ASSERT(arity == 5u);
            check_dependencies();
        }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T,
                typename a6T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6,
                const TaskAttributes& attr) :
            TaskInterface(attr), result_(result), func_(func), arg1_(a1), arg2_(a2),
            arg3_(a3), arg4_(a4), arg5_(a5), arg6_(a6), arg7_(), arg8_(), arg9_()
        {
            MADNESS_ASSERT(arity == 6u);
            check_dependencies();
        }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T,
                typename a6T, typename a7T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6,
                const a7T& a7, const TaskAttributes& attr) :
            TaskInterface(attr), result_(result), func_(func), arg1_(a1), arg2_(a2),
            arg3_(a3), arg4_(a4), arg5_(a5), arg6_(a6), arg7_(a7), arg8_(), arg9_()
        {
            MADNESS_ASSERT(arity == 7u);
            check_dependencies();
        }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T,
                typename a6T, typename a7T, typename a8T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7,
                const a8T& a8, const TaskAttributes& attr) :
            TaskInterface(attr), result_(result), func_(func), arg1_(a1), arg2_(a2),
            arg3_(a3), arg4_(a4), arg5_(a5), arg6_(a6), arg7_(a7), arg8_(a8), arg9_()
        {
            MADNESS_ASSERT(arity == 8u);
            check_dependencies();
        }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T,
                typename a6T, typename a7T, typename a8T, typename a9T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6,
                const a7T& a7, const a8T& a8, const a9T& a9, const TaskAttributes& attr) :
            TaskInterface(attr), result_(result), func_(func), arg1_(a1), arg2_(a2),
            arg3_(a3), arg4_(a4), arg5_(a5), arg6_(a6), arg7_(a7), arg8_(a8), arg9_(a9)
        {
            MADNESS_ASSERT(arity == 9u);
            check_dependencies();
        }

        TaskFn(const futureT& result, functionT func, const TaskAttributes& attr,
                archive::BufferInputArchive& input_arch) :
            TaskInterface(attr), result_(result), func_(func), arg1_(input_arch),
            arg2_(input_arch), arg3_(input_arch), arg4_(input_arch), arg5_(input_arch),
            arg6_(input_arch), arg7_(input_arch), arg8_(input_arch), arg9_(input_arch)
        {
            // No need to check dependencies since the arguments are from an archive
        }

        virtual ~TaskFn() { }

        const futureT& result() const { return result_; }


#ifdef HAVE_INTEL_TBB
        virtual tbb::task* execute() {
            detail::run_function(result_, func_, arg1_, arg2_, arg3_, arg4_,
                    arg5_, arg6_, arg7_, arg8_, arg9_);
            return nullptr;
        }
#else
      protected:
        virtual void run(const TaskThreadEnv& env) {
            detail::run_function(result_, func_, arg1_, arg2_, arg3_, arg4_,
                    arg5_, arg6_, arg7_, arg8_, arg9_);
        }
#endif // HAVE_INTEL_TBB

    }; // class TaskFn

} // namespace madness


#endif // MADNESS_WORLD_TASKFN_H__INCLUDED
