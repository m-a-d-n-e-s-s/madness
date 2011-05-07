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


#ifndef MADNESS_WORLD_TASK_FN_BASE_H__INCLUDED
#define MADNESS_WORLD_TASK_FN_BASE_H__INCLUDED

/// \file task_fn_base.h
/// \brief Defines task functions interface.

#include <world/worlddep.h>
#include <world/worldthread.h>
#include <world/worldfut.h>
#include <world/typestuff.h>

namespace madness {

    // Forward decls
    class World;

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
            this->run(world);
        }

        /// Send this task to \c dest and submit it to the task queue there.

        /// \param dest The process where the task should be sent to run
        /// \throw MadnessException If the process is not stealable
        /// \note If this task is in the task queue, it must be deleted and set
        /// to NULL or the task will be run twice and the result future will be
        /// set twice.
        virtual void send(ProcessID) const {
            MADNESS_EXCEPTION("Task send must be implemented by the task function.", false);
        }

        World* get_world() const { return const_cast<World*>(world); }

        virtual ~TaskInterface() { if (completion) completion->notify(); }

    }; // class TaskInterface



    /// Serialization container for sending tasks to remote nodes

    /// This is for internal use only. You should not use this class directly.
    /// \tparam refT The remote reference type for task result future
    /// \tparam functionT The task function type
    template <typename refT, typename functionT>
    struct TaskHandlerInfo {
        refT ref;               ///< Remote reference for a task result future
        functionT func;         ///< A task function
        TaskAttributes attr;    ///< Task attributes
        TaskHandlerInfo(const refT& ref, functionT func, const TaskAttributes& attr)
                : ref(ref), func(func), attr(attr) {}
        TaskHandlerInfo() : ref(), func(), attr() {}

        /// Serialization of object

        /// \tparam Archive The serialization archive type
        /// \param ar The serialization archive
        template <typename Archive>
        void serialize(const Archive& ar) {
            serialize_internal<functionT>(ar);
        }

    private:

        template <typename T>
        struct is_func_ptr {
            static const bool value = (std::is_function<T>::value && std::is_pointer<T>::value)
            || std::is_member_function_pointer<T>::value;
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

    namespace detail {

        // Remove Future, const, volatile, and reference qualifiers from the type
        template <typename T>
        struct remove_fcvr{
            typedef typename remove_future< typename std::remove_cv<
                    typename std::remove_reference<T>::type >::type >::type type;
        };

    }  // namespace detail


    /// Task function base class

    /// This holds the information that is common to all tasks: the function and
    /// the result. It is also the interface passed to \c WorldTaskQueue .
    /// \tparam fnT The task function
    template <typename fnT>
    class TaskFnBase : public TaskInterface {
    public:

        typedef fnT functionT;                          ///< The task function type
        typedef typename detail::remove_fcvr<typename result_of<fnT>::type>::type resultT;
                                                        ///< The result type of the function
        typedef Future<resultT> futureT;                ///< The future result type for the task function
        typedef typename futureT::remote_refT refT;     ///< The remote reference type for the task result future
        typedef TaskHandlerInfo<refT,functionT> infoT;  ///< The task info type for serialization

    protected:
        futureT result_;            ///< The task Future result
        const functionT func_;      ///< The task function

        /// Extract task info from AM argument list.

        /// \param args The AM arguments
        static infoT get_task_info(const AmArg& args) {
            infoT info;
            args & info;
            return info;
        }

    public:
        /// Construct Task base

        /// \param res The Future result that will hold the result of the task
        /// \param func The task function to be run
        /// \param attr The task attributes
        TaskFnBase(futureT res, functionT func, const TaskAttributes& attr) :
            TaskInterface(attr), result_(res), func_(func)
        { }

        /// Construct from task info generated on a remote node
        TaskFnBase(const infoT& info) :
            TaskInterface(info.attr), result_(info.ref), func_(info.func)
        { }

        /// Construct task info object

        /// \return The task info for this task to be sent to another node
        infoT task_info() const {
            return infoT(result_.remote_ref(* TaskInterface::get_world()), func_, *this);
        }

        /// Access the result future for this task.

        /// This function is here so the task queue can take a copy of the result
        /// \return A reference to the result
        const futureT& result() const { return result_; }

        /// Virtual destructor
        virtual ~TaskFnBase() { }
    }; // class TaskFnBase

} // namespace madness


#endif // MADNESS_WORLD_TASK_FN_BASE_H__INCLUDED
