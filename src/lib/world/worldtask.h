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

  
#ifndef WORLDTASK_H
#define WORLDTASK_H

/// \file worldtask.h
/// \brief Defines TaskInterface and implements WorldTaskQueue and associated stuff.

// To do:
// a) Redo this ontop of serializable tasks which will remote much of the clutter 
//    due to multiple length argument lists.
// b) Stealing which pretty much presume a) has been done

namespace madness {

    /// Contains attributes of a task

    /// \c generator_hint : Setting this to true hints that a task
    /// will produce additional tasks and gives the task receiving
    /// preferential scheduling in order to maximize available
    /// parallelism.  The default value is false, which does not
    /// preclude a task from producing other tasks and results
    /// in normal priority for the scheduling.
    ///
    /// \c stealable : Setting this to true indicates that a task
    /// may be migrated to another process for dynamic load balancing.
    /// The default value is false.
    class TaskAttributes {
        unsigned long flags;
        static const unsigned long one = 1ul;
        static const unsigned long GENERATOR = one<<0;
        static const unsigned long STEALABLE = one<<1;
        
    public:
        TaskAttributes() : flags(0) {};
        
        TaskAttributes(bool generator_hint, bool stealable=false) 
            : flags((generator_hint?GENERATOR:0) | (stealable?STEALABLE:0))
        {};
        
        bool is_generator() const {return flags&GENERATOR;};
        
        bool is_stealable() const {return flags&STEALABLE;};
        
        void set_generator_hint (bool generator_hint) {
            if (generator_hint) flags |= GENERATOR;
            else flags &= ~GENERATOR;
        };
        
        void set_stealable (bool stealable) {
            if (stealable) flags |= STEALABLE;
            else flags &= ~STEALABLE;
        };
        
        template <typename Archive>
        void serialize(Archive& ar) {
            ar & flags;
        }
    };

    /// Returns a const reference to a task attributes with generator hint set true
    const TaskAttributes& task_attr_generator();
    
    class WorldTaskQueue;
    class TaskInterface;
    
    void task_ready_callback_function(WorldTaskQueue*, TaskInterface*);
    
    /// All tasks must be derived from this public interface
    class TaskInterface : public DependencyInterface{
        friend class WorldTaskQueue;
    private:
        TaskAttributes attr;
        class TaskReadyCallback : public CallbackInterface {
            WorldTaskQueue* taskq;
            TaskInterface* task;
        public:
            TaskReadyCallback() : taskq(0), task(0) {};
            void register_callback(WorldTaskQueue* taskq, TaskInterface* task) {
                this->taskq = taskq;
                this->task = task;
                static_cast<DependencyInterface*>(task)->register_callback(this);
            };
            void notify() {
                if (taskq) ::madness::task_ready_callback_function(taskq,task);
                taskq = 0;
            };
        } callback;   //< Embedded here for efficiency (managed by taskq)
        
    public:
        /// Create a new task with ndepend dependencies (default 0) and given attributes
        TaskInterface(int ndepend=0, const TaskAttributes& attr = TaskAttributes())
            : DependencyInterface(ndepend)
            , attr(attr)
        {};

        /// Create a new task with zero dependencies and given attributes
        TaskInterface(const TaskAttributes& attr)
            : DependencyInterface(0)
            , attr(attr)
        {};

        bool is_stealable() const {return attr.is_stealable();};

        bool is_generator() const {return attr.is_generator();};

        /// Runs the task ... derived classes MUST implement this.
        virtual void run(World& world) = 0;

        /// Virtual base class destructor so that derived class destructor is called
        virtual ~TaskInterface() {};
    };

    class WorldTaskQueue;

    template <typename functionT> class TaskFunction;
    template <typename memfunT> class TaskMemfun;

    /// Queue to manage and run tasks.
    
    /// Note that even though the queue is presently single threaded,
    /// the invocation of poll() and/or fiber scheduling can still
    /// cause multiple entries through the task layer.  If modifying
    /// this class (indeed, any World class) be careful about keeping
    /// modifications to the data structure atomic w.r.t. AM
    /// operations and the execution of multiple fibers.
    ///
    /// The pending q for tasks that need probing to have their dependencies
    /// make progress has been removed.  It can be added back.  Similarly,
    /// rather keeping track of high-priority tasks we presently have 
    /// only one priority (if it is urgent, send an active message instead).
    class WorldTaskQueue : private NO_DEFAULTS {
        friend class TaskInterface;
    private:
        //typedef std::list< TaskInterface* > listT;
        typedef DQueue< TaskInterface* > listT;
        listT ready;       //< Tasks that are ready to run

        World& world;      //< The communication context
        const int maxfiber;//< Maximum no. of active fibers
        long nregistered;  //< Counts registered tasks with simple dependencies
        const ProcessID me;
        int suspended;      //< If non-zero, task processing is suspended
        long maxinq;        //< Maximum number of enqueued tasks

    public:
        WorldTaskQueue(World& world, int maxfiber=1) 
            : world(world)
            , maxfiber(maxfiber)
            , nregistered(0)
            , me(world.mpi.rank())
            , suspended(0) 
            , maxinq(100000)
        {};

        /// Returns the number of ready tasks 
        long nready() const {
            return ready.size();
        };

        /// Returns the number of pending (not-ready) tasks
        long npending() const {
            return nregistered;
        };

        /// Returns the total number of enqueued tasks (ready + not ready)
        long size() const {
            return nready() + npending();
        };

        /// Returns the maximum number of enqueued tasks before forced task execution
        long max_enqueued() const {
            return maxinq;
        };
        
        /// Sets the maximum number of enqueued tasks before forced task execution
        void set_max_enqueued(long value) {
            maxinq = value;
        };
        
        /// Add a new local task.  The task queue will eventually delete it.

        /// The task pointer (t) is assumed to have been created with
        /// \c new and when the task is eventually run the queue
        /// will call the task's destructor using \c delete.
        ///
        /// If the task has no outstanding dependencies it is immediately
        /// put into the ready queue.
        ///
        /// If the task has outstanding dependencies then it is
        /// assumed that other activities will be calling task->dec()
        /// to decrement the dependency count.  When this count goes
        /// to zero the callback will be invoked automatically to
        /// insert the task into the ready queue.  All that will be
        /// done here is registering the callback
        inline void add(TaskInterface* t) {
            nregistered++;
            if (t->probe()) add_ready_task(t);
            else t->callback.register_callback(this,t);
            World::poll_all(true); // Adding tasks is useful work so increase polling interval
        };


        /// Invoke "resultT (*function)(void)" as a local task

        /// A future is returned to eventually hold the result of the task.
        /// Future<void> is an empty class that may be ignored.
        template <typename functionT>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))> add(functionT function, const TaskAttributes& attr=TaskAttributes()) {
            return add(me,function,attr);
        }


        /// Invoke "resultT (*function)(void)" as a task, local or remote

        /// A future is returned to eventually hold the result of the task.
        /// Future<void> is an empty class that may be ignored.
        template <typename functionT>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))> add(ProcessID where, functionT function, const TaskAttributes& attr=TaskAttributes()) {
            Future<REMFUTURE(FUNCTION_RETURNT(functionT))> result;
            if (where == me) 
                add(new TaskFunction<functionT>(result, function,attr));
            else 
                TaskFunction<functionT>::sender(world, where, result, function, attr);
            return result;
        }
        

        /// Invoke "resultT (*function)(arg1T)" as a local task
        template <typename functionT, typename arg1T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))> add(functionT function, const arg1T& arg1, const TaskAttributes& attr=TaskAttributes()) {
            return add(me, function, arg1, attr);
        }


        /// Invoke "resultT (*function)(arg1T)" as a task, local or remote

        /// A future is returned to eventually hold the result of the
        /// task.  Future<void> is an empty class that may be ignored.
        ///
        /// Arguments must be (de)serializable and must of course make
        /// sense at the remote destination.  Fundamental types,
        /// simple STL containers, and pointers to World,
        /// WorldContainer, and user-defined types derived from
        /// WorldObject<> are automatically handled.  Anything else is
        /// your problem.
        ///
        /// An argument that is a future may be used to carry
        /// dependencies for local tasks.  An unready future cannot be
        /// used as an argument for a remote tasks --- i.e., remote
        /// tasks must be ready to execute (you can work around this
        /// by making a local task to submit the remote task once
        /// everything is ready).
        template <typename functionT, typename arg1T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))> add(ProcessID where, functionT function, 
                                                const arg1T& arg1, const TaskAttributes& attr=TaskAttributes()) {
            Future<REMFUTURE(FUNCTION_RETURNT(functionT))> result;
            if (where == me) {
                add(new TaskFunction<functionT>(result, function, arg1, attr));
            }
            else {
                MADNESS_ASSERT(future_probe(arg1));  // No dependencies allowed for remote tasks
                TaskFunction<functionT>::sender(world, where, result, function, arg1, attr);
            }
            return result;
        }
        

        /// Invoke "resultT (*function)(arg1T,arg2T)" as a local task
        template <typename functionT, typename arg1T, typename arg2T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))> add(functionT function, 
                                                const arg1T& arg1, const arg2T& arg2, const TaskAttributes& attr=TaskAttributes()) {
            return add(me,function,arg1,arg2,attr);
        }


        /// Invoke "resultT (*function)(arg1T,arg2T)" as a task, local or remote

        /// A future is returned to eventually hold the result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))> add(ProcessID where, functionT function, 
                                                const arg1T& arg1, const arg2T& arg2, const TaskAttributes& attr=TaskAttributes()) {
            Future<REMFUTURE(FUNCTION_RETURNT(functionT))> result;
            if (where == me) {
                add(new TaskFunction<functionT>(result, function, arg1, arg2, attr));
            }
            else {
                MADNESS_ASSERT(future_probe(arg1));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg2));  // No dependencies allowed for remote tasks
                TaskFunction<functionT>::sender(world, where, result, function, arg1, arg2, attr);
            }
            return result;
        }
        

        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T)" as a local task
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))> add(functionT function, 
                                                const arg1T& arg1, const arg2T& arg2, 
                                                const arg3T& arg3, const TaskAttributes& attr=TaskAttributes()) {
            return add(me,function,arg1,arg2,arg3,attr);
        }


        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T)" as a task, local or remote

        /// A future is returned to eventually hold the result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))> add(ProcessID where, functionT function, 
                                                const arg1T& arg1, const arg2T& arg2,
                                                const arg3T& arg3, const TaskAttributes& attr=TaskAttributes()) {
            Future<REMFUTURE(FUNCTION_RETURNT(functionT))> result;
            if (where == me) {
                add(new TaskFunction<functionT>(result, function, arg1, arg2, arg3, attr));
            }
            else {
                MADNESS_ASSERT(future_probe(arg1));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg2));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg3));  // No dependencies allowed for remote tasks
                TaskFunction<functionT>::sender(world, where, result, function, arg1, arg2, arg3, attr);
            }
            return result;
        }

        
        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T,arg4T)" as a task, local or remote

        /// A future is returned to eventually hold the result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))> add(ProcessID where, functionT function, 
                                                const arg1T& arg1, const arg2T& arg2,
                                                const arg3T& arg3, const arg4T& arg4, const TaskAttributes& attr=TaskAttributes()) {
            Future<REMFUTURE(FUNCTION_RETURNT(functionT))> result;
            if (where == me) {
                add(new TaskFunction<functionT>(result, function, arg1, arg2, arg3, arg4, attr));
            }
            else {
                MADNESS_ASSERT(future_probe(arg1));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg2));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg3));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg4));  // No dependencies allowed for remote tasks
                TaskFunction<functionT>::sender(world, where, result, function, arg1, arg2, arg3, arg4, attr);
            }
            return result;
        }


        
        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T,arg4T,arg5T)" as a task, local or remote

        /// A future is returned to eventually hold the result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))> add(ProcessID where, functionT function, 
                                                const arg1T& arg1, const arg2T& arg2,
                                                const arg3T& arg3, const arg4T& arg4,
                                                const arg5T& arg5, const TaskAttributes& attr=TaskAttributes()) {
            Future<REMFUTURE(FUNCTION_RETURNT(functionT))> result;
            if (where == me) {
                add(new TaskFunction<functionT>(result, function, arg1, arg2, arg3, arg4, arg5, attr));
            }
            else {
                MADNESS_ASSERT(future_probe(arg1));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg2));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg3));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg4));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg5));  // No dependencies allowed for remote tasks
                TaskFunction<functionT>::sender(world, where, result, function, arg1, arg2, arg3, arg4, arg5, attr);
            }
            return result;
        }
        
        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T)" as a task, local or remote

        /// A future is returned to eventually hold the result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))> add(ProcessID where, functionT function, 
                                                const arg1T& arg1, const arg2T& arg2,
                                                const arg3T& arg3, const arg4T& arg4,
                                                const arg5T& arg5, const arg6T& arg6, const TaskAttributes& attr=TaskAttributes()) {
            Future<REMFUTURE(FUNCTION_RETURNT(functionT))> result;
            if (where == me) {
                add(new TaskFunction<functionT>(result, function, arg1, arg2, arg3, arg4, arg5, arg6, attr));
            }
            else {
                MADNESS_ASSERT(future_probe(arg1));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg2));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg3));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg4));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg5));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg6));  // No dependencies allowed for remote tasks
                TaskFunction<functionT>::sender(world, where, result, function, arg1, arg2, arg3, arg4, arg5, arg6, attr);
            }
            return result;
        }

        /// Invoke "resultT (*function)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T)" as a task, local or remote

        /// A future is returned to eventually hold the result of the
        /// task.  Future<void> is an empty class that may be ignored.
        template <typename functionT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
        Future<REMFUTURE(FUNCTION_RETURNT(functionT))> add(ProcessID where, functionT function, 
                                                const arg1T& arg1, const arg2T& arg2,
                                                const arg3T& arg3, const arg4T& arg4,
                                                const arg5T& arg5, const arg6T& arg6, const arg7T& arg7, const TaskAttributes& attr=TaskAttributes()) {
            Future<REMFUTURE(FUNCTION_RETURNT(functionT))> result;
            if (where == me) {
                add(new TaskFunction<functionT>(result, function, arg1, arg2, arg3, arg4, arg5, arg6, arg7, attr));
            }
            else {
                MADNESS_ASSERT(future_probe(arg1));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg2));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg3));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg4));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg5));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg6));  // No dependencies allowed for remote tasks
                MADNESS_ASSERT(future_probe(arg7));  // No dependencies allowed for remote tasks
                TaskFunction<functionT>::sender(world, where, result, function, arg1, arg2, arg3, arg4, arg5, arg6, arg7, attr);
            }
            return result;
        }

        
        /// Invoke "resultT (obj.*memfun)()" as a local task
        template <typename memfunT>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> add(MEMFUN_OBJT(memfunT)& obj, memfunT memfun, const TaskAttributes& attr=TaskAttributes()) 
        {
            Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> result;
            add(new TaskMemfun<memfunT>(result,obj,memfun,attr));
            return result;
        }


        /// Invoke "resultT (obj.*memfun)(arg1T)" as a local task
        template <typename memfunT, typename arg1T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> add(MEMFUN_OBJT(memfunT)& obj, 
                                            memfunT memfun,
                                            const arg1T& arg1, const TaskAttributes& attr=TaskAttributes())
        {
            Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> result;
            add(new TaskMemfun<memfunT>(result,obj,memfun,arg1,attr));
            return result;
        }


        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T)" as a local task
        template <typename memfunT, typename arg1T, typename arg2T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> add(MEMFUN_OBJT(memfunT)& obj, 
                                            memfunT memfun,
                                            const arg1T& arg1, const arg2T& arg2, const TaskAttributes& attr=TaskAttributes())
        {
            Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> result;
            add(new TaskMemfun<memfunT>(result,obj,memfun,arg1,arg2,attr));
            return result;
        }


        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3)" as a local task
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> add(MEMFUN_OBJT(memfunT)& obj, 
                                            memfunT memfun,
                                            const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const TaskAttributes& attr=TaskAttributes())
        {
            Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> result;
            add(new TaskMemfun<memfunT>(result,obj,memfun,arg1,arg2,arg3,attr));
            return result;
        }


        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3,arg4)" as a local task
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> add(MEMFUN_OBJT(memfunT)& obj, 
                                            memfunT memfun,
                                            const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const TaskAttributes& attr=TaskAttributes())
        {
            Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> result;
            add(new TaskMemfun<memfunT>(result,obj,memfun,arg1,arg2,arg3,arg4,attr));
            return result;
        }



        /// Invoke "resultT (obj.*memfun)(arg1T,arg2T,arg3,arg4,arg5)" as a local task
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> add(MEMFUN_OBJT(memfunT)& obj, 
                                            memfunT memfun,
                                            const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg4T& arg4, const arg5T& arg5, const TaskAttributes& attr=TaskAttributes())
        {
            Future<REMFUTURE(MEMFUN_RETURNT(memfunT))> result;
            add(new TaskMemfun<memfunT>(result,obj,memfun,arg1,arg2,arg3,arg4,arg5,attr));
            return result;
        }


        struct TaskQFenceProbe {
            mutable WorldTaskQueue* q;
            TaskQFenceProbe(WorldTaskQueue* q) : q(q) {};
            bool operator()() const {
                return (q->nregistered==0 && q->ready.empty());
            };
        };

        /// Returns after all local tasks have completed AND am locally fenced
        void fence() {
            // Critical here is to ensure BOTH taskq is drained and
            // the AM locally fenced when we return.  Thus, drain the
            // queue, fence the AM, and if there are any new tasks we
            // have to repeat until all is quiet.  Do NOT poll
            // once this condition has been established.
            do {
                World::await(TaskQFenceProbe(this));
                world.am.fence();
            } while (nregistered || !ready.empty());
        };

        
        /// Private:  Call back for tasks that have satisfied their dependencies.
        inline void add_ready_task(TaskInterface* t) {
            if (t->is_generator()) {
                ready.push_front(t);
                print("WorldTaskQueue: just added generator task at front of q");
            }
            else {
                ready.push_back(t);
            }
            nregistered--;
            MADNESS_ASSERT(nregistered>=0);
        };


        /// Suspend processing of tasks
        inline void suspend() {
            suspended++;
        };


        /// Resume processing of tasks
        inline void resume() {
            suspended--;
            if (suspended < 0) {
                print("TaskQueue: suspended was negative --- mismatched suspend/resume pairing");
                suspended = 0;
            };
        };


        /// Run the next ready task, returning true if one was run
        inline bool run_next_ready_task() {
            if (suspended || ready.empty()) {
                return false;
            }
            else {
                TaskInterface *p = ready.front();
                ready.pop_front();

		// NB THIS ENFORCES SEQUENTIALLY CONSISTENT EXECUTION
		// ... can be relaxed at expense of possible unbounded memory consumption
                suspend();
                p->run(world);
                resume();

                delete p;
                return true;
            }
        };
   
    };


    // Internal: Convenience for serializing
    template <typename refT, typename functionT>
    struct TaskHandlerInfo {
        refT ref;
        functionT func;
        TaskAttributes attr;
        TaskHandlerInfo(const refT& ref, functionT func, const TaskAttributes& attr) 
            : ref(ref), func(func),attr(attr) {};
        TaskHandlerInfo() {};
        template <typename Archive> 
        void serialize(const Archive& ar) {
            ar & archive::wrap_opaque(*this);
        }
    };

    // Internal: Common functionality for TaskFunction and TaskMemfun classes
    class TaskFunctionBase : public TaskInterface {
    protected:
    public:

        TaskFunctionBase(const TaskAttributes& attributes) 
            : TaskInterface(attributes) 
        {};

        // Register non-ready future as a dependency
        template <typename T>
        inline void check_dependency(Future<T>& fut) {
            if (!fut.probe()) {
                inc();
                fut.register_callback(this);
            }
        }
    };

    // Task wrapping "resultT (*function)()"
    template <typename resultT>
    struct TaskFunction<resultT (*)()> : public TaskFunctionBase {
        typedef resultT (*functionT)();
        typedef Future<REMFUTURE(resultT)> futureT;
        typedef RemoteReference< FutureImpl<REMFUTURE(resultT)> > refT;

        static void handler(World& world, ProcessID src, const AmArg& arg) {
            TaskHandlerInfo<refT,functionT> info = arg;
            world.taskq.add(new TaskFunction<functionT>(futureT(info.ref),info.func,info.attr));
        };
        
        static void sender(World& world, ProcessID dest, Future<REMFUTURE(resultT)>& result, functionT func, const TaskAttributes& attr) {
            world.am.send(dest, handler, TaskHandlerInfo<refT,functionT>(result.remote_ref(world), func, attr));
        };

        futureT result; 
        const functionT func;
        TaskFunction(const futureT& result, functionT func, const TaskAttributes& attr) 
            : TaskFunctionBase(attr)
            , result(result)
            , func(func) 
            {};
        
        void run(World& world) {
            result.set(func());
        };
    };


    // Task wrapping "resultT (*function)(arg1)"
    template <typename resultT, typename arg1_type>
    struct TaskFunction<resultT (*)(arg1_type)> : public TaskFunctionBase {
        typedef resultT (*functionT)(arg1_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef Future<REMFUTURE(resultT)> futureT;
        typedef RemoteReference< FutureImpl<REMFUTURE(resultT)> > refT;

        static void handler(World& world, ProcessID src, void* buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg*) buf;
            TaskHandlerInfo<refT,functionT> info;
            arg1T arg1;
            arg->unstuff(nbyte, info, arg1);
            world.taskq.add(new TaskFunction<functionT>(futureT(info.ref),info.func,arg1,info.attr));
        };

        template <typename a1T>
        static void sender(World& world, ProcessID dest, const futureT& result, functionT func,
                           const a1T& a1, const TaskAttributes& attr) {
            LongAmArg* arg = new_long_am_arg();
            size_t nbyte = arg->stuff(TaskHandlerInfo<refT,functionT>(result.remote_ref(world), func, attr), 
                                      static_cast<arg1T>(a1));
            world.am.send_long_managed(dest, TaskFunction<functionT>::handler, arg, nbyte);
        }

        futureT result;
        const functionT func;
        Future<arg1T> arg1;

        template <typename a1T>
        TaskFunction(const futureT& result, functionT func, const a1T& a1, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), func(func), arg1(a1) {
            check_dependency(arg1);
        }

        void run(World& world) {
            result.set(func(arg1));
        };
    };


    // Task wrapping "resultT (*function)(arg1,arg2)"
    template <typename resultT, typename arg1_type, typename arg2_type>
    struct TaskFunction<resultT (*)(arg1_type,arg2_type)> : public TaskFunctionBase {
        typedef resultT (*functionT)(arg1_type,arg2_type);

        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef Future<REMFUTURE(resultT)> futureT;
        typedef RemoteReference< FutureImpl<REMFUTURE(resultT)> > refT;

        static void handler(World& world, ProcessID src, void* buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg*) buf;
            TaskHandlerInfo<refT,functionT> info;
            arg1T arg1;
            arg2T arg2;
            arg->unstuff(nbyte, info, arg1, arg2);
            world.taskq.add(new TaskFunction<functionT>(futureT(info.ref),info.func,arg1,arg2,info.attr));
        };

        template <typename a1T, typename a2T>
        static void sender(World& world, ProcessID dest, const futureT& result, functionT func,
                           const a1T& arg1, const a2T& arg2, const TaskAttributes& attr) {
            LongAmArg* arg = new_long_am_arg();
            size_t nbyte = arg->stuff(TaskHandlerInfo<refT,functionT>(result.remote_ref(world), func, attr), 
                                      static_cast<arg1T>(arg1), static_cast<arg2T>(arg2));
            world.am.send_long_managed(dest, TaskFunction<functionT>::handler, arg, nbyte);
        }

        futureT result;
        const functionT func;
        Future<arg1T> arg1;
        Future<arg2T> arg2;

        template <typename a1T, typename a2T>
        TaskFunction(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), func(func), arg1(a1), arg2(a2) {
            check_dependency(arg1);
            check_dependency(arg2);
        }

        void run(World& world) {
            result.set(func(arg1,arg2));
        };
    };

    // Task wrapping "resultT (*function)(arg1,arg2,arg3)"
    template <typename resultT, typename arg1_type, typename arg2_type, typename arg3_type>
    struct TaskFunction<resultT (*)(arg1_type,arg2_type,arg3_type)> : public TaskFunctionBase {
        typedef resultT (*functionT)(arg1_type,arg2_type,arg3_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef Future<REMFUTURE(resultT)> futureT;
        typedef RemoteReference< FutureImpl<REMFUTURE(resultT)> > refT;

        static void handler(World& world, ProcessID src, void* buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg*) buf;
            TaskHandlerInfo<refT,functionT> info;
            arg1T arg1;
            arg2T arg2;
            arg3T arg3;
            arg->unstuff(nbyte, info, arg1, arg2, arg3);
            world.taskq.add(new TaskFunction<functionT>(futureT(info.ref),info.func,arg1,arg2,arg3,info.attr));
        };

        template <typename a1T, typename a2T, typename a3T>
        static void sender(World& world, ProcessID dest, const futureT& result, functionT func,
                           const a1T& arg1, const a2T& arg2, const a3T& arg3, const TaskAttributes& attr) {
            LongAmArg* arg = new_long_am_arg();
            size_t nbyte = arg->stuff(TaskHandlerInfo<refT,functionT>(result.remote_ref(world), func, attr), 
                                      static_cast<arg1T>(arg1), static_cast<arg2T>(arg2), static_cast<arg3T>(arg3));
            world.am.send_long_managed(dest, TaskFunction<functionT>::handler, arg, nbyte);
        }

        futureT result;
        const functionT func;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;

        template <typename a1T, typename a2T, typename a3T>
        TaskFunction(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const a3T& a3, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), func(func), arg1(a1), arg2(a2), arg3(a3) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
        }

        void run(World& world) {
            result.set(func(arg1,arg2,arg3));
        };
    };


    // Task wrapping "resultT (*function)(arg1,arg2,arg3,arg4)"
    template <typename resultT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type>
    struct TaskFunction<resultT (*)(arg1_type,arg2_type,arg3_type,arg4_type)> : public TaskFunctionBase {
        typedef resultT (*functionT)(arg1_type,arg2_type,arg3_type,arg4_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef Future<REMFUTURE(resultT)> futureT;
        typedef RemoteReference< FutureImpl<REMFUTURE(resultT)> > refT;

        static void handler(World& world, ProcessID src, void* buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg*) buf;
            TaskHandlerInfo<refT,functionT> info;
            arg1T arg1;
            arg2T arg2;
            arg3T arg3;
            arg4T arg4;
            arg->unstuff(nbyte, info, arg1, arg2, arg3, arg4);
            world.taskq.add(new TaskFunction<functionT>(futureT(info.ref),info.func,arg1,arg2,arg3,arg4,info.attr));
        };

        template <typename a1T, typename a2T, typename a3T, typename a4T>
        static void sender(World& world, ProcessID dest, const futureT& result, functionT func,
                           const a1T& arg1, const a2T& arg2, const a3T& arg3, const a4T& arg4, const TaskAttributes& attr) {
            LongAmArg* arg = new_long_am_arg();
            size_t nbyte = arg->stuff(TaskHandlerInfo<refT,functionT>(result.remote_ref(world), func, attr), 
                                      static_cast<arg1T>(arg1), static_cast<arg2T>(arg2), static_cast<arg3T>(arg3), static_cast<arg4T>(arg4));
            world.am.send_long_managed(dest, TaskFunction<functionT>::handler, arg, nbyte);
        }

        futureT result;
        const functionT func;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;

        template <typename a1T, typename a2T, typename a3T, typename a4T>
        TaskFunction(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), func(func), arg1(a1), arg2(a2), arg3(a3), arg4(a4) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
        }

        void run(World& world) {
            result.set(func(arg1,arg2,arg3,arg4));
        };
    };


    // Task wrapping "resultT (*function)(arg1,arg2,arg3,arg4,arg5)"
    template <typename resultT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type, typename arg5_type>
    struct TaskFunction<resultT (*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type)> : public TaskFunctionBase {
        typedef resultT (*functionT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T;
        typedef Future<REMFUTURE(resultT)> futureT;
        typedef RemoteReference< FutureImpl<REMFUTURE(resultT)> > refT;

        static void handler(World& world, ProcessID src, void* buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg*) buf;
            TaskHandlerInfo<refT,functionT> info;
            arg1T arg1;
            arg2T arg2;
            arg3T arg3;
            arg4T arg4;
            arg5T arg5;
            arg->unstuff(nbyte, info, arg1, arg2, arg3, arg4, arg5);
            world.taskq.add(new TaskFunction<functionT>(futureT(info.ref),info.func,arg1,arg2,arg3,arg4,arg5,info.attr));
        };

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T>
        static void sender(World& world, ProcessID dest, const futureT& result, functionT func,
                           const a1T& arg1, const a2T& arg2, const a3T& arg3, const a4T& arg4, const a5T& arg5, const TaskAttributes& attr) {
            LongAmArg* arg = new_long_am_arg();
            size_t nbyte = arg->stuff(TaskHandlerInfo<refT,functionT>(result.remote_ref(world), func, attr), 
                                      static_cast<arg1T>(arg1), static_cast<arg2T>(arg2), static_cast<arg3T>(arg3), static_cast<arg4T>(arg4), static_cast<arg5T>(arg5));
            world.am.send_long_managed(dest, TaskFunction<functionT>::handler, arg, nbyte);
        }

        futureT result;
        const functionT func;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T>
        TaskFunction(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), func(func), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
        }

        void run(World& world) {
            result.set(func(arg1,arg2,arg3,arg4,arg5));
        };
    };

    // Task wrapping "resultT (*function)(arg1,arg2,arg3,arg4,arg5,arg6)"
    template <typename resultT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type, typename arg5_type, typename arg6_type>
    struct TaskFunction<resultT (*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type)> : public TaskFunctionBase {
        typedef resultT (*functionT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T;
        typedef REMFUTURE(REMCONST(REMREF(arg6_type))) arg6T;
        typedef Future<REMFUTURE(resultT)> futureT;
        typedef RemoteReference< FutureImpl<REMFUTURE(resultT)> > refT;

        static void handler(World& world, ProcessID src, void* buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg*) buf;
            TaskHandlerInfo<refT,functionT> info;
            arg1T arg1;
            arg2T arg2;
            arg3T arg3;
            arg4T arg4;
            arg5T arg5;
            arg6T arg6;
            arg->unstuff(nbyte, info, arg1, arg2, arg3, arg4, arg5, arg6);
            world.taskq.add(new TaskFunction<functionT>(futureT(info.ref),info.func,arg1,arg2,arg3,arg4,arg5,arg6,info.attr));
        };

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T>
        static void sender(World& world, ProcessID dest, const futureT& result, functionT func,
                           const a1T& arg1, const a2T& arg2, const a3T& arg3, const a4T& arg4, const a5T& arg5, const a6T& arg6, const TaskAttributes& attr) {
            LongAmArg* arg = new_long_am_arg();
            size_t nbyte = arg->stuff(TaskHandlerInfo<refT,functionT>(result.remote_ref(world), func, attr), 
                                      static_cast<arg1T>(arg1), static_cast<arg2T>(arg2), static_cast<arg3T>(arg3), static_cast<arg4T>(arg4), static_cast<arg5T>(arg5), static_cast<arg6T>(arg6));
            world.am.send_long_managed(dest, TaskFunction<functionT>::handler, arg, nbyte);
        }

        futureT result;
        const functionT func;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;
        Future<arg6T> arg6;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T>
        TaskFunction(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), func(func), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5), arg6(a6) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
            check_dependency(arg6);
        }

        void run(World& world) {
            result.set(func(arg1,arg2,arg3,arg4,arg5,arg6));
        };
    };


    // Task wrapping "resultT (*function)(arg1,arg2,arg3,arg4,arg5,arg6,arg7)"
    template <typename resultT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type, typename arg5_type, typename arg6_type, typename arg7_type>
    struct TaskFunction<resultT (*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type,arg7_type)> : public TaskFunctionBase {
        typedef resultT (*functionT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type,arg7_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T;
        typedef REMFUTURE(REMCONST(REMREF(arg6_type))) arg6T;
        typedef REMFUTURE(REMCONST(REMREF(arg7_type))) arg7T;
        typedef Future<REMFUTURE(resultT)> futureT;
        typedef RemoteReference< FutureImpl<REMFUTURE(resultT)> > refT;

        static void handler(World& world, ProcessID src, void* buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg*) buf;
            TaskHandlerInfo<refT,functionT> info;
            arg1T arg1;
            arg2T arg2;
            arg3T arg3;
            arg4T arg4;
            arg5T arg5;
            arg6T arg6;
            arg7T arg7;
            arg->unstuff(nbyte, info, arg1, arg2, arg3, arg4, arg5, arg6, arg7);
            world.taskq.add(new TaskFunction<functionT>(futureT(info.ref),info.func,arg1,arg2,arg3,arg4,arg5,arg6,arg7,info.attr));
        };

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T, typename a7T>
        static void sender(World& world, ProcessID dest, const futureT& result, functionT func,
                           const a1T& arg1, const a2T& arg2, const a3T& arg3, const a4T& arg4, const a5T& arg5, const a6T& arg6, const a7T& arg7, const TaskAttributes& attr) {
            LongAmArg* arg = new_long_am_arg();
            size_t nbyte = arg->stuff(TaskHandlerInfo<refT,functionT>(result.remote_ref(world), func, attr), 
                                      static_cast<arg1T>(arg1), static_cast<arg2T>(arg2), static_cast<arg3T>(arg3), static_cast<arg4T>(arg4), static_cast<arg5T>(arg5), static_cast<arg6T>(arg6), static_cast<arg7T>(arg7));
            world.am.send_long_managed(dest, TaskFunction<functionT>::handler, arg, nbyte);
        }

        futureT result;
        const functionT func;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;
        Future<arg6T> arg6;
        Future<arg7T> arg7;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T, typename a7T>
        TaskFunction(const futureT& result, functionT func, const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), func(func), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5), arg6(a6), arg7(a7) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
            check_dependency(arg6);
            check_dependency(arg7);
        }

        void run(World& world) {
            result.set(func(arg1,arg2,arg3,arg4,arg5,arg6,arg7));
        };
    };


    // Task wrapping "resultT (obj.*function)()"
    template <typename resultT, typename objT>
    struct TaskMemfun<resultT (objT::*)()> : public TaskFunctionBase {
        typedef resultT (objT::*memfunT)();
        typedef Future<REMFUTURE(resultT)> futureT;

        futureT result;
        objT& obj;
        const memfunT memfun;

        TaskMemfun(const futureT& result, objT& obj, memfunT memfun, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun) {}

        void run(World& world) {
            result.set((obj.*memfun)());
        };
    };

    // Task wrapping "resultT (obj.*function)(arg1)"
    template <typename resultT, typename objT, typename arg1_type>
    struct TaskMemfun<resultT (objT::*)(arg1_type)> : public TaskFunctionBase {
        typedef resultT (objT::*memfunT)(arg1_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef Future<REMFUTURE(resultT)> futureT;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;

        template <typename a1T>
            TaskMemfun(const futureT& result, objT& obj, memfunT memfun, const a1T& a1, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1) {
            check_dependency(arg1);
        }

        void run(World& world) {
            result.set((obj.*memfun)(arg1));
        };
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2)"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type>
    struct TaskMemfun<resultT (objT::*)(arg1_type,arg2_type)> : public TaskFunctionBase {
        typedef resultT (objT::*memfunT)(arg1_type,arg2_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef Future<REMFUTURE(resultT)> futureT;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;

        template <typename a1T, typename a2T>
            TaskMemfun(const futureT& result, objT& obj, memfunT memfun, const a1T& a1, const a2T& a2, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2) {
            check_dependency(arg1);
            check_dependency(arg2);
        }

        void run(World& world) {
            result.set((obj.*memfun)(arg1,arg2));
        };
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3)"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type>
    struct TaskMemfun<resultT (objT::*)(arg1_type,arg2_type,arg3_type)> : public TaskFunctionBase {
        typedef resultT (objT::*memfunT)(arg1_type,arg2_type,arg3_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef Future<REMFUTURE(resultT)> futureT;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;

        template <typename a1T, typename a2T, typename a3T>
            TaskMemfun(const futureT& result, objT& obj, memfunT memfun, 
                       const a1T& a1, const a2T& a2, const a3T& a3, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
        }

        void run(World& world) {
            result.set((obj.*memfun)(arg1,arg2,arg3));
        };
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3,arg4)"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type>
    struct TaskMemfun<resultT (objT::*)(arg1_type,arg2_type,arg3_type,arg4_type)> : public TaskFunctionBase {
        typedef resultT (objT::*memfunT)(arg1_type,arg2_type,arg3_type,arg4_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef Future<REMFUTURE(resultT)> futureT;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;

        template <typename a1T, typename a2T, typename a3T, typename a4T>
            TaskMemfun(const futureT& result, objT& obj, memfunT memfun, 
                       const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3), arg4(a4) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
        }

        void run(World& world) {
            result.set((obj.*memfun)(arg1,arg2,arg3,arg4));
        };
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3,arg4,arg5)"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type, typename arg5_type>
    struct TaskMemfun<resultT (objT::*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type)> : public TaskFunctionBase {
        typedef resultT (objT::*memfunT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T;
        typedef Future<REMFUTURE(resultT)> futureT;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T>
            TaskMemfun(const futureT& result, objT& obj, memfunT memfun, 
                       const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, a5T& a5, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
        }

        void run(World& world) {
            result.set((obj.*memfun)(arg1,arg2,arg3,arg4,arg5));
        };
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3,arg4,arg5,arg6)"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type, typename arg5_type, typename arg6_type>
    struct TaskMemfun<resultT (objT::*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type)> : public TaskFunctionBase {
        typedef resultT (objT::*memfunT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type);
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg6_type))) arg6T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T;
        typedef Future<REMFUTURE(resultT)> futureT;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;
        Future<arg6T> arg6;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T>
            TaskMemfun(const futureT& result, objT& obj, memfunT memfun, 
                       const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, a5T& a5, a6T& a6, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5), arg6(a6) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
            check_dependency(arg6);
        }

        void run(World& world) {
            result.set((obj.*memfun)(arg1,arg2,arg3,arg4,arg5,arg6));
        };
    };

    //
    // Same as above but now const
    //

    // Task wrapping "resultT (obj.*function)() const"
    template <typename resultT, typename objT>
    struct TaskMemfun<resultT (objT::*)()const> : public TaskFunctionBase {
        typedef resultT (objT::*memfunT)()const;
        typedef Future<REMFUTURE(resultT)> futureT;

        futureT result;
        objT& obj;
        const memfunT memfun;

        TaskMemfun(const futureT& result, objT& obj, memfunT memfun, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun) {};

        void run(World& world) {
            result.set((obj.*memfun)());
        };
    };

    // Task wrapping "resultT (obj.*function)(arg1)const"
    template <typename resultT, typename objT, typename arg1_type>
    struct TaskMemfun<resultT (objT::*)(arg1_type)const> : public TaskFunctionBase {
        typedef resultT (objT::*memfunT)(arg1_type)const;
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef Future<REMFUTURE(resultT)> futureT;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;

        template <typename a1T>
            TaskMemfun(const futureT& result, objT& obj, memfunT memfun, const a1T& a1, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1) {
            check_dependency(arg1);
        }

        void run(World& world) {
            result.set((obj.*memfun)(arg1));
        };
    };


    // Task wrapping "resultT (obj.*function)(arg1,arg2)const"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type>
    struct TaskMemfun<resultT (objT::*)(arg1_type,arg2_type)const> : public TaskFunctionBase {
        typedef resultT (objT::*memfunT)(arg1_type,arg2_type)const;
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef Future<REMFUTURE(resultT)> futureT;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;

        template <typename a1T, typename a2T>
            TaskMemfun(const futureT& result, objT& obj, memfunT memfun, 
                       const a1T& a1, const a2T& a2, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2) {
            check_dependency(arg1);
            check_dependency(arg2);
        }

        void run(World& world) {
            result.set((obj.*memfun)(arg1,arg2));
        };
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3)const"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type>
    struct TaskMemfun<resultT (objT::*)(arg1_type,arg2_type,arg3_type)const> : public TaskFunctionBase {
        typedef resultT (objT::*memfunT)(arg1_type,arg2_type,arg3_type)const;
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef Future<REMFUTURE(resultT)> futureT;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;

        template <typename a1T, typename a2T, typename a3T>
            TaskMemfun(const futureT& result, objT& obj, memfunT memfun, 
                       const a1T& a1, const a2T& a2, const a3T& a3, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
        }

        void run(World& world) {
            result.set((obj.*memfun)(arg1,arg2,arg3));
        };
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3,arg4)const"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type>
    struct TaskMemfun<resultT (objT::*)(arg1_type,arg2_type,arg3_type,arg4_type)const> : public TaskFunctionBase {
        typedef resultT (objT::*memfunT)(arg1_type,arg2_type,arg3_type,arg4_type)const;
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef Future<REMFUTURE(resultT)> futureT;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;

        template <typename a1T, typename a2T, typename a3T, typename a4T>
            TaskMemfun(const futureT& result, objT& obj, memfunT memfun, 
                       const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3), arg4(a4) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
        }

        void run(World& world) {
            result.set((obj.*memfun)(arg1,arg2,arg3,arg4));
        };
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3,arg4,arg5)const"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type, typename arg5_type>
    struct TaskMemfun<resultT (objT::*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type)const> : public TaskFunctionBase {
        typedef resultT (objT::*memfunT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type)const;
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T;
        typedef Future<REMFUTURE(resultT)> futureT;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T>
            TaskMemfun(const futureT& result, objT& obj, memfunT memfun, 
                       const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, a5T& a5, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
        }

        void run(World& world) {
            result.set((obj.*memfun)(arg1,arg2,arg3,arg4,arg5));
        };
    };

    // Task wrapping "resultT (obj.*function)(arg1,arg2,arg3,arg4,arg5,arg6)const"
    template <typename resultT, typename objT, typename arg1_type, typename arg2_type, typename arg3_type, typename arg4_type, typename arg5_type, typename arg6_type>
    struct TaskMemfun<resultT (objT::*)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type)const> : public TaskFunctionBase {
        typedef resultT (objT::*memfunT)(arg1_type,arg2_type,arg3_type,arg4_type,arg5_type,arg6_type)const;
        typedef REMFUTURE(REMCONST(REMREF(arg1_type))) arg1T;
        typedef REMFUTURE(REMCONST(REMREF(arg2_type))) arg2T;
        typedef REMFUTURE(REMCONST(REMREF(arg3_type))) arg3T;
        typedef REMFUTURE(REMCONST(REMREF(arg4_type))) arg4T;
        typedef REMFUTURE(REMCONST(REMREF(arg5_type))) arg5T; 
        typedef REMFUTURE(REMCONST(REMREF(arg6_type))) arg6T;
        typedef Future<REMFUTURE(resultT)> futureT;

        futureT result;
        objT& obj;
        const memfunT memfun;
        Future<arg1T> arg1;
        Future<arg2T> arg2;
        Future<arg3T> arg3;
        Future<arg4T> arg4;
        Future<arg5T> arg5;
        Future<arg5T> arg6;

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T, typename a6T>
            TaskMemfun(const futureT& result, objT& obj, memfunT memfun, 
                       const a1T& a1, const a2T& a2, const a3T& a3, const a4T& a4, a5T& a5, a6T& a6, const TaskAttributes& attr) 
            : TaskFunctionBase(attr), result(result), obj(obj), memfun(memfun), arg1(a1), arg2(a2), arg3(a3), arg4(a4), arg5(a5), arg6(a6) {
            check_dependency(arg1);
            check_dependency(arg2);
            check_dependency(arg3);
            check_dependency(arg4);
            check_dependency(arg5);
            check_dependency(arg6);
        }

        void run(World& world) {
            result.set((obj.*memfun)(arg1,arg2,arg3,arg4,arg5,arg6));
        };
    };

}


#endif
