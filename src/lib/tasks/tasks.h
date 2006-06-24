#ifndef TASKS_H_
#define TASKS_H_

/// \file tasks.h
/// \brief Implements TaskQueue, TaskInterface and other Task classes
#include <list>
#include <unistd.h>


namespace madness {
    // more efficient than the current polling loop(s)
    // would be a linking of input arguments
    // and tasks so that assigning an input argument
    // decrements a dependency count in the task,
    // and a zero count triggers an automatic move
    // to the ready queue.  not hard to do ... but
    // since it does not change the public nterface
    // we can do it later.

    /// The minimal interface required for all tasks in the task queue
    class TaskInterface {
    public:
        virtual bool probe() const = 0;
        virtual void run() = 0;
    };

    /// A task with one input argument.

    /// The operation can be a function or functor with this interface
    /// \code
    ///      void op(const T& in);
    /// \endcode
    template <typename OpT, typename T>
    class Task1in : public TaskInterface {
    private:
        OpT op;
        T arg;
    public:
        Task1in(OpT op, const T& arg) : op(op), arg(arg) {};

        bool probe() const {
            return arg.probe();
        };

        void run() {
            op(arg.get());
        };
    };

    /// A task with one output argument.

    /// The operation can be a function or functor with this interface
    /// \code
    ///      void op(SAV<T>& out);
    /// \endcode
    template <typename OpT, typename T>
    class Task1out : public TaskInterface {
    private:
        OpT op;
        T arg;
    public:
        Task1out(OpT op, T& arg) : op(op), arg(arg) {};

        bool probe() const {
            return true;
        };

        void run() {
            op(arg);
        };
    };

    /// A task with one input and one output argument.

    /// The operation can be a function or functor with this interface
    /// \code
    ///      void op(const inT& in, SAV<outT>& out);
    /// \endcode
    template <typename OpT, typename inT, typename outT>
    class Task1in1out : public TaskInterface {
    private:
        OpT op;
        inT inarg;
        outT outarg;
    public:
        Task1in1out(OpT op, inT& inarg, outT& outarg)
                : op(op), inarg(inarg), outarg(outarg) {};

        bool probe() const {
            return inarg.probe();
        };

        void run() {
            op(inarg.get(),outarg);
        };
    };

    /// A task with fixed-size arrays of input and output arguments

    /// The operation can be a function or functor with this interface
    /// \code
    ///      void op(const SAV<inT> in[n], SAV<outT> out[m]);
    /// \endcode
    template <typename OpT, typename inT, std::size_t n, typename outT, std::size_t m>
    class TaskAinAout : public TaskInterface {
    private:
        OpT op;
        inT inarg[n];
        outT outarg[m];
    public:
        TaskAinAout(OpT op, inT in[n], outT out[m]) : op(op) {
            for (int i=0; i<n; i++) inarg[i] = in[i];
            for (int i=0; i<m; i++)outarg[i] = out[i];
        };

        bool probe() const {
            for (int i=0; i<n; i++) if (!inarg[i].probe()) return false;
            return true;
        };

        void run() {
            op(inarg,outarg);
        };
    };

    class TaskQueue {
    private:
        // Is there a more efficient choice than list?
        std::list< TaskInterface* > ready;
        std::list< TaskInterface* > pending;

    public:
        /// Add a new task ... the task queue takes ownership of the pointer.
        void add(TaskInterface* t) {
            if (t->probe()) ready.push_back(t);
            else pending.push_back(t);
        };

        /// Probe pending tasks and move the first ready one to ready queue.

        /// Returns true if a ready task was found, false otherwise.
        inline bool probe() {
            if (!pending.empty()) {
                for (std::list<TaskInterface *>::iterator p = pending.begin();
                        p != pending.end(); ++p) {
                    TaskInterface *tp = *p;
                    if (tp->probe()) {
                        ready.push_back(tp);
                        pending.erase(p);
                        return true;
                    }
                }
            }
            return false;
        };

        /// Runs the next ready task if there is one, returns true if one was run
        inline bool run_next_ready_task() {
            if (ready.empty()) {
                return false;
            }
            else {
                TaskInterface *p = ready.front();
                p->run();
                ready.pop_front();
                return true;
            }
            return false;
        };

        // Need a probe all to ensure progress of multistep stuff

        // Need hooks to permit use of MPI_Testany instead of busy wait

        /// Runs until all tasks have been completed
        void wait() {
            while (!(pending.empty() && ready.empty())) {
                probe();
                if (!run_next_ready_task()) yield();
            }
        };
        
        void yield() {
            usleep(10);
        };
       };

    extern TaskQueue globalq;  // In tasks.cc
}
#endif /*TASKS_H_*/
