#ifndef TASKS_H_
#define TASKS_H_

/// \file tasks.h
/// \brief Implements TaskQueue, TaskInterface and other Task classes
#include <list>
#include <unistd.h>

#include <serialize/vecar.h>
using madness::archive::VectorInputArchive;
using madness::archive::VectorOutputArchive;


#include <misc/communicator.h>
// for default communicator and AM functionality


namespace madness {
    static void do_local_fence();
    /// The minimal interface required for all tasks inserted into the task queue
    class TaskInterface {
    public:
        /// Returns true if the task is ready to be run
        virtual bool probe() const = 0;
        /// Runs the task
        virtual void run() = 0;
        /// Virtual base class destructor so that derived class destructor is called
        virtual ~TaskInterface() {};
    };
    
    void task_generic_handler(Communicator& comm, ProcessID src, const AMArg& arg);
    void task_generic_broadcast_handler(Communicator& comm, ProcessID src, const AMArg& arg);
    
    typedef void (*GenericOpT)(Communicator& comm, ProcessID src, VectorInputArchive& ar);
    
    /// A generic local/remotely created task with arguments serialized in vector archive
    class TaskGeneric : public TaskInterface {
        mutable std::vector<unsigned char> v;
        const ProcessID src;
        const int tag;
        GenericOpT op;
        mutable bool ready;
        mutable MPI::Request req;
    public:
        /// This constructor is invoked locally
        TaskGeneric(GenericOpT op, const std::vector<unsigned char>& v) 
        : v(v), src(madness::comm_default->rank()), tag(1), op(op), ready(true) 
        {};
    
        /// This constructor is invoked by AM handler in response to remote request
        TaskGeneric(size_t size, ProcessID src, int tag, GenericOpT op) 
        : v(size), src(src), tag(tag), op(op), ready(size<=0) {
            if (size>0) req = madness::comm_default->Irecv(&v[0],size,src,tag);
        };
        
        /// This constructor is invoked by AM handler in response to a broadcast task
        TaskGeneric(ProcessID root, GenericOpT op,  const std::vector<unsigned char>& v) 
        : v(v), src(root), tag(1), op(op), ready(true) 
        {};        
        
        bool probe() const {return (ready || (ready=req.Test()));};
        
        // Good optimizations would be to make a form of VectorInputArchive
        // that deserializes into references rather than values and
        // to enable Tensor to wrap some memory that it does not manage.
        // These would eliminate most data motion on the receiver end.
        void run() {
            //print((int)v[0],(int)v[1],(int)v[2],(int)v[3],(int)v[4],(int)v[5],(int)v[6],(int)v[7]);
            VectorInputArchive ar(v); op(*madness::comm_default,src,ar);};

        ~TaskGeneric() {};
    };
    
    /// A local/remotely created task with arguments via the "active message" argument
    class TaskAM : public TaskInterface {
        AMArg arg;
        const ProcessID src;
        am_handlerT op;
    public:
        /// This constructor is invoked locally
        TaskAM(am_handlerT op, const AMArg& arg)
        : arg(arg), src(madness::comm_default->rank()), op(op) 
        {};
    
        /// This constructor is invoked by AM handler in response to remote request
        TaskAM(ProcessID src, am_handlerT op, const AMArg& arg) 
        : arg(arg), src(src), op(op) 
        {};     
        
        bool probe() const {return true;};
        
        void run() {op(*madness::comm_default,src,arg);};

        ~TaskAM() {};
    };
    
    
    /// Queue to manage and run tasks.  Currently single threaded, but not for long.
    
    /// Note that even though the queue is presently single threaded, the invocation of
    /// am_poll() can still cause multiple entries through the task layer.  If modifying
    /// this be careful about keeping modifications to the data structure atomic
    /// w.r.t. am_operations.
    class TaskQueue {
    private:
        // Is there a more efficient choice than list?
        std::list< TaskInterface* > ready;
        std::list< TaskInterface* > pending;
        HandleManager<GenericOpT> hgen;

    public:
        TaskQueue() {};
        
        /// Add a new local task.  The task queue takes ownership of the pointer.
        void add_local(TaskInterface* t) {
            if (t->probe()) ready.push_back(t);
            else pending.push_back(t);
            //print(madness::comm_default->rank(),"added task",(void *) t,"ntask=",pending.size(),ready.size());
            madness::comm_default->am_poll();
        };
        
        /// Add a new local high-priority task.  The task queue takes ownership of the pointer.
        void add_local_hp(TaskInterface* t) {
            if (t->probe()) ready.push_front(t);
            else pending.push_front(t);
            //print(madness::comm_default->rank(),"added task",(void *) t,"ntask=",pending.size(),ready.size());
            madness::comm_default->am_poll();
        };
        
        /// Add a new local/remote task passing arguments via vector archive
        
        /// The routine op should be of type GenericOpT and registered in the TaskQueue
        void add_generic(ProcessID dest, GenericOpT op, const std::vector<unsigned char>& v) {
            if (dest < 0 || dest == madness::comm_default->rank()) {
                add_local(new TaskGeneric(op, v));
            }
            else {
                long handle = hgen.get_handle(op);
                // While we are single threaded and wait for the data send to complete,
                // a single tag is OK.  Otherwise, we must use unique tags for each task.
                madness::comm_default->am_send(dest, task_generic_handler, AMArg(v.size(),
                     TASK_GENERIC_DATA_TAG, handle));
                if (v.size()) {
                    MPI::Request req = madness::comm_default->Isend(&v[0],v.size(),
                        MPI::BYTE, dest, TASK_GENERIC_DATA_TAG);
                    madness::comm_default->am_wait(req);
                }
            }
        };
        
        /// Add a new local/remote task calling a routine with the AM handler interface
        
        /// The routine op should be of type am_handlerT and registered in the communicator.
        void add_am(ProcessID dest, am_handlerT op, const AMArg& arg) {
            if (dest < 0 || dest == madness::comm_default->rank()) {
                add_local(new TaskAM(op, arg));
            }
            else {
                long handle = madness::comm_default->_am_handle_manager.get_handle(op);
                madness::comm_default->am_send(dest, handle, arg, false, -1);
            }
        };


        /// Broadcast a generic task to all other nodes, \em excluding the invoking node
        void add_generic_broadcast(GenericOpT op, const std::vector<unsigned char>& v) {
            add_generic_broadcast(hgen.get_handle(op),v);
        };
        
        /// Broadcast an AM handler interface task to all other nodes, \em excluding the invoking node
        void add_am_broadcast(am_handlerT op, const AMArg& arg) {
            long handle = madness::comm_default->_am_handle_manager.get_handle(op);
            madness::comm_default->am_broadcast(handle, arg, false, madness::comm_default->rank());
        };
           
        /// Broadcast a generic task to all other nodes, \em excluding the invoking node
        
        /// The root argument is for internal use only.
        void add_generic_broadcast(long handle, const std::vector<unsigned char>& v, 
                        ProcessID root=-1) {
            // ?? why can this not just invoke am_broadcast and then send data?
            long nproc = madness::comm_default->nproc();
            long rank = madness::comm_default->rank();
            if (root == -1) root = rank;
            ProcessID me=(rank+nproc-root)%nproc;
            ProcessID child0=(me<<1)+1+root, child1=(me<<1)+2+root;
            if (child0 >= nproc && child0<(nproc+root)) child0 -= nproc;
            if (child1 >= nproc && child1<(nproc+root)) child1 -= nproc;

            MPI::Request req0, req1;
            AMArg arg(v.size(), TASK_GENERIC_DATA_TAG, handle, root);
            
            if (child0<nproc) {
                madness::comm_default->am_send(child0, task_generic_broadcast_handler,arg);
                if (v.size()) req0 = madness::comm_default->Isend(&v[0],v.size(),
                     MPI::BYTE, child0, TASK_GENERIC_DATA_TAG);
            }
            if (child1<nproc) {
                madness::comm_default->am_send(child1, task_generic_broadcast_handler,arg); 
                if (v.size()) req1 = madness::comm_default->Isend(&v[0],v.size(),
                     MPI::BYTE, child1, TASK_GENERIC_DATA_TAG);
            }
            if (v.size()) {
                if (child0<nproc) madness::comm_default->am_wait(req0);
                if (child1<nproc) madness::comm_default->am_wait(req1);
            }     
            //print("in bcast",(int)v[0],(int)v[1],(int)v[2],(int)v[3],(int)v[4],(int)v[5],(int)v[6],(int)v[7]);
            
        };

        /// Register user handler for local/remote generic task interface
        long register_generic_op(GenericOpT op) {
            return hgen.insert(op);
        };
        
        GenericOpT get_generic_pointer(long handle) {
            return hgen.get_pointer(handle);
        };
        
        /// Probe pending tasks and move the first ready one to ready queue.
        bool probe() {
            madness::comm_default->am_poll();
            if (!pending.empty()) {
                madness::comm_default->am_suspend();
                for (std::list<TaskInterface *>::iterator p = pending.begin();
                        p != pending.end(); ++p) {
                    TaskInterface *tp = *p;
                    if (tp->probe()) {                      
                        ready.push_back(tp);  // Must modify to preserve HP status
                        pending.erase(p);
                        madness::comm_default->am_resume();
                        return true;
                    }
                }
                madness::comm_default->am_resume();
            }
            return false;
        };
        
        /// Wait for an MPI request to complete while processing tasks and AM
        void wait(MPI::Request& req) {
            madness::comm_default->am_poll();
            while (!req.Test()) {
                probe();
                if (!run_next_ready_task()) yield();
            }
        };


        /// Returns after all local tasks have completed
        void local_fence() {
            madness::comm_default->am_poll();
            while (!(pending.empty() && ready.empty())) {
                probe();
                if (!run_next_ready_task()) yield();
            }
        };
        
        /// Returns after completion of all tasks & AM communicator wide (collective)
        void global_fence() {
            do {
                local_fence();
                yield();
                local_fence();
            } while (madness::comm_default->am_ndiff_spmd(bind_mem_fun(this,&TaskQueue::wait),do_local_fence));
        };    

        
    private:
        /// Runs the next ready task if there is one, returns true if one was run
        bool run_next_ready_task() {
            madness::comm_default->am_poll();
            if (ready.empty()) {
                return false;
            }
            else {
                madness::comm_default->am_suspend();
                TaskInterface *p = ready.front();
                ready.pop_front();
                //print(madness::comm_default->rank(),"executing task",(void *) p);
                p->run();
                delete p;
                madness::comm_default->am_resume();
                madness::comm_default->am_poll();
                return true;
            }
            return false;
        };

        // Need a probe all to ensure progress of multistep stuff

        // Need hooks to permit use of MPI_Testany instead of busy wait
        
        inline void yield() {
            usleep(10);
        };
    };

    extern TaskQueue taskq;  // Initialized in tasks.cc
    
    static void do_local_fence() {taskq.local_fence();};
    
}
#endif /*TASKS_H_*/
