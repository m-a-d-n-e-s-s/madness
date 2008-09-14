#include <iostream>
#include <pthread.h>
#include <unistd.h>
#include <stdlib.h>

#include <world/worldthread.h>

using namespace std;
using namespace madness;

MADATOMIC_INT sum = MADATOMIC_INITIALIZE(0);
MADATOMIC_INT ndone = MADATOMIC_INITIALIZE(0);

void* doit(void *args) {
    for (int j=0; j<1000; j++) {
        for (int i=0; i<10000; i++) {
            MADATOMIC_INT_INC(&sum);
        }
        sched_yield();
    }
    MADATOMIC_INT_INC(&ndone);

    return 0;
}

void* doit_mutex(void *args) {
    Mutex* p = (Mutex *) args;
    for (int j=0; j<10; j++) {
        for (int i=0; i<100; i++) {
            p->lock();
            sum++;
            p->unlock();
        }
    }

    p->lock();
    ndone++;
    p->unlock();

    return 0;
}

/// Forward decls
void* pool_thread_main(void *v);
class Pool;

class TaskInterface {
public:
    TaskInterface() {}

    /// Actually execute the task
    virtual void run() = 0;
};



/// A singleton pool of threads for dynamic execution of tasks.
    
/// There is only one instance that is created the first time that you
/// obtain a pointer to it via \c Pool::instance().
class Pool {
private:
    struct pool_thread_args *thread_args;
    Thread *threads;            ///< Array of threads
    TaskInterface **tasks;		///< Circular buffer of pointers to tasks
    int nthreads;		///< No. of threads
    int maxtasks;		///< Max no. of tasks
    volatile int task_count;    ///< No. of pending tasks
    volatile int task_tail;     ///< Next task to do 
    volatile int task_head;	///< Place to add new task 

    Mutex task_mutex; ///< Protects all task info 
    
    static Pool* instance_ptr;
    
    Pool(int nthread, int maxtask);
    Pool(const Pool&) {};           // Verboten
    void operator=(const Pool&) {}; // Verboten

    int default_nthread();
    
    public:
        /// Return a pointer to the only instance constructing as necessary
        inline static Pool* instance(int nthread=-1, int maxtask=-1) {
            if (!instance_ptr) instance_ptr = new Pool(nthread, maxtask);
            return instance_ptr;
        };
        
    bool add_task(TaskInterface* task);
        
    ~Pool() {};
        
    friend void* pool_thread_main(void *);

    protected:
    void thread_main(Thread *thread);
};


Pool* Pool::instance_ptr = 0;

struct pool_thread_args {
    Pool *pool;
    Thread *thread;
};

/// The constructor is private to enforce the singleton model

/// Note that the main (invoking) thread is in addition to the
/// threads in the pool, so that setting POOL_NTHREAD to 0 will
/// make tasks sequentially executed by the main thread.
///
/// If nthreads < 0 use the value of the environment variable POOL_NTHREAD
/// .                or if that is not defined run with 0 threads.
/// If maxtasks < 0 use a default limit of 100000 outstanding tasks.
Pool::Pool(int nthread=-1, int maxtask=-1) : nthreads(nthread), maxtasks(maxtask) {
    int i;
    if (nthreads < 0) nthreads = default_nthread();
    if (maxtasks < 0) maxtasks = 100000;
    
    task_count = 0;
    task_head = 0;
    task_tail = 0;
    
    try {
        tasks = new TaskInterface*[maxtasks];
        if (nthreads > 0) threads = new Thread[nthreads];
        else threads = 0;
        if (nthreads > 0) thread_args = new struct pool_thread_args[nthreads];
        else thread_args = 0;
    }
    catch (...) {
        throw "memory allocation failed";
    }
    
    for (i=0; i<nthreads; i++) {
        threads[i].num = i;
        thread_args[i].pool = this;
        thread_args[i].thread = threads+i;
        threads[i].start(pool_thread_main, thread_args+i);
    }
    
    cout << "Pool: initialized with " << nthreads << " threads and " 
         << maxtasks << " tasks." << endl;
};

/// Get number of threads from the environment
int Pool::default_nthread() {
    int nthread;
    char *cnthread = getenv("POOL_NTHREAD");
    
    if (cnthread) {
        if (sscanf(cnthread, "%d", &nthread) != 1) 
            throw "POOL_NTHREAD is not an integer";
        else
            return nthread;
    }
    return 0;
};


/// Add a new task to the pool
bool Pool::add_task(TaskInterface* task) {
    task_mutex.lock();    // BEGIN CRITICAL SECTION
    
    if (task_head==task_tail && task_count>0) { // The queue is full 
        if (task_count != maxtasks) throw "task count confused in full queue";
        else throw "task queue is full";
    }
    
    tasks[task_head] = task;
    task_head = (task_head+1) % maxtasks;
    task_count++;
    task_mutex.unlock();  // END CRITICAL SECTION
    
    return true;
};

void Pool::thread_main(Thread *thread) {
    while (1) {
        TaskInterface *task = 0;

        task_mutex.lock();  // BEGIN CRITICAL SECTION
        if (task_count != 0) {
            if (task_count < 0) throw "task_count <= 0 inside worker thread";
        
            task_count--;
            task = tasks[task_tail];
            task_tail = (task_tail+1)%maxtasks;
        }
        task_mutex.unlock();  // END CRITICAL SECTION
        
        if (task) {
            task->run();          // What we are here to do
            delete task;
        }
    }
};


/// Little wrapper for passing bound member function as thread main program
void* pool_thread_main(void *v) {
    struct pool_thread_args *p = (struct pool_thread_args *) v;
    p->pool->thread_main(p->thread);
    return 0;
}

class Greet : public TaskInterface {
public:
    void run() {
        std::cout << "HI\n";
    }
};

class Adder : public TaskInterface {
public:
    void run() {
        MADATOMIC_INT_INC(&sum);
        if (sum >= 100000) MADATOMIC_INT_INC(&ndone);
    }
};

int main() {
    const int nthread = 4;
    Thread threads[nthread];

    sum = ndone = 0;
    for (int i=0; i<nthread; i++) threads[i].start(doit,0);
    while (MADATOMIC_INT_GET(&ndone) != nthread) sleep(1);
    cout << "SUM " << MADATOMIC_INT_GET(&sum) << endl;

    Mutex p;
    sum = ndone = 0;
    for (int i=0; i<nthread; i++) threads[i].start(doit_mutex,&p);
    while (MADATOMIC_INT_GET(&ndone) != nthread) sleep(1);
    cout << "SUM " << MADATOMIC_INT_GET(&sum) << endl;

    sum = ndone = 0;
    Pool *pool = Pool::instance();
    pool->add_task(new Greet());
    for (int i=0; i<100000; i++) pool->add_task(new Adder());
    while (!MADATOMIC_INT_GET(&ndone)) sleep(1);
    cout << "SUM " << MADATOMIC_INT_GET(&sum) << endl;

    return 0;
};
