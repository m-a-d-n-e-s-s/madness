#include <iostream>
#include <pthread.h>
#include <unistd.h>

using namespace std;

#include <world/madatomic.h>


/// Mutex using pthread operations
class Mutex {
private:
    pthread_mutex_t mutex;
    
    /// Copy constructor is forbidden
    Mutex(const Mutex& m) {};

    /// Assignment is forbidden
    void operator=(const Mutex& m) {};
    
public:
    /// Make and initialize a mutex ... initial state is unlocked
    Mutex() {pthread_mutex_init(&mutex, 0);}
    
    /// Acquire the mutex waiting if necessary
    inline void lock() {
        if (pthread_mutex_lock(&mutex)) throw "failed acquiring mutex";
    };
    
    /// Free a mutex owned by this thread
    inline void unlock() {
        if (pthread_mutex_unlock(&mutex)) throw "failed releasing mutex";
    };
    
    /// Return a pointer to the pthread mutex for use by a condition variable
    inline pthread_mutex_t* ptr() {return &mutex;};
    
    ~Mutex() {pthread_mutex_destroy(&mutex);};
};

#ifdef MUTEX_ATOMIC
/// Mutex using spin locks and atomic operations

/// This is provided for the Cray and other machines on which 
/// the pthread operations are desperately slow ... however, 
/// expect pathalogical slow down if you oversubscribe the processors.
/// Also cannot use these Mutexes with pthread condition variables
class MutexAtomic {
private:
    // Can add padding before and after to avoid false
    // sharing on cache-line
    MADATOMIC_INT flag;
    
    /// Copy constructor is forbidden
    MutexAtomic(const Mutex& m) {};

    /// Assignment is forbidden
    void operator=(const Mutex& m) {};
    
public:
    /// Make and initialize a mutex ... initial state is unlocked
    MutexAtomic() : flag(1) {MADATOMIC_FENCE;};
    
    /// Acquire the mutex waiting if necessary
    inline void lock() {
        while (1) {
            if (MADATOMIC_INT_COMPARE_AND_SWAP(&flag,1L,0L)) return;
            backoff(64);
            //sched_yield();
        }
        //              while (1) {
        //                if (MADATOMIC_INT_DEC_AND_TEST(&flag)) return;
        //               MADATOMIC_INT_INC(&flag);
        //           }
    };
    
    /// Free a mutex owned by this thread
    inline void unlock() {MADATOMIC_INT_SET(&flag,1L);};
    //        inline void unlock() {MADATOMIC_INT_INC(&flag);};
    
    ~MutexAtomic() {};
};
#endif

/// Simplified wrapper to hide pthread junk
class Thread {
public:
    pthread_t id;               ///< posix thread id
    int num;                    ///< integer thread id
    
    /// Default constructor ... must invoke \c start() to actually begin the thread.
    Thread() {};
    
    /// Create a thread and start it running f(args)
    Thread(void* (*f) (void *), void* args=0) {
        this->start(f,args);
    };
    
    /// Start a thread made with default constructor executing \c f(args)
    void start(void* (*f) (void *), void* args) {
        pthread_attr_t attr;
        // Want detached thread with kernel scheduling (so can use multiple cpus)
        pthread_attr_init(&attr);    
        // pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_DETACHED);
        // pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM); 
        if (pthread_create(&id, &attr, f, args))
            throw "failed creating thread";

        pthread_attr_destroy(&attr);
    }
        
    /// A thread can call this to terminate its execution
    static void exit() {pthread_exit(0);}
};

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

int main() {
    const int nthread = 4;
    Thread threads[nthread];

    for (int i=0; i<nthread; i++) threads[i].start(doit,0);

    while (MADATOMIC_INT_GET(&ndone) != nthread) sleep(1);

    cout << "SUM " << MADATOMIC_INT_GET(&sum) << endl;

    return 0;
};
