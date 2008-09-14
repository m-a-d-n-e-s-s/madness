#ifndef MAD_WORLDTHREAD_H
#define MAD_WORLDTHREAD_H

#include <pthread.h>
#include <unistd.h>
#include <world/madatomic.h>

namespace madness {

#define WORLD_MUTEX_ATOMIC

#ifdef WORLD_MUTEX_ATOMIC
    /// Mutex using spin locks and atomic operations
    
    /// Should be much faster than pthread operations for infrequent busy
    /// waiting ... however, expect pathalogical slowdown for busy
    /// waiting if you oversubscribe the processors.  Also cannot use
    /// these Mutexes with pthread condition variables.
    class Mutex {
    private:
        MADATOMIC_INT flag;
        
        /// Copy constructor is forbidden
        Mutex(const Mutex& m) {}
        
        /// Assignment is forbidden
        void operator=(const Mutex& m) {}
        
    public:
        /// Make and initialize a mutex ... initial state is unlocked
        Mutex() {
            MADATOMIC_INT_SET(&flag,1L);
        }
        
        /// Acquire the mutex waiting if necessary
        inline void lock() {
            while (1) {
                if (MADATOMIC_INT_DEC_AND_TEST(&flag)) return;
                MADATOMIC_INT_INC(&flag);
                //if (yield) sched_yield();
            }
        }
        
        /// Free a mutex owned by this thread
        inline void unlock() {
            MADATOMIC_INT_INC(&flag);
        }
        
        inline pthread_mutex_t* ptr() {throw "there is no pthread_mutex";}
        
        ~Mutex() {};
    };
    
#else
    
    /// Mutex using pthread operations
    class Mutex {
    private:
        pthread_mutex_t mutex;
        
        /// Copy constructor is forbidden
        Mutex(const Mutex& m) {}
        
        /// Assignment is forbidden
        void operator=(const Mutex& m) {}
        
    public:
        /// Make and initialize a mutex ... initial state is unlocked
        Mutex() {pthread_mutex_init(&mutex, 0);}
        
        /// Acquire the mutex waiting if necessary
        inline void lock() {
            if (pthread_mutex_lock(&mutex)) throw "failed acquiring mutex";
        }
        
        /// Free a mutex owned by this thread
        inline void unlock() {
            if (pthread_mutex_unlock(&mutex)) throw "failed releasing mutex";
        }
        
        /// Return a pointer to the pthread mutex for use by a condition variable
        inline pthread_mutex_t* ptr() {return &mutex;}

        ~Mutex() {pthread_mutex_destroy(&mutex);};
    };
    
#endif

    
    /// Simplified thread wrapper to hide pthread complexity
    class Thread {
    public:
        pthread_t id;               ///< posix thread id
        int num;                    ///< madness integer thread id
        
        /// Default constructor ... must invoke \c start() to actually begin the thread.
        Thread() {};
        
        /// Create a thread and start it running f(args)
        Thread(void* (*f) (void *), void* args=0) {
            this->start(f,args);
        };
        
        /// Start a thread made with default constructor executing \c f(args)
        void start(void* (*f) (void *), void* args=0) {
            pthread_attr_t attr;
            // Want detached thread with kernel scheduling so can use multiple cpus
            // RJH ... why not checking success/failure????
            pthread_attr_init(&attr);    
            pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_DETACHED);
            pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM); 
            if (pthread_create(&id, &attr, f, args))
                throw "failed creating thread";
            
            pthread_attr_destroy(&attr);
        }
        
        /// A thread can call this to terminate its execution
        static void exit() {pthread_exit(0);}
    };
    
}

#endif
