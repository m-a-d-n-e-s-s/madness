#include <madness_config.h>
#include <world/worldmutex.h>
#include <world/worldtime.h>
#include <iostream>
#include <vector>
#include <utility>

using namespace std;
using namespace madness;

template <typename T>
class SimpleQ : public Spinlock {
    const int len;
    volatile int ninq;
    T* q;

    SimpleQ(const SimpleQ& other); // Verboten

public:
    SimpleQ(int len=128) 
        : len(len)
        , ninq(0)
        , q(new T[len]) 
    {}

    ~SimpleQ() {
        delete q;
    }

    void push(const T& item) {
        lock();
        MADNESS_ASSERT(ninq <= len);
        q[ninq++] = item;
        unlock();
    }

    bool empty() const {
        return (ninq==0);
    }

    int size() const {
        return ninq;
    }

    int capacity() const {
        return len;
    }

    void clear() {
        ninq = 0;
    }

    const T& operator[](int i) const {
        return q[i];
    }

    // Atomically apply op() to all items in FIFO order and clear the queue
    template <typename opT>
    int process(const opT& op) {
        lock();
        int n = ninq;
        for (int i=0; i<ninq; i++) op(q[i]);
        clear();
        unlock();
        return n;
    }
};


template <typename T>
class DQueue {
    int n __attribute__((aligned(64)));        ///< Number of elements in the buffer
    int len;                                    ///< Current capacity
    T* buf;                                    ///< Actual buffer
    int front;               ///< Index of element at front of buffer
    int back;                ///< Index of element at back of buffer

    void grow() {
        if (len != n) MADNESS_EXCEPTION("assertion failure in dqueue::grow", len);
        int oldlen = len;
        if (len < 32768)
            len = 65536;
        else if (len <= 1048576)
            len *= 2;
        else
            len += 1048576;

        T* nbuf = new T[len];
        int lo = len/2 - oldlen/2;
        for (int i=front; i<int(oldlen); i++,lo++) {
            nbuf[lo] = buf[i];
        }
        if (front > 0) {
            for (int i=0; i<=back; i++,lo++) {
                nbuf[lo] = buf[i];
            }
        }
        front = len/2 - oldlen/2;
        back = front + n - 1;

        delete buf;

        buf = nbuf;
    }

    void sanity_check() const {
        int num = back - front + 1;
        if (num < 0) num += len;
        if (num==len && n==0) num=0;
        if (num==0 && n==len) num=len;
        //if (n != num) print("size",len,"front",front,"back",back,"n",n,"num",num);
        MADNESS_ASSERT(n == num);
    }

    DQueue(const DQueue& other); // Verboten
    
public:
    DQueue(int hint=32768) 
        : n(0)
        , len(hint>2 ? hint : 2)
        , buf(new T[len])
        , front(len/2)
        , back(front-1) {}
    
    virtual ~DQueue() {
        delete buf;
    }
    
    /// Insert value at front of queue
    void push_front(const T& value) {
        if (n == len) grow();
        n++;
        
        front--;
        if (front < 0) front = len - 1;
        buf[front] = value;
    }
    
    /// Insert element at back of queue
    void push_back(const T& value) {
        if (n == len) grow();
        n++;
        
        back++;
        if (back >= len) back = 0;
        buf[back] = value;
    }

    /// Pop value off the front of queue
    std::pair<bool,T> pop_front() {
        if (n) {
            n--;
            T result = buf[front++];
            if (front >= len) front = 0;
            return std::make_pair(true,result);
        }
        else {
            return std::make_pair(false,T());
        }
    }

    int size() const {
        return n;
    }

    int capacity() const {
        return len;
    }

    bool empty() const {
        return n==0;
    }
};

pthread_key_t idkey; // Must be external

static __thread int my_thread_id;

static inline void set_thread_id(int id) {
    //pthread_setspecific(idkey, (void*)(id));
    my_thread_id = id;
}

static inline int get_thread_id() {
    //long value = long(pthread_getspecific(idkey));
    //return int(value);
    return my_thread_id;
}

struct PoolTaskInterface {
    virtual void run() = 0;
};

struct InterThreadMessage {
    int n;
    PoolTaskInterface** buf;
    volatile int* count;

    InterThreadMessage() {}

    InterThreadMessage(int n, PoolTaskInterface** buf, volatile int* count) 
        : n(n), buf(buf), count(count)
    {}
};

class ThreadPool {

    struct ThreadPoolData {
        SimpleQ<InterThreadMessage> itmq;
        DQueue<PoolTaskInterface*> taskq;
        unsigned long ntask;
        
        ThreadPoolData() : ntask(0) {}
    };

    friend class PoolTaskInterface;
    typedef SimpleQ<InterThreadMessage> ITMQ;
    typedef DQueue<PoolTaskInterface*> TASKQ;

    static const int MAXNTHREAD = 32;
    static int nthread;
    static AtomicInt done;
    static ThreadPoolData* data[MAXNTHREAD];
    
    static void process_messages(ThreadPoolData* t) {
        ITMQ& itmq = t->itmq;
        TASKQ& taskq = t->taskq;
        itmq.lock();
        for (int i=0; i<itmq.size(); i++) {
            const InterThreadMessage& msg = itmq[i];
            int n = std::min(taskq.size()/2,msg.n);
            PoolTaskInterface** buf = msg.buf;
            for (int j=0; j<n; j++) buf[j] = taskq.pop_front().second;
            __asm__ __volatile__ ("" : : : "memory");
            *msg.count = n;
        }
        itmq.clear();
        itmq.unlock();
    }
    
    // Steal from target thread into taskq
    static void steal(int target, ThreadPoolData* t) {
        static const int n = 2048;
        PoolTaskInterface* buf[n];
        TASKQ& taskq = t->taskq;

        volatile int count = -1;
        __asm__ __volatile__ ("" : : : "memory");
        
        InterThreadMessage msg(n, buf, &count);
        data[target]->itmq.push(msg);
        while (count == -1) process_messages(t);

        int num = count;
        for (int i=0; i<num; i++) taskq.push_back(buf[i]);

        //std::cout << get_thread_id() << " from " << target << " got " << count << " " << taskq.size() << " " << t->ntask << endl;
    }

    static bool forever() {return false;}

    static void* main(void* args) {
        int id = long(args);
        set_thread_id(id);

        data[id] = new ThreadPoolData;

        done++;

        // Need to wait for everyone to be ready
        while (done != nthread); // Need to relax here

        work(forever);

        return 0;
    }

    static void start_thread(int id) {
        pthread_attr_t attr;
        // Want detached thread with kernel scheduling so can use multiple cpus
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_DETACHED);
#ifndef HAVE_IBMBGP	    
        pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
#endif
        pthread_t ptid;
        int result = pthread_create(&ptid, &attr, &ThreadPool::main, (void *)(id));
        if (result) MADNESS_EXCEPTION("failed creating thread", result);
        pthread_attr_destroy(&attr);
    }


    ThreadPool() {}; // Verboten;
    
public:

    static void initialize(int numthread) {
        MADNESS_ASSERT(numthread <= MAXNTHREAD);
        nthread = numthread;
        
        pthread_key_create(&idkey, NULL); // To store thread-specific id

        set_thread_id(0); // Main thread is zero
        data[0] = new ThreadPoolData;
        
        // Thread start up first binds thread to core, then allocates
        // queue structures to accomodate NUMA heap.  Done counts
        // threads that are ready
        done = 1;
        for (int i=1; i<nthread; i++) start_thread(i);
        while (done != nthread);
    }

    template <typename opT>
    static void work(const opT& terminated) {
        const int id = get_thread_id();
        ThreadPoolData* t = data[id];
        TASKQ& taskq = t->taskq;
        unsigned long& ntask = t->ntask;

        int target = id;
        while (!terminated()) {
            process_messages(t);

            std::pair<bool,PoolTaskInterface*> p = taskq.pop_front();
            if (p.first) {
                p.second->run();
                delete p.second;
                ntask++;
            }
            else if (nthread > 1) {
                do {
                    target++;
                    if (target >= nthread) target = 0;
                } while (target == id);
                steal(target, t);
            }
        }
    }

    static void add(PoolTaskInterface* task) {
        int id = get_thread_id();
        data[id]->taskq.push_back(task);
        process_messages(data[id]);
    }

    static int size() {return nthread;}

    static unsigned long ntaskdone(int id) {return data[id]->ntask;}
};

int ThreadPool::nthread;
ThreadPool::ThreadPoolData* ThreadPool::data[ThreadPool::MAXNTHREAD];
AtomicInt ThreadPool::done;

const int NGEN=100;
const int NTASK=100000;
AtomicInt cnt;

class Task : public PoolTaskInterface {
    const int gen;
    bool done;
public:
    Task(int gen) : gen(gen), done(false) {}

    void run() {
        MADNESS_ASSERT(!done);
        done = true;
        cnt++;
        if (gen > 1) ThreadPool::add(new Task(gen-1));
    }

    static bool finished() {return cnt==(NGEN*NTASK);}
};
    
int main() {
    SimpleQ<int> q;

    MADNESS_ASSERT(q.size() == 0 && q.empty());
    q.push(0); MADNESS_ASSERT(q.size() == 1 && !q.empty());
    q.push(1); MADNESS_ASSERT(q.size() == 2 && !q.empty());
    q.push(2); MADNESS_ASSERT(q.size() == 3 && !q.empty());

    q.lock();
    for (int i=0; i<q.size(); i++) MADNESS_ASSERT(q[i] == i);
    q.clear(); MADNESS_ASSERT(q.size() == 0 && q.empty());
    q.unlock();

    DQueue<int> dq;

    MADNESS_ASSERT(dq.size() == 0 && dq.empty());
    dq.push_back(3); MADNESS_ASSERT(dq.size() == 1 && !dq.empty());
    dq.push_back(4); MADNESS_ASSERT(dq.size() == 2 && !dq.empty());
    dq.push_back(5); MADNESS_ASSERT(dq.size() == 3 && !dq.empty());
    dq.push_front(2); MADNESS_ASSERT(dq.size() == 4 && !dq.empty());
    dq.push_front(1); MADNESS_ASSERT(dq.size() == 5 && !dq.empty());
    dq.push_front(0); MADNESS_ASSERT(dq.size() == 6 && !dq.empty());

    int j=0;
    while(1) {
        std::pair<bool,int> item = dq.pop_front();
        if (item.first) {
            MADNESS_ASSERT(item.second == j++);
            MADNESS_ASSERT(dq.size() == 6-j);
        }
        else {
            MADNESS_ASSERT(dq.size()==0 && dq.empty());
            break;
        }
    }

    ThreadPool::initialize(2);

    cnt = 0;
    for (int i=0; i<NTASK; i++)
        ThreadPool::add(new Task(NGEN));

    ThreadPool::work(Task::finished);

    for (int i=0; i<ThreadPool::size(); i++) cout << i << " " << ThreadPool::ntaskdone(i) << endl;

    return 0;
}
