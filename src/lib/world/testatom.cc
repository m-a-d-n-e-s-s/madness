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



class Greet : public PoolTaskInterface {
public:
    void run() {
        std::cout << "HI\n";
    }
};

class Adder : public PoolTaskInterface {
public:
    void run() {
        MADATOMIC_INT_INC(&sum);
        if (sum >= 100000000) MADATOMIC_INT_INC(&ndone);
    }
};

int main() {
    const int nthread = 4;
    Thread threads[nthread];

    try {
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
        ThreadPool::add(new Greet());
        for (int i=0; i<100000000; i++) ThreadPool::add(new Adder());
        while (!MADATOMIC_INT_GET(&ndone)) sleep(1);
        cout << "SUM " << MADATOMIC_INT_GET(&sum) << endl;
    }
    catch (const char* e) {
        cout << "string exception: " << e << endl;
    }
    catch (char const* e) {
        cout << "string exception: " << e << endl;
    }

    return 0;
};
