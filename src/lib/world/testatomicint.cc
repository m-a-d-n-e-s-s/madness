#include <world/worldthread.h>
#include <iostream>
#include <pthread.h>
#include <unistd.h>
#include <stdlib.h>

#include <world/atomicint.h>

using namespace std;
using namespace madness;

AtomicInt sum;
AtomicInt ndone;

void* doit(void *args) {
    for (int j=0; j<1000; j++) {
        for (int i=0; i<100000; i++) {
            sum++;
        }
        sched_yield();
    }
    ndone++;

    return 0;
}


class Greet : public PoolTaskInterface {
public:
    void run(const TaskThreadEnv& env) {
        std::cout << "HI\n";
    }

    
};

const int NDO = 10000000;

class Adder : public PoolTaskInterface {
public:
    void run(const TaskThreadEnv& env) {
        if (sum >= NDO) {
            ndone++;
        }
        else {
            sum++;
        }
    }
};

int main() {
    const int nthread = 3;
    Thread threads[nthread];

    try {
        sum = ndone = 0;
        for (int i=0; i<nthread; i++) threads[i].start(doit,0);
        while (ndone != nthread) sleep(1);
        cout << "SUM " << sum << endl;

        sum = ndone = 0;
        ThreadPool::add(new Greet());
        for (int i=0; i<(NDO+1000); i++) ThreadPool::add(new Adder());
        while (!ndone) {
            sleep(1);
            cout << int(ndone) << " " << int(sum) << endl;
        }
        sleep(1);
        cout << "SUM " << int(sum) << endl;
    }
    catch (const char* e) {
        cout << "string exception: " << e << endl;
    }
    catch (char const* e) {
        cout << "string exception: " << e << endl;
    }

    return 0;
};
