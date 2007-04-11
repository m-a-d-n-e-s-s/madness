#include <list>
#include <iostream>
#include <time.h>
#include <sys/time.h>

using namespace std;

double cpu_time() {
    return double(clock())/CLOCKS_PER_SEC;
}


int main() {
    double used;
    int** p = new int*[1000000];

    used = cpu_time();
    for (int i=0; i<1000000; i++) p[i] = new int(i);
    used = (cpu_time() - used);
    cout << "time to new 1M integers " << used << endl;

    used = cpu_time();
    for (int i=0; i<1000000; i++) delete p[i];
    used = (cpu_time() - used);
    cout << "time to del 1M integers " << used << endl;

        
    list<int> a;

    used = cpu_time();
    for (int i=0; i<1000000; i++) {
        a.push_back(i);
    }
    used = cpu_time() - used;
    cout << "time to push 1M integers in list " << used << endl;

    used = cpu_time();
    for (int i=0; i<1000000; i++) {
        a.pop_front();
    }
    used = cpu_time() - used;
    cout << "time to pop  1M integers in list " << used << endl;


    list<int*> b;

    used = cpu_time();
    for (int i=0; i<1000000; i++) {
        b.push_back(new int(i));
    }
    used = cpu_time() - used;
    cout << "time to new+push 1M integers in list " << used << endl;

    used = cpu_time();
    for (int i=0; i<1000000; i++) {
        delete b.front();
        b.pop_front();
    }
    used = cpu_time() - used;
    cout << "time to del+pop  1M integers in list " << used << endl;

    return 0;
}


