#include <iostream>
#include <unistd.h>
#include <pthread.h>
#include <papi.h>

using namespace std;

const int NUMEVENTS = 1;
int events[NUMEVENTS] = {PAPI_FP_OPS};

long long values1[NUMEVENTS];

void* f1(void* p) {
    if (PAPI_start_counters(events, NUMEVENTS) != PAPI_OK)
        throw "Could not start papi counters";

    double sum = 0.0;
    for (int i=0; i<100000; i++) {
        sum += 3.0;
    }

    if (PAPI_stop_counters(values1, NUMEVENTS) != PAPI_OK)
        throw "Could not stop papi counters";

    double* dp = (double*)(p);
    *dp = sum;
    return 0;
}

int main() {
    if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT)
        throw "Could not init PAPI";
    if (PAPI_thread_init(pthread_self) != PAPI_OK)
        throw "Could not init PAPI thread API";
    
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
    
    double sum1;
    pthread_t t1;
    pthread_create(&t1, &attr, f1, (void *)(&sum1));
    usleep(1000000);
    cout << "sum1 " << sum1 << " " << values1[0] << endl;
    return 0;
}
