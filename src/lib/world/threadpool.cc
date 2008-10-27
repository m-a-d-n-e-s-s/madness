#include <world/worldthread.h>

namespace madness {
    int ThreadBase::cpulo[3];
    int ThreadBase::cpuhi[3];
    bool ThreadBase::bind[3];
    
    ThreadPool* ThreadPool::instance_ptr = 0;
}

