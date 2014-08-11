#include <stdio.h>
#include <assert.h>

#include <sys/time.h>
#include <sys/types.h>
#include <sys/sysctl.h>

int main(int argc, char * argv[])
{
    {
        int mib[2], maxproc;
        size_t len;
        mib[0] = CTL_KERN;
        mib[1] = KERN_MAXPROC;
        len = sizeof(maxproc);
        sysctl(mib, 2, &maxproc, &len, NULL, 0);
        printf("maxproc = %d \n", maxproc);
    }
    {
        int ncpu;
        size_t len = sizeof(ncpu);
        int mib[2] = {CTL_HW,HW_NCPU};
        //mib[0] = CTL_HW;
        //mib[1] = HW_NCPU;
        sysctl(mib, 2, &ncpu, &len, NULL, 0);
        printf("ncpu = %d \n", ncpu);
    }
    {
        struct clockinfo c;
        int mib[2] = {CTL_KERN, KERN_CLOCKRATE};
        size_t len = sizeof(c);
        sysctl(mib, 2, &c, &len, NULL, 0);
        printf("c.hz = %d \n", c.hz);
    }
    {
        long freq;
        size_t len = sizeof(freq);
        sysctlbyname("hw.cpufrequency", &freq, &len, NULL, 0);
        printf("freq = %ld \n", freq);
    }
    return 0;
}
