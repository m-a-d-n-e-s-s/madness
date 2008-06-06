#include <unistd.h>
#include <cstdio>

namespace madness {
    
    /// Rename gmon.out for each process by ordering process termination
    
    /// Invoke with id and nproc as rank and size in MPI::COMM_WORLD
    void gprofexit(int id, int nproc) {
        char buf[64];
        if (id == 0) {
            for (int p=nproc-1; p>0; p--) {
                while (access("gmon.out",F_OK)) usleep(1000);
                sprintf(buf,"gmon.out.%d",id+1);
                if (rename("gmon.out",buf)) 
                    fprintf(stderr,"gprofexit: failed renaming gmon.out to %s", buf);
            }
        }
        else if ((id+1) < nproc) {
            // Wait for process id+1 to commence writing
            sprintf(buf,"gmon.out.%d",id+1);
            while (access(buf,F_OK)) usleep(10000);
        }
    }

}
