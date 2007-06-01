#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#ifdef SEEK_CUR
#undef SEEK_CUR
#endif
#ifdef SEEK_SET
#undef SEEK_SET
#endif
#ifdef SEEK_END
#undef SEEK_END
#endif
#include <mpi.h>

using namespace std;
namespace madness {

    extern "C" void xterm_debug_breakpoint() {
        std::cout << "xterm_debug_breakpoint" << std::endl;
    }

#ifndef HAS_XTERM_DEBUG
    void xterm_debug(const char* path, const char* display) {}
#else
    void xterm_debug(const char* path, const char* display) {
        int rank = MPI::COMM_WORLD.Get_rank();
        pid_t child;
        const char *argv[20], *xterm = "/usr/bin/xterm";
        char title[256], pid[256], geometry[256];
        int ix=(rank/3)%3;
        int iy=rank%3;
        sprintf(title, "Debugging process %d ", rank);
        sprintf(pid, "%d", getpid());
        sprintf(geometry,"%dx%d+%d+%d",80,24,ix*500,iy*280);

        if (path == 0) path = "test1";
        if (display == 0) display = getenv("DISPLAY");
        if (display == 0) return ;

        argv[0] = xterm;
        argv[1] = "-T";
        argv[2] = title;
        argv[3] = "-display";
        argv[4] = display;
        argv[5] = "-fn";
        argv[6] = "6x10";
        argv[7] = "-geometry";
        argv[8] = geometry;
        argv[9] = "-e";
        argv[10] = "gdb";
        argv[11] = "-q";
        argv[12] = path;
        argv[13] = pid;
        argv[14] = 0;
        if (rank == 0) {
            int i;
            printf("\n Starting xterms with debugger using command\n\n    ");
            for (i = 0; argv[i]; i++) printf("%s ", argv[i]);
            printf("\n\n");
            fflush(stdout);
        }

        child = fork();

        if (child < 0) {
            printf("debug: fork failed?\n\n");
        } else if (child > 0) {
            sleep(20);			/* Release cpu while debugger starts*/
            xterm_debug_breakpoint();
        } else {
            execv(xterm, (char*const*) argv);
            perror("");
            printf("util_debug: execv of xterm with debugger failed\n\n");
        }
    }
#endif
}

