#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <misc/communicator.h>

using namespace std;
namespace madness {
    
    void xterm_debug(const Communicator& comm, 
		     const char* path, const char* display)
    {
	pid_t child;
	const char *argv[20], *xterm="/usr/bin/xterm";
	char title[256], pid[256];
	
	sprintf(title, "Debugging process %d ", comm.rank());
	sprintf(pid, "%d", getpid());
	
	if (path == 0) path = "test1";
	if (display == 0) display = getenv("DISPLAY");
	if (display == 0) return;
	
	argv[0] = xterm;
	argv[1] = "-T";
	argv[2] = title;
	argv[3] = "-display";
	argv[4] = display;
	argv[5] = "-e";
	argv[6] = "gdb";
	argv[7] = path;
	argv[8] = pid;
	argv[9] = 0;
	if (comm.rank() == 0) {
	    int i;
	    printf("\n Starting xterms with debugger using command\n\n    ");
	    for (i=0; argv[i]; i++) printf("'%s' ", argv[i]);
	    printf("\n\n");
	    fflush(stdout);
	}
	
	child = fork();
	
	if (child < 0) {
	    printf("debug: fork failed?\n\n");
	}
	else if (child > 0) {
	    sleep(5);			/* Release cpu while debugger starts*/
	}
	else {
	    execv(xterm, (char*const*) argv);
	    perror("");
	    printf("util_debug: execv of xterm with debugger failed\n\n");
	}
    }
    
}

