/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


  $Id$
*/

#include <madness/madness_config.h>
#include <madness/world/safempi.h>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>

using namespace std;
namespace madness {

    extern "C" void xterm_debug_breakpoint() {
        std::cout << "xterm_debug_breakpoint" << std::endl;
    }

#if defined(HAVE_XTERM) && defined(HAVE_FORK) && (defined(HAVE_GDB) || defined(HAVE_LLDB)) && defined(HAVE_SLEEP)
    void xterm_debug(const char* path, const char* display) {
        int rank = SafeMPI::COMM_WORLD.Get_rank();
        pid_t child;
        const char *argv[20];
        const char xterm[] = XTERM_EXECUTABLE;
        const std::size_t bufsize = 255;
        char title[bufsize], pid[bufsize], geometry[bufsize];
        int ix=(rank/3)%3;
        int iy=rank%3;
        snprintf(title,bufsize,title, "Debugging process %d ", rank);
        snprintf(pid,bufsize, "%d", getpid());
        snprintf(geometry,bufsize,"%dx%d+%d+%d",80,24,ix*500,iy*280);

        if (path == 0) path = "test1";
        if (display == 0) display = getenv("DISPLAY");
        if (display == 0) return ;

        int nargs = 0;
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
#if defined(HAVE_GDB)
        argv[10] = GDB_EXECUTABLE;
        argv[11] = "-q";
        argv[12] = path;
        argv[13] = pid;
        nargs = 14;
#elif defined(HAVE_LLDB)
        argv[10] = LLDB_EXECUTABLE;
        argv[11] = "-p";
        argv[12] = pid;
        nargs = 13;
#endif
        argv[nargs] = 0;
        if (rank == 0) {
            int i;
            printf("\n Starting xterms with debugger using command\n\n    ");
            for (i = 0; argv[i]; ++i) printf("%s ", argv[i]);
            printf("\n\n");
            fflush(stdout);
        }

        child = fork();

        if (child < 0) {
            printf("debug: fork failed?\n\n");
        }
        else if (child > 0) {
            sleep(20);			/* Release cpu while debugger starts*/
            xterm_debug_breakpoint();
        }
        else {
            execv(xterm, (char*const*) argv);
            perror("");
            printf("util_debug: execv of xterm with debugger failed\n\n");
        }
    }

#else
    void xterm_debug(const char* path, const char* display) {}
#endif
}

