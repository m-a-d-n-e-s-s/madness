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
#include <unistd.h>
#include <cstdio>

namespace madness {

    /// Rename gmon.out for each process by ordering process termination

    /// Invoke with id and nproc as rank and size in MPI::COMM_WORLD
    void gprofexit(int id, int nproc) {
        const std::size_t bufsize=64;
        char buf[bufsize];
        if (id == 0) {
            for (int p=nproc-1; p>0; --p) {
                while (access("gmon.out",F_OK)) usleep(1000);
                snprintf(buf,bufsize,"gmon.out.%d",id+1);
                if (rename("gmon.out",buf))
                    fprintf(stderr,"gprofexit: failed renaming gmon.out to %s", buf);
            }
        }
        else if ((id+1) < nproc) {
            // Wait for process id+1 to commence writing
            snprintf(buf,bufsize,"gmon.out.%d",id+1);
            while (access(buf,F_OK)) usleep(10000);
        }
    }

}
