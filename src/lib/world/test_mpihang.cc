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


#include <mpi.h>
#include <iostream>

using namespace std;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int np = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    if (np != 2) MADNESS_EXCEPTION("2 only", np);

    int me = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    int other = me? 0 : 1;

    int a=0, b=-1;
    MPI_Request rsend;
    MPI_Isend(&a, sizeof(a), MPI_BYTE, other, 1, MPI_COMM_WORLD, &rsend);
    MPI_Request rrecv;
    MPI_Irecv(&b, sizeof(b), MPI_BYTE, other, 1, MPI_COMM_WORLD, &rrecv);

    MPI_Status status;
    int done = 0;
    while (!done) MPI_Request_get_status(rsend, &done, &status);
    while (!done) MPI_Request_get_status(rrecv, &done, &status);
    MPI_Test(&rsend, &done, &status);
    MPI_Test(&rrecv, &done, &status);

    cout << me << " got " << b << endl;

    MPI_Finalize();
    return 0;
}
