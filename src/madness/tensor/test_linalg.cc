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

  

/// \file ../tensor/test.cc
/// \brief Test code for LAPACK, Tensor+LAPACK, etc.

#include <madness/tensor/tensor_lapack.h>
#include <iostream>
#include <madness/madness_config.h>

using namespace madness;



int
main(int argc, char* argv[]) {

//vama    const int required = MADNESS_MPI_THREAD_LEVEL;
//vama
//vama
//vama#ifdef MADNESS_HAS_ELEMENTAL
//vama    SafeMPI::Init_thread(argc, argv, required);
//vama    const int myrank = mpi::CommRank( mpi::COMM_WORLD );
//vama#else
//vama    const int myrank = 0;
//vama#endif

    bool testok = test_tensor_lapack();
    std::cout << "Test " << (testok ? "passed" : "did not pass") << std::endl;

    return int(!testok);
}

