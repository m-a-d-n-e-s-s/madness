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

  

/// \file ../linalg/test.cc
/// \brief Test code for LAPACK, Tensor+LAPACK, etc.

#include <linalg/tensor_lapack.h>
#include <iostream>
#include <madness_config.h>

#ifdef MADNESS_HAS_ELEMENTAL
//#if 0
#include <mra/mra.h>
#include "elemental.hpp"
using namespace elem;
using namespace SafeMPI;
#endif

//#include <mra/mra.h>

using namespace madness;



int
main(int argc, char* argv[]) {

//vama#if defined(SERIALIZE_MPI) && !defined(MADNESS_HAS_ELEMENTAL)
//vama    int required = MPI_THREAD_SERIALIZED;
//vama#else
//vama    int required = MPI_THREAD_MULTIPLE;
//vama#endif
//vama
//vama
//vama#ifdef MADNESS_HAS_ELEMENTAL
//vama    SafeMPI::Init_thread(argc, argv, required);
//vama    const int myrank = mpi::CommRank( mpi::COMM_WORLD );
//vama#else
//vama    const int myrank = 0;
//vama#endif

#ifdef MADNESS_HAS_ELEMENTAL
    initialize( argc, argv );
    World world(SafeMPI::COMM_WORLD);
    const int myrank = world.rank();
#else
    const int myrank = 0;
#endif


    bool testok = test_tensor_lapack();
    if ( myrank==0 ) std::cout << "Test " << (testok ? "passed" : "did not pass") << std::endl;





#ifdef MADNESS_HAS_ELEMENTAL
    finalize();
#endif

    return 0;
}

