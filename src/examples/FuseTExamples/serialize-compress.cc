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
*/


/// \file examples/hello.cc
/// \brief Simplest example program for MADNESS
/// \defgroup hellowworldmad Hello world MADNESS style
/// \ingroup examples
///
/// Simplest program that initializes the MADNESS parallel runtime
/// using initialize(), makes a madness::World object, prints
/// a greeting, and then cleans up.
///
/// To initialize the MADNESS numerical environment you also need
/// \c startup(world,argc,argv) and should include mra/mra.h rather
/// than world/MADworld.h .

#include <madness/world/MADworld.h>
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/constants.h>
#include <madness/mra/FuseT/CopyOp.h>
#include <madness/mra/FuseT/CompressOp.h>
#include <madness/mra/FuseT/NothingOp.h>
#include <madness/mra/FuseT/OpExecutor.h>
using namespace madness;;

static const double L		= 20;
static const long	k		= 8;
static const double thresh	= 1e-6;
static const double c		= 2.0;
static const double alpha	= 1.9; // Exponent

static double uinitial(const coord_3d& r) 
{
	const double x=r[0], y=r[1], z=r[2];
	return exp(-alpha*(2*x*x + y*y + z*z)) * pow(constants::pi/alpha, -1.5);
}

int main(int argc, char** argv) 
{
    initialize(argc,argv);
    World world(SafeMPI::COMM_WORLD);

    startup(world, argc, argv);

	real_function_3d u1 = real_factory_3d(world).f(uinitial);

    world.gop.fence();

	printf ("======================================\n");
	printf ("Before u1.compress()\n");
	u1.compress();

	world.gop.fence();

    if (world.rank() == 0) print("u0: norm", u1.norm2()," trace", u1.trace());

	//finalize();
    return 0;
}


