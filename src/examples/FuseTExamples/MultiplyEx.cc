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
#include <madness/mra/FuseT/ReconstructOp.h>
#include <madness/mra/FuseT/NothingOp.h>
#include <madness/mra/FuseT/FuseT.h>
#include <madness/mra/FuseT/FusedExecutor.h>
#include <madness/mra/FuseT/MultiplyOp.h>
#include <madness/mra/FuseT/OpExecutor.h>

using namespace madness;;

static const double L		= 20;
static const long	k		= 8;
static const double thresh	= 1e-6;
//static const double thresh	= 1e-12;
static const double c		= 2.0;
static const double alpha	= 1.9; // Exponent

inline static double uinitial(const coord_3d& r) 
{
  const double x=r[0], y=r[1], z=r[2];
    return -2.0/(sqrt(x*x+1.2*y*y+z*z+1e-8));

}
inline static double uinitial1(const coord_3d& r) 
{
  const double x=r[0], y=r[1], z=r[2];
    return -3.0/(sqrt(x*x/0.7+y*y*2.8+z*z*1.2+1e-8));

}
inline static double uinitial12(const coord_3d& r) 
{
  return uinitial(r)*uinitial1(r);
}

int main(int argc, char** argv) 
{
  initialize(argc,argv);
  World world(SafeMPI::COMM_WORLD);

  startup(world, argc, argv);

  FunctionDefaults<3>::set_defaults(world);
  FunctionDefaults<3>::set_k(k);
  FunctionDefaults<3>::set_thresh(thresh);
  FunctionDefaults<3>::set_refine(true);
  FunctionDefaults<3>::set_initial_level(5);
  FunctionDefaults<3>::set_truncate_mode(1);
  FunctionDefaults<3>::set_cubic_cell(-L/2, L/2);
  
  
  real_function_3d u0 = real_factory_3d(world).f(uinitial);
  real_function_3d u1 = real_factory_3d(world).f(uinitial1);
  real_function_3d u12 = real_factory_3d(world).f(uinitial12);
  //u0.truncate();
  //u1.truncate();

  double u12_norm	= u12.norm2();
  double u12_trace	= u12.trace();

  if (world.rank() == 0) print("[Analytical Product] Initial norm", u12_norm,"trace", u12_trace);

  real_function_3d result_factory = real_factory_3d(world);
  real_function_3d result_fuset(result_factory);
  real_function_3d result_factory1 = real_factory_3d(world);
  real_function_3d result_fuset1(result_factory);

  MultiplyOp<double,3> op1("Multiply",&result_fuset, &u0, &u1, 0.0);
  MultiplyOp<double,3> op2("Multiply",&result_fuset1, &u0, &u1, 0.0);
  OpExecutor<double,3> exe(world);
  exe.execute(&op1, false);

  vector<PrimitiveOp<double,3>*> sequence;
  sequence.push_back(&op2);
  FuseT<double,3> odag(sequence);
  FusedOpSequence<double,3> fsequence = odag.getFusedOpSequence();
  FusedExecutor<double,3> fexecuter(world, &fsequence);
  fexecuter.execute();


  double result_fuset_norm	= result_fuset.norm2();
  double result_fuset_trace	= result_fuset.trace();
  double result_fuset1_norm	= result_fuset1.norm2();
  double result_fuset1_trace= result_fuset1.trace();
  if (world.rank() == 0) print("[Result Fuset] Norm", result_fuset_norm,"trace", result_fuset_trace);
  if (world.rank() == 0) print("[Result Fuset] Norm", result_fuset1_norm,"trace", result_fuset1_trace);

  finalize();
  return 0;
}
