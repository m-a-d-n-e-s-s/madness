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

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/constants.h>
#include <madness/mra/FuseT/InnerOp.h>
#include <madness/mra/FuseT/AddOp.h>
#include <madness/mra/FuseT/MultiplyOp.h>
#include <madness/mra/FuseT/CompressOp.h>
#include <madness/mra/FuseT/FuseT.h>
#include <madness/mra/FuseT/FusedExecutor.h>
//#include <madness/mra/FuseT/OpExecutor.h>

using namespace madness;;

static const double L		= 20;
static const long	k		= 8;
static const double thresh	= 1e-6;
static const double c		= 2.0;
static const double alpha	= 1.9; // Exponent

#define PI 3.1415926535897932385
#define LO 0.0000000000
#define HI 4.0000000000

static double sin_amp		= 1.0;
static double cos_amp		= 1.0;
static double sin_freq		= 1.0;
static double cos_freq		= 1.0;
static double sigma_x		= 1.0;
static double sigma_y		= 1.0;
static double sigma_z		= 1.0;
static double center_x		= 0.0;
static double center_y		= 0.0;
static double center_z		= 0.0;
static double gaussian_amp	= 1.0;
static double sigma_sq_x	= sigma_x*sigma_x;
static double sigma_sq_y	= sigma_y*sigma_y;
static double sigma_sq_z	= sigma_z*sigma_z;

double rtclock();

static double uinitial(const coord_3d& r) 
{
	const double x=r[0], y=r[1], z=r[2];
    return -2.0/(sqrt(x*x+y*y+z*z+1e-8));
}

static double random_function(const coord_3d& r) 
{
	const double x=r[0], y=r[1], z=r[2];

	const double dx = x - center_x;
	const double dy = y - center_y;
	const double dz = z - center_z;

	const double periodic_part = sin_amp * sin(sin_freq*(dx+dy+dz)) 
									+ cos_amp * cos(cos_freq*(dx+dy+dz));

	const double x_comp = dx*dx/sigma_sq_x;
	const double y_comp = dy*dy/sigma_sq_y;
	const double z_comp = dz*dz/sigma_sq_z;

	const double gaussian_part = -gaussian_amp/exp(sqrt(x_comp+y_comp+z_comp));

	return gaussian_part*gaussian_part;
}

static double get_rand() {
	double r3 = LO + static_cast<double>(rand())/(static_cast<double>(RAND_MAX/(HI-LO)));
	return r3;
}

static void randomizer()
{
	sin_amp = get_rand();
	cos_amp = get_rand();
	sin_freq = get_rand();
	cos_freq = get_rand();
	sigma_x = get_rand();
	sigma_y = get_rand();
	sigma_z = get_rand();
	center_x = get_rand()*L/(2.0*HI);
	center_y = get_rand()*L/(2.0*HI);
	center_z = get_rand()*L/(2.0*HI);
	gaussian_amp = get_rand();
	sigma_sq_x = sigma_x*sigma_x;
	sigma_sq_y = sigma_y*sigma_y;
	sigma_sq_z = sigma_z*sigma_z;
}

int main(int argc, char** argv) 
{
  initialize(argc,argv);
  World world(SafeMPI::COMM_WORLD);

  startup(world, argc, argv);

  FunctionDefaults<3>::set_k(k);
  FunctionDefaults<3>::set_thresh(thresh);
  FunctionDefaults<3>::set_refine(true);
  FunctionDefaults<3>::set_autorefine(false);
  FunctionDefaults<3>::set_cubic_cell(-L, L);
  FunctionDefaults<3>::set_max_refine_level(14);

  real_function_3d f1 = real_factory_3d(world).f(random_function);
  real_function_3d f2 = real_factory_3d(world).f(random_function);
  real_function_3d f3 = real_factory_3d(world).f(random_function);
  real_function_3d f4 = real_factory_3d(world).f(random_function);
  real_function_3d f5 = real_factory_3d(world).f(random_function);
	
  real_function_3d i1_factory = real_factory_3d(world);
  real_function_3d i2_factory = real_factory_3d(world);
  real_function_3d i3_factory = real_factory_3d(world);
  real_function_3d i4_factory = real_factory_3d(world);
  real_function_3d i5_factory = real_factory_3d(world);
  real_function_3d i6_factory = real_factory_3d(world);
  real_function_3d i7_factory = real_factory_3d(world);
  real_function_3d i8_factory = real_factory_3d(world);
  real_function_3d i9_factory = real_factory_3d(world);
  real_function_3d i1(i1_factory);
  real_function_3d i2(i2_factory);
  real_function_3d i3(i3_factory);
  real_function_3d i4(i4_factory);
  real_function_3d i5(i5_factory);
  real_function_3d i6(i6_factory);
  real_function_3d i7(i7_factory);
  real_function_3d i8(i8_factory);
  real_function_3d i9(i9_factory);

  f1.truncate();
  f2.truncate();
  f3.truncate();
  f4.truncate();
  f5.truncate();

// i1 = f1 * f2;
// i4 = f2 * f3;
// compress--- i1, f3, f4, f5
// i3 = f5 - i1 // minus?
// i2 = f3 + f4
// double r = inner (i2, i3);

  double clkbegin, clkend;
  if (world.rank() == 0) print ("=====================================================");
  if (world.rank() == 0) print ("   FuseT                                             ");
  if (world.rank() == 0) print ("=====================================================");

  world.gop.fence();
  clkbegin = rtclock();

  MultiplyOp<double,3> mul_op_1("Multiply",&i1,&f1,&f2,0.0);	// i1
  MultiplyOp<double,3> mul_op_2("Multiply",&i4,&f2,&f3,0.0);	// i4
  CompressOp<double,3> compress_op_1("Compress",&i5,&i1);		// i5 <-- i1
  CompressOp<double,3> compress_op_2("Compress",&i6,&f3);		// i6 <-- f3
  CompressOp<double,3> compress_op_3("Compress",&i7,&f4);		// i7 <-- f4
  CompressOp<double,3> compress_op_4("Compress",&i8,&f5);		// i8 <-- f5
  AddOp<double,3>      add_op_1("Add",&i3,&i8,&i5);				// f5(i8) + i1(i5)
  AddOp<double,3>      add_op_2("Add",&i2,&i6,&i7);				// f3(i6) + f4(i7)
  InnerOp<double,3>    inner_op_1("Inner",&i9,&i2,&i3);

  if (world.rank() == 0) print ("==after Ops...=================================================");
  vector<PrimitiveOp<double,3>*> sequence;
  sequence.push_back(&mul_op_1);
  sequence.push_back(&mul_op_2);
  sequence.push_back(&compress_op_1);
  sequence.push_back(&compress_op_2);
  sequence.push_back(&compress_op_3);
  sequence.push_back(&compress_op_4);
  sequence.push_back(&add_op_1);
  sequence.push_back(&add_op_2);
  sequence.push_back(&inner_op_1);

  FuseT<double,3>			odag(sequence);
  odag.processSequence();
  FusedOpSequence<double,3> fsequence = odag.getFusedOpSequence();
  FusedExecutor<double,3>	fexecuter(world, &fsequence);
  if (world.rank() == 0) print ("==before exe================================================");
  fexecuter.execute();

  clkend = rtclock() - clkbegin;
  if (world.rank() == 0) printf("Running Time: %f\n", clkend);

  if (world.rank() == 0) print ("=====================================================");
  if (world.rank() == 0) print ("   MADNESS                                           ");
  if (world.rank() == 0) print ("=====================================================");
  world.gop.fence();

  real_function_3d j1_factory = real_factory_3d(world);
  real_function_3d j2_factory = real_factory_3d(world);
  real_function_3d j3_factory = real_factory_3d(world);
  real_function_3d j4_factory = real_factory_3d(world);
  real_function_3d j1(j1_factory);
  real_function_3d j2(j4_factory);
  real_function_3d j3(j1_factory);
  real_function_3d j4(j4_factory);

  clkbegin = rtclock();

  j1 = f1 * f2;
  j4 = f2 * f3;

  j1.compress();
  f3.compress();
  f4.compress();
  f5.compress();
  
  j3 = f5 + j1;
  j2 = f3 + f4;

  double r = j2.inner(j3);

  clkend = rtclock() - clkbegin;
  if (world.rank() == 0) printf("Running Time: %f\n", clkend);

  if (world.rank() == 0) printf("[MADNESS]r: %f, [FuseT]i9: %f\n",r, inner_op_1._sum);
  world.gop.fence();

  finalize();
  return 0;
}

double rtclock()
{
	struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0)
	printf("Error return from gettimeofday: %d", stat);
    return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

