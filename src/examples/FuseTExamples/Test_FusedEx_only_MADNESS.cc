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
//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/world/MADworld.h>
#include <madness/world/world_object.h>
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/constants.h>
#include <madness/mra/FuseT/CopyOp.h>
#include <madness/mra/FuseT/InnerOp.h>
#include <madness/mra/FuseT/CompressOp.h>
#include <madness/mra/FuseT/ReconstructOp.h>
#include <madness/mra/FuseT/FusedExecutor.h>
#include <madness/mra/FuseT/OpExecutor.h>

#include <madness/mra/FuseT/FuseT.h>

#define FUSET_UNFUSED
#define FUSET_FUSED
#define MADNESS_VER
#define SIZE_FUNC	8

using namespace madness;

double rtclock();

static const double L = 20;     // Half box size
static const long k = 8;        // wavelet order
//static const double thresh = 1e-6; // precision
static const double thresh = 1e-12; // precision
static const double c = 2.0;       //
static const double tstep = 0.1;
static const double alpha = 1.9; // Exponent
static const double VVV = 0.2;  // Vp constant value

// Initial Gaussian with exponent alpha
static double uinitial(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return exp(-alpha*(2*x*x+3.2*y*y+1.7*z*z))*pow(constants::pi/alpha,-1.5);
}

static double uinitial2(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return exp(-alpha*(5*x*x+y*y+z*z))*pow(constants::pi/alpha,-1.5);
}

static double Vp(const coord_3d& r) {
    return VVV;
}

// Exact solution at time t
class uexact : public FunctionFunctorInterface<double,3> {
    double t;
public:
    uexact(double t) : t(t) {}

    double operator()(const coord_3d& r) const {
        const double x=r[0], y=r[1], z=r[2];
        double rsq = (x*x+y*y+z*z);

        return exp(VVV*t)*exp(-rsq*alpha/(1.0+4.0*alpha*t*c)) * pow(alpha/((1+4*alpha*t*c)*constants::pi),1.5);
    }
};

// Functor to compute exp(f) where f is a madness function
template<typename T, int NDIM>
struct unaryexp {
    void operator()(const Key<NDIM>& key, Tensor<T>& t) const {
        UNARY_OPTIMIZED_ITERATOR(T, t, *_p0 = exp(*_p0););
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};

int main(int argc, char** argv) 
{
	int i, j;
	double clkbegin, clkend;

    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);

    startup(world, argc, argv);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_cubic_cell(-L, L);
  //  FunctionDefaults<3>::set_max_refine_level(4);

	if (world.rank() == 0) print("================================================");
	if (world.rank() == 0) print(" Initializing Functions");
	if (world.rank() == 0) print("================================================");

	real_function_3d f_0 = real_factory_3d(world).f(uinitial);
	real_function_3d f_1 = real_factory_3d(world).f(uinitial2);
	real_function_3d f_2 = real_factory_3d(world).f(uinitial);
	real_function_3d g_0 = real_factory_3d(world).f(uinitial2);
	real_function_3d g_1 = real_factory_3d(world).f(uinitial);
	real_function_3d g_2 = real_factory_3d(world).f(uinitial2);

	f_0._treeName = "f_0";
	f_1._treeName = "f_1";
	f_2._treeName = "f_2";
	g_0._treeName = "g_0";
	g_1._treeName = "g_1";
	g_2._treeName = "g_2";

	f_0.truncate();
	f_1.truncate();
	f_2.truncate();
	g_0.truncate();
	g_1.truncate();
	g_2.truncate();

//
//
//
//
//
	if (world.rank() == 0) print("================================================");
	if (world.rank() == 0) print("By MADNESS");
	if (world.rank() == 0) print("================================================");

	if (world.rank() == 0) print("== [Running Operations] ================================");
	world.gop.fence();
	clkbegin = rtclock();
	f_0.compress();
	f_1.compress();
	f_2.compress();
	g_0.compress();
	g_1.compress();
	g_2.compress();

	double MAD_results[9];
	
	MAD_results[0] = f_0.inner(g_0);
	MAD_results[1] = f_0.inner(g_1);
	MAD_results[2] = f_0.inner(g_2);
	MAD_results[3] = f_1.inner(g_0);
	MAD_results[4] = f_1.inner(g_1);
	MAD_results[5] = f_1.inner(g_2);
	MAD_results[6] = f_2.inner(g_0);
	MAD_results[7] = f_2.inner(g_1);
	MAD_results[8] = f_2.inner(g_2);

	clkend = rtclock() - clkbegin;
	if (world.rank() == 0) printf("Running Time: %f\n", clkend);
	if (world.rank() == 0)
	{
		for (i=0; i<9; i++)
			printf ("Inner-Product [%d] = %f\n", i, MAD_results[i]);
	}

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

