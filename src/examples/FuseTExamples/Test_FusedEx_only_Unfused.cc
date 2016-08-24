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

	real_function_3d rf_0_factory = real_factory_3d(world);
	real_function_3d rf_0(rf_0_factory);
	real_function_3d rf_1_factory = real_factory_3d(world);
	real_function_3d rf_1(rf_1_factory);
	real_function_3d rf_2_factory = real_factory_3d(world);
	real_function_3d rf_2(rf_2_factory);
	real_function_3d rg_0_factory = real_factory_3d(world);
	real_function_3d rg_0(rg_0_factory);
	real_function_3d rg_1_factory = real_factory_3d(world);
	real_function_3d rg_1(rg_1_factory);
	real_function_3d rg_2_factory = real_factory_3d(world);
	real_function_3d rg_2(rg_2_factory);
	
	real_function_3d rfg_0_factory = real_factory_3d(world);
	real_function_3d rfg_0(rfg_0_factory);
	real_function_3d rfg_1_factory = real_factory_3d(world);
	real_function_3d rfg_1(rfg_1_factory);
	real_function_3d rfg_2_factory = real_factory_3d(world);
	real_function_3d rfg_2(rfg_2_factory);
	real_function_3d rfg_3_factory = real_factory_3d(world);
	real_function_3d rfg_3(rfg_3_factory);
	real_function_3d rfg_4_factory = real_factory_3d(world);
	real_function_3d rfg_4(rfg_4_factory);
	real_function_3d rfg_5_factory = real_factory_3d(world);
	real_function_3d rfg_5(rfg_5_factory);
	real_function_3d rfg_6_factory = real_factory_3d(world);
	real_function_3d rfg_6(rfg_6_factory);
	real_function_3d rfg_7_factory = real_factory_3d(world);
	real_function_3d rfg_7(rfg_7_factory);
	real_function_3d rfg_8_factory = real_factory_3d(world);
	real_function_3d rfg_8(rfg_8_factory);

	// Creating Operators
	CompressOp<double,3> op_compress_f_0("Compress", &rf_0, &f_0);
	CompressOp<double,3> op_compress_f_1("Compress", &rf_1, &f_1);
	CompressOp<double,3> op_compress_f_2("Compress", &rf_2, &f_2);
	CompressOp<double,3> op_compress_g_0("Compress", &rg_0, &g_0);
	CompressOp<double,3> op_compress_g_1("Compress", &rg_1, &g_1);
	CompressOp<double,3> op_compress_g_2("Compress", &rg_2, &g_2);

	InnerOp<double,3> op_inner_f_g_0_0("Inner", &rfg_0, &rf_0, &rg_0);
	InnerOp<double,3> op_inner_f_g_0_1("Inner", &rfg_1, &rf_0, &rg_1);
	InnerOp<double,3> op_inner_f_g_0_2("Inner", &rfg_2, &rf_0, &rg_2);
	InnerOp<double,3> op_inner_f_g_1_0("Inner", &rfg_3, &rf_1, &rg_0);
	InnerOp<double,3> op_inner_f_g_1_1("Inner", &rfg_4, &rf_1, &rg_1);
	InnerOp<double,3> op_inner_f_g_1_2("Inner", &rfg_5, &rf_1, &rg_2);
	InnerOp<double,3> op_inner_f_g_2_0("Inner", &rfg_6, &rf_2, &rg_0);
	InnerOp<double,3> op_inner_f_g_2_1("Inner", &rfg_7, &rf_2, &rg_1);
	InnerOp<double,3> op_inner_f_g_2_2("Inner", &rfg_8, &rf_2, &rg_2);

	// OpExecutor!!!!!!
	OpExecutor<double,3> exe(world);
	if (world.rank() == 0) print("================================================");
	if (world.rank() == 0) print("== By FuseT, UNFUSED");
	if (world.rank() == 0) print("================================================");
	clkbegin = rtclock();

	// Compress
	if (world.rank() == 0) print("== [Running Compress Operations] ========================");
	world.gop.fence();

	exe.execute(&op_compress_f_0,true);
	exe.execute(&op_compress_f_1,true);
	exe.execute(&op_compress_f_2,true);
	exe.execute(&op_compress_g_0,true);
	exe.execute(&op_compress_g_1,true);
	exe.execute(&op_compress_g_2,true);

	// Inner
	if (world.rank() == 0) print("== [Running Inner Operations] ========================");
	world.gop.fence();

	exe.execute(&op_inner_f_g_0_0,true);
	exe.execute(&op_inner_f_g_0_1,true);
	exe.execute(&op_inner_f_g_0_2,true);
	exe.execute(&op_inner_f_g_1_0,true);
	exe.execute(&op_inner_f_g_1_1,true);
	exe.execute(&op_inner_f_g_1_2,true);
	exe.execute(&op_inner_f_g_2_0,true);
	exe.execute(&op_inner_f_g_2_1,true);
	exe.execute(&op_inner_f_g_2_2,true);

	//
	//
	//
	clkend = rtclock() - clkbegin;
	if (world.rank() == 0) printf("Running Time: %f\n", clkend);

	if (world.rank() == 0)
	{
		printf ("Inner-Product [%d] = %f\n", 0, op_inner_f_g_0_0._sum);
		printf ("Inner-Product [%d] = %f\n", 1, op_inner_f_g_0_1._sum);
		printf ("Inner-Product [%d] = %f\n", 2, op_inner_f_g_0_2._sum);
		printf ("Inner-Product [%d] = %f\n", 3, op_inner_f_g_1_0._sum);
		printf ("Inner-Product [%d] = %f\n", 4, op_inner_f_g_1_1._sum);
		printf ("Inner-Product [%d] = %f\n", 5, op_inner_f_g_1_2._sum);
		printf ("Inner-Product [%d] = %f\n", 6, op_inner_f_g_2_0._sum);
		printf ("Inner-Product [%d] = %f\n", 7, op_inner_f_g_2_1._sum);
		printf ("Inner-Product [%d] = %f\n", 8, op_inner_f_g_2_2._sum);
	}
		
//
//
//
//
//

//
//
//
//
//
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

