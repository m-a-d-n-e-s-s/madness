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

//#define FUSET_UNFUSED
#define FUSET_FUSED
#define MADNESS_VER

using namespace madness;

#define FUNC_SIZE_N	3
#define FUNC_SIZE_M	3
double rtclock();

static const double L = 20;     // Half box size
static const long k = 8;        // wavelet order
static const double thresh = 1e-6; // precision
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
    //FunctionDefaults<3>::set_max_refine_level(4);

    if (world.rank() == 0) print("Creating (N + M) functions");
	real_function_3d u[FUNC_SIZE_N];
	real_function_3d g[FUNC_SIZE_M];
	
	for (i=0; i<FUNC_SIZE_N; i++)
	{ 
		u[i] = real_factory_3d(world).f(uinitial); 
		u[i].truncate();
	}
	for (i=0; i<FUNC_SIZE_M; i++)
	{
		g[i] = real_factory_3d(world).f(uinitial);
		g[i].truncate();
	}

	real_function_3d  u_result_factory[FUNC_SIZE_N];
	real_function_3d  g_result_factory[FUNC_SIZE_M];
	real_function_3d  u_result_factory_fused[FUNC_SIZE_N];
	real_function_3d  g_result_factory_fused[FUNC_SIZE_M];
	real_function_3d* u_result[FUNC_SIZE_N];
	real_function_3d* g_result[FUNC_SIZE_M];
	real_function_3d* u_result_fused[FUNC_SIZE_N];
	real_function_3d* g_result_fused[FUNC_SIZE_M];

	for (i=0; i<FUNC_SIZE_N; i++) {
		u_result_factory[i]			= real_factory_3d(world);
		u_result_factory_fused[i]	= real_factory_3d(world);
		u_result[i]					= new real_function_3d(u_result_factory[i]);
		u_result_fused[i]			= new real_function_3d(u_result_factory_fused[i]);
	}

	for (i=0; i<FUNC_SIZE_M; i++) {
		g_result_factory[i]			= real_factory_3d(world);
		g_result_factory_fused[i]	= real_factory_3d(world);
		g_result[i]					= new real_function_3d(g_result_factory[i]);
		g_result_fused[i]			= new real_function_3d(g_result_factory_fused[i]);
	}

	real_function_3d  ug_result_factory[FUNC_SIZE_N][FUNC_SIZE_M];
	real_function_3d  ug_result_factory_fused[FUNC_SIZE_N][FUNC_SIZE_M];
	real_function_3d* ug_result[FUNC_SIZE_N][FUNC_SIZE_M];
	real_function_3d* ug_result_fused[FUNC_SIZE_N][FUNC_SIZE_M];

	for (i=0; i<FUNC_SIZE_N; i++)
		for (j=0; j<FUNC_SIZE_M; j++)
		{
			ug_result_factory[i][j]			= real_factory_3d(world);
			ug_result_factory_fused[i][j]	= real_factory_3d(world);
			ug_result[i][j]					= new real_function_3d(ug_result_factory[i][j]);
			ug_result_fused[i][j]			= new real_function_3d(ug_result_factory_fused[i][j]);
		}

#ifdef FUSET_UNFUSED
	if (world.rank() == 0) print("================================================");
	if (world.rank() == 0) print("By FuseT, UnFUSED");
	if (world.rank() == 0) print("================================================");
	world.gop.fence();

	OpExecutor<double,3>	exe(world);
	CompressOp<double,3>*	op_u[FUNC_SIZE_N];
	CompressOp<double,3>*	op_g[FUNC_SIZE_M];
	InnerOp<double,3>*		op_ug[FUNC_SIZE_N][FUNC_SIZE_M];

	clkbegin = rtclock();

	for (i=0; i<FUNC_SIZE_N; i++) op_u[i] = new CompressOp<double,3>("CompressOp",u_result[i], &u[i]);
	for (i=0; i<FUNC_SIZE_M; i++) op_g[i] = new CompressOp<double,3>("CompressOp",g_result[i], &g[i]);

	for (i=0; i<FUNC_SIZE_N; i++) exe.execute(op_u[i], true);
	for (i=0; i<FUNC_SIZE_M; i++) exe.execute(op_g[i], true);

	for (i=0; i<FUNC_SIZE_N; i++)
		for (j=0; j<FUNC_SIZE_M; j++)
			op_ug[i][j] = new InnerOp<double,3>("InnerOp",ug_result[i][j],u_result[i],g_result[j]);	

	for (i=0; i<FUNC_SIZE_N; i++)
		for (j=0; j<FUNC_SIZE_M; j++)
			exe.execute(op_ug[i][j], true);

	clkend = rtclock() - clkbegin;
	if (world.rank() == 0) print ("================================================");
	world.gop.fence();
	for (i=0; i<FUNC_SIZE_N; i++)
		u[i].verify_tree();
	for (i=0; i<FUNC_SIZE_M; i++)
		g[i].verify_tree();
	if (world.rank() == 0) print ("================================================");
	world.gop.fence();

	if (world.rank() == 0) print ("================================================");
	if (world.rank() == 0) printf("Running Time: %f\n", clkend);
	if (world.rank() == 0) print ("================================================");
	world.gop.fence();

/*
	for (i=0; i<FUNC_SIZE_N; i++)
		for (j=0; j<FUNC_SIZE_M; j++)
			if (world.rank() == 0) print("Inner-Product Operators [",i,"][",j,"]", op_ug[i][j]->_sum);
	world.gop.fence();
*/
#endif 


#ifdef FUSET_FUSED
	if (world.rank() == 0) print("================================================");
	if (world.rank() == 0) print("By FuseT, FUSED");
	if (world.rank() == 0) print("================================================");

	CompressOp<double,3>*			op_u_fused[FUNC_SIZE_N];
	CompressOp<double,3>*			op_g_fused[FUNC_SIZE_M];
	InnerOp<double,3>*				op_ug_fused[FUNC_SIZE_N][FUNC_SIZE_M];
	vector<PrimitiveOp<double,3>*>	sequence;
	FuseT<double,3>*				odag;

	for (i=0; i<FUNC_SIZE_N; i++)	op_u_fused[i] = new CompressOp<double,3>("CompressOp-u",u_result_fused[i], &u[i]);
	for (i=0; i<FUNC_SIZE_M; i++)	op_g_fused[i] = new CompressOp<double,3>("CompressOp-g",g_result_fused[i], &g[i]);

	for (i=0; i<FUNC_SIZE_N; i++)	sequence.push_back(op_u_fused[i]);
	for (i=0; i<FUNC_SIZE_M; i++)	sequence.push_back(op_g_fused[i]);

	for (i=0; i<FUNC_SIZE_N; i++)
		for (j=0; j<FUNC_SIZE_M; j++)
			op_ug_fused[i][j] = new InnerOp<double,3>("InnerOp",ug_result_fused[i][j],u_result_fused[i],g_result_fused[j]);	

	for (i=0; i<FUNC_SIZE_N; i++)
		for (j=0; j<FUNC_SIZE_M; j++)
			sequence.push_back(op_ug_fused[i][j]);

	clkbegin = rtclock();

	odag = new FuseT<double,3>(sequence);
	odag->processSequence();
	
	FusedOpSequence<double,3> fsequence = odag->getFusedOpSequence();
	FusedExecutor<double,3> fexecutor(world, &fsequence);
	fexecutor.execute();
    world.gop.fence();

	clkend = rtclock() - clkbegin;
	if (world.rank() == 0) print ("================================================");
	if (world.rank() == 0) printf("Running Time: %f\n", clkend);
	if (world.rank() == 0) print ("================================================");
	world.gop.fence();
/*
	if (world.rank() == 0) print("After Inner-Product Operators - u*g");
	for (i=0; i<FUNC_SIZE_N; i++)
		for (j=0; j<FUNC_SIZE_M; j++)
			if (world.rank() == 0) print("Inner-Product Operators [",i,"][",j,"]", op_ug_fused[i][j]->_sum);
	world.gop.fence();
*/
#endif



#ifdef MADNESS_VER
	if (world.rank() == 0) print("================================================");
	if (world.rank() == 0) print("By MADNESS");
	if (world.rank() == 0) print("================================================");

	double omg[FUNC_SIZE_N][FUNC_SIZE_M];

	clkbegin = rtclock();

	for (i=0; i<FUNC_SIZE_N; i++) u[i].compress();
	for (i=0; i<FUNC_SIZE_M; i++) g[i].compress();

	for (i=0; i<FUNC_SIZE_N; i++)
		for (j=0; j<FUNC_SIZE_M; j++)
			omg[i][j] = u[i].inner(g[j]);

	clkend = rtclock() - clkbegin;
	if (world.rank() == 0) print ("================================================");
	if (world.rank() == 0) printf("Running Time: %f\n", clkend);
	if (world.rank() == 0) print ("================================================");
	world.gop.fence();

/*
	for (i=0; i<FUNC_SIZE_N; i++)
		for (j=0; j<FUNC_SIZE_M; j++)		
			if (world.rank() == 0) print("Inner-Product Operators [",i,"][",j,"]", omg[i][j]);
	world.gop.fence();
*/
#endif

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

