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
//#define FUSET_FUSED
//#define MADNESS_VER

/*!
\file heat2.cc
\brief Example Green function for the 3D heat equation with a linear term
\defgroup heatex2 Evolve in time 3D heat equation with a linear term
\ingroup examples

The source is <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local/trunk/src/apps/examples/heat2.cc>here</a>.

\par Points of interest
  - application of a function of a function to exponentiate the potential
  - use of a functor to compute the solution at an arbitrary future time
  - convolution with the Green's function


\par Background

This adds to the complexity of the other \ref exampleheat "heat equation example"
by including a linear term.  Specifically, we solve
\f[
  \frac{\partial u(x,t)}{\partial t} = c \nabla^2 u(x,t) + V_p(x,t) u(x,t)
\f]
If \f$ V_p = 0 \f$ time evolution operator is
\f[
  G_0(x,t) = \frac{1}{\sqrt{4 \pi c t}} \exp \frac{-x^2}{4 c t}
\f]
For non-zero \f$ V_p \f$ the time evolution is performed using the Trotter splitting
\f[
  G(x,t) = G_0(x,t/2) * \exp(V_p t) * G_0(x,t/2) + O(t^3)
\f]
In order to form an exact solution for testing, we choose \f$ V_p(x,t)=\mbox{constant} \f$
but the solution method is not limited to this choice.

*/

using namespace madness;

#define FUNC_SIZE	10
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
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);

    startup(world, argc, argv);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_cubic_cell(-L, L);
    FunctionDefaults<3>::set_max_refine_level(4);

	// Step: 1
	// for i = 0  to N-1
	//	u[i].create();
    if (world.rank() == 0) print("Creating functions");
	real_function_3d u[FUNC_SIZE];
	real_function_3d u_result_factory[FUNC_SIZE];
	real_function_3d* u_result[FUNC_SIZE];
	real_function_3d* u_result_fused[FUNC_SIZE];
	real_function_3d g[FUNC_SIZE];
	real_function_3d g_result_factory[FUNC_SIZE];
	real_function_3d* g_result[FUNC_SIZE];
	real_function_3d* g_result_fused[FUNC_SIZE];

	//real_function_3d temp_factory = real_factory_3d(world);
	//real_function_3d temp(temp_factory);
	real_function_3d temp_factory[FUNC_SIZE][FUNC_SIZE];
	real_function_3d* temp[FUNC_SIZE][FUNC_SIZE];

	for (i=0; i<FUNC_SIZE; i++)
	{
		u[i] = real_factory_3d(world).f(uinitial);
		u_result_factory[i] = real_factory_3d(world);


		g[i] = real_factory_3d(world).f(uinitial2);
		g_result_factory[i] = real_factory_3d(world);
	}


	for (i=0; i<FUNC_SIZE; i++)
	{
		for (j=0; j<FUNC_SIZE; j++)
		{
			temp_factory[i][j] = real_factory_3d(world);
			temp[i][j] = new real_function_3d(temp_factory[i][j]);
		}
	}


    //real_function_3d result_factory = real_factory_3d(world);
    //real_function_3d result(result_factory);
	for (i=0; i<FUNC_SIZE; i++)
	{
		// input
		u[i]._treeName = "u1";
		u[i].truncate();
		
		// result	
		u_result[i] = new real_function_3d(u_result_factory[i]);
		u_result_fused[i] = new real_function_3d(u_result_factory[i]);
	}

    if (world.rank() == 0) print("Creating g[0,..,N-1]");
	for (i=0; i<FUNC_SIZE; i++)
	{
		// input
		g[i]._treeName = "u1";
		g[i].truncate();

		// result
		g_result[i] = new real_function_3d(g_result_factory[i]);
		g_result_fused[i] = new real_function_3d(g_result_factory[i]);
	}

	double clkbegin, clkend;

#ifdef FUSET_UNFUSED
	if (world.rank() == 0) print("================================================");
	if (world.rank() == 0) print("By FuseT, UnFUSED");
	if (world.rank() == 0) print("================================================");
		
	clkbegin = rtclock();

	CompressOp<double,3>* op_u[FUNC_SIZE];
	CompressOp<double,3>* op_g[FUNC_SIZE];

	if (world.rank() == 0) print("Creating Compress Operators - u & g");
	for (i=0; i<FUNC_SIZE; i++)
	{
		op_u[i] = new CompressOp<double,3>("CompressOp-u",u_result[i], &u[i]);
		op_g[i] = new CompressOp<double,3>("CompressOp-g",g_result[i], &g[i]);
	}

	OpExecutor<double,3> exe(world);
	if (world.rank() == 0) print("Running Compress Operators - u");
	for (i=0; i<FUNC_SIZE; i++)
	{
		exe.execute(op_u[i], true);
	}

	if (world.rank() == 0) print("Running Compress Operators - g");
	for (i=0; i<FUNC_SIZE; i++)
	{
		exe.execute(op_g[i], true);
	}

	if (world.rank() == 0) print("Creating Inner-Product Operators - u*g");
	InnerOp<double,3>* op_ug[FUNC_SIZE][FUNC_SIZE];
	for (i=0; i<FUNC_SIZE; i++)
	{
		for (j=0; j<FUNC_SIZE; j++)
		{
			//op_ug[i][j] = new InnerOp<double,3>("InnerOp-u*g",&temp,u_result[i],g_result[j]);	
			op_ug[i][j] = new InnerOp<double,3>("InnerOp-u*g",temp[i][j],u_result[i],g_result[j]);	
		}
	}

	if (world.rank() == 0) print("Running Inner-Product Operators - u*g");
	for (i=0; i<FUNC_SIZE; i++)
	{
		for (j=0; j<FUNC_SIZE; j++)
		{
			exe.execute(op_ug[i][j], true);
		}
	}

	clkend = rtclock() - clkbegin;

	if (world.rank() == 0) printf("Running Time: %f\n", clkend);
	if (world.rank() == 0) print("After Inner-Product Operators - u*g");

	for (i=0; i<FUNC_SIZE; i++)
		for (j=0; j<FUNC_SIZE; j++)
			if (world.rank() == 0) print("Inner-Product Operators [",i,"][",j,"]", op_ug[i][j]->_sum);
	// Correctness by using MADNESS
#endif 


#ifdef FUSET_FUSED
	if (world.rank() == 0) print("================================================");
	if (world.rank() == 0) print("By FuseT, FUSED");
	if (world.rank() == 0) print("================================================");

	clkbegin = rtclock();
	CompressOp<double,3>* op_u_fused[FUNC_SIZE];
	CompressOp<double,3>* op_g_fused[FUNC_SIZE];
	for (i=0; i<FUNC_SIZE; i++)
	{
		op_u_fused[i] = new CompressOp<double,3>("CompressOp-u",u_result_fused[i], &u[i]);
		op_g_fused[i] = new CompressOp<double,3>("CompressOp-g",g_result_fused[i], &g[i]);
	}

	real_function_3d temp_factory_fused = real_factory_3d(world);
	real_function_3d temp_fused(temp_factory_fused);
	InnerOp<double,3>* op_ug_fused[FUNC_SIZE][FUNC_SIZE];
	for (i=0; i<FUNC_SIZE; i++)
	{
		for (j=0; j<FUNC_SIZE; j++)
		{
			op_ug_fused[i][j] = new InnerOp<double,3>("InnerOp-u*g",&temp_fused,u_result_fused[i],g_result_fused[j]);	
		}
	}

	vector<PrimitiveOp<double,3>*> sequence;
	for (i=0; i<FUNC_SIZE; i++)
	{
		sequence.push_back(op_u_fused[i]);
		sequence.push_back(op_g_fused[i]);
	}
	for (i=0; i<FUNC_SIZE; i++)
	{
		for (j=0; j<FUNC_SIZE; j++)
		{
			sequence.push_back(op_ug_fused[i][j]);
		}
	}

	FuseT<double,3> odag(sequence);
	odag.processSequence();
	
	FusedOpSequence<double,3> fsequence = odag.getFusedOpSequence();
	FusedExecutor<double,3> fexecutor(world, &fsequence);
	fexecutor.execute();
    world.gop.fence();

	clkend = rtclock() - clkbegin;
	if (world.rank() == 0) printf("Running Time: %f\n", clkend);

	if (world.rank() == 0) print("After Inner-Product Operators - u*g");
	for (i=0; i<FUNC_SIZE; i++)
	{
		for (j=0; j<FUNC_SIZE; j++)
		{
			if (world.rank() == 0) print("Inner-Product Operators [",i,"][",j,"]", op_ug_fused[i][j]->_sum);
		}
	}
	// Correctness by using MADNESS
#endif
/*
	for (i=0; i<FUNC_SIZE; i++)
	{
		for (j=0; j<FUNC_SIZE; j++)
		{
			if (world.rank() == 0) printf("[%d][%d] - [%d][%d]: %f ", i, j, i, j, op_ug_fused[i][j]->_sum - op_ug[i][j]->_sum);
		}
		if (world.rank() == 0) printf("\n");
	}
*/
#ifdef MADNESS_VER
	if (world.rank() == 0) print("================================================");
	if (world.rank() == 0) print("By MADNESS");
	if (world.rank() == 0) print("================================================");

	clkbegin = rtclock();
	for (i=0; i<FUNC_SIZE; i++)
	{
		u[i].compress();
		g[i].compress();
	}

	double omg[FUNC_SIZE][FUNC_SIZE];
	for (i=0; i<FUNC_SIZE; i++)
	{
		for (j=0; j<FUNC_SIZE; j++)
		{
			omg[i][j] = u[i].inner(g[j]);
		}
	}

	clkend = rtclock() - clkbegin;
	if (world.rank() == 0) printf("Running Time: %f\n", clkend);

	for (i=0; i<FUNC_SIZE; i++)
		for (j=0; j<FUNC_SIZE; j++)		
			if (world.rank() == 0) print("Inner-Product Operators [",i,"][",j,"]", omg[i][j]);
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

