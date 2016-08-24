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
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/constants.h>
#include <madness/mra/FuseT/InnerOp.h>
#include <madness/mra/FuseT/CompressOp.h>
#include <madness/mra/FuseT/OpExecutor.h>
#include <madness/mra/FuseT/FusedExecutor.h>
#include <madness/mra/FuseT/FuseT.h>
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

static const double L = 20;     // Half box size
static const long k = 8;        // wavelet order
static const double thresh = 1e-6; // precision
static const double c = 2.0;       //
static const double tstep = 0.1;
static const double alpha = 1.9; // Exponent
static const double VVV = 0.2;  // Vp constant value

#define FUNC_SIZE	2
#define FUNC_SIZE_M	2

double rtclock();

// Initial Gaussian with exponent alpha
static double uinitial(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return exp(-alpha*(2*x*3.2*x+y*y+1.7*z*z))*pow(constants::pi/alpha,-1.5);
}

static double guess(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return 6.0*exp(-2.0*sqrt(x*x+y*y+z*z+1e-4));
}

static double ghaly(const coord_3d& r) {
	std::srand(time(NULL));
	const double randVal = std::rand()/1000000000.0;
    const double x=r[0], y=r[1], z=r[2];
    return 3.0*exp(-2.0*sqrt(x*x + randVal*randVal + y*y + z*z + 1e-4));
}
static double uinitial1(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return exp(-alpha*(2*x*x+1.4*y*y+z*z))*pow(constants::pi/alpha,-1.5);
};

static double Vp(const coord_3d& r) {
    return VVV;
}

class alpha_functor : public FunctionFunctorInterface<double,3> {
private:
    double coeff;
public:
    alpha_functor(double coeff=1.0) : coeff(coeff) {}

    virtual double operator()(const coord_3d& r) const {
        const double x=r[0], y=r[1], z=r[2];
        return (coeff * (x*x + y*y + z*z) * sin(x*x + y*y + z*z));
    }
};
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
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);

    startup(world, argc, argv);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_cubic_cell(-L, L);
	FunctionDefaults<3>::set_max_refine_level(14);

	if (world.rank() == 0) print ("====================================================");
    if (world.rank() == 0) printf("   Initializing Functions\n");
    if (world.rank() == 0) printf("     %d Functions, %d Functions\n", FUNC_SIZE, FUNC_SIZE_M);
	if (world.rank() == 0) print ("====================================================");
    world.gop.fence();

	// 2 * N Functions
	real_function_3d  h[FUNC_SIZE];
	real_function_3d  g[FUNC_SIZE_M];

	// N * N Results Functions by Inner-Product
	real_function_3d  temp_factory_h[FUNC_SIZE];
	real_function_3d  temp_factory_g[FUNC_SIZE_M];
	real_function_3d* temp_h[FUNC_SIZE];
	real_function_3d* temp_g[FUNC_SIZE_M];

	real_function_3d  temp_factory_hg[FUNC_SIZE][FUNC_SIZE_M];
	real_function_3d* temp_hg[FUNC_SIZE][FUNC_SIZE_M];

	int i, j;
	double clkbegin, clkend;
	clkbegin = rtclock();
	for (i=0; i<FUNC_SIZE; i++) {
		//if (world.rank() == 0) print (" Creating h...", i);
		h[i]	= real_factory_3d(world).f(guess);
	}

	for (i=0; i<FUNC_SIZE_M; i++) {
		//if (world.rank() == 0) print (" Creating g...", i);
		g[i]	= real_factory_3d(world).f(uinitial);
	}

	//FunctionDefaults<3>::set_max_refine_level(30);
	for (i=0; i<FUNC_SIZE; i++) 
	{
		temp_factory_h[i]	= real_factory_3d(world);
		temp_h[i]			= new real_function_3d(temp_factory_h[i]);
	}
		
	for (j=0; j<FUNC_SIZE_M; j++)
	{
		temp_factory_g[j]	= real_factory_3d(world);
		temp_g[j]			= new real_function_3d(temp_factory_g[j]);
	}
	
	for (i=0; i<FUNC_SIZE; i++)
		for (j=0; j<FUNC_SIZE_M; j++)
		{
			temp_factory_hg[i][j]	= real_factory_3d(world);
			temp_hg[i][j]			= new real_function_3d(temp_factory_hg[i][j]);
		}


	for (i=0; i<FUNC_SIZE; i++)
	{
		h[i].truncate();
	//	h[i].compress();
	}
	for (i=0; i<FUNC_SIZE_M; i++)
	{
		g[i].truncate();
	//	g[i].compress();
	}

	double result_h_norm;
	double result_g_norm;
	double result_h_trace;
	double result_g_trace;

	clkend = rtclock() - clkbegin;
	if (world.rank() == 0) printf("Running Time: %f\n", clkend);
	world.gop.fence();


	if (world.rank() == 0) print ("====================================================");
	if (world.rank() == 0) print ("==      FUSET-UN-FUSED       ==========================");
	if (world.rank() == 0) print ("====================================================");
	world.gop.fence();

	CompressOp<double,3>* inner_op_h[FUNC_SIZE];
	CompressOp<double,3>* inner_op_g[FUNC_SIZE_M];
	clkbegin = rtclock();
	for (i=0; i<FUNC_SIZE; i++)
		inner_op_h[i] = new CompressOp<double,3>("Inner",temp_h[i],&h[i]);
		
	for (j=0; j<FUNC_SIZE_M; j++)
		inner_op_g[j] = new CompressOp<double,3>("Inner",temp_g[j],&g[j]);

	//
	if (world.rank() == 0) print ("======= OpExecutor =========================================");
	world.gop.fence();
	OpExecutor<double,3> exe(world);
	for (i=0; i<FUNC_SIZE; i++)
		exe.execute(inner_op_h[i], true);

	for (j=0; j<FUNC_SIZE_M; j++)
		exe.execute(inner_op_g[j], true);


//	InnerOp<double,3>* inner_op_hg[FUNC_SIZE][FUNC_SIZE_M];
/*	
	if (world.rank() == 0) print ("======= OpExecutor =========================================");
	for (i=0; i<FUNC_SIZE; i++)
		for (j=0; j<FUNC_SIZE_M; j++)
			inner_op_hg[i][j] = new InnerOp<double,3>("Inner", temp_hg[i][j], temp_h[i], temp_g[j]);

	for (i=0; i<FUNC_SIZE; i++)
		for (j=0; j<FUNC_SIZE_M; j++)
			exe.execute(inner_op_hg[i][j],true);
*/
/*

	vector<PrimitiveOp<double,3>*> sequence;
	for (i=0; i<FUNC_SIZE; i++)
		sequence.push_back(inner_op_h[i]);

	for (j=0; j<FUNC_SIZE_M; j++)
		sequence.push_back(inner_op_g[j]);

	FuseT<double,3> odag(sequence);
	odag.processSequence();

	FusedOpSequence<double,3> fsequence = odag.getFusedOpSequence();
	FusedExecutor<double,3> fexecuter(world, &fsequence);
	fexecuter.execute();
*/
	clkend = rtclock() - clkbegin;
	if (world.rank() == 0) printf("Running Time: %f\n", clkend);
	world.gop.fence();

	for (i=0; i<FUNC_SIZE; i++)
	{
		result_h_norm = temp_h[i]->norm2();
		result_h_trace = temp_h[i]->trace();
		if (world.rank() == 0) print (i,"norm:", result_h_norm, " trace", result_h_trace);
	}

	for (j=0; j<FUNC_SIZE_M; j++)	
	{
		result_g_norm = temp_g[j]->norm2();
		result_g_trace = temp_g[j]->trace();
		if (world.rank() == 0) print (i,"norm:", result_g_norm, " trace", result_g_trace);
	}

	for (i=0; i<FUNC_SIZE; i++)
		for (j=0; j<FUNC_SIZE_M; j++)
	//		if (world.rank() == 0) printf ("%d:%d = %f\n", i, j, inner_op_hg[i][j]->_sum);
	world.gop.fence();
//
//
//
	if (world.rank() == 0) print ("====================================================");
	if (world.rank() == 0) print ("==      MADNESS       ==============================");
	if (world.rank() == 0) print ("====================================================");
	world.gop.fence();


	clkbegin = rtclock();
	
	for (i=0; i<FUNC_SIZE; i++)
		h[i].compress();		

	for (j=0; j<FUNC_SIZE_M; j++)
		g[j].compress();

	clkend = rtclock() - clkbegin;
	if (world.rank() == 0) printf("Running Time: %f\n", clkend);
	world.gop.fence();

	for (i=0; i<FUNC_SIZE; i++)
	{
		result_h_norm = h[i].norm2();
		result_h_trace = h[i].trace();
		if (world.rank() == 0) print (i,"norm:", result_h_norm, " trace", result_h_trace);
	}

	for (j=0; j<FUNC_SIZE_M; j++)	
	{
		result_g_norm = g[j].norm2();
		result_g_trace = g[j].trace();
		if (world.rank() == 0) print (i,"norm:", result_g_norm, " trace", result_g_trace);
	}

	double hello[FUNC_SIZE][FUNC_SIZE_M];
	for (i=0; i<FUNC_SIZE; i++)
		for (j=0; j<FUNC_SIZE_M; j++)
			hello[i][j] = h[i].inner(g[j]);

	for (i=0; i<FUNC_SIZE; i++)
		for (j=0; j<FUNC_SIZE_M; j++)
			if (world.rank() == 0) printf ("%d:%d = %f\n", i, j, hello[i][j]);


	world.gop.fence();

    finalize();    
    return 0;
}

double rtclock()
{struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0)
	printf("Error return from gettimeofday: %d", stat);
    return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

