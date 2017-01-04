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

using namespace madness;

static const double L = 20;     // Half box size
static const long k = 8;        // wavelet order
static const double thresh = 1e-6; // precision   // w/o diff. and 1e-12 -> 64 x 64
static const double c = 2.0;       //
static const double tstep = 0.1;
static const double alpha = 1.9; // Exponent
static const double VVV = 0.2;  // Vp constant value

#define FUNC_SIZE	64
#define FUNC_SIZE_M	64
#define MIN_NODES	4
#define SCALE		MIN_NODES/4

double rtclock();

// Initial Gaussian with exponent alpha
static double uinitial(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return exp(-alpha*(2*x*3.2*x+y*y+1.7*z*z))*pow(constants::pi/alpha,-1.5);
}

static double uinitial2(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
	srand(time(0));
    return exp(-alpha*(5*x*x+y*y+z*z))*pow(constants::pi/alpha,-1.5);
}
static double uinitial1(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return exp(-alpha*(2*x*x+1.4*y*y+z*z))*pow(constants::pi/alpha,-1.5);
};

static double ghaly(const coord_3d& r) {
	std::srand(time(NULL));
	const double randVal = std::rand()/1000000000.0;
    const double x=r[0], y=r[1], z=r[2];
    return 3.0*exp(-2.0*sqrt(x*x + randVal*randVal + y*y + z*z + 1e-4));
}
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
    if (world.rank() == 0) printf("     Max-refine-level: %d Functions\n", 14);
	if (world.rank() == 0) print ("====================================================");
    world.gop.fence();

	// N and M Functions
	real_function_3d  h[FUNC_SIZE];
	real_function_3d  g[FUNC_SIZE_M];

	// N and M (Result) Functions
	real_function_3d  temp_factory_h[FUNC_SIZE];
	real_function_3d  temp_factory_g[FUNC_SIZE_M];
	real_function_3d* temp_h[FUNC_SIZE];
	real_function_3d* temp_g[FUNC_SIZE_M];

	// (SCALE) * N * M Results Functions by Inner-Product
	real_function_3d  temp_factory[SCALE][FUNC_SIZE][FUNC_SIZE_M];
	real_function_3d* temp[SCALE][FUNC_SIZE][FUNC_SIZE_M];

	int i, j, k;
	double clkbegin, clkend;
	clkbegin = rtclock();

	for (i=0; i<FUNC_SIZE; i++) 
	{
		h[i]				= real_factory_3d(world).f(ghaly);
		//h[i]				= real_factory_3d(world).f(uinitial);
		temp_factory_h[i]	= real_factory_3d(world);
		temp_h[i]			= new real_function_3d(temp_factory_h[i]);
	}

	for (i=0; i<FUNC_SIZE_M; i++)
	{
		g[i]				= real_factory_3d(world).f(ghaly);
		//g[i]				= real_factory_3d(world).f(uinitial2);
		temp_factory_g[i]	= real_factory_3d(world);
		temp_g[i]			= new real_function_3d(temp_factory_g[i]);
	}

	for (k=0; k<SCALE; k++) {
		for (i=0; i<FUNC_SIZE; i++) {
			for (j=0; j<FUNC_SIZE_M; j++) {
				temp_factory[k][i][j]	= real_factory_3d(world);
				temp[k][i][j]			= new real_function_3d(temp_factory[k][i][j]);
			}
		}
	}

	for (i=0; i<FUNC_SIZE; i++)		h[i].truncate();
	for (i=0; i<FUNC_SIZE_M; i++)	g[i].truncate();

	clkend = rtclock() - clkbegin;
	if (world.rank() == 0) printf("[Initializing Functions] Running Time: %f\n", clkend);

	if (world.rank() == 0) print ("====================================================");
	if (world.rank() == 0) print ("==      FUSET-FUSED       ==========================");
	if (world.rank() == 0) print ("====================================================");
	world.gop.fence();

	clkbegin = rtclock();

	CompressOp<double,3>*			c_op_h[FUNC_SIZE];
	CompressOp<double,3>*			c_op_g[FUNC_SIZE_M];
	InnerOp<double,3>*				inner_op_ug[SCALE][FUNC_SIZE][FUNC_SIZE_M];
	vector<PrimitiveOp<double,3>*>	sequence;

	for (i=0; i<FUNC_SIZE; i++) 	c_op_h[i] = new CompressOp<double,3>("Compress",temp_h[i],&h[i]);
	for (j=0; j<FUNC_SIZE_M; j++)	c_op_g[j] = new CompressOp<double,3>("Compress",temp_g[j],&g[j]);

	for (k=0; k<SCALE; k++)
		for (i=0; i<FUNC_SIZE; i++)
			for (j=0; j<FUNC_SIZE_M; j++)
				inner_op_ug[k][i][j] = new InnerOp<double,3>("Inner",temp[i][j],temp_h[i],temp_g[j]);

	for (i=0; i<FUNC_SIZE; i++) 	sequence.push_back(c_op_h[i]);
	for (j=0; j<FUNC_SIZE_M; j++)	sequence.push_back(c_op_g[j]);

	for (k=0; k<SCALE; k++)
		for (i=0; i<FUNC_SIZE; i++)
			for (j=0; j<FUNC_SIZE_M; j++)
				sequence.push_back(inner_op_ug[k][i][j]);

	FuseT<double,3> odag(sequence);
	odag.processSequence();

	FusedOpSequence<double,3> fsequence = odag.getFusedOpSequence();
	FusedExecutor<double,3> fexecuter(world, &fsequence);
	fexecuter.execute();

	clkend = rtclock() - clkbegin;
	if (world.rank() == 0) printf("[the fused version by FuseT] Running Time: %f\n", clkend);
	world.gop.fence();

	for (k=0; k<SCALE; k++)
		for (i=0; i<FUNC_SIZE; i++)
			for (j=0; j<FUNC_SIZE_M; j++)	
				if (world.rank() == 0) printf ("[%d](%d,%d): %f\n",k,i,j,inner_op_ug[k][i][j]->_sum); 
//
//
//
	if (world.rank() == 0) print ("====================================================");
	if (world.rank() == 0) print ("==      MADNESS       ==============================");
	if (world.rank() == 0) print ("====================================================");
	world.gop.fence();

	clkbegin = rtclock();
	double resultInner[SCALE][FUNC_SIZE][FUNC_SIZE_M] = {0.0,};
	
	for (i=0; i<FUNC_SIZE; i++)
		h[i].compress();
	for (i=0; i<FUNC_SIZE_M; i++)
		g[i].compress();

	for (k=0; k<SCALE; k++)
		for (i=0; i<FUNC_SIZE; i++)
			for (j=0; j<FUNC_SIZE_M; j++)
				resultInner[k][i][j] = h[i].inner(g[j]);

	clkend = rtclock() - clkbegin;
	if (world.rank() == 0) printf("[MADNESS] Running Time: %f\n", clkend);
	world.gop.fence();

	for (k=0; k<SCALE; k++)
		for (i=0; i<FUNC_SIZE; i++)
			for (j=0; j<FUNC_SIZE_M; j++)	
				if (world.rank() == 0) printf ("[%d](%d,%d): %f\n",k,i,j,resultInner[k][i][j]); 

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

