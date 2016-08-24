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
#include <madness/tensor/distributed_matrix.h>
#include <madness/tensor/cblas.h>
#include <madness/mra/FuseT/MatrixInnerOp.h>
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

static const double L		= 20;     // Half box size
static const long	k		= 8;        // wavelet order
static const double thresh	= 1e-6; // precision   // w/o diff. and 1e-12 -> 64 x 64
static const double c		= 2.0;       //
static const double tstep	= 0.1;
static const double alpha	= 1.9; // Exponent
static const double VVV		= 0.2;  // Vp constant value

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

#define FUNC_SIZE	32
#define FUNC_SIZE_M	32

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

static double random_function(const coord_3d& r) {
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


typedef DistributedMatrix<double> distmatT;
typedef Function<double,3> functionT;
typedef std::vector<functionT> vecfuncT;

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
	//FunctionDefaults<3>::set_max_refine_level(8);


	if (world.rank() == 0) print ("====================================================");
    if (world.rank() == 0) printf("   Initializing Functions\n");
    if (world.rank() == 0) printf("     %d Functions, %d Functions\n", FUNC_SIZE, FUNC_SIZE_M);
	if (world.rank() == 0) print ("====================================================");
    world.gop.fence();

	// 2 * N Functions
	real_function_3d  h[FUNC_SIZE];
	real_function_3d  g[FUNC_SIZE_M];
	real_function_3d  output[FUNC_SIZE*FUNC_SIZE_M];

	real_function_3d  temp_factory_h[FUNC_SIZE];
	real_function_3d  temp_factory_g[FUNC_SIZE_M];
	real_function_3d* temp_h[FUNC_SIZE];
	real_function_3d* temp_g[FUNC_SIZE_M];

	// N * N Results Functions by Inner-Product
	real_function_3d  temp_factory[FUNC_SIZE][FUNC_SIZE_M];
	real_function_3d* temp[FUNC_SIZE][FUNC_SIZE_M];

	int i, j;
	double clkbegin, clkend;
	clkbegin = rtclock();

	for (i=0; i<FUNC_SIZE; i++) 
	{
		randomizer();
		h[i]				= real_factory_3d(world).f(random_function);
		temp_factory_h[i]	= real_factory_3d(world);
		temp_h[i]			= new real_function_3d(temp_factory_h[i]);
	}

	for (i=0; i<FUNC_SIZE_M; i++)
	{
		randomizer();
		g[i]				= real_factory_3d(world).f(random_function);
//		temp_factory_g[i]	= real_factory_3d(world);
//		temp_g[i]			= new real_function_3d(temp_factory_g[i]);
	}

	for (i=0; i<FUNC_SIZE; i++)
		for (j=0; j<FUNC_SIZE_M; j++)
			output[i*FUNC_SIZE + j] = h[i]*g[j];
/*
	for (i=0; i<FUNC_SIZE; i++) 
		for (j=0; j<FUNC_SIZE_M; j++)
		{
			temp_factory[i][j]	= real_factory_3d(world);
			temp[i][j]			= new real_function_3d(temp_factory[i][j]);
		}
*/
	clkend = rtclock() - clkbegin;
	if (world.rank() == 0) printf("Running Time: %f\n", clkend);
	if (world.rank() == 0) print ("====================================================");
	if (world.rank() == 0) print ("==      FUSET-UNFUSED       ========================");
	if (world.rank() == 0) print ("====================================================");
	world.gop.fence();

	clkbegin = rtclock();
	for (i=0; i<FUNC_SIZE; i++)
		for (j=0; j<FUNC_SIZE_M; j++)
			output[i*FUNC_SIZE + j].compress();

	clkend = rtclock() - clkbegin;
	if (world.rank() == 0) printf("Running Time-- compress(): %f\n", clkend);
	world.gop.fence();
/*
	for (i=0; i<FUNC_SIZE; i++)
		h[i].compress();

	for (j=0; j<FUNC_SIZE_M; j++)
		g[j].compress();
*/
	vecfuncT fs;
	vecfuncT gs;

	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
		fs.push_back(output[i]);

	for (i=FUNC_SIZE*FUNC_SIZE_M/2; i<FUNC_SIZE*FUNC_SIZE_M; i++)
		gs.push_back(output[i]);

/*
	for (i=0; i<FUNC_SIZE; i++)
		fs.push_back(h[i]);

	for (j=0; j<FUNC_SIZE_M; j++)
		gs.push_back(g[j]);
*/
	clkbegin = rtclock();
	MatrixInnerOp<double,3>* matrix_inner_op = new MatrixInnerOp<double, 3>("MatrixInner", temp_h[0], fs, gs, false);

	OpExecutor<double,3> exe(world);
	exe.execute(matrix_inner_op, false);

	clkend = rtclock() - clkbegin;
	if (world.rank() == 0) printf("Running Time: %f\n", clkend);
	world.gop.fence();


	if (world.rank() == 0)
	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
		for (j=0; j<FUNC_SIZE_M*FUNC_SIZE_M/2; j++)	
		{	
		//	printf ("(%d,%d): %.12f\n", i, j, (*matrix_inner_op->_r)(i, j));
			printf ("(%d,%d): %f\n", i, j, (*matrix_inner_op->_r)(i, j));
		}
	world.gop.fence();

//
//
//
	if (world.rank() == 0) print ("====================================================");
	if (world.rank() == 0) print ("==      MADNESS - individual inner      ============");
	if (world.rank() == 0) print ("====================================================");
	world.gop.fence();

	vecfuncT v_f;
	vecfuncT v_g;
/*
	for (i=0; i<FUNC_SIZE; i++)
		v_f.push_back(h[i]);

	for (j=0; j<FUNC_SIZE_M; j++)
		v_g.push_back(g[j]);
*/
	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
		v_f.push_back(output[i]);

	for (i=FUNC_SIZE*FUNC_SIZE_M/2; i<FUNC_SIZE*FUNC_SIZE_M; i++)
		v_g.push_back(output[i]);

	clkbegin = rtclock();
	Tensor<double> ghaly = matrix_inner(world, v_f, v_g);
	//Tensor<double> ghaly = matrix_inner_old(world, v_f, v_g);

/*
	double resultInner[FUNC_SIZE][FUNC_SIZE_M] = {0.0, };
	for (i=0; i<FUNC_SIZE; i++)
		for (j=0; j<FUNC_SIZE_M; j++)
			resultInner[i][j] = h[i].inner(g[j]);
*/	

	if (world.rank() == 0) { 
		clkend = rtclock() - clkbegin;
		printf("Running Time: %f\n", clkend);
	}
	world.gop.fence();

	if (world.rank() == 0)
	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++) {
		for (j=0; j<FUNC_SIZE_M*FUNC_SIZE/2; j++){
			//printf ("matrix_inner_old: r(%d,%d): %.12f\n", i, j, ghaly(i, j));
			printf ("matrix_inner_old: r(%d,%d): %f\n", i, j, ghaly(i, j));
		}
	}
	world.gop.fence();

//
//
//
/*
	if (world.rank() == 0) print ("====================================================");
	if (world.rank() == 0) print ("==      MADNESS       ==============================");
	if (world.rank() == 0) print ("====================================================");
	world.gop.fence();

	clkbegin = rtclock();
	double resultInner[FUNC_SIZE][FUNC_SIZE_M] = {0.0,};
	
	for (i=0; i<FUNC_SIZE; i++)
		h[i].compress();
	for (i=0; i<FUNC_SIZE_M; i++)
		g[i].compress();

	for (i=0; i<FUNC_SIZE; i++)
		for (j=0; j<FUNC_SIZE_M; j++)
			resultInner[i][j] = h[i].inner(g[j]);

	clkend = rtclock() - clkbegin;
	if (world.rank() == 0) printf("Running Time: %f\n", clkend);
	world.gop.fence();
	
	for (i=0; i<FUNC_SIZE; i++)
		for (j=0; j<FUNC_SIZE_M; j++)	
		{
			if (world.rank() == 0) printf ("%d:%d = %f\n", i, j, resultInner[i][j]); 
		}
	world.gop.fence();
*/

	if (world.rank() == 0) print ("====================================================");
	if (world.rank() == 0) print ("==      GEMM          ==============================");
	if (world.rank() == 0) print ("====================================================");
	world.gop.fence();

	double* A;
	double* B;
	double* C;

	A = (double*)malloc(sizeof(double)*3*4);
	B = (double*)malloc(sizeof(double)*4*4);
	C = (double*)malloc(sizeof(double)*3*4);

	//
	for (i=0; i<4*3; i++)
		A[i] = i*1.0;

	//
	for (i=0; i<4; i++)
		B[i] = 1.0;

	for (i=4; i<8; i++)
		B[i] = 2.0;

	for (i=8; i<12; i++)
		B[i] = 3.0;

	for (i=12; i<16; i++)
		B[i] = 4.0;

	//
	for (i=0; i<3*4; i++)
		C[i] = 0.0;

	//
	cblas::gemm(cblas::CBLAS_TRANSPOSE::Trans, cblas::CBLAS_TRANSPOSE::NoTrans, 3, 4, 4, 1, A, 4, B, 4, 1, C, 3);

	if (world.rank() == 0)
	{
		for (i=0; i<3*4; i++) {
			printf ("%f ", C[i]);
		}
	}

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

