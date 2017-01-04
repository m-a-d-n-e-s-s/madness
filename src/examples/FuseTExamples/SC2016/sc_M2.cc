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
#include <madness/mra/FuseT/ReconstructOp.h>
#include <madness/mra/FuseT/OpExecutor.h>
#include <madness/mra/FuseT/DerivativeOp.h>
#include <madness/mra/FuseT/FusedExecutor.h>
#include <madness/mra/FuseT/FuseT.h>
#include <apps/chem/SCFOperators.h>
#include <apps/chem/SCF.h>

using namespace madness;

static const double L		= 20;     // Half box size
static const long	k		= 8;        // wavelet order
//static const double thresh	= 1e-6; // precision   // w/o diff. and 1e-12 -> 64 x 64
static const double c		= 2.0;       //
static const double tstep	= 0.1;
static const double alpha	= 1.9; // Exponent
static const double VVV		= 0.2;  // Vp constant value

#define VMRA_1
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

static double uinitial(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return exp(-alpha*(2*x*x+3.2*y*y+1.7*z*z))*pow(constants::pi/alpha,-1.5);
}
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

typedef Function<double,3> functionT;
typedef std::vector<functionT> vecfuncT;

struct write_test_input {

    double eprec=FunctionDefaults<3>::get_thresh()*0.1;

    std::string filename_;
    write_test_input(std::string mol="lih") : filename_("test_input") {
        std::ofstream of(filename_);
        of << "dft\n";
        of << "xc hf\n";
        of << "no_orient\n";
        of << "k 8\n";
        of << "protocol 1.e-5 \n";
        of << "nuclear_corrfac  slater 2.0\n";
        of << "end\n";

        if (mol=="lih") {
            of << "geometry\n";
            of << "eprec " << eprec << std::endl;
            of << "Li 0.0    0.0 0.0\n";
            of << "H  1.4375 0.0 0.0\n";
            of << "end\n";
        } else if (mol=="hf") {
            double eprec=1.e-5;
            of << "geometry\n";
            of << "eprec " << eprec << std::endl;
            of << "F  0.1    0.0 0.2\n";
            of << "H  1.4375 0.0 0.0\n";
            of << "end\n";
        }
        of.close();
    }

    ~write_test_input() {
        std::remove(filename_.c_str());
    }

    std::string filename() const {return filename_;}
};

int main(int argc, char** argv) 
{
	// input
	// (1) M			-- M and M functions
	// (2) thresh		-- threshold
	// (3) max-refine	-- max-refine-level
	// (4) type			-- 0: all, 1: FuseT, 2: vmra, 3: OpExecutor

	// M1. Kinetic Energy Matrix Calculation : vmra.h vs FusedExecutor (Reconstruct + DerivativeOp + CompressOp + InnerMatrixOp) 
	int		max_refine_level	= 14; //
	double	thresh				= 1e-12; // precision   // w/o diff. and 1e-12 -> 64 x 64
	int		FUNC_SIZE			= 4;
	int		FUNC_SIZE_M			= 4;
	int		type				= 0;

	if (argc == 5)
	{
		FUNC_SIZE			= atoi(argv[1]);
		FUNC_SIZE_M			= FUNC_SIZE;
		max_refine_level	= atoi(argv[3]);
		thresh				= atof(argv[2]);
		type				= atoi(argv[4]);
	}

    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);

    startup(world, argc, argv);

	if (world.rank() == 0) print ("====================================================");
	if (world.rank() == 0) printf("  Micro Benchmark #2 \n");
	if (world.rank() == 0) printf("  %d functions based on %d and %d random functions\n", FUNC_SIZE*FUNC_SIZE_M, FUNC_SIZE, FUNC_SIZE_M);
	if (world.rank() == 0) printf("  threshold: %13.4g, max_refine_level: %d\n", thresh, max_refine_level);
	if (world.rank() == 0) print ("====================================================");
	world.gop.fence();

	// Setting FunctionDefaults
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_cubic_cell(-L, L);
	FunctionDefaults<3>::set_max_refine_level(max_refine_level);


	if (world.rank() == 0) print ("====================================================");
    if (world.rank() == 0) print ("==   Initializing Functions   ======================");
	if (world.rank() == 0) print ("====================================================");
    world.gop.fence();

	// 2 * N Functions
	real_function_3d  h[FUNC_SIZE];
	real_function_3d  g[FUNC_SIZE_M];
	real_function_3d  output[FUNC_SIZE*FUNC_SIZE_M];

	real_function_3d  result_factory(world);
	real_function_3d  result(result_factory);

	int i, j;
	double clkbegin, clkend;
	clkbegin = rtclock();

	// M functions
	for (i=0; i<FUNC_SIZE; i++) {
		randomizer();
		h[i]				= real_factory_3d(world).f(random_function);
	}

	// N functions
	for (i=0; i<FUNC_SIZE_M; i++) {
		randomizer();
		g[i]				= real_factory_3d(world).f(random_function);
	}

	// M*N functions
	for (i=0; i<FUNC_SIZE; i++)
		for (j=0; j<FUNC_SIZE_M; j++)
			output[i*FUNC_SIZE_M + j] = h[i]*g[j];

	// compressing M*N functions
	for (i=0; i<FUNC_SIZE; i++)
		for (j=0; j<FUNC_SIZE_M; j++)
			output[i*FUNC_SIZE_M + j].compress();

#ifdef VMRA
	// M1. Kinetic Energy Matrix Calculation : vmra.h vs FusedExecutor (Reconstruct + DerivativeOp + CompressOp + InnerMatrixOp) 
	if (world.rank() == 0) print ("====================================================");
	if (world.rank() == 0) print ("=== Kinetic Energy Matrix Calculation ==============");
	if (world.rank() == 0) print ("=== vmra.h =========================================");
	if (world.rank() == 0) print ("====================================================");
	world.gop.fence();

	distmatT r;
	vecfuncT v_f;
	vecfuncT v_g;

    write_test_input test_input;
	SCF temp(world,test_input.filename().c_str());

	clkbegin = rtclock();
	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
		v_f.push_back(output[i]);

	for (i=FUNC_SIZE*FUNC_SIZE_M/2; i<FUNC_SIZE*FUNC_SIZE_M; i++)
		v_g.push_back(output[i]);

	r = temp.kinetic_energy_matrix(world, v_f, v_g);

	clkend = rtclock() - clkbegin;
	if (world.rank() == 0)	printf("Running Time: %f\n", clkend);
	world.gop.fence();
#endif
//
//
//
	if (world.rank() == 0) print ("====================================================");
	if (world.rank() == 0) print ("=== FusedExecutor (Reconstruct + Derivative + ======");
	if (world.rank() == 0) print ("===                Compress + MatrixInner)    ======");
	if (world.rank() == 0) print ("====================================================");
	world.gop.fence();
/*
	distmatT SCF::kinetic_energy_matrix(World & world, const vecfuncT & vbra, const vecfuncT & vket) const {
         PROFILE_MEMBER_FUNC(SCF);
         MADNESS_ASSERT(vbra.size() == vket.size());
         int n = vbra.size();
         distmatT r = column_distributed_matrix<double>(world, n, n);
         reconstruct(world, vbra);
         reconstruct(world, vket);
     }
*/
	// Results for ReconstructOp
	real_function_3d  reconstruct_factory_h[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d  reconstruct_factory_g[FUNC_SIZE_M*FUNC_SIZE/2];
	real_function_3d* reconstruct_h[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d* reconstruct_g[FUNC_SIZE_M*FUNC_SIZE/2];

	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
	{
		reconstruct_factory_h[i]	= real_factory_3d(world);
		reconstruct_factory_g[i]	= real_factory_3d(world);
		reconstruct_h[i]			= new real_function_3d(reconstruct_factory_h[i]);
		reconstruct_g[i]			= new real_function_3d(reconstruct_factory_g[i]);
	}

	// Results for DerivativeOp
	real_function_3d  derivative_factory_h_x[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d  derivative_factory_h_y[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d  derivative_factory_h_z[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d  derivative_factory_g_x[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d  derivative_factory_g_y[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d  derivative_factory_g_z[FUNC_SIZE*FUNC_SIZE_M/2];

	real_function_3d* derivative_h_x[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d* derivative_h_y[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d* derivative_h_z[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d* derivative_g_x[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d* derivative_g_y[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d* derivative_g_z[FUNC_SIZE*FUNC_SIZE_M/2];

	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
	{
		derivative_factory_h_x[i] = real_factory_3d(world);
		derivative_factory_h_y[i] = real_factory_3d(world);
		derivative_factory_h_z[i] = real_factory_3d(world);
		derivative_factory_g_x[i] = real_factory_3d(world);
		derivative_factory_g_y[i] = real_factory_3d(world);
		derivative_factory_g_z[i] = real_factory_3d(world);

		derivative_h_x[i]	= new real_function_3d(derivative_factory_h_x[i]);
		derivative_h_y[i]	= new real_function_3d(derivative_factory_h_y[i]);
		derivative_h_z[i]	= new real_function_3d(derivative_factory_h_z[i]);
		derivative_g_x[i]	= new real_function_3d(derivative_factory_g_x[i]);
		derivative_g_y[i]	= new real_function_3d(derivative_factory_g_y[i]);
		derivative_g_z[i]	= new real_function_3d(derivative_factory_g_z[i]);
	}

	real_derivative_3d D_h_x	= free_space_derivative<double,3>(world,0);
	real_derivative_3d D_h_y	= free_space_derivative<double,3>(world,1);
	real_derivative_3d D_h_z	= free_space_derivative<double,3>(world,2);
	//real_derivative_3d D_g_x	= free_space_derivative<double,3>(world,0);
	//real_derivative_3d D_g_y	= free_space_derivative<double,3>(world,1);
	//real_derivative_3d D_g_z	= free_space_derivative<double,3>(world,2);

	// Results for CompressOp
	real_function_3d  compress_factory_h_x[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d  compress_factory_h_y[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d  compress_factory_h_z[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d  compress_factory_g_x[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d  compress_factory_g_y[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d  compress_factory_g_z[FUNC_SIZE*FUNC_SIZE_M/2];

	real_function_3d* compress_h_x[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d* compress_h_y[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d* compress_h_z[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d* compress_g_x[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d* compress_g_y[FUNC_SIZE*FUNC_SIZE_M/2];
	real_function_3d* compress_g_z[FUNC_SIZE*FUNC_SIZE_M/2];

	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
	{
		compress_factory_h_x[i] = real_factory_3d(world);
		compress_factory_h_y[i] = real_factory_3d(world);
		compress_factory_h_z[i] = real_factory_3d(world);
		compress_factory_g_x[i] = real_factory_3d(world);
		compress_factory_g_y[i] = real_factory_3d(world);
		compress_factory_g_z[i] = real_factory_3d(world);

		compress_h_x[i] = new real_function_3d(compress_factory_h_x[i]);
		compress_h_y[i] = new real_function_3d(compress_factory_h_y[i]);
		compress_h_z[i] = new real_function_3d(compress_factory_h_z[i]);
		compress_g_x[i] = new real_function_3d(compress_factory_g_x[i]);
		compress_g_y[i] = new real_function_3d(compress_factory_g_y[i]);
		compress_g_z[i] = new real_function_3d(compress_factory_g_z[i]);
	}

	// Results for MatrixInnerOp
	real_function_3d  matrixinner_factory_x(world);
	real_function_3d  matrixinner_factory_y(world);
	real_function_3d  matrixinner_factory_z(world);
	real_function_3d  matrixinner_x(matrixinner_factory_x);
	real_function_3d  matrixinner_y(matrixinner_factory_y);
	real_function_3d  matrixinner_z(matrixinner_factory_z);

	clkbegin = rtclock();

	// Reconstruct Op
	ReconstructOp<double,3>* reconstruct_op_h[FUNC_SIZE*FUNC_SIZE_M/2];	// vbra
	ReconstructOp<double,3>* reconstruct_op_g[FUNC_SIZE*FUNC_SIZE_M/2];	// bket

	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
	{
		reconstruct_op_h[i]	= new ReconstructOp<double,3>("ReconstructOp", reconstruct_h[i], &output[i]);
		reconstruct_op_g[i]	= new ReconstructOp<double,3>("ReconstructOp", reconstruct_g[i], &output[i + (FUNC_SIZE*FUNC_SIZE_M/2)]);
	}
	
/*
	vecfuncT dvx_bra = apply(world, *(gradop[0]), vbra, false);
	vecfuncT dvy_bra = apply(world, *(gradop[1]), vbra, false);
	vecfuncT dvz_bra = apply(world, *(gradop[2]), vbra, false);
	vecfuncT dvx_ket = apply(world, *(gradop[0]), vket, false);
	vecfuncT dvy_ket = apply(world, *(gradop[1]), vket, false);
	vecfuncT dvz_ket = apply(world, *(gradop[2]), vket, false);
	world.gop.fence();
*/

	// Derivative Op
	DerivativeOp<double,3>* derivative_op_x_b[FUNC_SIZE*FUNC_SIZE_M/2];
	DerivativeOp<double,3>* derivative_op_y_b[FUNC_SIZE*FUNC_SIZE_M/2];
	DerivativeOp<double,3>* derivative_op_z_b[FUNC_SIZE*FUNC_SIZE_M/2];
	DerivativeOp<double,3>* derivative_op_x_k[FUNC_SIZE*FUNC_SIZE_M/2];
	DerivativeOp<double,3>* derivative_op_y_k[FUNC_SIZE*FUNC_SIZE_M/2];
	DerivativeOp<double,3>* derivative_op_z_k[FUNC_SIZE*FUNC_SIZE_M/2];

	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
	{
		derivative_op_x_b[i] = new DerivativeOp<double,3>("Derivative00",derivative_h_x[i],reconstruct_h[i], world,&D_h_x);
		derivative_op_y_b[i] = new DerivativeOp<double,3>("Derivative01",derivative_h_y[i],derivative_h_x[i],world,&D_h_y);
		derivative_op_z_b[i] = new DerivativeOp<double,3>("Derivative02",derivative_h_z[i],derivative_h_y[i],world,&D_h_z);
	}

	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
	{
		derivative_op_x_k[i] = new DerivativeOp<double,3>("Derivative10",derivative_g_x[i],reconstruct_g[i], world,&D_h_x);
		derivative_op_y_k[i] = new DerivativeOp<double,3>("Derivative11",derivative_g_y[i],derivative_g_x[i],world,&D_h_y);
		derivative_op_z_k[i] = new DerivativeOp<double,3>("Derivative12",derivative_g_z[i],derivative_g_y[i],world,&D_h_z);
	}


	// Compress Op
	CompressOp<double,3>* compress_op_x_b[FUNC_SIZE];
	CompressOp<double,3>* compress_op_y_b[FUNC_SIZE];
	CompressOp<double,3>* compress_op_z_b[FUNC_SIZE];
	CompressOp<double,3>* compress_op_x_k[FUNC_SIZE];
	CompressOp<double,3>* compress_op_y_k[FUNC_SIZE];
	CompressOp<double,3>* compress_op_z_k[FUNC_SIZE];

	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
	{
		compress_op_x_b[i] = new CompressOp<double,3>("CompressOp",compress_h_x[i],derivative_h_x[i]);
		compress_op_y_b[i] = new CompressOp<double,3>("CompressOp",compress_h_y[i],derivative_h_y[i]);
		compress_op_z_b[i] = new CompressOp<double,3>("CompressOp",compress_h_z[i],derivative_h_z[i]);
		compress_op_x_k[i] = new CompressOp<double,3>("CompressOp",compress_g_x[i],derivative_g_x[i]);
		compress_op_y_k[i] = new CompressOp<double,3>("CompressOp",compress_g_y[i],derivative_g_y[i]);
		compress_op_z_k[i] = new CompressOp<double,3>("CompressOp",compress_g_z[i],derivative_g_z[i]);
	}

/*
	r += matrix_inner(r.distribution(), dvx_bra, dvx_ket, true);
	r += matrix_inner(r.distribution(), dvy_bra, dvy_ket, true);
	r += matrix_inner(r.distribution(), dvz_bra, dvz_ket, true);
	r *= 0.5;
	return r;
*/
	
	// MatrixInner Op
	vecfuncT h_x;
	vecfuncT h_y;
	vecfuncT h_z;
	vecfuncT g_x;
	vecfuncT g_y;
	vecfuncT g_z;

	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
	{
		h_x.push_back(*compress_h_x[i]);	
		h_y.push_back(*compress_h_y[i]);	
		h_z.push_back(*compress_h_z[i]);	
		g_x.push_back(*compress_g_x[i]);	
		g_y.push_back(*compress_g_y[i]);	
		g_z.push_back(*compress_g_z[i]);	
	}

	MatrixInnerOp<double,3>* matrixinner_op_a = new MatrixInnerOp<double,3>("MatrixInner", &matrixinner_x, h_x, g_x, true, false);
	MatrixInnerOp<double,3>* matrixinner_op_b = new MatrixInnerOp<double,3>("MatrixInner", &matrixinner_y, h_y, g_y, true, false);
	MatrixInnerOp<double,3>* matrixinner_op_c = new MatrixInnerOp<double,3>("MatrixInner", &matrixinner_z, h_z, g_z, true, false);

	clkend = rtclock() - clkbegin;
	if (world.rank() == 0)	printf("Running Time: %f\n", clkend);
	world.gop.fence();
//
//
	clkbegin = rtclock();
	vector<PrimitiveOp<double,3>*>	sequence;

	// Pushing ReconstructOp
	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
		sequence.push_back(reconstruct_op_h[i]);
	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
		sequence.push_back(reconstruct_op_g[i]);

	// Pushing DerivativeOp

	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
	{
		sequence.push_back(derivative_op_x_k[i]);
		//sequence.push_back(derivative_op_x_b[i]);
	}
	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
	{
		sequence.push_back(derivative_op_y_k[i]);
		sequence.push_back(derivative_op_z_k[i]);
	}
	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
	{
		//sequence.push_back(derivative_op_y_b[i]);
		//sequence.push_back(derivative_op_z_b[i]);
	}
/*
	// Pushing CompressOp
	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
		sequence.push_back(compress_op_x_b[i]);
	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
		sequence.push_back(compress_op_y_b[i]);
	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
		sequence.push_back(compress_op_z_b[i]);
	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
		sequence.push_back(compress_op_x_k[i]);
	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
		sequence.push_back(compress_op_y_k[i]);
	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++)
		sequence.push_back(compress_op_z_k[i]);

	// Pushing MatrixInnerOp
	//MatrixInnerOp<double,3>* matrixinner_op_a = new MatrixInnerOp<double,3>("MatrixInner", &matrixinner_x, h_x, g_x, true);
	//MatrixInnerOp<double,3>* matrixinner_op_b = new MatrixInnerOp<double,3>("MatrixInner", &matrixinner_y, h_y, g_y, true);
	//MatrixInnerOp<double,3>* matrixinner_op_c = new MatrixInnerOp<double,3>("MatrixInner", &matrixinner_z, h_z, g_z, true);

	sequence.push_back(matrixinner_op_a);
	sequence.push_back(matrixinner_op_b);
	sequence.push_back(matrixinner_op_c);
*/


	// Processing a sequence of Operators
	FuseT<double,3> odag(sequence);
	odag.processSequence();

	FusedOpSequence<double,3> fsequence = odag.getFusedOpSequence();
	FusedExecutor<double,3> fexecuter(world, &fsequence);
	fexecuter.execute();
		
	clkend = rtclock() - clkbegin;
	if (world.rank() == 0) printf ("Done!\n");
	if (world.rank() == 0)	printf("Running Time: %f\n", clkend);
	world.gop.fence();
/*
	if (world.rank() == 0)
	for (i=0; i<FUNC_SIZE*FUNC_SIZE_M/2; i++) {
		for (j=0; j<FUNC_SIZE_M*FUNC_SIZE/2; j++){
			printf ("inner fused: r(%d,%d): %f\n", i, j, inner_op_ug[i][j]->_sum);
		}
	}
	world.gop.fence();
*/

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

