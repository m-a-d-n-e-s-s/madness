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
#include <madness/mra/FuseT/CopyOp.h>
#include <madness/mra/FuseT/CompressOp.h>
#include <madness/mra/FuseT/MatrixInnerOp.h>
#include <madness/mra/FuseT/OpExecutor.h>

using namespace madness;

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
    return exp(-alpha*(2*x*x+y*y+z*z))*pow(constants::pi/alpha,-1.5);
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
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);

    startup(world, argc, argv);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_cubic_cell(-L, L);

    real_function_3d u0 = real_factory_3d(world).f(uinitial);
    u0.truncate();

    double u0_norm = u0.norm2();
    double u0_trace = u0.trace();
	
    if (world.rank() == 0) print("Initial norm", u0_norm,"trace", u0_trace);
    world.gop.fence();

    // Make exponential of Vp
    real_function_3d result_factory = real_factory_3d(world);
    real_function_3d result(result_factory);

    double result_init_norm = result.norm2();
    double result_init_trace = result.trace();
	
    if (world.rank() == 0) print("Initial Result norm", result_init_norm,"trace", result_init_trace);
    world.gop.fence();

	CopyOp<double,3> op1("Copy",&result,&u0);
    OpExecutor<double,3> exe(world);
    exe.execute(&op1, false);
    world.gop.fence();

    u0_norm = u0.norm2();
    u0_trace = u0.trace();
	double result_norm = result.norm2();
	double result_trace = result.trace();
    
    if (world.rank() == 0) print("u0 norm", u0_norm," u0 trace", u0_trace);
    world.gop.fence();
    if (world.rank() == 0) print("Result norm", result_norm," result trace", result_trace);
    world.gop.fence();
   
	finalize(); 
    return 0;
}

