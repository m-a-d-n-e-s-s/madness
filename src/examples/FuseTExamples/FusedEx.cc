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
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/constants.h>
#include <madness/mra/FuseT/CopyOp.h>
#include <madness/mra/FuseT/InnerOp.h>
#include <madness/mra/FuseT/CompressOp.h>
#include <madness/mra/FuseT/ReconstructOp.h>
#include <madness/mra/FuseT/FusedExecutor.h>
#include <madness/mra/FuseT/OpExecutor.h>
#include <madness/mra/FuseT/DerivativeOp.h>
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
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);

    startup(world, argc, argv);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_cubic_cell(-L, L);
    FunctionDefaults<3>::set_max_refine_level(4);
    real_function_3d u0 = real_factory_3d(world).f(uinitial);
    u0._treeName = "u0";
    u0.truncate();
    real_function_3d u1 = real_factory_3d(world).f(uinitial2);
    u1._treeName = "u1";
    u1.truncate();

    double u0_norm = u0.norm2();
    double u0_trace = u0.trace();
    double u1_norm = u1.norm2();
    double u1_trace = u1.trace();
	
    if (world.rank() == 0) print("Initial norm", u0_norm,"trace", u0_trace);
    if (world.rank() == 0) print("Initial norm", u1_norm,"trace", u1_trace);
    u0.compress();
    u1.compress();
    double in = u0.inner(u1);
    u1.reconstruct();
    world.gop.fence();

    // Make exponential of Vp
    real_function_3d rreconstruct_factory = real_factory_3d(world);
    real_function_3d rreconstruct(rreconstruct_factory);
    rreconstruct._treeName="RReconstruct";
    real_function_3d rcopy_factory = real_factory_3d(world);
    real_function_3d rcopy(rcopy_factory);
    rcopy._treeName="RCopy";
    real_function_3d rcmp_factory = real_factory_3d(world);
    real_function_3d rcmp(rcmp_factory);
    rcmp._treeName ="RCompress";
    real_function_3d rcmp1_factory = real_factory_3d(world);
    real_function_3d rcmp1(rcmp1_factory);
    rcmp1._treeName ="RCompress1";
    real_function_3d rinner_factory = real_factory_3d(world);
    real_function_3d rinner(rinner_factory);
    rinner._treeName ="RInner";

    
    ReconstructOp<double,3> op0("Reconstruct",&rreconstruct,&u0);
    CopyOp<double,3> op1("Copy",&rcopy,&rreconstruct);
    CompressOp<double,3> op2("Compress_rcopy",&rcmp,&rcopy);
    CompressOp<double,3> op3("Compress_u1",&rcmp1,&u1);
    InnerOp<double,3> op4("Inner",&rinner,&rcmp,&rcmp1);
    //InnerOp<double,3> op4("Inner",&rinner,&u0,&u1);
    
    bool fused = true;
    bool printInfo = true;
    
    if(fused){
	vector<PrimitiveOp<double,3>* > sequence;
	sequence.push_back(&op0);
	sequence.push_back(&op1);
	sequence.push_back(&op2);
	sequence.push_back(&op3);
	sequence.push_back(&op4);
	FuseT<double,3> odag(sequence);
	odag.processSequence();

	if(printInfo && world.rank() == 0){
	  odag.printOpsAndTrees();
	  odag.printValidSequences();
      }
	FusedOpSequence<double,3> fsequence = odag.getFusedOpSequence();
	FusedExecutor<double,3> fexecuter(world, &fsequence);
	fexecuter.execute();

	double rreconstruct_norm =rreconstruct.norm2();
	double rcopy_norm =rcopy.norm2();
	double rcmp_norm =rcmp.norm2();
	double rcmp1_norm =rcmp1.norm2();
	double rreconstruct_trace  =rreconstruct.trace();
	double rcopy_trace =rcopy.trace();
	double rcmp_trace =rcmp.trace();
	double rcmp1_trace =rcmp1.trace();

	//double in = u0.inner(u1);
	if(world.rank() == 0){
	    cout<<endl<<"Fused Execution Results :"<<endl;
	    cout<<"RReconstruct Norm : "<<rreconstruct_norm<<" and trace : "<<rreconstruct_trace<<endl;
	    cout<<"RCopy Norm : "<<rcopy_norm<<" and trace : "<<rcopy_trace<<endl;
	    cout<<"RCmp Norm : "<<rcmp_norm<<" and trace : "<<rcmp_trace<<endl;
	    cout<<"RCmp1 Norm : "<<rcmp1_norm<<" and trace : "<<rcmp1_trace<<endl;
	    cout<<"Madness Inner Product "<< in<<". FuseT Inner Product is "<<op4._sum<<endl;
	}
	
	
    }else{

	OpExecutor<double,3> exe(world);
	exe.execute(&op2, false);
	world.gop.fence();
    
    }

//	double result_norm = result.norm2();
//	double result_trace = u0.trace();
    
    //if (world.rank() == 0) print("Result norm", result_norm," result trace", result_trace);
    world.gop.fence();
    finalize();
    return 0;
}

