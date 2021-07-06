//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/nonlinsol.h>




//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

#include <madness/mra/lbdeux.h>
#include <madness/mra/qmprop.h>

#include <madness/misc/misc.h>
#include <madness/misc/ran.h>

#include <madness/tensor/systolic.h>
#include <madness/tensor/solvers.h>
#include <madness/tensor/elem.h>


#include <chem/xcfunctional.h>
#include <chem/SCFOperators.h>


using namespace madness;


typedef std::shared_ptr< WorldDCPmapInterface< Key<3> > > pmapT;
typedef Vector<double,3> coordT;
typedef std::shared_ptr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef std::vector<functionT> vecfuncT;
typedef std::pair<vecfuncT,vecfuncT> pairvecfuncT;
typedef std::vector<pairvecfuncT> subspaceT;
typedef Tensor<double> tensorT;
typedef FunctionFactory<double,3> factoryT;
typedef SeparatedConvolution<double,3> operatorT;
typedef std::shared_ptr<operatorT> poperatorT;
typedef Function<std::complex<double>,3> complex_functionT;
typedef std::vector<complex_functionT> cvecfuncT;
typedef Convolution1D<double_complex> complex_operatorT;


static const double R = 1.4;    // bond length
static const double L = 64.0*R; // box size
static const long k = 8;        // wavelet order
static const double thresh = 1e-5; // precision

//XCfunctional xc;
std::vector< std::shared_ptr<real_derivative_3d> > gradop;

//functionT make_dft_potential(World & world, XCfunctional& xc, const vecfuncT& vf,
//        int ispin, XCfunctional::xc_contribution what)
//{
//	return multiop_values<double, xc_potential, 3>(xc_potential(xc, ispin, what), vf);
//}
//
//double make_dft_energy(World & world, XCfunctional& xc, const vecfuncT& vf, int ispin)
//{
//	functionT vlda = multiop_values<double, xc_functional, 3>(xc_functional(xc), vf);
//	return vlda.trace();
//}

static double guess(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return (exp(-sqrt(x*x+y*y+(z-R/2)*(z-R/2)+1e-8))+
            exp(-sqrt(x*x+y*y+(z+R/2)*(z+R/2)+1e-8)));
}

static double V(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return -1.0/sqrt(x*x+y*y+(z-R/2)*(z-R/2)+1e-8)+
           -1.0/sqrt(x*x+y*y+(z+R/2)*(z+R/2)+1e-8);
}

//template <class solverT>
double iterate_ground(World& world, NonlinearSolver& solver,
		real_function_3d& V, real_function_3d& psi,
		double& eps) {
	real_convolution_3d op = BSHOperator3D(world, sqrt(-2*eps), 0.001, 1e-6);
	real_function_3d Vpsi = (V*psi);
	Vpsi.scale(-2.0).truncate();
	real_function_3d tmp = op(Vpsi).truncate();
	double norm = tmp.norm2();
	real_function_3d r = tmp-psi;
	double rnorm = r.norm2();
	double eps_new = eps - 0.5*inner(Vpsi,r)/(norm*norm);
	if (world.rank() == 0) {
		print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
	}
	psi = solver.update(psi, r);
	psi.scale(1.0/psi.norm2());
	eps = eps_new;
	return rnorm;
}

int main(int argc, char** argv) {
	initialize(argc, argv);
	World world(SafeMPI::COMM_WORLD);

	startup(world,argc,argv);
	std::cout.precision(6);

	FunctionDefaults<3>::set_k(k);
	FunctionDefaults<3>::set_thresh(thresh);
	FunctionDefaults<3>::set_refine(true);
	FunctionDefaults<3>::set_initial_level(5);
	FunctionDefaults<3>::set_truncate_mode(1);
	FunctionDefaults<3>::set_cubic_cell(-L/2, L/2);

	if (world.rank() == 0) print("\n  Solving for the KS aux. wave function\n");
	functionT Vnuc = real_factory_3d(world).f(V);
	functionT psi  = real_factory_3d(world).f(guess);
	psi.truncate();
	psi.scale(1.0/psi.norm2());

	std::string xc_data;
	xc_data="lda";
	//xc_data="GGA_X_PBE 1.";
	//xc_data="GGA_C_PBE 1.";
	//xc_data="GGA_X_B88 1.";
//	XCfunctional xc;
//	xc.initialize(xc_data, false, world);

	poperatorT coulop;
	coulop = poperatorT(CoulombOperatorPtr(world, 1e-10, thresh));
	gradop = gradient_operator<double,3>(world);
	real_convolution_3d op = CoulombOperator(world, 0.001, 1e-6);

	double eps = -0.5;
	{
	    NonlinearSolver solver;
	    for (int iter=0; iter<20; iter++) {
	        functionT rho = square(psi).truncate();

	        rho.reconstruct();
//	        vecfuncT delrho;
//	        vecfuncT vf(XCfunctional::number_xc_args);
//
//	        vf[XCfunctional::enum_rhoa]=rho;
//	        // ADD SIGMA
//	        if (xc.is_gga()) {
//	            for(int axis = 0;axis < 3;++axis){
//	                Derivative<double,3> D = free_space_derivative<double,3>(world, axis);
//	                delrho.push_back(D(rho));
//	            }
//	            functionT saa = delrho[0]*delrho[0]+delrho[1]*delrho[1]+delrho[2]*delrho[2];
//	            vf[XCfunctional::enum_saa]=(saa); // sigma_aa
//	            if (vf.size()) {
//	                reconstruct(world, vf);
//	                refine_to_common_level(world,vf); // Ugly but temporary (I hope!)
//	            }
//	        }
//	        //double exc = make_dft_energy(world,xc, vf, 0);
//	        //print("exc=",exc );
//
//	        real_function_3d vxc =  make_dft_potential(world,xc, vf, 0, XCfunctional::potential_rho);
//
//	        if (xc.is_gga()) {
//	            functionT vsigaa = make_dft_potential(world,xc, vf, 0, XCfunctional::potential_same_spin).truncate();
//
//	            for (int axis=0; axis<1; axis++) {
//	                Derivative<double,3> D = free_space_derivative<double,3>(world, axis);
//	                real_function_3d  gradn = D(rho);
//	                real_function_3d  ddel = vsigaa*gradn;
//	                ddel.scale(2.0);
//
//	                real_function_3d vxc2=D(ddel).truncate();
//	                vxc = vxc - vxc2;
//	            }
//	        }
	        XCOperator<double,3> xc(world,xc_data,false,rho,rho);
	        real_function_3d vxc=xc.make_xc_potential();


	        if(iter == 1000){
	            coord_3d r(0.0);
	            for (int i=0; i<201; i++) {
	                r[2] = -L/2. + L*i/200.0;
	                print(r[2], rho(r));
	                //print(r[2], vsigaa(r));
	                //print(r[0], vxc(r)*r[0]*2.);
	            }
	        }
	        real_function_3d potential = Vnuc + 2* op(rho).truncate() +vxc.truncate();

	        double err  = iterate_ground(world, solver, potential, psi, eps);
	        if (err < thresh) break;
	    }
	}
	double kinetic_energy = 0.0;
	for (int axis=0; axis<3; axis++) {
	    real_derivative_3d D = free_space_derivative<double,3>(world, axis);
	    real_function_3d dpsi = D(psi);
	    kinetic_energy += inner(dpsi,dpsi);
	}

	real_function_3d rho = square(psi);
	double two_electron_energy = 2*inner(op(rho),rho); // <u|rho> = <phi phi | 1/r12 | phi phi>
	rho.reconstruct();

//	vecfuncT vf;
//	vf.push_back(rho);
//	double exc=make_dft_energy(world,xc, vf, 0); // Exc
    XCOperator<double,3> xc(world,xc_data,false,rho,rho);
    real_function_3d vxc=xc.make_xc_potential();
    double exc=xc.compute_xc_energy();

	double nuclear_repulsion_energy = 1.0/R;
	double nuclear_attraction_energy = 2.0*inner(Vnuc,rho); // <V|rho> = <phi|V|phi>
	double total_energy = kinetic_energy + two_electron_energy +
	        nuclear_attraction_energy + exc + nuclear_repulsion_energy;
	double virial = (nuclear_attraction_energy + two_electron_energy ) / kinetic_energy;

	if (world.rank() == 0) {
	    print("");
	    print("            Kinetic energy ", kinetic_energy);
	    print(" Nuclear attraction energy ", nuclear_attraction_energy);
	    print(" Nuclear repulsion energy  ", nuclear_repulsion_energy);
	    print("       Two-electron energy ", two_electron_energy);
	    print("                 XC energy ", exc);
	    print("              Total energy ", total_energy);
	    print("                    Virial ", virial);
	}


	world.gop.fence();
	finalize();
	return 0;
}
