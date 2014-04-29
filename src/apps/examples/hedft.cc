#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/operator.h>
#include <examples/nonlinsol.h>




#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/operator.h>

#include <mra/lbdeux.h>
#include <mra/qmprop.h>

#include <misc/misc.h>
#include <misc/ran.h>

#include <tensor/systolic.h>
#include <tensor/solvers.h>
#include <tensor/elem.h>


#include <moldft/xcfunctional.h>


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



static const double L = 50.0; // box size
static const long k = 8;        // wavelet order
static const double thresh = 1e-8; // precision

std::vector< std::shared_ptr<real_derivative_3d> > gradop;
XCfunctional xc;

double make_dft_energy(World & world, const vecfuncT& vf, int ispin)
{
	functionT vlda = multiop_values<double, xc_functional, 3>(xc_functional(xc, ispin), vf);
	return vlda.trace();
}

functionT make_dft_potential(World & world, const vecfuncT& vf, int ispin, int what)
{
	return multiop_values<double, xc_potential, 3>(xc_potential(xc, ispin, what), vf);
}

static double guess(const coord_3d& r) {
	const double x=r[0], y=r[1], z=r[2];
	return (2.0*exp(-sqrt(x*x+y*y+z*z+1e-8)));
}

static double V(const coord_3d& r) {
	const double x=r[0], y=r[1], z=r[2];
	return  -2.0/sqrt(x*x+y*y+z*z+1e-8);
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
	xc.initialize(xc_data, false);

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
			vecfuncT delrho;
                        vecfuncT vf;

			vf.push_back(rho);
// ADD SIGMA
			if (xc.is_gga()) {
                              for(int axis = 0;axis < 3;++axis){
                                       Derivative<double,3> D = free_space_derivative<double,3>(world, axis);
                                       delrho.push_back(D(rho));
                              }
                              functionT saa = delrho[0]*delrho[0]+delrho[1]*delrho[1]+delrho[2]*delrho[2];
                              vf.push_back(saa); // sigma_aa
                              if (vf.size()) {
                                    reconstruct(world, vf);
                                    rho.refine_to_common_level(vf); // Ugly but temporary (I hope!)
                              }
                        }
                        //double exc = make_dft_energy(world, vf, 0);
                        //print("exc=",exc );

                        real_function_3d vxc =  make_dft_potential(world, vf, 0, 0);

			if (xc.is_gga()) {
                                functionT vsigaa = make_dft_potential(world, vf, 0, 1).truncate();

                                for (int axis=0; axis<1; axis++) {
                                        Derivative<double,3> D = free_space_derivative<double,3>(world, axis);
                                        real_function_3d  gradn = D(rho);
                                        real_function_3d  ddel = vsigaa*gradn;
                                        ddel.scale(2.0);

                                        real_function_3d vxc2=D(ddel).truncate();
                                        vxc = vxc - vxc2;
                                }
                        }


                        if(iter == 0){
			   coord_3d r(0.0);
			   for (int i=0; i<201; i++) {
			      //r[2] = -L/2. + L*i/200.0;
			      r[0] =  20.*i/400.0;
			      //print(r[2], vf[1](r));
			      //print(r[2], vsigaa(r));
		//	      print(r[0], vxc(r)*r[0]);
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
	vecfuncT vf;
	vf.push_back(rho);
	double exc=make_dft_energy(world, vf, 0); // Exc
	double nuclear_attraction_energy = 2.0*inner(Vnuc,rho); // <V|rho> = <phi|V|phi>
	double total_energy = kinetic_energy + two_electron_energy + 
		nuclear_attraction_energy + exc;
	double virial = (nuclear_attraction_energy + two_electron_energy ) / kinetic_energy;

	if (world.rank() == 0) {
		print("");
		print("            Kinetic energy ", kinetic_energy);
		print(" Nuclear attraction energy ", nuclear_attraction_energy);
		print("       Two-electron energy ", two_electron_energy);
		print("                 XC energy ", exc);
		print("              Total energy ", total_energy);
		print("                    Virial ", virial);
	}


	world.gop.fence();
	finalize();
	return 0;
}
