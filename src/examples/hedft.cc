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

#include <madness/mra/legendre.h>


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



static const double length = 50.0; // box size
static const long k = 8;        // wavelet order
static const double thresh = 1e-8; // precision

std::vector< std::shared_ptr<real_derivative_3d> > gradop;
//XCfunctional xc;

//double make_dft_energy(World & world, const vecfuncT& vf, int ispin)
//{
//	functionT vlda = multiop_values<double, xc_functional, 3>(xc_functional(xc), vf);
//	return vlda.trace();
//}
//
//functionT make_dft_potential(World & world, const vecfuncT& vf, int ispin,
//        XCfunctional::xc_contribution what)
//{
//	return multiop_values<double, xc_potential, 3>(xc_potential(xc, ispin, what), vf);
//}

static double guess(const coord_3d& r) {
	const double x=r[0], y=r[1], z=r[2];
	return (2.0*exp(-sqrt(x*x+y*y+z*z+1e-8)));
}

static double V(const coord_3d& r) {
	const double x=r[0], y=r[1], z=r[2];
	return  -2.0/sqrt(x*x+y*y+z*z+1e-8);
}


/*static double guess_density(const coord_3d& r) {
	const double x=r[0], y=r[1], z=r[2];
	return (4.0*exp(-2.0*sqrt(x*x+y*y+z*z+1e-6)))/(12.566358048);
}

//d/dx**2+d/dy**2+d/dz**2
static double guess_gradient(const coord_3d& r) {
	const double x=r[0], y=r[1], z=r[2];
        double dx=(-8.0*x*exp(-2.0*sqrt(x*x+y*y+z*z+1e-6)))/(sqrt(z*z+y*y+x*x+1e-6)*12.566358048);
        double dy=(-8.0*y*exp(-2.0*sqrt(x*x+y*y+z*z+1e-6)))/(sqrt(z*z+y*y+x*x+1e-6)*12.566358048);
        double dz=(-8.0*z*exp(-2.0*sqrt(x*x+y*y+z*z+1e-6)))/(sqrt(z*z+y*y+x*x+1e-6)*12.566358048);
	return dx*dx+dy*dy+dz*dz;
}*/


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
	std::cout.precision(12);

	FunctionDefaults<3>::set_k(k);
	FunctionDefaults<3>::set_thresh(thresh);
	FunctionDefaults<3>::set_refine(true);
	FunctionDefaults<3>::set_initial_level(5);
	FunctionDefaults<3>::set_truncate_mode(1);
	FunctionDefaults<3>::set_cubic_cell(-length/2, length/2);

	if (world.rank() == 0) print("\n  Solving for the KS aux. wave function\n");
	functionT Vnuc = real_factory_3d(world).f(V);
	functionT psi  = real_factory_3d(world).f(guess);
	psi.truncate();
	psi.scale(1.0/psi.norm2());


	std::string xc_data;
	//xc_data="LDA_X 1.0 LDA_C_VWN 1.0";
	xc_data="GGA_X_PBE 1.0 GGA_C_PBE 1.0";
	//xc_data="GGA_X_PBE 1.";
	//xc_data="GGA_C_PBE 1.";
	//xc_data="GGA_X_B88 1.";
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
	        //functionT rho  = real_factory_3d(world).f(guess_density);
	        //rho.truncate();
	        rho.reconstruct();

//	        vecfuncT delrho;
//	        vecfuncT vf(XCfunctional::number_xc_args);    // number of possible intermediates
//
////	        vf.push_back(rho);
//	        vf[XCfunctional::enum_rhoa]=rho;
//
//	        // ADD SIGMA
//	        functionT saa;
//	        if (xc.is_gga()) {
//	            for(int axis = 0;axis < 3;++axis){
//	                Derivative<double,3> D = free_space_derivative<double,3>(world, axis);
//	                delrho.push_back(D(rho));
//	            }
//	            saa = delrho[0]*delrho[0]+delrho[1]*delrho[1]+delrho[2]*delrho[2];
//	            //saa  = real_factory_3d(world).f(guess_gradient);
//	            //saa.truncate();
//	            //saa.reconstruct();
//
//	            vf[XCfunctional::enum_saa]=saa;
//
//	            reconstruct(world, vf);
//	            refine_to_common_level(world,vf); // Ugly but temporary (I hope!)
//	        }
//	        //double exc = make_dft_energy(world, vf, 0);
//	        //print("exc=",exc );
//
//	        real_function_3d vxc =  make_dft_potential(world, vf, 0, XCfunctional::potential_rho);
//
//
//	        functionT vsigaa;
//	        if (xc.is_gga()) {
//	            vsigaa = make_dft_potential(world, vf, 0, XCfunctional::potential_same_spin).truncate();
//
//	            for (int axis=0; axis<3; axis++) {
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

	        real_function_3d potential = Vnuc + 2* op(rho).truncate() +vxc.truncate();

	        //deactivate test printing for now
	        if(iter == 0 && false){
	            std::ofstream file;
	            std::ostringstream iter_str;
	            iter_str << iter;

	            file.open ("iter"+iter_str.str()+".txt");
	            file.precision(12);
	            coord_3d r(0.0);
	            /*double av, avrho, avpot;*/
	            //double dx=0.1;
	            //int imax=length/dx;
	            int imax=1024;
	            for (int i=0; i<=imax; i++) {
	                r[0] = -length/2. + length*i/imax;
	                /*av=0.0;avrho=0.0;avpot=0.0;
			      for (int j=0; j<=imax; j++) {
			         r[1] = -L/2. + L*j/imax;
			         for (int k=0; k<=imax; k++) {
			            r[0] = -L/2. + L*k/imax;

			            //print(r[2], vf[1](r));
			            //print(r[2], vsigaa(r));
			            //print(r[0], vxc(r)*r[0]);
			            //print(r[0], vxc(r)*r[0], r[2], vf[1](r), r[2], vsigaa(r));
			            //if (xc.is_gga())
			            //   {file << r[0] << "\t" << vxc(r)*r[0] << "\t" << r[2] << "\t" << vf[1](r) << "\t" << r[2] << "\t" << vsigaa(r) << "\n";}
			            //else
			            //{file << r[0] << "\t" << r[1] << "\t" << r[2] << "\t" << vxc(r) << "\t" << vxc(r)*r[0] << "\t" << vxc(r)*r[1] << "\t" << vxc(r)*r[2] << "\n" ;}

			            av+=vxc(r);
			            avrho+=rho(r);
			            avpot+=potential(r);
			         }
			      }*/
	                r[2]=0.01;r[1]=0.01;
	                /*file << r[2] << "\t" << vxc(r) << "\t" << vxc_trunc(r) << "\t" << av/((imax+1)*(imax+1)) << "\t";
                              file << potential(r) << "\t" << avpot/((imax+1)*(imax+1)) << "\t" << Vnuc(r) << "\t"  << rho(r) << "\t" << avrho/((imax+1)*(imax+1));
                              file << "\t" << psi(r) << "\n" ;*/

	                //functionT rho_nt = square(psi);
	                //functionT rho_rec = square(psi).truncate();
	                //rho_rec.reconstruct();
	                functionT rho_nt = real_factory_3d(world);
	                functionT rho_rec = real_factory_3d(world);

	                file << r[0] << "\t" << vxc(r) << "\t" << potential(r) << "\t" << Vnuc(r) << "\t";
	                file << rho(r) << "\t" << rho_nt(r) << "\t"  << rho_rec(r) << "\t" << psi(r);

	                file << "\n";
	            }
	            file.close();

	            ///////////////////////////////////////////////////////////
	            /*{
                              int npt_plot=101;
                              tensorT plot_cell;
                              //plot_cell = tensorT(3L,2L);
                              std::vector<long> npt(3,npt_plot);

                              if (plot_cell.size() == 0)
                                  plot_cell = copy(FunctionDefaults<3>::get_cell());

                                  plotdx(rho, "density.dx", plot_cell, npt, true);
                                  plotdx(Vnuc, "vnuc.dx", plot_cell, npt, true);

                           }*/
	            ///////////////////////////////////////////////////////////

	        }

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
//	vecfuncT vf(13);
//	vf.push_back(rho);
//    vf[XCfunctional::enum_rhoa]=rho;

//	vecfuncT delrho;
    XCOperator<double,3> xc(world,xc_data,false,rho,rho);
    real_function_3d vxc=xc.make_xc_potential();

//	if (xc.is_gga()) {
//	    for(int axis = 0;axis < 3;++axis){
//	        Derivative<double,3> D = free_space_derivative<double,3>(world, axis);
//	        delrho.push_back(D(rho));
//	    }
//	    functionT saa = delrho[0]*delrho[0]+delrho[1]*delrho[1]+delrho[2]*delrho[2];
////	    vf.push_back(saa); // sigma_aa
//	    vf[XCfunctional::enum_saa]=saa;
//
//	    reconstruct(world, vf);
//	    refine_to_common_level(world,vf); // Ugly but temporary (I hope!)
//	}
//    double exc=make_dft_energy(world, vf, 0); // Exc
	double exc=xc.compute_xc_energy(); // Exc
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



        //print out scaling functions
        /*{
        std::ofstream file;
 	file.open ("scalingc.txt");
	file.precision(12);
        double phi[k];
        double dx=0.01;
        int l=1;
        int imax=2*l/dx;
        double x=-l;
        for (int n=0;n<=64;n++) {
        x=-1;
        for (int i=0; i<=imax; i++) {
          legendre_polynomials(x,k,phi);
          file << x+n*2*l;
          for (int j=0;j<k;j++) {
            file << "\t" << phi[j];
          }
          file << "\n";
          x+=dx;
        }
        }
        file.close();
        }*/

	world.gop.fence();
	finalize();
	return 0;
}



