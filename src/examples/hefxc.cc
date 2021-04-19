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



static const double Length = 50.0; // box size
static const long k = 8;        // wavelet order
static const double thresh = 1e-6; // precision

std::vector< std::shared_ptr<real_derivative_3d> > gradop;

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
//
//functionT make_dft_kernel(World & world, const vecfuncT& vf, int ispin,
//        int what)
//{
//    // xc_kernel has been disables and replaced by xc_kernel_apply, because
//    // it's more stable for GGA's
////        return multiop_values<double, xc_kernel, 3>(xc_kernel(xc, ispin, what), vf);
//}


static double guess(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return (2.0*exp(-sqrt(x*x+y*y+z*z+1e-8)));
}

// static double V(const coord_3d& r) {
//     const double x=r[0], y=r[1], z=r[2];
//     return  -2.0/sqrt(x*x+y*y+z*z+1e-8);
// }


// static double guess_density(const coord_3d& r) {
//     const double x=r[0], y=r[1], z=r[2];
//     return (4.0*exp(-2.0*sqrt(x*x+y*y+z*z+1e-6)))/(12.566358048);
// }


// //d/dx**2+d/dy**2+d/dz**2
// static double guess_gradient(const coord_3d& r) {
//     const double x=r[0], y=r[1], z=r[2];
//     double dx=(-8.0*x*exp(-2.0*sqrt(x*x+y*y+z*z+1e-6)))/(sqrt(z*z+y*y+x*x+1e-6)*12.566358048);
//     double dy=(-8.0*y*exp(-2.0*sqrt(x*x+y*y+z*z+1e-6)))/(sqrt(z*z+y*y+x*x+1e-6)*12.566358048);
//     double dz=(-8.0*z*exp(-2.0*sqrt(x*x+y*y+z*z+1e-6)))/(sqrt(z*z+y*y+x*x+1e-6)*12.566358048);
//     return dx*dx+dy*dy+dz*dz;
// }



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
    FunctionDefaults<3>::set_cubic_cell(-Length/2, Length/2);

    if (world.rank() == 0) print("\n  Solving for the KS aux. wave function\n");
    functionT psi  = real_factory_3d(world).f(guess);
    psi.truncate();
    psi.scale(1.0/psi.norm2());


    std::string xc_data;
    //xc_data="GGA_X_PBE 1.0 RHOTOL 1e-7 RHOMIN 0.0 SIGTOL 0.0 SIGMIN 0.0";

    //xc_data="GGA_X_PBE 1.0 RHOTOL 1e-7 RHOMIN 0.0 SIGTOL 1.e-10 SIGMIN 1.e-10";
    //xc_data="GGA_X_PBE 1.0 RHOTOL 1e-7 RHOMIN 0.0 SIGTOL 1.e-10 SIGMIN 1.e-10";

    //	xc_data="GGA_X_B88 1.0 RHOTOL 1e-5 RHOMIN 1e-4 SIGTOL 1e-3 SIGMIN 1e-5";
    //xc_data="GGA_X_B88 1.0 RHOTOL 1e-7 RHOMIN 0.0 SIGTOL 0.0 SIGMIN 0.0";
    //robert suggestions
    //xc_data="LDA_X 1.0 LDA_C_VWN 1.0";
    //xc_data="GGA_C_LYP 1.0 GGA_X_B88 1.0 RHOTOL 1e-7 RHOMIN 0.0 SIGTOL 0.0 SIGMIN 0.0";
    //xc_data="GGA_X_PBE 1.0 GGA_C_PBE 1.0 RHOTOL 0.0 RHOMIN 0.0 SIGTOL 1e-8 SIGMIN 1e-8";
    //xc_data="GGA_X_B88 1.0 RHOTOL 1e-7 RHOMIN 0.0 SIGTOL 0.0 SIGMIN 0.0";
    xc_data="GGA_C_LYP 1.0 GGA_X_B88 1.0 RHOTOL 1e-7 RHOMIN 0.0 SIGTOL 1e-10 SIGMIN 1e-10";
    //xc_data="GGA_X_PBE 1.";
    //xc_data="GGA_C_PBE 1.";
    //xc_data="GGA_X_B88 1.";
    //	xc.initialize(xc_data, false, world);
    XCOperator<double,3> xc(world,xc_data,false,psi,psi);

    gradop = gradient_operator<double,3>(world);

    {
        functionT rho = square(psi).truncate();
        rho.reconstruct();

        vecfuncT delrho;
        vecfuncT vf;

        vf.push_back(rho);

        // ADD SIGMA
//        functionT saa;
//        if (xc.xc->is_gga()) {
//            for(int axis = 0;axis < 3;++axis){
//                Derivative<double,3> D = free_space_derivative<double,3>(world, axis);
//                delrho.push_back(D(rho));
//            }
//            saa = delrho[0]*delrho[0]+delrho[1]*delrho[1]+delrho[2]*delrho[2];
//
//            vf.push_back(saa); // sigma_aa
//            if (vf.size()) {
//                reconstruct(world, vf);
//                //                       rho.refine_to_common_level(vf); // Ugly but temporary (I hope!)
//                refine_to_common_level(world,vf); // Ugly but temporary (I hope!)
//            }
//        }
//        double exc = make_dft_energy(world, vf, 0);
        double exc = xc.compute_xc_energy();
        print("exc=",exc );

//        real_function_3d  vxco = make_dft_potential(world, vf, 0, XCfunctional::potential_rho); //.truncate();
//#if 1
//        if (xc.is_gga() ) {
//            // get Vsigma_aa (if it is the case and Vsigma_bb)
//            functionT vsigaa = make_dft_potential(world, vf, 0, XCfunctional::potential_same_spin); //.truncate();
//
//            for (int axis=0; axis<3; axis++) {
//                functionT gradn = delrho[axis];
//                functionT ddel = vsigaa*gradn;
//                ddel.scale(4.0);
//                Derivative<double,3> D = free_space_derivative<double,3>(world, axis);
//                functionT vxc2=D(ddel);
//                vxco = vxco - vxc2;//.truncate();
//            }
//        } //is gga
//#endif

        real_function_3d vxco=xc.make_xc_potential();
//        real_function_3d vxc = real_factory_3d(world);
//        vxc.scale(0.0);
//
//
//        real_function_3d d1 = make_dft_kernel(world, vf, 0, 0); //.truncate();
//        d1.scale(constants::pi);
//        functionT vxc0 =  d1*rho;
//        vxc0.scale(2.);
//        vxc = vxc + vxc0;
//
//
//#if 1
//        if (xc.is_gga() ) {
//            // get Vsigma_aa (if it is the case and Vsigma_bb)
//            functionT d2 = make_dft_kernel(world, vf, 0, 1);
//            d2.scale(constants::pi);
//            //
//            functionT vxc1 = d2 * saa ;//  d2e/drds *s
//            vxc1.scale(1.*1.*4.);
//            vxc = vxc + vxc1;
//            //
//            functionT d3 = make_dft_kernel(world, vf, 0, 2);
//            d3.scale(constants::pi);
//            //
//            functionT d4 = make_dft_potential(world, vf, 0, XCfunctional::potential_same_spin);
//            d4.scale(constants::pi);
//            d4.scale(1.*2);
//            //
//            //
//            functionT fxct4 = d2 * rho;
//            fxct4.scale(1.*1.*2); //2 because rho
//            functionT fxct5 = d3 * saa ;
//            fxct5.scale(1.*1.*4); //4 cause saa
//            functionT fxct6 = fxct5 + fxct4;
//
//            for (int axis=0; axis<3; axis++) {
//                functionT gradn = delrho[axis ];
//
//                functionT ddel = (fxct6 + d4) * gradn;
//                ddel.scale(1.*1.*2);
//
//                Derivative<double,3> D = free_space_derivative<double,3>(world, axis);
//                functionT fxc2=D(ddel);
//                vxc = vxc + fxc2 ;//.truncate();
//                //
//            }
//        } //isgga
//
//        //  vxc.scale(4);
//#endif
        real_function_3d vxc=xc.apply_xc_kernel(rho);

        std::ofstream file;
        file.open ("fxc.txt");
        file.precision(12);
#if 0
        std::ofstream fd1;
        fd1.open ("fd1.txt");
        fd1.precision(12);
        std::ofstream fd2;
        fd2.open ("fd2.txt");
        fd2.precision(12);
        std::ofstream fd3;
        fd3.open ("fd3.txt");
        fd3.precision(12);
        std::ofstream fd4;
        fd4.open ("fd4.txt");
        fd4.precision(12);
#endif
        coord_3d r(0.0);
        /*double av, avrho, avpot;*/
        //double dx=0.1;
        //int imax=Length/dx;
        int imax=1024;
        for (int i=0; i<=imax; i++) {
            r[0] = -Length/2. + Length*i/imax;
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

            file << r[0] << "\t" << vxc(r) << "\t"<< "\t" << vxco(r) << "\t";
            //file << rho(r) ;
            //    fd1 << r[0] << "\t" << d1(r) << "\n" ;
            //    fd2 << r[0] << "\t" << d2(r) << "\n" ;
            //    fd3 << r[0] << "\t" << d3(r) << "\n" ;
            //    fd4 << r[0] << "\t" << d4(r) << "\n" ;
            /*
                   file << r[0] << "\t" << vxc(r) << "\t" << potential(r) << "\t" << Vnuc(r) << "\t";
                   file << rho(r) << "\t" << rho_nt(r) << "\t"  << rho_rec(r) << "\t" << psi(r);

             */

            file << "\n";
        }
        file.close();
#if 0
        fd1.close();
        fd2.close();
        fd3.close();
        fd4.close();
#endif

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



    world.gop.fence();
    finalize();
    return 0;
}



