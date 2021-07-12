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
*/

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness.h>
#include <chem/SCFOperators.h>

using namespace madness;

static double alpha=10.0;   // exponent of the density

// slater function centered at the origin
double slater(const coord_3d& xyz) {
    return exp(-alpha*xyz.normf());
}
// slater function centered at the origin
double slater2(const coord_3d& xyz) {
    double s=slater(xyz);
    return s*s;
}

bool check_err(double err, double thresh, std::string msg) {
    if (fabs(err)>thresh) {
        print("\nfailing test:",msg,"\n");
        return true;
    }
    return false;
}


/// test the Slater exchange potential and the kernel

/// The Slater exchange is a simple function of rho,
/// \f[
///    E_\mathrm{xc} = \int\epsilon_\mathrm{xc}[\rho]\rho(r) d^3r
/// \f]
/// with
/// \f[
///    v_\mathrm{xc} = \epsilon_\mathrm{xc} + \rho\frac{\partial \epsilon}{\partial \rho}
///                  = 4/3\epsilon_\mathrm{xc}
/// \f]
/// and
/// \f[
///    k_\mathrm{xc} = \frac{\partial^2\epsilon}{\partial \rho^2} = 4/3*1/3 * \epsilon_\mathrm{xc}
/// \f]
int test_slater_exchange(World& world) {

    const double thresh=FunctionDefaults<3>::get_thresh();
    real_function_3d dens=real_factory_3d(world).f(slater2).truncate_on_project();

    // construct the XC operator
    const bool spin_polarized=false;
    const std::string xc_data="LDA_X";
    XCOperator<double,3> xc(world,xc_data,spin_polarized,copy(dens),copy(dens));

    double energy=xc.compute_xc_energy();
    print("xc energy:",energy);

    const real_function_3d potential=xc.make_xc_potential();
    const real_function_3d vphi=potential*dens;
    double energy2=inner(dens,potential);
    energy2=energy2 * 2.0/(4./3.); // fac 2 for closed shell
    print("xc energy via potential:",energy2);
    double ratio=energy/energy2;
    print("ratio ",ratio);
    double err=std::abs(ratio-1.0);
    if (check_err(err,thresh,"dft potential error")) return 1;

    const real_function_3d vphiphi=xc.apply_xc_kernel(dens);
    double energy3=inner(dens,vphiphi);
    energy3=energy3*4.0/(4./3.)/(1./3.);// 2 fac 2 for closed shell
    print("xc energy via kernel:",energy3);
    ratio=energy/energy3;
    print("ratio ",ratio);
    if (check_err(err,thresh,"dft kernel error")) return 1;
//    plot_plane(world,kernel,dens,vphiphi,"kernel");


    return 0;

}

int main(int argc, char** argv) {
    madness::initialize(argc, argv);

    madness::World world(SafeMPI::COMM_WORLD);
    world.gop.fence();
    startup(world,argc,argv);

    FunctionDefaults<3>::set_thresh(1.e-6);
    FunctionDefaults<3>::set_cubic_cell(-10, 10);


    int result=0;

    result+=test_slater_exchange(world);

    if (world.rank()==0) {
        if (result==0) print("\ntests passed\n");
        else print("\ntests failed\n");
    }
    madness::finalize();
    return result;
}
