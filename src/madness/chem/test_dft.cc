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

#include <madness.h>
#include<madness/chem/SCFOperators.h>
#include<madness/chem/xcfunctional.h>

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
    err=std::abs(ratio-1.0);
    if (check_err(err,thresh,"dft kernel error")) return 1;
//    plot_plane(world,kernel,dens,vphiphi,"kernel");


    return 0;

}

/// the exact-exchange fraction of a hybrid must not depend on how it was named

/// A hybrid can be requested either through one of the hardcoded aliases or by
/// its libxc name. Both must end up with the same admixture: the fraction is a
/// property of the functional, queried from libxc, not of the input line.
/// Before this was queried, the libxc-name form silently ran with no exact
/// exchange at all.
int test_hybrid_coefficients(World& world) {

    if (world.rank()==0) print("\nentering test_hybrid_coefficients");

    // {input line, expected exact-exchange fraction}
    std::vector<std::pair<std::string,double> > cases = {
            {"LDA",              0.0},
            {"PBE",              0.0},
            {"PBE0",             0.25},
            {"HYB_GGA_XC_PBEH",  0.25},
            {"B3LYP",            0.2},
            {"HYB_GGA_XC_B3LYP", 0.2},
            {"HF",               1.0}
    };

    int result=0;
    for (const auto& c : cases) {
        XCfunctional xcfunc;
        xcfunc.initialize(c.first, false, world);
        const double coeff=xcfunc.hf_exchange_coefficient();
        if (world.rank()==0) printf("  %-20s hf_coeff = %6.4f (expected %6.4f)\n",
                c.first.c_str(), coeff, c.second);
        if (check_err(coeff-c.second, 1.e-10, "hf exchange coefficient of "+c.first)) result=1;
    }
    return result;
}


/// the meta-gga machinery, checked against two exact identities

/// (a) int(tau) must equal the kinetic energy. tau = 1/2 sum_i |grad psi_i|^2, so
///     this pins the factor of 1/2 that libxc has required since version 2.0.0
///     and would catch a tau built from the wrong quantity.
/// (b) the non-multiplicative operator must satisfy its own weak form,
///     <phi|v_tau|psi> = 1/2 int (de/dtau) grad(phi).grad(psi).
///     That is the integration by parts the nested form
///     -1/2 sum_x D_x(v_tau D_x psi) is supposed to realize, and it is also the
///     self-adjointness that make_fock_matrix relies on.
int test_meta_gga(World& world) {

    if (world.rank()==0) print("\nentering test_meta_gga");

    const double thresh=FunctionDefaults<3>::get_thresh();
    int result=0;

    // a normalized s-type orbital, and the closed-shell density that goes with it
    real_function_3d psi=real_factory_3d(world).f(slater);
    psi.scale(1.0/psi.norm2());
    real_function_3d dens=2.0*psi*psi;          // closed shell

    XCOperator<double,3> xc(world,"MGGA_X_TPSS",false,copy(dens),copy(dens));
    MADNESS_CHECK(xc.has_tau_term());
    // per-spin occupation: madness stores one occupied set per spin channel with
    // occupation 1, and tau is stored per spin -- the closed-shell factor 2 is
    // applied by the consumer, as the int(tau) check below does explicitly
    Tensor<double> occ(1l); occ(0l)=1.0;
    xc.set_tau(vecfuncT(1,psi),occ);

    // (a) int(tau) == T. tau is stored per spin, so the total is twice the alpha
    //     value for a closed-shell system.
    const double tau_integral=2.0*xc.get_tau(0).trace();
    Kinetic<double,3> T(world);
    const double kinetic=2.0*T(psi,psi);        // factor 2 for closed shell
    // relative: the test orbital is deliberately cuspy, so its kinetic energy is
    // O(100) and an absolute tolerance would say nothing
    if (world.rank()==0) print("  int(tau)",tau_integral," kinetic",kinetic,
                               " rel err",(tau_integral-kinetic)/kinetic);
    if (check_err((tau_integral-kinetic)/kinetic,thresh*10.0,
                  "int(tau) != kinetic energy")) result=1;

    // (b) the weak form of the operator. make_xc_potential() is what computes vtau.
    const real_function_3d pot=xc.make_xc_potential();
    const real_function_3d vtau=xc.get_vtau();
    MADNESS_CHECK(vtau.is_initialized());

    real_function_3d phi=real_factory_3d(world).f(slater2);
    phi.scale(1.0/phi.norm2());

    for (const real_function_3d& ket : {psi, phi}) {
        const double lhs=inner(phi,xc.apply_tau_term(vecfuncT(1,ket))[0]);
        double rhs=0.0;
        for (int axis=0; axis<3; ++axis) {
            real_derivative_3d D(world,axis);
            rhs+=0.5*inner(vtau*D(phi),D(ket));
        }
        if (world.rank()==0) print("  <phi|v_tau|ket>",lhs," weak form",rhs);
        if (check_err(lhs-rhs,std::max(1.e-8,std::fabs(rhs)*1.e-3),
                      "meta-gga operator does not match its weak form")) result=1;
    }

    // (c) an unoccupied orbital must not contribute to tau. The SCF hands
    //     set_tau() its orbital vectors, which are sized nmo and therefore carry
    //     the virtuals; an unweighted sum of |grad psi|^2 over them inflates tau,
    //     and with it the meta-gga energy and potential, without changing the
    //     density -- so nothing else in this file would notice.
    {
        real_function_3d virt=real_factory_3d(world).f(slater2);
        virt.scale(1.0/virt.norm2());
        XCOperator<double,3> xc2(world,"MGGA_X_TPSS",false,copy(dens),copy(dens));
        Tensor<double> occ2(2l); occ2(0l)=1.0; occ2(1l)=0.0;   // one occupied, one virtual
        xc2.set_tau({psi,virt},occ2);
        // relative, like the int(tau) check: tau is O(100) here, so the two
        // independently truncated operators differ by truncation noise. Including
        // the virtual would add 0.5|grad virt|^2, which is of the same order as
        // tau itself -- orders of magnitude above this tolerance.
        const double dtau=(xc2.get_tau(0)-xc.get_tau(0)).norm2()/xc.get_tau(0).norm2();
        if (world.rank()==0) print("  ||tau(occ+virtual) - tau(occ)|| / ||tau||",dtau);
        if (check_err(dtau,thresh*10.0,"a virtual orbital changed tau")) result=1;
    }

    return result;
}


/// pointwise check of the spin-polarized de/dtau against a finite difference
///
/// Nothing else exercises the polarized meta-gga unpacking: the closed-shell path
/// uses a different libxc entry point (nspin=1), a different index into vxc()'s
/// result vector, and a different stride. This works on raw tensors, so it isolates
/// the functional layer from everything MRA does.
int test_meta_gga_dedtau_polarized(World& world) {
    if (world.rank()==0) print("\nentering test_meta_gga_dedtau_polarized");
    int result=0;

    XCfunctional xcfunc;
    xcfunc.initialize("MGGA_X_TPSS MGGA_C_TPSS",true,world);   // spin polarized
    MADNESS_CHECK(xcfunc.needs_tau());

    // eight sample points, sweeping from nearly spin-compensated to strongly
    // polarized -- the regime an open-shell tail lives in, where the minority
    // channel is orders of magnitude below the majority one
    const long n=2;
    std::vector<double> rhoa={0.30,0.30,0.10,0.10,0.050,0.050,0.020,0.020};
    std::vector<double> rhob={0.25,0.15,0.05,0.01,1.e-3,1.e-4,1.e-5,1.e-6};
    const double za=0.4, zb=0.7;        // zeta = grad(ln rho), one direction

    auto build=[&](double dtaua, double dtaub) {
        std::vector<Tensor<double> > t(XCfunctional::number_xc_args);
        auto mk=[&](int which, const std::vector<double>& v) {
            t[which]=Tensor<double>(n,n,n);
            double* p=t[which].ptr();
            for (long j=0; j<n*n*n; ++j) p[j]=v[j];
        };
        std::vector<double> ta(n*n*n), tb(n*n*n);
        for (long j=0; j<n*n*n; ++j) {
            // sigma_ss = rho_s^2 zeta_s^2, and tau_s >= tau_W = sigma_ss/(8 rho_s);
            // sit safely above the bound so the clamp in make_libxc_args is inactive
            ta[j]=1.5*(rhoa[j]*rhoa[j]*za*za)/(8.0*rhoa[j])+dtaua;
            tb[j]=1.5*(rhob[j]*rhob[j]*zb*zb)/(8.0*rhob[j])+dtaub;
        }
        mk(XCfunctional::enum_rhoa,rhoa);   mk(XCfunctional::enum_rhob,rhob);
        mk(XCfunctional::enum_taua,ta);     mk(XCfunctional::enum_taub,tb);
        // no chi: make_libxc_args contracts the zeta components below
        std::vector<double> zax(n*n*n,za), zbx(n*n*n,zb), zero(n*n*n,0.0);
        mk(XCfunctional::enum_zetaa_x,zax); mk(XCfunctional::enum_zetaa_y,zero);
        mk(XCfunctional::enum_zetaa_z,zero);
        mk(XCfunctional::enum_zetab_x,zbx); mk(XCfunctional::enum_zetab_y,zero);
        mk(XCfunctional::enum_zetab_z,zero);
        return t;
    };

    const std::vector<Tensor<double> > t0=build(0.0,0.0);
    const Tensor<double> va=xcfunc.vxc(t0,0)[7];
    const Tensor<double> vb=xcfunc.vxc(t0,1)[7];

    for (int spin=0; spin<2; ++spin) {
        for (long j=0; j<n*n*n; ++j) {
            // step scaled to tau itself: tau_beta spans several decades here
            const double tau=(spin==0 ? 1.5*rhoa[j]*za*za/8.0 : 1.5*rhob[j]*zb*zb/8.0);
            const double h=1.e-4*std::max(tau,1.e-8);
            const Tensor<double> ep=xcfunc.exc(spin==0 ? build(h,0.0) : build(0.0,h));
            const Tensor<double> em=xcfunc.exc(spin==0 ? build(-h,0.0) : build(0.0,-h));
            const double fd=(ep.ptr()[j]-em.ptr()[j])/(2.0*h);
            const double an=(spin==0 ? va.ptr()[j] : vb.ptr()[j]);
            // de/dtau is screened on the same-spin density (ggatol), so a point
            // below that floor is legitimately zero and carries no information
            const double dens=(spin==0 ? rhoa[j] : rhob[j]);
            const bool screened=(dens<xcfunc.get_ggatol());
            if (world.rank()==0)
                printf("  spin %d  rho %8.1e  de/dtau %13.6e  finite diff %13.6e%s\n",
                       spin,dens,an,fd,screened ? "   (screened)" : "");
            if (screened) continue;
            if (check_err(an-fd,std::max(1.e-8,std::fabs(fd)*1.e-4),
                          "polarized de/dtau does not match the finite difference")) result=1;
        }
    }
    return result;
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
    result+=test_hybrid_coefficients(world);
    result+=test_meta_gga(world);
    result+=test_meta_gga_dedtau_polarized(world);

    if (world.rank()==0) {
        if (result==0) print("\ntests passed\n");
        else print("\ntests failed\n");
    }
    madness::finalize();
    return result;
}
