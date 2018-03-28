#include "iterate.h"
#include "potential.h"
#include "energies.h"
#include "orthonorm.h"
#include "states.h"
#include "loadbalance.h"
#include "densities.h"

// Makes BSH Operators
std::vector<poperatorT> BSHoperators(World& world, const real_tensorT& E_q,
                                                   const double Es,
                                                   const double tol,
                                                   const double k_fq)
{
    double ntol = FunctionDefaults<3>::get_thresh() * 0.01;
    int n       = E_q.dim(0);

    std::vector<poperatorT> GBSH(n);
    double lo = ntol * 0.1;

    for (int i = 0; i < n; i++) {
        double E_iq = E_q(i);
        E_iq -= Es;
        if (E_iq > 0.0) {if (world.rank() == 0) {print(E_iq, Es);} E_iq = -0.05;}
        GBSH[i] = poperatorT(BSHOperatorPtr3D(world, sqrt(-(1.0/k_fq) * E_iq), lo, ntol));
    }
    world.gop.fence();
    return GBSH;
}


// Iteration routine
void iterate(World& world, const real_functionT& U_q,
                           comp_vecfuncT& psi_qu,
                           comp_vecfuncT& psi_qd,
                           real_tensorT& E_q,
                           double& delta_psi,
                           const real_functionT& rho,
                           real_functionT& rho_q,
                           const int iter,
                           const double chi,
                           const double tol,
                           const int spinorbit,
                           const int meff,
                           const int avg_wav,
                           const double prec,
                           const real_tensorT b,
                           const real_tensorT bp,
                           const double k_fq)
{
    double ntol = FunctionDefaults<3>::get_thresh() * prec;
    int nt      = psi_qu.size();
    double Es   = 0.0;
    const double_complex one( 1.0,0.0);

    complex_derivative_3d Dx(world, 0);
    complex_derivative_3d Dy(world, 1);
    complex_derivative_3d Dz(world, 2);
    world.gop.fence();

    comp_vecfuncT dpsi_qu_dx = apply(world, Dx, psi_qu, false);
    world.gop.fence();

    truncate(world, dpsi_qu_dx, ntol);
    comp_vecfuncT dpsi_qu_dy = apply(world, Dy, psi_qu, false);
    world.gop.fence();

    truncate(world, dpsi_qu_dy, ntol);
    comp_vecfuncT dpsi_qu_dz = apply(world, Dz, psi_qu, false);
    world.gop.fence();

    truncate(world, dpsi_qu_dz, ntol);

    comp_vecfuncT dpsi_qd_dx = apply(world, Dx, psi_qd, false);
    world.gop.fence();

    truncate(world, dpsi_qd_dx, ntol);
    comp_vecfuncT dpsi_qd_dy = apply(world, Dy, psi_qd, false);
    world.gop.fence();

    truncate(world, dpsi_qd_dy, ntol);
    comp_vecfuncT dpsi_qd_dz = apply(world, Dz, psi_qd, false);
    world.gop.fence();

    truncate(world, dpsi_qd_dz, ntol);

    if (world.rank() == 0) {print("       ortho-normalize: ");}
    world.gop.fence();

    comp_vecfuncT Upsi_qu = mul(world, U_q, psi_qu, false);
    world.gop.fence();
    comp_vecfuncT Upsi_qd = mul(world, U_q, psi_qd, false);
    world.gop.fence();

    if (spinorbit == 1) {
        if (world.rank() == 0) {print("           Calculate Uso ... ");}
        real_functionT W_q = b[4] * rho + bp[4] * rho_q;
        world.gop.fence();
        Uso(world, Upsi_qu, Upsi_qd, W_q, dpsi_qu_dx, dpsi_qd_dx, dpsi_qu_dy, dpsi_qd_dy, dpsi_qu_dz, dpsi_qd_dz);
        world.gop.fence();
        truncate2(world, Upsi_qu, Upsi_qd, prec);
        world.gop.fence();
    }

    if (meff == 1) {
        if (world.rank() == 0) {print("           Calculate Umeff ... ");}
        real_functionT B_q = b[1] * rho - bp[1] * rho_q;
        world.gop.fence();
        Umeff(world, Upsi_qu, Upsi_qd, B_q, dpsi_qu_dx, dpsi_qu_dy, dpsi_qu_dz, dpsi_qd_dx, dpsi_qd_dy, dpsi_qd_dz, prec);
        world.gop.fence();
        truncate2(world, Upsi_qu, Upsi_qd, prec);
        world.gop.fence();
    }

    if (world.rank() == 0) {print("           Do orthonorm ... ");}
    orthonorm(world, psi_qu, psi_qd, Upsi_qu, Upsi_qd, E_q, Es, tol, U_q, dpsi_qu_dx, dpsi_qd_dx, dpsi_qu_dy, dpsi_qd_dy,
                     dpsi_qu_dz, dpsi_qd_dz, spinorbit, meff, prec, k_fq);
    truncate2(world, Upsi_qu, Upsi_qd, prec);
    world.gop.fence();

    world.gop.fence();
    dpsi_qu_dx.clear(); dpsi_qd_dx.clear(); dpsi_qu_dy.clear(); dpsi_qd_dy.clear();
    dpsi_qu_dz.clear(); dpsi_qd_dz.clear();
    world.gop.fence();

    if (world.rank() == 0) {print("       Apply BSH operator ... ");}
    std::vector<poperatorT> ops = BSHoperators(world, E_q, Es, tol, k_fq);
    world.gop.fence();

    loadbalance_v2(world, Upsi_qu, Upsi_qd);
    world.gop.fence();

    //if (world.rank() == 0 && details == 1) {std::cout << "             to Upsi_qu ... " << std::endl;}
    comp_vecfuncT phi_qu = apply(world, ops, Upsi_qu);
    world.gop.fence();
    truncate1(world, phi_qu, prec);
    world.gop.fence();

    //if (world.rank() == 0 && details == 1) {std::cout << "             to Upsi_qd ... " << std::endl;}
    comp_vecfuncT phi_qd = apply(world, ops, Upsi_qd);
    world.gop.fence();
    truncate1(world, phi_qd, prec);
    world.gop.fence();

    Upsi_qu.clear(); ops.clear(); Upsi_qd.clear();
    world.gop.fence();

    //if (world.rank() == 0 && details == 1) {print("       Update states ... ");}
    normalize_ud(world, phi_qu, phi_qd);
    world.gop.fence();

    std::vector<double> rnorm_qu = norm2s(world, sub(world, psi_qu, phi_qu));
    std::vector<double> rnorm_qd = norm2s(world, sub(world, psi_qd, phi_qd));
    std::vector<double> rnorm_q(nt);
    world.gop.fence();

    compress(world, phi_qu, false);
    compress(world, phi_qd, false);
    compress(world, psi_qu, false);
    compress(world, psi_qd, false);
    world.gop.fence();

    if (avg_wav == 1) {
        gaxpy(world, 1.0-chi, psi_qu, chi, phi_qu, false);
        gaxpy(world, 1.0-chi, psi_qd, chi, phi_qd, false);
    }
    else {
        psi_qu = copy(world, phi_qu);
        psi_qd = copy(world, phi_qd);
    }

    world.gop.fence();
    phi_qu.clear(); phi_qd.clear();
    for (int i = 0; i < nt; i++) {rnorm_q[i] = std::sqrt(rnorm_qu[i] * rnorm_qu[i]
                                                       + rnorm_qd[i] * rnorm_qd[i]);}
    world.gop.fence();

    for (int i = 0; i < nt; i++) {delta_psi = std::max(delta_psi, rnorm_q[i]);}

    normalize_ud(world, psi_qu, psi_qd);
    world.gop.fence();
    rho_q = ndensity(world, psi_qu, psi_qd);
    world.gop.fence();

    truncate2(world, psi_qu, psi_qd, prec);
    world.gop.fence();
    normalize_ud(world, psi_qu, psi_qd);
    world.gop.fence();
}




