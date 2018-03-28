#include "orthonorm.h"
#include "potential.h"
#include "energies.h"
#include "states.h"

void orthonorm( World& world, comp_vecfuncT& psi_qu,
                              comp_vecfuncT& psi_qd,
                              comp_vecfuncT& Upsi_qu,
                              comp_vecfuncT& Upsi_qd,
                              real_tensorT& E_q,
                              double& Es,
                              const double tol,
                              const real_functionT& U_q,
                              const comp_vecfuncT& dpsi_qu_dx,
                              const comp_vecfuncT& dpsi_qd_dx,
                              const comp_vecfuncT& dpsi_qu_dy,
                              const comp_vecfuncT& dpsi_qd_dy,
                              const comp_vecfuncT& dpsi_qu_dz,
                              const comp_vecfuncT& dpsi_qd_dz,
                              const int spinorbit,
                              const int meff,
                              const double prec,
                              const double k_fq)
{
    int nst = psi_qu.size();              // number of states

    comp_tensorT C_q;
    comp_tensorT H_q(nst, nst);           // hamiltonian
    comp_tensorT S_q(nst, nst);           // overlap matrix

    //if (world.rank() == 0 && details == 1) {print("             Compute K-matrix ... ");}
    comp_tensorT K_q = Kmatrix(world, psi_qu, psi_qd, dpsi_qu_dx, dpsi_qd_dx, dpsi_qu_dy,
                                      dpsi_qd_dy, dpsi_qu_dz, dpsi_qd_dz, prec, k_fq);
    world.gop.fence();

    //if (world.rank() == 0 && details == 1) {print("             Compute S-matrix and H-matrix... ");}
    S_q =       matrix_inner(world, psi_qu, psi_qu,  true) + matrix_inner(world, psi_qd, psi_qd,  true);
    H_q = K_q + matrix_inner(world, psi_qu, Upsi_qu, true) + matrix_inner(world, psi_qd, Upsi_qd, true);
    world.gop.fence();

    //if (world.rank() == 0 && details == 1) {print("             Compute sygv ... ");}
    sygv(H_q, S_q, 1, C_q, E_q);
    world.gop.fence();

    //if (world.rank() == 0 && details == 1) {print("             transform ... ");}
    if ((spinorbit != 1) && (meff != 1)) {
        E_q = E_q(Slice(0, nst - 1));
        psi_qu = transform(world, psi_qu, C_q(_, Slice(0, nst - 1)), 1.e-3 * tol, false);
        world.gop.fence();
        psi_qd = transform(world, psi_qd, C_q(_, Slice(0, nst - 1)), 1.e-3 * tol, false);
        world.gop.fence();
        Upsi_qu = mul(world, U_q, psi_qu, false);
        world.gop.fence();
        Upsi_qd = mul(world, U_q, psi_qd, false);
    }
    else {
        E_q = E_q(Slice(0, nst - 1));
        psi_qu  = transform(world, psi_qu,  C_q(_, Slice(0, nst - 1)), 1.e-3 * tol, false);
        world.gop.fence();
        psi_qd  = transform(world, psi_qd,  C_q(_, Slice(0, nst - 1)), 1.e-3 * tol, false);
        world.gop.fence();
        Upsi_qu = transform(world, Upsi_qu, C_q(_, Slice(0, nst - 1)), 1.e-3 * tol, false);
        world.gop.fence();
        Upsi_qd = transform(world, Upsi_qd, C_q(_, Slice(0, nst - 1)), 1.e-3 * tol, false);
    }
    world.gop.fence();

    double maxen = 0.0;
    double minen = 0.0;
    for (unsigned int i = 0; i < psi_qu.size(); i++) {
        maxen = std::max(maxen, E_q[i]);
        minen = std::min(minen, E_q[i]);
    }

    Es = std::max(15.0, 2.0 * maxen);
    if (maxen > 0) {if(world.rank() == 0) {print("              Energy Es = ", maxen);}}

    truncate2(world, psi_qu,  psi_qd,  prec);
    truncate2(world, Upsi_qu, Upsi_qd, prec);
    world.gop.fence();
    normalize_ud(world, psi_qu, psi_qd);
    world.gop.fence();

    compress(world, Upsi_qu, false);
    compress(world, Upsi_qd, false);
    compress(world, psi_qu, false);
    compress(world, psi_qd, false);
    world.gop.fence();

    gaxpy(world, 1.0, Upsi_qu, -1.0 * Es, psi_qu, false);
    world.gop.fence();
    gaxpy(world, 1.0, Upsi_qd, -1.0 * Es, psi_qd, false);
    world.gop.fence();

    scale(world, Upsi_qu, -1.0/k_fq);
    scale(world, Upsi_qd, -1.0/k_fq);
    world.gop.fence();
}




