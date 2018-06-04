#include "orthonorm.h"
#include "potential.h"
#include "energies.h"
#include "states.h"

void orthonorm( World& world, comp_vecfuncT& psi_qu,
                              comp_vecfuncT& psi_qd,
                              comp_vecfuncT& Upsi_qu,
                              comp_vecfuncT& Upsi_qd,
                              real_tensorT& E_q,
                              double& Es1,
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

    comp_tensorT K_q = Kmatrix(world, psi_qu, psi_qd, dpsi_qu_dx, dpsi_qd_dx, dpsi_qu_dy,
                                      dpsi_qd_dy, dpsi_qu_dz, dpsi_qd_dz, prec, k_fq);
    world.gop.fence();

    S_q =       matrix_inner(world, psi_qu, psi_qu,  true) + matrix_inner(world, psi_qd, psi_qd,  true);
    H_q = K_q + matrix_inner(world, psi_qu, Upsi_qu, true) + matrix_inner(world, psi_qd, Upsi_qd, true);
    world.gop.fence();

    sygv(H_q, S_q, 1, C_q, E_q);
    world.gop.fence();

    comp_tensorT C_q1 = 1.0 * C_q;
    comp_tensorT C_q2 = 1.0 * C_q;
    comp_tensorT C_q3 = 1.0 * C_q;

    if ((spinorbit != 1) && (meff != 1)) {
        E_q = E_q(Slice(0, nst - 1));
        psi_qu = transform(world, psi_qu, C_q(_, Slice(0, nst - 1)), 1.e-3 * tol, false);
        //world.gop.fence();
        psi_qd = transform(world, psi_qd, C_q1(_, Slice(0, nst - 1)), 1.e-3 * tol, false);
        world.gop.fence();
        Upsi_qu = mul(world, U_q, psi_qu, false);
        world.gop.fence();
        Upsi_qd = mul(world, U_q, psi_qd, false);
    }
    else {
        E_q = E_q(Slice(0, nst - 1));
        psi_qu  = transform(world, psi_qu,  C_q(_, Slice(0, nst - 1)), 1.e-3 * tol, false);
        //world.gop.fence();
        psi_qd  = transform(world, psi_qd,  C_q1(_, Slice(0, nst - 1)), 1.e-3 * tol, false);
        //world.gop.fence();
        Upsi_qu = transform(world, Upsi_qu, C_q2(_, Slice(0, nst - 1)), 1.e-3 * tol, false);
        //world.gop.fence();
        Upsi_qd = transform(world, Upsi_qd, C_q3(_, Slice(0, nst - 1)), 1.e-3 * tol, false);
    }
    world.gop.fence();

    double maxen = 0.0;
    double minen = 0.0;
    for (unsigned int i = 0; i < psi_qu.size(); i++) {
        maxen = std::max(maxen, E_q[i]);
        minen = std::min(minen, E_q[i]);
    }

    Es1 = std::max(15.0, 2.0 * maxen);
    double Es2 = Es1;
    if (maxen > 0) {if(world.rank() == 0) {print("              Energy Es = ", Es1);}}

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

    gaxpy(world, 1.0, Upsi_qu, -1.0 * Es1, psi_qu, false);
    //world.gop.fence();
    gaxpy(world, 1.0, Upsi_qd, -1.0 * Es2, psi_qd, false);
    world.gop.fence();

    scale(world, Upsi_qu, -1.0/k_fq);
    scale(world, Upsi_qd, -1.0/k_fq);
    world.gop.fence();
}




