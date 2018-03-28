#include "densities.h"
#include "loadbalance.h"
#include "states.h"
#include "potential.h"
#include "energies.h"

// Calculates matrix of kinetic energies
comp_tensorT Kmatrix(World& world, const comp_vecfuncT& psi_qu,
                                   const comp_vecfuncT& psi_qd,
                                   const comp_vecfuncT& dpsi_qu_dx,
                                   const comp_vecfuncT& dpsi_qd_dx,
                                   const comp_vecfuncT& dpsi_qu_dy,
                                   const comp_vecfuncT& dpsi_qd_dy,
                                   const comp_vecfuncT& dpsi_qu_dz,
                                   const comp_vecfuncT& dpsi_qd_dz,
                                   const double prec,
                                   const double k_fq)

{
    double ntol = FunctionDefaults<3>::get_thresh() * prec;
    const double_complex one(1.0,0.0);
    comp_vecfuncT ddpsi_qu, ddpsi_qd;

    reconstruct(world, psi_qu);
    reconstruct(world, psi_qd);
    world.gop.fence();

    complex_derivative_3d Dx(world, 0);
    complex_derivative_3d Dy(world, 1);
    complex_derivative_3d Dz(world, 2);
    world.gop.fence();

    ddpsi_qu = apply(world, Dx, dpsi_qu_dx, false);
    world.gop.fence();

    truncate(world, ddpsi_qu, ntol);

    comp_vecfuncT ddpsi_qu_dy = apply(world, Dy, dpsi_qu_dy, false);
    comp_vecfuncT ddpsi_qu_dz = apply(world, Dz, dpsi_qu_dz, false);
    world.gop.fence();

    gaxpy(world, one, ddpsi_qu, double_complex(1,0), ddpsi_qu_dy, false);
    world.gop.fence();
    gaxpy(world, one, ddpsi_qu, double_complex(1,0), ddpsi_qu_dz, false);
    world.gop.fence();
    truncate(world, ddpsi_qu, ntol);
    world.gop.fence();

    ddpsi_qu_dy.clear(),  ddpsi_qu_dz.clear();
    comp_tensorT K1_q = matrix_inner(world, psi_qu, ddpsi_qu);
    world.gop.fence();

    ddpsi_qd = apply(world, Dx, dpsi_qd_dx, false);
    world.gop.fence();
    truncate(world, ddpsi_qd, ntol);

    comp_vecfuncT ddpsi_qd_dy = apply(world, Dy, dpsi_qd_dy, false);
    comp_vecfuncT ddpsi_qd_dz = apply(world, Dz, dpsi_qd_dz, false);
    world.gop.fence();

    gaxpy(world, one, ddpsi_qd, double_complex(1,0), ddpsi_qd_dy, false);
    world.gop.fence();
    gaxpy(world, one, ddpsi_qd, double_complex(1,0), ddpsi_qd_dz, false);
    world.gop.fence();
    truncate(world, ddpsi_qd, ntol);
    world.gop.fence();

    ddpsi_qd_dy.clear(), ddpsi_qd_dz.clear();
    comp_tensorT K2_q = matrix_inner(world, psi_qd, ddpsi_qd);
    world.gop.fence();

    comp_tensorT K_q  = -k_fq * (K1_q + K2_q);
    world.gop.fence();
    return K_q;
}


// Make a table of kinetic and potential energies for each state as output info
void Energies(World& world, const real_functionT& U_q,
                            const comp_vecfuncT&  psi_qu,
                            const comp_vecfuncT&  psi_qd,
                            const real_functionT& rho,
                            const real_functionT& rho_q,
                            real_tensorT& Etot,
                            const int iter,
                            const comp_tensorT& K_q)
{
    reconstruct(world, psi_qu);
    reconstruct(world, psi_qd);

    int nv = psi_qu.size();

    real_tensorT Ekin(nv);
    for (int i = 0; i < nv; i++) {Ekin[i] += K_q(i,i).real();}
    world.gop.fence();

    if (iter == 0) {
        real_tensorT Epot(nv);
        world.gop.fence();

        for (int i = 0; i < nv; i++) {
            double_complex epot = inner(psi_qu[i], U_q * psi_qu[i])
                                + inner(psi_qd[i], U_q * psi_qd[i]); // need to vectorize
            Epot[i] = epot.real();
        }
        Etot = Ekin + Epot;
    }

    if (world.rank() == 0) {
        std::ostringstream filename1;
        filename1 << "En_0" << std::setw(6) << std::setfill('0') << iter << ".txt";
        FILE * fid1 = fopen(filename1.str().c_str(), "a");

        fprintf(fid1,"\n");

        for (int i = 0; i < nv; i++) {
            fprintf(fid1, "%u,\t%f,\t%f,\t%f\n", i, Ekin[i], Etot[i] - Ekin[i], Etot[i]);
        }

        if (world.rank() == 0) {fprintf(fid1,"\n");}
        if (world.rank() == 0) {fclose(fid1);}
    }
    world.gop.fence();
}


// Calculate total binding energy for output info
double Ebinding(World& world, const real_functionT& rho_p,
                              const real_functionT& rho_n,
                              const real_functionT& tau_p,
                              const real_functionT& tau_n,
                              const real_functionT& lap_p,
                              const real_functionT& lap_n,
                              const real_functionT& dJ_p,
                              const real_functionT& dJ_n,
                              const real_functionT& U_c,
                              const comp_vecfuncT&  psi_pu,
                              const comp_vecfuncT&  psi_pd,
                              const int iter,
                              const int spinorbit,
                              const real_tensorT b,
                              const real_tensorT bp,
                              const double alpha,
                              const double k_fn,
                              const double k_fp)

{
    real_functionT rho = rho_p + rho_n;
    real_functionT U01 =   0.5e0 * b[0] * rho * rho;
    real_functionT U02 = - 0.5e0 * bp[0] * (rho_p * rho_p + rho_n * rho_n);
    double E0 = U01.trace() + U02.trace();

    real_functionT U11 =  b[1] * (rho_p + rho_n) * (tau_p + tau_n);
    real_functionT U12 = -1.e0 * (bp[1] * (rho_p * tau_p) + bp[1] * (rho_n * tau_n));
    double E1 = U11.trace() + U12.trace();

    real_functionT U21 = - 0.5e0 * b[2] * (rho_p + rho_n) * (lap_p + lap_n);
    real_functionT U22 = + 0.5e0 * bp[2] * rho_p * lap_p + 0.5e0 * bp[2] * rho_n * lap_n;
    double E2 = U21.trace() + U22.trace();

    real_functionT power = real_factoryT(world);
    power = 0.0 * power + alpha;
    real_functionT ho2  = binary_op(rho, power, rho_power());

    real_functionT U31 =  (1.e0/3.e0) * b[3] * ho2 * rho * rho;
    real_functionT U32 = -(1.e0/3.e0) * bp[3] * ho2 * (rho_p * rho_p + rho_n * rho_n);
    double E3 = U31.trace() + U32.trace();

    double E4 = 0.0;
    if (spinorbit == 1) {
        real_functionT U41 = -b[4] * rho * (dJ_n + dJ_p);
        real_functionT U42 = -bp[4] * (rho_p * dJ_p + rho_n * dJ_n);
        E4 = U41.trace() + U42.trace();
    }

    double Ec = 0.0;
    world.gop.fence();
    for (unsigned int i = 0; i < psi_pu.size(); i++) {
        double_complex ec = 0.5e0 * inner(U_c * psi_pu[i], psi_pu[i])
                          + 0.5e0 * inner(U_c * psi_pd[i], psi_pd[i]);
        Ec += ec.real();
    }
    world.gop.fence();

    real_functionT U_ex = binary_op(rho_p, rho_p, uex());
    U_ex = 0.75e0 * e2 * rho_p * rho_p * U_ex;
    double Eex = U_ex.trace();

    double Ekin  = k_fp * (tau_p.trace()) + k_fn * (tau_n.trace());
    double Etotal = E0 + E1 + E2 + E3 + E4 + Ec + Eex + Ekin;

    if (iter%10 == 0 || iter == 0) {
        if (world.rank() == 0) {
            std::ostringstream filename1;
            filename1 << "En_0" << std::setw(6) << std::setfill('0') << iter << ".txt";
            FILE * fid1 = fopen(filename1.str().c_str(), "w");
            fprintf(fid1,"\n");
            fprintf(fid1,"%s \t%u\n","Iteration:",iter);
            fprintf(fid1,"%s\n","Energies integrated from density functional: ");
            fprintf(fid1,"\n");
            fprintf(fid1,"%s \t%8.4f,\t%s \t%8.4f,\t%s \t%8.4f\n","E0 part:",E0,"E1 part:",E1,"E2 part:",E2);
            fprintf(fid1,"%s \t%8.4f,\t%s \t%8.4f\n","E3 part:",E3,"E4 part:",E4);
            fprintf(fid1,"%s \t%8.4f,\t%s \t%8.4f,\t%s \t%8.4f\n","Coulomb:",Ec + Eex,"Direct:",Ec,"Exchange:",Eex);
            fprintf(fid1,"%s \t%8.4f\n","Kinetic:",Ekin);
            fprintf(fid1,"%s \t%8.4f\n","Total  :",E0 + E1 + E2 + E3 + E4 + Ec + Eex + Ekin);
            fclose(fid1);
        }
    }
    world.gop.fence();
    return Etotal;
}


