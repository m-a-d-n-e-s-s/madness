#ifndef ENERGIES_H
#define ENERGIES_H

#include "input.h"

extern comp_tensorT Kmatrix(World& world, const comp_vecfuncT& psi_qu,
                                   const comp_vecfuncT& psi_qd,
                                   const comp_vecfuncT& dpsi_qu_dx,
                                   const comp_vecfuncT& dpsi_qd_dx,
                                   const comp_vecfuncT& dpsi_qu_dy,
                                   const comp_vecfuncT& dpsi_qd_dy,
                                   const comp_vecfuncT& dpsi_qu_dz,
                                   const comp_vecfuncT& dpsi_qd_dz,
                                   const double prec,
                                   const double k_fq);


extern void Energies(World& world, const real_functionT& U_q,
                            const comp_vecfuncT&  psi_qu,
                            const comp_vecfuncT&  psi_qd,
                            const real_functionT& rho,
                            const real_functionT& rho_q,
                            real_tensorT& Etot,
                            const int iter,
                            const comp_tensorT& K_q);



extern double Ebinding(World& world, const real_functionT& rho_p,
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
                              const double k_fp);


#endif
