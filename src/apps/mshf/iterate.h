#ifndef ITERATE_H
#define ITERATE_H

#include "input.h"

// Makes BSH Operators
extern std::vector<poperatorT> BSHoperators(World& world, const real_tensorT& E_q,
                                                   const double Es,
                                                   const double tol,
                                                   const double k_fq);


// Iteration routine
extern void iterate(World& world, const real_functionT& U_q,
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
                           const double k_fq);

#endif
