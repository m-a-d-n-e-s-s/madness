#ifndef ORTHONORM_H
#define ORTHONORM_H

#include "input.h"

// Orthonormalizes function-vector computing eigenvectors
extern void orthonorm( World& world, comp_vecfuncT& psi_qu,
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
                              const double k_fq);



#endif

