#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include "input.h"

extern real_functionT laplacian(World& world, real_functionT& rho_q, const double brad, const double prec);

extern real_functionT laplacian1(World& world, const comp_vecfuncT& dpsi_qu_dx, const comp_vecfuncT& dpsi_qd_dx, 
                                               const comp_vecfuncT& dpsi_qu_dy, const comp_vecfuncT& dpsi_qd_dy, 
                                               const comp_vecfuncT& dpsi_qu_dz, const comp_vecfuncT& dpsi_qd_dz, 
                                               const comp_vecfuncT& psi_qu, const comp_vecfuncT& psi_qd);

extern real_functionT laplacian2(World& world, const comp_vecfuncT& dpsi_qu_dx, const comp_vecfuncT& dpsi_qd_dx,
                                               const comp_vecfuncT& dpsi_qu_dy, const comp_vecfuncT& dpsi_qd_dy,
                                               const comp_vecfuncT& dpsi_qu_dz, const comp_vecfuncT& dpsi_qd_dz,
                                               const comp_vecfuncT& psi_qu, const comp_vecfuncT& psi_qd,
                                               const double prec);

#endif

