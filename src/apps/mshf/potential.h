#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "input.h"


// Calculates Fgamma for Coulomb exchange potential WITH screening
struct Fgamma
{
    void operator()(const Key<3>& key, real_tensor FGAMMA,
                    const real_tensor& rho_p,
                    const real_tensor& lambda)
    const {
        ITERATOR(FGAMMA, double length = lambda(IND);
                         double den    = rho_p(IND);
                         double gamma  = length/(std::pow(3.0 * M_PI * M_PI * den, 1.0/3.0));
                         double t1     = atan(2.0/gamma);
                         double ln1    = log(1.0 + 4.0/(gamma*gamma));
                         double g2     = gamma * gamma;
                         if (den < 1.e-10) FGAMMA(IND) = 0.0;
                         else FGAMMA(IND) = 1.0 - (4.0/3.0) * gamma * t1 + 0.5 * g2 * ln1 + (1.0/6.0) * g2 * (1 - 0.25 * g2 * ln1););
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};


// Calculates Coulomb exchange potential via Slater approximation
struct uex
{
    void operator()(const Key<3>& key,
                    real_tensor C_ex,
                    const real_tensor& rho,
                    const real_tensor& rho_u)
    const {
        ITERATOR(C_ex, double d = rho_u(IND);
                       double p = rho(IND);
                       double coeff = 1.0/3.0;
                       double pre = (-1.0) * pow(3.0/M_PI, 1.0/3.0);
                       if (p <= 0.0 || d <= 0.0) C_ex(IND) = 0.0;
                       else C_ex(IND) = pre * p * std::pow(p, coeff - 2.0););
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};


// Calculates effective-mass potential 
extern void Umeff(World& world, comp_vecfuncT& Upsi_qu,
                                comp_vecfuncT& Upsi_qd,
                                const real_functionT& B_q,
                                const comp_vecfuncT& dpsi_qu_dx,
                                const comp_vecfuncT& dpsi_qu_dy,
                                const comp_vecfuncT& dpsi_qu_dz,
                                const comp_vecfuncT& dpsi_qd_dx,
                                const comp_vecfuncT& dpsi_qd_dy,
                                const comp_vecfuncT& dpsi_qd_dz,
                                const double prec);


// Calculates spin-orbit potential
extern void Uso(World& world, comp_vecfuncT& Upsi_qu,
                              comp_vecfuncT& Upsi_qd,
                              const real_functionT& W_q,
                              const comp_vecfuncT& dpsi_qu_dx,
                              const comp_vecfuncT& dpsi_qd_dx,
                              const comp_vecfuncT& dpsi_qu_dy,
                              const comp_vecfuncT& dpsi_qd_dy,
                              const comp_vecfuncT& dpsi_qu_dz,
                              const comp_vecfuncT& dpsi_qd_dz);


// Makes Skyrme Potential
extern void Potential(World& world, const comp_vecfuncT& psi_pu,
                             const comp_vecfuncT& psi_pd, 
                             const comp_vecfuncT& psi_nu,
                             const comp_vecfuncT& psi_nd,
                             real_tensorT& energy_p,
                             real_tensorT& energy_n,
                             real_functionT& U_p,
                             real_functionT& U_n, 
                             real_functionT& rho_p,
                             real_functionT& rho_n,
                             real_functionT& rho, 
                             const int iter,
                             double& BE,
                             const double delta_psi,
                             const double tol,
                             const double brad, 
                             real_functionT& U_p_old,
                             real_functionT& U_n_old, 
                             real_functionT& lap_p_old,
                             real_functionT& lap_n_old,
                             const double L,
                             const int jellium,
                             const int spinorbit,
                             const int screening,
                             const double screenl,
                             const double avg_pot,
                             const double avg_lap,
                             const int lap_comp, 
                             const double prec,
                             const int vtk_output,
                             const int txt_output,
                             const real_tensorT b,
                             const real_tensorT bp,
                             const double alpha,
                             const double k_fn,
                             const double k_fp,
                             const int timing);


#endif
