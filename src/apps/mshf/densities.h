#ifndef DENSITIES_H
#define DENSITIES_H

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

// Calculates rho^(alpha+1) for Skyrme potential
struct rho_power
{
    void operator()(const Key<3>& key,
                    real_tensor RHOPOWER,
                    const real_tensor& rho,
                    const real_tensor& power)
    const {
        ITERATOR(RHOPOWER, double p   = rho(IND);
                           double pow = power(IND);
                           if (p <= 0.0) RHOPOWER(IND) = 0.0;
                           else RHOPOWER(IND) = std::pow(p, pow););
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};



// Calculates absolute of a real function
struct abs_func
{
    void operator()(const Key<3>& key, real_tensor ABS,
                    const real_tensor& rho,
                    const real_tensor& rho_u)
    const {
        ITERATOR(ABS, double d = rho(IND);
                      if( d < 0.0 ) d *= -1.0;
                      ABS(IND) = std::abs(d););
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};


// Calculated number density rho
extern real_functionT ndensity(World& world, const comp_vecfuncT& psi_qu,
                                      const comp_vecfuncT& psi_qd);


// Calculates kinetic density tau
extern real_functionT kdensity(World& world, const comp_vecfuncT& dpsi_qu_dx,
                                        const comp_vecfuncT& dpsi_qd_dx,
                                        const comp_vecfuncT& dpsi_qu_dy,
                                        const comp_vecfuncT& dpsi_qd_dy,
                                        const comp_vecfuncT& dpsi_qu_dz,
                                        const comp_vecfuncT& dpsi_qd_dz);

// Calculates spin-orbit density
extern real_functionT so_density(World& world, const comp_vecfuncT& dpsi_qu_dx,
                                        const comp_vecfuncT& dpsi_qd_dx,
                                        const comp_vecfuncT& dpsi_qu_dy,
                                        const comp_vecfuncT& dpsi_qd_dy,
                                        const comp_vecfuncT& dpsi_qu_dz,
                                        const comp_vecfuncT& dpsi_qd_dz,
                                        const double prec);

#endif

