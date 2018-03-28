#ifndef OUTPUT_H
#define OUTPUT_H

#include "input.h"

// Make output for plotting
extern void output(World& world,  const real_functionT& rho_p,
                           const real_functionT& rho_n,
                           const real_functionT& tau,
                           const real_functionT& lap_p,
                           const real_functionT& lap_n,
                           const real_functionT& U,
                           const real_functionT& Uc,
                           const real_functionT& Uex,
                           const int iter,
                           const double L,
                           const int vtk_output,
                           const int txt_output);


#endif
