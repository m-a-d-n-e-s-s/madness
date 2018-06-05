#ifndef GROUND_H
#define GROUND_H

#include "input.h"

extern void ground_state(World& world, const int A,
                                       const int Z,
                                       const double Box,
                                       const int initial,
                                       const int boundary,
                                       const int jellium,
                                       const int spinorbit,
                                       const int meff,
                                       const int screening,
                                       const double screenl,
                                       const int avg_pot,
                                       const int avg_lap,
                                       const int avg_wav,
                                       const int lap_comp,
                                       const double prec,
                                       const int vtk_output,
                                       const int txt_output,
                                       const int use_infile,
                                       double knumber,
                                       double thresh,
                                       double chi,
                                       double tol,
                                       double brad,
                                       int IO_nodes,
                                       const real_tensorT b,
                                       const real_tensorT bp,
                                       const double alpha,
                                       const double k_fn,
                                       const double k_fp,
                                       const int project,
                                       const int timing);

#endif
