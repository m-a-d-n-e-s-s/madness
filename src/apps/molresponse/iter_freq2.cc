
// Copyright 2021 Adrian Hurtado
#include <math.h>

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include "../chem/NWChem.h" // For nwchem interface
#include "../chem/SCFOperators.h"
#include "../chem/molecule.h"
#include "Plot_VTK.h"
#include "TDDFT.h"
#include "chem/potentialmanager.h"
#include "chem/projector.h" // For easy calculation of (1 - \hat{\rho}^0)
#include "madness/mra/funcdefaults.h"
#include "molresponse/basic_operators.h"
#include "molresponse/density.h"
#include "molresponse/global_functions.h"
#include "molresponse/property.h"
#include "molresponse/response_functions.h"
#include "molresponse/timer.h"
#include <madness/world/worldmem.h>

// Iterate Frequency Response
void TDDFT::iterate_freq2(World& world) {
  // Variables needed to iterate
  size_t iteration = 0; // Iteration counter
  QProjector<double, 3> projector(world, ground_orbitals);
  size_t n = r_params.num_orbitals(); // Number of ground state orbitals
  size_t m = r_params.n_states();     // Number of excited states

  real_function_3d v_xc;  // For TDDFT
  bool converged = false; // Converged flag
  const double dconv = std::max(FunctionDefaults<3>::get_thresh(), r_params.dconv());

  Tensor<double> bsh_residualsX(m);
  Tensor<double> bsh_residualsY(m);
  vecfuncT rho_omega_old(m);

  // initialize DFT XC functional operator
  XCOperator<double, 3> xc = create_XCOperator(world, ground_orbitals, r_params.xc());

  // create X space residuals
  X_space residuals(world, m, n);
  X_space old_Chi(world, m, n);
  X_space newChi(world, m, n);
  // Create the X space
  // vector of Xvectors
  std::vector<X_vector> Xvector;
  std::vector<X_vector> Xresidual;
  for (size_t b = 0; b < m; b++) {
    Xvector.push_back(X_vector(Chi, b));
    Xresidual.push_back(X_vector(residuals, b));
  }
  // If DFT, initialize the XCOperator<double,3>

  NonLinearXsolver kain_x_space;
  size_t nkain = m; // (r_params.omega() != 0.0) ? 2 * m : m;
  for (size_t b = 0; b < nkain; b++) {
    kain_x_space.push_back(
        XNonlinearSolver<X_vector, double, X_space_allocator>(X_space_allocator(world, n), false));
    if (r_params.kain())
      kain_x_space[b].set_maxsub(r_params.maxsub());
  }
  //
  double omega_n = r_params.omega();
  omega_n = abs(omega_n);
  omega[0] = omega_n;
  // We compute with positive frequencies
  print("Warning input frequency is assumed to be positive");
  print("Computing at positive frequency omega = ", omega_n);
  double x_shifts(0);
  double y_shifts(0);
  // if less negative orbital energy + frequency is positive or greater than 0
  if ((ground_energies[n - 1] + omega_n) >= 0.0) {
    // Calculate minimum shift needed such that \eps + \omega + shift < 0
    print("*** we are shifting just so you know!!!");
    x_shifts = -.05 - (omega_n + ground_energies[n - 1]);
  }
  std::vector<poperatorT> bsh_x_ops = make_bsh_operators_response(world, x_shifts, omega_n);
  std::vector<poperatorT> bsh_y_ops;

  // Negate omega to make this next set of BSH operators \eps - omega
  if (omega_n != 0.0) {
    omega_n = -omega_n;
    std::vector<poperatorT> bsh_x_ops = make_bsh_operators_response(world, y_shifts, omega_n);
    omega_n = -omega_n;
  }
  // create couloumb operator
  // Now to iterate
  // while (iteration < r_params.maxiter() and !converged) {

  for (int iter = 0; iter < r_params.maxiter(); ++iter) {
    // Basic output
    if (r_params.print_level() >= 1) {
      molresponse::start_timer(world);
      if (world.rank() == 0)
        printf("\n   Iteration %d at time %.1fs\n", static_cast<int>(iteration), wall_time());
      if (world.rank() == 0)
        print(" -------------------------------");
    }
    // compute rho_omega
    vecfuncT rho_omega(m);
    molresponse::start_timer(world);
    if (r_params.omega() == 0.0) {
      rho_omega = transition_density(world, ground_orbitals, Chi.X, Chi.X);
    } else {
      rho_omega = transition_density(world, ground_orbitals, Chi.X, Chi.Y);
    }
    molresponse::end_timer(world, "Make density omega");
    print_meminfo(world.rank(), "Make density omega");

    if (iter > 0) {
      if (world.rank() == 0 and (r_params.print_level() > 2))
        print("residuals", bsh_residualsX, bsh_residualsY);
    }

    if (iter < 2 || (iter % 10) == 0) {
      molresponse::start_timer(world);
      loadbal(world, rho_omega, Chi, old_Chi);
      molresponse::end_timer(world, "Load balancing");
      print_meminfo(world.rank(), "Load balancing");
    }

    vector<double> density_residuals;

    if (iter > 0) {
      density_residuals = norm2s(world, (rho_omega - rho_omega_old));
      if (world.rank() == 0 and (r_params.print_level() > 2))
        print("delta rho", density_residuals, "residuals", bsh_residualsX, bsh_residualsY);
    }

    // If omega = 0.0, x = y
    if (r_params.omega() == 0.0)
      Chi.Y = Chi.X.copy();
    old_Chi = Chi.copy();

    update_x_space_response(world,
                            old_Chi,
                            Chi,
                            newChi,
                            xc,
                            bsh_x_ops,
                            bsh_y_ops,
                            projector,
                            x_shifts,
                            omega_n,
                            kain_x_space,
                            Xvector,
                            Xresidual,
                            bsh_residualsX,
                            bsh_residualsY,
                            iteration);

    // apply bsh
    // x_norm is the norm of the res vector
    // Check convergence
    if (std::max(bsh_residualsX.absmax(), bsh_residualsY.absmax()) < dconv * 5.0 and
        iteration > 0) {
      if (r_params.print_level() >= 1)
        molresponse::end_timer(world, "This iteration:");
      if (world.rank() == 0)
        print("\n   Converged!");
      converged = true;
      break;
    }
    // Update counter
    iteration += 1;
    // X_space PQ(P, rhs_y);

    Tensor<double> G = polarizability();
    // Polarizability Tensor
    print("Polarizability Tensor");
    print(G);
    // Save
    if (r_params.save()) {
      molresponse::start_timer(world);
      save(world, r_params.save_file());
      if (r_params.print_level() >= 1)
        molresponse::end_timer(world, "Save:");
    }
    // Basic output
    if (r_params.print_level() >= 1)
      molresponse::end_timer(world, " This iteration:");
    // plot orbitals
    if (r_params.plot_all_orbitals()) {
      PlotGroundandResponseOrbitals(world, iteration, Chi.X, Chi.Y, r_params, g_params);
    }
    /*
for (size_t b = 0; b < m; b++) {
Xvector[b] = (X_vector(world, 0));
Xresidual[b] = (X_vector(world, 0));
}
    */

  } // while converged
}