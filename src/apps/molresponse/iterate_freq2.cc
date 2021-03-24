

// Copyright 2021 Adrian Hurtado
#include <math.h>

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include "../chem/SCFOperators.h"
#include "../chem/molecule.h"
#include "NWChem.h"  // For nwchem interface
#include "Plot_VTK.h"
#include "TDDFT.h"
#include "chem/potentialmanager.h"
#include "chem/projector.h"  // For easy calculation of (1 - \hat{\rho}^0)
#include "madness/mra/funcdefaults.h"
#include "molresponse/basic_operators.h"
#include "molresponse/density.h"
#include "molresponse/global_functions.h"
#include "molresponse/property.h"
#include "molresponse/response_functions.h"
#include "molresponse/timer.h"

// Iterate Frequency Response
void TDDFT::IterateFrequencyResponse2(World& world,
                                    X_space Chi
                                    ) {
  // Variables needed to iterate
  size_t iteration = 0;  // Iteration counter
  QProjector<double, 3> projector(
      world, Gparams.orbitals);     // Projector to project out ground state
  size_t n = Gparams.num_orbitals;  // Number of ground state orbitals
  size_t m = Rparams.states;        // Number of excited states
  Tensor<double> x_norms(m);
  // Holds the norms of x function residuals (for convergence)
  Tensor<double> y_norms(m);
  // Holds the norms of y function residuals (for convergence)

  // Holds wave function corrections
  response_space x_differences(world, m, n);
  // Holds wave function corrections
  response_space y_differences(world, m, n);
  response_space x_residuals(world, m, n);
  response_space y_residuals(world, m, n);
  // response functions
  response_space old_x_response(world, m, n);
  response_space old_y_response(world, m, n);
  real_function_3d v_xc;   // For TDDFT
  bool converged = false;  // Converged flag

  // initialize DFT XC functional operator
  XCOperator xc = create_xcoperator(world, Gparams.orbitals, Rparams.xc);

  /*
   * X space refers to X and Y vector spaces |X,Y>
   * X vector is a single |X_b,Y_b> b is a single response state
   * For kain we need a vector of X_vectors
   */

  // create X space residuals
  X_space residuals(world, m, n);
  // Create the X space
  X_space X(x_response, y_response);
  // vector of Xvectors
  std::vector<X_vector> Xvector;
  std::vector<X_vector> Xresidual;
  for (size_t b = 0; b < m; b++) {
    Xvector.push_back(X_vector(X, b));
    Xresidual.push_back(X_vector(residuals, b));
  }
  // If DFT, initialize the XCOperator
  std::vector<XNonlinearSolver<X_vector, double, X_space_allocator>>
      kain_x_space;
  size_t nkain = m;  // (Rparams.omega != 0.0) ? 2 * m : m;
  for (size_t b = 0; b < nkain; b++) {
    kain_x_space.push_back(
        XNonlinearSolver<X_vector, double, X_space_allocator>(
            X_space_allocator(world, n), false));
    if (Rparams.kain) kain_x_space[b].set_maxsub(Rparams.maxsub);
  }

  double omega_n = Rparams.omega;
  omega_n = abs(omega_n);
  omega[0] = omega_n;
  // We compute with positive frequencies
  print("Warning input frequency is assumed to be positive");
  print("Computing at positive frequency omega = ", omega_n);
  double x_shifts{0};
  double y_shifts{0};
  // if less negative orbital energy + frequency is positive or greater than 0
  print("Ground State orbitals");
  print(Gparams.energies);
  if ((Gparams.energies[n - 1] + omega_n) >= 0.0) {
    // Calculate minimum shift needed such that \eps + \omega + shift < 0
    print("*** we are shifting just so you know!!!");
    x_shifts = -(omega_n + Gparams.energies[n - 1]);
  }

  // Construct BSH operators
  std::vector<std::shared_ptr<real_convolution_3d>> bsh_x_operators =
      CreateBSHOperatorPropertyVector(
          world, x_shifts, Gparams.energies, omega_n, .001, 1e-6);
  std::vector<std::shared_ptr<real_convolution_3d>> bsh_y_operators;

  // Negate omega to make this next set of BSH operators \eps - omega
  if (omega_n != 0.0) {
    omega_n = -omega_n;
    bsh_y_operators = CreateBSHOperatorPropertyVector(
        world, y_shifts, Gparams.energies, omega_n, .001, 1e-6);
    omega_n = -omega_n;
  }
  // create couloumb operator
  // Now to iterate
  while (iteration < Rparams.max_iter and !converged) {
    molresponse::start_timer(world);
    // Basic output
    if (Rparams.print_level >= 1) {
      if (world.rank() == 0)
        printf("\n   Iteration %d at time %.1fs\n",
               static_cast<int>(iteration),
               wall_time());
      if (world.rank() == 0) print(" -------------------------------");
    }

    // If omega = 0.0, x = y
    if (Rparams.omega == 0.0) y_response = x_response.copy();
    // Save current to old
    // deep copy of response functions
    old_x_response = x_response.copy();
    old_y_response = y_response.copy();
    if (Rparams.print_level == 3) {
      print("old x norms in iteration after copy  : ", iteration);
      print(old_x_response.norm2());

      print("old y norms in iteration after copy: ", iteration);
      print(old_y_response.norm2());
    }

    // print level 3
    if (Rparams.print_level >= 3) {
      print(
          "x norms in iteration before Iterate XY and after computing "
          "rho_omega "
          ": ",
          iteration,
          " norm : ",
          x_response.norm2());
      print(x_response.norm2());
    }
    IterateXY(world,
              x_response,
              y_response,
              rhs_x,
              rhs_y,
              xc,
              x_shifts,
              Gparams,
              Rparams,
              bsh_x_operators,
              bsh_y_operators,
              ham_no_diag,
              iteration);
    // Get the difference between old and new
    //
    if (Rparams.print_level == 3) {
      print("x norms in iteration after Iterate XY : ", iteration);
      print(x_response.norm2());

      print("y norms in iteration after IterateXY: ", iteration);
      print(y_response.norm2());
    }
    //
    // I need to compute a residual in this new space
    x_differences = old_x_response - x_response;

    if (omega_n != 0.0) y_differences = old_y_response - y_response;

    // Next calculate 2-norm of these vectors of differences
    // Remember: the entire vector is one state
    for (size_t i = 0; i < m; i++) x_norms(i) = norm2(world, x_differences[i]);
    if (omega_n != 0.0) {
      for (size_t i = 0; i < m; i++)
        y_norms(i) = norm2(world, y_differences[i]);
    }

    // Basic output
    if (Rparams.print_level >= 0 and world.rank() == 0) {
      if (omega_n != 0.0) {
        std::cout << "res " << iteration << " X :";
        for (size_t i(0); i < m; i++) {
          std::cout << x_norms[i] << "  ";
        }
        std::cout << " Y :";
        for (size_t i(0); i < m; i++) {
          std::cout << y_norms[i] << "  ";
        }
        std::cout << endl;
      } else {
        print("resX ", iteration, " :", x_norms);
      }
    }

    if (Rparams.kain) {
      /*
      if (omega_n == 0) {
        rho_omega =
            transition_density(world, Gparams.orbitals, x_response, x_response);
      } else {
        rho_omega =
            transition_density(world, Gparams.orbitals, x_response, y_response);
      }
      */

      X = X_space(x_response, y_response);
      residuals = X_space(x_differences, y_differences);

      // seperate X_space vectors into individual vectors
      for (size_t b = 0; b < m; b++) {
        Xvector[b] = (X_vector(X, b));
        Xresidual[b] = (X_vector(residuals, b));
      }

      // Add y functions to bottom of x functions
      // (for KAIN)

      molresponse::start_timer(world);
      for (size_t b = 0; b < nkain; b++) {
        X_vector kain_X = kain_x_space[b].update(
            Xvector[b], Xresidual[b], FunctionDefaults<3>::get_thresh(), 3.0);
        x_response[b].assign(kain_X.X[0].begin(), kain_X.X[0].end());
        y_response[b].assign(kain_X.Y[0].begin(), kain_X.Y[0].end());
      }
      molresponse::end_timer(world, " KAIN update:");
    }
    if (iteration > 0) {
      for (size_t b = 0; b < m; b++) {
        do_step_restriction(
            world, old_x_response[b], x_response[b], "x_response");
        if (omega_n != 0.0) {
          do_step_restriction(
              world, old_y_response[b], y_response[b], "y_response");
        }
      }
    }
    // print x norms
    x_response.truncate_rf();
    if (omega_n == 0.0) y_response = x_response.copy();
    if (omega_n != 0.0) y_response.truncate_rf();

    if (Rparams.print_level >= 1) {
      print("x norms in iteration after truncate: ", iteration);
      print(x_response.norm2());

      print("y norms in iteration after truncate: ", iteration);
      print(y_response.norm2());
    }

    // Check convergence
    if (std::max(x_norms.absmax(), y_norms.absmax()) < Rparams.dconv and
        iteration > 0) {
      if (Rparams.print_level >= 1)
        molresponse::end_timer(world, "This iteration:");
      if (world.rank() == 0) print("\n   Converged!");
      converged = true;
      break;
    }
    // Update counter
    iteration += 1;

    Tensor<double> G(m, m);
    response_space grp(world, m, m);

    for (size_t i(0); i < m; i++) {
      for (size_t j(0); j < m; j++) {
        grp[i][j] =
            dot(world, P[i], x_response[j]) + dot(world, Q[i], y_response[j]);
        G(i, j) = grp[i][j].trace();
        G(i, j) = -2 * G(i, j);
      }
    }
    print("polarizability tensor");
    print(G);
    // Save
    if (Rparams.save) {
      molresponse::start_timer(world);
      save(world, Rparams.save_file);
      if (Rparams.print_level >= 1) molresponse::end_timer(world, "Save:");
    }
    // Basic output
    if (Rparams.print_level >= 1)
      molresponse::end_timer(world, " This iteration:");
    // plot orbitals
    if (Rparams.plot_all_orbitals) {
      PlotGroundandResponseOrbitals(
          world, iteration, x_response, y_response, Rparams, Gparams);
    }
    X = X_space(world, m, n);
    for (size_t b = 0; b < m; b++) {
      Xvector[b] = (X_vector(world, 0));
      Xresidual[b] = (X_vector(world, 0));
    }
  }
}
