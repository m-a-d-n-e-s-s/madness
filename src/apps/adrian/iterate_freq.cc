
#include <math.h>

#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include "../chem/SCFOperators.h"
#include "../chem/molecule.h"
#include "NWChem.h"  // For nwchem interface
#include "Plot_VTK.h"
#include "TDDFT.h"
#include "TDHF_Basic_Operators2.h"
#include "adrian/ResponseFunction2.h"
#include "adrian/density.h"
#include "adrian/global_functions.h"
#include "adrian/property.h"
#include "adrian/timer.h"
#include "chem/potentialmanager.h"
#include "chem/projector.h"  // For easy calculation of (1 - \hat{\rho}^0)
#include "madness/mra/funcdefaults.h"

// Iterate Frequency Response
void TDHF::IterateFrequencyResponse(World& world,
                                    response_space& rhs_x,
                                    response_space& rhs_y) {
  // Variables needed to iterate
  int iteration = 0;  // Iteration counter
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

  X_space residuals(world, m, n);
  X_space X(x_response, y_response);

  std::vector<X_vector> Xvector;
  std::vector<X_vector> Xresidual;

  for (size_t b = 0; b < m; b++) {
    Xvector.push_back(X_vector(X, b));
    Xresidual.push_back(X_vector(residuals, b));
  }
  // If DFT, initialize the XCOperator
  XCOperator xc = create_xcoperator(world, Gparams.orbitals, Rparams.xc);

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
  double x_shifts = 0;
  double y_shifts = 0;
  // if less negative orbital energy + frequency is positive or greater than 0
  if ((Gparams.energies[n - 1] + omega_n) >= 0.0) {
    // Calculate minimum shift needed such that \eps + \omega + shift < 0
    print("*** we are shifting just so you know!!!");
    x_shifts = -(.05 + omega_n + Gparams.energies[n - 1]);
  }

  // Construct BSH operators
  std::vector<std::shared_ptr<real_convolution_3d>> bsh_x_operators =
      CreateBSHOperatorPropertyVector(world,
                                      x_shifts,
                                      Gparams.energies,
                                      omega_n,
                                      Rparams.small,
                                      FunctionDefaults<3>::get_thresh());
  std::vector<std::shared_ptr<real_convolution_3d>> bsh_y_operators;

  // Negate omega to make this next set of BSH operators \eps - omega
  if (omega_n != 0.0) {
    omega_n = -omega_n;
    bsh_y_operators =
        CreateBSHOperatorPropertyVector(world,
                                        y_shifts,
                                        Gparams.energies,
                                        omega_n,
                                        Rparams.small,
                                        FunctionDefaults<3>::get_thresh());
    omega_n = -omega_n;
  }
  // create couloumb operator
  real_convolution_3d op =
      CoulombOperator(world, Rparams.small, FunctionDefaults<3>::get_thresh());
  // Two ways single vector or vector vector style
  // here I create the orbital products for elctron interaction terms
  response_space orbital_products(world, n, n);

  for (size_t k = 0; k < n; k++) {
    // important to do orb[i]*all orbs
    orbital_products[k] =
        apply(world, op, mul(world, Gparams.orbitals[k], Gparams.orbitals));
  }
  orbital_products.truncate_rf();
  print("orbital_products norms");
  print(orbital_products.norm2());

  vector_real_function_3d rho_omega;
  // Now to iterate
  while (iteration < Rparams.max_iter and !converged) {
    start_timer(world);
    // Basic output
    if (Rparams.print_level >= 1) {
      if (world.rank() == 0)
        printf("\n   Iteration %d at time %.1fs\n", iteration, wall_time());
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

    if (omega_n == 0) {
      rho_omega =
          transition_density(world, Gparams.orbitals, x_response, x_response);
    } else {
      rho_omega =
          transition_density(world, Gparams.orbitals, x_response, y_response);
    }
    // print level 3
    if (Rparams.print_level == 3) {
      print(
          "x norms in iteration before Iterate XY and after computing "
          "rho_omega "
          ": ",
          iteration);
      print(x_response.norm2());

      print(
          "y norms in iteration before IterateXY and after computing "
          "rho_omega",
          iteration);
      print(y_response.norm2());
    }

    IterateXY(world,
              rho_omega,
              orbital_products,
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
        std::cout << "resX " << iteration << ":" << std::endl;
        print(x_norms);
        std::cout << endl;
        std::cout << "resY " << iteration << ":" << std::endl;
        print(y_norms);
      } else {
        std::cout << "resX " << iteration << ":" << std::endl;
        print(x_norms);
      }
      std::cout << endl;
    }

    // Check convergence
    if (std::max(x_norms.absmax(), y_norms.absmax()) < Rparams.dconv and
        iteration > 0) {
      if (Rparams.print_level >= 1) end_timer(world, "This iteration:");
      if (world.rank() == 0) print("\n   Converged!");
      converged = true;
      break;
    }
    // KAIN solver update
    // Returns next set of components
    // If not kain, save the new components
    // print x norms
    /*
    print("x norms in iteration before kain: ", iteration);
    print(x_response.norm2());

    print("y norms in iteration before kain: ", iteration);
    print(y_response.norm2());
    */
    if (Rparams.kain) {
      if (omega_n == 0) {
        rho_omega =
            transition_density(world, Gparams.orbitals, x_response, x_response);
      } else {
        rho_omega =
            transition_density(world, Gparams.orbitals, x_response, y_response);
      }

      X = X_space(x_response, y_response);
      residuals = X_space(x_differences, y_differences);

      // seperate X_space vectors into individual vectors
      for (size_t b = 0; b < m; b++) {
        Xvector[b] = (X_vector(X, b));
        Xresidual[b] = (X_vector(residuals, b));
        /*
        for (real_function_3d fx : Xvector[b].X[0]) {
          print("norm xvector x before kain ", fx.norm2());
        }
        for (real_function_3d fy : Xvector[b].Y[0]) {
          print("norm xvector y before kain ", fy.norm2());
        }
        for (real_function_3d fx : Xresidual[b].X[0]) {
          print("norm xvector x before kain ", fx.norm2());
        }
        for (real_function_3d fy : Xresidual[b].Y[0]) {
          print("norm xvector y before kain ", fy.norm2());
        }
        */
      }

      // Add y functions to bottom of x functions
      // (for KAIN)

      start_timer(world);
      for (size_t b = 0; b < nkain; b++) {
        X_vector kain_X = kain_x_space[b].update(
            Xvector[b], Xresidual[b], FunctionDefaults<3>::get_thresh(), 3.0);

        /*
        for (real_function_3d fx : kain_X.X[0]) {
          print("norm xvector x after kain ", fx.norm2());
        }
        for (real_function_3d fx : kain_X.Y[0]) {
          print("norm xvector y after kain ", fx.norm2());
        }
        */
        x_response[b].assign(kain_X.X[0].begin(), kain_X.X[0].end());
        y_response[b].assign(kain_X.Y[0].begin(), kain_X.Y[0].end());
      }
      end_timer(world, " KAIN update:");

      /*
            print("x norms in iteration after kain: ", iteration);
            print(x_response.norm2());

            print("y norms in iteration after kain: ", iteration);
            print(y_response.norm2());
            print("x norms in iteration after transfer: ", iteration);
            print(x_response.norm2());

            print("y norms in iteration after transfer: ", iteration);
            print(y_response.norm2());
            */
    }

    // print x norms
    // Apply mask
    /*
    for (int i = 0; i < m; i++) x_response[i] = mask * x_response[i];
    if (omega_n != 0.0) {
      for (int i = 0; i < m; i++) y_response[i] = mask * y_response[i];
    }
    */
    // print x norms
    print("x norms in iteration after mask: ", iteration);
    print(x_response.norm2());

    print("y norms in iteration after mask: ", iteration);
    print(y_response.norm2());
    // Update counter
    iteration += 1;

    // Done with the iteration.. truncate
    x_response.truncate_rf();
    if (omega_n != 0.0) x_response.truncate_rf();
    /*
        print("x norms in iteration after truncation: ", iteration);
        print(x_response.norm2());

        print("y norms in iteration after truncation: ", iteration);
        print(y_response.norm2());
        */
    // Save
    if (Rparams.save) {
      start_timer(world);
      save(world, Rparams.save_file);
      if (Rparams.print_level >= 1) end_timer(world, "Save:");
    }
    // Basic output
    if (Rparams.print_level >= 1) end_timer(world, " This iteration:");
    // plot orbitals
    if (Rparams.plot_all_orbitals) {
      PlotGroundandResponseOrbitals(
          world, iteration, x_response, y_response, Rparams, Gparams);
    }
  }
  X = X_space(world, m, n);
  for (size_t b = 0; b < m; b++) {
    Xvector[b] = (X_vector(world, 0));
    Xresidual[b] = (X_vector(world, 0));
  }
}
