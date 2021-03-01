/*
 *
 *   Small class to hold all the functionality needed in creating/manipulating
 *   the  term in the response equations. (Gamma is the perturbed two
 *   electron piece.)
 *
 */

#ifndef MADNESS_APPS_TDHF_RESPONSEPOTENTIAL_INCLUDE
#define MADNESS_APPS_TDHF_RESPONSEPOTENTIAL_INCLUDE

#include "../chem/SCFOperators.h"
#include "../chem/projector.h"
#include "molresponse/response_functions.h"
#include "molresponse/x_space.h"
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <memory>
#include <vector>

namespace madness {

class ResponsePotential {
  // Member Variables
private:
  World& world; // MPI communicator
  std::vector<real_function_3d>
      ground_orbitals;             // Pointer to ground state orbitals
  real_function_3d ground_rho;     // Pointer to ground state density
  real_function_3d v_nuc;          // Pointer to nuclear potential
  QProjector<double, 3> projector; // Will project out ground state
  const double small;              // Smallest lengthscale for coulomb op.
  const double thresh;             // Truncation threshold for coulomb op.
  real_convolution_3d op;          // Coulomb operator
  response_space potential;        // For storing the ground state
  const bool is_dft;         // Flag for calc. dft potential   (default = false)
  const std::string xc_name; // Name of xc functional          (default = "")
  const bool store_potential; // Flag for storing the potential (default = true)

public:
  // Member Functions
private:
  // Calculates ground state density
  real_function_3d calc_density(std::vector<real_function_3d>& orbs) {
    std::vector<real_function_3d> vsq = square(world, orbs);
    compress(world, vsq);
    real_function_3d rho = real_factory_3d(world);
    rho.compress();
    for (unsigned int i = 0; i < vsq.size(); i++) {
      rho.gaxpy(1.0, vsq[i], 1.0, false);
    }
    world.gop.fence();
    return rho;
  }

  // For TDA, y will be zero
  std::vector<real_function_3d> perturbed_density(response_space& x,
                                                  response_space& y) {
    // Get sizes
    size_t m = x.size();
    size_t n = x[0].size();

    // Return container
    std::vector<real_function_3d> densities =
        zero_functions<double, 3>(world, m);

    // Run over virtual...
    for (size_t i = 0; i < m; i++) {
      // Run over occupied...
      for (size_t j = 0; j < n; j++) {
        densities[i] += ground_orbitals[j] * (x[i][j] + y[i][j]);
      }
    }
    // Done!
    return densities;
  }

  // For DFT, uses the XCOperator to construct perturubed vxc
  std::vector<real_function_3d> create_perturbed_vxc(response_space& x) {
    // Create XCoperator
    XCOperator xc(world, xc_name, false, ground_rho,
                  ground_rho); // Assumes closed shell

    // Need a blank ResponseFunction for y
    response_space y(world, x.size(), x[0].size());

    // Get transition density
    std::vector<real_function_3d> drho = perturbed_density(x, y);

    // Return container
    std::vector<real_function_3d> vxc;

    // Finally create vxc
    for (unsigned int i = 0; i < x.size(); i++) {
      vxc.push_back(xc.apply_xc_kernel(drho[i]));
    }

    return vxc;
  }

  // For DFT, uses the XCOperator to construct ground vxc
  real_function_3d create_ground_vxc() {
    // Create XCoperator
    XCOperator xc(world, xc_name, false, ground_rho,
                  ground_rho); // Assumes closed shell

    // Return its potential
    return xc.make_xc_potential();
  }

public:
  // ResponsePotential constructor
  ResponsePotential(World& world, std::vector<real_function_3d>& ground,
                    real_function_3d& v_nuc, double small, double thresh,
                  size_t   r_states, int g_states, bool is_dft = false,
                    std::string xc = "", bool store = true)
      : world(world), ground_orbitals(ground), ground_rho(calc_density(ground)),
        v_nuc(v_nuc), projector(QProjector<double, 3>(world, ground_orbitals)),
        small(small), thresh(thresh), op(CoulombOperator(world, small, thresh)),
        potential(response_space()), is_dft(is_dft), xc_name(xc),
        store_potential(store) {}

  /// Calculates coulomb like terms in the potential
  /// @param[in]    x        current perturbed x orbitals
  /// @param[inout] gammaK   perturbed potential applied to ground state
  /// orbitals
  /// @param[inout] groundJ  ground state potential applied to perturbed
  /// orbitals
  void coulomb_terms(response_space& x, response_space& gammaK,
                     response_space& groundJ) {
    // Calculate intermediaries
    // If store_potential is true, only calculate it once
    size_t num_states = x.num_states;
    size_t num_orbitals = x.num_orbitals;

    if (!store_potential or potential.num_states == 0) {
      // Clear any old data
      potential.clear();

      // Calculate potential
      for (unsigned int i = 0; i < ground_orbitals.size(); i++) {
        std::vector<real_function_3d> psif =
            mul_sparse(world, ground_orbitals[i], ground_orbitals, thresh);
        truncate(world, psif);
        psif = apply(world, op, psif);
        truncate(world, psif);
        potential.push_back(psif);
      }
    }

    // Calc. gammaK
    // (Will hold either exchange or v_xc)
    gammaK = response_space(world, num_states, num_orbitals);
    if (is_dft) {
      // DFT, need d^2/drho^2 E[rho]
      // Create vxc
      std::vector<real_function_3d> vxc = create_perturbed_vxc(x);

      // Apply vxc
      for (unsigned int i = 0; i < x.num_states; i++) {
        for (unsigned int j = 0; j < x.num_orbitals; j++) {
          // Negative counters the subtraction of of 2*J-K
          // (K needs to be added for DFT)
          gammaK[i][j] = -1.0 * vxc[i] * ground_orbitals[j];
        }
      }
    } else {
      // Hartree-Fock perturbed exchange
      gammaK.compress_rf();
      for (size_t p =0; p < x.num_orbitals; p++) {
        for (size_t k = 0; k < x.num_states; k++) {
          for (unsigned int i = 0; i < x.num_orbitals; i++) {
            real_function_3d t1 = potential[i][p] * x[k][i];
            t1.compress();
            gammaK[k][p].gaxpy(1.0, t1, 1.0, true);
          }
        }
      }
    }

    // Calc. groundJ
    // Scaling by 2 (from spin integration) here
    potential.compress_rf();
    real_function_3d coulomb = 2.0 * potential[0][0];
    coulomb.compress();
    for (unsigned int i = 1; i < x.num_orbitals; i++) {
      coulomb.gaxpy(1.0, potential[i][i], 2.0, true);
    }
    // Add in nuclear potential
    coulomb = coulomb + v_nuc;
    groundJ = x * coulomb;

    // Project out ground states from gammaK
    for (unsigned int i = 0; i < gammaK.num_states; i++) {
      gammaK[i] = projector(gammaK[i]);
    }
    // Project out ground states from groundJ
    for (unsigned int i = 0; i < groundJ.num_states; i++) {
      groundJ[i] = projector(groundJ[i]);
    }
  }

  /// Calculates exchange like terms in the potential
  /// @param[in]    x        current perturbed orbitals
  /// @param[inout] gammaJ   perturbed potential applied to ground state
  /// orbitals
  /// @param[inout] groundK  ground state potential applied to perturbed
  /// orbitals
  void exchange_terms(response_space& x, response_space& gammaJ,
                      response_space& groundK) {
    // Clear inputs
    gammaJ.clear();
    groundK.clear();

    // Going to calculate each transition density
    // and then use it to calculare correct pieces of
    // gammaJ and groundK before moving on
    gammaJ = response_space(world, x.size(), x[0].size());
    groundK = response_space(world, x.size(), x[0].size());
    x.compress_rf();
    gammaJ.compress_rf();
    groundK.compress_rf();

    if (is_dft) {
      // Doing DFT
      // Still need derivative of the coulomb operator
      real_function_3d rho = real_function_3d(world);
      for (size_t k = 0; k < x.num_states; k++) {
        // Get transition density
        rho = dot(world, x[k], ground_orbitals);
        rho = apply(op, rho);

        // Post multiply by ground states
        for (size_t p =0; p<x.num_orbitals; p++) {
          gammaJ[k][p] = rho * ground_orbitals[p];
        }
      }

      // Now do v_xc
      // Negative counters the subtraction of of 2*J-K
      // (K needs to be added for DFT)
      real_function_3d vxc = create_ground_vxc();
      groundK = x * vxc;
      groundK = groundK * -1.0;
    } else {
      // Doing Hartree-Fock
      for (size_t k = 0; k < x.num_states; k++) {
        for (size_t p =0; x.num_orbitals; p++) {
          for (unsigned int i = 0; i < x.num_orbitals; i++) {
            // Get the transition density
            real_function_3d rho = x[k][p] * ground_orbitals[i];

            // Apply coulomb op
            rho = apply(op, rho);

            real_function_3d t1 = rho * ground_orbitals[p];
            t1.compress();
            real_function_3d t2 = rho * ground_orbitals[i];
            t2.compress();

            // Add pieces to gammaJ
            gammaJ[k][p].gaxpy(1.0, t1, 1.0, true);

            // Add pieces to groundK
            groundK[k][p].gaxpy(1.0, t2, 1.0, true);
          }
        }
      }
    }

    // Spin integration for singlet states yields this 2
    gammaJ = gammaJ * 2.0;

    // Project out ground states from gammaK
    for (unsigned int i = 0; i < gammaJ.num_states; i++) {
      gammaJ[i] = projector(gammaJ[i]);
    }
    // Project out ground states from groundJ
    for (unsigned int i = 0; i < groundK.num_states; i++) {
      groundK[i] = projector(groundK[i]);
    }
  }
};

} // End namespace madness

#endif

// Dueces
