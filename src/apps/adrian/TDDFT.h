/*!
   \file TDHF.h
   \brief Header file for the TDHF class, which iteratively solves the linear
   response HF equations in the Tamm-Dancoff approximation. \ingroup response
   \addtogroup response

   \par Introduction

   Building on the work presented in the paper from Yanai:
      Yanai, Fann, Beylkin, Harrison; Phys. Chem. Chem. Phys., 2015, 17,
   31405-31416

   Solving equation 37 from Yanai (latex formatting):
   \f$
      \~{x}_p(r) = -2[-\nabla^2 - 2(\epsilon_p^0 + \omega)]^{-1} [\hat{V}^0
   x_p(r) + (1 - \hat{\rho}^0) \Gamma_p(r)] \f$ with \f$ \Gamma_p(r) = \{
   \frac{\partial \hat{g}}{\partial \rho}[\rho^0] \times (\sum_i^{occ} x_i(r)
   \phi_i^\dagger(r'))\} \phi_p(r) \f$

   He lists 12 steps to solve these equations:
      1.  Obtain ground state orbitals {\phi_p} and energies {\epsilon_p}
      2.  Compute a representation of \frac{ \partial^2 E_{xc}}{\partial
   \rho^2}[\rho^0]
      3.  Create guess response functions

    [ 4.  Compute transition density (sum of products of occupied orbitals with
   guess response functions) [ 5.  Obtain \f$\Gamma_p(r)\f$ for current density
    [ 6.  Compute \f$\hat{V}^0 x_p^{(k)}\f$ (does contain HF potential, but rest
   is static as its the ground state values) [ 7.  Obtain initial eigenvalues
   \f$\omega^k\f$ from a matrix diagonalization of [     \f$   A x = S x \omega
   \f$ [     where S is the overlap matrix of the response functions, and A has
   the form [     \f$   A_{ij} = \sum_p \int dr x_p^{(i}}(1 -
   \hat{\rho}^0)[(\hat{F}^0 - \epsilon_p^0) x_p^{(j)}(r) + [ \Gamma_p^{(j)}(r)
   \phi_p(r)] \f$ [     \f$   S_{ij} = \sum_p \int dr x_p^{(i)}(r)
   x_p^{(j)}(r)\f$ [ 8.  Rotate the gamma and potential functions according to
   eigenvectors of the Hamiltonian. [ 9.  Apply BSH integral operator to the
   integral equations (eq. 37)
      10. Repeat steps 4-9 until the residual is within your tolerance
*/

#ifndef MADNESS_APPS_TDHF_H_INCLUDED
#define MADNESS_APPS_TDHF_H_INCLUDED

#include <madness/constants.h>
#include <madness/mra/mra.h>
#include <madness/mra/nonlinsol.h>  // The kain solver
#include <madness/mra/operator.h>
#include <math.h>
#include <stdio.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <random>
#include <vector>

#include "../chem/SCFOperators.h"
#include "../chem/molecule.h"
#include "../chem/xcfunctional.h"
#include "GroundParameters.h"
#include "ResponseFunction2.h"
#include "ResponseParameters.h"
#include "ResponsePotential.h"

using namespace madness;

// Functor from SCF.cc (it wasn't linking right, no idea why, so just copied and
// renamed here) A copy of a MADNESS functor to compute the cartesian moment x^i
// * y^j * z^k (i, j, k integer and >= 0)
class BS_MomentFunctor : public FunctionFunctorInterface<double, 3> {
 private:
  const int i, j, k;

 public:
  BS_MomentFunctor(int i, int j, int k) : i(i), j(j), k(k) {}
  BS_MomentFunctor(const std::vector<int> &x) : i(x[0]), j(x[1]), k(x[2]) {}
  double operator()(const Vector<double, 3> &r) const {
    double xi = 1.0, yj = 1.0, zk = 1.0;
    for (int p = 0; p < i; ++p) xi *= r[0];
    for (int p = 0; p < j; ++p) yj *= r[1];
    for (int p = 0; p < k; ++p) zk *= r[2];
    return xi * yj * zk;
  }
};

/// an N-dimensional real-valued Gaussian function

/// the function looks like
/// \[
/// f(r) = x^i y^j .. z^k exp(-alpha r^2)
/// \]
template <std::size_t NDIM>
class GaussianGuess : public FunctionFunctorInterface<double, NDIM> {
  typedef Vector<double, NDIM> coordT;

 public:
  /// ctor

  /// @param[in]  origin  the origin of the Gauss function
  /// @param[in]  alpha   the exponent exp(-alpha r^2)
  /// @param[in]  ijk     the monomial x^i y^j z^k exp(-alpha r^2) (for NDIM)
  GaussianGuess(const coordT &origin, const double alpha,
                const std::vector<int> ijk = std::vector<int>(NDIM))
      : origin(origin), exponent(alpha), ijk(ijk) {}

  coordT origin;
  double exponent;       ///< exponent of the guess
  std::vector<int> ijk;  ///< cartesian exponents

  double operator()(const coordT &xyz) const {
    double arg = 0.0, prefac = 1.0;
    for (std::size_t i = 0; i < NDIM; ++i) {
      arg += (xyz[i] - origin[i]) * (xyz[i] - origin[i]);
      prefac *= pow(xyz[i], ijk[i]);
    }
    const double e = exponent * arg;
    return prefac * exp(-e);
  }
};

/// Given a molecule and ground state orbitals, solve the response equations
/// in the Tamm-Danchoff approximation.
class TDHF {
 public:
  // ResponseParameter object to hold all user input variables
  ResponseParameters Rparams;

 private:
  // Member variables

  // GroundParameter object to hold all variables needed from
  // ground state calculation. Read from an archive
  GroundParameters Gparams;

  // Tensors for holding energies
  // residuals, and shifts
  Tensor<double> omega;        // Energies of response functions
  Tensor<double> e_residuals;  // Residuals of energies

  // Information that is inferred from input file
  std::vector<real_function_3d>
      act_orbitals;  // Ground state orbitals being used in calculation
  Tensor<double>
      act_ground_energies;  // Ground state energies being used for calculation
  Tensor<double> hamiltonian;  // Ground state hamiltonian tensor
  Tensor<double> ham_no_diag;  // Ground state ham. without diagonal (Used when
                               // localized orbitals are given)
  std::vector<int> active;     // The labels of orbitals selected as "active"
  unsigned int act_num_orbitals;  // Number of ground state orbitals being used
                                  // in calculation

  // XCfunction object for DFT calculations
  XCfunctional xcf;

  // Mask function to handle boundary conditions
  real_function_3d mask;

  // Functions
  real_function_3d stored_v_nuc;   // Stored nuclear potential from ground state
  real_function_3d stored_v_coul;  // Stored coulomb potential from ground state

  ResponseFunction
      x_response;  // Excited states to be solved for.
                   //    Note on storage: The response functions are calculated
                   //    by calculating each transition of occupied --> virtual,
                   //    and thus the actual response function is a sum of of
                   //    all contributions to a specific virtual.

  ResponseFunction
      y_response;  // De-excitation states to be solved for.
                   //    Note on storage: The response functions are calculated
                   //    by calculating each transition of occupied --> virtual,
                   //    and thus the actual response function is a sum of of
                   //    all contributions to a specific virtual.

  ResponseFunction stored_potential;  // The ground state potential, stored only
                                      // if store_potential is true (default is
                                      // false). Holds the integrals
                                      //   \int dr \frac{\phi_i^\dagger
                                      //   phi_j}{\left| r - r' \right|}

 public:
  // Collective constructor for response uses contents of file \c filename and
  // broadcasts to all nodes
  TDHF(World &world,             // MADNESS world object
       const char *input_file);  // Input file

  // Collective constructor for Response uses contens of steream \c input and
  // broadcasts to all nodes
  TDHF(World &world,                          // MADNESS world object
       std::shared_ptr<std::istream> input);  // Pointer to input stream

  // Saves a response calculation
  void save(World &world);

  // Loads a response calculation
  void load(World &world, std::string name);

  // Normalizes in the response sense
  void normalize(World &world, ResponseFunction &f);

  // Normalizes in the response sense
  void normalize(World &world, ResponseFunction &f, ResponseFunction &g);

  // Prints norms of the given vector
  void print_norms(World &world, ResponseFunction function);

  // Prints molecule geometry
  void print_molecule(World &world);

  // Returns a set of vector of vector of real_function_3d of proper size,
  // initialized to zero
  ResponseFunction response_zero_functions(World &world, int m, int n);

  // Returns a list of solid harmonics
  std::map<std::vector<int>, real_function_3d> solid_harmonics(World &world,
                                                               int n);

  // Returns initial response functions
  ResponseFunction create_trial_functions(
      World &world, int k, std::vector<real_function_3d> &orbitals,
      int print_level);

  // Returns dipole operator * molecular orbitals
  ResponseFunction dipole_guess(World &world,
                                std::vector<real_function_3d> orbitals);

  // Returns the derivative of the coulomb operator, applied to ground state
  // orbitals
  ResponseFunction CreateCoulombDerivative(
      World &world,
      ResponseFunction &f,                 // response functions
      std::vector<real_function_3d> &phi,  // orbitals
      double small, double thresh);

  // Returns the derivative of the exchange operator, applied to the ground
  // state orbitals This is the function for TDA only
  ResponseFunction CreateExchangeDerivative(World &world, ResponseFunction &f,
                                            std::vector<real_function_3d> &phi,
                                            double small, double thresh);

  // Returns the diagonal (letter A) elements of response matrix
  ResponseFunction CreateA(World &world, ResponseFunction &fe,
                            ResponseFunction &gamma, ResponseFunction &Vf,
                            ResponseFunction &f,
                            std::vector<real_function_3d> &ground_orbitals,
                            Tensor<double> &hamiltonian, int print_level,
                            std::string xy);

  // Returns the off diagonal (letter B) elements of response matrix
  ResponseFunction CreateB(World &world, ResponseFunction &f,
                            ResponseFunction &g,
                            std::vector<real_function_3d> &orbitals,
                            double small, double thresh, int print_level);

  // Returns gamma (the perturbed 2 electron piece)
  ResponseFunction CreateGamma(World &world, ResponseFunction &f,
                               ResponseFunction &g,
                               std::vector<real_function_3d> &phi, double small,
                               double thresh, int print_level, std::string xy);

  ResponseFunction CreateH(World &world, ResponseFunction &f,
                               std::vector<real_function_3d> &phi, double small,
                               double thresh, int print_level, std::string xy);


  // Returns the coulomb potential of the ground state
  // Note: No post multiplication involved here
  real_function_3d Coulomb(World &world);

  // Returns the result of ground state exchange applied to response functions
  ResponseFunction exchange(World &world, ResponseFunction &f);

  // Returns the ground state potential applied to response functions
  ResponseFunction CreatePotential(World &world, ResponseFunction &f,
                                    XCOperator xc, int print_level,
                                    std::string xy);

  // Returns a tensor, where entry (i,j) = inner(a[i], b[j]).sum()
  Tensor<double> expectation(World &world, ResponseFunction &a,
                             ResponseFunction &b);

  // Returns the ground state fock operator applied to response functions
  ResponseFunction CreateFock(World &world, ResponseFunction &Vf,
                               ResponseFunction &f, int print_level,
                               std::string xy);

  // Returns the hamiltonian matrix, equation 45 from the paper
  Tensor<double> create_response_matrix(
      World &world, ResponseFunction &fe, ResponseFunction &gamma,
      ResponseFunction &V, ResponseFunction &f,
      std::vector<real_function_3d> &ground_orbitals, Tensor<double> &energies,
      int print_level, std::string xy);

  // Constructs full response matrix of
  // [ A  B ] [ X ] = w [ X ]
  // [-B -A ] [ Y ]     [ Y ]
  Tensor<double> CreateFullResponseMatrix(
      World &world, ResponseFunction &x_b, ResponseFunction &Vx,
      ResponseFunction &B_x, ResponseFunction &fe_x, ResponseFunction &x,
      ResponseFunction &y_b, ResponseFunction &Vy, ResponseFunction &B_y,
      ResponseFunction &fe_y, ResponseFunction &y,
      std::vector<real_function_3d> &ground_orbitals,
      Tensor<double> &ground_ham, double small, double thresh, int print_level);

  // Returns the shift needed for each orbital to make sure
  // -2.0 * (ground_state_energy + excited_state_energy) is positive
  Tensor<double> create_shift(World &world, Tensor<double> &ground,
                              Tensor<double> &omega, int print_level,
                              std::string xy);

  // Returns the shift needed for each orbital to make sure
  // (ground_state_energy + excited_state_energy + shift) = target
  Tensor<double> create_shift_target(World &world, Tensor<double> &ground,
                                     Tensor<double> &omega, double target,
                                     int print_level, std::string xy);

  // Returns the given shift applied to the given potentials
  ResponseFunction apply_shift(World &world, Tensor<double> &shifts,
                               ResponseFunction &V, ResponseFunction &f);

  // Returns a vector of BSH operators
  std::vector<std::vector<std::shared_ptr<real_convolution_3d>>>
  create_bsh_operators(World &world, Tensor<double> &shift,
                       Tensor<double> &ground, Tensor<double> &omega,
                       double small, double thresh);

  // Returns the second order update to the energy
  Tensor<double> calculate_energy_update(World &world, ResponseFunction &gamma,
                                         ResponseFunction &f_residuals,
                                         ResponseFunction &new_f,
                                         int print_level, std::string xy);

  // Returns response functions that have been orthonormalized via
  // modified Gram-Schmidt. Note: This is specifically designed for
  // response functions only
  ResponseFunction gram_schmidt(World &world, ResponseFunction &f);

  void gram_schmidt(World &world, ResponseFunction &f, ResponseFunction &g);

  // Returns the max norm of the given vector of functions
  double calculate_max_residual(World &world, ResponseFunction &f);

  // Selects the 'active' orbitals from ground state orbitals to be used in the
  // calculation (based on energy distance from the HOMO.) Function needs
  // knowledge of Gparams.orbitals and Gparams.ground_energies. Function sets
  // act_orbitals and num_act_orbitals.
  void select_active_subspace(World &world);

  // Selects from a list of functions and energies the k functions with the
  // lowest energy
  ResponseFunction select_functions(World &world, ResponseFunction &f,
                                    Tensor<double> &energies, int k,
                                    int print_level);

  // Calculates the exponentiation of a matrix through first order (I think)
  Tensor<double> matrix_exponential(const Tensor<double> &A);

  // Computes the unitary transformation that diagonalizes the fock matrix
  Tensor<double> get_fock_transformation(World &world, Tensor<double> &overlap,
                                         Tensor<double> &fock,
                                         Tensor<double> &evals,
                                         const double thresh_degenerate);

  // Sorts the given Tensor of energies
  Tensor<int> sort_eigenvalues(World &world, Tensor<double> &vals,
                               Tensor<double> &vecs);

  // Diagonalizes the fock matrix, taking care of degerate states
  Tensor<double> diag_fock_matrix(World &world, Tensor<double> &fock,
                                  ResponseFunction &psi, ResponseFunction &Vpsi,
                                  ResponseFunction &gamma, ResponseFunction &fe,
                                  Tensor<double> &evals,
                                  Tensor<double> &overlap, const double thresh);

  // Transforms the given matrix of functions according to the given
  // transformation matrix. Used to update orbitals / potentials
  ResponseFunction transform(World &world, ResponseFunction &f,
                             Tensor<double> &U);

  // If using a larger subspace to diagonalize in, this will put everything in
  // the right spot
  void augment(World &world, Tensor<double> &S_x, Tensor<double> &A_x,
               ResponseFunction &x_gamma, ResponseFunction &x_response,
               ResponseFunction &V_x_response, ResponseFunction &x_fe,
               Tensor<double> &old_S, Tensor<double> &old_A,
               ResponseFunction &old_x_gamma, ResponseFunction &old_x_resopnse,
               ResponseFunction &old_V_x_response, ResponseFunction &old_x_fe,
               int print_level);

  // If using a larger subspace to diagonalize in, this will put everything in
  // the right spot
  void augment_full(
      World &world, Tensor<double> &S, Tensor<double> &A, ResponseFunction &B_x,
      ResponseFunction &x_gamma, ResponseFunction &x_response,
      ResponseFunction &V_x_response, ResponseFunction &x_fe,
      ResponseFunction &B_y, ResponseFunction &y_gamma,
      ResponseFunction &y_response, ResponseFunction &V_y_response,
      ResponseFunction &y_fe, Tensor<double> &old_S, Tensor<double> &old_A,
      ResponseFunction &old_B_x, ResponseFunction &old_x_gamma,
      ResponseFunction &old_x_response, ResponseFunction &old_V_x_response,
      ResponseFunction &old_x_fe, ResponseFunction &old_B_y,
      ResponseFunction &old_y_gamma, ResponseFunction &old_y_response,
      ResponseFunction &old_V_y_response, ResponseFunction &old_y_fe,
      int print_level);

  // If using a larger subspace to diagonalize in, after diagonalization this
  // will put everything in the right spot
  void unaugment(World &world, int m, int iter, Tensor<double> &omega,
                 Tensor<double> &S_x, Tensor<double> &A_x,
                 ResponseFunction &x_gamma, ResponseFunction &x_response,
                 ResponseFunction &V_x_response, ResponseFunction &x_fe,
                 Tensor<double> &old_S, Tensor<double> &old_A,
                 ResponseFunction &old_x_gamma,
                 ResponseFunction &old_x_resopnse,
                 ResponseFunction &old_V_x_response, ResponseFunction &old_x_fe,
                 int print_level);

  // If using a larger subspace to diagonalize in, after diagonalization this
  // will put everything in the right spot
  void unaugment_full(
      World &world, int m, int iter, Tensor<double> &U, Tensor<double> &omega,
      Tensor<double> &S, Tensor<double> &A, ResponseFunction &x_gamma,
      ResponseFunction &x_response, ResponseFunction &V_x_response,
      ResponseFunction &x_fe, ResponseFunction &B_x, ResponseFunction &y_gamma,
      ResponseFunction &y_response, ResponseFunction &V_y_response,
      ResponseFunction &y_fe, ResponseFunction &B_y, Tensor<double> &old_S,
      Tensor<double> &old_A, ResponseFunction &old_x_gamma,
      ResponseFunction &old_x_response, ResponseFunction &old_V_x_response,
      ResponseFunction &old_x_fe, ResponseFunction &old_B_x,
      ResponseFunction &old_y_gamma, ResponseFunction &old_y_response,
      ResponseFunction &old_V_y_response, ResponseFunction &old_y_fe,
      ResponseFunction &old_B_y, int print_level);

  // Diagonalize the full response matrix, taking care of degenerate states
  Tensor<double> DiagonalizeFullResponseMatrix(
      World &world, Tensor<double> &S, Tensor<double> &A, ResponseFunction &x,
      ResponseFunction &Vx, ResponseFunction &x_g, ResponseFunction &x_fe,
      ResponseFunction &B_x, ResponseFunction &y, ResponseFunction &Vy,
      ResponseFunction &y_g, ResponseFunction &y_fe, ResponseFunction &B_y,
      Tensor<double> &omega, const double thresh, int print_level);

  // Similar to what robert did above in "get_fock_transformation"
  Tensor<double> GetFullResponseTransformation(World &world,
                                                  Tensor<double> &S,
                                                  Tensor<double> &A,
                                                  Tensor<double> &evals,
                                                  const double thresh);

  // Sorts the given Tensor and response functions in place
  void sort(World &world, Tensor<double> &vals, ResponseFunction &f);

  // Creates the XCOperator object and initializes it with correct parameters
  XCOperator create_xcoperator(World &world,
                               std::vector<real_function_3d> orbitals,
                               std::string xc);

  // Uses an XCOperator to construct v_xc for the ground state density
  // Returns d^2/d rho^2 E_xc[rho]
  std::vector<real_function_3d> create_fxc(
      World &world, std::vector<real_function_3d> &orbitals,
      ResponseFunction &f, ResponseFunction &g);

  std::vector<real_function_3d> CreateXCDerivative(
      World &world, std::vector<real_function_3d> &orbitals,
      ResponseFunction &f );

  // Iterates the trial functions until covergence or it runs out of iterations
  void Iterate(World &world);

  // Constructs and prints a more detailed analysis of response functions
  // Uses member variables
  void analysis(World &world);

  // Simplified iterate scheme for guesses
  void iterate_guess(World &world, ResponseFunction &guesses);

  // Create and diagonalize the CIS matrix for improved initial guess
  ResponseFunction diagonalize_CIS_guess(
      World &world, std::vector<real_function_3d> &virtuals,
      Tensor<double> &omega, std::vector<real_function_3d> &orbitals,
      Tensor<double> &energies, double small, double thresh, int print_level);

  // Adds random noise to function f
  ResponseFunction add_randomness(World &world, ResponseFunction &f,
                                  double magnitude);

  // Creates the transition density
  std::vector<real_function_3d> transition_density(
      World &world, std::vector<real_function_3d> &orbitals,
      ResponseFunction &x, ResponseFunction &y);

  // Creates the ground state hamiltonian for the orbitals in the active
  // subspace (aka the orbitals in tda_act_orbitals)
  Tensor<double> CreateGroundHamiltonian(World &world,
                                           std::vector<real_function_3d> f,
                                           int print_level);

  // Sets the different k/thresh levels
  template <std::size_t NDIM>
  void set_protocol(World &world, double thresh);

  // Verifies that correct order of polynomial is in use for all
  void check_k(World &world, double thresh, int k);

  // Creates random guess functions semi-intelligently(?)
  ResponseFunction create_random_guess(World &world, int m, int n,
                                       std::vector<real_function_3d> &grounds,
                                       Molecule &molecule);

  // Creates random guess functions semi-intelligently(?)
  std::vector<real_function_3d> create_random_guess(
      World &world, int m, std::vector<real_function_3d> &grounds,
      Molecule &molecule);

  // Creates an initial guess using NWChem outputs from a ground state
  // calculation Requires:
  //    1. nwchem output file (named as "base_name.out")
  //    2. nwchem movecs file (named as "base_name.movecs")
  ResponseFunction create_nwchem_guess(World &world, int m);

  // Creates potentials using the ResponsePotential object
  // Potentials are modified in place
  void create_all_potentials(World &world, ResponseFunction &x,
                             ResponseFunction &x_gamma, ResponseFunction &x_V0,
                             ResponsePotential &potentials, int print_level);

  // Solves the response equations for response states
  void solve(World &world);

  // Iterates the response functions until converged or out of iterations
  void iterate_polarizability(World &world, ResponseFunction &dipoles);

  // Calculates polarizability according to
  // alpha_ij(\omega) = -sum_{m occ} <psi_m(0)|r_i|psi_mj(1)(\omega)> +
  // <psi_mj(1)(-\omega)|r_i|psi_m(0)>
  void polarizability(World &world, Tensor<double> polar);

  // Solves the response equations for the polarizability
  void solve_polarizability(World &world);
};

#endif

// Deuces
