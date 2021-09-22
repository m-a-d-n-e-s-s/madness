// Copyright 2021 Adrian Hurtado

#ifndef SRC_APPS_MOLRESPONSE_TDDFT_H_
#define SRC_APPS_MOLRESPONSE_TDDFT_H_

#include <chem/SCFOperators.h>
#include <chem/molecule.h>
#include <chem/xcfunctional.h>
#include <madness/constants.h>
#include <madness/mra/mra.h>
#include <madness/mra/nonlinsol.h>  // The kain solver
#include <madness/mra/operator.h>
#include <madness/tensor/distributed_matrix.h>
#include <madness/tensor/solvers.h>
#include <math.h>
#include <molresponse/basic_operators.h>
#include <molresponse/density.h>
#include <molresponse/ground_parameters.h>
#include <molresponse/load_balance.h>
#include <molresponse/property.h>
#include <molresponse/response_functions.h>
#include <molresponse/response_parameters.h>
#include <molresponse/response_potential.h>
#include <molresponse/timer.h>
#include <molresponse/x_space.h>
#include <stdio.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

// Functor from SCF.cc (it wasn't linking right, no idea why, so just copied and
// renamed here) A copy of a MADNESS functor to compute the cartesian moment x^i
// * y^j * z^k (i, j, k integer and >= 0)
class BS_MomentFunctor : public FunctionFunctorInterface<double, 3> {
 private:
  const int i, j, k;

 public:
  BS_MomentFunctor(int i, int j, int k) : i(i), j(j), k(k) {}
  explicit BS_MomentFunctor(const std::vector<int>& x)
      : i(x[0]), j(x[1]), k(x[2]) {}
  double operator()(const Vector<double, 3>& r) const {
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
  GaussianGuess(const coordT& origin,
                const double alpha,
                const std::vector<int> ijk = std::vector<int>(NDIM))
      : origin(origin), exponent(alpha), ijk(ijk) {}

  coordT origin;
  double exponent;       ///< exponent of the guess
  std::vector<int> ijk;  ///< cartesian exponents

  double operator()(const coordT& xyz) const {
    double arg = 0.0, prefac = 1.0;
    for (std::size_t i = 0; i < NDIM; ++i) {
      arg += (xyz[i] - origin[i]) * (xyz[i] - origin[i]);
      prefac *= pow(xyz[i], ijk[i]);
    }
    const double e = exponent * arg;
    return prefac * exp(-e);
  }
};

class ResidualResponseVectors {
 public:
  response_space x;
  response_space y;
  ResidualResponseVectors(World& world, size_t m, int n) {
    x = response_space(world, m, n);
    y = response_space(world, m, n);
  }
};
typedef Vector<double, 3> coordT;
typedef std::shared_ptr<FunctionFunctorInterface<double, 3>> functorT;
typedef Function<double, 3> functionT;
typedef std::vector<functionT> vecfuncT;
typedef std::pair<vecfuncT, vecfuncT> pairvecfuncT;
typedef std::vector<pairvecfuncT> subspaceT;
typedef Tensor<double> tensorT;
typedef DistributedMatrix<double> distmatT;
typedef FunctionFactory<double, 3> factoryT;
typedef SeparatedConvolution<double, 3> operatorT;
typedef std::shared_ptr<operatorT> poperatorT;
typedef Function<std::complex<double>, 3> complex_functionT;
typedef std::vector<complex_functionT> cvecfuncT;
typedef Convolution1D<double_complex> complex_operatorT;
typedef std::vector<XNonlinearSolver<X_vector, double, X_space_allocator>>
    NonLinearXsolver;

class TDDFT {
 public:
  // ResponseParameter object to hold all user input variables
  std::shared_ptr<PotentialManager> potentialmanager;
  // CalculationParameters param;
  functionT mask;

  ResponseParameters r_params;
  GroundParameters g_params;
  Molecule molecule;

  // Tensors for holding energies
  // residuals, and shifts
  Tensor<double> omega;        // Energies of response functions
  Tensor<double> e_residuals;  // Residuals of energies

  // Information that is inferred from input file
  // Ground state orbitals being used in calculation
  std::vector<real_function_3d> act_orbitals;
  std::vector<real_function_3d> ground_orbitals;
  Tensor<double> ground_energies;  // Ground state hamiltonian tensor
  // Ground state energies being used for calculation
  Tensor<double> act_ground_energies;
  Tensor<double> hamiltonian;  // Ground state hamiltonian tensor
  Tensor<double> ham_no_diag;  // Ground state ham. without diagonal (Used when
                               // localized orbitals are given)
  std::vector<int> active;     // The labels of orbitals selected as "active"
  unsigned int act_num_orbitals;  // Number of ground state orbitals being used
                                  // in calculation

  // XCfunction object for DFT calculations
  XCfunctional xcf;

  // Mask function to handle boundary conditions

  // Functions
  real_function_3d stored_v_nuc;   // Stored nuclear potential from ground state
  real_function_3d stored_v_coul;  // Stored coulomb potential from ground state

  X_space Chi;
  X_space PQ;
  density_vector rho;
  functionT rho0;
  vecfuncT rho_omega;

  response_space stored_potential;  // The ground state potential, stored only
                                    // if store_potential is true (default is
                                    // false). Holds the integrals
                                    //   \int dr \frac{\phi_i^\dagger
                                    //   phi_j}{\left| r - r' \right|}
  PropertyBase p;                   // for frequency calculations

  // Get the response Function
  X_space& GetXspace();
  X_space& GetPQspace();
  response_space& GetPVector();
  response_space& GetQVector();
  ResponseParameters GetResponseParameters();
  GroundParameters GetGroundParameters();
  PropertyBase GetPropertyObject();
  // Get Frequencies Omega
  Tensor<double> GetFrequencyOmega();

  poperatorT coulop;
  std::vector<std::shared_ptr<real_derivative_3d>> gradop;
  double vtol;

  // Member variables
 public:
  // Collective constructor for response uses contents of file \c filename and

  TDDFT(World& world, density_vector& rho);
  // Saves a response calculation
  void save(World& world, std::string name);

  // Loads a response calculation
  void load(World& world, std::string name);
  // Initial load balance using vnuc
  void initial_load_bal(World& world);
  void loadbal(World& world, vecfuncT rho_omega, X_space Chi, X_space Chi_old);
  void orbital_load_balance(World& world,
                            vecfuncT& psi0,
                            vecfuncT& psi0_copy,
                            X_space& Chi,
                            X_space& Chi_copy);
  // Normalizes in the response sense
  void normalize(World& world, response_space& f);

  // Normalizes in the response sense
  void normalize(World& world, response_space& f, response_space& g);
  // Normalize X_space xx-yy=1
  void normalize(World& world, X_space& Chi);

  // Prints norms of the given vector
  void print_norms(World& world, response_space function);

  // Returns a set of vector of vector of real_function_3d of proper size,
  // initialized to zero
  response_space response_zero_functions(World& world, size_t m, int n);

  // Returns a list of solid harmonics
  std::map<std::vector<int>, real_function_3d> solid_harmonics(World& world,
                                                               int n);
  // returns a map of real form spherical harmonics  n=level...returns (n+1)^2
  // functions
  std::map<std::vector<int>, real_function_3d> simple_spherical_harmonics(
      World& world,
      int n);
  // returns a map of real form spherical harmonics  n=level...returns (n+1)^2
  // functions
  std::vector<real_function_3d> createDipoleFunctionMap(World& world);
  // Returns initial response functions
  X_space create_trial_functions(World& world,
                                 size_t k,
                                 std::vector<real_function_3d>& orbitals,
                                 size_t print_level);
  X_space create_trial_functions2(World& world,
                                  std::vector<real_function_3d>& orbitals,
                                  size_t print_level);

  response_space PropertyRHS(World& world, PropertyBase& p) const;
  // Returns the derivative of the coulomb operator, applied to ground state
  // orbitals

  // Returns the diagonal (letter A) elements of response matrix
  X_space compute_gamma_full(World& world,
                             X_space& Chi,
                             XCOperator<double, 3> xc);
  X_space compute_gamma_static(World& world,
                               X_space& Chi,
                               XCOperator<double, 3> xc);
  X_space compute_gamma_tda(World& world,
                            X_space& Chi,
                            XCOperator<double, 3> xc);
  // Note: No post multiplication involved here
  real_function_3d Coulomb(World& world);

  // Returns the result of ground state exchange applied to response functions
  response_space exchange(World& world, response_space& f);

  // Returns the ground state potential applied to response functions
  void make_nuclear_potential(World& world);
  X_space compute_V0X(World& world,
                      X_space& Chi,
                      XCOperator<double, 3> xc,
                      bool compute_Y);
  X_space compute_F0X(World& world,
                      X_space& Chi,
                      XCOperator<double, 3> xc,
                      bool compute_Y);

  // Returns a tensor, where entry (i,j) = inner(a[i], b[j]).sum()
  Tensor<double> expectation(World& world,
                             const response_space& a,
                             const response_space& b);
  Tensor<double> expectation2(World& world,
                              const response_space& a,
                              const response_space& b);
  void PrintRFExpectation(World& world,
                          response_space f,
                          response_space g,
                          std::string fname,
                          std::string gname);
  void PrintResponseVectorNorms(World& world,
                                response_space f,
                                std::string fname);
  // Returns the ground state fock operator applied to response functions
  void xy_from_XVector(response_space& x,
                       response_space& y,
                       std::vector<X_vector>& Xvectors);

  void vector_stats(const std::vector<double>& v,
                    double& rms,
                    double& maxabsval) const;

  void vector_stats_new(const Tensor<double> v,
                        double& rms,
                        double& maxabsval) const;

  double do_step_restriction(World& world,
                             const vecfuncT& x,
                             vecfuncT& x_new,
                             std::string spin) const;
  double do_step_restriction(World& world,
                             const vecfuncT& x,
                             vecfuncT& x_new,
                             std::string spin,
                             double maxrotn) const;

  double do_step_restriction(World& world,
                             const vecfuncT& x,
                             const vecfuncT& y,
                             vecfuncT& x_new,
                             vecfuncT& y_new,
                             std::string spin) const;

  X_space Compute_Theta_X(World& world,
                          X_space& Chi,
                          XCOperator<double, 3> xc,
                          std::string calc_type);
  X_space Compute_Lambda_X(World& world,
                           X_space& Chi,
                           XCOperator<double, 3> xc,
                           std::string calc_type);
  // Returns the hamiltonian matrix, equation 45 from the paper
  // -2.0 * (ground_state_energy + excited_state_energy) is positive
  Tensor<double> create_shift(World& world,
                              Tensor<double>& ground,
                              Tensor<double>& omega,
                              size_t print_level,
                              std::string xy);

  // Returns the shift needed for each orbital to make sure
  // (ground_state_energy + excited_state_energy + shift) = target
  Tensor<double> create_shift_target(World& world,
                                     Tensor<double>& ground,
                                     Tensor<double>& omega,
                                     double target,
                                     size_t print_level,
                                     std::string xy);

  // Returns the given shift applied to the given potentials
  response_space apply_shift(World& world,
                             Tensor<double>& shifts,
                             response_space& V,
                             response_space& f);
  // single shift value
  response_space apply_shift(World& world,
                             double& shift,
                             response_space& V,
                             response_space& f);

  // Returns a vector of BSH operators
  std::vector<std::vector<std::shared_ptr<real_convolution_3d>>>
  create_bsh_operators(World& world,
                       Tensor<double>& shift,
                       Tensor<double>& ground,
                       Tensor<double>& omega,
                       double lo,
                       double thresh);

  std::vector<poperatorT> make_bsh_operators_response(
      World& world,
      double& shift,
      double& omega) const;  // Returns a vector of BSH operators
  std::vector<std::vector<std::shared_ptr<real_convolution_3d>>>
  CreateBSHOperatorPropertyVector(World& world,
                                  Tensor<double>& shift,
                                  Tensor<double>& ground,
                                  Tensor<double>& omega,
                                  double lo,
                                  double thresh);
  // here omega and shifts are doubles
  std::vector<std::shared_ptr<real_convolution_3d>>
  CreateBSHOperatorPropertyVector(World& world,
                                  double& shift,
                                  Tensor<double>& ground,
                                  double& omega,
                                  double lo,
                                  double thresh);

  void update_x_space_response(World& world,
                               X_space& old_Chi,
                               X_space& Chi,
                               X_space& residual,
                               XCOperator<double, 3>& xc,
                               std::vector<poperatorT>& bsh_x_ops,
                               std::vector<poperatorT>& bsh_y_ops,
                               QProjector<double, 3>& projector,
                               double& x_shifts,
                               double& omega_n,
                               NonLinearXsolver& kain_x_space,
                               std::vector<X_vector> Xvector,
                               std::vector<X_vector> Xresidual,
                               Tensor<double>& bsh_residualsX,
                               Tensor<double>& bsh_residualsY,
                               size_t iteration,
                               Tensor<double>& maxrotn);

  X_space bsh_update_response(World& world,
                              X_space& theta_X,
                              std::vector<poperatorT>& bsh_x_ops,
                              std::vector<poperatorT>& bsh_y_ops,
                              QProjector<double, 3>& projector,
                              double& x_shifts);

  X_space compute_residual(World& world,
                           X_space& old_Chi,
                           X_space& temp,
                           Tensor<double>& bsh_residualsX,
                           Tensor<double>& bsh_residualsY,
                           std::string calc_type);

  void print_residual_norms(World& world,
                            X_space& old_Chi,
                            bool compute_y,
                            size_t iteration);
  X_space bsh_update_excited(World& world,
                             X_space& theta_X,
                             QProjector<double, 3>& projector,
                             std::vector<bool>& converged);

  void update_x_space_excited(World& world,
                              X_space& old_Chi,
                              X_space& Chi,
                              X_space& old_Lambda_X,
                              X_space& residuals,
                              XCOperator<double, 3>& xc,
                              QProjector<double, 3>& projector,
                              Tensor<double>& omega,
                              NonLinearXsolver& kain_x_space,
                              std::vector<X_vector>& Xvector,
                              std::vector<X_vector>& Xresidual,
                              Tensor<double>& energy_residuals,
                              Tensor<double>& old_energy,
                              Tensor<double>& bsh_residualsX,
                              Tensor<double>& bsh_residualsY,
                              Tensor<double>& S,
                              Tensor<double>& old_S,
                              Tensor<double>& A,
                              Tensor<double>& old_A,
                              std::vector<bool>& converged,
                              size_t iteration,
                              Tensor<double>& maxrotn);
  void compute_new_omegas_transform(World& world,
                                    X_space& old_Chi,
                                    X_space& Chi,
                                    X_space& old_Lambda_X,
                                    X_space& Lambda_X,
                                    Tensor<double>& omega,
                                    Tensor<double>& old_energy,
                                    Tensor<double>& S,
                                    Tensor<double>& old_S,
                                    Tensor<double>& A,
                                    Tensor<double>& old_A,
                                    Tensor<double>& energy_residuals,
                                    size_t iteration);

  X_space compute_residual_excited(World& world,
                                   X_space& old_Chi,
                                   X_space& Chi,
                                   XCOperator<double, 3>& xc,
                                   QProjector<double, 3>& projector,
                                   Tensor<double>& bsh_residualsX,
                                   Tensor<double>& bsh_residualsY,
                                   std::vector<bool>& converged);
  void kain_x_space_update(World& world,
                           X_space& temp,
                           X_space& res,
                           NonLinearXsolver& kain_x_space,
                           std::vector<X_vector>& Xvector,
                           std::vector<X_vector>& Xresidual);
  void x_space_step_restriction(World& world,
                                X_space& old_Chi,
                                X_space& temp,
                                bool restrict_y,
                                Tensor<double>& maxrotn);
  // Returns the second order update to the energy
  Tensor<double> calculate_energy_update(World& world,
                                         response_space& gamma,
                                         response_space& f_residuals,
                                         response_space& new_f,
                                         size_t print_level,
                                         std::string xy);

  // Returns response functions that have been orthonormalized via
  // modified Gram-Schmidt. Note: This is specifically designed for
  // response functions only
  vecfuncT make_density(World& world, X_space& Chi, std::string calc_type);
  response_space gram_schmidt(World& world, response_space& f);

  void gram_schmidt(World& world, response_space& f, response_space& g);

  // Returns the max norm of the given vector of functions
  double calculate_max_residual(World& world, response_space& f);

  // Selects the 'active' orbitals from ground state orbitals to be used in
  // the calculation (based on energy distance from the HOMO.) Function needs
  // knowledge of g_params.orbitals and g_params.ground_energies. Function
  // sets
  void select_active_subspace(World& world);
  // Selects from a list of functions and energies the k functions with the
  // lowest energy
  response_space select_functions(World& world,
                                  response_space& f,
                                  Tensor<double>& energies,
                                  size_t k,
                                  size_t print_level);

  // Calculates the exponentiation of a matrix through first order (I think)
  Tensor<double> matrix_exponential(const Tensor<double>& A);

  // Computes the unitary transformation that diagonalizes the fock matrix
  Tensor<double> get_fock_transformation(World& world,
                                         Tensor<double>& overlap,
                                         Tensor<double>& fock,
                                         Tensor<double>& evals,
                                         const double thresh_degenerate);

  // Sorts the given Tensor of energies
  Tensor<int> sort_eigenvalues(World& world,
                               Tensor<double>& vals,
                               Tensor<double>& vecs);

  // Diagonalizes the fock matrix, taking care of degerate states

  Tensor<double> diagonalizeFockMatrix(World& world,
                                       X_space& Chi,
                                       X_space& Lambda_X,
                                       Tensor<double>& evals,
                                       Tensor<double>& A,
                                       Tensor<double>& S,
                                       const double thresh);
  // Transforms the given matrix of functions according to the given
  // transformation matrix. Used to update orbitals / potentials
  response_space transform(World& world, response_space& f, Tensor<double>& U);

  // If using a larger subspace to diagonalize in, this will put everything in

  void augment(World& world,
               X_space& Chi,
               X_space& old_Chi,
               X_space& Lambda_X,
               X_space& last_Lambda_X,
               Tensor<double>& S_x,
               Tensor<double>& A_x,
               Tensor<double>& old_S,
               Tensor<double>& old_A,
               size_t print_level);
  // If using a larger subspace to diagonalize in, this will put everything in
  // the right spot
  void augment_full(World& world,
                    X_space& Chi,
                    X_space& old_Chi,
                    X_space& Lambda_X,
                    X_space& last_Lambda_X,
                    Tensor<double>& S_x,
                    Tensor<double>& A_x,
                    Tensor<double>& old_S,
                    Tensor<double>& old_A,
                    size_t print_level);

  // If using a larger subspace to diagonalize in, after diagonalization this
  // will put everything in the right spot
  void unaugment(World& world,
                 X_space& Chi,
                 X_space& old_Chi,
                 X_space& Lambda_X,
                 X_space& last_Lambda_X,
                 Tensor<double>& omega,
                 Tensor<double>& S_x,
                 Tensor<double>& A_x,
                 Tensor<double>& old_S,
                 Tensor<double>& old_A,
                 size_t num_states,
                 size_t iter,
                 size_t print_level);
  // If using a larger subspace to diagonalize in, after diagonalization this
  // will put everything in the right spot
  void unaugment_full(World& world,
                      X_space& Chi,
                      X_space& old_Chi,
                      X_space& Lambda_X,
                      X_space& last_Lambda_X,
                      Tensor<double>& omega,
                      Tensor<double>& S_x,
                      Tensor<double>& A_x,
                      Tensor<double>& old_S,
                      Tensor<double>& old_A,
                      size_t num_states,
                      size_t iter,
                      size_t print_level);
  Tensor<double> diagonalizeFullResponseMatrix(World& world,
                                               X_space& Chi,
                                               X_space& Lambda_X,
                                               Tensor<double>& omega,
                                               Tensor<double>& S,
                                               Tensor<double>& A,
                                               const double thresh,
                                               size_t print_level);

  // Similar to what robert did above in "get_fock_transformation"
  Tensor<double> GetFullResponseTransformation(World& world,
                                               Tensor<double>& S,
                                               Tensor<double>& A,
                                               Tensor<double>& evals,
                                               const double thresh);

  // Sorts the given Tensor and response functions in place
  void sort(World& world, Tensor<double>& vals, response_space& f);

  void deflateGuesses(World& world,
                      X_space& Chi,
                      X_space& Lambda_X,
                      Tensor<double>& S,
                      Tensor<double>& omega,
                      size_t& iteration,
                      size_t& m);
  void deflateTDA(World& world,
                  X_space& Chi,
                  X_space& old_Chi,
                  X_space& Lambda_X,
                  X_space& old_Lambda_X,
                  Tensor<double>& S,
                  Tensor<double> old_S,
                  Tensor<double> old_A,
                  Tensor<double>& omega,
                  size_t& iteration,
                  size_t& m);

  void deflateFull(World& world,
                   X_space& Chi,
                   X_space& old_Chi,
                   X_space& Lambda_X,
                   X_space& old_Lambda_X,
                   Tensor<double>& S,
                   Tensor<double> old_S,
                   Tensor<double> old_A,
                   Tensor<double>& omega,
                   size_t& iteration,
                   size_t& m);
  // Creates the XCOperator<double,3>  object and initializes it with correct
  // parameters
  XCOperator<double, 3> create_XCOperator(
      World& world,
      std::vector<real_function_3d> orbitals,
      std::string xc);

  // Uses an XCOperator<double,3>  to construct v_xc for the ground state
  // density Returns d^2/d rho^2 E_xc[rho]
  std::vector<real_function_3d> create_fxc(
      World& world,
      std::vector<real_function_3d>& orbitals,
      response_space& f,
      response_space& g);

  std::vector<real_function_3d> GetWxcOnFDensities(
      World& world,
      const std::vector<real_function_3d>& orbitals,
      const response_space& f);
  std::vector<real_function_3d> GetConjugateWxcOnFDensities(
      World& world,
      const std::vector<real_function_3d>& orbitals,
      const response_space& f);

  std::vector<real_function_3d> CreateXCDerivative(
      World& world,
      const std::vector<real_function_3d>& orbitals,
      const response_space& f);

  // Iterates the trial functions until covergence or it runs out of
  // iterations
  void iterate_excited(World& world, X_space& Chi);

  // Uses member variables
  void analysis(World& world, X_space& Chi);

  // Simplified iterate scheme for guesses
  void iterate_guess(World& world, X_space& guesses);

  // Create and diagonalize the CIS matrix for improved initial guess
  response_space diagonalize_CIS_guess(World& world,
                                       std::vector<real_function_3d>& virtuals,
                                       Tensor<double>& omega,
                                       std::vector<real_function_3d>& orbitals,
                                       Tensor<double>& energies,
                                       double lo,
                                       double thresh,
                                       size_t print_level);

  // Adds random noise to function f
  response_space add_randomness(World& world,
                                response_space& f,
                                double magnitude);

  // Creates the transition density
  functionT make_ground_density(World& world, const vecfuncT& v);
  // Creates the transition density
  std::vector<real_function_3d> transition_density(
      World& world,
      std::vector<real_function_3d>& orbitals,
      response_space& x,
      response_space& y);
  std::vector<real_function_3d> transition_densityTDA(
      World& world,
      std::vector<real_function_3d> const& orbitals,
      response_space& x);
  // Get transition density from f and orbitals
  std::vector<real_function_3d> GetTransitionDensities(
      World& world,
      const std::vector<real_function_3d>& orbitals,
      const response_space& f);

  std::vector<real_function_3d> GetConjugateTransitionDensities(
      World& world,
      const std::vector<real_function_3d>& orbitals,
      const response_space& f);
  // Creates the ground state hamiltonian for the orbitals in the active
  // subspace (aka the orbitals in tda_act_orbitals)
  Tensor<double> CreateGroundHamiltonian(World& world,
                                         std::vector<real_function_3d> f,
                                         size_t print_level);

  // Sets the different k/thresh levels
  template <std::size_t NDIM>
  void set_protocol(World& world, double thresh);

  // Verifies that correct order of polynomial is in use for all
  void check_k(World& world, double thresh, size_t k);

  // Creates random guess functions semi-intelligently(?)
  X_space create_random_guess(World& world,
                              size_t m,
                              size_t n,
                              vector_real_function_3d& grounds,
                              Molecule& molecule);

  // Creates random guess functions semi-intelligently(?)
  std::vector<real_function_3d> create_random_guess(
      World& world,
      size_t m,
      std::vector<real_function_3d>& grounds,
      Molecule& molecule);

  // Creates an initial guess using NWChem outputs from a ground state
  // calculation Requires:
  //    1. nwchem output file (named as "base_name.out")
  //    2. nwchem movecs file (named as "base_name.movecs")
  X_space create_nwchem_guess(World& world, size_t m);

  // Creates potentials using the ResponsePotential object
  // Potentials are modified in place
  void create_all_potentials(World& world,
                             response_space& x,
                             response_space& x_gamma,
                             response_space& x_V0,
                             ResponsePotential& potentials,
                             size_t print_level);

  // Solves the response equations for response states
  void solve_excited_states(World& world);

  // Iterates the response functions until converged or out of iterations

  void iterate_freq2(World& world);
  // Calculates polarizability according to
  // alpha_ij(\omega) = -sum_{m occ} <psi_m(0)|r_i|psi_mj(1)(\omega)> +
  // <psi_mj(1)(-\omega)|r_i|psi_m(0)>
  Tensor<double> polarizability();
  void PrintPolarizabilityAnalysis(World& world,
                                   const Tensor<double> polar_tensor);

  class plotCoords {
   public:
    coord_3d lo, hi;

    plotCoords() {
      lo[0] = 0.0;
      lo[1] = 0.0;
      lo[2] = 0.0;
      hi[0] = 0.0;
      hi[1] = 0.0;
      hi[2] = 0.0;
    }
    plotCoords(size_t direction, double Lp) {
      lo[0] = 0.0;
      lo[1] = 0.0;
      lo[2] = 0.0;
      hi[0] = 0.0;
      hi[1] = 0.0;
      hi[2] = 0.0;
      lo[0] = -Lp;
      hi[0] = Lp;
    }
  };
  plotCoords SetPlotCoord(size_t i, double Lp);

  void PlotGroundandResponseOrbitals(World& world,
                                     size_t iteration,
                                     response_space& x_response,
                                     response_space& y_response,
                                     ResponseParameters const& r_params,
                                     GroundParameters const& g_params);
  void compute_and_print_polarizability(World& world,
                                        X_space& Chi,
                                        X_space& PQ,
                                        std::string message);
  void plot_excited_states(World& world,
                           size_t iteration,
                           response_space& x_response,
                           response_space& y_response,
                           ResponseParameters const& r_params,
                           GroundParameters const& g_params);
  // Solves the response equations for the polarizability
  void solve_response_states(World& world);
};
#endif  // SRC_APPS_MOLRESPONSE_TDDFT_H_

// Deuces
