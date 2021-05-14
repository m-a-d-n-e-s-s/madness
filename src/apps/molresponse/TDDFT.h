// Copyright 2021 Adrian Hurtado

#ifndef SRC_APPS_MOLRESPONSE_TDDFT_H_
#define SRC_APPS_MOLRESPONSE_TDDFT_H_

#include <chem/molecule.h>
#include <chem/xcfunctional.h>
#include <madness/constants.h>
#include <madness/mra/mra.h>
#include <madness/mra/nonlinsol.h>  // The kain solver
#include <madness/mra/operator.h>
#include <madness/tensor/distributed_matrix.h>
#include <madness/tensor/solvers.h>
#include <math.h>
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

#include "../chem/SCFOperators.h"
#include "../chem/molecule.h"
#include "../chem/xcfunctional.h"
#include "molresponse/basic_operators.h"
#include "molresponse/ground_parameters.h"
#include "molresponse/load_balance.h"
#include "molresponse/property.h"
#include "molresponse/response_functions.h"
#include "molresponse/response_parameters.h"
#include "molresponse/response_potential.h"
#include "molresponse/timer.h"

// Functor from SCF.cc (it wasn't linking right, no idea why, so just copied and
// renamed here) A copy of a MADNESS functor to compute the cartesian moment x^i
// * y^j * z^k (i, j, k integer and >= 0)
class BS_MomentFunctor : public FunctionFunctorInterface<double, 3> {
 private:
  const int i, j, k;

 public:
  BS_MomentFunctor(int i, int j, int k) : i(i), j(j), k(k) {}
  explicit BS_MomentFunctor(const std::vector<int>& x) : i(x[0]), j(x[1]), k(x[2]) {}
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
  GaussianGuess(const coordT& origin, const double alpha, const std::vector<int> ijk = std::vector<int>(NDIM))
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

struct GammaResponseFunctions {
  response_space gamma;
  response_space gamma_conjugate;
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

class TDDFT {
 public:
  // ResponseParameter object to hold all user input variables
  std::shared_ptr<PotentialManager> potentialmanager;
  Molecule molecule;
  // CalculationParameters param;
  XCfunctional xc;
  functionT mask;

  ResponseParameters r_params;

  GroundParameters g_params;

  // Tensors for holding energies
  // residuals, and shifts
  Tensor<double> omega;        // Energies of response functions
  Tensor<double> e_residuals;  // Residuals of energies

  // Information that is inferred from input file
  std::vector<real_function_3d> act_orbitals;  // Ground state orbitals being used in calculation
  Tensor<double> act_ground_energies;          // Ground state energies being used for calculation
  Tensor<double> hamiltonian;                  // Ground state hamiltonian tensor
  Tensor<double> ham_no_diag;                  // Ground state ham. without diagonal (Used when
                                               // localized orbitals are given)
  std::vector<int> active;                     // The labels of orbitals selected as "active"
  unsigned int act_num_orbitals;               // Number of ground state orbitals being used
                                               // in calculation

  // XCfunction object for DFT calculations
  XCfunctional xcf;

  // Mask function to handle boundary conditions

  // Functions
  real_function_3d stored_v_nuc;   // Stored nuclear potential from ground state
  real_function_3d stored_v_coul;  // Stored coulomb potential from ground state

  X_space Chi;
  X_space PQ;

  response_space stored_potential;  // The ground state potential, stored only
                                    // if store_potential is true (default is
                                    // false). Holds the integrals
                                    //   \int dr \frac{\phi_i^\dagger
                                    //   phi_j}{\left| r - r' \right|}
  Property p;                       // for frequency calculations

  // Get the response Function
  X_space& GetXspace();
  X_space& GetPQspace();
  response_space& GetPVector();
  response_space& GetQVector();
  ResponseParameters GetResponseParameters();
  GroundParameters GetGroundParameters();
  Property GetPropertyObject();
  // Get Frequencies Omega
  Tensor<double> GetFrequencyOmega();

  // Member variables
 public:
  // Collective constructor for response uses contents of file \c filename and
  // broadcasts to all nodes
  TDDFT(World& world,                   // MADNESS world object
        const char* filename);          // Input file
  TDDFT(World& world,                   // MADNESS world object
        const std::string& inputfile);  // Input file

  // Collective constructor for Response uses contens of steream \c input and
  // broadcasts to all nodes
  TDDFT(World& world,                          // MADNESS world object
        std::shared_ptr<std::istream> input);  // Pointer to input stream

  TDDFT(World& world, ResponseParameters rparams, GroundParameters gparams);
  // Saves a response calculation
  void save(World& world, std::string name);

  // Loads a response calculation
  void load(World& world, std::string name);

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
  std::map<std::vector<int>, real_function_3d> solid_harmonics(World& world, int n);
  // returns a map of real form spherical harmonics  n=level...returns (n+1)^2
  // functions
  std::map<std::vector<int>, real_function_3d> simple_spherical_harmonics(World& world, int n);
  // returns a map of real form spherical harmonics  n=level...returns (n+1)^2
  // functions
  std::vector<real_function_3d> createDipoleFunctionMap(World& world);
  // Returns initial response functions
  response_space create_trial_functions(World& world,
                                        size_t k,
                                        std::vector<real_function_3d>& orbitals,
                                        size_t print_level);
  response_space create_trial_functions2(World& world, std::vector<real_function_3d>& orbitals, size_t print_level);

  response_space PropertyRHS(World& world, Property& p) const;
  // Returns the derivative of the coulomb operator, applied to ground state
  // orbitals
  response_space CreateCoulombDerivativeRF(World& world,
                                           const response_space& f,                   // response functions
                                           const std::vector<real_function_3d>& phi,  // orbitals
                                           double small,
                                           double thresh);

  response_space CreateCoulombDerivativeRFDagger(World& world,
                                                 const response_space& f,
                                                 const std::vector<real_function_3d>& phi,
                                                 double small,
                                                 double thresh);

  // Returns the derivative of the exchange operator, applied to the ground
  // state orbitals This is the function for TDA only
  response_space CreateExchangeDerivativeRF(World& world,
                                            const response_space& f,
                                            const std::vector<real_function_3d>& phi,
                                            double small,
                                            double thresh);

  response_space CreateExchangeDerivativeRFDagger(World& world,
                                                  const response_space& f,
                                                  const std::vector<real_function_3d>& phi,
                                                  double small,
                                                  double thresh);

  response_space CreateXCDerivativeRF(World& world,
                                      const response_space& f,
                                      const std::vector<real_function_3d>& phi,
                                      double small,
                                      double thresh);
  response_space CreateXCDerivativeRFDagger(World& world,
                                            const response_space& f,
                                            const std::vector<real_function_3d>& phi,
                                            double small,
                                            double thresh);

  // Returns the diagonal (letter A) elements of response matrix
  response_space createAf(World& world,
                          response_space& Vf,
                          response_space& F0_f,
                          response_space& Epsilonf,
                          response_space& Hf,
                          response_space& f,
                          std::vector<real_function_3d>& orbitals,
                          size_t print_level,
                          std::string xy);

  // Returns the off diagonal (letter B) elements of response matrix
  response_space createBf(World& world,
                          response_space& Gf,
                          std::vector<real_function_3d>& orbitals,
                          size_t print_level);

  // Returns gamma (the perturbed 2 electron piece)
  response_space CreateGamma(World& world,
                             response_space& f,
                             response_space& g,
                             std::vector<real_function_3d>& phi,
                             double small,
                             double thresh,
                             size_t print_level,
                             std::string xy);

  response_space ComputeHf(World& world,
                           const response_space& f,
                           const std::vector<real_function_3d>& orbitals,
                           double small,
                           double thresh,
                           size_t print_level,
                           std::string xy);

  response_space ComputeGf(World& world,
                           const response_space& f,
                           const std::vector<real_function_3d>& orbitals,
                           double small,
                           double thresh,
                           size_t print_level,
                           std::string xy);
  GammaResponseFunctions ComputeGammaFunctions(World& world,
                                               response_space& x,
                                               response_space& y,
                                               XCOperator<double, 3> xc,
                                               const GroundParameters& Gparams,
                                               const ResponseParameters& Rparams,
                                               bool compute_Y);
  X_space ComputeGammaFull(World& world, X_space& Chi, XCOperator<double, 3> xc);
  X_space ComputeGammaStatic(World& world, X_space& Chi, XCOperator<double, 3> xc);
  X_space ComputeGammaTDA(World& world, X_space& Chi, XCOperator<double, 3> xc);
  // Note: No post multiplication involved here
  real_function_3d Coulomb(World& world);

  // Returns the result of ground state exchange applied to response functions
  response_space exchange(World& world, response_space& f);

  // Returns the ground state potential applied to response functions
  response_space CreatePotential(World& world,
                                 response_space& f,
                                 XCOperator<double, 3> xc,
                                 size_t print_level,
                                 std::string xy);

  void computeElectronResponse(World& world,
                               ElectronResponseFunctions& I,
                               response_space& x,
                               response_space& y,
                               std::vector<real_function_3d>& orbitals,
                               XCOperator<double, 3> xc,
                               Tensor<double>& hamiltonian,
                               Tensor<double>& ham_no_diag,
                               double small,
                               double thresh,
                               size_t print_level,
                               std::string xy);

  // Returns a tensor, where entry (i,j) = inner(a[i], b[j]).sum()
  Tensor<double> expectation(World& world, const response_space& a, const response_space& b);
  Tensor<double> expectation2(World& world, const response_space& a, const response_space& b);
  void PrintRFExpectation(World& world, response_space f, response_space g, std::string fname, std::string gname);
  void PrintResponseVectorNorms(World& world, response_space f, std::string fname);
  // Returns the ground state fock operator applied to response functions
  response_space CreateFock(World& world, response_space& Vf, response_space& f, size_t print_level, std::string xy);

  Zfunctions ComputeZFunctions(World& world,
                               const std::vector<real_function_3d> rho_omega,
                               response_space orbital_products,
                               response_space& x,
                               response_space& y,
                               XCOperator<double, 3> xc,
                               double x_shifts,
                               double y_shifts,
                               const GroundParameters& Gparams,
                               const ResponseParameters& Rparams,
                               Tensor<double> ham_no_diagonal);
  void xy_from_XVector(response_space& x, response_space& y, std::vector<X_vector>& Xvectors);

  void vector_stats(const std::vector<double>& v, double& rms, double& maxabsval) const;

  double do_step_restriction(World& world, const vecfuncT& x, vecfuncT& x_new, std::string spin) const;

  X_space Compute_Theta_X(World& world, X_space& Chi, XCOperator<double, 3> xc, bool compute_Y);
  X_space Compute_Lambda_X(World& world, X_space& Chi, XCOperator<double, 3> xc, bool compute_Y);
  // Returns the hamiltonian matrix, equation 45 from the paper
  Tensor<double> CreateResponseMatrix(World& world,
                                      response_space& x,
                                      ElectronResponseFunctions& I,
                                      std::vector<real_function_3d>& ground_orbitals,
                                      size_t print_level,
                                      std::string xy);

  // Constructs full response matrix of
  // [ A  B ] [ X ] = w [ X ]
  // [-B -A ] [ Y ]     [ Y ]

  Tensor<double> CreateFullResponseMatrix(World& world,
                                          response_space& x,  // x response functions
                                          response_space& y,  // y response functions
                                          ElectronResponseFunctions& I,
                                          std::vector<real_function_3d>& ground_orbitals,  // ground state orbitals
                                          double small,
                                          double thresh,
                                          size_t print_level);
  // Returns the shift needed for each orbital to make sure
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
  response_space apply_shift(World& world, Tensor<double>& shifts, response_space& V, response_space& f);
  // single shift value
  response_space apply_shift(World& world, double& shift, response_space& V, response_space& f);

  // Returns a vector of BSH operators
  std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> create_bsh_operators(World& world,
                                                                                      Tensor<double>& shift,
                                                                                      Tensor<double>& ground,
                                                                                      Tensor<double>& omega,
                                                                                      double small,
                                                                                      double thresh);

  // Returns a vector of BSH operators
  std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> CreateBSHOperatorPropertyVector(World& world,
                                                                                                 Tensor<double>& shift,
                                                                                                 Tensor<double>& ground,
                                                                                                 Tensor<double>& omega,
                                                                                                 double small,
                                                                                                 double thresh);
  // here omega and shifts are doubles
  std::vector<std::shared_ptr<real_convolution_3d>> CreateBSHOperatorPropertyVector(World& world,
                                                                                    double& shift,
                                                                                    Tensor<double>& ground,
                                                                                    double& omega,
                                                                                    double small,
                                                                                    double thresh);
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
  response_space gram_schmidt(World& world, response_space& f);

  void gram_schmidt(World& world, response_space& f, response_space& g);

  // Returns the max norm of the given vector of functions
  double calculate_max_residual(World& world, response_space& f);

  // Selects the 'active' orbitals from ground state orbitals to be used in
  // the calculation (based on energy distance from the HOMO.) Function needs
  // knowledge of Gparams.orbitals and Gparams.ground_energies. Function sets
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
  Tensor<int> sort_eigenvalues(World& world, Tensor<double>& vals, Tensor<double>& vecs);

  // Diagonalizes the fock matrix, taking care of degerate states
  Tensor<double> diagonalizeFockMatrix(World& world,
                                       Tensor<double>& fock,
                                       response_space& psi,
                                       ElectronResponseFunctions& I,
                                       Tensor<double>& evals,
                                       Tensor<double>& overlap,
                                       const double thresh);

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
  // the right spot
  void augment(World& world,
               Tensor<double>& S_x,
               Tensor<double>& A_x,
               ElectronResponseFunctions& Current,
               ElectronResponseFunctions& Last,
               response_space& x_response,
               Tensor<double>& old_S,
               Tensor<double>& old_A,
               response_space& old_x_resopnse,
               size_t print_level);

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
                 size_t num_states,
                 size_t iter,
                 Tensor<double>& omega,
                 Tensor<double>& S_x,
                 Tensor<double>& A_x,
                 ElectronResponseFunctions& Current,
                 ElectronResponseFunctions& Last,
                 response_space& x_response,
                 Tensor<double>& old_S,
                 Tensor<double>& old_A,
                 response_space& old_x_response,
                 size_t print_level);

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
                      size_t m,
                      size_t iter,
                      Tensor<double>& U,
                      Tensor<double>& omega,
                      Tensor<double>& S,
                      Tensor<double>& A,
                      response_space& x_gamma,
                      response_space& x_response,
                      response_space& V_x_response,
                      response_space& x_fe,
                      response_space& B_x,
                      response_space& y_gamma,
                      response_space& y_response,
                      response_space& V_y_response,
                      response_space& y_fe,
                      response_space& B_y,
                      Tensor<double>& old_S,
                      Tensor<double>& old_A,
                      response_space& old_x_gamma,
                      response_space& old_x_response,
                      response_space& old_V_x_response,
                      response_space& old_x_fe,
                      response_space& old_B_x,
                      response_space& old_y_gamma,
                      response_space& old_y_response,
                      response_space& old_V_y_response,
                      response_space& old_y_fe,
                      response_space& old_B_y,
                      size_t print_level);

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
  // Diagonalize the full response matrix, taking care of degenerate states
  Tensor<double> diagonalizeFullResponseMatrix(World& world,
                                               Tensor<double>& S,
                                               Tensor<double>& A,
                                               response_space& x,
                                               response_space& y,
                                               ElectronResponseFunctions& I,
                                               Tensor<double>& omega,
                                               const double thresh,
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

  void deflateTDA(World& world,
                  Tensor<double>& S,
                  Tensor<double> old_S,
                  Tensor<double> old_A,
                  response_space& x_response,
                  response_space& old_x_response,
                  ElectronResponseFunctions& ElectronResponses,
                  ElectronResponseFunctions& OldElectronResponses,
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
  XCOperator<double, 3> create_XCOperator(World& world, std::vector<real_function_3d> orbitals, std::string xc);

  // Uses an XCOperator<double,3>  to construct v_xc for the ground state density
  // Returns d^2/d rho^2 E_xc[rho]
  std::vector<real_function_3d> create_fxc(World& world,
                                           std::vector<real_function_3d>& orbitals,
                                           response_space& f,
                                           response_space& g);

  std::vector<real_function_3d> GetWxcOnFDensities(World& world,
                                                   const std::vector<real_function_3d>& orbitals,
                                                   const response_space& f);
  std::vector<real_function_3d> GetConjugateWxcOnFDensities(World& world,
                                                            const std::vector<real_function_3d>& orbitals,
                                                            const response_space& f);

  std::vector<real_function_3d> CreateXCDerivative(World& world,
                                                   const std::vector<real_function_3d>& orbitals,
                                                   const response_space& f);

  // Iterates the trial functions until covergence or it runs out of
  // iterations
  void Iterate(World& world);
  void Iterate(World& world, X_space& Chi);

  // Constructs and prints a more detailed analysis of response functions
  // Uses member variables
  void analysis(World& world);
  void analysis(World& world, X_space& Chi);

  // Simplified iterate scheme for guesses
  void IterateGuess(World& world, response_space& guesses);

  // Create and diagonalize the CIS matrix for improved initial guess
  response_space diagonalize_CIS_guess(World& world,
                                       std::vector<real_function_3d>& virtuals,
                                       Tensor<double>& omega,
                                       std::vector<real_function_3d>& orbitals,
                                       Tensor<double>& energies,
                                       double small,
                                       double thresh,
                                       size_t print_level);

  // Adds random noise to function f
  response_space add_randomness(World& world, response_space& f, double magnitude);

  // Creates the transition density
  std::vector<real_function_3d> transition_density(World& world,
                                                   std::vector<real_function_3d> const& orbitals,
                                                   response_space& x,
                                                   response_space& y);
  std::vector<real_function_3d> transition_densityTDA(World& world,
                                                      std::vector<real_function_3d> const& orbitals,
                                                      response_space& x);
  // Get transition density from f and orbitals
  std::vector<real_function_3d> GetTransitionDensities(World& world,
                                                       const std::vector<real_function_3d>& orbitals,
                                                       const response_space& f);

  std::vector<real_function_3d> GetConjugateTransitionDensities(World& world,
                                                                const std::vector<real_function_3d>& orbitals,
                                                                const response_space& f);
  // Creates the ground state hamiltonian for the orbitals in the active
  // subspace (aka the orbitals in tda_act_orbitals)
  Tensor<double> CreateGroundHamiltonian(World& world, std::vector<real_function_3d> f, size_t print_level);

  // Sets the different k/thresh levels
  template <std::size_t NDIM>
  void set_protocol(World& world, double thresh);

  // Verifies that correct order of polynomial is in use for all
  void check_k(World& world, double thresh, size_t k);

  // Creates random guess functions semi-intelligently(?)
  response_space create_random_guess(World& world,
                                     size_t m,
                                     size_t n,
                                     std::vector<real_function_3d>& grounds,
                                     Molecule& molecule);

  // Creates random guess functions semi-intelligently(?)
  std::vector<real_function_3d> create_random_guess(World& world,
                                                    size_t m,
                                                    std::vector<real_function_3d>& grounds,
                                                    Molecule& molecule);

  // Creates an initial guess using NWChem outputs from a ground state
  // calculation Requires:
  //    1. nwchem output file (named as "base_name.out")
  //    2. nwchem movecs file (named as "base_name.movecs")
  response_space create_nwchem_guess(World& world, size_t m);

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
  //
  void IterateFrequencyResponse(World& world, response_space& rhs_x, response_space& rhs_y);

  void iterate_freq_2(World& world);
  // Calculates polarizability according to
  // alpha_ij(\omega) = -sum_{m occ} <psi_m(0)|r_i|psi_mj(1)(\omega)> +
  // <psi_mj(1)(-\omega)|r_i|psi_m(0)>
  void polarizability(World& world, Tensor<double> polar);
  void PrintPolarizabilityAnalysis(World& world, const Tensor<double> polar_tensor, const Tensor<double> omega);

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
                                     ResponseParameters const& Rparams,
                                     GroundParameters const& Gparams);
  void plot_excited_states(World& world,
                           size_t iteration,
                           response_space& x_response,
                           response_space& y_response,
                           ResponseParameters const& Rparams,
                           GroundParameters const& Gparams);
  // Solves the response equations for the polarizability
  void compute_freq_response(World& world, std::string property, X_space& Chi, X_space& PQ);
};
#endif  // SRC_APPS_MOLRESPONSE_TDDFT_H_

// Deuces
