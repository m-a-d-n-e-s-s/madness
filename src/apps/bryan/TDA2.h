/*
   Iteratively solves the linear response HF equations 

   Building on the work presented in the paper from Yanai:
      Yanai, Fann, Beylkin, Harrison; Phys. Chem. Chem. Phys., 2015, 17, 31405-31416

   Solving equation 37 from Yanai (latex formatting):

      \~{x}_p(r) = -2[-\nabla^2 - 2(\epsilon_p^0 + \omega)]^{-1} [\hat{V}^0 x_p(r) + (1 - \hat{\rho}^0) \Gamma_p(r)]

      with

      \Gamma_p(r) = \{ \frac{\partial \hat{g}}{\partial \rho}[\rho^0] \times (\sum_i^{occ} x_i(r) \phi_i^\dagger(r'))\} \phi_p(r)

 
   He lists 12 steps to solve these equations:
      1.  Obtain ground state orbitals {\phi_p} and energies {\epsilon_p}
      2.  Compute a representation of \frac{ \partial^2 E_{xc}}{\partial \rho^2}[\rho^0]	
      3.  Create guess response functions
    
    [ 4.  Compute transition density (sum of products of occupied orbitals with guess response functions)
    [ 5.  Obtain {\Gamma_p(r)} for current density
    [ 6.  Compute \hat{V}^0 x_p^{(k)} (does contain HF potential, but rest is static as its the ground state values)
    [ 7.  In the first iteration only, obtain initial eigenvalues \omega^k from a matrix diagonalization of
    [        A x = S x \omega
    [     where S is the overlap matrix of the response functions, and A has the form
    [        A_{ij} = \sum_p \int dr x_p^{(i}}(1 - \hat{\rho}^0)[(\hat{F}^0 - \epsilon_p^0) x_p^{(j)}(r) +
    [                 \Gamma_p^{(j)}(r) \phi_p(r)]
    [        S_{ij} = \sum_p \int dr x_p^{(i)}(r) x_p^{(j)}(r)
    [ 8.  Compute the projection (1 - \hat{\rho}^0) \Gamma_p^{(k)}
    [ 9.  Apply BSH integral operator to the integral equations (eq. 37)
    [ 10. Correct trial response functions according to
             \delta \omega^{(k)} = - \frac{ \sum_p\left< \hat{V}^0 x_p^{(k)}(r) + (1 - \hat{\rho}^0) \Gamma_p^{(k)}(r)\right|
                                                \left. x_p^{(k)} - \~{x}_p^{(k)} \right> }
                                          { \sum_p \left| \left| \~{x}_p^{(k)} \right| \right|^2 }
    [ 11. Apply Gram-Schmidt orthonormalization to the response functions

      12. Repeat steps 4-11 until the residual is within your tolerance
*/

#ifndef MADNESS_APPS_TDA_H_INCLUDED
#define MADNESS_APPS_TDA_H_INCLUDED

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/constants.h>
#include <vector>
#include <math.h> 
#include <stdio.h> 
#include <iomanip>
#include <complex>
#include <cmath>
#include "../chem/molecule.h"


using namespace madness;

// Functor from SCF.cc (it wasn't linking right, no idea why, so just copied and renamed here)
/// A copy of a MADNESS functor to compute the cartesian moment x^i * y^j * z^k (i, j, k integer and >= 0)
class BS_MomentFunctor : public FunctionFunctorInterface<double,3> {
private:
    const int i, j, k;
public:
    BS_MomentFunctor(int i, int j, int k) : i(i), j(j), k(k) {}
    BS_MomentFunctor(const std::vector<int>& x) : i(x[0]), j(x[1]), k(x[2]) {}
    double operator()(const Vector<double,3>& r) const {
        double xi=1.0, yj=1.0, zk=1.0;
        for (int p=0; p<i; ++p) xi *= r[0];
        for (int p=0; p<j; ++p) yj *= r[1];
        for (int p=0; p<k; ++p) zk *= r[2];
        return xi*yj*zk;
    }
};

/// an N-dimensional real-valued Gaussian function

/// the function looks like
/// \[
/// f(r) = x^i y^j .. z^k exp(-alpha r^2)
/// \]
template<std::size_t NDIM>
class GaussianGuess : public FunctionFunctorInterface<double,NDIM> {
    typedef Vector<double,NDIM> coordT;

public:

    /// ctor

    /// @param[in]  origin  the origin of the Gauss function
    /// @param[in]  alpha   the exponent exp(-alpha r^2)
    /// @param[in]  ijk     the monomial x^i y^j z^k exp(-alpha r^2) (for NDIM)
    GaussianGuess(const coordT& origin, const double alpha,
            const std::vector<int> ijk=std::vector<int>(NDIM))
            : origin(origin), exponent(alpha), ijk(ijk) {
    }

    coordT origin;
    double exponent;        ///< exponent of the guess
    std::vector<int> ijk;   ///< cartesian exponents

    double operator()(const coordT& xyz) const {
        double arg=0.0, prefac=1.0;
        for (std::size_t i=0; i<NDIM;++i) {
            arg+=(xyz[i]-origin[i])*(xyz[i]-origin[i]);
            prefac*=pow(xyz[i],ijk[i]);
        }
        const double e=exponent*arg;
        return prefac*exp(-e);
    }
};

/// Given a molecule and ground state orbitals, solve the response equations
/// in the Tamm-Danchoff approximation.
class TDA 
{
   private:
      // Member variables

      // For timers
      std::vector<double> tda_ttt;
      std::vector<double> tda_sss;

      // Tensors for holding energies 
      // residuals, and shifts
      Tensor<double> tda_omega;        // Energies of response functions
      Tensor<double> tda_e_residuals;  // Residuals of energies

      // Information that is loaded from input file
      std::vector<real_function_3d> tda_orbitals;         // Ground state orbitals
      std::vector<real_function_3d> tda_act_orbitals;     // Ground state orbitals being used in calculation
      bool tda_spinrestricted;                            // Indicates open or closed shell
      Tensor<double> tda_ground_energies;                 // Ground state energies
      Tensor<double> tda_act_ground_energies;             // Ground state energies being used for calculation
      unsigned int tda_num_orbitals;                      // Number of ground state orbitals
      unsigned int tda_act_num_orbitals;                  // Number of ground state orbitals being used in calculation
      Tensor<double> tda_occ;                             // Occupied orbital occupation numbers
      double tda_L;                                       // Box size
      int tda_order;                                      // Wavelet order
      Molecule tda_molecule;                              // The molecule object from ground state calculation
      double tda_thresh;                                  // Derived from 'order' (thresh = 10^(-order + 2))
      double tda_small;                                   // Just chose 1e-4 (arbitrary, smallest distance to resolve)

      // Functions
      std::vector<std::vector<real_function_3d>> tda_x_response;   // Excited states to be solved for. 
                                                                   //    Note on storage: The response functions are calculated
                                                                   //    by calculating each transition of occupied --> virtual,
                                                                   //    and thus the actual response function is a sum of of all
                                                                   //    contributions to a specific virtual.

      // User input variables
      double tda_energy_threshold;
      double tda_range;
      int tda_max_iterations;
      int tda_num_excited; 
      int tda_print_level;    // Controls the amount and style of printing. Higher values print more
                              //   Values |   What gets printed
                              //   ----------------------------
                              //     0    |   Only convergence information from
                              //          |   each iteration. Default level.
                              //   ----------------------------
                              //     1    |   Changes style to print out each step
                              //          |   in the calculation, along with timings
                              //   ----------------------------
                              //     2    |   Debug level. Prints EVERYTHING!!!
   public:
      // Start a timer
      void start_timer(World & world);

      // Needed to do timers correctly
      double pop(std::vector<double> & v);

      // Stop a timer
      Tensor<double> end_timer(World & world);

      // Constructor, everything is initialized inside here
      TDA(World & world,           // MADNESS world object
          char* input_file,        // File with converged ground-state orbitals
                                   //    Note: You need to leave off the .00000 for this to work
          int k,                   // Number of desired response functions
          int print_level);        // Specifies how much printing should be done (see above)

      // Normalizes in the response sense
      void normalize(World & world,
                     std::vector<std::vector<real_function_3d>> & f); 

      // Prints norms of the given vector 
      void print_norms(World & world,
                       std::vector<std::vector<real_function_3d>> function); 

      // Prints relevant MADNESS parameters
      void print_madness_params(World & world);

      // Prints molecule geometry
      void print_molecule(World & world);

      // Prints response information
      void print_response_params(World & world);

      // Returns a set of vector of vector of real_function_3d of proper size, initialized to zero
      std::vector<std::vector<real_function_3d>> tda_zero_functions(World & world,
                                                                    int m,
                                                                    int n);

      // Returns a list of symmetry related functions for correct
      // pointgroup of the provided molecule
      std::vector<real_function_3d> symmetry(World & world);
     
      // Returns initial response functions
      std::vector<std::vector<real_function_3d>> create_trial_functions(World & world,
                                                                        int k,
                                                                        Tensor<double> energies,
                                                                        std::vector<real_function_3d> & orbitals,
                                                                        int print_level);

      // Returns the derivative of the coulomb operator, applied to ground state orbitals
      std::vector<std::vector<real_function_3d>> create_coulomb_derivative(World & world,
                                                                           std::vector<std::vector<real_function_3d>> & f,
                                                                           std::vector<real_function_3d> & orbitals,
                                                                           double small,
                                                                           double thresh);

      // Overloaded to get default values correct
      std::vector<std::vector<real_function_3d>> create_coulomb_derivative(World & world);

      // Returns the derivative of the exchange operator, applied to the ground state orbitals
      std::vector<std::vector<real_function_3d>> create_exchange_derivative(World & world,
                                                                            std::vector<std::vector<real_function_3d>> & f,
                                                                            std::vector<real_function_3d> & orbitals,
                                                                            double small,
                                                                            double thresh);

      // Overloaded to get default values correct
      std::vector<std::vector<real_function_3d>> create_exchange_derivative(World & world);

      // Returns gamma (the perturbed 2 electron piece)
      std::vector<std::vector<real_function_3d>> create_gamma(World & world,
                                                              std::vector<std::vector<real_function_3d>> & f,
                                                              std::vector<real_function_3d> & orbitals,
                                                              double small,
                                                              double thresh,
                                                              int print_level);

      // Overloaded to get default values correct 
      std::vector<std::vector<real_function_3d>> create_gamma(World & world);

      // Returns the coulomb potential of the ground state
      // Note: No post multiplication involved here
      real_function_3d coulomb(World & world);

      // Returns the result of ground state exchange applied to response functions
      std::vector<std::vector<real_function_3d>> exchange(World & world,
                                                          std::vector<std::vector<real_function_3d>> & f);

      // Overloaded to get defaults right 
      std::vector<std::vector<real_function_3d>> exchange(World & world);

      // Returns the ground state potential applied to response functions
      std::vector<std::vector<real_function_3d>> create_potential(World & world,
                                                                  std::vector<std::vector<real_function_3d>> & f,
                                                                  int print_level);

      // Returns the ground state potential applied to response functions
      std::vector<std::vector<real_function_3d>> create_potential(World & world);
      
      // Returns a tensor, where entry (i,j) = inner(a[i], b[i,j])
      Tensor<double> v_m_inner(World & world,
                               std::vector<real_function_3d> & a,
                               std::vector<std::vector<real_function_3d>> & b);

      // Returns a tensor, where entry (i,j) = inner(a[i], b[j]).sum()
      Tensor<double> expectation(World & world,
                                 std::vector<std::vector<real_function_3d>> & a,
                                 std::vector<std::vector<real_function_3d>> & b);
 
      // Returns the overlap matrix of the given response functions
      Tensor<double> create_S(World & world,
                              std::vector<std::vector<real_function_3d>> & f,
                              int print_level);

      // Overloaded to get default parameters correct 
      Tensor<double> create_S(World & world);

      // Returns the ground state fock operator applied to response functions
      std::vector<std::vector<real_function_3d>> create_fock(World & world,
                                                             std::vector<std::vector<real_function_3d>> & V,
                                                             std::vector<std::vector<real_function_3d>> & f,
                                                             int print_level);
                                                             
      // Returns the hamiltonian matrix, equation 45 from the paper
      Tensor<double> create_hamiltonian(World & world,
                                        std::vector<std::vector<real_function_3d>> & gamma,
                                        std::vector<std::vector<real_function_3d>> & V,
                                        std::vector<std::vector<real_function_3d>> & f,
                                        std::vector<real_function_3d> & ground_orbitals,
                                        std::vector<real_function_3d> & full_ground_orbitals,
                                        Tensor<double> & energies,
                                        int print_level);

      // Overloaded to get default parameters correct
      Tensor<double> create_hamiltonian(World & world,
                                        std::vector<std::vector<real_function_3d>> & gamma,
                                        std::vector<std::vector<real_function_3d>> & V);

      // Returns the shift needed for each orbital to make sure
      // -2.0 * (ground_state_energy + excited_state_energy) is positive
      Tensor<double> create_shift(World & world,
                                  int m,
                                  int n,
                                  Tensor<double> & ground,
                                  Tensor<double> & omega,
                                  int print_level);

      // Overloaded for parameter defaults
      Tensor<double> create_shift(World & world);

      // Returns the given shift applied to the given potentials
      std::vector<std::vector<real_function_3d>> apply_shift(World & world,
                                                             Tensor<double> & shifts,
                                                             std::vector<std::vector<real_function_3d>> & V,
                                                             std::vector<std::vector<real_function_3d>> & f);

      // Overloaded for default parameters 
      std::vector<std::vector<real_function_3d>> apply_shift(World & world,
                                                             Tensor<double> & shifts,
                                                             std::vector<std::vector<real_function_3d>> & V);


      // Returns a vector of BSH operators
      std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> create_bsh_operators(World & world,
                                                                                          Tensor<double> & shift,
                                                                                          Tensor<double> & ground,
                                                                                          Tensor<double> & omega,
                                                                                          double small,
                                                                                          double thresh);
      // Overloaded for parameter defaults 
      std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> create_bsh_operators(World & world,
                                                                                          Tensor<double> & shift);

      // Returns the second order update to the energy
      Tensor<double> calculate_energy_update(World & world,
                                             std::vector<std::vector<real_function_3d>> & gamma,
                                             std::vector<std::vector<real_function_3d>> & f_residuals,
                                             std::vector<std::vector<real_function_3d>> & new_f,
                                             int print_level);

      // Returns response functions that have been orthonormalized via
      // modified Gram-Schmidt. Note: This is specifically designed for
      // response functions only
      std::vector<std::vector<real_function_3d>> gram_schmidt(World & world,
                                                              std::vector<std::vector<real_function_3d>> & f);

      // Returns the max norm of the given vector of functions
      double calculate_max_residual(World & world, 
                                    std::vector<std::vector<real_function_3d>> & f);

      // Selects the 'active' orbitals from ground state orbitals to be used in the calculation (based on energy distance
      // from the HOMO.) Function has knowledge of tda_orbitals and tda_ground_energies. Function sets tda_act_orbitals and
      // tda_num_act_orbitals.
      void select_active_subspace(World & world);

      // Selects from a list of functions and energies the k functions with the lowest 
      // energy
      std::vector<std::vector<real_function_3d>> select_trial_functions(World & world,
                                                                        std::vector<std::vector<real_function_3d>> & f,
                                                                        Tensor<double> & energies,
                                                                        int k,
                                                                        int print_level); 

      // Calculates the exponentiation of a matrix through first order (I think)
      Tensor<double> matrix_exponential(const Tensor<double> & A);

      // Computes the unitary transformation that diagonalizes the fock matrix
      Tensor<double> get_fock_transformation(World & world,
                                             const Tensor<double> & overlap,
                                             Tensor<double> & fock,
                                             Tensor<double> & evals,
                                             const double thresh_degenerate);

      // Diagonalizes the fock matrix, taking care of degerate states
      Tensor<double> diag_fock_matrix(World & world,
                                      Tensor<double> & fock,
                                      std::vector<std::vector<real_function_3d>> & psi,
                                      std::vector<std::vector<real_function_3d>> & Vpsi,                     
                                      std::vector<std::vector<real_function_3d>> & gamma,                     
                                      Tensor<double> & evals,
                                      Tensor<double> & overlap,
                                      const double thresh); 

      // Transforms the given matrix of functions according to the given
      // transformation matrix. Used to update orbitals / potentials
      std::vector<std::vector<real_function_3d>> transform(World & world,
                                                           std::vector<std::vector<real_function_3d>> & f,
                                                           Tensor<double> & U);

      // Need to calculate: \sum_{j \neq k} x_p^{(j)} \Omega_{jk}
      // and add to RHS before BSH to allow correct for the rotated
      // potentials. Omega here is simply the Hamiltonian matrix.
      std::vector<std::vector<real_function_3d>> rotation_correction_term(World & world,
                                                                          std::vector<std::vector<real_function_3d>> f,
                                                                          Tensor<double> A,
                                                                          int print_level);

      // Sorts the given Tensor and vector of functions in place
      Tensor<int> sort(World & world,
                       Tensor<double> & vals,
                       Tensor<double> & vals_residuals,
                       std::vector<std::vector<real_function_3d>> & f, 
                       std::vector<std::vector<real_function_3d>> & f_diff); 

      // Prints iterate headers, used in default printing
      void print_iterate_headers(World & world);

      // Iterates the trial functions until covergence or it runs out of iterations
      void iterate(World & world);

      // Diagonalizes the given functions
      void diagonalize_guess(World & world,
                             std::vector<std::vector<real_function_3d>> & f,
                             Tensor<double> & omega,
                             std::vector<real_function_3d> & orbitals,
                             std::vector<real_function_3d> & full_orbitals,
                             Tensor<double> & energies,
                             double thresh,
                             double small,
                             int print_level);

      // Calculate polarizability tensor from response orbitals
      Tensor<double> polarizability(World & world,
                                    std::vector<std::vector<std::vector<real_function_3d>>> & f); 

      // Adds random noise to function f
      std::vector<std::vector<real_function_3d>> add_randomness(World & world,
                                                                std::vector<std::vector<real_function_3d>> & f);

      // Solves the response equations
      void solve(World & world);
};
#endif

// Dueces
