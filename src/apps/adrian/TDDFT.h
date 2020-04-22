
#ifndef MADNESS_APPS_TDDFT_H_INCLUDED
#define MADNESS_APPS_TDDFT_H_INCLUDED

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/constants.h>
#include <madness/mra/nonlinsol.h> // The kain solver
#include <vector>
#include <math.h>
#include <stdio.h>
#include <iomanip>
#include <complex>
#include <cmath>
#include <random>
#include <algorithm>
#include "../chem/molecule.h"
#include "../chem/SCFOperators.h"
#include "../chem/xcfunctional.h"
#include "ResponseParameters.h"
#include "GroundParameters.h"
#include "ResponseFunction2.h"
#include "ResponsePotential.h"

using namespace madness;

class BS_MomentFunctor : public FunctionFunctorInterface<double, 3>
{
private:
    const int i, j, k;

public:
    BS_MomentFunctor(int i, int j, int k) : i(i), j(j), k(k) {}
    BS_MomentFunctor(const std::vector<int> &x) : i(x[0]), j(x[1]), k(x[2]) {}
    double operator()(const Vector<double, 3> &r) const
    {
        double xi = 1.0, yj = 1.0, zk = 1.0;
        for (int p = 0; p < i; ++p)
            xi *= r[0];
        for (int p = 0; p < j; ++p)
            yj *= r[1];
        for (int p = 0; p < k; ++p)
            zk *= r[2];
        return xi * yj * zk;
    }
};
// Functor for moment...not sure what this is for yet

/// an N-dimensional real-valued Gaussian function

/// the function looks like
/// \[
/// f(r) = x^i y^j .. z^k exp(-alpha r^2)
/// j\]
template <std::size_t NDIM>
class GaussianGuess : public FunctionFunctorInterface<double, NDIM>
{
    typedef Vector<double, NDIM> coordT;

public:
    /// ctor

    /// @param[in]  origin  the origin of the Gauss function
    /// @param[in]  alpha   the exponent exp(-alpha r^2)
    /// @param[in]  ijk     the monomial x^i y^j z^k exp(-alpha r^2) (for NDIM)
    GaussianGuess(const coordT &origin, const double alpha,
                  const std::vector<int> ijk = std::vector<int>(NDIM))
        : origin(origin), exponent(alpha), ijk(ijk)
    {
    }

    coordT origin;
    double exponent;      ///< exponent of the guess
    std::vector<int> ijk; ///< cartesian exponents

    double operator()(const coordT &xyz) const
    {
        double arg = 0.0, prefac = 1.0;
        for (std::size_t i = 0; i < NDIM; ++i)
        {
            arg += (xyz[i] - origin[i]) * (xyz[i] - origin[i]);
            prefac *= pow(xyz[i], ijk[i]);
        }
        const double e = exponent * arg;
        return prefac * exp(-e);
    }
};
// returns a gaussian guess

class TDDFT
{
    // Member Variables
    typedef std::vector<real_function_3d> vecOrbs;
    typedef ResponseFunction RF;

public:
    ResponseParameters Rparams;

private:
    GroundParameters Gparams;
    // vector
    vecOrbs act_orbitals; // {phi_i}
    Tensor<double> omega;
    Tensor<double> e_residuals;
    Tensor<double> act_ground_energies; //{epsilon_i}
    Tensor<double> hamiltonian;         //---not sure what this means yet
    Tensor<double> ham_no_diag;         //---no diagonal

    //XCFunctional object of DFT
    //Functions
    real_function_3d mask;         // BC
    real_function_3d stored_v_nuc; // vnuc
    //Response Functions variables
    RF x_response;
    RF y_response;
    RF stored_potential;

public:
    //Constructors
    TDDFT(World &world, const char *input_file);
    TDDFT(World &world, std::shared_ptr<std::istream> input);

    // save load
    void save(World &world);
    void load(World &world, std::string name);

    //normalize
    void normalize(World &world, RF &f);
    void normalize(World &world, RF &f, RF &g);

    //prints
    void print_norms(World &world, RF function);
    void print_molecule(World &world);

    // Return a list of solid harmonics---???
    std::map<std::vector<int>, real_function_3d> solid_harmonics(World &world, int n);

    // Initial Response functions
    RF response_zero_functions(World &world, int m, int n); //zero functions
    RF create_trial_functions(World *world, int k, vecOrbs &orbitals, int print_level);
    RF dipole_guess(World &world, vecOrbs orbitals);
    //Create Derivatives  dJ/dn dVxc/dn
    RF create_coulomb_derivative(World &world, RF &f, vecOrbs &orbitals, double small, double thresh);
    // Only for TDA
    RF create_exchange_derivative(World &world, RF &f, vecOrbs &orbitals, double small, double thresh);

    // Diagonal Elements of response matrix
    RF create_A(World &world, RF &fe, RF &gamma, RF &V, RF &f, vecOrbs ground_orbitals,
                Tensor<double> &hamiltonian, int print_level, std::string xy);

    // Diagonal Elements of response matrix  off Diag B
    RF create_B(World &world, RF &f, RF &g, vecOrbs &orbitals, double small, double thresh, int print_level);

    RF create_gamma(World &world, RF &f, RF &g, vecOrbs &orbitals, double small, double thresh, int print_level, std::string xy);

    // Returns the coulomb potential of the ground state
    real_function_3d coulomb(World &world);

    // Returns the result of the ground state
    RF exchange(World &world, RF &f);

    // GS exchange applied to Response Function
    RF create_potential(World &world, RF &f, XCOperator xc, int print_level, std::string xy);

    // Return a tensor, where entry (i,j)=inner(a[i],b[j]).sum()
    Tensor<double> expection(World &world, RF *a, RF *b);

    //Returns the GS Operator appied to response function
    RF create_fock(World & world, RF &V, RF &f, int print_level, std::string xy);

    // Returns the hamiltonian matrix, equation 45 from the paper
    Tensor<double> create_response_matrix(World & world, RF &fe, RF &gamma, RF &V, RF & f, vecOrbs &ground_orbitals, Tensor<double> & energies, int print_level, std::string xy);

    // Construct the full response matrix of
    // [  A     B ] [ X ] = w [ X ]
    // [ -B   -A ] [ Y ] =    [ Y ]

        Tensor<double> create_full_response_matrix(World & world, 
                                                 ResponseFunction & x_b,
                                                 ResponseFunction & Vx,
                                                 ResponseFunction & B_x,
                                                 ResponseFunction & fe_x,
                                                 ResponseFunction & x,
                                                 ResponseFunction & y_b,
                                                 ResponseFunction & Vy,
                                                 ResponseFunction & B_y,
                                                 ResponseFunction & fe_y,
                                                 ResponseFunction & y,
                                                 std::vector<real_function_3d> & ground_orbitals,
                                                 Tensor<double> & ground_ham,
                                                 double small,
                                                 double thresh,
                                                 int print_level);
    

    
    
};

#endif