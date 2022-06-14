//
// Created by adrianhurtado on 1/24/22.
//

#ifndef MADNESS_RESPONSEBASE_HPP
#define MADNESS_RESPONSEBASE_HPP

#include <functional>
#include <numeric>
#include <utility>
#include <vector>

#include "global_functions.h"
#include "load_balance.h"
#include "madness/mra/functypedefs.h"
#include "madness/mra/mra.h"
#include "madness/tensor/tensor.h"
#include "timer.h"
#include "x_space.h"


using namespace madness;

class response_timing {
    std::map<std::string, std::vector<double>> wall_time_data;
    std::map<std::string, std::vector<double>> cpu_time_data;
    int iter;

public:
    response_timing();

    void to_json(json& j);

    void print_data();

    void add_data(std::map<std::string, std::pair<double, double>> values);
};
class ResponseTester;

struct residuals {

    X_space residual;
    Tensor<double> x;
    Tensor<double> y;
};


using gamma_orbitals = std::tuple<X_space, vector_real_function_3d, vector_real_function_3d>;

class ResponseBase {
public:
    friend ResponseTester;
    ResponseBase(World& world, const CalcParams& params);
    void solve(World& world);
    virtual void initialize(World& world) = 0;
    virtual void iterate(World& world) = 0;
    //virtual void iterate();
    CalcParams get_parameter() const { return {ground_calc, molecule, r_params}; }
    vector_real_function_3d get_orbitals() const { return ground_orbitals; }
    void output_json();

    json j_molresponse{};
    response_timing time_data;
    mutable std::map<std::string, std::pair<double, double>> iter_timing;

protected:
    // Given molecule returns the nuclear potential of the molecule
    ResponseParameters r_params;
    Molecule molecule;
    GroundStateCalculation ground_calc;
    bool converged = false;


    XCfunctional xcf;
    real_function_3d mask;

    std::shared_ptr<PotentialManager> potential_manager;
    // shared pointers to Operators
    poperatorT coulop;// shared pointer to seperated convolution operator
    std::vector<std::shared_ptr<real_derivative_3d>> gradop;

    // Stored functions
    mutable real_function_3d stored_v_nuc; // Stored nuclear potential from ground state
    mutable real_function_3d stored_v_coul;// Stored coulomb potential from ground state

    // Ground state orbitals and energies
    vector_real_function_3d ground_orbitals{};
    Tensor<double> ground_energies;

    // Information that is inferred from input file
    // Ground state orbitals being used in calculation
    Tensor<double> hamiltonian;
    Tensor<double> ham_no_diag;
    // Tensors for holding energies
    // residuals, and shifts

    Tensor<double> e_residuals;// Residuals of energies

    // Mask function to handle boundary conditions

    functionT ground_density;// ground state density

    mutable response_space stored_potential;// The ground state potential, stored only
                                            // if store_potential is true (default is

    double vtol{};

    X_space Chi;


    /// Sets the Function protocol dependent on the truncation threshold.
    /// Sets the polynomial order of basis functions k
    /// Then creates shared coulomb operator, gradient operator, ground density
    /// AS well as the ground state density
    /// \param world
    /// \param thresh
    void set_protocol(World& world, double thresh) {
        int k;
        // Allow for imprecise conversion of threshold
        if (thresh >= 0.9e-2) k = 4;
        else if (thresh >= 0.9e-4)
            k = 6;
        else if (thresh >= 0.9e-6)
            k = 8;
        else if (thresh >= 0.9e-8)
            k = 10;
        else
            k = 12;

        // k defaults to make sense with thresh, override by providing k in
        // input file
        if (r_params.k() == -1) {
            FunctionDefaults<3>::set_k(k);
        } else {
            FunctionDefaults<3>::set_k(r_params.k());
        }

        // Set Function Defaults
        FunctionDefaults<3>::set_thresh(thresh);
        FunctionDefaults<3>::set_refine(true);
        FunctionDefaults<3>::set_initial_level(2);


        FunctionDefaults<3>::set_autorefine(false);
        FunctionDefaults<3>::set_apply_randomize(false);
        FunctionDefaults<3>::set_project_randomize(false);
        GaussianConvolution1DCache<double>::map.clear();

        double safety = 0.1;
        vtol = FunctionDefaults<3>::get_thresh() * safety;
        coulop = poperatorT(CoulombOperatorPtr(world, r_params.lo(), thresh));
        gradop = gradient_operator<double, 3>(world);
        potential_manager = std::make_shared<PotentialManager>(molecule, "a");
        potential_manager->make_nuclear_potential(world);
        // GaussianConvolution1DCache<double>::map.clear();//(TODO:molresponse-What
        // Create the masking function
        mask = real_function_3d(real_factory_3d(world).f(mask3).initial_level(4).norefine());
        ground_density = make_ground_density(world);
        // Basic print
        if (world.rank() == 0) {
            print("\nSolving NDIM=", 3, " with thresh", thresh, "    k",
                  FunctionDefaults<3>::get_k(), "  dconv", std::max(thresh, r_params.dconv()),
                  "\n");
        }
    }

    virtual void check_k(World& world, double thresh, int k);
    functionT make_ground_density(World& world) const;
    std::pair<Tensor<double>, Tensor<double>> ComputeHamiltonianPair(World& world) const;
    real_function_3d Coulomb(World& world) const;
    XCOperator<double, 3> make_xc_operator(World& world) const;
    virtual void save(World& world, const std::string& name) = 0;
    virtual void load(World& world, const std::string& name) = 0;
    vecfuncT make_density(World& world, const X_space& chi) const;

    void load_balance_chi(World& world);
    vector<poperatorT> make_bsh_operators_response(World& world, double& shift,
                                                   double& omega) const;


    X_space compute_residual(World& world, X_space& old_Chi, X_space& temp,
                             Tensor<double>& bsh_residualsX, Tensor<double>& bsh_residualsY,
                             std::string calc_type);
    X_space kain_x_space_update(World& world, const X_space& temp, const X_space& res,
                                NonLinearXsolver& kain_x_space, vector<X_vector>& Xvector,
                                vector<X_vector>& Xresidual);
    void x_space_step_restriction(World& world, X_space& old_Chi, X_space& temp, bool restrict_y,
                                  const double& maxrotn);
    void vector_stats(const vector<double>& v, double& rms, double& maxabsval) const;
    double do_step_restriction(World& world, const vector_real_function_3d& x,
                               vector_real_function_3d& x_new, std::string spin) const;
    double do_step_restriction(World& world, const vector_real_function_3d& x,
                               vector_real_function_3d& x_new, std::string spin,
                               double maxrotn) const;
    double do_step_restriction(World& world, const vector_real_function_3d& x,
                               const vector_real_function_3d& y, vector_real_function_3d& x_new,
                               vector_real_function_3d& y_new, std::string spin) const;
    void PlotGroundandResponseOrbitals(World& world, size_t iteration, response_space& x_response,
                                       response_space& y_response,
                                       const ResponseParameters& r_params,
                                       const GroundStateCalculation& g_params);
    void vector_stats_new(const Tensor<double> v, double& rms, double& maxabsval) const;

    static gamma_orbitals orbital_load_balance(World& world, const gamma_orbitals&,
                                               const double load_balance);
    X_space compute_gamma_tda(World& world, const gamma_orbitals& density,
                              const XCOperator<double, 3>& xc) const;
    X_space compute_gamma_static(World& world, const gamma_orbitals&,
                                 const XCOperator<double, 3>& xc) const;
    X_space compute_gamma_full(World& world, const gamma_orbitals&,
                               const XCOperator<double, 3>& xc) const;
    X_space compute_V0X(World& world, const X_space& X, const XCOperator<double, 3>& xc,
                        bool compute_Y) const;
    X_space compute_lambda_X(World& world, const X_space& chi, XCOperator<double, 3>& xc,
                             const std::string& calc_type) const;
    X_space compute_theta_X(World& world, const X_space& chi, XCOperator<double, 3> xc,
                            const std::string& calc_type) const;
    X_space compute_F0X(World& world, const X_space& X, const XCOperator<double, 3>& xc,
                        bool compute_Y) const;
    void analyze_vectors(World& world, const vecfuncT& x, const std::string& response_state);
    vecfuncT project_ao_basis(World& world, const AtomicBasisSet& aobasis);


    vecfuncT project_ao_basis_only(World& world, const AtomicBasisSet& aobasis,
                                   const Molecule& molecule);
    void converged_to_json(json& j);
    residuals compute_residual(World& world, X_space& old_Chi, X_space& temp,
                               std::string calc_type);
    std::tuple<X_space, X_space, X_space> compute_response_potentials(
            World& world, const X_space& chi, XCOperator<double, 3>& xc,
            const std::string& calc_type) const;
};


// Some helper functions
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Check k given response parameters
/// \param world
/// \param Chi
/// \param thresh
/// \param k
void check_k(World& world, X_space& Chi, double thresh, int k);


response_space add_randomness(World& world, const response_space& f, double magnitude);
void normalize(World& world, response_space& f);
void normalize(World& world, X_space& Chi);


static double kronecker(size_t l, size_t n) {
    if (l == n) return 1.0;
    return 0.0;
}

std::map<std::vector<int>, real_function_3d> solid_harmonics(World& world, int n);

/***
 * @brief Prints the norms of the functions of a response space
 *
 *
 *
 * @param world
 * @param f
 */
void print_norms(World& world, const response_space& f);
// Returns a list of solid harmonics such that:
// solid_harm.size() * num_ground_orbs > 2 * num. resp. components
vector_real_function_3d make_xyz_functions(World& world);
// Selects from a list of functions and energies the k functions with the
// lowest energy
response_space select_functions(World& world, response_space f, Tensor<double>& energies, size_t k,
                                size_t print_level);
// Sorts the given tensor of eigenvalues and
// response functions
void sort(World& world, Tensor<double>& vals, response_space& f);
void sort(World& world, Tensor<double>& vals, X_space& f);


// Specialized for response calculations that returns orthonormalized
// functions
response_space gram_schmidt(World& world, const response_space& f);
// Specialized for response calculations that returns orthonormalized
// functions
void gram_schmidt(World& world, response_space& f, response_space& g);
/// Computes the transition density between set of two response functions x and y.
/// Uses std::transform to iterate between x and y vectors
/// \param world
/// \param orbitals
/// \param x
/// \param y
/// \return
vector_real_function_3d transition_density(World& world, const vector_real_function_3d& orbitals,
                                           const response_space& x, const response_space& y);

vector_real_function_3d transition_densityTDA(World& world, const vector_real_function_3d& orbitals,
                                              const response_space& x);

response_space transform(World& world, const response_space& f, const Tensor<double>& U);

X_space transform(World &world, const X_space &x, const Tensor<double> &U);

// result(i,j) = inner(a[i],b[j]).sum()
Tensor<double> expectation(World& world, const response_space& A, const response_space& B);


class ResponseTester {

public:
    void load_calc(World& world, ResponseBase* p, double thresh) {
        p->set_protocol(world, thresh);
        p->load(world, p->r_params.restart_file());
        p->check_k(world, thresh, FunctionDefaults<3>::get_k());
    }
    X_space compute_gamma_full(World& world, ResponseBase* p, double thresh) {
        XCOperator<double, 3> xc = p->make_xc_operator(world);
        X_space gamma =
                p->compute_gamma_full(world, {p->Chi, p->ground_orbitals, p->ground_orbitals}, xc);
        return gamma;
    }
    X_space compute_lambda_X(World& world, ResponseBase* p, double thresh) {
        XCOperator<double, 3> xc = p->make_xc_operator(world);
        X_space gamma = p->compute_lambda_X(world, p->Chi, xc, p->r_params.calc_type());
        return gamma;
    }
    std::pair<X_space, X_space> compute_VFOX(World& world, ResponseBase* p, bool compute_y) {
        XCOperator<double, 3> xc = p->make_xc_operator(world);
        X_space V = p->compute_V0X(world, p->Chi, xc, compute_y);
        X_space F = p->compute_F0X(world, p->Chi, xc, compute_y);
        return {V, F};
    }
};
#endif// MADNESS_RESPONSEBASE_HPP
