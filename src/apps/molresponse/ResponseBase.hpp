//
// Created by adrianhurtado on 1/24/22.
//

#ifndef MADNESS_RESPONSEBASE_HPP
#define MADNESS_RESPONSEBASE_HPP

#include <functional>
#include <numeric>
#include <utility>
#include <vector>

#include <madness/mra/nonlinsol.h>
#include <global_functions.h>
#include "madness/mra/functypedefs.h"
#include "madness/mra/mra.h"
#include<madness/chem/SCF.h>
#include "madness/tensor/tensor.h"
#include "madness/tensor/tensor_json.hpp"
#include "timer.h"
#include "x_space.h"


using namespace madness;


class inner_strategy {

public:
    virtual ~inner_strategy() = default;
    virtual Tensor<double> compute_inner(const X_space &x, const X_space &y) const = 0;
};

class Context {

private:
    std::unique_ptr<inner_strategy> strategy_;

public:
    explicit Context(std::unique_ptr<inner_strategy> &&strategy = {})
        : strategy_(std::move(strategy)) {}
    void set_strategy(std::unique_ptr<inner_strategy> &&strategy) {
        strategy_ = std::move(strategy);
    }
    void print_inner(const X_space &x, const X_space &y) const {
        if (strategy_) {
            std::cout << "Context: Computing inner using the strategy (not sure how it'll do it)\n";
            auto result = strategy_->compute_inner(x, y);
            std::cout << result << "\n";
        } else {
            std::cout << "Context: Strategy isn't set\n";
        }
    }

    Tensor<double> inner(const X_space &x, const X_space &y) const {
        if (strategy_) {
            return strategy_->compute_inner(x, y);
        } else {
            throw madness::MadnessException("Inner product Stratgey isn't set",
                                            "Need to set a strategy", 2, 455, "inner",
                                            "ResponseBase.hpp");
        }
    }
};

class full_inner_product : public inner_strategy {
public:
    Tensor<double> compute_inner(const X_space &x, const X_space &y) const override {
        return inner(x, y);
    }
};

class static_inner_product : public inner_strategy {
public:
    [[nodiscard]] Tensor<double> compute_inner(const X_space &x, const X_space &y) const override {
        return response_space_inner(x.x, y.x);
    }
};
typedef std::vector<XNonlinearSolver<vector_real_function_3d, double, response_matrix_allocator>>
        response_solver;
typedef std::vector<XNonlinearSolver<real_function_3d, double, response_function_allocator>>
        response_function_solver;

class response_timing {
    std::map<std::string, std::vector<double>> wall_time_data;
    std::map<std::string, std::vector<double>> cpu_time_data;
    int iter;

public:
    response_timing();

    void to_json(json &j);

    void print_data();

    void add_data(std::map<std::string, std::pair<double, double>> values);
};
class response_data {
    std::map<std::string, std::vector<Tensor<double>>> function_data;
    int iter;

public:
    response_data();

    void to_json(json &j);

    void add_data(std::map<std::string, Tensor<double>> values);
};

class ResponseTester;

struct residuals {
    X_space residual;
    Tensor<double> residual_norms;
};


using gamma_orbitals = std::tuple<X_space, vector_real_function_3d>;

class ResponseBase {
public:
    friend ResponseTester;

    ResponseBase(World &world, const CalcParams &params);

    void solve(World &world);

    virtual void initialize(World &world) = 0;

    virtual void iterate(World &world) = 0;

    //virtual void iterate();
    auto get_parameter() const -> CalcParams { return {ground_calc, molecule, r_params}; }

    auto get_orbitals() const -> vector_real_function_3d { return ground_orbitals; }

    auto get_chi() const -> X_space { return Chi.copy(); };

    void output_json();

    json j_molresponse{};
    response_timing time_data;
    response_data function_data;
    mutable std::map<std::string, std::pair<double, double>> iter_timing;
    mutable std::map<std::string, Tensor<double>> iter_function_data;

    Context response_context;

protected:
    // Given molecule returns the nuclear potential of the molecule
    ResponseParameters r_params;
    Molecule molecule;
    GroundStateCalculation ground_calc;
    bool all_done = false;


    XCfunctional xcf;
    real_function_3d mask;

    std::shared_ptr<PotentialManager> potential_manager;
    // shared pointers to Operators
    poperatorT shared_coulomb_operator;// shared pointer to seperated convolution operator
    std::vector<poperatorT> coul_ops;
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
    void set_protocol(World &world, double thresh) {
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
        shared_coulomb_operator = poperatorT(CoulombOperatorPtr(world, r_params.lo(), thresh));
        coul_ops.clear();
        for (int i = 0; i < r_params.num_states(); i++) {
            coul_ops.push_back(poperatorT(
                    CoulombOperatorPtr(world, r_params.lo(), FunctionDefaults<3>::get_thresh())));
        }
        gradop = gradient_operator<double, 3>(world);
        potential_manager = std::make_shared<PotentialManager>(molecule, "a");
        potential_manager->make_nuclear_potential(world);
        // GaussianConvolution1DCache<double>::map.clear();//(TODO:molresponse-What
        // Create the masking function
        mask = real_function_3d(real_factory_3d(world).f(mask3).initial_level(4).norefine());

        ground_density = make_ground_density(world);
        ground_density.truncate(FunctionDefaults<3>::get_thresh());
        // Basic print
        if (world.rank() == 0) {
            print("\nSolving NDIM=", 3, " with thresh", thresh, "    k",
                  FunctionDefaults<3>::get_k(), "  dconv", std::max(thresh, r_params.dconv()),
                  "\n");
        }
    }

    virtual void check_k(World &world, double thresh, int k);

    auto make_ground_density(World &world) const -> functionT;

    auto ComputeHamiltonianPair(World &world) const -> std::pair<Tensor<double>, Tensor<double>>;

    auto Coulomb(World &world) const -> real_function_3d;

    auto make_xc_operator(World &world) const -> XCOperator<double, 3>;

    virtual void save(World &world, const std::string &name) = 0;

    virtual void load(World &world, const std::string &name) = 0;

    auto make_density(World &world, const X_space &chi) const -> vecfuncT;

    void load_balance_chi(World &world);

    auto make_bsh_operators_response(World &world, double &shift, const double omega) const
            -> vector<poperatorT>;


    auto kain_x_space_update(World &world, const X_space &chi, const X_space &residual_chi,
                             response_solver &kain_x_space) -> X_space;

    void x_space_step_restriction(World &world, const X_space &old_Chi, X_space &temp,
                                  bool restrict_y, const double &max_bsh_rotation);

    void plotResponseOrbitals(World &world, size_t iteration, const response_space &x_response,
                              const response_space &y_response,
                              const ResponseParameters &responseParameters,
                              const GroundStateCalculation &g_params);

    static auto orbital_load_balance(World &world, const gamma_orbitals &, double load_balance)
            -> gamma_orbitals;

    auto compute_gamma_tda(World &world, const gamma_orbitals &density,
                           const XCOperator<double, 3> &xc) const -> X_space;

    auto compute_gamma_static(World &world, const gamma_orbitals &,
                              const XCOperator<double, 3> &xc) const -> X_space;

    auto compute_gamma_full(World &world, const gamma_orbitals &,
                            const XCOperator<double, 3> &xc) const -> X_space;

    auto compute_V0X(World &world, const X_space &X, const XCOperator<double, 3> &xc,
                     bool compute_Y) const -> X_space;

    auto compute_lambda_X(World &world, const X_space &chi, XCOperator<double, 3> &xc,
                          const std::string &calc_type) const -> X_space;

    auto compute_theta_X(World &world, const X_space &chi, const XCOperator<double, 3> &xc,
                         const std::string &calc_type) const -> X_space;

    auto compute_F0X(World &world, const X_space &X, const XCOperator<double, 3> &xc,
                     bool compute_Y) const -> X_space;

    void analyze_vectors(World &world, const vecfuncT &x, const std::string &response_state);

    auto project_ao_basis(World &world, const AtomicBasisSet &aobasis) -> vecfuncT;


    static auto project_ao_basis_only(World &world, const AtomicBasisSet &aobasis,
                                      const Molecule &mol) -> vecfuncT;

    void converged_to_json(json &j);

    auto update_residual(World &world, const X_space &chi, const X_space &g_chi,
                         const std::string &calc_type, const Tensor<double> &old_residuals)
            -> residuals;

    auto compute_response_potentials(World &world, const X_space &chi, XCOperator<double, 3> &xc,
                                     const std::string &calc_type) const
            -> std::tuple<X_space, X_space, X_space>;

    // compute exchange |i><i|J|p>
    auto exchangeHF(const vecfuncT &ket, const vecfuncT &bra, const vecfuncT &vf) const
            -> vecfuncT {
        World &world = ket[0].world();
        auto n = bra.size();
        auto nf = ket.size();
        double tol = FunctionDefaults<3>::get_thresh();/// Important this is
        double mul_tol = 0.0;
        const double lo = r_params.lo();


        std::shared_ptr<real_convolution_3d> poisson;
        /// consistent with Coulomb
        vecfuncT Kf = zero_functions_compressed<double, 3>(world, nf, true);

        reconstruct(world, bra);
        reconstruct(world, ket);
        reconstruct(world, vf);

        // i-j sym
        for (int i = 0; i < n; ++i) {
            // for each |i> <i|phi>
            vecfuncT psi_f = mul_sparse(world, bra[i], vf, mul_tol, true);/// was vtol
            truncate(world, psi_f, tol, true);
            // apply to vector of products <i|phi>..<i|1> <i|2>...<i|N>
            psi_f = apply(world, *shared_coulomb_operator, psi_f);
            truncate(world, psi_f, tol, true);
            // multiply by ket i  <i|phi>|i>: <i|1>|i> <i|2>|i> <i|2>|i>
            psi_f = mul_sparse(world, ket[i], psi_f, mul_tol, true);/// was vtol
            /// Generalized A*X+y for vectors of functions ---- a[i] = alpha*a[i] +
            // 1*Kf+occ[i]*psi_f
            gaxpy(world, double(1.0), Kf, double(1.0), psi_f);
        }
        truncate(world, Kf, tol, true);
        return Kf;
    }

    static void print_inner(World &world, const std::string &name, const X_space &left,
                            const X_space &right);

    void function_data_to_json(json &j_mol_in, size_t iter, const Tensor<double> &x_norms,
                               const Tensor<double> &x_abs_norms, const Tensor<double> &rho_norms,
                               const Tensor<double> &rho_res_norms);
    X_space compute_TX(World &world, const X_space &X, bool compute_Y) const;
    vecfuncT update_density(World &world, const X_space &chi, const vecfuncT &old_density) const;
};


// Some helper functions
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Check k given response parameters
/// \param world
/// \param Chi
/// \param thresh
/// \param k
void check_k(World &world, X_space &Chi, double thresh, int k);


auto add_randomness(World &world, const response_space &f, double magnitude) -> response_space;

void normalize(World &world, response_space &f);

void normalize(World &world, X_space &Chi);


static auto kronecker(size_t l, size_t n) -> double {
    if (l == n) return 1.0;
    return 0.0;
}

auto solid_harmonics(World &world, int n) -> std::map<std::vector<int>, real_function_3d>;

/***
 * @brief Prints the norms of the functions of a response space
 *
 *
 *
 * @param world
 * @param f
 */
void print_norms(World &world, const response_space &f);

// Returns a list of solid harmonics such that:
// solid_harm.size() * num_ground_orbs > 2 * num. resp. components
auto make_xyz_functions(World &world) -> vector_real_function_3d;

// Selects from a list of functions and energies the k functions with the
// lowest energy
auto select_functions(World &world, response_space f, Tensor<double> &energies, size_t k,
                      size_t print_level) -> response_space;

// Sorts the given tensor of eigenvalues and
// response functions
void sort(World &world, Tensor<double> &vals, response_space &f);

void sort(World &world, Tensor<double> &vals, X_space &f);


// Specialized for response calculations that returns orthonormalized
// functions
auto gram_schmidt(World &world, const response_space &f) -> response_space;

/// Computes the transition density between set of two response functions x and y.
/// Uses std::transform to iterate between x and y vectors
/// \param world
/// \param orbitals
/// \param x
/// \param y
/// \return
auto transition_density(World &world, const vector_real_function_3d &orbitals,
                        const response_space &x, const response_space &y)
        -> vector_real_function_3d;

auto transition_densityTDA(World &world, const vector_real_function_3d &orbitals,
                           const response_space &x) -> vector_real_function_3d;

auto transform(World &world, const response_space &f, const Tensor<double> &U) -> response_space;

auto transform(World &world, const X_space &x, const Tensor<double> &U) -> X_space;

// result(i,j) = inner(a[i],b[j]).sum()
auto expectation(World &world, const response_space &A, const response_space &B) -> Tensor<double>;

void inner_to_json(World &world, const std::string &name, const Tensor<double> &m_val,
                   std::map<std::string, Tensor<double>> &data);
class ResponseTester {

public:
    static void load_calc(World &world, ResponseBase *p, double thresh) {
        p->set_protocol(world, thresh);
        p->load(world, p->r_params.restart_file());
        p->check_k(world, thresh, FunctionDefaults<3>::get_k());
    }

    static X_space compute_gamma_full(World &world, ResponseBase *p, double thresh) {
        XCOperator<double, 3> xc = p->make_xc_operator(world);
        X_space gamma = p->compute_gamma_full(world, {p->Chi, p->ground_orbitals}, xc);
        return gamma;
    }

    X_space compute_lambda_X(World &world, ResponseBase *p, double thresh) {
        XCOperator<double, 3> xc = p->make_xc_operator(world);
        X_space gamma = p->compute_lambda_X(world, p->Chi, xc, p->r_params.calc_type());
        return gamma;
    }

    std::pair<X_space, X_space> compute_VFOX(World &world, ResponseBase *p, bool compute_y) {
        XCOperator<double, 3> xc = p->make_xc_operator(world);
        X_space V = p->compute_V0X(world, p->Chi, xc, compute_y);
        X_space F = p->compute_F0X(world, p->Chi, xc, compute_y);
        return {V, F};
    }
};


#endif// MADNESS_RESPONSEBASE_HPP
