//
// Created by adrianhurtado on 1/24/22.
//

#ifndef MADNESS_RESPONSEBASE_HPP
#define MADNESS_RESPONSEBASE_HPP

#include "x_space.h"
#include "global_functions.h"
#include "load_balance.h"
#include "timer.h"
#include <utility>

class ResponseBase {
public:
    ResponseBase(World& world, const CalcParams& params);
    void solve(World& world);
    virtual void initialize(World &world)=0;
    //virtual void iterate();
    CalcParams get_parameter() const { return {ground_calc, molecule, r_params}; }
    vector_real_function_3d get_orbitals() const { return ground_orbitals; }

protected:
    // Given molecule returns the nuclear potential of the molecule
    ResponseParameters r_params;
    Molecule molecule;
    GroundStateCalculation ground_calc;



    XCfunctional xcf;
    static real_function_3d mask;

    std::shared_ptr<PotentialManager> potential_manager;
    // shared pointers to Operators
    poperatorT coulop;// shared pointer to seperated convolution operator
    std::vector<std::shared_ptr<real_derivative_3d>> gradop{};

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
    Tensor<double> omega;      // Energies of response functions
    Tensor<double> e_residuals;// Residuals of energies

    // Mask function to handle boundary conditions

    functionT ground_density;// ground state density
    vecfuncT rho_omega{};      // response density

    mutable response_space stored_potential;// The ground state potential, stored only
                                            // if store_potential is true (default is

    double vtol{};
    json j_molresponse{};

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

    void check_k(World& world, double thresh, int k);
    functionT make_ground_density(World& world) const;
    std::pair<Tensor<double>,Tensor<double>> ComputeHamiltonianPair(World& world) const;
    real_function_3d Coulomb(World& world) const;
    XCOperator<double, 3> make_xc_operator(World& world) const;
    void save(World& world, const std::string& name);
    void load(World& world, const std::string& name);
    vecfuncT make_density(World& world);

    vector<real_function_3d> transition_density(World& world, vector<real_function_3d>& orbitals,
                                                response_space& x, response_space& y);
    std::vector<real_function_3d> transition_densityTDA(World& world,
                                                        const response_space& x)const ;
    void load_balance(World& world);
    vector<poperatorT> make_bsh_operators_response(World& world, double& shift,
                                                   double& omega) const;
    X_space Compute_Theta_X(World& world, X_space& Chi, XCOperator<double, 3> xc,
                            std::string calc_type);
    X_space compute_F0X(World& world, X_space& Chi, XCOperator<double, 3> xc, bool compute_Y);
    X_space compute_V0X(World& world, X_space& X, XCOperator<double, 3> xc, bool compute_Y);
    void orbital_load_balance(World& world, vecfuncT& psi0, vecfuncT& psi0_copy, X_space& X,
                              X_space& Chi_copy);
    X_space compute_gamma_tda(World& world, X_space& X, XCOperator<double, 3> xc);
    X_space compute_gamma_full(World& world, X_space& X, const XCOperator<double, 3>& xc);
    X_space compute_gamma_static(World& world, X_space& X, XCOperator<double, 3> xc);
    X_space compute_residual(World& world, X_space& old_Chi, X_space& temp,
                             Tensor<double>& bsh_residualsX, Tensor<double>& bsh_residualsY,
                             std::string calc_type);
    X_space kain_x_space_update(World& world, const X_space& temp, const X_space& res,
                                NonLinearXsolver& kain_x_space, vector<X_vector>& Xvector,
                                vector<X_vector>& Xresidual);
    void x_space_step_restriction(World& world, X_space& old_Chi, X_space& temp, bool restrict_y,
                                  Tensor<double>& maxrotn);
    void vector_stats(const vector<double>& v, double& rms, double& maxabsval) const;
    double do_step_restriction(World& world, const vecfuncT& x, vecfuncT& x_new,
                               std::string spin) const;
    double do_step_restriction(World& world, const vecfuncT& x, vecfuncT& x_new, std::string spin,
                               double maxrotn) const;
    double do_step_restriction(World& world, const vecfuncT& x, const vecfuncT& y, vecfuncT& x_new,
                               vecfuncT& y_new, std::string spin) const;
    void PlotGroundandResponseOrbitals(World& world, size_t iteration, response_space& x_response,
                                       response_space& y_response,
                                       const ResponseParameters& r_params,
                                       const GroundStateCalculation& g_params);
    void vector_stats_new(const Tensor<double> v, double& rms, double& maxabsval) const;
};

// Some helper functions
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Check k given response parameters
/// \param world
/// \param Chi
/// \param thresh
/// \param k
void check_k(World &world, X_space &Chi, double thresh, int k);


response_space add_randomness(World &world, const response_space& f, double magnitude);
void normalize(World &world, response_space &f);
void normalize(World &world, X_space &Chi);


static double kronecker(size_t l, size_t n) {
    if (l == n) return 1.0;
    return 0.0;
}

std::map<std::vector<int>, real_function_3d> solid_harmonics(World &world, int n);

/***
 * @brief Prints the norms of the functions of a response space
 *
 *
 *
 * @param world
 * @param f
 */
void print_norms(World &world, const response_space & f) {

    Tensor<double> norms(f.size(), f[0].size());
    // Calc the norms
    long i=0;
    for (const auto & fi:f){
        for(const auto & fij:fi){
            norms(i++) = fij.norm2();
        }
    }
    // Print em in a smart way
    if (world.rank() == 0) print(norms);
}
// Returns a list of solid harmonics such that:
// solid_harm.size() * num_ground_orbs > 2 * num. resp. components
vector_real_function_3d make_xyz_functions(World &world) {
    // Container to return

    // Create the basic x, y, z, constant and zero
    real_function_3d x =
            real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(std::vector<int>{1, 0, 0})));
    real_function_3d y =
            real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 1, 0})));
    real_function_3d z =
            real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 0, 1})));

    std::vector<real_function_3d> funcs={x,y,z};
    return funcs;
}
// Selects from a list of functions and energies the k functions with the
// lowest energy
response_space select_functions(World &world,
                                       response_space &f,
                                       Tensor<double> &energies,
                                       size_t k,
                                       size_t print_level) {
    // Container for result
    response_space answer;

    // Debugging output
    if (print_level >= 1) {
        if (world.rank() == 0) print("\n   Selecting the", k, "lowest excitation energy components.\n");
    }

    // Get rid of extra functions and save
    // the first k
    while (f.size() > k) f.pop_back();
    answer = f;
    answer.truncate_rf();

    // Get rid of extra energies and save
    // the first k
    energies = energies(Slice(0, k - 1));

    // Basic output
    if (print_level >= 1) {
        if (world.rank() == 0) print("   The selected components have excitation energies:");
        if (world.rank() == 0) print(energies);
    }

    // Done
    return answer;
}
// Sorts the given tensor of eigenvalues and
// response functions
void sort(World &world, Tensor<double> &vals, response_space &f) {
    // Get relevant sizes
    size_t k = vals.size();

    // Copy everything...
    response_space f_copy(f);
    Tensor<double> vals_copy = copy(vals);
    Tensor<double> vals_copy2 = copy(vals);

    // Now sort vals_copy
    std::sort(vals_copy.ptr(), vals_copy.ptr() + vals_copy.size());

    // Now sort the rest of the things, using the sorted energy list
    // to find the correct indices
    for (size_t i = 0; i < k; i++) {
        // Find matching index in sorted vals_copy
        size_t j = 0;
        while (fabs(vals_copy(i) - vals_copy2(j)) > 1e-8 && j < k) j++;

        // Put corresponding function, difference function, value residual and
        // value in the correct place
        f[i] = f_copy[j];
        vals(i) = vals_copy(i);

        // Change the value of vals_copy2[j] to help deal with duplicates?
        vals_copy2(j) = 10000.0;
    }
}

void sort(World &world, Tensor<double> &vals, X_space &f) {
    // Get relevant sizes
    size_t k = vals.size();

    // Copy everything...
    X_space f_copy(f);
    Tensor<double> vals_copy = copy(vals);
    Tensor<double> vals_copy2 = copy(vals);

    // Now sort vals_copy
    std::sort(vals_copy.ptr(), vals_copy.ptr() + vals_copy.size());

    // Now sort the rest of the things, using the sorted energy list
    // to find the correct indices
    for (size_t i = 0; i < k; i++) {
        // Find matching index in sorted vals_copy
        size_t j = 0;
        while (fabs(vals_copy(i) - vals_copy2(j)) > 1e-8 && j < k) j++;

        // Put corresponding function, difference function, value residual and
        // value in the correct place
        f.X[i] = f_copy.X[j];
        f.Y[i] = f_copy.Y[j];

        vals(i) = vals_copy(i);

        // Change the value of vals_copy2[j] to help deal with duplicates?
        vals_copy2(j) = 10000.0;
    }
}
// Specialized for response calculations that returns orthonormalized
// functions
response_space gram_schmidt(World &world, const response_space &f) {
    // Sizes inferred
    size_t m = f.size();

    // Return container
    response_space result = f.copy();

    // Orthogonalize
    for (size_t j = 0; j < m; j++) {
        // Need to normalize the row
        double norm = norm2(world, result[j]);

        // Now scale each entry
        scale(world, result[j], 1.0 / norm);

        // Project out from the rest of the vectors
        for (size_t k = j + 1; k < m; k++) {
            // Temp function to hold the sum
            // of inner products
            // vmra.h function, line 627
            double temp = inner(result[j], result[k]);

            // Now subtract
            gaxpy(world, 1.0, result[k], -temp, result[j]);
        }
    }
    result.truncate_rf();

    // Done
    return result;
}

#endif// MADNESS_RESPONSEBASE_HPP
