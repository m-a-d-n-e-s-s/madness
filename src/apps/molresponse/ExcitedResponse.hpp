//
// Created by adrianhurtado on 1/28/22.
//

#ifndef MADNESS_EXCITEDRESPONSE_HPP
#define MADNESS_EXCITEDRESPONSE_HPP
#include "ResponseBase.hpp"

struct ExcitedSpace {

    X_space chi;
    X_space l_chi;
};


class ExcitedResponse : public ResponseBase {

public:
    ExcitedResponse(World& world, const CalcParams& params) : ResponseBase(world, params) {}
    void initialize(World& world) override;
    void iterate_trial(World& world, X_space& trial);
    friend class ExcitedTester;

private:
    Tensor<double> omega;
    void iterate(World& world) override;

    X_space make_random_trial(World& world, size_t m) const;
    X_space make_nwchem_trial(World& world, size_t m) const;

    X_space create_trial_functions2(World& world) const;
    X_space create_trial_functions(World& world, size_t k) const;
    X_space create_virtual_ao_guess(World& world) const;


    void deflateTDA(World& world, X_space& Chi, X_space& old_Chi, X_space& Lambda_X,
                    X_space& old_Lambda_X, Tensor<double>& S, Tensor<double> old_S,
                    Tensor<double> old_A, Tensor<double>& omega, size_t& iteration, size_t& m);
    void deflateFull(World& world, X_space& Chi, X_space& old_Chi, X_space& Lambda_X,
                     X_space& old_Lambda_X, Tensor<double>& S, Tensor<double> old_S,
                     Tensor<double> old_A, Tensor<double>& omega, size_t& iteration, size_t& m);
    void augment(World& world, X_space& Chi, X_space& old_Chi, X_space& Lambda_X,
                 X_space& last_Lambda_X, Tensor<double>& S, Tensor<double>& A,
                 Tensor<double>& old_S, Tensor<double>& old_A, size_t print_level);
    void augment_full(World& world, X_space& Chi, X_space& old_Chi, X_space& Lambda_X,
                      X_space& last_Lambda_X, Tensor<double>& S, Tensor<double>& A,
                      Tensor<double>& old_S, Tensor<double>& old_A, size_t print_level);
    void unaugment(World& world, X_space& Chi, X_space& old_Chi, X_space& Lambda_X,
                   X_space& last_Lambda_X, Tensor<double>& omega, Tensor<double>& S_x,
                   Tensor<double>& A_x, Tensor<double>& old_S, Tensor<double>& old_A,
                   size_t num_states, size_t iter, size_t print_level);
    void unaugment_full(World& world, X_space& Chi, X_space& old_Chi, X_space& Lambda_X,
                        X_space& last_Lambda_X, Tensor<double>& omega, Tensor<double>& S_x,
                        Tensor<double>& A_x, Tensor<double>& old_S, Tensor<double>& old_A,
                        size_t num_states, size_t iter, size_t print_level);
    Tensor<double> diagonalizeFullResponseMatrix(World& world, X_space& Chi, X_space& Lambda_X,
                                                 Tensor<double>& omega, Tensor<double>& S,
                                                 Tensor<double>& A, const double thresh,
                                                 size_t print_level);
    Tensor<double> GetFullResponseTransformation(World& world, Tensor<double>& S, Tensor<double>& A,
                                                 Tensor<double>& evals,
                                                 const double thresh_degenerate);
    void deflateGuesses(World& world, X_space& Chi, X_space& Lambda_X, Tensor<double>& S,
                        Tensor<double>& frequencies, size_t& iteration, size_t& m) const;

    Tensor<double> diagonalizeFockMatrix(World& world, X_space& Chi, X_space& Lambda_X,
                                         Tensor<double>& evals, Tensor<double>& A,
                                         Tensor<double>& S, const double thresh) const;
    /// compute the unitary transformation that diagonalizes the fock matrix

    /// @param[in]  world   the world
    /// @param[in]  overlap the overlap matrix of the orbitals
    /// @param[inout]       fock    the fock matrix; diagonal upon exit
    /// @param[out] evals   the orbital energies
    /// @param[in]  thresh_degenerate       threshold for orbitals being
    /// degenerate
    /// @return             the unitary matrix U: U^T F U = evals
    Tensor<double> get_fock_transformation(World& world, Tensor<double>& overlap,
                                           Tensor<double>& fock, Tensor<double>& evals,
                                           const double thresh_degenerate) const;

    // Sorts the given Tensor of energies
    Tensor<int> sort_eigenvalues(World& world, Tensor<double>& vals, Tensor<double>& vecs) const;
    // Construct the Hamiltonian
    // Returns the shift needed to make sure that
    // -2.0 * (ground_state_energy + excited_state_energy)
    // is negative. Please note: The same shift needs to
    // be applied to the potential.
    Tensor<double> create_shift(World& world, const Tensor<double>& ground,
                                const Tensor<double>& omega, std::string xy) const;


    // Returns the given shift applied to the given potential
    response_space apply_shift(World& world, const Tensor<double>& shifts, const response_space& V,
                               const response_space& f);


    // Function to make a vector of BSH operators using ground and excited
    // state energies
    std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> create_bsh_operators(
            World& world, const Tensor<double>& shift, const Tensor<double>& ground,
            const Tensor<double>& omega, const double lo, const double thresh) const;

    static void excited_to_json(json& j_mol_in, size_t iter, const Tensor<double>& omega);

    /**
 * @brief Computes the BSH Update for an excited state calculation.  Passes in
 * omega and computes the necessary shifts in the potential, computes BSH
 * operators and applys BSH operator
 *
 * \f$ \chi^m=-2\hat{G} * \Theta\chi \f$
 *
 * @param world
 * @param theta_X
 * @param projector
 * @param converged
 * @return X_space
 */
    X_space bsh_update_excited(World& world, const Tensor<double>& omega, X_space& theta_X,
                               QProjector<double, 3>& projector);
    void analysis(World& world, const X_space& chi);
    void save(World& world, const std::string& name) override;
    void load(World& world, const std::string& name) override;

    std::pair<Tensor<double>, Tensor<double>> excited_eig(World& world, Tensor<double>& S,
                                                          Tensor<double>& A,
                                                          const double thresh_degenerate);

    std::tuple<X_space, X_space, X_space, X_space> rotate_excited_vectors(
            World& world, const Tensor<double>& U, const X_space& chi, const X_space& l_chi,
            const X_space& v0_chi, const X_space& gamma_chi);
    std::tuple<Tensor<double>, X_space, X_space, X_space, X_space> rotate_excited_space(
            World& world, X_space& chi, X_space& lchi, X_space& v_chi, X_space& gamma_chi);
    std::tuple<Tensor<double>, X_space, X_space, residuals>
    update_response(World &world, X_space &Chi, XCOperator<double, 3> &xc,
                    QProjector<double, 3> &projector, response_solver &kain_x_space,
                    response_matrix &Xvector, response_matrix &Xresidual, size_t iter,
                    const double &maxrotn, const Tensor<double> old_residuals,
                    const X_space &xres_old);
    X_space create_response_guess(World& world) const;
    std::tuple<Tensor<double>, Tensor<double>, Tensor<double>> reduce_subspace(
            World& world, Tensor<double>& S, Tensor<double>& A, const double thresh_degenerate);
};

class ExcitedTester {
private:
public:
    ExcitedTester(World& world, ExcitedResponse& calc, double thresh) {
        print("Setting Function Defaults");
        calc.set_protocol(world, thresh);
        calc.check_k(world, thresh, FunctionDefaults<3>::get_k());
    }
    X_space test_ao_guess(World& world, ExcitedResponse& calc);
};


#endif//MADNESS_EXCITEDRESPONSE_HPP
