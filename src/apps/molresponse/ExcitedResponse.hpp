//
// Created by adrianhurtado on 1/28/22.
//

#ifndef MADNESS_EXCITEDRESPONSE_HPP
#define MADNESS_EXCITEDRESPONSE_HPP
#include "ResponseBase.hpp"


class ExcitedResponse : public ResponseBase {

public:
    ExcitedResponse(World& world, const CalcParams& params) : ResponseBase(world, params) {}
    void initialize(World& world) override;
    void iterate_trial(World& world, X_space& trial);

private:
    Tensor<double> omega;
    void iterate(World& world) override;

    X_space make_random_trial(World& world, size_t m) const;
    X_space make_nwchem_trial(World& world, size_t m) const;

    X_space create_trial_functions2(World& world) const;
    X_space create_trial_functions(World& world, size_t k) const;


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
                        Tensor<double>& omega, size_t& iteration, size_t& m) const;

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

    void excited_to_json(json& j_mol_in, size_t iter, const Tensor<double>& res_X,
                         const Tensor<double>& res_Y, const Tensor<double>& density_res,
                         const Tensor<double>& omega);

    void update_x_space_excited(World& world, X_space& old_Chi, X_space& Chi, X_space& old_Lambda_X,
                                X_space& res, XCOperator<double, 3>& xc,
                                QProjector<double, 3>& projector, Tensor<double>& omega_n,
                                NonLinearXsolver& kain_x_space, std::vector<X_vector>& Xvector,
                                std::vector<X_vector>& Xresidual, Tensor<double>& energy_residuals,
                                Tensor<double>& old_energy, Tensor<double>& bsh_residualsX,
                                Tensor<double>& bsh_residualsY, Tensor<double>& S,
                                Tensor<double>& old_S, Tensor<double>& A, Tensor<double>& old_A,
                                size_t iter, Tensor<double>& maxrotn);

    // Load Balancing
    void compute_new_omegas_transform(World& world, X_space& old_Chi, X_space& Chi,
                                      X_space& old_Lambda_X, X_space& Lambda_X,
                                      Tensor<double>& omega, Tensor<double>& old_energy,
                                      Tensor<double>& S, Tensor<double>& old_S, Tensor<double>& A,
                                      Tensor<double>& old_A, Tensor<double>& energy_residuals,
                                      size_t iter);

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
};


#endif//MADNESS_EXCITEDRESPONSE_HPP
