//
// Created by adrianhurtado on 2/3/22.
//

#ifndef MADNESS_FREQUENCYRESPONSE_HPP
#define MADNESS_FREQUENCYRESPONSE_HPP
#include "ResponseBase.hpp"

using path = std::filesystem::path;

class FrequencyResponse;


using RHS_Generator = std::function<X_space(World &, ResponseBase &)>;

response_space vector_to_PQ(World &world, const vector_real_function_3d &rhs_operators, const vector_real_function_3d &ground_orbitals);
X_space nuclear_generator(World &world, ResponseBase &calc);
X_space dipole_generator(World &world, ResponseBase &calc);
// using RHS_Generator = std::function<X_space(World&, FrequencyResponse&)>;

// Create a quadratic response class

class QuadraticResponse : public ResponseBase
{

    // A quadratic response class needs X_space vectors and one ground state.
    // It will compute the 3rd order response property at the 3 given frequencies
    // Beta(omegaA;omegaB,OmegaC)=tr(xA,vBC)+tr(muA,pBC)+tr(muA,qBC)
    // Where xA, xB, xC are the response functions at the frequencies omegaA, omegaB, omegaC
    // And 2nd order right hand side perturbation vector vBC(xB,xC) and
    // pBC and qBC are the homogeneous components of the 2nd order density matrix response
    // made entirely from first order vectors xB, xC


public:
    QuadraticResponse(World &world, const CalcParams &params, RHS_Generator rhs_generator) : ResponseBase(world, params), generator(std::move(rhs_generator))
    {
        FunctionDefaults<3>::set_cubic_cell(-r_params.L(), r_params.L());
        FunctionDefaults<3>::set_truncate_mode(1);
        auto thresh = r_params.protocol()[0];

        int k;
        // Allow for imprecise conversion of threshold
        if (thresh >= 0.9e-2)
            k = 4;
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
        if (r_params.k() == -1)
        {
            FunctionDefaults<3>::set_k(k);
        }
        else
        {
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
        gradop = gradient_operator<double, 3>(world);
        potential_manager = std::make_shared<PotentialManager>(molecule, "a");
        potential_manager->make_nuclear_potential(world);
        // GaussianConvolution1DCache<double>::map.clear();//(TODO:molresponse-What
        // Create the masking function
        mask = real_function_3d(real_factory_3d(world).f(mask3).initial_level(4).norefine());

        ground_density = make_ground_density(world);
        ground_density.truncate(FunctionDefaults<3>::get_thresh());
        // Basic print
        if (world.rank() == 0)
        {
            print("\nSolving NDIM=", 3, " with thresh", thresh, "    k", FunctionDefaults<3>::get_k(), "  dconv", std::max(thresh, r_params.dconv()), "\n");
        }
        if (world.rank() == 0)
        {
            print("Successfully set protocol");
        }

        // ground state orbitals change
        bool redo = false;
        // Verify ground state orbitals have correct k
        if (FunctionDefaults<3>::get_k() != ground_orbitals[0].k())
        {
            // Re-read orbitals from the archive (assuming
            // the archive has orbitals stored at a higher
            if (world.rank() == 0)
            {
                print("check k: ground orbitals");
            }
            // k value than what was previously computed
            ground_calc.read(world);
            if (world.rank() == 0)
            {
                print("check k: read ground orbitals");
            }
            // k value than what was previously computed
            reconstruct(world, ground_orbitals);
            if (world.rank() == 0)
            {
                print("check k: reconstruct ground orbitals");
            }
            // Reset correct k (its set in g_params.read)
            FunctionDefaults<3>::set_k(k);
            // Project each ground state to correct k
            for (auto &orbital : ground_orbitals)
            {
                orbital = project(orbital, FunctionDefaults<3>::get_k(), thresh, false);
            }
            world.gop.fence();
            if (world.rank() == 0)
            {
                print("check k: project ground orbitals");
            }
            // Clean up a bit
            truncate(world, ground_orbitals);
            if (world.rank() == 0)
            {
                print("check k: truncate ground orbitals");
            }
            // Now that ground orbitals have correct k lets make the ground density
            // again
            ground_density = make_ground_density(world);
            ground_density.truncate(FunctionDefaults<3>::get_thresh());
            if (world.rank() == 0)
            {
                print("check k: make ground density");
            }
            // Ground state orbitals changed, clear old hamiltonian
            redo = true;
        }
        // Recalculate ground state hamiltonian here
        if (redo or !hamiltonian.has_data())
        {
            if (world.rank() == 0)
            {
                print("check k: re-do hamiltonian");
            }
            auto [HAM, HAM_NO_DIAG] = ComputeHamiltonianPair(world);
            if (world.rank() == 0)
            {
                print("check k: output hamiltonian");
            }
            // TODO this doesn't seem right...
            hamiltonian = HAM;
            ham_no_diag = HAM_NO_DIAG;
        }

        // If we stored the potential, check that too
        if (r_params.store_potential())
        {
            if (FunctionDefaults<3>::get_k() != stored_potential[0][0].k())
            {
                // Project the potential into correct k
                for (auto &potential_vector : stored_potential)
                {
                    reconstruct(world, potential_vector);
                    for (auto &vi : potential_vector)
                    {
                        vi = project(vi, FunctionDefaults<3>::get_k(), thresh, false);
                    }
                    world.gop.fence();
                }
            }
            if (FunctionDefaults<3>::get_k() != stored_v_coul.k())
                stored_v_coul = project(stored_v_coul, FunctionDefaults<3>::get_k(), thresh, false);
            if (FunctionDefaults<3>::get_k() != stored_v_nuc.k())
                stored_v_nuc = project(stored_v_nuc, FunctionDefaults<3>::get_k(), thresh, false);
        }
        // project the mask
        if (FunctionDefaults<3>::get_k() != mask.k())
        {
            mask = project(mask, FunctionDefaults<3>::get_k(), thresh, false);
            if (world.rank() == 0)
            {
                print("check k: project mask");
            }
        }
        if (world.rank() == 0)
        {
            print("check k: project Chi");
        }
        // Make sure everything is done before leaving
        world.gop.fence();
    }


    void set_x_data(World &world, const std::array<double, 3> &freqABC, const std::array<path, 3> &restart_files)
    {
        this->frequencies = freqABC;
        if (freqABC.size() != 3)
        {
            throw std::runtime_error("Quadratic response requires 3 freqABC");
        }
        // print ABC frequencies
        print("ABC frequencies: ", freqABC[0], " ", freqABC[1], " ", freqABC[2]);

        for (size_t i = 0; i < freqABC.size(); i++)
        {
            auto omega = freqABC[i];

            if (omega == 0.0)
            {
                frequency_contexts[i].set_strategy(std::make_unique<static_inner_product>(), std::make_unique<J1StrategyStable>(), std::make_unique<K1StrategyStatic>(),
                                                   std::make_unique<VXC1StrategyStandard>(), std::make_unique<StaticDensityStrategy>(), std::make_unique<LoadFrequencyXSpace>(), r_params);
            }
            else
            {
                frequency_contexts[i].set_strategy(std::make_unique<full_inner_product>(), std::make_unique<J1StrategyStable>(), std::make_unique<K1StrategyFull>(),
                                                   std::make_unique<VXC1StrategyStandard>(), std::make_unique<FullDensityStrategy>(), std::make_unique<LoadFrequencyXSpace>(), r_params);
            }

            x_data[i] = frequency_contexts[i].load_x_space(world, restart_files[i], r_params, omega);
            ::check_k(world, x_data[i].first, FunctionDefaults<3>::get_thresh(), r_params.k());
        }
    };

    void initialize(World &world) override;

    void load(World &world, const std::string &name) override;

    void save(World &world, const std::string &name) override;
    void iterate(World &world) override;

    Tensor<double> compute_beta(World &world);
    std::pair<Tensor<double>, std::vector<std::string>> compute_beta_v2(World &world, const double &omega_b, const double &omega_c);


private:
    std::vector<int> index_A;
    std::vector<int> index_B;
    std::vector<int> index_C;
    std::vector<std::string> bc_directions;
    std::vector<std::string> a{"X", "Y", "Z"};
    bool indicies_set;

    std::map<int, std::string> xyz = {{0, "X"}, {1, "Y"}, {2, "Z"}};


    std::array<Context, 3> frequency_contexts;
    std::array<double, 3> frequencies;
    std::array<XData, 3> x_data;
    std::pair<X_space, X_space> setup_XBC(World &world, const double &omega_b, const double &omega_c);
    RHS_Generator generator;
    std::pair<X_space, X_space> dipole_perturbation(World &world, const X_space &left, const X_space &right) const;
    X_space compute_g1_term(World &world, const X_space &left, const X_space &right, const X_space &apply) const;
    X_space compute_coulomb_term(World &world, const X_space &A, const X_space &B, const X_space &D) const;
    X_space compute_exchange_term(World &world, const X_space &A, const X_space &B, const X_space &x_apply) const;
    std::tuple<X_space, X_space, X_space, X_space> compute_zeta_response_vectors(World &world, const X_space &B, const X_space &C);
    std::pair<X_space, X_space> compute_first_order_fock_matrix_terms(World &world, const X_space &A, const X_space &phi0, const X_space &B) const;
    std::pair<X_space, X_space> compute_first_order_fock_matrix_terms_v2(World &world, const X_space &B, const X_space &C, const X_space &g1b, const X_space &g1c, const X_space &VB, const X_space &VC,
                                                                         const X_space &phi0) const;
    std::pair<Tensor<double>, std::vector<std::string>> compute_beta_tensor(World &world, const X_space &AB_left, const X_space &AB_right, const X_space &BA_left, const X_space &BA_right,
                                                                            const X_space &XA, const X_space &VBC);
    X_space compute_second_order_perturbation_terms(World &world, const X_space &B, const X_space &C, const X_space &zeta_bc_x, const X_space &zeta_bc_y, const X_space &zeta_cb_x,
                                                    const X_space &zeta_cb_y, const X_space &phi0);
    X_space compute_second_order_perturbation_terms_v2(World &world, const X_space &B, const X_space &C, const X_space &zeta_bc_x, const X_space &zeta_bc_y, const X_space &zeta_cb_x,
                                                       const X_space &zeta_cb_y, const X_space &phi0);
    std::tuple<X_space, X_space, X_space, X_space, X_space, X_space> compute_beta_exchange(World &world, const X_space &B, const X_space &C, const X_space &zeta_bc_left, const X_space &zeta_bc_right,
                                                                                           const X_space &zeta_cb_left, const X_space &zeta_cb_right, const X_space &phi0);
    std::tuple<X_space, X_space, X_space, X_space, X_space, X_space> compute_beta_coulomb(World &world, const X_space &B, const X_space &C, const X_space &zeta_bc_left, const X_space &zeta_bc_right,
                                                                                          const X_space &zeta_cb_left, const X_space &zeta_cb_right, const X_space &phi0);
};


// Create a quadratic response class

class PODResponse : public ResponseBase
{

    // A quadratic response class needs X_space vectors and one ground state.
    // It will compute the 3rd order response property at the 3 given frequencies
    // Beta(omegaA;omegaB,OmegaC)=tr(xA,vBC)+tr(muA,pBC)+tr(muA,qBC)
    // Where xA, xB, xC are the response functions at the frequencies omegaA, omegaB, omegaC
    // And 2nd order right hand side perturbation vector vBC(xB,xC) and
    // pBC and qBC are the homogeneous components of the 2nd order density matrix response
    // made entirely from first order vectors xB, xC
    // vBC(xB,xC) = Q*(
public:
    PODResponse(World &world, const CalcParams &params, RHS_Generator rhs_generator) : ResponseBase(world, params), generator(std::move(rhs_generator))
    {
        FunctionDefaults<3>::set_cubic_cell(-r_params.L(), r_params.L());
        FunctionDefaults<3>::set_truncate_mode(1);
        auto thresh = r_params.protocol()[0];

        int k;
        // Allow for imprecise conversion of threshold
        if (thresh >= 0.9e-2)
            k = 4;
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
        if (r_params.k() == -1)
        {
            FunctionDefaults<3>::set_k(k);
        }
        else
        {
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
        gradop = gradient_operator<double, 3>(world);
        potential_manager = std::make_shared<PotentialManager>(molecule, "a");
        potential_manager->make_nuclear_potential(world);
        // GaussianConvolution1DCache<double>::map.clear();//(TODO:molresponse-What
        // Create the masking function
        mask = real_function_3d(real_factory_3d(world).f(mask3).initial_level(4).norefine());

        ground_density = make_ground_density(world);
        ground_density.truncate(FunctionDefaults<3>::get_thresh());
        // Basic print
        if (world.rank() == 0)
        {
            print("\nSolving NDIM=", 3, " with thresh", thresh, "    k", FunctionDefaults<3>::get_k(), "  dconv", std::max(thresh, r_params.dconv()), "\n");
        }
        if (world.rank() == 0)
        {
            print("Successfully set protocol");
        }

        // ground state orbitals change
        bool redo = false;
        // Verify ground state orbitals have correct k
        if (FunctionDefaults<3>::get_k() != ground_orbitals[0].k())
        {
            // Re-read orbitals from the archive (assuming
            // the archive has orbitals stored at a higher
            if (world.rank() == 0)
            {
                print("check k: ground orbitals");
            }
            // k value than what was previously computed
            ground_calc.read(world);
            if (world.rank() == 0)
            {
                print("check k: read ground orbitals");
            }
            // k value than what was previously computed
            reconstruct(world, ground_orbitals);
            if (world.rank() == 0)
            {
                print("check k: reconstruct ground orbitals");
            }
            // Reset correct k (its set in g_params.read)
            FunctionDefaults<3>::set_k(k);
            // Project each ground state to correct k
            for (auto &orbital : ground_orbitals)
            {
                orbital = project(orbital, FunctionDefaults<3>::get_k(), thresh, false);
            }
            world.gop.fence();
            if (world.rank() == 0)
            {
                print("check k: project ground orbitals");
            }
            // Clean up a bit
            truncate(world, ground_orbitals);
            if (world.rank() == 0)
            {
                print("check k: truncate ground orbitals");
            }
            // Now that ground orbitals have correct k lets make the ground density
            // again
            ground_density = make_ground_density(world);
            ground_density.truncate(FunctionDefaults<3>::get_thresh());
            if (world.rank() == 0)
            {
                print("check k: make ground density");
            }
            // Ground state orbitals changed, clear old hamiltonian
            redo = true;
        }
        // Recalculate ground state hamiltonian here
        if (redo or !hamiltonian.has_data())
        {
            if (world.rank() == 0)
            {
                print("check k: re-do hamiltonian");
            }
            auto [HAM, HAM_NO_DIAG] = ComputeHamiltonianPair(world);
            if (world.rank() == 0)
            {
                print("check k: output hamiltonian");
            }
            // TODO this doesn't seem right...
            hamiltonian = HAM;
            ham_no_diag = HAM_NO_DIAG;
        }

        // If we stored the potential, check that too
        if (r_params.store_potential())
        {
            if (FunctionDefaults<3>::get_k() != stored_potential[0][0].k())
            {
                // Project the potential into correct k
                for (auto &potential_vector : stored_potential)
                {
                    reconstruct(world, potential_vector);
                    for (auto &vi : potential_vector)
                    {
                        vi = project(vi, FunctionDefaults<3>::get_k(), thresh, false);
                    }
                    world.gop.fence();
                }
            }
            if (FunctionDefaults<3>::get_k() != stored_v_coul.k())
                stored_v_coul = project(stored_v_coul, FunctionDefaults<3>::get_k(), thresh, false);
            if (FunctionDefaults<3>::get_k() != stored_v_nuc.k())
                stored_v_nuc = project(stored_v_nuc, FunctionDefaults<3>::get_k(), thresh, false);
        }
        // project the mask
        if (FunctionDefaults<3>::get_k() != mask.k())
        {
            mask = project(mask, FunctionDefaults<3>::get_k(), thresh, false);
            if (world.rank() == 0)
            {
                print("check k: project mask");
            }
        }
        if (world.rank() == 0)
        {
            print("check k: project Chi");
        }
        // Make sure everything is done before leaving
        world.gop.fence();
    }


    void set_x_data(World &world, const std::vector<double> &freqABC, const std::vector<path> &restart_files)
    {
        this->frequencies = freqABC;

        frequency_contexts.resize(freqABC.size());
        x_data.resize(freqABC.size());
        for (size_t i = 0; i < freqABC.size(); i++)
        {
            auto omega = freqABC[i];

            if (omega == 0.0)
            {
                frequency_contexts[i].set_strategy(std::make_unique<static_inner_product>(), std::make_unique<J1StrategyStable>(), std::make_unique<K1StrategyStatic>(),
                                                   std::make_unique<VXC1StrategyStandard>(), std::make_unique<StaticDensityStrategy>(), std::make_unique<LoadFrequencyXSpace>(), r_params);
            }
            else
            {
                frequency_contexts[i].set_strategy(std::make_unique<full_inner_product>(), std::make_unique<J1StrategyStable>(), std::make_unique<K1StrategyFull>(),
                                                   std::make_unique<VXC1StrategyStandard>(), std::make_unique<FullDensityStrategy>(), std::make_unique<LoadFrequencyXSpace>(), r_params);
            }

            // print omega before load_x_space
            print("omega: ", omega);
            print("restart_files[i]: ", restart_files[i]);

            x_data[i] = frequency_contexts[i].load_x_space(world, restart_files[i].string(), r_params, omega);
            ::check_k(world, x_data[i].first, FunctionDefaults<3>::get_thresh(), r_params.k());
        }
    };

    void initialize(World &world) override;

    void load(World &world, const std::string &name) override;

    void save(World &world, const std::string &name) override;
    void iterate(World &world) override;

    void compute_pod_modes(World &world);


    void compute_pod_modes_2(World &world);

private:
    std::vector<Context> frequency_contexts;
    std::vector<double> frequencies;
    std::vector<XData> x_data;
    RHS_Generator generator;
};


class FrequencyResponse : public ResponseBase
{

public:
    FrequencyResponse(World &world, const CalcParams &params, double frequency, RHS_Generator rhs) : ResponseBase(world, params), omega{frequency}, generator{std::move(rhs)}, PQ{}
    {
        if (omega == 0.0)
        {
            response_context.set_strategy(std::make_unique<static_inner_product>(), std::make_unique<J1StrategyStable>(), std::make_unique<K1StrategyStatic>(),
                                          std::make_unique<VXC1StrategyStandard>(), std::make_unique<StaticDensityStrategy>(), std::make_unique<LoadFrequencyXSpace>(), r_params);
        }
        else
        {
            response_context.set_strategy(std::make_unique<full_inner_product>(), std::make_unique<J1StrategyStable>(), std::make_unique<K1StrategyFull>(), std::make_unique<VXC1StrategyStandard>(),
                                          std::make_unique<FullDensityStrategy>(), std::make_unique<LoadFrequencyXSpace>(), r_params);
        }
        PQ = generator(world, *this);
    }
    void initialize(World &world) override;

    void load(World &world, const std::string &name) override;

    void check_k(World &world, double thresh, int k) override
    {
        ResponseBase::check_k(world, thresh, k);
        ::check_k(world, PQ, thresh, k);
    }

    X_space PQ;

    RHS_Generator generator;

private:
    double omega;
    void iterate(World &world) override;
    X_space bsh_update_response(World &world, X_space &theta_X, vector<poperatorT> &bsh_x_ops, vector<poperatorT> &bsh_y_ops, QProjector<double, 3> &projector, double &x_shifts);
    static void frequency_to_json(json &j_mol_in, size_t iter, const Tensor<double> &polar_ij, const Tensor<double> &res_polar_ij);
    void save(World &world, const std::string &name) override;
    std::tuple<X_space, residuals, vector_real_function_3d> update_response(World &world, X_space &chi, XCOperator<double, 3> &xc, std::vector<poperatorT> &bsh_x_ops,
                                                                            std::vector<poperatorT> &bsh_y_ops, QProjector<double, 3> &projector, double &x_shifts, double &omega_n,
                                                                            response_solver &kain_x_space, size_t iteration, const double &max_rotation, const vector_real_function_3d &rho_old,
                                                                            const Tensor<double> &old_residuals, const X_space &xres_old);
};


#endif // MADNESS_FREQUENCYRESPONSE_HPP
