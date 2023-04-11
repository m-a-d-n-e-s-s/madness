//
// Created by adrianhurtado on 2/3/22.
//

#ifndef MADNESS_FREQUENCYRESPONSE_HPP
#define MADNESS_FREQUENCYRESPONSE_HPP
#include "ResponseBase.hpp"


class FrequencyResponse;


using RHS_Generator = std::function<X_space(World &, FrequencyResponse &)>;

response_space vector_to_PQ(World &world,
                            const vector_real_function_3d &rhs_operators,
                            const vector_real_function_3d &ground_orbitals);
X_space nuclear_generator(World &world, FrequencyResponse &calc);
X_space dipole_generator(World &world, FrequencyResponse &calc);
//using RHS_Generator = std::function<X_space(World&, FrequencyResponse&)>;

class FrequencyResponse : public ResponseBase {

public:
    FrequencyResponse(World &world, const CalcParams &params, double frequency,
                      RHS_Generator rhs)
        : ResponseBase(world, params), omega{frequency},
          generator{std::move(rhs)}, PQ{} {
        if (omega == 0.0) {
            response_context.set_strategy(
                    std::make_unique<static_inner_product>(),
                    std::make_unique<J1StrategyStable>(),
                    std::make_unique<K1StrategyStatic>(),
                    std::make_unique<VXC1StrategyStandard>(),
                    std::make_unique<StaticDensityStrategy>());
        } else {
            response_context.set_strategy(
                    std::make_unique<full_inner_product>(),
                    std::make_unique<J1StrategyStable>(),
                    std::make_unique<K1StrategyFull>(),
                    std::make_unique<VXC1StrategyStandard>(),
                    std::make_unique<FullDensityStrategy>());
        }
        PQ = generator(world, *this);
    }
    void initialize(World &world) override;

    void load(World &world, const std::string &name) override;

    void check_k(World &world, double thresh, int k) override {
        ResponseBase::check_k(world, thresh, k);
        ::check_k(world, PQ, thresh, k);
    }

    X_space PQ;

    RHS_Generator generator;

private:
    double omega;
    void iterate(World &world) override;
    X_space bsh_update_response(World &world, X_space &theta_X,
                                vector<poperatorT> &bsh_x_ops,
                                vector<poperatorT> &bsh_y_ops,
                                QProjector<double, 3> &projector,
                                double &x_shifts);
    static void frequency_to_json(json &j_mol_in, size_t iter,
                                  const Tensor<double> &polar_ij,
                                  const Tensor<double> &res_polar_ij);
    static void compute_and_print_polarizability(World &world, X_space &Chi,
                                                 X_space &pq,
                                                 std::string message);
    void save(World &world, const std::string &name) override;
    std::tuple<X_space, residuals, vector_real_function_3d> update_response(
            World &world, X_space &chi, XCOperator<double, 3> &xc,
            std::vector<poperatorT> &bsh_x_ops,
            std::vector<poperatorT> &bsh_y_ops,
            QProjector<double, 3> &projector, double &x_shifts, double &omega_n,
            response_solver &kain_x_space, size_t iteration,
            const double &max_rotation, const vector_real_function_3d &rho_old,
            const Tensor<double> &old_residuals, const X_space &xres_old);
};


#endif//MADNESS_FREQUENCYRESPONSE_HPP
