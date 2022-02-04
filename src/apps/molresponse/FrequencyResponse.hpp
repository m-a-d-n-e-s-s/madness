//
// Created by adrianhurtado on 2/3/22.
//

#ifndef MADNESS_FREQUENCYRESPONSE_HPP
#define MADNESS_FREQUENCYRESPONSE_HPP
#include "ResponseBase.hpp"


class FrequencyResponse;

using RHS_Generator = std::function<X_space(World&, FrequencyResponse&)>;

class FrequencyResponse : public ResponseBase {

public:
    FrequencyResponse(World& world, const CalcParams& params) : ResponseBase(world, params) {}
    void initialize(World& world) override;

private:
    X_space PQ;
    RHS_Generator rhs_generator;
    void generate_PQ(World& world) { PQ = rhs_generator(world, *this); }
    void iterate(World& world) override;
    void check_k(World& world, double thresh, int k) override {
        ResponseBase::check_k(world, thresh, k);
        ::check_k(world, PQ, thresh, k);
    }
};

response_space vector_to_PQ(World& world, const vector_real_function_3d& p,
                            const vector_real_function_3d& ground_orbitals, double lo) {

    response_space rhs(world, p.size(), ground_orbitals.size());

    reconstruct(world, ground_orbitals);

    QProjector<double, 3> Qhat(world, ground_orbitals);

    std::vector<real_function_3d> orbitals = ground_orbitals;

    auto f = [&](auto property) {
      auto phat_phi = mul(world, property, ground_orbitals, lo);
      truncate(world, phat_phi);
      // rhs[i].truncate_vec();

      // project rhs vectors for state
      phat_phi = Qhat(phat_phi);
      world.gop.fence();
      return phat_phi;
    };
    std::transform(p.begin(), p.end(), rhs.begin(), f);
    return rhs;
}

X_space dipole_PQ(World& world, FrequencyResponse calc) {

    auto [gc,molecule,r_params]=calc.get_parameter();
    vector_real_function_3d dipole_vectors(3);
    size_t i = 0;
    for (auto& d: dipole_vectors) {

        std::vector<int> f(3, 0);
        f[i++] = 1;
        d = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(f)));
    }
    truncate(world, dipole_vectors, true);
    X_space PQ;
    PQ.X=vector_to_PQ(world,dipole_vectors,calc.get_orbitals(),r_params.lo());
    PQ.Y=PQ.X;
    return PQ;
}


X_space nuclear_derivative_PQ(World& world, FrequencyResponse calc) {
    auto [gc,molecule,r_params]=calc.get_parameter();
    auto num_operators = size_t(molecule.natom() * 3);
    auto nuclear_vector = vecfuncT(num_operators);

    for (size_t atom = 0; atom < molecule.natom(); ++atom) {
        for (size_t axis = 0; axis < 3; ++axis) {
            FunctorT func(new MolecularDerivativeFunctor(molecule, atom, axis));
            nuclear_vector.at(atom * 3 + axis) = FunctionT(
                    FactoryT(world).functor(func).nofence().truncate_on_project().truncate_mode(0));
        }
    }

    X_space PQ;
    PQ.X=vector_to_PQ(world,nuclear_vector,calc.get_orbitals(),r_params.lo());
    PQ.Y=PQ.X;
    return PQ;

}



#endif//MADNESS_FREQUENCYRESPONSE_HPP
