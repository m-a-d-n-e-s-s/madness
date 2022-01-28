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

private:
    X_space make_random_trial(World& world, size_t m) const;
    X_space make_nwchem_trial(World& world, size_t m) const;

    X_space create_trial_functions2(World& world);
    X_space create_trial_functions(World& world, size_t k);

    void iterate_trial(World & world, X_space & trial)const ;


};


#endif//MADNESS_EXCITEDRESPONSE_HPP
