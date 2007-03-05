#ifndef TWOSCALE_H
#define TWOSCALE_H

#include <madness_config.h>
#include <tensor/tensor.h>
#include <world/world.h>

namespace madness {
    extern void load_coeffs(World& world);
    extern bool two_scale_coefficients(int k,
                                           Tensor<double>* h0, Tensor<double>* h1,
                                           Tensor<double>* g0, Tensor<double>* g1);
    extern bool two_scale_hg(int k, Tensor<double>* hg);
    extern bool test_two_scale_coefficients();

    extern bool autoc(int k, Tensor<double>* c);
    extern bool test_autoc();
}

#endif
