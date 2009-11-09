#ifndef MADNESS_QM_PROP_H
#define MADNESS_QM_PROP_H

/// \file qmprop.h
/// \brief Prototypes for qm propagator

namespace madness {
    template <int NDIM>
    SeparatedConvolution<double_complex,NDIM>
    qm_free_particle_propagator(World& world, int k, double bandlimit, double timestep);

    template <int NDIM>
    SeparatedConvolution<double_complex,NDIM>*
    qm_free_particle_propagatorPtr(World& world, int k, double bandlimit, double timestep);
}


#endif
