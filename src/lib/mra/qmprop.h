#ifndef MADNESS_MRA_QMPROP_H__INCLUDED
#define MADNESS_MRA_QMPROP_H__INCLUDED

/// \file qmprop.h
/// \brief Prototypes for qm propagator

namespace madness {
    Convolution1D<double_complex>*
    qm_1d_free_particle_propagator(int k, double bandlimit, double timestep, double width);

    template <int NDIM>
    SeparatedConvolution<double_complex,NDIM>
    qm_free_particle_propagator(World& world, int k, double bandlimit, double timestep);

    template <int NDIM>
    SeparatedConvolution<double_complex,NDIM>*
    qm_free_particle_propagatorPtr(World& world, int k, double bandlimit, double timestep);
}


#endif
