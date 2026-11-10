/// \file magnetic_kernels.h
/// \brief All magnetic Green's-function kernels and the factory; see magnetic_kernels.md

#ifndef MADNESS_MRA_MAGNETIC_KERNELS_H__INCLUDED
#define MADNESS_MRA_MAGNETIC_KERNELS_H__INCLUDED

#include <madness/mra/magnetic_kernel.h>
#include <madness/mra/magnetic_kernel_chain.h>
#include <madness/mra/magnetic_kernel_ft.h>
#include <madness/mra/magnetic_kernel_mehler.h>

namespace madness {

/// the names accepted by make_magnetic_kernel
inline std::vector<std::string> magnetic_kernel_names() { return {"ft", "chain_symmetric", "chain_landau", "mehler"}; }

/// kernel of kind p.kind; requires FunctionDefaults<3> (cubic cell, k, thresh) to be set
inline std::shared_ptr<MagneticKernel> make_magnetic_kernel(World& world, const MagneticKernelParameters& p) {
    if (p.kind == "ft") return std::make_shared<MagneticKernelFT>(world, p);
    if (p.kind == "chain_symmetric") return std::make_shared<MagneticKernelChain>(world, p, MagneticGauge::symmetric);
    if (p.kind == "chain_landau") return std::make_shared<MagneticKernelChain>(world, p, MagneticGauge::landau);
    if (p.kind == "mehler") return std::make_shared<MagneticKernelMehler>(world, p);
    MADNESS_EXCEPTION(("make_magnetic_kernel: unknown kernel " + p.kind).c_str(), 1);
}

}  // namespace madness

#endif  // MADNESS_MRA_MAGNETIC_KERNELS_H__INCLUDED
