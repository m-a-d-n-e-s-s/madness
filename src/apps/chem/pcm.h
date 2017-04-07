/*
 This file is part of MADNESS.

 Copyright (C) 2007,2010 Oak Ridge National Laboratory

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

 For more information please contact:

 Robert J. Harrison
 Oak Ridge National Laboratory
 One Bethel Valley Road
 P.O. Box 2008, MS-6367

 email: harrisonrj@ornl.gov
 tel:   865-241-3937
 fax:   865-572-0680
 */

#ifndef SRC_APPS_CHEM_PCM_H_
#define SRC_APPS_CHEM_PCM_H_

#include <memory>

#include <madness/tensor/tensor.h>
#include <madness/world/vector.h>
#include <madness/mra/functypedefs.h>
#include <chem/molecule.h>

#ifdef MADNESS_HAS_PCM
#include "PCMSolver/pcmsolver.h"
#endif

namespace madness {

/// interface class to the PCMSolver library

/// PCMSolver, an API for the Polarizable Continuum Model electrostatic problem.
/// Main authors: R. Di Remigio, L. Frediani, K. Mozgawa
class PCM {

public:

    /// default ctor
    PCM() : pcm_context(0) {}

    /// constructor for the interface to PCM

    /// @param[in]  world   the world
    /// @param[in]  mol     the molecule (coordinates and charges of the nuclei)
    /// @param[in]  pcm_data    pcm input data as read from the input file
    /// @param[in]  verbose print the PCM header
    PCM(World& world, const Molecule& mol, const std::string pcm_data,
            const bool verbose);

    /// compute the potential induced by the surrounding solvent

    /// @param[in]  coulomb_potential   the (positive) potential caused by the electron density
    /// @param[in]  dynamic     compute the contribution to a response kernel
    /// @return     the pcm potential with correct sign: J - K + V_nuc + V_pcm
    real_function_3d compute_pcm_potential(const real_function_3d& coulomb_potential,
            const bool dynamic=false) const;

    /// compute the PCM energy based on the most recent call of compute_pcm_potential
    double compute_pcm_energy() const;

private:
#ifdef MADNESS_HAS_PCM
    /// the main pcmsolver object
    std::shared_ptr<pcmsolver_context_t> pcm_context;

    /// default input generator
    PCMInput pcmsolver_input() const;
#else
    void* pcm_context;
#endif

    /// compute the molecular electrostatic potential from the nuclei
    Tensor<double> nuclear_mep(int nr_nuclei, const Tensor<double>& charges,
                         const Tensor<double>& coordinates, const int grid_size,
                         const Tensor<double>& grid) const;

    /// symmetry info for PCM, needs to be memory controlled by the PCM class
    Tensor<int> symmetry_info;

    /// molecular coordinates, needs to be memory controlled by the PCM class
    Tensor<double> coordinates;

    /// nuclear charges of the molecule, needs to be memory controlled by the PCM class
    Tensor<double> charges;

    /// file name for the total (nuclear + electronic) molecular electronic potential
    std::string mep_lbl;

    /// file name for the total (nuclear + electronic) apparent surface charge
    std::string asc_lbl;


};

} /* namespace madness */

#endif /* SRC_APPS_CHEM_PCM_H_ */
