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

#include <map>
#include <memory>
#include <utility>

#include <madness/tensor/tensor.h>
#include <madness/world/vector.h>
#include <madness/mra/functypedefs.h>
#include <madness/mra/QCCalculationParametersBase.h>
#include <madness/mra/commandlineparser.h>
#include<madness/chem/molecule.h>
#include<madness/chem/CalculationParameters.h>

#ifdef MADNESS_HAS_PCM
#include <PCMSolver/pcmsolver.h>
#include <PCMSolver/PCMInput.h>
#endif

namespace madness {

/// input parameters for the polarizable continuum model -- the `pcm` data group

/// Declared outside MADNESS_HAS_PCM on purpose: ParameterManager and madqc hold and
/// print this group whether or not PCMSolver was compiled in, and only the conversion
/// to PCMSolver's PCMInput needs the library.
///
/// The group refines `pcm_data` of the `dft` group, which stays the on/off switch:
/// `pcm_data <solvent>` on its own is a complete specification, and the solvent it
/// names arrives here as a *derived* value that an explicit `solvent` in this group
/// overrides.
///
/// Only the thirteen PCMInput fields that PCMSolver's host-struct reader actually
/// consults are exposed. `patch_level`, `coarsity`, `min_distance`, `der_order` and
/// `equation_type` are accepted by the struct but never read by
/// Input::reader(const PCMInput&), and upstream documents no keyword for them.
class PCMParameters : public QCCalculationParametersBase {
public:
    static constexpr char const* tag = "pcm";

    PCMParameters();

    PCMParameters(World& world, const commandlineparser& parser) : PCMParameters() {
        read_input_and_commandline_options(world, parser, tag);
    }

    std::string get_tag() const override { return std::string(tag); }

    /// adopt the solvent named by the dft group's `pcm_data`, derive from it, validate

    /// The solvent is set as a *derived* value, so an explicit `solvent` in the `pcm`
    /// group wins; `epsilon` and `probe_radius` are then derived from the solvent in
    /// turn, and only stand as given for `solvent explicit`. Validation that
    /// allowed_values cannot express -- buffer capacities and cross-key consistency --
    /// happens here, as OEP_Parameters does for its density thresholds.
    void set_derived_values(const CalculationParameters& cparam);

    std::string solvent() const { return get<std::string>("solvent"); }
    std::string solver_type() const { return get<std::string>("solver_type"); }
    double correction() const { return get<double>("correction"); }
    std::string cavity_type() const { return get<std::string>("cavity_type"); }
    double area() const { return get<double>("area"); }
    std::string radii_set() const { return get<std::string>("radii_set"); }
    bool scaling() const { return get<bool>("scaling"); }
    double min_radius() const { return get<double>("min_radius"); }
    std::string restart_name() const { return get<std::string>("restart_name"); }
    double probe_radius() const { return get<double>("probe_radius"); }
    double epsilon() const { return get<double>("epsilon"); }

    /// whether PCMSolver should read its own input file instead of our parameters

    /// Set from `pcm_data pcmsolver_reader_own` in the dft group; not a key of this
    /// group, because in that mode none of the keys below are consulted.
    bool reader_own() const { return reader_own_; }

    /// the solvents PCMSolver tabulates

    /// Maps the canonical name to {static permittivity, probe radius / Angstrom};
    /// values from PCMSolver's src/utils/Solvent.cpp. Two uses: deriving `epsilon`
    /// and `probe_radius` from `solvent`, and spelling out the `solvent explicit`
    /// equivalent for the names that do not fit PCMInput::solvent.
    static const std::map<std::string, std::pair<double, double> >& solvent_data();

    /// longest solvent name PCMInput::solvent (char[16]) can hold with its NUL
    static constexpr std::size_t max_solvent_name_length = 15;

#ifdef MADNESS_HAS_PCM
    /// fill PCMSolver's host input struct from these parameters
    PCMInput to_pcmsolver_input() const;
#endif

private:
    bool reader_own_ = false;
};

/// interface class to the PCMSolver library

/// PCMSolver, an API for the Polarizable Continuum Model electrostatic problem.
/// Main authors: R. Di Remigio, L. Frediani, K. Mozgawa
class PCM {

public:

    /// default ctor
    PCM() = default;

    /// constructor for the interface to PCM

    /// @param[in]  world   the world
    /// @param[in]  mol     the molecule (coordinates and charges of the nuclei)
    /// @param[in]  param   the `pcm` data group, already reconciled with `pcm_data`
    /// @param[in]  verbose print the PCM header and echo the pcm parameter group
    PCM(World& world, const Molecule& mol, const PCMParameters& param,
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
