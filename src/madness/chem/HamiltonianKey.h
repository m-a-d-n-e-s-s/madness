/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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

/// \file HamiltonianKey.h
/// \brief what operator a set of orbitals solves, in one place

#ifndef MADNESS_CHEM_HAMILTONIANKEY_H__INCLUDED
#define MADNESS_CHEM_HAMILTONIANKEY_H__INCLUDED

#include <madness/mra/mra.h>  // before pcm.h, which is not self-contained
#include <madness/chem/CalculationParameters.h>
#include <madness/chem/molecule.h>
#include <madness/chem/pcm.h>
#include <nlohmann/json.hpp>

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace madness {

/// everything besides the geometry that changes the operator the SCF solves
///
/// Orbitals converged for one key are a good guess for another, but their
/// convergence claim is about a different problem -- and unlike a k or thresh
/// mismatch that cannot be repaired by reprojecting. The restart planner, the
/// archive header and the madqc results file all compare through this one
/// struct, so they agree on what "the same Hamiltonian" means.
///
/// Every field has a value meaning "not recorded" (empty string, 0, empty
/// vector, -1): an old archive or a seeding tool that did not know it. Not
/// recorded on either side is never a mismatch.
///
/// Deliberately absent:
///  - geometry: compared separately, with a tolerance (compare_geometry)
///  - localize: selects which orbitals of the same solution are stored, not the
///    operator; it lives in the archive header next to the key
///  - dispersion: shifts the energy without touching the orbitals, so it is a
///    property of the results, not of the archive
struct HamiltonianKey {
    std::string xc;               ///< exchange-correlation functional
    double eprec = 0.0;           ///< molecular smoothing parameter
    std::string ncf;              ///< nuclear correlation factor, "type:a"; nemo only
    std::vector<double> field;    ///< external electric field
    std::string core_type;        ///< "none" or "mcp"
    int psp_calc = -1;            ///< pseudopotentials on all atoms: 0/1, -1 unknown
    std::string pcm;              ///< "none", or the solvent model, see make_hamiltonian_key

    /// why \p requested is a different operator, or an empty string if it is not
    ///
    /// The phrase goes into the restart plan's `why`, so it reads as a clause:
    /// "archive was written at eprec 1e-04, this run uses 1e-06".
    std::string mismatch(const HamiltonianKey& requested) const {
        auto fmt = [](const double d) {
            std::stringstream ss;
            ss << std::scientific << std::setprecision(0) << d;
            return ss.str();
        };
        auto differ = [](const std::string& a, const std::string& b) {
            return not a.empty() and not b.empty() and a != b;
        };
        if (eprec != 0.0 and requested.eprec != 0.0 and
                std::abs(eprec / requested.eprec - 1.0) > 1.e-10)
            return "archive was written at eprec " + fmt(eprec) + ", this run uses " +
                   fmt(requested.eprec);
        if (differ(xc, requested.xc))
            return "archive was written with xc '" + xc + "', this run uses '" +
                   requested.xc + "'";
        if (differ(ncf, requested.ncf))
            return "archive was written with the nuclear correlation factor '" + ncf +
                   "', this run uses '" + requested.ncf + "'";
        if (not field.empty() and not requested.field.empty()) {
            bool same = field.size() == requested.field.size();
            for (std::size_t i = 0; same and i < field.size(); ++i)
                same = std::abs(field[i] - requested.field[i]) <= 1.e-12;
            if (not same) return "archive was written in a different external field";
        }
        if (differ(core_type, requested.core_type))
            return "archive was written with core_type '" + core_type +
                   "', this run uses '" + requested.core_type + "'";
        if (psp_calc >= 0 and requested.psp_calc >= 0 and psp_calc != requested.psp_calc)
            return std::string("archive was written ") + (psp_calc ? "with" : "without") +
                   " pseudopotentials, this run " + (requested.psp_calc ? "uses" : "does not use") +
                   " them";
        if (differ(pcm, requested.pcm))
            return "archive was written with solvent '" + pcm + "', this run uses '" +
                   requested.pcm + "'";
        return std::string();
    }

    bool operator==(const HamiltonianKey& other) const {
        return xc == other.xc and eprec == other.eprec and ncf == other.ncf and
               field == other.field and core_type == other.core_type and
               psp_calc == other.psp_calc and pcm == other.pcm;
    }

    nlohmann::json to_json() const {
        return {{"xc", xc},           {"eprec", eprec},         {"ncf", ncf},
                {"field", field},     {"core_type", core_type}, {"psp_calc", psp_calc},
                {"pcm", pcm}};
    }

    static HamiltonianKey from_json(const nlohmann::json& j) {
        HamiltonianKey k;
        k.xc = j.value("xc", std::string());
        k.eprec = j.value("eprec", 0.0);
        k.ncf = j.value("ncf", std::string());
        k.field = j.value("field", std::vector<double>());
        k.core_type = j.value("core_type", std::string());
        k.psp_calc = j.value("psp_calc", -1);
        k.pcm = j.value("pcm", std::string());
        return k;
    }
};

/// the key a calculation with these parameters solves
///
/// @param[in] ncf the nuclear correlation factor as SCF::restart_ncf spells it,
///                empty for moldft
inline HamiltonianKey make_hamiltonian_key(const CalculationParameters& cparam,
                                           const Molecule& molecule,
                                           PCMParameters pcm_param,
                                           const std::string& ncf = "") {
    HamiltonianKey k;
    k.xc = cparam.xc();
    k.eprec = molecule.parameters.eprec();
    k.ncf = ncf;
    k.field = molecule.parameters.field();
    k.core_type = molecule.parameters.core_type();
    k.psp_calc = molecule.parameters.psp_calc() ? 1 : 0;
    // pcm_data is the on/off switch and names the solvent; the pcm group refines
    // it (SCF::SCF, PCMParameters::set_derived_values). In PCMSolver's reader-own
    // mode the model lives in a file we do not parse, so only its name is known.
    const std::string pcm_data = cparam.pcm_data();
    if (pcm_data == "none" or pcm_data == "pcmsolver_reader_own") {
        k.pcm = pcm_data;
    } else {
        // derive the solvent's epsilon/probe the way SCF's ctor does, so a key
        // built from raw input parameters equals the one the engine writes
        pcm_param.set_derived_values(cparam);
        std::stringstream ss;
        ss << pcm_param.solvent() << ":" << pcm_param.solver_type()
           << ":eps=" << pcm_param.epsilon() << ":probe=" << pcm_param.probe_radius();
        k.pcm = ss.str();
    }
    return k;
}

} // namespace madness

#endif // MADNESS_CHEM_HAMILTONIANKEY_H__INCLUDED
