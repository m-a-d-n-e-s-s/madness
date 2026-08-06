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

/// \file Restart.h
/// \brief the header of a restartdata archive, in one place

#ifndef MADNESS_CHEM_RESTART_H__INCLUDED
#define MADNESS_CHEM_RESTART_H__INCLUDED

#include <madness/chem/molecule.h>
#include <madness/mra/mra.h>

#include <optional>
#include <string>

namespace madness {

/// what the functions in a restartdata archive actually ARE
///
/// moldft stores psi; nemo stores the regularized F = psi/R; znemo the complex
/// variant. Both engines write the same filename with the same version tag, so
/// without this the two are indistinguishable and loading one as the other is
/// silently wrong.
///
/// The numeric values are part of the archive format: they are written as an int.
/// Add new members at the end and NEVER renumber an existing one, or archives
/// already on disk change meaning.
enum class Representation : int {
    unknown = 0,   ///< not recorded; a version-4 archive, which predates the field
    mo      = 1,   ///< moldft orbitals psi
    nemo    = 2,   ///< nemo's regularized orbitals F = psi/R
    znemo   = 3,   ///< complex regularized orbitals
};

inline std::string to_string(const Representation r) {
    switch (r) {
        case Representation::mo:    return "mo";
        case Representation::nemo:  return "nemo";
        case Representation::znemo: return "znemo";
        default:                    return "unknown";
    }
}

/// map a stored int back to a Representation, tolerating values we do not know
inline Representation representation_from_int(const int i) {
    switch (i) {
        case 1: return Representation::mo;
        case 2: return Representation::nemo;
        case 3: return Representation::znemo;
        default: return Representation::unknown;   // includes a newer writer's value
    }
}

/// the header of a `<prefix>.restartdata` archive
///
/// This is the single definition of that header. Before it existed the layout
/// was open-coded in four places -- SCF::save_mos, SCF::load_mos,
/// MolecularOrbitals::{save,read}_restartdata, and two response tools -- which
/// had to be kept in field order by hand.
///
/// The orbitals follow the header in the archive, one block per spin (see
/// MolecularOrbitals::{save,load}_mos). Everything needed to decide whether an
/// archive is usable lives *in the header*, so it can be inspected without
/// reading a single MRA coefficient -- see peek_restartdata().
///
/// ## Versioning
///
/// Version 5 appends to version 4 rather than reordering it, so the v4 layout is
/// a strict prefix of v5. read() therefore dispatches on the stored version and
/// reads the tail only when it is there; v4 archives written by any earlier
/// MADNESS remain loadable, and the fields they lack take the defaults below.
/// write() always emits the current version.
///
/// If you add a field: append it, bump CURRENT_VERSION, and give it a default
/// that means "unknown" so a v4 archive stays meaningful. Do NOT reorder.
struct RestartMetadata {

    static constexpr unsigned int CURRENT_VERSION = 5;

    /// version of the archive this was read from, or CURRENT_VERSION for a fresh one
    unsigned int version = CURRENT_VERSION;

    // --- version 4 fields, in archive order -------------------------------
    double current_energy = 1.e10;
    bool spin_restricted = true;
    double L = 0.0;
    int k = 0;
    Molecule molecule;
    std::string xc;
    std::string localize;
    double converged_for_thresh = 1.e10;

    // --- appended in version 5 --------------------------------------------

    /// density convergence the orbitals were converged to. Not stored by v4, so
    /// a v4 archive reads back as "unconverged" on this axis rather than
    /// claiming a convergence it never recorded.
    double converged_for_dconv = 1.e10;

    /// what the stored functions are; see Representation
    Representation representation = Representation::unknown;

    /// nuclear correlation factor, e.g. "slater:2.0". Meaningful for nemo/znemo
    /// only; an empty string means none was recorded.
    std::string ncf;

    /// molecular smoothing parameter. Changes the nuclear potential, so orbitals
    /// from a different eprec solve a different Hamiltonian and are not the same
    /// solution -- unlike k or thresh, this cannot be repaired by reprojecting.
    double eprec = 0.0;

    /// MADNESS version that wrote the archive, for provenance in bug reports
    std::string madness_version;

    // NB: no truncate_mode, k-per-function, thresh-per-function or tree_state
    // here. FunctionImpl::store/load already serialize those with every orbital
    // (mra/funcimpl.h), so the loaded functions carry them -- which is how
    // SCF::load_mos can query amo[0].k() and amo[0].thresh() after loading.
    // Duplicating them in the header would create two sources of truth in one
    // file. Only what is needed BEFORE any orbital is read belongs here.

    /// write the header at the current archive position, always at CURRENT_VERSION
    template <typename Archive>
    void write(Archive& ar) const {
        const unsigned int v = CURRENT_VERSION;
        ar & v;
        ar & current_energy & spin_restricted;
        ar & L & k & molecule & xc & localize & converged_for_thresh;
        // version 5 tail. The representation goes out as an int; see the warning
        // on Representation about never renumbering.
        const int rep = static_cast<int>(representation);
        ar & converged_for_dconv & rep & ncf & eprec & madness_version;
    }

    /// read the header from the current archive position
    ///
    /// Accepts version 4 and 5. Throws on anything else, since the field order
    /// would be unknown and the archive position left wrong for the orbitals.
    template <typename Archive>
    void read(Archive& ar) {
        ar & version;
        MADNESS_CHECK_THROW(version == 4 or version == 5,
                "unsupported restartdata version: only 4 and 5 can be read");

        ar & current_energy & spin_restricted;
        ar & L & k & molecule & xc & localize & converged_for_thresh;

        if (version >= 5) {
            int rep = 0;
            ar & converged_for_dconv & rep & ncf & eprec & madness_version;
            representation = representation_from_int(rep);
        }
        // else: the v5 members keep their defaults, which all read as "unknown"
    }

    /// true if these orbitals are at least as converged as the request
    bool is_converged_to(const double thresh, const double dconv) const {
        return converged_for_thresh <= thresh and converged_for_dconv <= dconv;
    }

    /// true if the stored functions are what the asking engine expects
    ///
    /// A v4 archive records no representation; treat that as compatible rather
    /// than rejecting every pre-existing archive, and rely on the geometry and
    /// convergence checks instead.
    bool representation_matches(const Representation wanted) const {
        return representation == Representation::unknown or representation == wanted;
    }

    std::string print_to_string() const {
        std::stringstream ss;
        ss << "restartdata v" << version
           << "  representation " << madness::to_string(representation)
           << (ncf.empty() ? std::string() : "  ncf " + ncf)
           << "  k " << k << "  L " << L
           << "  converged to thresh " << converged_for_thresh
           << " dconv " << converged_for_dconv
           << "  energy " << current_energy;
        return ss.str();
    }
};

/// read ONLY the header of a restartdata archive, without loading any orbitals
///
/// This is what makes "is the archive on disk good enough?" answerable cheaply:
/// SCF::load_mos pulls every MRA function off disk before any of the header can
/// be looked at.
///
/// @param[in] world     the world; the parallel archive broadcasts, so all ranks
///                      see the same result and the return value is collective
/// @param[in] filename  archive name WITHOUT the chunk suffix, e.g. "mad.restartdata"
/// @return    the header, or nullopt if the archive is absent or unreadable
inline std::optional<RestartMetadata>
peek_restartdata(World& world, const std::string& filename) {
    try {
        archive::ParallelInputArchive<archive::BinaryFstreamInputArchive>
            ar(world, filename.c_str());
        RestartMetadata meta;
        meta.read(ar);
        return meta;
    } catch (...) {
        // absent, truncated, or a version we cannot parse. The caller decides
        // whether that is benign (nothing to restart from) or alarming (a file
        // exists but cannot be read) -- see restartdata_exists().
        return std::nullopt;
    }
}

/// true if the first chunk of a restartdata archive is present on disk
///
/// Lets a caller tell "no archive" from "archive present but unreadable", which
/// want different reactions: the first is normal, the second is a corrupt file
/// or a format mismatch and should be reported rather than silently recomputed.
inline bool restartdata_exists(const std::string& filename) {
    return std::filesystem::exists(filename + ".00000");
}

} // namespace madness

#endif // MADNESS_CHEM_RESTART_H__INCLUDED
