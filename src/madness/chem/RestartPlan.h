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

/// \file RestartPlan.h
/// \brief decide once, per geometry, where the initial orbitals come from

#ifndef MADNESS_CHEM_RESTARTPLAN_H__INCLUDED
#define MADNESS_CHEM_RESTARTPLAN_H__INCLUDED

#include <madness/chem/CalculationParameters.h>
#include <madness/chem/Restart.h>

#include <optional>
#include <string>
#include <vector>

namespace madness {

/// what the user asked for with the `restart` keyval
///
/// `automatic` is spelled "auto" in input files -- `auto` is a C++ keyword.
enum class RestartMode {
    automatic,   ///< look at what is on disk and choose (the default)
    none,        ///< ignore everything on disk, start from the initial guess
    iterate,     ///< read restartdata and keep iterating
    read_only,   ///< read restartdata and do not iterate, whatever its precision
    ao,          ///< read the AO projections (restartaodata) and iterate
    nwchem,      ///< read an NWChem movecs file and iterate
};

/// the six spellings accepted by the `restart` keyval, for allowed_values
inline std::vector<std::string> restart_mode_names() {
    return {"auto", "none", "iterate", "read_only", "ao", "nwchem"};
}

inline std::string to_string(const RestartMode m) {
    switch (m) {
        case RestartMode::automatic: return "auto";
        case RestartMode::none:      return "none";
        case RestartMode::iterate:   return "iterate";
        case RestartMode::read_only: return "read_only";
        case RestartMode::ao:        return "ao";
        case RestartMode::nwchem:    return "nwchem";
    }
    return "auto";
}

/// parse the `restart` keyval
///
/// Throws on an unknown spelling. The keyval carries allowed_values, so
/// QCParameter normally rejects a bad string before this is ever reached; this
/// check is the backstop for values set programmatically.
inline RestartMode restart_mode_from_string(const std::string& s) {
    if (s == "auto")      return RestartMode::automatic;
    if (s == "none")      return RestartMode::none;
    if (s == "iterate")   return RestartMode::iterate;
    if (s == "read_only") return RestartMode::read_only;
    if (s == "ao")        return RestartMode::ao;
    if (s == "nwchem")    return RestartMode::nwchem;
    MADNESS_EXCEPTION("unknown restart mode", 1);
}

/// where the initial orbitals will come from
enum class RestartSource : int {
    initial_guess = 0,   ///< atomic guess, i.e. no restart at all
    restartdata   = 1,   ///< <prefix>.restartdata, MRA orbitals
    restartao     = 2,   ///< <prefix>.restartaodata, AO expansion coefficients
    nwchem        = 3,   ///< an NWChem movecs file
};

inline std::string to_string(const RestartSource s) {
    switch (s) {
        case RestartSource::restartdata: return "restartdata";
        case RestartSource::restartao:   return "restartaodata";
        case RestartSource::nwchem:      return "nwchem";
        default:                         return "initial_guess";
    }
}

/// what a look at the disk found -- pure disk facts, no engine opinion in here
///
/// Separated from the decision so that plan_restart() stays a pure function of
/// its arguments and the whole decision table can be tested without touching a
/// filesystem or constructing an MRA function.
struct RestartSources {

    /// <prefix>.restartdata.00000 exists
    bool restartdata_present = false;

    /// its header parsed. Set only when restartdata_present; a present file with
    /// no metadata is a corrupt or future-version archive, which is worth saying
    /// out loud rather than quietly recomputing.
    std::optional<RestartMetadata> meta;

    /// <prefix>.restartaodata exists
    bool restartao_present = false;

    /// an NWChem file was named in the input
    bool nwfile_named = false;
};

/// which restart sources the asking engine can actually read
///
/// Not every engine can read every source, and pretending otherwise turns a
/// clear "not implemented" into a confusing "file not found". moldft can read
/// all of them; nemo can read neither the AO projections (`SCF::restart_aos` is
/// reachable only from `SCF::get_initial_orbitals`, which nemo never calls, and
/// the file holds <AO|psi> where nemo needs <AO|F>) nor NWChem movecs.
struct RestartCapabilities {
    bool ao = false;       ///< can read <prefix>.restartaodata
    bool nwchem = false;   ///< can read an NWChem movecs file

    static RestartCapabilities all() { return {true, true}; }
    static RestartCapabilities restartdata_only() { return {false, false}; }
};

/// how the requested geometry relates to the one in an archive
enum class GeometryMatch {
    same,                    ///< same atoms at the same places
    displaced,               ///< same atoms, moved -- e.g. a geometry optimization step
    different_composition,   ///< different atoms altogether
};

/// compare two molecules, atom by atom and in order
///
/// Order matters: the orbitals in an archive are expanded about the atoms in
/// the order the archive's molecule lists them, so a permuted molecule is not
/// the same molecule as far as restart is concerned.
inline GeometryMatch compare_geometry(const Molecule& archive, const Molecule& requested,
                                      const double tol = 1.e-8) {
    if (archive.natom() != requested.natom()) return GeometryMatch::different_composition;
    for (std::size_t i = 0; i < requested.natom(); ++i) {
        if (archive.get_atomic_number(i) != requested.get_atomic_number(i))
            return GeometryMatch::different_composition;
    }
    for (std::size_t i = 0; i < requested.natom(); ++i) {
        const Atom a = archive.get_atom(i);
        const Atom b = requested.get_atom(i);
        if (std::abs(a.x - b.x) > tol or std::abs(a.y - b.y) > tol or
            std::abs(a.z - b.z) > tol)
            return GeometryMatch::displaced;
    }
    return GeometryMatch::same;
}

/// the decision: one source, one starting rung, one reason
struct RestartPlan {

    RestartMode mode = RestartMode::automatic;   ///< what was asked for
    RestartSource source = RestartSource::initial_guess;

    /// false means "return what is on disk without solving anything"
    bool iterate = true;

    /// rung of `protocol()` to start at; rungs below it are already covered
    std::size_t protocol_start = 0;

    /// energy from the archive, meaningful only when iterate==false
    double stale_energy = 1.e10;

    /// true if this decision should be logged as a WARNING rather than as info
    ///
    /// Two cases: a file was there but could not be used (silently recomputing
    /// would hide a disk fault or a wrong prefix), and read_only handing back an
    /// energy that does not meet the requested precision.
    bool warn = false;

    /// one line, for the log and for the results json
    std::string why;

    /// true if orbitals have to be read from disk before anything else happens
    bool needs_load() const { return source != RestartSource::initial_guess; }

    std::string print_to_string() const {
        return "restart " + madness::to_string(mode) + ": from " +
               madness::to_string(source) + ", " +
               (iterate ? "starting at protocol rung " + std::to_string(protocol_start)
                        : std::string("no iterations")) +
               " -- " + why;
    }

    /// serialization, so the plan can be decided on one rank and broadcast
    template <typename Archive>
    void serialize(Archive& ar) {
        int m = static_cast<int>(mode);
        int s = static_cast<int>(source);
        ar & m & s & iterate & protocol_start & stale_energy & warn & why;
        mode = static_cast<RestartMode>(m);
        source = static_cast<RestartSource>(s);
    }
};

/// format a threshold compactly for a log line ("1e-06", not "0.000001")
inline std::string format_thresh(const double t) {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(0) << t;
    return ss.str();
}

/// index of the first protocol rung that is tighter than an achieved precision
///
/// nullopt means the ladder holds nothing tighter, i.e. the archive already
/// covers every rung. The 0.999 slack absorbs the round trip of a threshold
/// through an archive; without it a run converged to 1e-6 can fail to recognize
/// the 1e-6 rung as covered.
inline std::optional<std::size_t>
first_rung_tighter_than(const std::vector<double>& protocol, const double achieved) {
    for (std::size_t i = 0; i < protocol.size(); ++i)
        if (protocol[i] < achieved * 0.999) return i;
    return std::nullopt;
}

/// decide where the initial orbitals come from -- the whole decision table
///
/// Pure: no filesystem, no MPI, no MRA. That is the point of the split, and it
/// is what makes the table testable line by line (test_restart.cc).
///
/// @param[in] mode       the user's `restart` keyval
/// @param[in] disk       what survey_restart_sources() found
/// @param[in] can        which sources the asking engine can read
/// @param[in] protocol   the precision ladder, `CalculationParameters::protocol()`
/// @param[in] user_dconv the user's `dconv`
/// @param[in] requested  the geometry this calculation is for
/// @param[in] wanted     the representation the asking engine stores (mo/nemo/znemo)
/// @param[in] eprec      the requested molecular smoothing parameter; 0 skips the
///                       check, which is what a caller that does not know it passes
/// @param[in] xc         the requested exchange-correlation functional; an empty
///                       string skips the check
/// @param[in] ncf        the requested nuclear correlation factor, e.g. "slater:2.0";
///                       empty for an engine that has none, and skips the check
inline RestartPlan plan_restart(const RestartMode mode, const RestartSources& disk,
                                 const RestartCapabilities& can,
                                 const std::vector<double>& protocol,
                                 const double user_dconv, const Molecule& requested,
                                 const Representation wanted,
                                 const double eprec = 0.0,
                                 const std::string& xc = "",
                                 const std::string& ncf = "") {

    MADNESS_CHECK_THROW(not protocol.empty(), "empty protocol in plan_restart");
    const std::size_t last = protocol.size() - 1;

    const bool ao_available = disk.restartao_present and can.ao;
    const bool nwchem_available = disk.nwfile_named and can.nwchem;

    // the precision this run has to reach. dconv cannot be demanded tighter
    // than the final rung represents, which is what SCFProtocol also does.
    const double target_thresh = protocol.back();
    const double target_dconv = std::max(target_thresh, user_dconv);

    RestartPlan plan;
    plan.mode = mode;

    // does the archive solve the same Hamiltonian this run is asking about?
    //
    // eprec, xc and the nuclear correlation factor all change the operator, not
    // just its representation: orbitals from a different one are a perfectly good
    // guess, but their convergence claim is about another problem, and -- unlike a
    // k or thresh mismatch -- that cannot be repaired by reprojecting. An unset
    // value on either side (0.0 / "") means "not recorded" -- v4 archives, the
    // seeding tools, and engines that have no ncf -- which is not evidence of a
    // mismatch. Returns an empty string when the Hamiltonians agree, otherwise the
    // phrase that goes into `why`.
    auto hamiltonian_mismatch = [&](const RestartMetadata& meta) {
        if (meta.eprec != 0.0 and eprec != 0.0 and
                std::abs(meta.eprec / eprec - 1.0) > 1.e-10)
            return "archive was written at eprec " + format_thresh(meta.eprec) +
                   ", this run uses " + format_thresh(eprec);
        if (not meta.xc.empty() and not xc.empty() and meta.xc != xc)
            return "archive was written with xc '" + meta.xc + "', this run uses '" +
                   xc + "'";
        if (not meta.ncf.empty() and not ncf.empty() and meta.ncf != ncf)
            return "archive was written with the nuclear correlation factor '" +
                   meta.ncf + "', this run uses '" + ncf + "'";
        return std::string();
    };

    // ---- fall back to the initial guess, recording why ---------------------
    auto give_up = [&](const std::string& why) {
        plan.source = RestartSource::initial_guess;
        plan.iterate = true;
        plan.protocol_start = 0;
        plan.why = why;
        return plan;
    };

    // ---- use restartdata, iterating from wherever it left off --------------
    auto continue_from_archive = [&](const RestartMetadata& meta, const std::string& why) {
        plan.source = RestartSource::restartdata;
        plan.stale_energy = meta.current_energy;
        const auto rung = first_rung_tighter_than(protocol, meta.converged_for_thresh);
        plan.iterate = true;
        plan.protocol_start = rung.value_or(last);
        plan.why = why;
        return plan;
    };

    if (mode == RestartMode::none) {
        // Contradictory input, and worth saying so here: with `nwfile` set, the AO
        // basis itself is read from the NWChem file (SCF's ctor), so there is no
        // atomic guess left to fall back to and SCF::initial_guess asserts with a
        // message that explains nothing.
        MADNESS_CHECK_THROW(not nwchem_available,
                "restart none together with `nwfile`: the AO basis is read from the "
                "nwchem file, so there is no atomic guess to fall back to. Use "
                "restart=nwchem, or drop the `nwfile` keyval.");
        return give_up("restart none: ignoring anything on disk");
    }

    // ---- explicitly requested sources: never fall back silently ------------
    //
    // A silent fallback here burns hours: the user asked to resume a long run
    // and would get a fresh one instead, indistinguishable from the outside
    // until the wall clock says so.
    if (mode == RestartMode::nwchem) {
        MADNESS_CHECK_THROW(can.nwchem,
                "restart nwchem: this calculation cannot read NWChem orbitals");
        MADNESS_CHECK_THROW(disk.nwfile_named,
                "restart nwchem: no nwchem file given -- set the `nwfile` keyval");
        plan.source = RestartSource::nwchem;
        plan.iterate = true;
        plan.protocol_start = 0;
        plan.why = "restart nwchem: as requested";
        return plan;
    }

    if (mode == RestartMode::ao) {
        MADNESS_CHECK_THROW(can.ao,
                "restart ao: this calculation cannot read AO projections");
        MADNESS_CHECK_THROW(disk.restartao_present,
                "restart ao: no restartaodata file found");
        plan.source = RestartSource::restartao;
        plan.iterate = true;
        plan.protocol_start = 0;
        plan.why = "restart ao: as requested";
        return plan;
    }

    if (mode == RestartMode::iterate or mode == RestartMode::read_only) {
        MADNESS_CHECK_THROW(disk.restartdata_present,
                "restart iterate/read_only: no restartdata archive found");
        MADNESS_CHECK_THROW(disk.meta.has_value(),
                "restart iterate/read_only: the restartdata archive could not be read");
        const RestartMetadata& meta = disk.meta.value();
        MADNESS_CHECK_THROW(meta.representation_matches(wanted),
                "restart iterate/read_only: the restartdata archive holds a different "
                "kind of orbital than this calculation uses");

        if (mode == RestartMode::read_only) {
            MADNESS_CHECK_THROW(
                    compare_geometry(meta.molecule, requested) == GeometryMatch::same,
                    "restart read_only: archive geometry does not match the requested geometry");
            // The user asserted these orbitals are the answer. Respect that even
            plan.source = RestartSource::restartdata;
            plan.iterate = false;
            plan.protocol_start = last;
            plan.stale_energy = meta.current_energy;
            const std::string other = hamiltonian_mismatch(meta);
            if (not other.empty()) {
                // Not overridden -- the user asked for these orbitals and gets
                // them -- but handing back an energy for a different Hamiltonian
                // without saying so is how a wrong number ends up in a table.
                plan.warn = true;
                plan.why = "restart read_only: " + other +
                           " -- returning the archive's energy for a DIFFERENT "
                           "Hamiltonian, as requested";
                return plan;
            }
            const bool good = meta.is_converged_to(target_thresh, target_dconv);
            plan.warn = not good;
            plan.why = good
                    ? "restart read_only: archive is converged to the requested precision"
                    : "restart read_only: archive is converged only to thresh " +
                      format_thresh(meta.converged_for_thresh) + " dconv " +
                      format_thresh(meta.converged_for_dconv) + ", NOT to the requested " +
                      format_thresh(target_thresh) + "/" + format_thresh(target_dconv) +
                      " -- returning its energy anyway, as requested";
            return plan;
        }
        // iterate: honour the request even if the archive already looks converged,
        // and in that case re-verify at the final rung rather than doing nothing.
        return continue_from_archive(meta, "restart iterate: as requested");
    }

    // ---- automatic ---------------------------------------------------------
    MADNESS_CHECK(mode == RestartMode::automatic);

    if (not disk.restartdata_present) {
        // nwchem before the AO projections, which reverses the order the old
        // precedence ladder used (SCF::get_initial_orbitals: restartdata, ao,
        // NWChem, guess). Naming a file in the input is a statement of intent; a
        // leftover restartaodata is not. In practice the two rarely coexist --
        // save_mos deliberately writes no restartaodata when nwfile is set,
        // because the AO basis then comes from the nwchem file.
        if (nwchem_available) {
            plan.source = RestartSource::nwchem;
            plan.why = "restart auto: no restartdata, but an nwchem file was given";
            return plan;
        }
        if (ao_available) {
            plan.source = RestartSource::restartao;
            plan.why = "restart auto: no restartdata, using the AO projections";
            return plan;
        }
        return give_up("restart auto: nothing on disk");
    }

    // a file is there. Anything below this point that rejects it is a surprise
    // and must be reported, not swallowed.
    auto reject_archive = [&](const std::string& why) {
        plan.warn = true;
        if (ao_available) {
            plan.source = RestartSource::restartao;
            plan.iterate = true;
            plan.protocol_start = 0;
            plan.why = why + "; using the AO projections instead";
            return plan;
        }
        return give_up(why + "; starting from the initial guess instead");
    };

    if (not disk.meta.has_value())
        return reject_archive("restart auto: restartdata exists but its header could "
                              "not be read (truncated, or written by a newer MADNESS)");

    const RestartMetadata& meta = disk.meta.value();

    if (not meta.representation_matches(wanted))
        return reject_archive("restart auto: restartdata holds '" +
                              madness::to_string(meta.representation) +
                              "' orbitals, this calculation wants '" +
                              madness::to_string(wanted) + "'");

    const GeometryMatch geom = compare_geometry(meta.molecule, requested);
    if (geom == GeometryMatch::different_composition) {
        // Not a fault: a different molecule in the same directory. The AO file
        // is indexed by the same atoms, so it is no better.
        plan.warn = false;
        return give_up("restart auto: restartdata is for a different molecule");
    }
    if (geom == GeometryMatch::displaced) {
        // A geometry optimization step. MRA orbitals at the old geometry are
        // wrong where it matters most -- at the nuclei -- so go through the AO
        // expansion, which is re-centred on the new positions.
        if (ao_available) {
            plan.source = RestartSource::restartao;
            plan.iterate = true;
            plan.protocol_start = 0;
            plan.why = "restart auto: geometry moved, using the AO projections";
            return plan;
        }
        return give_up("restart auto: geometry moved and no AO projections available");
    }

    // A different eprec, functional or nuclear correlation factor is a different
    // Hamiltonian: keep the orbitals as a guess, but throw away the archive's
    // convergence claim, which is about another problem. Without this an
    // `xc=lda` run in a directory holding a converged `xc=hf` archive skips the
    // SCF entirely and reports the HF energy.
    const std::string other = hamiltonian_mismatch(meta);
    if (not other.empty())
        return continue_from_archive(meta, "restart auto: " + other +
                                           " -- a different Hamiltonian, so re-converging");

    if (meta.is_converged_to(target_thresh, target_dconv)) {
        plan.source = RestartSource::restartdata;
        plan.iterate = false;
        plan.protocol_start = last;
        plan.stale_energy = meta.current_energy;
        plan.why = "restart auto: archive is converged to thresh " +
                   format_thresh(target_thresh) + " and dconv " +
                   format_thresh(target_dconv);
        return plan;
    }

    return continue_from_archive(meta, "restart auto: archive converged only to thresh " +
                                       format_thresh(meta.converged_for_thresh) +
                                       " dconv " + format_thresh(meta.converged_for_dconv) +
                                       ", continuing");
}

/// look at the disk: what restart data is there?
///
/// The existence tests run on rank 0 and are broadcast, so ranks cannot diverge
/// on them. peek_restartdata() is collective by construction (the parallel
/// archive broadcasts what it reads), so every rank passes through it together
/// or none does.
inline RestartSources survey_restart_sources(World& world, const std::string& prefix,
                                             const bool nwfile_named) {
    RestartSources disk;
    disk.nwfile_named = nwfile_named;

    int flags[2] = {0, 0};
    if (world.rank() == 0) {
        flags[0] = restartdata_exists(prefix + ".restartdata") ? 1 : 0;
        flags[1] = std::filesystem::exists(prefix + ".restartaodata") ? 1 : 0;
    }
    world.gop.broadcast(flags, 2, 0);
    disk.restartdata_present = (flags[0] == 1);
    disk.restartao_present = (flags[1] == 1);

    if (disk.restartdata_present) disk.meta = peek_restartdata(world, prefix + ".restartdata");
    return disk;
}

/// survey the disk and decide, identically on every rank
///
/// Agreement between ranks is the point. Ranks that disagree about where the
/// orbitals come from diverge on collective Function construction and hang -- a
/// hang is a far worse failure than a wrong answer, because it leaves nothing to
/// debug. Two things secure it:
///
///  * every input to plan_restart() is already rank-invariant -- the existence
///    tests are broadcast by survey_restart_sources(), the header is broadcast by
///    the parallel archive, and the parameters and molecule are replicated -- so
///    all ranks run the decision and any throw from an explicitly requested but
///    missing source is collective rather than a rank-0 abort into a hang;
///  * the result is then broadcast anyway, so the "identical inputs, identical
///    output" argument does not have to stay true for correctness.
///
/// @param[in] can which sources the asking engine can read; see RestartCapabilities
/// @param[in] ncf the nuclear correlation factor this run uses (SCF::restart_ncf);
///                empty for an engine that has none
inline RestartPlan make_restart_plan(World& world, const RestartMode mode,
                                     const CalculationParameters& param,
                                     const Molecule& requested,
                                     const Representation wanted,
                                     const RestartCapabilities& can,
                                     const std::string& ncf = "") {

    const RestartSources disk =
            survey_restart_sources(world, param.prefix(), param.nwfile() != "none");

    RestartPlan plan = plan_restart(mode, disk, can, param.protocol(), param.dconv(),
                                    requested, wanted, requested.parameters.eprec(),
                                    param.xc(), ncf);
    world.gop.broadcast_serializable(plan, 0);

    if (world.rank() == 0 and param.print_level() > 1) {
        if (plan.warn) print("WARNING:", plan.why);
        else print(plan.print_to_string());
        if (disk.meta.has_value() and param.print_level() > 2)
            print("   archive:", disk.meta.value().print_to_string());
    }
    return plan;
}

} // namespace madness

#endif // MADNESS_CHEM_RESTARTPLAN_H__INCLUDED
