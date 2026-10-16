/// \file test_restart.cc
/// \brief the restart decision table, checked line by line
///
/// plan_restart() is a pure function of (mode, what is on disk, what the engine
/// can read, protocol, dconv, geometry, representation). That is the whole
/// reason it was factored out of SCF and Nemo: the decision table can be walked
/// exhaustively here with no MRA data, no SCF, no archive and no filesystem.

#include <madness/chem/RestartPlan.h>
#include <madness/world/MADworld.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

using namespace madness;

namespace {

int failures = 0;
int checks = 0;

void check(const bool ok, const std::string& what) {
    ++checks;
    if (not ok) {
        ++failures;
        print("FAILED:", what);
    }
}

/// LiH, the geometry every archive in this test claims to hold
Molecule lih() {
    Molecule m;
    m.add_atom(0.0, 0.0, 0.0, 3.0, 3);
    m.add_atom(0.0, 0.0, 3.0, 1.0, 1);
    return m;
}

/// LiH with the H pushed out by 0.1 a0 -- a geometry optimization step
Molecule lih_displaced() {
    Molecule m;
    m.add_atom(0.0, 0.0, 0.0, 3.0, 3);
    m.add_atom(0.0, 0.0, 3.1, 1.0, 1);
    return m;
}

/// BeH2: a different molecule in the same directory
Molecule beh2() {
    Molecule m;
    m.add_atom(0.0, 0.0, 0.0, 4.0, 4);
    m.add_atom(0.0, 0.0, 2.5, 1.0, 1);
    m.add_atom(0.0, 0.0, -2.5, 1.0, 1);
    return m;
}

/// a version-5 archive header for LiH, converged as told
RestartMetadata converged_archive(const double thresh, const double dconv,
                                  const Representation rep = Representation::mo,
                                  const Molecule& mol = lih()) {
    RestartMetadata meta;
    meta.version = 5;
    meta.molecule = mol;
    meta.k = 6;
    meta.L = 20.0;
    meta.current_energy = -8.0;
    meta.converged_for_thresh = thresh;
    meta.converged_for_dconv = dconv;
    meta.representation = rep;
    return meta;
}

/// a version-4 archive: no representation and no converged_for_dconv recorded
RestartMetadata v4_archive(const double thresh) {
    RestartMetadata meta = converged_archive(thresh, 1.e10, Representation::unknown);
    meta.version = 4;
    return meta;
}

RestartSources with_archive(const RestartMetadata& meta, const bool ao = false,
                            const bool nw = false) {
    RestartSources disk;
    disk.restartdata_present = true;
    disk.meta = meta;
    disk.restartao_present = ao;
    disk.nwfile_named = nw;
    return disk;
}

RestartSources nothing(const bool ao = false, const bool nw = false) {
    RestartSources disk;
    disk.restartao_present = ao;
    disk.nwfile_named = nw;
    return disk;
}

const std::vector<double> ladder = {1.e-4, 1.e-6};
const double user_dconv = 1.e-4;

/// one row of the table
void expect(const std::string& name, const RestartMode mode, const RestartSources& disk,
            const RestartCapabilities& can, const RestartSource source,
            const bool iterate, const std::size_t start,
            const Molecule& requested = lih(),
            const Representation wanted = Representation::mo) {
    try {
        const RestartPlan plan =
                plan_restart(mode, disk, can, ladder, user_dconv, requested, wanted);
        check(plan.source == source, name + ": source is " +
                madness::to_string(plan.source) + ", expected " + madness::to_string(source));
        check(plan.iterate == iterate, name + ": iterate is " +
                std::to_string(plan.iterate) + ", expected " + std::to_string(iterate));
        if (iterate)
            check(plan.protocol_start == start, name + ": starts at rung " +
                    std::to_string(plan.protocol_start) + ", expected " + std::to_string(start));
    } catch (const std::exception& e) {
        check(false, name + ": threw unexpectedly");
    }
}

/// a row that must be refused rather than answered
void expect_throw(const std::string& name, const RestartMode mode,
                  const RestartSources& disk, const RestartCapabilities& can) {
    bool threw = false;
    try {
        plan_restart(mode, disk, can, ladder, user_dconv, lih(), Representation::mo);
    } catch (const MadnessException&) {
        threw = true;
    }
    check(threw, name + ": should have thrown, but returned a plan");
}

void test_helpers() {
    // first_rung_tighter_than: nullopt means the ladder is fully covered
    check(first_rung_tighter_than(ladder, 1.e10).value() == 0, "first rung: nothing achieved");
    check(first_rung_tighter_than(ladder, 1.e-4).value() == 1, "first rung: 1e-4 achieved");
    check(not first_rung_tighter_than(ladder, 1.e-6).has_value(), "first rung: 1e-6 achieved");
    // the 0.999 slack: a threshold that made a round trip through an archive and
    // came back a hair large must still count as covering its own rung
    check(not first_rung_tighter_than(ladder, 1.e-6 * 1.0005).has_value(),
          "first rung: 1e-6 achieved, within slack");

    check(compare_geometry(lih(), lih()) == GeometryMatch::same, "geometry: same");
    check(compare_geometry(lih(), lih_displaced()) == GeometryMatch::displaced,
          "geometry: displaced");
    check(compare_geometry(lih(), beh2()) == GeometryMatch::different_composition,
          "geometry: different composition");
    // the atom ORDER is part of the identity: the orbitals in an archive are
    // expanded about the atoms in the order that archive lists them
    Molecule hli;
    hli.add_atom(0.0, 0.0, 3.0, 1.0, 1);
    hli.add_atom(0.0, 0.0, 0.0, 3.0, 3);
    check(compare_geometry(lih(), hli) == GeometryMatch::different_composition,
          "geometry: permuted atoms are a different molecule");

    // --restart_info should accept whatever the user types, including the chunk
    // suffix that tab completion appends
    check(restartdata_basename("he") == "he", "basename: bare prefix");
    check(restartdata_basename("he.restartdata") == "he", "basename: with .restartdata");
    check(restartdata_basename("he.restartdata.00000") == "he", "basename: with chunk suffix");
    check(restartdata_basename("dir/he.restartdata.00000") == "dir/he",
          "basename: path is preserved");
    check(restartdata_basename("restartdata") == "restartdata",
          "basename: a bare name that IS the suffix is left alone");

    check(restart_mode_from_string("auto") == RestartMode::automatic, "parse auto");
    check(restart_mode_from_string("read_only") == RestartMode::read_only, "parse read_only");
    for (const auto& name : restart_mode_names())
        check(madness::to_string(restart_mode_from_string(name)) == name,
              "mode round trip: " + name);
}

/// the keyval's allowed_values and this enum are two lists that must agree
///
/// They cannot be one list: CalculationParameters.h is included BY RestartPlan.h,
/// so the dependency cannot run the other way. This walks every spelling the
/// parser knows through the keyval and back, which catches a mode added to the
/// enum but not to the keyval. The reverse drift -- a spelling the keyval accepts
/// but the parser does not -- surfaces immediately and loudly as "unknown restart
/// mode" the first time anyone uses it.
void test_keyval_agrees_with_enum() {
    for (const auto& name : restart_mode_names()) {
        CalculationParameters param;
        bool ok = true;
        try {
            param.set_user_defined_value("restart", name);
        } catch (const std::exception&) {
            ok = false;
        }
        check(ok, "keyval accepts restart=" + name);
        if (ok) check(param.restart() == name, "keyval stores restart=" + name);
        if (ok) check(madness::to_string(restart_mode_from_string(param.restart())) == name,
                      "keyval round trips restart=" + name);
    }
    // and the list is actually constraining
    CalculationParameters param;
    bool rejected = false;
    try {
        param.set_user_defined_value("restart", std::string("1"));
    } catch (const std::exception&) {
        rejected = true;
    }
    check(rejected, "keyval rejects the retired boolean spelling `restart 1`");

    check(CalculationParameters().restart() == "auto", "restart defaults to auto");
}

void test_automatic() {
    const auto moldft = RestartCapabilities::all();
    const auto nemo = RestartCapabilities::restartdata_only();
    const auto A = RestartMode::automatic;

    expect("auto, empty directory", A, nothing(), moldft,
           RestartSource::initial_guess, true, 0);

    // fully converged: hand back what is on disk, do not solve anything
    expect("auto, converged archive", A, with_archive(converged_archive(1.e-6, 1.e-4)), moldft,
           RestartSource::restartdata, false, 1);

    // converged only through the first rung: resume at the second, which is the
    // whole point -- rung 0 is not repeated
    expect("auto, archive at 1e-4", A, with_archive(converged_archive(1.e-4, 1.e-3)), moldft,
           RestartSource::restartdata, true, 1);

    // right thresh, but never converged in dconv (this is every v4 archive)
    expect("auto, thresh met but dconv unknown", A, with_archive(v4_archive(1.e-6)), moldft,
           RestartSource::restartdata, true, 1);

    // dconv is relaxed to the final rung: a user asking for dconv 1e-12 on a
    // 1e-6 ladder cannot demand tighter than the representation supports
    {
        const RestartPlan plan = plan_restart(A, with_archive(converged_archive(1.e-6, 1.e-6)),
                moldft, ladder, 1.e-12, lih(), Representation::mo);
        check(not plan.iterate, "auto: dconv is relaxed to the final rung");
    }

    // a file that cannot be parsed is a damaged archive, not an empty directory
    {
        RestartSources disk = nothing();
        disk.restartdata_present = true;   // present, but meta stayed nullopt
        const RestartPlan plan = plan_restart(A, disk, moldft, ladder, user_dconv,
                                              lih(), Representation::mo);
        check(plan.source == RestartSource::initial_guess, "auto, corrupt: falls back");
        check(plan.warn, "auto, corrupt: reports the unusable file");

        disk.restartao_present = true;
        const RestartPlan plan2 = plan_restart(A, disk, moldft, ladder, user_dconv,
                                               lih(), Representation::mo);
        check(plan2.source == RestartSource::restartao, "auto, corrupt: prefers the AO file");
        check(plan2.warn, "auto, corrupt with AO: still reports it");
    }

    // nemo orbitals are F = psi/R; moldft must not read them as psi
    expect("auto, nemo archive for moldft", A,
           with_archive(converged_archive(1.e-6, 1.e-4, Representation::nemo), true), moldft,
           RestartSource::restartao, true, 0);
    expect("auto, moldft archive for nemo", A,
           with_archive(converged_archive(1.e-6, 1.e-4, Representation::mo)), nemo,
           RestartSource::initial_guess, true, 0, lih(), Representation::nemo);

    // a v4 archive records no representation. Rejecting all of them would break
    // every archive written before the field existed, so "unknown" is accepted.
    expect("auto, v4 archive is representation-agnostic", A,
           with_archive(v4_archive(1.e-6)), nemo,
           RestartSource::restartdata, true, 1, lih(), Representation::nemo);

    // a geometry optimization step: through the AO expansion, which is
    // re-centred on the new nuclei, or from scratch if there is none
    expect("auto, geometry moved", A,
           with_archive(converged_archive(1.e-6, 1.e-4), true), moldft,
           RestartSource::restartao, true, 0, lih_displaced());
    expect("auto, geometry moved, no AO file", A,
           with_archive(converged_archive(1.e-6, 1.e-4)), moldft,
           RestartSource::initial_guess, true, 0, lih_displaced());

    // a different molecule in the same directory. The AO file is indexed by the
    // same atoms, so it is no better -- and this is not a fault to warn about.
    {
        const RestartPlan plan = plan_restart(A,
                with_archive(converged_archive(1.e-6, 1.e-4), true), moldft, ladder,
                user_dconv, beh2(), Representation::mo);
        check(plan.source == RestartSource::initial_guess, "auto, other molecule: fresh start");
        check(not plan.warn, "auto, other molecule: not a fault");
    }

    // eprec changes the nuclear potential, so a converged archive at a different
    // eprec is a good guess for a different Hamiltonian -- iterate, do not trust
    // its convergence. Reprojecting cannot repair this the way it repairs k.
    {
        RestartMetadata meta = converged_archive(1.e-6, 1.e-4);
        meta.eprec = 1.e-4;
        const RestartPlan same = plan_restart(A, with_archive(meta), moldft, ladder,
                user_dconv, lih(), Representation::mo, HamiltonianKey{.eprec = 1.e-4});
        check(not same.iterate, "auto: matching eprec keeps the convergence claim");

        const RestartPlan diff = plan_restart(A, with_archive(meta), moldft, ladder,
                user_dconv, lih(), Representation::mo, HamiltonianKey{.eprec = 1.e-6});
        check(diff.iterate, "auto: differing eprec re-converges");
        check(diff.source == RestartSource::restartdata,
              "auto: differing eprec still uses the orbitals as a guess");

        // 0 on either side means "not recorded", which is not a mismatch
        const RestartPlan unknown = plan_restart(A, with_archive(meta), moldft, ladder,
                user_dconv, lih(), Representation::mo, HamiltonianKey{});
        check(not unknown.iterate, "auto: unknown requested eprec is not a mismatch");
        const RestartPlan unknown2 = plan_restart(A,
                with_archive(converged_archive(1.e-6, 1.e-4)), moldft, ladder, user_dconv,
                lih(), Representation::mo, HamiltonianKey{.eprec = 1.e-6});
        check(not unknown2.iterate, "auto: unrecorded archive eprec is not a mismatch");
    }

    // The functional is the Hamiltonian just as eprec is: an `xc=lda` run in a
    // directory holding a converged `xc=hf` archive must NOT skip the SCF and
    // report the HF energy. The orbitals are still the best guess available.
    {
        RestartMetadata meta = converged_archive(1.e-6, 1.e-4);
        meta.xc = "hf";
        const RestartPlan same = plan_restart(A, with_archive(meta), moldft, ladder,
                user_dconv, lih(), Representation::mo, HamiltonianKey{.xc = "hf"});
        check(not same.iterate, "auto: matching xc keeps the convergence claim");

        const RestartPlan diff = plan_restart(A, with_archive(meta), moldft, ladder,
                user_dconv, lih(), Representation::mo, HamiltonianKey{.xc = "lda"});
        check(diff.iterate, "auto: a different functional re-converges");
        check(diff.source == RestartSource::restartdata,
              "auto: a different functional still uses the orbitals as a guess");
        // the convergence CLAIM is dropped, but not the precision the orbitals are
        // represented at: they are accurate functions of the wrong operator, so
        // there is nothing to gain from reprojecting them back down to rung 0
        check(diff.protocol_start == 1,
              "auto: a different functional resumes at the archive's precision");

        // an empty string on either side is "not recorded", not a mismatch
        const RestartPlan unknown = plan_restart(A, with_archive(meta), moldft, ladder,
                user_dconv, lih(), Representation::mo, HamiltonianKey{});
        check(not unknown.iterate, "auto: unknown requested xc is not a mismatch");
        const RestartPlan unknown2 = plan_restart(A,
                with_archive(converged_archive(1.e-6, 1.e-4)), moldft, ladder, user_dconv,
                lih(), Representation::mo, HamiltonianKey{.xc = "lda"});
        check(not unknown2.iterate, "auto: unrecorded archive xc is not a mismatch");
    }

    // the same for the nuclear correlation factor: two nemo runs in one directory
    // that differ only in `ncf` solve different Hamiltonians under the same
    // representation tag, so the representation check cannot catch them
    {
        RestartMetadata meta = converged_archive(1.e-6, 1.e-4, Representation::nemo);
        meta.ncf = "slater:2.000000";
        const RestartPlan same = plan_restart(A, with_archive(meta), nemo, ladder,
                user_dconv, lih(), Representation::nemo, HamiltonianKey{.ncf = "slater:2.000000"});
        check(not same.iterate, "auto: matching ncf keeps the convergence claim");

        const RestartPlan diff = plan_restart(A, with_archive(meta), nemo, ladder,
                user_dconv, lih(), Representation::nemo, HamiltonianKey{.ncf = "slater:1.000000"});
        check(diff.iterate, "auto: a different ncf re-converges");
        check(diff.source == RestartSource::restartdata,
              "auto: a different ncf still uses the orbitals as a guess");
    }

    // nwchem, when there is nothing better and the engine can read it
    expect("auto, nwchem file named", A, nothing(false, true), moldft,
           RestartSource::nwchem, true, 0);
    expect("auto, nwchem named but unreadable by nemo", A, nothing(false, true), nemo,
           RestartSource::initial_guess, true, 0, lih(), Representation::nemo);
    // restartdata wins over both
    expect("auto, restartdata beats nwchem and ao", A,
           with_archive(converged_archive(1.e-4, 1.e-3), true, true), moldft,
           RestartSource::restartdata, true, 1);

    // An archive holding fewer orbitals than requested (virtuals added on
    // restart) is a guess for the missing ones, however converged it is; one
    // holding as many or more is used as it stands.
    {
        auto short_archive = [](const std::size_t nmo) {
            RestartSources disk = with_archive(converged_archive(1.e-6, 1.e-4));
            disk.nmo_alpha = nmo;
            return disk;
        };
        auto plan_for = [&](const RestartMode mode, const std::size_t nmo) {
            return plan_restart(mode, short_archive(nmo), moldft, ladder, user_dconv, lih(),
                                Representation::mo, HamiltonianKey{}, 7);
        };
        const RestartPlan fewer = plan_for(A, 5);
        check(fewer.source == RestartSource::restartdata, "auto, short archive: reads the archive");
        check(fewer.iterate, "auto, short archive: iterates");
        check(fewer.protocol_start == 1, "auto, short archive: at the archive's rung");
        check(fewer.archive_nmo_alpha == 5, "auto, short archive: records the archive's count");
        check(fewer.why.find("missing virtuals") != std::string::npos,
              "auto, short archive: says why");
        check(not plan_for(A, 7).iterate, "auto, archive with the requested count: no iterations");
        check(not plan_for(A, 9).iterate, "auto, archive with more orbitals: no iterations");
        check(plan_for(A, 0).iterate == false, "auto, archive count unknown: not treated as short");

        const RestartPlan iterate = plan_for(RestartMode::iterate, 5);
        check(iterate.source == RestartSource::restartdata and iterate.iterate and
              iterate.protocol_start == 1, "iterate, short archive: continues from the archive");

        bool threw = false;
        try {
            plan_for(RestartMode::read_only, 5);
        } catch (const MadnessException&) { threw = true; }
        check(threw, "read_only, short archive: refused");
        check(not plan_for(RestartMode::read_only, 7).iterate, "read_only, full archive: accepted");
    }

    // what the plan tells freeze_occupied about the archive
    {
        const RestartPlan good = plan_restart(A, with_archive(converged_archive(1.e-4, 1.e-3)), moldft,
                                              ladder, user_dconv, lih(), Representation::mo, HamiltonianKey{.xc = "hf"});
        check(good.archive_converged and good.archive_same_hamiltonian,
              "auto, converged archive: converged for this Hamiltonian");
        const RestartPlan never = plan_restart(A, with_archive(converged_archive(1.e10, 1.e10)), moldft,
                                               ladder, user_dconv, lih(), Representation::mo);
        check(never.source == RestartSource::restartdata and not never.archive_converged,
              "auto, archive that never converged: read, but not converged");
        RestartMetadata lda = converged_archive(1.e-6, 1.e-4);
        lda.xc = "lda";
        const RestartPlan other = plan_restart(A, with_archive(lda), moldft, ladder, user_dconv, lih(),
                                               Representation::mo, HamiltonianKey{.xc = "hf"});
        check(other.source == RestartSource::restartdata and other.iterate and
              other.archive_converged and not other.archive_same_hamiltonian,
              "auto, xc mismatch: re-converged, flagged as another Hamiltonian");
    }
}

void test_explicit_modes() {
    const auto moldft = RestartCapabilities::all();
    const auto nemo = RestartCapabilities::restartdata_only();

    // none ignores everything, however good it is
    expect("none, converged archive present", RestartMode::none,
           with_archive(converged_archive(1.e-6, 1.e-4), true), moldft,
           RestartSource::initial_guess, true, 0);
    // ... except that `none` plus an nwfile is contradictory, since with an nwfile
    // the AO basis itself comes from that file and there is no atomic guess left
    expect_throw("none together with nwfile", RestartMode::none, nothing(false, true), moldft);
    // for nemo the nwfile is ignored entirely, so there is nothing to contradict
    expect("none with nwfile, engine ignores nwchem", RestartMode::none,
           nothing(false, true), nemo, RestartSource::initial_guess, true, 0, lih(),
           Representation::nemo);

    // read_only respects the user even when the archive is worse than requested:
    // warn and hand back the stale energy rather than second-guessing
    expect("read_only, converged archive", RestartMode::read_only,
           with_archive(converged_archive(1.e-6, 1.e-4)), moldft,
           RestartSource::restartdata, false, 1);
    {
        const RestartPlan plan = plan_restart(RestartMode::read_only,
                with_archive(converged_archive(1.e-3, 1.e-2)), moldft, ladder, user_dconv,
                lih(), Representation::mo);
        check(not plan.iterate, "read_only, unconverged archive: still no iterations");
        check(plan.stale_energy == -8.0, "read_only: carries the archive's energy");
        check(plan.warn, "read_only, unconverged archive: warns");
        check(plan.why.find("NOT to the requested") != std::string::npos,
              "read_only, unconverged archive: says so");
    }
    {
        const RestartPlan plan = plan_restart(RestartMode::read_only,
                with_archive(converged_archive(1.e-6, 1.e-4)), moldft, ladder, user_dconv,
                lih(), Representation::mo);
        check(not plan.warn, "read_only, converged archive: no warning");
    }

    // iterate is honoured even on an archive that already looks done -- the
    // final rung is re-verified rather than skipped
    expect("iterate, converged archive", RestartMode::iterate,
           with_archive(converged_archive(1.e-6, 1.e-4)), moldft,
           RestartSource::restartdata, true, 1);
    expect("iterate, archive at 1e-4", RestartMode::iterate,
           with_archive(converged_archive(1.e-4, 1.e-3)), moldft,
           RestartSource::restartdata, true, 1);

    expect("ao, file present", RestartMode::ao, nothing(true), moldft,
           RestartSource::restartao, true, 0);
    expect("nwchem, file named", RestartMode::nwchem, nothing(false, true), moldft,
           RestartSource::nwchem, true, 0);

    // An explicitly named source that is not there must stop the run. Falling
    // back silently burns hours: someone resuming a long calculation would get a
    // fresh one instead, and only the wall clock would say so.
    expect_throw("iterate, nothing on disk", RestartMode::iterate, nothing(), moldft);
    expect_throw("read_only, nothing on disk", RestartMode::read_only, nothing(), moldft);
    expect_throw("ao, no AO file", RestartMode::ao, nothing(), moldft);
    expect_throw("nwchem, no file named", RestartMode::nwchem, nothing(), moldft);
    expect_throw("ao, engine cannot read it", RestartMode::ao, nothing(true), nemo);
    expect_throw("nwchem, engine cannot read it", RestartMode::nwchem,
                 nothing(false, true), nemo);
    {
        // present but unparseable is also a refusal, not a fallback
        RestartSources disk = nothing();
        disk.restartdata_present = true;
        expect_throw("iterate, corrupt archive", RestartMode::iterate, disk, moldft);
    }
    // an explicit mode will not read the wrong kind of orbital either
    {
        RestartSources disk =
                with_archive(converged_archive(1.e-6, 1.e-4, Representation::nemo));
        bool threw = false;
        try {
            plan_restart(RestartMode::iterate, disk, moldft, ladder, user_dconv, lih(),
                         Representation::mo);
        } catch (const MadnessException&) { threw = true; }
        check(threw, "iterate, nemo archive for moldft: refused");
    }
}

void test_single_rung_ladder() {
    // A one-rung protocol is common (every deck that writes `protocol [1e-4]`).
    // protocol_start must stay in range for it.
    const std::vector<double> one = {1.e-4};
    const auto moldft = RestartCapabilities::all();
    const RestartPlan a = plan_restart(RestartMode::automatic,
            with_archive(converged_archive(1.e-4, 1.e-4)), moldft, one, 1.e-4, lih(),
            Representation::mo);
    check(a.protocol_start == 0, "single rung, converged: rung 0");
    check(not a.iterate, "single rung, converged: no iterations");

    const RestartPlan b = plan_restart(RestartMode::automatic,
            with_archive(converged_archive(1.e-2, 1.e-2)), moldft, one, 1.e-4, lih(),
            Representation::mo);
    check(b.protocol_start == 0, "single rung, unconverged: rung 0");
    check(b.iterate, "single rung, unconverged: iterate");
}

void test_serialization() {
    // The plan is broadcast, so it has to survive a round trip through an
    // archive. A field forgotten in serialize() would let ranks disagree about
    // where the orbitals come from, and ranks that disagree hang.
    RestartPlan plan;
    plan.mode = RestartMode::read_only;
    plan.source = RestartSource::restartao;
    plan.iterate = false;
    plan.protocol_start = 3;
    plan.stale_energy = -8.25;
    plan.warn = true;
    plan.why = "a reason";
    plan.archive_nmo_alpha = 5;
    plan.archive_converged = true;
    plan.archive_same_hamiltonian = true;

    const std::size_t bufsize = 4096;
    std::vector<unsigned char> buf(bufsize);
    {
        archive::BufferOutputArchive ar(buf.data(), bufsize);
        ar & plan;
    }
    RestartPlan copy;
    {
        archive::BufferInputArchive ar(buf.data(), bufsize);
        ar & copy;
    }
    check(copy.mode == plan.mode, "serialize: mode");
    check(copy.source == plan.source, "serialize: source");
    check(copy.iterate == plan.iterate, "serialize: iterate");
    check(copy.protocol_start == plan.protocol_start, "serialize: protocol_start");
    check(copy.stale_energy == plan.stale_energy, "serialize: stale_energy");
    check(copy.warn == plan.warn, "serialize: warn");
    check(copy.why == plan.why, "serialize: why");
    check(copy.archive_nmo_alpha == plan.archive_nmo_alpha, "serialize: archive_nmo_alpha");
    check(copy.archive_converged == plan.archive_converged, "serialize: archive_converged");
    check(copy.archive_same_hamiltonian == plan.archive_same_hamiltonian,
          "serialize: archive_same_hamiltonian");
}

/// the operator key: one comparison for the planner, the header and the results file
void test_hamiltonian_key() {
    HamiltonianKey a;
    a.xc = "hf"; a.eprec = 1.e-4; a.field = {0.0, 0.0, 0.0}; a.core_type = "none";
    a.psp_calc = 0; a.pcm = "none"; a.ncf = "slater:1.0";

    check(a.mismatch(a).empty(), "key: equal to itself");
    check(a.mismatch(HamiltonianKey{}).empty(), "key: an unrecorded request is not a mismatch");
    check(HamiltonianKey{}.mismatch(a).empty(), "key: an unrecorded archive is not a mismatch");

    auto differs = [&](const std::string& name, auto&& change) {
        HamiltonianKey b = a;
        change(b);
        check(not a.mismatch(b).empty(), "key: " + name + " is a different Hamiltonian");
    };
    differs("xc", [](HamiltonianKey& b) { b.xc = "lda"; });
    differs("eprec", [](HamiltonianKey& b) { b.eprec = 1.e-6; });
    differs("ncf", [](HamiltonianKey& b) { b.ncf = "slater:2.0"; });
    differs("field", [](HamiltonianKey& b) { b.field = {0.0, 0.0, 1.e-3}; });
    differs("core_type", [](HamiltonianKey& b) { b.core_type = "mcp"; });
    differs("psp_calc", [](HamiltonianKey& b) { b.psp_calc = 1; });
    differs("pcm", [](HamiltonianKey& b) { b.pcm = "water:iefpcm:eps=78.39:probe=1.385"; });

    const HamiltonianKey back = HamiltonianKey::from_json(a.to_json());
    check(back == a, "key: json round trip");
    check(HamiltonianKey::from_json(nlohmann::json::object()) == HamiltonianKey{},
          "key: an empty json reads as unrecorded");

    // the planner sees the new fields through the header: a converged archive in a
    // different external field keeps its orbitals as a guess, but re-converges
    const RestartMode A = RestartMode::automatic;
    const RestartCapabilities moldft = RestartCapabilities::all();
    RestartMetadata meta = converged_archive(1.e-6, 1.e-4);
    meta.xc = "hf";
    HamiltonianKey requested;
    requested.xc = "hf";
    const RestartPlan same = plan_restart(A, with_archive(meta), moldft, ladder, user_dconv,
                                          lih(), Representation::mo, requested);
    check(not same.iterate, "key: a v5 archive leaves the new fields unrecorded");
}

/// localize selects which orbitals of the same solution are stored
void test_localize() {
    const RestartMode A = RestartMode::automatic;
    const RestartCapabilities moldft = RestartCapabilities::all();
    RestartMetadata meta = converged_archive(1.e-6, 1.e-4);
    meta.localize = "boys";
    auto plan_for = [&](const std::string& localize) {
        return plan_restart(A, with_archive(meta), moldft, ladder, user_dconv, lih(),
                            Representation::mo, HamiltonianKey{}, 0, localize);
    };
    check(not plan_for("boys").iterate, "localize: the same method reuses the archive");
    check(not plan_for("").iterate, "localize: an unrecorded request is not a mismatch");
    const RestartPlan other = plan_for("canon");
    check(other.iterate, "localize: a different method iterates");
    check(other.source == RestartSource::restartdata, "localize: from the archive's orbitals");
    check(other.protocol_start == ladder.size() - 1, "localize: at the final rung only");
    check(other.archive_same_hamiltonian, "localize: still the same Hamiltonian");
}

/// header v6: the id, the electron counts and the rest of the key survive a round
/// trip; a v5 header still reads, with everything v6 added left unrecorded
void test_header_v6(World& world) {
    const std::size_t bufsize = 1 << 16;
    std::vector<unsigned char> buf(bufsize);

    RestartMetadata meta = converged_archive(1.e-6, 1.e-4);
    meta.version = RestartMetadata::CURRENT_VERSION;
    meta.archive_id = new_archive_id(world);
    meta.origin = "scf";
    meta.nalpha = 2; meta.nbeta = 1; meta.nmo_beta = 3;
    HamiltonianKey key;
    key.xc = "lda"; key.eprec = 1.e-4; key.field = {0.0, 0.0, 1.e-3}; key.core_type = "none";
    key.psp_calc = 0; key.pcm = "none";
    meta.set_hamiltonian_key(key);
    {
        archive::BufferOutputArchive ar(buf.data(), bufsize);
        meta.write(ar);
    }
    RestartMetadata back;
    {
        archive::BufferInputArchive ar(buf.data(), bufsize);
        back.read(ar);
    }
    check(back.version == 6, "v6: version");
    check(back.archive_id == meta.archive_id and back.archive_id != 0, "v6: archive_id");
    check(back.origin == "scf", "v6: origin");
    check(back.nalpha == 2 and back.nbeta == 1 and back.nmo_beta == 3, "v6: electron counts");
    check(back.hamiltonian_key() == key, "v6: the whole Hamiltonian key");

    // a v5 writer: the v4 fields, then the v5 tail, and nothing after it
    {
        archive::BufferOutputArchive ar(buf.data(), bufsize);
        const unsigned int v = 5;
        const int rep = static_cast<int>(Representation::mo);
        ar & v & meta.current_energy & meta.spin_restricted & meta.L & meta.k & meta.molecule
           & meta.xc & meta.localize & meta.converged_for_thresh;
        ar & meta.converged_for_dconv & rep & meta.ncf & meta.eprec & meta.madness_version;
    }
    RestartMetadata v5;
    {
        archive::BufferInputArchive ar(buf.data(), bufsize);
        v5.read(ar);
    }
    check(v5.version == 5 and v5.archive_id == 0, "v5: reads, with no archive id");
    check(v5.nalpha == -1 and v5.psp_calc == -1 and v5.field.empty(),
          "v5: the v6 fields read as not recorded");
    check(v5.hamiltonian_key().xc == "lda", "v5: the fields it had are kept");

    const ArchiveId a = new_archive_id(world), b = new_archive_id(world);
    check(a != 0 and b != 0 and a != b, "new_archive_id: nonzero and fresh every time");
    check(archive_id_from_string(archive_id_to_string(a)) == a, "archive id: hex round trip");
    check(archive_id_from_string("not hex") == 0, "archive id: garbage reads as none");
}

/// the AO projections are used only if their header matches the atoms and the archive
void test_ao_header(World& world) {
    const RestartMode A = RestartMode::automatic;
    const RestartCapabilities moldft = RestartCapabilities::all();
    RestartMetadata meta = converged_archive(1.e-6, 1.e-4);
    meta.archive_id = 42;

    auto displaced_plan = [&](const std::optional<RestartAOHeader>& header) {
        RestartSources disk = with_archive(meta, true);
        disk.ao_header = header;
        return plan_restart(A, disk, moldft, ladder, user_dconv, lih_displaced(),
                            Representation::mo);
    };
    RestartAOHeader h;
    h.archive_id = 42;
    h.molecule = lih();
    check(displaced_plan(h).source == RestartSource::restartao, "ao: same id is used");
    check(displaced_plan(std::nullopt).source == RestartSource::restartao,
          "ao: a legacy file without a header is used");
    h.archive_id = 43;
    check(displaced_plan(h).source == RestartSource::initial_guess,
          "ao: projections of another archive are not used");
    h.archive_id = 42;
    h.molecule = beh2();
    check(displaced_plan(h).source == RestartSource::initial_guess,
          "ao: projections for other atoms are not used");

    // the plan names the archive it read
    const RestartPlan plan = plan_restart(A, with_archive(meta), moldft, ladder, user_dconv,
                                          lih(), Representation::mo);
    check(plan.archive_id == 42, "plan: carries the archive id");

    // present() tells a headed file from a legacy one without a failing read
    if (world.rank() == 0) {
        {
            archive::BinaryFstreamOutputArchive ar("test_restart_ao_new");
            RestartAOHeader w;
            w.archive_id = 7;
            w.molecule = lih();
            w.write(ar);
            ar << Tensor<double>(2, 2);
        }
        {
            archive::BinaryFstreamOutputArchive ar("test_restart_ao_old");
            ar << Tensor<double>(2, 2);
        }
        check(RestartAOHeader::present("test_restart_ao_new"), "ao: a headed file is recognized");
        check(not RestartAOHeader::present("test_restart_ao_old"), "ao: a legacy file is recognized");
        const auto peeked = RestartAOHeader::peek("test_restart_ao_new");
        check(peeked.has_value() and peeked->archive_id == 7, "ao: peek reads the id");
        std::filesystem::remove("test_restart_ao_new");
        std::filesystem::remove("test_restart_ao_old");
    }
    world.gop.fence();
}

/// writing under a temporary name and renaming leaves one consistent archive
void test_commit(World& world) {
    const std::string name = "test_restart_commit.restartdata";
    // stale parts from an earlier, wider write, and a previous archive
    if (world.rank() == 0) {
        for (const char* p : {".00000", ".00001", ".00002"}) std::ofstream(name + p) << "old";
    }
    world.gop.fence();
    {
        archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(world, (name + ".tmp").c_str(), 1);
        ar & 12345;
    }
    commit_parallel_archive(world, name + ".tmp", name);
    int value = 0;
    {
        archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(world, name.c_str());
        ar & value;
    }
    check(value == 12345, "commit: the new archive is in place");
    if (world.rank() == 0) {
        check(not std::filesystem::exists(name + ".tmp.00000"), "commit: no temporary left");
        check(not std::filesystem::exists(name + ".00001") and
              not std::filesystem::exists(name + ".00002"), "commit: stale parts removed");
        std::filesystem::remove(name + ".00000");
    }
    world.gop.fence();
}

}  // namespace

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);

    test_helpers();
    test_keyval_agrees_with_enum();
    test_automatic();
    test_explicit_modes();
    test_single_rung_ladder();
    test_serialization();
    test_hamiltonian_key();
    test_localize();
    test_header_v6(world);
    test_ao_header(world);
    test_commit(world);

    if (world.rank() == 0) {
        print("");
        print("test_restart:", checks - failures, "of", checks, "checks passed");
        print(failures == 0 ? "test_restart passed" : "test_restart FAILED");
    }

    finalize();
    return failures == 0 ? 0 : 1;
}
