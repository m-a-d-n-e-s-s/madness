#include "GroundState.hpp"

#include <madness/chem/SCFOperators.h>
#include <madness/chem/xcfunctional.h>
#include <madness/tensor/tensor_json.hpp>
#include <madness/world/MADworld.h>

#include "broadcast_json.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace molresponse_v3 {

using namespace madness;

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

GroundState::GroundState(World& world, std::shared_ptr<SCF> scf)
    : scf_(std::move(scf)),
      original_k_(FunctionDefaults<3>::get_k()),
      current_k_(FunctionDefaults<3>::get_k()) {
    // from_archive overwrites original_k_/current_k_ with the archive's k right
    // after this; the FunctionDefaults reads here are just safe initializers.
}

GroundState GroundState::from_archive(World& world,
                                       const std::string& archive_path,
                                       const Molecule& molecule) {
    // Step 1: Read archive header to get L, xc, nmo_alpha, spin_restricted
    auto header = read_archive_header(world, archive_path);

    // Step 2: Build CalculationParameters with values from the archive.
    // SCF::load_mos checks L == param.L(), so we must set L correctly.
    // SCF constructor calls set_derived_values(molecule) which computes
    // nalpha/nbeta from nuclear charge + nopen.
    CalculationParameters params;
    params.set_user_defined_value("l", header.L);
    params.set_user_defined_value("xc", header.xc);

    // Derive prefix from archive path: strip ".restartdata" suffix
    std::string prefix = archive_path;
    auto pos = prefix.rfind(".restartdata");
    if (pos != std::string::npos) {
        prefix = prefix.substr(0, pos);
    }
    params.set_user_defined_value("prefix", prefix);
    // Workaround for QCCalculationParametersBase::tostring<std::string>
    // case-folding: re-set the prefix value directly via the QCParameter
    // struct so the original path case is preserved (paths like
    // /home/User/Projects/... would otherwise lowercase, breaking
    // SCF::load_mos's archive lookup).
    params.get_parameter("prefix").set_user_defined_value(prefix);

    if (world.rank() == 0) {
        print("ARCHIVE_HEADER: spin_restricted=", header.spin_restricted,
              " nmo_alpha=", header.nmo_alpha, " L=", header.L,
              " k=", header.k, " xc=", header.xc);
    }

    // Electron bookkeeping (review io finding): set_derived_values computes
    // nalpha/nbeta from Z + charge + nopen, and SCF::load_mos reads EXACTLY
    // param.nmo_alpha() orbitals — with the neutral default a charged species
    // silently DROPS its HOMO (anion) or mis-derives nopen. The archive side
    // knows the truth: closed shell nelec = 2*nmo_alpha exactly; open shell
    // takes nalpha/nbeta from the moldft calc_info.json sibling.
    const double Z = molecule.total_nuclear_charge();
    if (header.spin_restricted) {
        const double charge = Z - 2.0 * static_cast<double>(header.nmo_alpha);
        params.set_user_defined_value("charge", charge);
        if (world.rank() == 0 && std::abs(charge) > 1e-6)
            print("ARCHIVE_HEADER: inferred total charge =", charge,
                  " (nelec =", 2 * header.nmo_alpha, ")");
    } else {
        // rank 0 reads the sibling calc_info.json; broadcast (na, nb, found).
        long na = -1, nb = -1;
        if (world.rank() == 0) {
            const auto ci = std::filesystem::path(prefix + ".calc_info.json");
            std::ifstream in(ci);
            if (in) {
                try {
                    nlohmann::json j;
                    in >> j;
                    na = j.value("calcinfo_nalpha", -1);
                    nb = j.value("calcinfo_nbeta", -1);
                } catch (const std::exception &) { na = nb = -1; }
            }
        }
        std::vector<long> nab{na, nb};
        world.gop.broadcast_serializable(nab, 0);
        na = nab[0]; nb = nab[1];
        if (na > 0 && nb >= 0) {
            const double charge = Z - static_cast<double>(na + nb);
            params.set_user_defined_value("charge", charge);
            params.set_user_defined_value("nopen", static_cast<int>(na - nb));
            if (world.rank() == 0)
                print("ARCHIVE_HEADER: open-shell from calc_info: nalpha=", na,
                      " nbeta=", nb, " nopen=", na - nb,
                      " inferred charge=", charge);
        } else {
            // Legacy fallback: neutral assumption (pre-fix behavior), loudly.
            int nelec = static_cast<int>(Z);
            int nopen = 2 * static_cast<int>(header.nmo_alpha) - nelec;
            if (nopen < 0) nopen = 0;
            params.set_user_defined_value("nopen", nopen);
            if (world.rank() == 0)
                print("ARCHIVE_HEADER: WARNING — no ", prefix,
                      ".calc_info.json next to the archive; ASSUMING NEUTRAL "
                      "(nopen=", nopen, "). A charged open-shell restart needs "
                      "that file (moldft writes it) or the loaded orbital set "
                      "will be wrong.");
        }
    }

    // Step 3: Construct SCF (calls set_derived_values internally)
    auto scf = std::make_shared<SCF>(world, params, molecule);

    // Step 4: Load orbitals via SCF's own archive reader
    scf->load_mos(world);

    // Step 5: Build GroundState
    GroundState gs(world, scf);
    gs.original_k_ = header.k;
    gs.current_k_ = header.k;
    return gs;
}

GroundState::ArchiveHeader
GroundState::read_archive_header(World& world,
                                  const std::string& archive_path) {
    // NB: this parse MIRRORS SCF::save_mos / load_mos's field order
    // (madness/chem/SCF.cc, archive version 4) — it re-reads the header
    // because SCF::load_mos offers no header-only entry point. If moldft's
    // restartdata format drifts, update BOTH in lockstep; the guard below
    // turns silent field-order drift into an actionable error. The header
    // values are broadcast by the parallel archive, so the throw is
    // collective.
    ArchiveHeader h;
    archive::ParallelInputArchive<archive::BinaryFstreamInputArchive>
        ar(world, archive_path.c_str());

    ar & h.version;
    if (h.version != 4) {
        throw std::runtime_error(
            "GroundState::read_archive_header: " + archive_path +
            " has archive version " + std::to_string(h.version) +
            ", but this reader mirrors moldft's version-4 layout "
            "(SCF::save_mos, madness/chem/SCF.cc). Regenerate the ground-state "
            "archive with a matching moldft, or update read_archive_header to "
            "the new field order.");
    }
    ar & h.energy;
    ar & h.spin_restricted;
    ar & h.L;
    ar & h.k;

    // Read molecule (needed to advance archive position, but we discard it
    // since the caller provides the canonical molecule)
    Molecule mol_discard;
    ar & mol_discard;
    ar & h.xc;
    ar & h.localize_method;
    ar & h.converged_for_thresh;

    // Read nmo_alpha (needed to infer nopen for open-shell)
    ar & h.nmo_alpha;

    return h;
}

// ---------------------------------------------------------------------------
// Response-specific preparation
// ---------------------------------------------------------------------------

void GroundState::prepare(World& world, double vtol,
                           const poperatorT& coulop,
                           const std::string& fock_json_file,
                           bool verbose) {
    auto target_k = FunctionDefaults<3>::get_k();
    auto thresh = FunctionDefaults<3>::get_thresh();

    // Reload pristine MOs when the orbitals we hold cannot represent the
    // target protocol: a k change (projection needs the original-k source),
    // OR a TIGHTER thresh than the orbitals currently carry — re-truncating
    // an already-truncated function cannot restore the discarded precision,
    // so a thresh-only climb without the reload silently drags the coarse
    // rung's error into the fine rung (review io finding, verified).
    const bool need_pristine =
        (current_k_ != target_k) ||
        (current_thresh_ >= 0.0 && thresh < current_thresh_);
    if (need_pristine) {
        // Reload pristine MOs at original_k_ from the archive, then project to
        // target. GroundState is only built via from_archive, so the
        // checkpoint to reload from always exists.
        scf_->load_mos(world);
        if (original_k_ != target_k) {
            reconstruct(world, scf_->amo);
            for (auto& orbital : scf_->amo) {
                orbital = project(orbital, target_k, thresh, true);
            }

            if (!is_spin_restricted()) {
                reconstruct(world, scf_->bmo);
                for (auto& orbital : scf_->bmo) {
                    orbital = project(orbital, target_k, thresh, true);
                }
            }
        }
        truncate(world, scf_->amo, thresh);
        if (!is_spin_restricted()) {
            truncate(world, scf_->bmo, thresh);
        }
        current_k_ = target_k;
    } else {
        truncate(world, scf_->amo, thresh);
        if (!is_spin_restricted()) {
            truncate(world, scf_->bmo, thresh);
        }
    }
    current_thresh_ = thresh;

    // Build QProjectors from current orbitals
    q_alpha_ = QProjector<double, 3>(scf_->get_amo());
    if (!is_spin_restricted()) {
        q_beta_ = QProjector<double, 3>(scf_->get_bmo());
    }

    // Initialize SCF's internal state for the current protocol.
    // set_protocol sets FunctionDefaults (k, thresh, cubic_cell),
    // vtol (used by apply_potential for mul_sparse), coulop, gradop, mask.
    // Without this, apply_potential and make_fock_matrix produce wrong results
    // because vtol is uninitialized and the box boundaries may be wrong.
    scf_->set_protocol<3>(world, thresh);

    // Build V_local for the response solver (multiplicative potential)
    build_v_local(world, vtol, coulop);

    // Build Fock matrices using SCF's own methods
    build_fock_matrices(world, vtol, fock_json_file);

    prepared_ = true;

    if (verbose && world.rank() == 0) {   // F2d-ii: suppressed in subworlds
        print_info();
    }
}

void GroundState::build_v_local(World& world, double vtol,
                                 const poperatorT& coulop) {
    // Density: alpha + beta (for restricted, brho = arho so total = 2*arho)
    auto arho = scf_->make_density(world, scf_->aocc, scf_->amo);
    functionT rho;
    if (is_spin_restricted()) {
        rho = 2.0 * arho;
    } else {
        auto brho = scf_->make_density(world, scf_->bocc, scf_->bmo);
        rho = arho + brho;
    }

    // Nuclear potential
    scf_->make_nuclear_potential(world);
    world.gop.fence();
    auto vnuc = scf_->potentialmanager->vnuclear();
    vnuc.truncate(vtol);

    // Coulomb potential from total density
    auto vcoul = apply(*coulop, rho);
    vcoul.truncate(vtol);
    rho.clear();

    // Assemble V_local = V_nuc + V_coul
    v_local_ = vnuc + vcoul;
    vnuc.clear();
    vcoul.clear();

    // Add XC potential if DFT (not pure HF)
    if (scf_->xc.is_dft() && scf_->xc.hf_exchange_coefficient() != 1.0) {
        XCOperator<double, 3> xc_op(world, scf_->param.xc(),
                                     is_spin_restricted(),
                                     arho, arho);
        v_local_ += xc_op.make_xc_potential();
    }
    arho.clear();
}

void GroundState::build_fock_matrices(World& world, double vtol,
                                       const std::string& fock_json_file) {
    auto thresh = FunctionDefaults<3>::get_thresh();
    auto current_k = FunctionDefaults<3>::get_k();

    // Try loading Fock matrices from JSON file first
    bool loaded_from_file = false;
    if (!fock_json_file.empty() && std::filesystem::exists(fock_json_file)) {
        auto fock_json = broadcast_json_file(world, fock_json_file);
        world.gop.fence();

        auto protocol_key = std::string("thresh: ") + std::to_string(thresh)
                          + std::string(" k: ") + std::to_string(current_k);
        if (fock_json.contains(protocol_key)) {
            focka_ = tensor_from_json<double>(fock_json[protocol_key]["focka"]);
            if (!is_spin_restricted() && fock_json[protocol_key].contains("fockb")) {
                fockb_ = tensor_from_json<double>(fock_json[protocol_key]["fockb"]);
            } else {
                fockb_ = copy(focka_);
            }
            loaded_from_file = true;
        }
    }

    if (!loaded_from_file) {
        // Build Fock matrices using SCF's own pipeline:
        //   1. apply_potential: builds V_eff*phi (V_local + V_xc + K)
        //   2. make_fock_matrix: T + <phi|V_eff|phi>, symmetrized

        // SCF::apply_potential needs vlocal = V_nuc + V_coul (without V_xc;
        // apply_potential adds V_xc and K internally).
        // Rebuild the bare vlocal for SCF's interface.
        auto arho = scf_->make_density(world, scf_->aocc, scf_->amo);
        functionT rho;
        if (is_spin_restricted()) {
            rho = 2.0 * arho;
        } else {
            auto brho = scf_->make_density(world, scf_->bocc, scf_->bmo);
            rho = arho + brho;
        }
        arho.clear();

        scf_->make_nuclear_potential(world);
        world.gop.fence();
        auto vnuc = scf_->potentialmanager->vnuclear();
        vnuc.truncate(vtol);

        // coulop already set by set_protocol in prepare()
        auto vcoul = apply(*scf_->coulop, rho);
        vcoul.truncate(vtol);
        rho.clear();

        functionT vlocal_bare = vnuc + vcoul;
        vnuc.clear();
        vcoul.clear();

        // Alpha Fock matrix
        double exca = 0.0, enla = 0.0, ekina = 0.0;
        auto Vpsia = scf_->apply_potential(
            world, scf_->aocc, scf_->amo, vlocal_bare, exca, enla, 0);
        focka_ = scf_->make_fock_matrix(
            world, scf_->amo, Vpsia, scf_->aocc, ekina);
        Vpsia.clear();

        // Beta Fock matrix
        if (!is_spin_restricted() && num_beta() > 0) {
            double excb = 0.0, enlb = 0.0, ekinb = 0.0;
            auto Vpsib = scf_->apply_potential(
                world, scf_->bocc, scf_->bmo, vlocal_bare, excb, enlb, 1);
            fockb_ = scf_->make_fock_matrix(
                world, scf_->bmo, Vpsib, scf_->bocc, ekinb);
            Vpsib.clear();
        } else {
            fockb_ = copy(focka_);
        }
    }

    // Build off-diagonal versions (zero the diagonal)
    focka_no_diag_ = copy(focka_);
    for (long i = 0; i < num_alpha(); i++) {
        focka_no_diag_(i, i) = 0.0;
    }

    fockb_no_diag_ = copy(fockb_);
    for (long i = 0; i < static_cast<long>(fockb_.dim(0)); i++) {
        fockb_no_diag_(i, i) = 0.0;
    }
}

// ---------------------------------------------------------------------------
// Cached operator accessors
// ---------------------------------------------------------------------------

const real_function_3d& GroundState::V_local() const {
    MADNESS_CHECK(prepared_);
    return v_local_;
}

const tensorT& GroundState::focka() const {
    MADNESS_CHECK(prepared_);
    return focka_;
}

const tensorT& GroundState::fockb() const {
    MADNESS_CHECK(prepared_);
    return fockb_;
}

const tensorT& GroundState::focka_no_diag() const {
    MADNESS_CHECK(prepared_);
    return focka_no_diag_;
}

const tensorT& GroundState::fockb_no_diag() const {
    MADNESS_CHECK(prepared_);
    return fockb_no_diag_;
}

const QProjector<double, 3>& GroundState::Q() const {
    MADNESS_CHECK(prepared_);
    return q_alpha_;
}

const QProjector<double, 3>& GroundState::Q_alpha() const {
    MADNESS_CHECK(prepared_);
    return q_alpha_;
}

const QProjector<double, 3>& GroundState::Q_beta() const {
    MADNESS_CHECK(prepared_);
    return q_beta_;
}

// ---------------------------------------------------------------------------
// Info
// ---------------------------------------------------------------------------

void GroundState::print_info() const {
    print("GROUND_STATE_INFO (v3)");
    print("  xc =", scf_->param.xc());
    print("  spin_restricted =", is_spin_restricted());
    print("  num_alpha =", num_alpha());
    if (!is_spin_restricted()) {
        print("  num_beta =", num_beta());
    }
    print("  box_l =", L());
    print("  k =", k());
    print("  orbital_energies_alpha =", scf_->aeps);
    if (!is_spin_restricted() && scf_->beps.size() > 0) {
        print("  orbital_energies_beta =", scf_->beps);
    }
    if (prepared_) {
        print("  focka shape =", focka_.dim(0), "x", focka_.dim(1));
        if (!is_spin_restricted()) {
            print("  fockb shape =", fockb_.dim(0), "x", fockb_.dim(1));
        }
    }
    print("  prepared =", prepared_);
    print("------------------------------------------------------------");
}

} // namespace molresponse_v3
