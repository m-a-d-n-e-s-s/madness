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

/// \file lcao_scan.cc
/// \brief the LCAO state scan (33_state_scan_interface.md)

#include <madness/chem/lcao_scan.h>
#include <madness/tensor/tensor_lapack.h>
#include <madness/world/print.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <map>

namespace madness {
namespace lcao {

namespace {

constexpr double same_energy = 1.e-6;       ///< Eh: solutions closer than this may be one state
constexpr double same_s2 = 1.e-4;           ///< with equal energies: degenerate partners (e.g. pi_x, pi_y)
constexpr double same_overlap = 0.99;       ///< |det C1_occ^T S C2_occ| per spin above this: one state
constexpr double retry_shift = 1.0;         ///< the retry of a start that does not converge (31_, stage a)
constexpr int retry_maxiter = 300;

/// the first nocc columns of C
Tensor<double> occupied(const Tensor<double>& C, const long nocc) {
    return nocc > 0 ? copy(C(_, Slice(0, nocc - 1))) : Tensor<double>();
}

/// C_occ C_occ^T of the first nocc columns
Tensor<double> occ_density(const Tensor<double>& C, const long nocc) {
    if (nocc == 0) return Tensor<double>(C.dim(0), C.dim(0));
    const Tensor<double> Co = occupied(C, nocc);
    return inner(Co, transpose(Co));
}

/// C with the columns i and j exchanged
Tensor<double> swap_columns(const Tensor<double>& C, const long i, const long j) {
    Tensor<double> D = copy(C);
    D(_, i) = C(_, j);
    D(_, j) = C(_, i);
    return D;
}

bool is_index(const std::string& s) {
    return not s.empty() and std::all_of(s.begin(), s.end(), [](const unsigned char c) { return std::isdigit(c); });
}

} // namespace


LCAOStateScan::LCAOStateScan(World& world, const Molecule& molecule, const AtomicBasisSet& aobasis, const int nalpha,
                             const int nbeta, const bool unrestricted, const LCAOParameters& param,
                             const std::string& minbasis, const double charge, const bool collective)
    : world_(world), molecule_(molecule), aobasis_(aobasis), nalpha_(nalpha), nbeta_(nbeta),
      unrestricted_(unrestricted), param_(param), minbasis_(minbasis), charge_(charge),
      scf_(world, molecule, aobasis, nalpha, nbeta, param, collective) {
    scf_.set_unrestricted(unrestricted);
    for (const std::string& s : param_.scan_starts())
        MADNESS_CHECK_THROW(s == "sad" or s == "core" or s == "swaps", "scan_starts: the starts are sad, core and swaps");
    MADNESS_CHECK_THROW(param_.state() == "lowest" or is_index(param_.state()),
                        "state: lowest, or the index of a state in the scan's listing");
    MADNESS_CHECK_THROW(param_.scan_max_states() > 0, "scan_max_states must be positive");
}


double LCAOStateScan::occupied_overlap(const Tensor<double>& C1, const Tensor<double>& C2, const long nocc) const {
    if (nocc == 0) return 1.0;
    const Tensor<double> M = inner(transpose(occupied(C1, nocc)), inner(scf_.overlap(), occupied(C2, nocc)));
    Tensor<double> U, s, VT;
    svd(M, U, s, VT);
    double d = 1.0;
    for (long i = 0; i < s.size(); ++i) d *= s(i);
    return d;
}


void LCAOStateScan::run_start(const std::string& label, const Tensor<double>& Pa, const Tensor<double>& Pb,
                              const bool mom, const Tensor<double>& Ca_occ, const Tensor<double>& Cb_occ) {
    const bool printme = world_.rank() == 0 and param_.print_level() > 0;
    SCFOptions opt = SCFOptions::from(param_);
    opt.econv = param_.scan_econv();
    opt.dconv = param_.scan_dconv();
    opt.print_level = 0;
    opt.mom = mom;
    scf_.iterate(Pa, Pb, opt, Ca_occ, Cb_occ);
    std::string how = label;
    int its = scf_.iterations();
    if (not scf_.converged()) {
        opt.level_shift = std::max(retry_shift, opt.level_shift);
        opt.damping = 0.0;
        opt.maxiter = std::max(retry_maxiter, opt.maxiter);
        scf_.iterate(Pa, Pb, opt, Ca_occ, Cb_occ);
        how += " (retry)";
        its += scf_.iterations();
    }
    if (not scf_.converged()) {
        failed_.push_back(label);
        if (printme) printf("  start %-30s NOT CONVERGED (%d iterations)\n", how.c_str(), its);
        return;
    }
    const LCAOSCF::Result r = scf_.result();
    for (LCAOState& s : states_) {
        if (std::abs(s.result.energies.total - r.energies.total) > same_energy) continue;
        if (occupied_overlap(s.result.Ca, r.Ca, nalpha_) < same_overlap) continue;
        if (nbeta_ > 0 and occupied_overlap(s.result.Cb, r.Cb, nbeta_) < same_overlap) continue;
        s.found_by.push_back(how);
        if (printme)
            printf("  start %-30s %4d its  E %18.10f  <S^2> %9.6f  -> #%d\n", how.c_str(), its,
                   r.energies.total, r.s2, s.id);
        return;
    }
    LCAOState s;
    s.result = r;
    s.found_by = {how};
    s.id = nextid_++;
    states_.push_back(s);
    if (printme)
        printf("  start %-30s %4d its  E %18.10f  <S^2> %9.6f  -> #%d (new)\n", how.c_str(), its, r.energies.total,
               r.s2, s.id);
}


void LCAOStateScan::run() {
    const bool printme = world_.rank() == 0 and param_.print_level() > 0;
    scf_.setup();
    const std::vector<std::string> starts = param_.scan_starts();
    const auto wanted = [&starts](const std::string& s) {
        return std::find(starts.begin(), starts.end(), s) != starts.end();
    };
    if (printme)
        printf("\nLCAO state scan: %s, %d alpha / %d beta electrons; starts:%s%s%s\n",
               scf_.restricted() ? "RHF" : "UHF", nalpha_, nbeta_, wanted("sad") ? " sad" : "",
               wanted("core") ? " core" : "", wanted("swaps") ? " swaps" : "");

    for (const std::string s : {"sad", "core"}) {
        if (not wanted(s)) continue;
        Tensor<double> Pa, Pb;
        const std::string used = scf_.start_densities(s, Pa, Pb);
        run_start(used == s ? s : s + " (as " + used + ")", Pa, Pb, false);
    }

    // occupation swaps of the states of the primary starts, per spin (RHF: of the spatial orbitals)
    if (wanted("swaps")) {
        struct Swap {
            const char* name;
            long hole, particle;      // HOMO - hole -> LUMO + particle
        };
        const Swap swaps[] = {{"H->L", 0, 0}, {"H-1->L", 1, 0}, {"H->L+1", 0, 1}};
        const bool rhf = scf_.restricted();
        const std::size_t nprimary = states_.size();
        for (std::size_t k = 0; k < nprimary; ++k) {
            const LCAOSCF::Result base = states_[k].result;
            const int id = states_[k].id;
            for (int spin = 0; spin < (rhf ? 1 : 2); ++spin) {
                const long nocc = (spin == 0) ? nalpha_ : nbeta_;
                const Tensor<double>& C = (spin == 0) ? base.Ca : base.Cb;
                if (nocc == 0) continue;
                for (const Swap& sw : swaps) {
                    const long from = nocc - 1 - sw.hole, to = nocc + sw.particle;
                    if (from < 0 or to >= C.dim(1)) continue;
                    Tensor<double> Ca = base.Ca, Cb = base.Cb;
                    (spin == 0 ? Ca : Cb) = swap_columns(C, from, to);
                    if (rhf) Cb = Ca;
                    const std::string label = std::string("swap ") + (rhf ? "" : (spin == 0 ? "a " : "b ")) + sw.name +
                                              " of #" + std::to_string(id);
                    run_start(label, occ_density(Ca, nalpha_), occ_density(Cb, nbeta_), true, occupied(Ca, nalpha_),
                              occupied(Cb, nbeta_));
                }
            }
        }
    }
    MADNESS_CHECK_THROW(not states_.empty(), "LCAO state scan: no start converged");

    // the listing: ascending energy, at most scan_max_states
    std::stable_sort(states_.begin(), states_.end(), [](const LCAOState& a, const LCAOState& b) {
        return a.result.energies.total < b.result.energies.total;
    });
    if (states_.size() > std::size_t(param_.scan_max_states())) states_.resize(param_.scan_max_states());

    // IAO charges and spin populations of each state (analytic overlaps, every rank)
    for (LCAOState& s : states_) {
        try {
            s.populations = population_analysis(molecule_, aobasis_, occupied(s.result.Ca, nalpha_),
                                                occupied(s.result.Cb, nbeta_), {"iao"}, minbasis_, minbasis_,
                                                charge_, 0);
        } catch (const std::exception& e) {
            s.populations = nlohmann::json{{"error", std::string(e.what())}};
        }
    }

    // the seed
    if (param_.state() == "lowest") {
        seed_ = 0;
    } else {
        seed_ = std::stol(param_.state());
        MADNESS_CHECK_THROW(seed_ < long(states_.size()), "state: the scan's listing has no state with this index");
    }
    scf_.restore(states_[seed_].result);
}


void LCAOStateScan::print() const {
    if (world_.rank() != 0) return;
    const bool rhf = scf_.restricted();
    std::map<int, long> index;      // discovery id -> listing index
    for (std::size_t i = 0; i < states_.size(); ++i) index[states_[i].id] = long(i);
    const auto relabel = [&index](std::string s) {     // "#id" -> "state <index>" in the found-by labels
        const std::size_t p = s.find('#');
        if (p == std::string::npos) return s;
        std::size_t q = p + 1;
        while (q < s.size() and std::isdigit(static_cast<unsigned char>(s[q]))) ++q;
        const int id = std::stoi(s.substr(p + 1, q - p - 1));
        const auto it = index.find(id);
        return s.substr(0, p) + (it == index.end() ? "a dropped state" : "state " + std::to_string(it->second)) +
               s.substr(q);
    };
    printf("\nLCAO state scan: %zu distinct states (%s, %ld basis functions, %d alpha / %d beta)\n", states_.size(),
           rhf ? "RHF" : "UHF", scf_.nbf(), nalpha_, nbeta_);
    printf("  converged to scan_econv %.0e, scan_dconv %.0e; one state: |dE| < %.0e Eh and |det C1_occ^T S C2_occ| "
           "> %.2f per spin\n", param_.scan_econv(), param_.scan_dconv(), same_energy, same_overlap);
    if (not failed_.empty()) {
        printf("  starts that did not converge:");
        for (const std::string& f : failed_) printf(" %s;", relabel(f).c_str());
        printf("\n");
    }
    printf("\n state        energy (Eh)    dE (mEh)      <S^2>   HOMO a/b (Eh)        LUMO a/b (Eh)        found by\n");
    const double e0 = states_.empty() ? 0.0 : states_[0].result.energies.total;
    for (std::size_t i = 0; i < states_.size(); ++i) {
        const LCAOSCF::Result& r = states_[i].result;
        const auto eps = [](const Tensor<double>& e, const long k) { return (k >= 0 and k < e.size()) ? e(k) : NAN; };
        std::string by;
        for (const std::string& f : states_[i].found_by) by += (by.empty() ? "" : ", ") + relabel(f);
        const long partner = degenerate_partner(i);
        if (partner >= 0) by = "degenerate with state " + std::to_string(partner) + "; " + by;
        printf("%5zu  %18.10f  %9.3f  %9.6f   %8.4f/%8.4f    %8.4f/%8.4f    %s%s\n", i, r.energies.total,
               (r.energies.total - e0) * 1e3, r.s2, eps(r.epsa, nalpha_ - 1), eps(r.epsb, nbeta_ - 1),
               eps(r.epsa, nalpha_), eps(r.epsb, nbeta_), by.c_str(), long(i) == seed_ ? "   <- seed" : "");
    }
    // IAO charges, and spin populations for open shells
    printf("\n IAO (minimal basis %s): charge%s per atom\n", minbasis_.c_str(), rhf ? "" : " / spin population");
    for (std::size_t i = 0; i < states_.size(); ++i) {
        const nlohmann::json& p = states_[i].populations;
        if (not p.contains("schemes") or not p["schemes"].contains("iao")) {
            printf("%5zu  (no populations: %s)\n", i, p.dump().c_str());
            continue;
        }
        const nlohmann::json& iao = p["schemes"]["iao"];
        const std::vector<std::string> atoms = p["atoms"].get<std::vector<std::string>>();
        const std::vector<double> q = iao["charges"].get<std::vector<double>>();
        // spin populations of every UHF state, also of broken-symmetry ones with nalpha == nbeta (where the JSON
        // has no "spin"): the alpha minus the beta electrons per atom
        std::vector<double> sp;
        if (not rhf) {
            const std::vector<double> ea = iao["electrons_alpha"].get<std::vector<double>>();
            const std::vector<double> eb = iao["electrons_beta"].get<std::vector<double>>();
            for (std::size_t a = 0; a < ea.size(); ++a) sp.push_back(ea[a] - eb[a]);
        }
        printf("%5zu ", i);
        for (std::size_t a = 0; a < atoms.size(); ++a) {
            if (a > 0 and a % 6 == 0) printf("\n      ");
            if (sp.empty())
                printf("  %s%zu %+.3f", atoms[a].c_str(), a, q[a]);
            else
                printf("  %s%zu %+.3f/%+.3f", atoms[a].c_str(), a, q[a], sp[a]);
        }
        printf("\n");
    }
    printf("\n");
}


long LCAOStateScan::degenerate_partner(const std::size_t i) const {
    for (std::size_t j = 0; j < i; ++j)
        if (std::abs(states_[i].result.energies.total - states_[j].result.energies.total) < same_energy and
            std::abs(states_[i].result.s2 - states_[j].result.s2) < same_s2)
            return long(j);
    return -1;
}


nlohmann::json LCAOStateScan::to_json() const {
    std::map<int, long> index;
    for (std::size_t i = 0; i < states_.size(); ++i) index[states_[i].id] = long(i);
    nlohmann::json j;
    j["method"] = scf_.restricted() ? "RHF" : "UHF";
    j["nalpha"] = nalpha_;
    j["nbeta"] = nbeta_;
    j["nbf"] = scf_.nbf();
    j["basis"] = param_.basis();
    j["starts"] = param_.scan_starts();
    j["scan_econv"] = param_.scan_econv();
    j["scan_dconv"] = param_.scan_dconv();
    j["seed"] = seed_;
    j["failed_starts"] = failed_;
    j["states"] = nlohmann::json::array();
    for (std::size_t i = 0; i < states_.size(); ++i) {
        const LCAOSCF::Result& r = states_[i].result;
        const auto eps = [](const Tensor<double>& e, const long k) { return (k >= 0 and k < e.size()) ? e(k) : 0.0; };
        nlohmann::json s;
        s["index"] = i;
        s["energy"] = r.energies.total;
        s["s2"] = r.s2;
        s["iterations"] = r.iterations;
        s["homo"] = {eps(r.epsa, nalpha_ - 1), eps(r.epsb, nbeta_ - 1)};
        s["lumo"] = {eps(r.epsa, nalpha_), eps(r.epsb, nbeta_)};
        s["found_by"] = states_[i].found_by;
        s["discovery_id"] = states_[i].id;
        s["degenerate_with"] = degenerate_partner(i);
        s["populations"] = states_[i].populations;
        j["states"].push_back(s);
    }
    return j;
}

} // namespace lcao
} // namespace madness
