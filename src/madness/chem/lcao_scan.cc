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
#include <madness/constants.h>
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
    for (const std::string& s : param_.scan_starts()) {
        MADNESS_CHECK_THROW(s == "sad" or s == "core" or s == "swaps" or s == "bs" or s == "flip",
                            "scan_starts: the starts are sad, core, swaps, bs and flip");
        MADNESS_CHECK_THROW(s != "bs" or (nalpha_ == nbeta_ and unrestricted_),
                            "scan_starts bs: needs nalpha == nbeta (nopen 0) and spin_restricted false");
        MADNESS_CHECK_THROW(s != "flip" or (unrestricted_ and nbeta_ > 0 and not param_.flip_atoms().empty()),
                            "scan_starts flip: needs spin_restricted false, beta electrons and flip_atoms");
    }
    for (const int a : param_.flip_atoms())
        MADNESS_CHECK_THROW(a >= 0 and std::size_t(a) < molecule_.natom(), "flip_atoms: an atom index out of range");
    MADNESS_CHECK_THROW(param_.state() == "lowest" or is_index(param_.state()),
                        "state: lowest, or the index of a state in the scan's listing");
    MADNESS_CHECK_THROW(param_.scan_max_states() > 0, "scan_max_states must be positive");
}


void LCAOStateScan::flip_start() {
    const bool printme = world_.rank() == 0 and param_.print_level() > 0;
    const std::vector<int> atoms = param_.flip_atoms();
    std::string label = "flip of atoms";
    for (const int a : atoms) label += " " + std::to_string(a);
    // the high-spin state, with the same integrals
    scf_.set_occupations(nalpha_ + 1, nbeta_ - 1);
    Tensor<double> Pa, Pb;
    scf_.start_densities("sad", Pa, Pb);
    SCFOptions opt = SCFOptions::from(param_);
    opt.econv = param_.scan_econv();
    opt.dconv = param_.scan_dconv();
    opt.print_level = 0;
    opt.level_shift = std::max(retry_shift, opt.level_shift);
    opt.maxiter = std::max(retry_maxiter, opt.maxiter);
    scf_.iterate(Pa, Pb, opt);
    const LCAOSCF::Result hs = scf_.result();
    scf_.set_occupations(nalpha_, nbeta_);
    if (printme)
        printf("  high-spin state (%d alpha / %d beta) for %s: %s, E %18.10f  <S^2> %9.6f\n", nalpha_ + 1,
               nbeta_ - 1, label.c_str(), hs.converged ? "converged" : "NOT CONVERGED", hs.energies.total, hs.s2);

    // the fragment's share of an orbital: its Mulliken population on the basis functions of flip_atoms
    const Tensor<double>& S = scf_.overlap();
    std::vector<bool> on(S.dim(0), false);
    for (const Shell& sh : scf_.shells())
        if (std::find(atoms.begin(), atoms.end(), sh.atom) != atoms.end())
            for (int c = 0; c < sh.ncart(); ++c) on[sh.offset + c] = true;
    const auto share = [&](const Tensor<double>& psi) {
        const Tensor<double> Spsi = inner(S, psi);
        double p = 0.0;
        for (long mu = 0; mu < psi.size(); ++mu)
            if (on[mu]) p += psi(mu) * Spsi(mu);
        return p;
    };
    // the two highest alpha orbitals of the high-spin state, rotated to separate their shares (a fine grid)
    const long n1 = nalpha_ - 1, n2 = nalpha_;
    const Tensor<double> phi1 = copy(hs.Ca(_, n1)), phi2 = copy(hs.Ca(_, n2));
    double best = -1.0, theta = 0.0;
    for (int k = 0; k < 720; ++k) {
        const double th = constants::pi * k / 720.0;
        const double d = share(std::cos(th) * phi1 + std::sin(th) * phi2) -
                         share(-std::sin(th) * phi1 + std::cos(th) * phi2);
        if (std::abs(d) > best) {
            best = std::abs(d);
            theta = th;
        }
    }
    Tensor<double> psi1 = std::cos(theta) * phi1 + std::sin(theta) * phi2;
    Tensor<double> psi2 = -std::sin(theta) * phi1 + std::cos(theta) * phi2;
    if (share(psi2) > share(psi1)) std::swap(psi1, psi2);      // psi1: the one on flip_atoms
    if (printme)
        printf("  flip: the two highest alpha orbitals have %.3f and %.3f of their population on the atoms\n",
               share(psi1), share(psi2));
    // alpha: the high-spin alpha orbitals without psi1; beta: the high-spin beta orbitals and psi1, orthonormalized
    const long n = S.dim(0);
    Tensor<double> Ca(n, nalpha_), Cb(n, nbeta_);
    for (long i = 0; i < nalpha_ - 1; ++i) Ca(_, i) = hs.Ca(_, i);
    Ca(_, nalpha_ - 1) = psi2;
    for (long i = 0; i < nbeta_ - 1; ++i) Cb(_, i) = hs.Cb(_, i);
    Cb(_, nbeta_ - 1) = psi1;
    const auto orthonormal = [&S](const Tensor<double>& C) {
        Tensor<double> s = inner(transpose(C), inner(S, C));
        Tensor<double> U, e;
        syev(0.5 * (s + transpose(s)), U, e);
        Tensor<double> Ui = copy(U);
        for (long j = 0; j < e.size(); ++j)
            for (long i = 0; i < Ui.dim(0); ++i) Ui(i, j) /= std::sqrt(e(j));
        return Tensor<double>(inner(C, inner(Ui, transpose(U))));
    };
    Ca = orthonormal(Ca);
    Cb = orthonormal(Cb);
    run_start(label, inner(Ca, transpose(Ca)), inner(Cb, transpose(Cb)), true, Ca, Cb);
}


StabilityBlock LCAOStateScan::followed_block() const {
    return scf_.restricted() ? StabilityBlock::rhf_rhf : StabilityBlock::uhf_uhf;
}


void LCAOStateScan::analyze(LCAOState& s) const {
    s.stability.clear();
    const std::vector<StabilityBlock> blocks = scf_.restricted()
            ? std::vector<StabilityBlock>{StabilityBlock::rhf_rhf, StabilityBlock::rhf_uhf}
            : std::vector<StabilityBlock>{StabilityBlock::uhf_uhf};
    for (const StabilityBlock b : blocks) {
        const OrbitalHessian H(scf_, s.result, nalpha_, nbeta_, b);
        s.stability.push_back(lowest_roots(H, param_.stability_roots()));
    }
    s.stable = true;
    for (const StabilityRoots& r : s.stability)
        if (r.block == followed_block() and not r.eigenvalues.empty() and r.eigenvalues[0] < -param_.stability_tol())
            s.stable = false;
}


void LCAOStateScan::follow(const std::size_t k) {
    const LCAOState s = states_[k];         // a copy: run_start may grow states_
    const StabilityRoots* root = nullptr;
    for (const StabilityRoots& r : s.stability)
        if (r.block == followed_block()) root = &r;
    if (root == nullptr or root->eigenvalues.empty()) return;
    const bool rhf = scf_.restricted();
    const auto densities = [&](const double t, Tensor<double>& Pa, Tensor<double>& Pb) {
        const Tensor<double> Ca = rotate_occupied(s.result.Ca, nalpha_, root->xa[0], t);
        Pa = inner(Ca, transpose(Ca));
        if (rhf) {
            Pb = Pa;
        } else if (nbeta_ > 0) {
            const Tensor<double> Cb = rotate_occupied(s.result.Cb, nbeta_, root->xb[0], t);
            Pb = inner(Cb, transpose(Cb));
        } else {
            Pb = Tensor<double>(Pa.dim(0), Pa.dim(1));
        }
    };
    // the step along the mode (the eigenvector has norm 1): the lowest energy of a few, then reconverge
    double tbest = 0.0, ebest = s.result.energies.total;
    for (const double t : {0.1, 0.2, 0.4, 0.7, 1.0}) {
        Tensor<double> Pa, Pb;
        densities(t, Pa, Pb);
        const double e = scf_.energy(Pa, Pb);
        if (e < ebest) {
            ebest = e;
            tbest = t;
        }
    }
    if (tbest == 0.0) tbest = 0.2;          // no lower energy on the line: reconverge from a small step anyway
    Tensor<double> Pa, Pb;
    densities(tbest, Pa, Pb);
    char label[128];
    snprintf(label, sizeof(label), "stability of #%d (%s %+.4f, step %.1f)", s.id, to_string(root->block).c_str(),
             root->eigenvalues[0], tbest);
    run_start(label, Pa, Pb, false);
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
    if (printme) {
        printf("\nLCAO state scan: %s, %d alpha / %d beta electrons; starts:", scf_.restricted() ? "RHF" : "UHF",
               nalpha_, nbeta_);
        for (const std::string& s : starts) printf(" %s", s.c_str());
        printf("\n");
    }

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
    // broken symmetry: the alpha/beta HOMO-LUMO mix of the states of the primary starts (alpha +45 degrees,
    // beta -45 degrees), reconverged; nalpha == nbeta, UHF
    if (wanted("bs")) {
        const std::size_t nprimary = std::min(states_.size(), std::size_t(2));
        for (std::size_t k = 0; k < nprimary; ++k) {
            const LCAOSCF::Result base = states_[k].result;
            if (base.Ca.dim(1) <= nalpha_) continue;
            const double c = std::sqrt(0.5);
            const auto mixed = [&](const Tensor<double>& C, const long nocc, const double sign) {
                Tensor<double> Co = occupied(C, nocc);
                Co(_, nocc - 1) = c * C(_, nocc - 1) + sign * c * C(_, nocc);
                return Co;
            };
            const Tensor<double> Coa = mixed(base.Ca, nalpha_, 1.0), Cob = mixed(base.Cb, nbeta_, -1.0);
            run_start("bs of #" + std::to_string(states_[k].id), inner(Coa, transpose(Coa)),
                      inner(Cob, transpose(Cob)), false);
        }
    }

    // a defined broken-symmetry state: the high-spin state (nalpha+1, nbeta-1), its two highest alpha orbitals
    // rotated so that one lies on flip_atoms as much as possible (Mulliken population of its basis functions),
    // that one moved to beta, then reconverged with the occupation held by maximum overlap
    if (wanted("flip")) flip_start();

    MADNESS_CHECK_THROW(not states_.empty(), "LCAO state scan: no start converged");

    // stability: every state, also those that following finds (states_ grows); an unstable one is followed
    if (param_.stability()) {
        std::size_t follows = 0;
        const std::size_t max_follows = 2 * states_.size() + 6;
        for (std::size_t k = 0; k < states_.size(); ++k) {
            analyze(states_[k]);
            if (printme) {
                printf("  stability of #%d:", states_[k].id);
                for (const StabilityRoots& r : states_[k].stability)
                    printf("  %s %+.5f%s", to_string(r.block).c_str(), r.eigenvalues.empty() ? NAN : r.eigenvalues[0],
                           r.converged ? "" : " (Davidson not converged)");
                printf("  -> %s\n", states_[k].stable ? "stable" : "unstable");
            }
            if (not states_[k].stable and follows < max_follows) {
                ++follows;
                follow(k);
            }
        }
    }

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
        while (seed_ < long(states_.size()) and not states_[seed_].stable) ++seed_;
        if (seed_ == long(states_.size())) {
            seed_ = 0;
            if (printme) madness::print("WARNING: no listed state is stable; the seed is the lowest state");
        }
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
    printf("\n state        energy (Eh)    dE (mEh)      <S^2>   HOMO a/b (Eh)        LUMO a/b (Eh)        "
           "lowest Hessian eigenvalue(s)        found by\n");
    const double e0 = states_.empty() ? 0.0 : states_[0].result.energies.total;
    for (std::size_t i = 0; i < states_.size(); ++i) {
        const LCAOSCF::Result& r = states_[i].result;
        const auto eps = [](const Tensor<double>& e, const long k) { return (k >= 0 and k < e.size()) ? e(k) : NAN; };
        std::string by;
        for (const std::string& f : states_[i].found_by) by += (by.empty() ? "" : ", ") + relabel(f);
        const long partner = degenerate_partner(i);
        if (partner >= 0) by = "degenerate with state " + std::to_string(partner) + "; " + by;
        std::string hess;
        for (const StabilityRoots& sr : states_[i].stability) {
            char buf[64];
            snprintf(buf, sizeof(buf), "%s%s %+.5f", hess.empty() ? "" : " ", to_string(sr.block).c_str(),
                     sr.eigenvalues.empty() ? NAN : sr.eigenvalues[0]);
            hess += buf;
        }
        if (hess.empty()) hess = "-";
        else hess += states_[i].stable ? "" : " UNSTABLE";
        printf("%5zu  %18.10f  %9.3f  %9.6f   %8.4f/%8.4f    %8.4f/%8.4f    %-34s %s%s\n", i, r.energies.total,
               (r.energies.total - e0) * 1e3, r.s2, eps(r.epsa, nalpha_ - 1), eps(r.epsb, nbeta_ - 1),
               eps(r.epsa, nalpha_), eps(r.epsb, nbeta_), hess.c_str(), by.c_str(),
               long(i) == seed_ ? "   <- seed" : "");
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
    j["stability"] = param_.stability();
    j["stability_tol"] = param_.stability_tol();
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
        s["stable"] = states_[i].stable;
        s["stability"] = nlohmann::json::array();
        for (const StabilityRoots& sr : states_[i].stability)
            s["stability"].push_back({{"block", to_string(sr.block)}, {"eigenvalues", sr.eigenvalues},
                                      {"converged", sr.converged}, {"products", sr.products}});
        s["populations"] = states_[i].populations;
        j["states"].push_back(s);
    }
    return j;
}

} // namespace lcao
} // namespace madness
