// seed_moldft_from_dalton — build a moldft ground-state SEED by projecting the
// occupied DALTON/molden MOs onto the MRA basis and writing a version-4
// `.restartdata` archive that moldft `--restart` resumes from. Ground-state
// analogue of tools/seed_from_dalton (which seeds the ES/response bundle).
//
// The write layout mirrors SCF::save_mos (chem/SCF.cc, version 4) EXACTLY so
// SCF::load_mos accepts it with no moldft change:
//   uint version; double energy; bool spin_restricted;
//   double L; int k; Molecule; string xc; string localize; double conv_thresh;
//   uint nmo; Tensor eps; Tensor occ; vector<int> set; Function mo[nmo];
// (closed-shell only: spin_restricted=true, no beta block.)
//
// Constraints / couplings:
//   * NP=1 only (matches seed_from_dalton; single-client projection).
//   * L MUST equal the target moldft run's box (load_mos THROWS otherwise);
//     pass --L to match the deck's `l` (default 200, the value used here).
//   * k/thresh mismatch is handled by load_mos (it re-projects), so seeding at
//     the first protocol rung (k6/1e-4) is fine even for a 1e-4,1e-6 climb.
//   * xc/localize are LOADED INTO the run's params by load_mos, so they must be
//     what the run intends (defaults hf / canon).
//   * MADNESS re-converges to its own MRA-HF minimum; the DALTON orbitals are a
//     starting guess, not frozen.
//
// Usage:
//   seed_moldft_from_dalton --molden=molden.inp --n-occ=N --out=PREFIX
//                           [--L=200] [--thresh=1e-4] [--energy=<E_h>]
//                           [--xc=hf] [--localize=canon] [--nio=1]
//   -> writes PREFIX.restartdata ; run:  moldft (with `restart 1`, prefix PREFIX)

#include "dalton_gto.hpp"

#include "../ResponseProtocol.hpp"

#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>
#include <madness/world/parallel_archive.h>
#include <madness/chem/molecule.h>

#include <cctype>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

namespace {

// Projects one DALTON MO (its AO-coefficient column) as a real MRA function.
// Same evaluation path as seed_from_dalton's DaltonResponseFunctor.
class DaltonMOFunctor : public madness::FunctionFunctorInterface<double, 3> {
    const DaltonMoldenBasis& basis;
    std::vector<double> c;                       // AO coeffs of this MO
    std::vector<madness::coord_3d> centers;
public:
    DaltonMOFunctor(const DaltonMoldenBasis& b, std::vector<double> w)
        : basis(b), c(std::move(w)) {
        for (const auto& sh : basis.shells) centers.push_back({sh.cx, sh.cy, sh.cz});
    }
    double operator()(const madness::coord_3d& r) const override {
        double val = 0.0, bf[9];
        for (size_t s = 0; s < basis.shells.size(); s++) {
            const auto& sh = basis.shells[s];
            sh.evaluate(r[0], r[1], r[2], bf);
            const int off = basis.ao_offsets[s];
            for (int k = 0; k < sh.n_ao; k++)
                val += c[static_cast<size_t>(off + k)] * bf[k];
        }
        return val;
    }
    std::vector<madness::coord_3d> special_points() const override { return centers; }
};

// Minimal symbol -> atomic number for the elements in the current study set.
// DALTON molden labels duplicate atoms with an index suffix (H1, H2, O, ...),
// so strip trailing non-alphabetic characters to recover the element.
int symbol_to_Z(const std::string& s) {
    std::string elem;
    for (char c : s) { if (std::isalpha(static_cast<unsigned char>(c))) elem += c; else break; }
    static const std::vector<std::pair<std::string, int>> tab = {
        {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5}, {"C", 6},
        {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10}, {"Na", 11}, {"Mg", 12},
        {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17}, {"Ar", 18}};
    for (const auto& [sym, z] : tab) if (sym == elem) return z;
    throw std::runtime_error("seed_moldft_from_dalton: unknown element symbol '" + s +
                             "' (parsed element '" + elem + "')");
}

}  // namespace

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    startup(world, argc, argv, true);
    commandlineparser parser(argc, argv);

    if (world.size() != 1) {
        if (world.rank() == 0)
            print("ERROR: seed_moldft_from_dalton must run on a single rank (NP=1).");
        finalize();
        return 2;
    }
    if (!parser.key_exists("molden") || !parser.key_exists("n-occ") ||
        !parser.key_exists("out")) {
        if (world.rank() == 0) {
            print("Usage: seed_moldft_from_dalton --molden=molden.inp --n-occ=N "
                  "--out=PREFIX");
            print("  [--L=200] [--thresh=1e-4] [--energy=<E_h>] [--xc=hf] "
                  "[--localize=canon] [--nio=1]");
        }
        finalize();
        return 2;
    }

    const std::string molden_path = parser.value_raw("molden");
    const std::string out_prefix  = parser.value_raw("out");
    const int    n_occ  = std::stoi(parser.value("n-occ"));
    const double L      = parser.key_exists("L") ? std::stod(parser.value("L")) : 200.0;
    const double thresh = parser.key_exists("thresh") ? std::stod(parser.value("thresh")) : 1e-4;
    const double energy = parser.key_exists("energy") ? std::stod(parser.value("energy")) : 0.0;
    const std::string xc  = parser.key_exists("xc") ? parser.value("xc") : "hf";
    const std::string loc = parser.key_exists("localize") ? parser.value("localize") : "canon";
    const int    nio    = parser.key_exists("nio") ? std::stoi(parser.value("nio")) : 1;
    const int    k      = default_k_for_thresh(thresh);

    {
        DaltonMoldenResult molden = read_molden(molden_path);
        const int n_ao = molden.n_ao;
        if (n_occ <= 0 || n_occ > molden.n_mo) {
            if (world.rank() == 0)
                print("ERROR: n_occ=", n_occ, " out of range (n_mo=", molden.n_mo, ")");
            finalize();
            return 2;
        }

        // Molecule from the molden geometry (Bohr); q = atn = Z (all-electron).
        Molecule molecule;
        for (size_t a = 0; a < molden.atom_symbols.size(); ++a) {
            const int Z = symbol_to_Z(molden.atom_symbols[a]);
            molecule.add_atom(molden.coords[a][0], molden.coords[a][1],
                              molden.coords[a][2], static_cast<double>(Z), Z);
        }

        // Discretization BEFORE projection (cell must match the run's L).
        Tensor<double> cell(3L, 2L);
        for (int i = 0; i < 3; i++) { cell(i, 0) = -L; cell(i, 1) = L; }
        FunctionDefaults<3>::set_cell(cell);
        FunctionDefaults<3>::set_k(k);
        FunctionDefaults<3>::set_thresh(thresh);

        // Project the n_occ occupied MOs (columns 0..n_occ-1 of the column-major
        // C[mu + n_ao*mo] coefficient matrix).
        std::vector<real_function_3d> amo;
        for (int i = 0; i < n_occ; i++) {
            std::vector<double> col(
                molden.mo_coeffs.begin() + static_cast<ptrdiff_t>(i * n_ao),
                molden.mo_coeffs.begin() + static_cast<ptrdiff_t>((i + 1) * n_ao));
            std::shared_ptr<madness::FunctionFunctorInterface<double, 3>> ff =
                std::make_shared<DaltonMOFunctor>(molden.basis, std::move(col));
            amo.push_back(FunctionFactory<double, 3>(world)
                              .functor(ff).thresh(thresh).truncate_on_project());
        }

        // DALTON MOs are orthonormal in the GTO metric, but independent MRA
        // projection leaves small mutual overlaps. moldft's .restartdata path
        // (SCF::load_mos) assumes an orthonormal set and does NOT re-orthonormalize
        // on load (a normal moldft save already writes one), so without this the
        // first density is invalid and the SCF collapses (a -106.9 Ha variational
        // collapse was observed on water without it). Symmetric (Löwdin)
        // orthonormalization cleans the set while preserving orbital character.
        double max_offdiag = 0.0;
        {
            Tensor<double> S = matrix_inner(world, amo, amo, true);
            for (int i = 0; i < n_occ; i++)
                for (int j = 0; j < n_occ; j++)
                    if (i != j) max_offdiag = std::max(max_offdiag, std::abs(S(i, j)));
        }
        amo = orthonormalize_symmetric(amo);

        // Orbital energies from the molden; occupations follow moldft's
        // ALPHA-occupation convention: aocc = 1.0 per occupied spatial orbital
        // (SCF::make_density uses aocc, then spin_restricted DOUBLES the density,
        // SCF.cc:522-527). The molden's chemist-convention occ=2.0 would double
        // the density -> 20 electrons, 2x every energy term, variational collapse
        // (observed -106.9 Ha on water before this fix).
        Tensor<double> aeps(static_cast<long>(n_occ)), aocc(static_cast<long>(n_occ));
        double smo_norm_dev = 0.0;
        for (int i = 0; i < n_occ; i++) {
            aeps(i) = (i < static_cast<int>(molden.mo_energies.size())) ? molden.mo_energies[i] : 0.0;
            aocc(i) = 1.0;   // moldft doubles for spin_restricted closed shell
            smo_norm_dev = std::max(smo_norm_dev, std::abs(amo[i].norm2() - 1.0));
        }
        std::vector<int> aset(static_cast<size_t>(n_occ), 0);

        // Write the version-4 restartdata (SCF::save_mos layout, closed-shell).
        const bool spin_restricted = true;
        const double conv_thresh = thresh;
        const std::string archivename = out_prefix + ".restartdata";
        archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(
            world, archivename.c_str(), nio);
        unsigned int version = 4;
        ar & version;
        ar & energy & spin_restricted;
        ar & L & k & molecule & xc & loc & conv_thresh;
        ar & static_cast<unsigned int>(amo.size());
        ar & aeps & aocc & aset;
        for (unsigned int i = 0; i < amo.size(); ++i) ar & amo[i];

        if (world.rank() == 0) {
            print("seed_moldft_from_dalton: wrote", n_occ, "occupied MOs ->", archivename);
            print("  molden   =", molden_path, "  n_ao =", n_ao, "  n_mo =", molden.n_mo);
            print("  L =", L, "  k =", k, "  thresh =", thresh, "  xc =", xc,
                  "  localize =", loc, "  energy =", energy);
            print("  atoms    =", static_cast<int>(molden.atom_symbols.size()),
                  "  max |S_ij| (i!=j) pre-orthonormalization =", max_offdiag,
                  "  |proj MO norm - 1| max (post) =", smo_norm_dev);
            print("  Run:  moldft  (input with `restart 1`, `prefix ", out_prefix,
                  "`, `l ", L, "`, `xc ", xc, "`) -> SCF resumes from seed.");
        }
    }  // functions destruct before finalize

    finalize();
    return 0;
}
