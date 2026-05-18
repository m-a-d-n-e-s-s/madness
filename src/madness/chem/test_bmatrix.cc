// test_bmatrix.cc — B matrix comparison: commutator approach vs LRF-B
//
// Runs a Nemo SCF from an input file, builds the Info struct via
// CCPotentials::update_info, then compares two ways to compute the B matrix:
//   (a) Commutator:  B = B^Ue + B^g + B^KffK + B^X  (exact 3D, thresh_6D-dependent)
//   (b) LRF-B:       orbital-weighted LRF of f12*phi_i*phi_j  (independent of thresh_6D)
//
// Loads MP2 pair amplitudes from CCPair archives in the working directory.
// Tries (in order):
//   iter_final.MP2_pair_u_{ij}   — converged pair amplitude
//   MP2_pair_u_{ij}_const        — CCPair with constant source term
// Falls back to u_ij = 0 if neither is found (smoke test).
//
// Usage:
//   test_bmatrix --input <inputfile>

#include <fstream>
#include <madness.h>
#include <madness/chem/nemo.h>
#include <madness/chem/CCPotentials.h>
#include <madness/chem/CCParameters.h>
#include <madness/chem/CCStructures.h>
#include <madness/chem/localizer.h>
#include <madness/chem/b_matrix.h>
#include <madness/chem/lowrankfunction.h>
#include <madness/mra/commandlineparser.h>

using namespace madness;

int main(int argc, char** argv) {
    World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);
    commandlineparser parser(argc, argv);

    int status = 0;
#ifdef ENABLE_GENTENSOR
    try {
        // ── SCF ──────────────────────────────────────────────────────────
        std::shared_ptr<Nemo> nemo(new Nemo(world, parser));
        nemo->value();

        auto cell = FunctionDefaults<3>::get_cell();
        FunctionDefaults<6>::set_cubic_cell(cell(0,0), cell(0,1));
        FunctionDefaults<6>::set_tensor_type(TT_2D);

        // ── CC setup ──────────────────────────────────────────────────────
        CCParameters ccparam(world, parser);
        FunctionDefaults<3>::set_thresh(ccparam.thresh_3D());
        FunctionDefaults<6>::set_thresh(ccparam.thresh_6D());
        FunctionDefaults<6>::set_k(FunctionDefaults<3>::get_k());
        FunctionDefaults<3>::set_truncate_mode(3);
        FunctionDefaults<6>::set_truncate_mode(3);

        Tensor<double> fmat = nemo->compute_fock_matrix(nemo->get_calc()->amo, nemo->get_calc()->aocc);
        long nfrozen = Localizer::determine_frozen_orbitals(fmat);
        ccparam.set_derived_value<long>("freeze", nfrozen);

        CCPotentials CCOPS(world, nemo, ccparam);
        CCOPS.reset_nemo(nemo);
        Info info = CCOPS.update_info(ccparam, nemo);

        const size_t freeze = info.parameters.freeze();
        const size_t nocc   = info.mo_ket.size();

        if (world.rank() == 0) {
            print_header1("B Matrix Comparison: Commutator vs LRF-B");
            print("nocc =", nocc, " freeze =", freeze);
            print("thresh_6D =", ccparam.thresh_6D());
            print("orbital energies:", info.orbital_energies);
        }

        // ── Hylleraas functional ──────────────────────────────────────────
        if (world.rank() == 0) print_header2("Hylleraas functional");

        auto load_raw_u = [&](const std::string& fname)
                -> std::pair<bool, real_function_6d> {
            int found = 0;
            if (world.rank() == 0)
                found = std::ifstream(fname + ".00000").good() ? 1 : 0;
            world.gop.broadcast(found, 0);
            if (!found) return {false, real_factory_6d(world)};
            if (world.rank() == 0) printf("  loading u from %s\n", fname.c_str());
            real_function_6d u = real_factory_6d(world);
            archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(world, fname.c_str(), 1);
            ar & u;
            return {true, u};
        };

        auto load_ccpair_u = [&](const std::string& fname, size_t i, size_t j)
                -> std::pair<bool, real_function_6d> {
            int found = 0;
            if (world.rank() == 0)
                found = std::ifstream(fname + ".00000").good() ? 1 : 0;
            world.gop.broadcast(found, 0);
            if (!found) return {false, real_factory_6d(world)};
            if (world.rank() == 0) printf("  loading CCPair from %s\n", fname.c_str());
            CCPair cc_pair(i, j, GROUND_STATE, CT_MP2);
            archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(world, fname.c_str(), 1);
            ar & cc_pair;
            if (!cc_pair.functions.empty()) {
                try {
                    real_function_6d u = cc_pair.function();
                    if (u.is_initialized()) return {true, u};
                } catch (...) {}
            }
            if (cc_pair.constant_part.is_initialized())
                return {true, cc_pair.constant_part};
            return {false, real_factory_6d(world)};
        };

        for (size_t i = freeze; i < nocc; ++i) {
            for (size_t j = i; j < nocc; ++j) {
                std::string ij = std::to_string(i) + std::to_string(j);
                real_function_6d u_ij = real_factory_6d(world);
                std::string loaded_from = "zero";

                auto [ok1, u1] = load_raw_u("iter_final.MP2_pair_u_" + ij);
                if (ok1) { u_ij = u1; loaded_from = "iter_final.MP2_pair_u_" + ij; }
                else {
                    auto [ok2, u2] = load_ccpair_u("MP2_pair_u_" + ij + "_const", i, j);
                    if (ok2) { u_ij = u2; loaded_from = "MP2_pair_u_" + ij + "_const"; }
                }
                if (loaded_from == "zero" && world.rank() == 0)
                    printf("  no archive found for pair (%zu,%zu), using u = 0\n", i, j);

                double e_hyl = BMatrix::compute_hylleraas_pair(world, info, i, j, u_ij);
                if (world.rank() == 0)
                    printf("  e_hyl(%zu,%zu) = %+.10f  (u from %s)\n",
                           i, j, e_hyl, loaded_from.c_str());
            }
        }

        // ── Build common operators ─────────────────────────────────────────
        auto fock_ptr = nemo->make_fock_operator();
        auto f12_sep  = CCConvolutionOperatorPtr<double,3>(world, OT_F12, ccparam)->get_op();
        auto& molecule = nemo->get_calc()->molecule;

        // ── (a) Commutator: exact 6D [K,f12] ─────────────────────────────
        if (world.rank() == 0)
            print_header2("(a) Commutator: B^Ue + B^g + B^KffK_6d + B^X  (exact 6D)");

        double t0_comm = wall_time();
        auto B_comm = BMatrix::compute_6d(world, info);
        double t_comm = wall_time() - t0_comm;

        if (world.rank() == 0) {
            for (size_t ii = freeze; ii < nocc; ++ii)
                for (size_t jj = ii; jj < nocc; ++jj)
                    printf("  B_6d(%zu,%zu,%zu,%zu) = %+.10f\n",
                           ii, jj, ii, jj,
                           B_comm(ii-freeze, jj-freeze, ii-freeze, jj-freeze));
            printf("  wall time: %.1f s\n", t_comm);
        }

        // ── (b) LRF-B: orbital-weighted LRF ──────────────────────────────
        if (world.rank() == 0)
            print_header2("(b) LRF-B: orbital-weighted f12*phi_i*phi_j, random grid");

        LowRankFunctionParameters lrf_params;
        lrf_params.set_derived_value<double>("radius", 5.0);
        lrf_params.set_derived_value<double>("volume_element", 0.1);
        lrf_params.set_derived_value<std::string>("gridtype", std::string("random"));
        lrf_params.set_derived_value<bool>("canonicalize", true);

        double t0_lrf = wall_time();
        std::vector<std::pair<size_t,size_t>> pairs;
        std::vector<double> B_lrf_vals;

        for (size_t ii = freeze; ii < nocc; ++ii) {
            for (size_t jj = ii; jj < nocc; ++jj) {
                LRFunctorF12<double,6> lrf_orb(f12_sep, info.mo_ket[ii], info.mo_ket[jj]);
                auto lrf_ij = LowRankFunctionFactory<double,6>(lrf_params, molecule)
                                  .project(lrf_orb, 1.e-3);
                if (world.rank() == 0)
                    printf("  LRF pair (%zu,%zu): rank = %d  l2_rel = %.4e\n",
                           ii, jj, (int)lrf_ij.rank()[0], lrf_ij.l2error(lrf_orb));
                double B_val = BMatrix::compute_via_lrf_pair(world, info, ii, jj, lrf_ij, fock_ptr);
                pairs.push_back({ii, jj});
                B_lrf_vals.push_back(B_val);
            }
        }
        double t_lrf = wall_time() - t0_lrf;

        // ── Summary table ─────────────────────────────────────────────────
        if (world.rank() == 0) {
            print_header2("Summary");
            printf("  thresh_6D = %.1e\n\n", ccparam.thresh_6D());
            printf("  %-26s  %14s  %14s  %8s\n",
                   "method", "B(i,j,i,j)", "vs commutator", "time(s)");
            printf("  %s\n", std::string(70, '-').c_str());
            for (size_t idx = 0; idx < pairs.size(); ++idx) {
                auto [ii, jj] = pairs[idx];
                double Bc = B_comm(ii-freeze, jj-freeze, ii-freeze, jj-freeze);
                printf("  commutator (%zu,%zu)        %+14.8f  %14s  %8.1f\n",
                       ii, jj, Bc, "---", t_comm);
                printf("  LRF-B      (%zu,%zu)        %+14.8f  %+14.8f  %8.1f\n",
                       ii, jj, B_lrf_vals[idx], B_lrf_vals[idx] - Bc, t_lrf);
            }
            if (status == 0) print("\n  PASSED");
        }

    } catch (std::exception& e) {
        if (world.rank() == 0) print("exception:", e.what());
        status = 1;
    }
#else
    if (world.rank() == 0)
        print("test_bmatrix requires compilation with ENABLE_GENTENSOR=1");
    status = 77; // skip
#endif
    finalize();
    return status;
}
