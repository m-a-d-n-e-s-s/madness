/*
 * test_bsh_apply_strategy.cc
 *
 * Tier-1 unit test for the BSH-apply backend chooser (BSHApplyStrategy.h). The chooser is
 * a pure function, so this encodes the measured regime grid (water10/C40/valinomycin/w213)
 * and the batch ramp as assertions -- the CI regression guard for the dispatch logic.
 * See docs/bsh_apply_macrotask.md "IMPLEMENTATION SPEC" (test plan, Tier 1).
 */

#include <madness/chem/BSHApplyStrategy.h>
#include <cstdio>
#include <string>

using namespace madness;

namespace {

int failures = 0;

const char* name(BSHBackend b) { return b == BSHBackend::Tile ? "tile" : "macrotask"; }

// per-node RAM / ranks-per-node, the realistic per-rank budget on the 1 TB test nodes
double budget_per_rank_gb(double node_ram_gb, long ranks_per_node) {
    return node_ram_gb / double(std::max<long>(1, ranks_per_node));
}

void expect_backend(const std::string& label, const BSHApplyContext& c, BSHBackend want) {
    BSHApplyPlan p = choose_bsh_apply_plan(c);
    const bool ok = (p.backend == want);
    std::printf("  [%s] %-9s expected=%-9s got=%-9s : %s\n",
                ok ? "PASS" : "FAIL", label.c_str(), name(want), name(p.backend),
                p.reason.c_str());
    if (!ok) ++failures;
}

void expect_batch(const std::string& label, const BSHApplyContext& c, double thresh, long want) {
    long got = choose_bsh_batch(c, thresh);
    const bool ok = (got == want);
    std::printf("  [%s] %-20s batch(thresh=%.0e) expected=%ld got=%ld\n",
                ok ? "PASS" : "FAIL", label.c_str(), thresh, want, got);
    if (!ok) ++failures;
}

void expect_true(const std::string& label, bool cond) {
    std::printf("  [%s] %s\n", cond ? "PASS" : "FAIL", label.c_str());
    if (!cond) ++failures;
}

}  // namespace

int main(int /*argc*/, char** /*argv*/) {
    std::printf("=== BSH-apply chooser: measured regime grid ===\n");
    const double b125 = budget_per_rank_gb(1000.0, 8);   // 1 TB node, 8 ranks/node
    const double b500 = budget_per_rank_gb(1000.0, 2);   // 1 TB node, 2 ranks/node (w213 layout)

    // single node -> tile (it fits and is faster: no gather, no inter-node comm)
    expect_backend("water10", {/*nmo*/50,  /*k*/10, /*nproc*/8,  /*n_nodes*/1, b125, 3}, BSHBackend::Tile);
    expect_backend("C40",     {/*nmo*/161, /*k*/10, /*nproc*/8,  /*n_nodes*/1, b125, 3}, BSHBackend::Tile);
    // multinode -> macrotask (wins/ties speed + bounds memory; w213 tile OOMs)
    expect_backend("valino",  {/*nmo*/300, /*k*/10, /*nproc*/64, /*n_nodes*/8, b125, 3}, BSHBackend::Macrotask);
    expect_backend("w213",    {/*nmo*/1065,/*k*/10, /*nproc*/60, /*n_nodes*/30,b500, 3}, BSHBackend::Macrotask);

    std::printf("=== edge cases ===\n");
    // single node, no budget given -> gate off -> tile
    expect_backend("1node/no-budget", {161, 10, 8, 1, /*budget*/0.0, 3}, BSHBackend::Tile);
    // single node but would not fit (tight RAM, large k) -> macrotask (OOM safety net)
    expect_backend("1node/OOM",       {200, 14, 8, 1, /*budget*/10.0, 3}, BSHBackend::Macrotask);
    // multinode always macrotask even for a tiny system
    expect_backend("2node/tiny",      {20,  6,  2, 2, b125, 3}, BSHBackend::Macrotask);

    std::printf("=== topology+protocol-aware batch ===\n");
    const BSHApplyContext c40b    {161, 10, 8,  1,  b125, 3};                       // single node
    const BSHApplyContext valinob {300, 10, 64, 8,  b125, 3};                       // 8 ranks/node (intra)
    const BSHApplyContext w213b    {1065,10, 60, 30, budget_per_rank_gb(384,2), 3}; // 2/node, real 384GB node
    const BSHApplyContext w213_1tb {1065,10, 60, 30, b500, 3};                      // 2/node, 1 TB node
    const BSHApplyContext w213_nob {1065,10, 60, 30, 0.0,  3};                      // 2/node, no budget
    // loose/medium: working set tiny -> batch up everywhere
    expect_batch("C40 loose",        c40b,     1.0e-4, 4);
    expect_batch("w213 medium",      w213b,    1.0e-6, 4);
    // tight, single-node or intra-node-multinode -> b=1 (memory / compute-monolith)
    expect_batch("C40 tight",        c40b,     1.0e-8, 1);
    expect_batch("valino tight",     valinob,  1.0e-8, 1);   // rpn=8 -> not inter-node-sparse
    // tight + inter-node-sparse (>=2 nodes, <=2 ranks/node) + budget -> widen, memory-bounded
    expect_batch("w213 tight",       w213b,    1.0e-8, 5);   // 0.10*192/3.44 -> 5
    expect_batch("w213 tight 1TB",   w213_1tb, 1.0e-8, 8);   // 0.10*500/3.44 -> 14, capped at B_HI=8
    // tight + inter-node-sparse but no budget -> cannot bound the working set -> b=1
    expect_batch("w213 tight no-bdg",w213_nob, 1.0e-8, 1);

    std::printf("=== estimator sanity (C40-like ~34 GB/rank) ===\n");
    const double w_c40 = estimate_tile_peak_per_rank_gb({161, 10, 8, 1, 0.0, 3});
    std::printf("  estimate(C40) = %.1f GB/rank\n", w_c40);
    expect_true("estimate(C40) in [25,45] GB/rank", w_c40 > 25.0 && w_c40 < 45.0);

    std::printf("%s: %d failure(s)\n", failures == 0 ? "PASSED" : "FAILED", failures);
    return failures == 0 ? 0 : 1;
}
