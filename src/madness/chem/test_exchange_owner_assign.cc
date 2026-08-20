/// \file test_exchange_owner_assign.cc
/// \brief tests the pure batch-splitting and owner-assignment helpers of the exchange operator

#include <madness/chem/exchangeoperator.h>
#include <madness/world/test_utilities.h>

#include <cmath>
#include <map>
#include <set>
#include <vector>

using namespace madness;

namespace {

/// the property that makes owner-pinning fetch-free: a task's owner stores one of its batches
bool owner_is_eligible(const long i, const long j, const long owner, const long n) {
    if (i == j) return owner == (i % n);
    return owner == (i % n) or owner == (j % n);
}

/// every triangular task assigned exactly once, and always to an eligible worker
bool assignment_is_valid(const std::map<std::pair<long,long>,long>& owner,
                         const long n, const long M) {
    if (owner.size() != std::size_t(M) * std::size_t(M + 1) / 2) return false;
    for (long i = 0; i < M; ++i) {
        for (long j = 0; j <= i; ++j) {
            auto it = owner.find({i, j});
            if (it == owner.end()) return false;
            if (it->second < 0 or it->second >= n) return false;
            if (not owner_is_eligible(i, j, it->second, n)) return false;
        }
    }
    return true;
}

std::vector<long> task_counts(const std::map<std::pair<long,long>,long>& owner, const long n) {
    std::vector<long> c(n, 0);
    for (const auto& [task, r] : owner) ++c[r];
    return c;
}

double max_load(const std::map<std::pair<long,long>,long>& owner, const long n,
                const std::vector<double>& cost) {
    std::vector<double> load(n, 0.0);
    for (const auto& [task, r] : owner) load[r] += cost[exchange_sym_tri(task.first, task.second)];
    return *std::max_element(load.begin(), load.end());
}

/// a cost model with the inhomogeneity screening actually produces: a few expensive
/// batches (compact core orbitals) among cheap ones, so a pair's cost is set by its operands
std::vector<double> hot_batch_cost(const long M, const std::set<long>& hot) {
    std::vector<double> cost(std::size_t(M) * std::size_t(M + 1) / 2, 0.0);
    auto weight = [&hot](long k) { return hot.count(k) ? 10.0 : 1.0; };
    for (long i = 0; i < M; ++i)
        for (long j = 0; j <= i; ++j)
            cost[exchange_sym_tri(i, j)] = weight(i) * weight(j);
    return cost;
}

}   // namespace

int test_batch_split() {
    test_output t1("testing exchange_sym_owner_split");

    t1.checkpoint(exchange_sym_owner_nbatch(0, 8, 1) == 0, "no batches for an empty vector");
    t1.checkpoint(exchange_sym_owner_split(0, 8, 1).empty(), "empty split for an empty vector");
    t1.checkpoint(exchange_sym_owner_nbatch(40, 8, 1) == 8, "level 1 gives one batch per rank");
    t1.checkpoint(exchange_sym_owner_nbatch(40, 8, 3) == 24, "level 3 gives three batches per rank");
    t1.checkpoint(exchange_sym_owner_nbatch(5, 8, 2) == 5, "clamped: never more batches than orbitals");
    t1.checkpoint(exchange_sym_owner_nbatch(40, 0, 0) == 1, "degenerate rank/level counts clamp to 1");

    bool ok = true;
    for (long n = 1; n <= 64; ++n) {
        for (long nsw : {1L, 2L, 3L, 8L, 40L}) {
            for (long level : {1L, 2L, 5L}) {
                const auto split = exchange_sym_owner_split(n, nsw, level);
                const long M = exchange_sym_owner_nbatch(n, nsw, level);
                if (long(split.size()) != M) ok = false;
                long covered = 0, mn = n, mx = 0, expect_begin = 0;
                for (const auto& b : split) {
                    if (b.begin != expect_begin) ok = false;      // contiguous, no gap or overlap
                    if (b.size() <= 0) ok = false;                // no empty batch
                    mn = std::min(mn, b.size());
                    mx = std::max(mx, b.size());
                    covered += b.size();
                    expect_begin = b.end;
                }
                if (covered != n) ok = false;                     // covers [0,n)
                if (mx - mn > 1) ok = false;                      // sizes differ by at most one
            }
        }
    }
    t1.checkpoint(ok, "splits are contiguous, cover [0,n), and are balanced to within one");
    return t1.end();
}

int test_row_owner_split() {
    test_output t1("testing exchange_row_owner_split");

    t1.checkpoint(exchange_row_owner_split(0, 8).empty(), "empty split for an empty vector");

    bool same = true, one_per_rank = true;
    for (long n = 1; n <= 64; ++n) {
        for (long nsw : {1L, 2L, 3L, 8L, 40L}) {
            const auto row = exchange_row_owner_split(n, nsw);
            if (row != exchange_sym_owner_split(n, nsw, 1)) same = false;
            if (long(row.size()) != std::min<long>(std::max<long>(1, nsw), n)) one_per_rank = false;
        }
    }
    // the contiguity, coverage and balance of these boundaries are covered by test_batch_split
    t1.checkpoint(same, "the same boundaries as the granularity-1 symmetric split");
    t1.checkpoint(one_per_rank, "exactly min(nsubworld, n) batches, which the column-to-rank assignment needs");
    return t1.end();
}

int test_row_owner_grid() {
    test_output t1("testing exchange_row_owner_grid");

    // deliberately ncolumn != nrow: vf and the bra/ket pair are independent lengths, and every
    // consumer so far has had them equal, so an implementation that conflated them would pass
    // every existing test
    const long ncol = 3, nrow = 8, nsw = 4;
    const auto grid = exchange_row_owner_grid(ncol, nrow, nsw);
    const auto columns = exchange_row_owner_split(ncol, nsw);
    const auto rows = exchange_row_owner_split(nrow, nsw);

    t1.checkpoint(grid.size() == columns.size() * rows.size(), "the grid is the full rectangle");

    // assert which split feeds which slot, not just the shape: a swap gives a grid of the same
    // size whenever the two lengths happen to agree, which is the case that has always been run
    bool sides_ok = true, pairs_unique = true;
    std::set<std::pair<long,long>> seen;
    for (const auto& [c, r] : grid) {
        if (std::find(columns.begin(), columns.end(), c) == columns.end()) sides_ok = false;
        if (std::find(rows.begin(), rows.end(), r) == rows.end()) sides_ok = false;
        if (not seen.insert({c.begin, r.begin}).second) pairs_unique = false;
    }
    t1.checkpoint(sides_ok, "the column slot comes from the vf split and the row slot from the bra split");
    t1.checkpoint(pairs_unique and seen.size() == grid.size(), "every column/row pair appears exactly once");

    // one batch per rank per dimension, so the shorter dimension is what limits the columns
    t1.checkpoint(long(columns.size()) == std::min(ncol, nsw) and long(rows.size()) == std::min(nrow, nsw),
                  "each dimension is split to one batch per rank, clamped by its own length");
    t1.checkpoint(exchange_row_owner_grid(0, nrow, nsw).empty() and
                  exchange_row_owner_grid(ncol, 0, nsw).empty(),
                  "an empty dimension gives an empty grid");
    return t1.end();
}

int test_row_owner_assignment() {
    test_output t1("testing exchange_row_owner_assign");

    // ncolumn != nrow throughout: the two dimensions are independent lengths
    bool one_rank_per_column = true, complete = true, in_range = true, spread_ok = true;
    for (long ncol : {1L, 2L, 3L, 5L, 8L}) {
        for (long nrow : {1L, 4L, 7L}) {
            for (long nworker : {1L, 2L, 4L, 8L}) {
                const auto owner = exchange_row_owner_assign(ncol, nrow, nworker);
                if (long(owner.size()) != ncol * nrow) complete = false;
                std::map<long, std::set<long>> ranks_of_column;
                std::vector<long> per_rank(nworker, 0);
                for (const auto& [task, rank] : owner) {
                    if (rank < 0 or rank >= nworker) in_range = false;
                    ranks_of_column[task.first].insert(rank);
                    ++per_rank[rank];
                }
                for (const auto& [column, ranks] : ranks_of_column)
                    if (ranks.size() != 1) one_rank_per_column = false;
                // each worker gets whole columns, so the counts differ by at most one column
                const long mx = *std::max_element(per_rank.begin(), per_rank.end());
                const long mn = *std::min_element(per_rank.begin(), per_rank.end());
                if (mx - mn > nrow) spread_ok = false;
            }
        }
    }
    t1.checkpoint(complete, "every column/row task is assigned exactly once");
    t1.checkpoint(in_range, "every owner is a valid worker");
    t1.checkpoint(one_rank_per_column,
                  "all tasks of a column share one worker, so its held batch is never fetched");
    t1.checkpoint(spread_ok, "workers differ by at most one column's worth of tasks");

    // distinct columns land on distinct workers while there are workers to spare -- that is what
    // makes each rank hold a different vf batch
    const auto owner = exchange_row_owner_assign(4, 3, 4);
    std::set<long> distinct;
    for (long c = 0; c < 4; ++c) distinct.insert(owner.at({c, 0}));
    t1.checkpoint(distinct.size() == 4, "distinct columns go to distinct workers when enough exist");

    t1.checkpoint(exchange_row_owner_assign(4, 3, 4) == exchange_row_owner_assign(4, 3, 4),
                  "deterministic, so every rank agrees without communicating");
    t1.checkpoint(exchange_row_owner_assign(0, 3, 4).empty() and
                  exchange_row_owner_assign(4, 0, 4).empty() and
                  exchange_row_owner_assign(4, 3, 0).empty(),
                  "degenerate dimensions or worker count give an empty assignment");
    return t1.end();
}

int test_owner_assignment() {
    test_output t1("testing exchange owner assignment");

    bool rr_valid = true, ca_valid = true, spread_ok = true;
    for (long n : {1L, 2L, 3L, 8L, 13L}) {
        for (long level : {1L, 2L, 3L}) {
            const long M = n * level;
            const auto rr = exchange_sym_round_robin_assign(n, M);
            if (not assignment_is_valid(rr, n, M)) rr_valid = false;
            const auto counts = task_counts(rr, n);
            const long spread = *std::max_element(counts.begin(), counts.end())
                              - *std::min_element(counts.begin(), counts.end());
            if (spread > 2) spread_ok = false;

            const auto cost = hot_batch_cost(M, {0, 1});
            if (not assignment_is_valid(exchange_sym_cost_aware_assign(n, M, cost), n, M)) ca_valid = false;
        }
    }
    t1.checkpoint(rr_valid, "round-robin assigns every task once to an eligible owner");
    t1.checkpoint(ca_valid, "cost-aware assigns every task once to an eligible owner");
    // phase 2 aims for a spread of 1 but may only move a task to a worker it is eligible
    // for, so it stops early on some (n, M); 2 is the bound observed over the sweep
    t1.checkpoint(spread_ok, "round-robin balances the task count to within two");

    t1.checkpoint(exchange_sym_round_robin_assign(0, 8).empty() and
                  exchange_sym_round_robin_assign(8, 0).empty(),
                  "degenerate worker or batch counts give an empty assignment");

    // an empty cost vector must still produce a usable assignment -- the first two exchange
    // calls have no measured cost yet
    t1.checkpoint(assignment_is_valid(exchange_sym_cost_aware_assign(8, 16, {}), 8, 16),
                  "cost-aware degrades gracefully with no cost information");

    t1.checkpoint(exchange_sym_round_robin_assign(8, 16) == exchange_sym_round_robin_assign(8, 16) and
                  exchange_sym_cost_aware_assign(8, 16, hot_batch_cost(16, {0, 1}))
                          == exchange_sym_cost_aware_assign(8, 16, hot_batch_cost(16, {0, 1})),
                  "both are deterministic, so every rank agrees without communicating");

    // the reason the cost-aware strategy exists: with inhomogeneous tasks, balancing the
    // count leaves the worker holding the expensive batches as the straggler
    const long n = 8, M = 16;
    const auto cost = hot_batch_cost(M, {0, 8});      // both hot batches round-robin onto worker 0
    const double rr_max = max_load(exchange_sym_round_robin_assign(n, M), n, cost);
    const double ca_max = max_load(exchange_sym_cost_aware_assign(n, M, cost), n, cost);
    print("max per-worker cost: round-robin", rr_max, " cost-aware", ca_max);
    t1.checkpoint(ca_max < rr_max, "cost-aware lowers the peak load on an inhomogeneous cost model");

    return t1.end();
}

int main(int argc, char** argv) {
    madness::initialize(argc, argv);
    int result = 0;
    result += test_batch_split();
    result += test_row_owner_split();
    result += test_row_owner_grid();
    result += test_row_owner_assignment();
    result += test_owner_assignment();
    if (result == 0) print("\n --> all tests \033[32m passed \033[0m \n");
    else print("\n --> all tests \033[31m failed \033[0m \n");
    madness::finalize();
    return result;
}
