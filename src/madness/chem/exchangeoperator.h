#ifndef SRC_APPS_CHEM_EXCHANGEOPERATOR_H_
#define SRC_APPS_CHEM_EXCHANGEOPERATOR_H_

#include<madness.h>
#include<madness/world/cloud.h>
#include<madness/mra/macrotaskq.h>
#include<madness/chem/SCFOperators.h>

#include <algorithm>
#include <fstream>
#include <list>
#include <map>
#include <utility>
#include <vector>

namespace madness {

// forward declaration
class SCF;
class Nemo;

/// Number of batches M for the owner-pinned symmetric exchange algorithm.

/// M = granularity_level * nsubworld, so the level is directly the number of batches per
/// rank and every level is selectable regardless of nsubworld parity. Clamped to [1, n]:
/// there cannot be more batches than orbitals, and a caller that hits the clamp is in the
/// small-problem regime where batches stop being one-per-rank.
///
/// Batch k is owned by rank (k mod nsubworld), so away from the clamp every rank owns
/// exactly `granularity_level` batches.
inline long exchange_sym_owner_nbatch(const std::size_t n, const long nsubworld,
                                      const long granularity_level) {
    if (n == 0) return 0;
    const long nsw = std::max<long>(1, nsubworld);
    const long level = std::max<long>(1, granularity_level);
    return std::min<long>(level * nsw, long(n));
}

/// Split a vector of length n into M = exchange_sym_owner_nbatch(...) contiguous batches.

/// Sizes differ by at most 1; the first (n mod M) batches get the extra element, which
/// avoids a runt batch of 1. Single source of truth for the symmetric owner-pinned batch
/// boundaries — the task grid, the cloud batch storage and the owner assignment all derive
/// from it, so they align by construction. Depends only on (n, nsubworld, granularity_level).
inline std::vector<Batch_1D> exchange_sym_owner_split(const std::size_t n, const long nsubworld,
                                                      const long granularity_level) {
    std::vector<Batch_1D> out;
    if (n == 0) return out;
    const long nbatch = std::max<long>(1, exchange_sym_owner_nbatch(n, nsubworld, granularity_level));
    const long bs_floor = long(n) / nbatch;
    const long rem = long(n) - bs_floor * nbatch;
    long begin = 0;
    for (long b = 0; b < nbatch; ++b) {
        const long sz = bs_floor + (b < rem ? 1 : 0);
        out.emplace_back(begin, begin + sz);
        begin += sz;
    }
    return out;
}

/// Batch boundaries for the asymmetric row/column split: one batch per rank.

/// Deliberately the same boundaries as exchange_sym_owner_split at granularity 1, so batches stored
/// for one dimension align with the partition. Named apart because "sym" reads wrong in the
/// asymmetric path, and because granularity is not a knob here -- the column-to-rank assignment
/// relies on there being exactly one batch per rank.
inline std::vector<Batch_1D> exchange_row_owner_split(const std::size_t n, const long nsubworld) {
    return exchange_sym_owner_split(n, nsubworld, 1);
}

/// The asymmetric task grid: every (column, row) pair over two independent splits.

/// The column dimension is the operand that stays put -- vf, whose batch the running rank holds --
/// and the row dimension is the one that rotates, carrying bra and ket, which share this range
/// because the operator sums over pairs. There is no symmetry to exploit, so the grid is the full
/// rectangle rather than a triangle, and the two dimensions are split independently: bra and ket
/// have the same length as each other but need not match vf's.
inline std::vector<std::pair<Batch_1D, Batch_1D>>
exchange_row_owner_grid(const std::size_t ncolumn, const std::size_t nrow, const long nsubworld) {
    const std::vector<Batch_1D> columns = exchange_row_owner_split(ncolumn, nsubworld);
    const std::vector<Batch_1D> rows = exchange_row_owner_split(nrow, nsubworld);
    std::vector<std::pair<Batch_1D, Batch_1D>> out;
    out.reserve(columns.size() * rows.size());
    for (const auto& c : columns)
        for (const auto& r : rows) out.emplace_back(c, r);
    return out;
}

/// Owner of every task in the asymmetric grid: all tasks of a column go to one worker.

/// The column operand is the one that stays put, so putting a whole column on one worker means the
/// batch it holds is never fetched, and that worker computes the *complete* result for those
/// columns -- there is no cross-rank reduction of results beyond the final gather. Handing columns
/// out in order is balanced by construction, the grid being a full rectangle, so no cost model is
/// needed here. Keyed by batch index, like exchange_sym_round_robin_assign.
inline std::map<std::pair<long,long>,long>
exchange_row_owner_assign(const long ncolumn, const long nrow, const long nworker) {
    std::map<std::pair<long,long>,long> owner;
    if (ncolumn <= 0 or nrow <= 0 or nworker <= 0) return owner;
    for (long c = 0; c < ncolumn; ++c)
        for (long r = 0; r < nrow; ++r) owner[{c, r}] = c % nworker;
    return owner;
}

/// Triangular index of a batch pair, collapsing (i,j) and (j,i): a*(a+1)/2 + b.

/// The indexing convention of the per-task cost vector, shared by whoever records a
/// measured cost and by exchange_sym_cost_aware_assign, which consumes it.
inline long exchange_sym_tri(const long i, const long j) {
    const long a = std::max(i, j), b = std::min(i, j);
    return a * (a + 1) / 2 + b;
}

/// Round-robin owner assignment for the symmetric algorithm's triangular task matrix.

/// Pure function of (n, M): n workers (= nsubworld), M batches, batch k owned by (k mod n).
/// Task (i,j) with 0 <= j <= i < M is eligible for worker t iff
///   i==j: (i mod n)==t                    -- a diagonal task has a single owner
///   i!=j: (i mod n)==t or (j mod n)==t    -- the worker owns at least one operand
/// Eligibility is what makes the assignment fetch-free: a task always lands on a worker
/// that already stores one of its two batches. Phase 1 round-robins eligible tasks across
/// workers; phase 2 then rebalances off-diagonal tasks (diagonals never move) toward a
/// task-count spread of 1.
///
/// Phase 2 is best effort, not a guarantee: it can only move a task to a worker the task
/// is eligible for, so it stops early when the hottest worker holds nothing the coldest
/// ones may take. The residual spread is 2 rather than 1 for some (n, M) — measured over
/// n <= 128 and 1 to 3 batches per rank, never worse than 2.
///
/// \return map (i,j) -> owner. Deterministic, so every rank computes the same assignment
///         without communicating.
inline std::map<std::pair<long,long>,long>
exchange_sym_round_robin_assign(const long n, const long M) {
    std::map<std::pair<long,long>,long> owner;
    if (n <= 0 or M <= 0) return owner;
    auto eligible = [n](long i, long j, long t) -> bool {
        if (i == j) return (i % n) == t;
        return ((i % n) == t) or ((j % n) == t);
    };
    // all tasks, ascending by (i,j)
    const std::size_t ntask = std::size_t(M) * std::size_t(M + 1) / 2;
    std::vector<std::pair<long,long>> tasks;
    tasks.reserve(ntask);
    for (long i = 0; i < M; ++i)
        for (long j = 0; j <= i; ++j)
            tasks.emplace_back(i, j);

    std::vector<std::vector<std::pair<long,long>>> T(n);

    // Phase 1 — round-robin placement. Each worker takes the first task, in ascending
    // (i,j), that is still free and eligible for it. `cursor[t]` remembers how far worker
    // t has scanned: everything before it is either taken or permanently ineligible for t,
    // and neither condition is ever undone, so the cursor only moves forward. That makes
    // the phase linear in the task count per worker rather than quadratic overall.
    std::vector<char> taken(ntask, 0);
    std::vector<std::size_t> cursor(n, 0);
    std::size_t placed = 0;
    long t = 0, misses = 0;
    while (placed < ntask and misses < n) {
        std::size_t& c = cursor[t];
        while (c < ntask and (taken[c] or not eligible(tasks[c].first, tasks[c].second, t))) ++c;
        if (c < ntask) {
            taken[c] = 1;
            T[t].push_back(tasks[c]);
            ++placed;
            ++c;
            misses = 0;
        } else {
            ++misses;
        }
        t = (t + 1) % n;
    }
    // any leftovers; every task is eligible for someone, so this should not trigger
    for (std::size_t idx = 0; idx < ntask; ++idx)
        if (not taken[idx]) T[tasks[idx].first % n].push_back(tasks[idx]);

    // Phase 2 — rebalance until the spread is <= 1 (the iteration cap is a backstop)
    for (long iter = 0; iter < 10000; ++iter) {
        long big = 0, small = 0;
        for (long p = 1; p < n; ++p) {
            if (long(T[p].size()) > long(T[big].size())) big = p;
            if (long(T[p].size()) < long(T[small].size())) small = p;
        }
        if (long(T[big].size()) - long(T[small].size()) <= 1) break;
        bool moved = false;
        std::vector<long> order(n);
        for (long p = 0; p < n; ++p) order[p] = p;
        std::sort(order.begin(), order.end(),
                  [&T](long a, long b){ return T[a].size() > T[b].size(); });
        for (long bi = 0; bi < n and not moved; ++bi) {
            for (long si = n - 1; si > bi and not moved; --si) {
                const long bt = order[bi], st = order[si];
                if (long(T[bt].size()) - long(T[st].size()) <= 1) continue;
                for (long idx = 0; idx < long(T[bt].size()); ++idx) {
                    const auto& tk = T[bt][idx];
                    if (tk.first == tk.second) continue;           // never move diagonals
                    if (not eligible(tk.first, tk.second, st)) continue;
                    T[st].push_back(tk);
                    T[bt].erase(T[bt].begin() + idx);
                    moved = true;
                    break;
                }
            }
        }
        if (not moved) break;   // accept the residual imbalance
    }

    for (long p = 0; p < n; ++p)
        for (const auto& tk : T[p])
            owner[tk] = p;
    return owner;
}

/// Cost-aware owner assignment for the symmetric algorithm's triangular task matrix.

/// Same ownership and eligibility as exchange_sym_round_robin_assign, so the fetch-free
/// invariant is preserved, but it balances per-worker total COST instead of task COUNT.
/// Screening makes the tasks strongly inhomogeneous for large molecules, and that is what
/// a count-based assignment cannot see. `cost` is indexed by exchange_sym_tri and holds a
/// relative per-task cost, in practice the previous call's measured wall time; entries
/// past its end count as zero, so an empty or short vector degrades to cost-blind.
/// Greedy largest-cost-first onto the less loaded eligible worker, then a bounded local
/// search that relieves the hottest worker.
///
/// \return map (i,j) -> owner. Deterministic given the same `cost`.
inline std::map<std::pair<long,long>,long>
exchange_sym_cost_aware_assign(const long n, const long M, const std::vector<double>& cost) {
    std::map<std::pair<long,long>,long> owner;
    if (n <= 0 or M <= 0) return owner;
    auto C = [&](long i, long j) -> double {
        const long t = exchange_sym_tri(i, j); return (t < long(cost.size())) ? cost[t] : 0.0;
    };
    std::vector<double> load(n, 0.0);
    // diagonals first: they have no choice of owner
    for (long i = 0; i < M; ++i) { const long r = i % n; load[r] += C(i,i); owner[{i,i}] = r; }
    // off-diagonals, largest cost first, ties broken by (i,j) to stay deterministic
    std::vector<std::pair<long,long>> off;
    off.reserve(std::size_t(M) * std::size_t(M > 0 ? M - 1 : 0) / 2);
    for (long i = 0; i < M; ++i) for (long j = 0; j < i; ++j) off.emplace_back(i, j);
    std::sort(off.begin(), off.end(), [&](const std::pair<long,long>& A, const std::pair<long,long>& B){
        const double ca = C(A.first,A.second), cb = C(B.first,B.second);
        return (ca != cb) ? (ca > cb) : (A < B);
    });
    for (const auto& [i,j] : off) {
        const long ri = i % n, rj = j % n;
        const long r = (ri == rj) ? ri : (load[ri] <= load[rj] ? ri : rj);
        owner[{i,j}] = r; load[r] += C(i,j);
    }
    // relieve the hottest worker by flipping its off-diagonals to their other eligible owner
    for (long pass = 0; pass < 50; ++pass) {
        long hot = 0; for (long p = 1; p < n; ++p) if (load[p] > load[hot]) hot = p;
        bool improved = false;
        for (const auto& [i,j] : off) {
            auto it = owner.find({i,j});
            if (it->second != hot) continue;
            const long ri = i % n, rj = j % n;
            if (ri == rj) continue;                          // no choice
            const long other = (hot == ri) ? rj : ri;
            const double c = C(i,j);
            if (load[other] + c < load[hot]) {               // strictly lowers the hot worker
                load[hot] -= c; load[other] += c; it->second = other; improved = true;
            }
        }
        if (not improved) break;
    }
    return owner;
}

/// One task's record for the exchange profiler.

/// Written only when MAD_EXCH_TASK_PROFILE is set. What it adds over the aggregate counters is
/// **attribution**: the counters say how many batches arrived from where, this says which task
/// waited and for how long, so a straggler can be identified rather than inferred.
///
/// It deliberately does not carry a per-stage breakdown of the compute (multiply / apply /
/// multiply). That would mean timing calls inside the numerical kernels, and the aggregate split
/// already exists in the operator's own timers.
struct ExchTaskProfile {
    long task_id = -1;
    long universe_rank = 0;                 ///< keys the output file: one per process
    unsigned long subworld_id = 0;
    int subworld_nrank = 0;
    double thresh = 0.0;                    ///< which protocol tier this task ran in
    long k = 0;
    bool diagonal = false;
    long row_begin = 0, row_end = 0, col_begin = 0, col_end = 0;
    double wall_start = 0.0, wall_end = 0.0;
    double wait_for_data_wall = 0.0;        ///< task entry until its operands are in hand
    double compute_wall = 0.0, compute_cpu = 0.0;
    /// wall inside each stage of the tile loop, accumulated over its rows. Honest without adding
    /// any fence: every stage below runs with fence=true, so each completes before the next is
    /// timed. What they do not cover -- building the per-row update vector, the compresses and the
    /// accumulating gaxpys -- shows up as the emitted `other` residual.
    double mul1_wall = 0.0, apply_wall = 0.0, mul2_wall = 0.0, truncate_wall = 0.0;
    int operand_source = -1;                ///< worst of its fetches: 0 resident, 1 ahead, 2 cold
    bool waited = false;                    ///< a cold fetch happened, so this task paid latency
    double peak_rss_gb = 0.0;

    void reset() { *this = ExchTaskProfile(); }
    /// keep the worst source, since that is the one that set the task's wait
    void observe_fetch_tier(const int tier) {
        if (tier > operand_source) operand_source = tier;
        if (tier == 2) waited = true;
    }
};

/// Is per-task exchange profiling on? Read once per process.
inline bool exch_task_profile_enabled() {
    static const bool on = (std::getenv("MAD_EXCH_TASK_PROFILE") != nullptr);
    return on;
}

/// Send a symmetric application down the general (bra, ket, vf) path? Read once per process.

/// A symmetric application computes the same numbers either way -- the general path just does the
/// whole rectangle where the symmetric one does a triangle and reuses each intermediate twice. So
/// with this set, moldft becomes a differential test of the general path on data whose answer is
/// already known, which is otherwise reachable only from nemo and molresponse. Debug only.
inline bool exch_force_general_path() {
    static const bool on = (std::getenv("MAD_EXCH_FORCE_GENERAL") != nullptr);
    return on;
}

/// Append one record to exch_taskprof.r<rank>.jsonl.

/// The stream stays open for the life of the process: one file per rank appended across every
/// application, which bounds the file count at the rank count however many iterations run, and
/// avoids a per-task open/close -- expensive metadata traffic on a parallel filesystem.
inline void exch_write_task_profile(const ExchTaskProfile& p) {
    static std::ofstream os;
    if (not os.is_open()) {
        os.open("exch_taskprof.r" + std::to_string(p.universe_rank) + ".jsonl", std::ios::app);
        if (not os.is_open()) return;
    }
    os << "{\"task\":" << p.task_id
       << ",\"rank\":" << p.universe_rank
       << ",\"subworld\":" << p.subworld_id
       << ",\"subworld_nrank\":" << p.subworld_nrank
       << ",\"thresh\":" << p.thresh
       << ",\"k\":" << p.k
       << ",\"diagonal\":" << (p.diagonal ? "true" : "false")
       << ",\"row\":[" << p.row_begin << "," << p.row_end << "]"
       << ",\"col\":[" << p.col_begin << "," << p.col_end << "]"
       << ",\"wall_start\":" << p.wall_start
       << ",\"wall\":" << (p.wall_end - p.wall_start)
       << ",\"wait_for_data\":" << p.wait_for_data_wall
       << ",\"compute_wall\":" << p.compute_wall
       << ",\"compute_cpu\":" << p.compute_cpu
       << ",\"compute_components_wall\":{"
       <<   "\"mul1\":" << p.mul1_wall << ",\"apply\":" << p.apply_wall
       <<   ",\"mul2\":" << p.mul2_wall << ",\"truncate\":" << p.truncate_wall
       <<   ",\"other\":" << (p.compute_wall - p.mul1_wall - p.apply_wall
                               - p.mul2_wall - p.truncate_wall)
       << "}"
       << ",\"operand_source\":" << p.operand_source
       << ",\"waited\":" << (p.waited ? "true" : "false")
       << ",\"peak_rss_gb\":" << p.peak_rss_gb
       << "}\n";
    os.flush();   // keep it readable mid-run; there is one write per task, not per node
}

/// Write the measured per-task cost matrix of one application, for offline inspection.

/// Behind the same flag as the rest of the profiler. One file per application per rank 0:
/// `Ccall<NNN>_k<K>.csv`, holding `i,j,cost` over the lower-triangular batch grid, with the call
/// index, wavelet order and batch count in a leading comment so a reader needs no other context.
/// This is the matrix the next application's placement is derived from, so dumping it is how a
/// straggler can be traced back to the cost that put it there.
inline void exch_write_cost_matrix(const long call_index, const long k, const long M,
                                   const std::vector<double>& cost) {
    if (M <= 0 or cost.empty()) return;
    char name[64];
    std::snprintf(name, sizeof(name), "Ccall%03ld_k%ld.csv", call_index, k);
    std::ofstream os(name);
    if (not os.is_open()) return;
    os << "# call=" << call_index << " k=" << k << " M=" << M << "\n";
    os << "i,j,cost\n";
    for (long i = 0; i < M; ++i)
        for (long j = 0; j <= i; ++j) {
            const long t = i * (i + 1) / 2 + j;
            if (t < long(cost.size())) os << i << "," << j << "," << cost[t] << "\n";
        }
}

/// which of the three exchange operand vectors a stored batch belongs to
enum ExchangeBatchDim { EXCHANGE_BATCH_VF = 0, EXCHANGE_BATCH_BRA = 1, EXCHANGE_BATCH_KET = 2 };

/// Per-invocation salt for the exchange batch record keys, from the ket identities.

/// Taken from the ket vector because that one is available in full on both sides — where
/// the batches are stored and where a task fetches them — so both derive the same record
/// keys without anyone communicating a manifest. Function implementation ids are handed out
/// in collective creation order, so the salt is identical on every rank and changes when
/// the operator is applied to freshly built functions.
///
/// \warning It keys on identity, not on content: two applications over the same function
///          objects produce the same salt, so a cache keyed by it must not outlive an
///          in-place mutation of those functions.
template<typename T, std::size_t NDIM>
inline long exchange_batch_salt(const std::vector<Function<T, NDIM>>& ket) {
    std::size_t k = 0x5a17ull;
    for (const auto& f : ket) hash_combine(k, hash_value(f.get_impl()->id()));
    return long(k);
}

/// Deterministic cloud record key for one stored batch: (salt, dimension, range).

/// Range-keyed, so the storing side and the fetching side agree on the key by construction.
inline long exchange_batch_record_key(const long salt, const int dim, const Batch_1D& r) {
    std::size_t k = std::size_t(salt);
    hash_combine(k, std::size_t(dim));
    hash_combine(k, std::size_t(r.begin));
    hash_combine(k, std::size_t(r.end));
    return long(k);
}

/// Bounded cache of fetched exchange batches, keyed by cloud record key.

/// A rank reuses the batches it owns across all of its tasks, so those are **pinned** and
/// never evicted; only the batches fetched from elsewhere are transient and bounded. An
/// earlier single bound over both let the churning transients evict the reusable owned
/// batches, which re-fetched them at every owned-batch switch. Pinning the owned ones costs
/// a fixed share of the data (orbitals per rank) and is independent of batch granularity.
///
/// Entries are held in a list, so a reference returned by find() or insert() stays valid
/// until that entry is evicted — promotion and insertion of other entries do not move it.
template<typename keyT, typename dataT>
class ExchangeBatchLRU {
public:
    /// how many non-owned entries may be resident; at least one is always allowed
    void set_transient_capacity(const std::size_t c) { transient_capacity_ = std::max<std::size_t>(1, c); }
    std::size_t transient_capacity() const { return transient_capacity_; }

    bool contains(const keyT& key) const {
        for (const auto& s : slots_) if (s.key == key) return true;
        return false;
    }

    /// \return the cached batch, promoted to most-recently-used, or nullptr if absent
    const dataT* find(const keyT& key) {
        for (auto it = slots_.begin(); it != slots_.end(); ++it) {
            if (it->key == key) {
                slots_.splice(slots_.begin(), slots_, it);
                return &slots_.front().data;
            }
        }
        return nullptr;
    }

    /// Insert as most-recently-used. `pinned` marks a batch this rank owns.
    const dataT& insert(const keyT& key, dataT&& data, const bool pinned) {
        slots_.push_front(Slot{key, std::move(data), pinned});
        while (n_transient() > transient_capacity_) {
            // drop the least-recently-used transient entry; pinned ones are skipped
            for (auto it = slots_.end(); it != slots_.begin(); ) {
                --it;
                if (not it->pinned) { slots_.erase(it); break; }
            }
        }
        return slots_.front().data;
    }

    /// drop every entry; the capacity setting survives
    void clear() { slots_.clear(); }

    std::size_t size() const { return slots_.size(); }

    std::size_t n_transient() const {
        std::size_t c = 0;
        for (const auto& s : slots_) if (not s.pinned) ++c;
        return c;
    }

private:
    struct Slot { keyT key; dataT data; bool pinned; };
    std::list<Slot> slots_;                  // front = most recently used
    std::size_t transient_capacity_ = 2;
};


/// One coefficient node in transit during the exchange finalize.

/// Carries the destination function index, the tree-node key and the whole node, which is
/// the same content the per-node active message would have sent -- just batched.
template <typename T, std::size_t NDIM>
struct FinalizeNodeRec {
    std::size_t f = 0;                  ///< index into the destination function vector
    Key<NDIM> key;                      ///< tree-node key
    FunctionNode<T, NDIM> node;         ///< the source node
    template <typename Archive>
    void serialize(Archive& ar) { ar & f & key & node; }
};

/// Receiving end of the exchange finalize, living in the world the transfer rides on.

/// Sources push bulk chunks here with task(), so deserializing and accumulating happens on a
/// worker and never on the communication thread. The arithmetic is the same as
/// FunctionNode::gaxpy_inplace; the accessor write lock is what serializes two senders that
/// reach the same key at once.
template <typename T, std::size_t NDIM>
class FinalizeReducer : public WorldObject<FinalizeReducer<T, NDIM>> {
public:
    typedef FunctionImpl<T, NDIM> implT;
    typedef FinalizeNodeRec<T, NDIM> recT;

    explicit FinalizeReducer(World& world)
        : WorldObject<FinalizeReducer<T, NDIM>>(world) {
        this->process_pending();
    }

    /// point at the destinations for the next drain; local, called by every rank
    void set_targets(std::vector<implT*> dests, T beta) {
        dests_ = std::move(dests);
        beta_ = beta;
    }

    void accumulate_chunk(const std::vector<recT>& recs) {
        for (const auto& r : recs) {
            // always on: this guards a raw pointer dereference, so it must not be a check that
            // release builds compile out. ASSERTION_TYPE=disable is a supported configuration
            MADNESS_CHECK_THROW(r.f < dests_.size() and dests_[r.f] != nullptr,
                                "exchange finalize: a chunk names a destination this rank has not "
                                "been given, so its reducer targets were never set");
            typename implT::dcT::accessor acc;
            dests_[r.f]->get_coeffs().insert(acc, r.key);   // get-or-create, locks the key
            acc->second.template gaxpy_inplace<T, T>(T(1.0), r.node, beta_);
        }
    }

private:
    std::vector<implT*> dests_;
    T beta_ = T(1.0);
};

/// Accumulate `src_vec` into `dest_vec` in bulk, one message per destination rank per chunk.

/// Replaces one active message per source tree node. The two vectors may live in different
/// worlds -- subworld or node into universe, or subworld into node -- so routing uses the
/// **destination** process map while the transfer rides `transport_world`, which must be the
/// destination's world. Completion is the caller's own fence.
///
/// \warning Collective on `transport_world`: every rank must call it once per finalize. The
///          fence below is a readiness barrier, so that no chunk can arrive at a rank that
///          has not yet repointed its reducer. A rank whose `src_vec` is empty still has to
///          reach it -- **do not add an early return for having nothing to send.**
template <typename T, std::size_t NDIM>
void coalesced_gaxpy(World& transport_world,
                     FinalizeReducer<T, NDIM>& reducer,
                     std::vector<Function<T, NDIM>>& dest_vec,
                     std::vector<Function<T, NDIM>>& src_vec,
                     const T beta,
                     const std::size_t chunk_entries) {
    typedef FunctionImpl<T, NDIM> implT;
    typedef FinalizeNodeRec<T, NDIM> recT;

    std::vector<implT*> dests(dest_vec.size(), nullptr);
    for (std::size_t f = 0; f < dest_vec.size(); ++f)
        dests[f] = dest_vec[f].get_impl().get();
    reducer.set_targets(dests, beta);
    transport_world.gop.fence();          // the readiness barrier -- see the warning above

    std::map<ProcessID, std::vector<recT>> buckets;
    const std::size_t nf = std::min(src_vec.size(), dest_vec.size());
    for (std::size_t f = 0; f < nf; ++f) {
        auto simpl = src_vec[f].get_impl();
        auto dimpl = dest_vec[f].get_impl();
        if (not simpl or not dimpl) continue;
        const implT& dref = *dimpl;
        for (auto it = simpl->get_coeffs().begin(); it != simpl->get_coeffs().end(); ++it) {
            const ProcessID owner = dref.get_coeffs().owner(it->first);
            std::vector<recT>& b = buckets[owner];
            b.push_back(recT{f, it->first, it->second});
            if (b.size() >= chunk_entries) {
                reducer.task(owner, &FinalizeReducer<T, NDIM>::accumulate_chunk, b);
                b.clear();
            }
        }
    }
    for (auto& kv : buckets)
        if (not kv.second.empty())
            reducer.task(kv.first, &FinalizeReducer<T, NDIM>::accumulate_chunk, kv.second);
}


template<typename T, std::size_t NDIM>
class Exchange<T,NDIM>::ExchangeImpl {
    typedef Function<T, NDIM> functionT;
    typedef std::vector<functionT> vecfuncT;

    static inline std::atomic<long> apply_timer;
    static inline std::atomic<long> mul2_timer;
    static inline std::atomic<long> mul1_timer; ///< timing
    static inline double elapsed_time;

    static void reset_timer() {
        mul1_timer = 0l;
        mul2_timer = 0l;
        apply_timer = 0l;
        elapsed_time = 0.0;
        MacroTaskExchangeSimple::reset_batch_cache_counters();
    }

public:
    nlohmann::json gather_timings(World& world) const {
        double t1 = double(mul1_timer) * 0.001;
        double t2 = double(apply_timer) * 0.001;
        double t3 = double(mul2_timer) * 0.001;
        // the operand batches a task found resident vs had to fetch from their owner. Summed
        // here rather than read off the local counters, so this reports the whole run like
        // every other number in this object
        double resident = double(MacroTaskExchangeSimple::batch_cache_hits());
        double fetched = double(MacroTaskExchangeSimple::batch_cache_misses());
        double ahead = double(MacroTaskExchangeSimple::batch_prefetch_hits());
        world.gop.sum(t1);
        world.gop.sum(t2);
        world.gop.sum(t3);
        world.gop.sum(resident);
        world.gop.sum(fetched);
        world.gop.sum(ahead);
        nlohmann::json j;
        j["multiply1"] = t1;
        j["apply"] = t2;
        j["multiply2"] = t3;
        j["total"] = elapsed_time;
        j["batch_cache_hits"] = long(resident);
        j["batch_cache_misses"] = long(fetched);
        j["batch_prefetch_hits"] = long(ahead);
        return j;
    }

    void print_timer(World& world) const {
        auto timings= gather_timings(world);
        if (world.rank() == 0) {
            printf(" cpu time spent in multiply1   %8.2fs\n", timings["multiply1"].template get<double>());
            printf(" cpu time spent in apply       %8.2fs\n", timings["apply"].template get<double>());
            printf(" cpu time spent in multiply2   %8.2fs\n", timings["multiply2"].template get<double>());
            printf(" total wall time               %8.2fs\n", timings["total"].template get<double>());
            // only the owner-pinned path fetches operand batches, so this stays quiet for
            // every other algorithm; zero fetches would mean that path never ran
            const long resident = timings["batch_cache_hits"].template get<long>();
            const long fetched = timings["batch_cache_misses"].template get<long>();
            const long ahead = timings["batch_prefetch_hits"].template get<long>();
            if (resident + fetched + ahead > 0)
                printf(" operand batches resident/ahead/fetched %6ld /%6ld /%6ld\n",
                       resident, ahead, fetched);
        }
    }


    typedef Exchange<T,NDIM>::ExchangeAlgorithm Algorithm;
    Algorithm algorithm_ = multiworld_efficient_row;
    MacroTaskInfo macro_task_info = MacroTaskInfo::preset("default");

    /// default ctor
    ExchangeImpl(World& world, const double lo, const double thresh) : world(world), lo(lo), thresh(thresh) {}

    /// ctor with a conventional calculation
    ExchangeImpl(World& world, const SCF *calc, const int ispin) ;

    /// ctor with a nemo calculation
    ExchangeImpl(World& world, const Nemo *nemo, const int ispin);

    /// set the bra and ket orbital spaces, and the occupation

    /// @param[in]	bra		bra space, must be provided as complex conjugate
    /// @param[in]	ket		ket space
    void set_bra_and_ket(const vecfuncT& bra, const vecfuncT& ket) {
        mo_bra = copy(world, bra);
        mo_ket = copy(world, ket);
    }

    std::string info() const {return "K";}

    static auto set_poisson(World& world, const double lo, const double econv = FunctionDefaults<3>::get_thresh()) {
        return std::shared_ptr<real_convolution_3d>(CoulombOperatorPtr(world, lo, econv));
    }

    /// apply the exchange operator on a vector of functions

    /// note that only one spin is used (either alpha or beta orbitals)
    /// @param[in]  vket       the orbitals |i> that the operator is applied on
    /// @return     a vector of orbitals  K| i>
    vecfuncT operator()(const vecfuncT& vket) const;

    bool is_symmetric() const { return symmetric_; }

    ExchangeImpl& set_taskq(std::shared_ptr<MacroTaskQ> taskq1) {
        this->taskq=taskq1;
        return *this;
    }

    ExchangeImpl& symmetric(const bool flag) {
        symmetric_ = flag;
        return *this;
    }

    ExchangeImpl& set_macro_task_info(const MacroTaskInfo& info) {
        macro_task_info = info;
        return *this;
    }

    ExchangeImpl& set_macro_task_info(const std::vector<std::string>& info) {
        macro_task_info.from_vector_of_strings(info);
        if (world.rank() == 0 && printdebug()) {
            print("set macrotaskinfo to");
            print(macro_task_info);
        }
        return *this;
    }

    ExchangeImpl& set_algorithm(const Algorithm& alg) {
        algorithm_ = alg;
        return *this;
    }

    ExchangeImpl& set_printlevel(const long& level) {
        printlevel=level;
        return *this;
    }

    ExchangeImpl& set_batch_granularity(const long level) {
        MADNESS_CHECK_THROW(level >= 1, "exchange batch granularity must be at least 1");
        batch_granularity_ = level;
        return *this;
    }

    ExchangeImpl& set_cost_aware_assignment(const bool flag) {
        cost_aware_assign_ = flag;
        return *this;
    }

    ExchangeImpl& set_accumulation_mode(const int mode) {
        MADNESS_CHECK_THROW(mode == 1 or mode == 2, "exchange accumulation mode must be 1 or 2");
        accumulation_mode_ = mode;
        return *this;
    }

    std::shared_ptr<MacroTaskQ> get_taskq() const {return taskq;}

    World& get_world() const {return world;}

    nlohmann::json get_statistics() const {return statistics;}

    /// return some statistics about the current settings
    nlohmann::json gather_statistics() const {
        nlohmann::json j;
        j["symmetric"] = symmetric_;
        j["lo"] = lo;
        j["thresh"] = thresh;
        j["mul_tol"] = mul_tol;
        j["printlevel"] = printlevel;
        j["algorithm"] = to_string(algorithm_);
        j["macro_task_info"] = macro_task_info.to_json();
        auto timings = gather_timings(world);
        j.update(timings);
        return j;
    }

private:

    /// exchange using macrotasks, i.e. apply K on a function in individual worlds
    vecfuncT K_macrotask_efficient(const vecfuncT& vket, const double mul_tol = 0.0) const;

    /// exchange using macrotasks, i.e. apply K on a function in individual worlds row-wise
    vecfuncT K_macrotask_efficient_row(const vecfuncT& vket, const double mul_tol = 0.0) const;

    /// computing the full square of the double sum (over vket and the K orbitals)
    vecfuncT K_small_memory(const vecfuncT& vket, const double mul_tol = 0.0) const;

    /// computing the upper triangle of the double sum (over vket and the K orbitals)
    vecfuncT K_large_memory(const vecfuncT& vket, const double mul_tol = 0.0) const;

    /// computing the upper triangle of the double sum (over vket and the K orbitals)
    static vecfuncT compute_K_tile(World& world, const vecfuncT& mo_bra, const vecfuncT& mo_ket,
                                   const vecfuncT& vket, std::shared_ptr<real_convolution_3d> poisson,
                                   const bool symmetric, const double mul_tol = 0.0);

    inline bool printdebug() const {return printlevel >= 10; }
    inline bool printprogress() const {return (printlevel>=4) and (not (printdebug()));}
    inline bool printtimings() const {return printlevel>=3;}
    inline bool printtimings_detail() const {return printlevel>=4;}

    World& world;
    std::shared_ptr<MacroTaskQ> taskq;
    bool symmetric_ = false;      /// is the exchange matrix symmetric? K phi_i = \sum_k \phi_k \int \phi_k \phi_i
    vecfuncT mo_bra, mo_ket;    ///< MOs for bra and ket
    double lo = 1.e-4;
    double thresh = FunctionDefaults<NDIM>::get_thresh();
    long printlevel = 0;
    double mul_tol = FunctionDefaults<NDIM>::get_thresh()*0.1;
    /// batches per rank in the owner-pinned symmetric partition; 1 is the coarsest split
    /// and the one with the lowest peak memory
    long batch_granularity_ = 1;
    /// how the tile results are gathered: 1 = subworld buffer then universe, 2 = also reduce
    /// within a node first (default; degrades to 1 on a single node)
    int accumulation_mode_ = 2;
    /// place tasks by measured cost instead of by counting them
    bool cost_aware_assign_ = true;

    mutable nlohmann::json statistics;  ///< statistics of the Cloud (timings, memory)  and of the parameters of this run

    class MacroTaskExchangeSimple : public MacroTaskOperationBase {

        long nresult;
        double lo = 1.e-4;
        double mul_tol = 1.e-7;
        bool symmetric = false;
        /// pin each task to a rank that owns one of its batches, over the owner-pinned split
        bool owner_pinned = false;
        long granularity_level = 1;
        /// Per-application salt for the batch record keys, from exchange_batch_salt.

        /// Carried on the task rather than re-derived where it is used: every rank builds the
        /// same task objects collectively, so a constructor argument already reaches every
        /// rank, and deriving it instead requires whichever operand vector it is taken from to
        /// be available in full wherever a key is formed. That is true only while the ket is
        /// passed unbatched.
        long batch_salt_ = 0;
        /// 1 = sum into a subworld buffer and drain that into the universe;
        /// 2 = additionally reduce within a node first, so only one rank per node scatters
        ///     across nodes. Degrades to 1 automatically when there is a single node.
        int accumulation_mode_ = 2;
        /// place tasks by their measured cost rather than by counting them
        bool cost_aware_ = true;
        /// (column batch offset, row batch offset) -> owning rank, filled by
        /// prepare_owner_assignment and read by owner_hint
        std::map<std::pair<long,long>,long> owner_map_;
        /// this process's rank in the universe, to recognise the batches it owns
        long universe_rank_ = 0;

        /// Batches fetched from the cloud, reused across the tasks that run in one subworld.

        /// Static because the tasks of a subworld are separate objects: the cache has to
        /// outlive any one of them to be reused at all. It is therefore scoped by hand to
        /// the subworld that filled it -- see ensure_cache_world, which drops it when the
        /// subworld changes so a batch from one subworld is never read in another.
        static inline ExchangeBatchLRU<long, vecfuncT> batch_cache_;
        static inline long cache_world_id_ = -1;
        static inline std::atomic<long> batch_cache_hits_;
        static inline std::atomic<long> batch_cache_misses_;
        static inline std::atomic<long> batch_prefetch_hits_;

        /// One batch requested ahead of the task that will read it.

        /// A strict double buffer: `next` is requested while this task computes, and every
        /// task promotes it to `current` before computing, so **at most two requests are ever
        /// in flight per rank**. That bound is load-bearing, not a tuning choice: an earlier
        /// design that allowed several concurrent requests corrupted the transport on large
        /// replies, because the outstanding replies could outnumber what the receive path was
        /// prepared for.
        struct PrefetchSlot {
            bool valid = false;
            long key = 0;
            Future<batch_bytesT> fut;
        };
        static inline PrefetchSlot prefetch_current_;   ///< promoted from the previous task
        static inline PrefetchSlot prefetch_next_;      ///< requested during this task

        /// What each task cost last time, to place them better this time.

        /// Screening makes the tiles strongly uneven for large molecules, and a placement that
        /// balances task *counts* cannot see that. The measured wall time of each tile is
        /// recorded here, summed across ranks after the call -- each tile ran on exactly one
        /// rank, so summing unions the contributions -- and used as the reference for the next
        /// call. Kept across calls and across protocol changes on purpose: only the relative
        /// cost matters, and its structure barely moves between them.
        static inline std::vector<double> cost_reference_;   ///< from the previous call
        static inline std::vector<double> cost_this_call_;    ///< rank-local, summed after the call
        static inline std::map<long,long> batch_begin_to_index_;   ///< batch offset -> index
        static inline long exchange_call_index_ = 0;

        /// This task's profile record. NOT static: it belongs to the single operator() call on
        /// this object, and the fetch writes into it through `this` during that call.
        mutable ExchTaskProfile prof_;
        static inline long prof_task_seq_ = 0;   ///< per-process task counter, for identity only

        /// Where this subworld's tile results are summed before they leave it, and where a
        /// node's subworlds are summed before that leaves the node. Static for the same reason
        /// the batch cache is -- each task batch is a separate object -- and so subject to the
        /// same two lifetime rules: the world-id guards below stop a stale read, and cleanup()
        /// releases them while their world is still alive. Neither alone is enough.
        static inline vecfuncT Kf_local_;
        static inline bool Kf_local_initialized_ = false;
        static inline long Kf_local_world_id_ = -1;
        static inline vecfuncT Kf_node_;
        static inline bool Kf_node_initialized_ = false;
        static inline long Kf_node_world_id_ = -1;
        /// Receiving endpoints for the two drains, one per world they transfer within.

        /// These are WorldObjects, so they are bound to a world exactly as a Function is, and the
        /// same lifetime rule applies: **release them while that world is still alive.** The node
        /// world is built per application (the queue owns it and a fresh queue is built per call),
        /// so a cached node reducer that survives the application is registered in a world that no
        /// longer exists. The world id is kept beside the pointer rather than read back out of it,
        /// because comparing `reducer->get_world().id()` is itself a read of the dead world.
        static inline std::shared_ptr<FinalizeReducer<T, NDIM>> universe_reducer_;
        static inline long universe_reducer_world_id_ = -1;
        static inline std::shared_ptr<FinalizeReducer<T, NDIM>> node_reducer_;
        static inline long node_reducer_world_id_ = -1;
        /// each drain happens once per rank, not once per task object
        static inline bool finalize_stage1_done_ = false;
        static inline bool finalize_universe_done_ = false;

        static void clear_local_caches() {
            batch_cache_.clear();
            // drop requests issued in a subworld that is no longer the one we are in
            prefetch_current_ = PrefetchSlot();
            prefetch_next_ = PrefetchSlot();
        }

        /// entries per chunk in coalesced_gaxpy, sized so one message stays modest at any k
        static std::size_t finalize_chunk_entries() {
            const std::size_t k = FunctionDefaults<NDIM>::get_k();
            const std::size_t per_node = std::size_t(1) << (2 * NDIM);   // (2k)^NDIM / k^NDIM bound
            return std::max<std::size_t>(1, (1u << 20) / (per_node * k * k * k * sizeof(T) + 1));
        }

        /// One reducer per transport world, rebuilt when that world changes. Collective:
        /// constructing a WorldObject is, so every rank of `world` reaches this together.
        static FinalizeReducer<T, NDIM>& get_universe_reducer(World& world) {
            if (not universe_reducer_ or universe_reducer_world_id_ != long(world.id())) {
                universe_reducer_ = std::make_shared<FinalizeReducer<T, NDIM>>(world);
                universe_reducer_world_id_ = long(world.id());
            }
            return *universe_reducer_;
        }

        static FinalizeReducer<T, NDIM>& get_node_reducer(World& world) {
            if (not node_reducer_ or node_reducer_world_id_ != long(world.id())) {
                node_reducer_ = std::make_shared<FinalizeReducer<T, NDIM>>(world);
                node_reducer_world_id_ = long(world.id());
            }
            return *node_reducer_;
        }

        /// drop cached batches when the subworld changes; they belong to the old one
        void ensure_cache_world(World& world) const {
            if (cache_world_id_ != long(world.id())) {
                clear_local_caches();
                cache_world_id_ = long(world.id());
            }
        }

        /// Fetch one owner-pinned batch, from the local cache if it is resident.

        /// A miss goes straight to the owning rank over the cloud's point-to-point batch
        /// path. Batches this rank owns are pinned, since every one of its tasks needs
        /// them; the others are transient and bounded.
        ///
        /// \return a reference into the cache, valid until that entry is evicted
        const vecfuncT& fetch_batch(World& world, Cloud& cloud, const long record) const {
            ensure_cache_world(world);
            if (const vecfuncT* resident = batch_cache_.find(record)) {
                ++batch_cache_hits_;
                if (profile_active()) prof_.observe_fetch_tier(0);
                return *resident;
            }
            const bool owned = (cloud.batch_owner(record) == universe_rank_);
            // requested one task ahead, so the transfer overlapped that task's compute
            for (PrefetchSlot* slot : {&prefetch_current_, &prefetch_next_}) {
                if (slot->valid and slot->key == record) {
                    vecfuncT data = cloud.template deserialize_batch_p2p<T, NDIM>(
                            world, slot->fut, record, /*cache_result=*/false);
                    *slot = PrefetchSlot();
                    ++batch_prefetch_hits_;
                    if (profile_active()) prof_.observe_fetch_tier(1);
                    return batch_cache_.insert(record, std::move(data), owned);
                }
            }
            ++batch_cache_misses_;
            if (profile_active()) prof_.observe_fetch_tier(2);
            vecfuncT data = cloud.template fetch_batch_p2p<T, NDIM>(world, record, /*cache_result=*/false);
            return batch_cache_.insert(record, std::move(data), owned);
        }

        /// custom partitioning for the exchange operator in exchangeoperator.h

        /// arguments are: result[i] += sum_k vket[k] \int 1/r vbra[k] f[i]
        /// with f and vbra being batched, result and vket being passed on as a whole
        class MacroTaskPartitionerExchange : public MacroTaskPartitioner {
        public:
            MacroTaskPartitionerExchange(const bool symmetric, const bool owner_pinned = false,
                                         const long granularity_level = 1)
                    : symmetric(symmetric), owner_pinned(owner_pinned),
                      granularity_level(granularity_level) {
                max_batch_size=30;
            }

            bool symmetric = false;
            /// build the grid from the owner-pinned split instead of the size-driven one,
            /// so a task's two batches coincide with batches some rank owns
            bool owner_pinned = false;
            long granularity_level = 1;

            partitionT do_partitioning(const std::size_t& vsize1, const std::size_t& vsize2,
                                       const std::string policy) const override {

                if (owner_pinned and not symmetric) {
                    // full grid over two independent splits; see exchange_row_owner_grid
                    partitionT result;
                    for (const auto& [column, row] : exchange_row_owner_grid(vsize1, vsize2,
                                                                            long(nsubworld))) {
                        Batch batch(column, row, _);
                        result.push_back(std::make_pair(batch, compute_priority(batch)));
                    }
                    return result;
                }

                if (owner_pinned) {
                    // lower-triangular grid over one granularity-aware split, shared with
                    // the owner assignment: input[0] = batch i (column), input[1] = batch j
                    // (row), j <= i. Owners are assigned by prepare_owner_assignment.
                    const std::vector<Batch_1D> batches =
                            exchange_sym_owner_split(vsize1, long(nsubworld), granularity_level);
                    partitionT result;
                    for (long i = 0; i < long(batches.size()); ++i) {
                        for (long j = 0; j <= i; ++j) {
                            Batch batch(batches[i], batches[j], _);
                            result.push_back(std::make_pair(batch, compute_priority(batch)));
                        }
                    }
                    return result;
                }

                partitionT partition1 = do_1d_partition(vsize1, policy);
                partitionT partition2 = do_1d_partition(vsize2, policy);
                partitionT result;
                for (auto i = partition1.begin(); i != partition1.end(); ++i) {
                    if (symmetric) {
                        for (auto j = i; j != partition1.end(); ++j) {
                            Batch batch(i->first.input[0], j->first.input[0], _);
                            double priority=compute_priority(batch);
                            result.push_back(std::make_pair(batch,priority));
                        }
                    } else {
                        for (auto j = partition2.begin(); j != partition2.end(); ++j) {
                            Batch batch(i->first.input[0], j->first.input[0], _);
                            double priority=compute_priority(batch);
                            result.push_back(std::make_pair(batch,priority));
                        }
                    }
                }
                return result;
            }

            /// compute the priority of this task for non-dumb scheduling

            /// \return the priority as double number (no limits)
            double compute_priority(const Batch& batch) const override {
                MADNESS_CHECK(batch.input.size() == 2);   // must be quadratic batches
                long nrow = batch.input[0].size();
                long ncol = batch.input[1].size();
                return double(nrow * ncol);
            }
        };

    public:
        MacroTaskExchangeSimple(const long nresult, const double lo, const double mul_tol,
                                const bool symmetric, const bool owner_pinned = false,
                                const long granularity_level = 1, const long universe_rank = 0,
                                const int accumulation_mode = 2, const bool cost_aware = true,
                                const long batch_salt = 0)
                : nresult(nresult), lo(lo), mul_tol(mul_tol), symmetric(symmetric),
                  owner_pinned(owner_pinned), granularity_level(granularity_level),
                  batch_salt_(batch_salt),
                  accumulation_mode_(accumulation_mode), cost_aware_(cost_aware),
                  universe_rank_(universe_rank) {
            partitioner.reset(new MacroTaskPartitionerExchange(symmetric, owner_pinned, granularity_level));
            name="MacroTaskExchangeSimple";
        }

        /// how often a task's operand batch was already resident, and how often it had to be
        /// fetched from the rank owning it. No fetches at all means the path never ran.
        static long batch_cache_hits() { return batch_cache_hits_; }
        static long batch_cache_misses() { return batch_cache_misses_; }
        static long batch_prefetch_hits() { return batch_prefetch_hits_; }
        /// per-task profiling on for this task?
        bool profile_active() const { return owner_pinned and exch_task_profile_enabled(); }
        static std::vector<double>& cost_this_call() { return cost_this_call_; }
        static long exchange_call_index() { return exchange_call_index_; }
        /// number of batches this application split into, which squares to the cost matrix
        static long cost_matrix_dimension() { return long(batch_begin_to_index_.size()); }
        /// make this call's measured costs the reference for the next one
        static void commit_cost_reference() { cost_reference_ = cost_this_call_; }
        static void reset_batch_cache_counters() {
            batch_cache_hits_ = 0l; batch_cache_misses_ = 0l; batch_prefetch_hits_ = 0l;
        }

        /// Assign every task to the rank that will own one of its two batches.

        /// Called by the macrotask queue after partitioning and before it asks for each
        /// task's owner. The batch boundaries come from the same split the partitioner
        /// used, so a task's (column, row) batch offsets identify a pair of batch indices,
        /// and exchange_sym_round_robin_assign turns that pair into an owner. Every rank
        /// runs this over the same partition and gets the same map without communicating.
        void prepare_owner_assignment(const MacroTaskPartitioner::partitionT& partition,
                                      const long nsubworld) {
            owner_map_.clear();
            if (not owner_pinned or nsubworld <= 0) return;

            if (not symmetric) {
                // Both dimensions are indexed from the partition rather than from nresult, which is
                // the column length only: the row dimension carries bra and ket and has its own.
                std::map<long,long> column_index, row_index;
                for (const auto& [task_batch, priority] : partition) {
                    MADNESS_CHECK_THROW(task_batch.input.size() >= 2,
                                        "owner-pinned exchange expects two-dimensional task batches");
                    column_index[task_batch.input[0].begin] = 0;   // renumbered below, in order
                    row_index[task_batch.input[1].begin] = 0;
                }
                long next = 0;
                for (auto& [begin, index] : column_index) index = next++;
                next = 0;
                for (auto& [begin, index] : row_index) index = next++;

                const auto assignment = exchange_row_owner_assign(long(column_index.size()),
                                                                  long(row_index.size()), nsubworld);
                for (const auto& [task_batch, priority] : partition) {
                    const long c = column_index[task_batch.input[0].begin];
                    const long r = row_index[task_batch.input[1].begin];
                    auto it = assignment.find({c, r});
                    owner_map_[{task_batch.input[0].begin, task_batch.input[1].begin}] =
                            (it != assignment.end()) ? it->second : (c % nsubworld);
                }
                // Cost-aware placement stays off here: its vector is indexed by exchange_sym_tri,
                // which is meaningless for a rectangle, and a column-per-worker grid is already
                // balanced. An empty vector is what stops the tiles recording into it.
                cost_this_call_.clear();
                batch_begin_to_index_.clear();
                return;
            }

            const std::vector<Batch_1D> split =
                    exchange_sym_owner_split(nresult, nsubworld, granularity_level);
            const long M = long(split.size());
            std::map<long,long> begin_to_index;
            for (long k = 0; k < M; ++k) begin_to_index[split[k].begin] = k;

            // Cost-aware placement needs a reference from a previous call, and the first call
            // is not representative of the ones that follow it -- cold caches and an initial
            // guess -- so it takes effect from the third call on, with the second call's
            // measurements as its reference.
            ++exchange_call_index_;
            batch_begin_to_index_ = begin_to_index;
            const std::size_t ntask = std::size_t(M) * std::size_t(M + 1) / 2;
            cost_this_call_.assign(ntask, 0.0);
            const bool use_cost = cost_aware_ and exchange_call_index_ >= 3
                    and cost_reference_.size() == ntask;
            const std::map<std::pair<long,long>,long> assignment =
                    use_cost ? exchange_sym_cost_aware_assign(nsubworld, M, cost_reference_)
                             : exchange_sym_round_robin_assign(nsubworld, M);

            for (const auto& [task_batch, priority] : partition) {
                MADNESS_CHECK_THROW(task_batch.input.size() >= 2,
                                    "owner-pinned exchange expects two-dimensional task batches");
                const long column_begin = task_batch.input[0].begin;
                const long row_begin = task_batch.input[1].begin;
                auto ic = begin_to_index.find(column_begin);
                auto jr = begin_to_index.find(row_begin);
                MADNESS_CHECK_THROW(ic != begin_to_index.end() and jr != begin_to_index.end(),
                                    "owner-pinned exchange: a task batch is not one of the split batches");
                // the assignment is keyed on the lower triangle, so order the pair
                const long i = std::max(ic->second, jr->second);
                const long j = std::min(ic->second, jr->second);
                auto it = assignment.find({i, j});
                owner_map_[{column_begin, row_begin}] =
                        (it != assignment.end()) ? it->second : (i % nsubworld);
            }
        }

        /// \return the rank this task is pinned to, or -1 to leave the choice to the queue
        long owner_hint(const Batch& task_batch, const long nsubworld) const override {
            if (not owner_pinned or nsubworld <= 0 or owner_map_.empty()) return -1;
            MADNESS_CHECK_THROW(task_batch.input.size() >= 2,
                                "owner-pinned exchange expects two-dimensional task batches");
            auto it = owner_map_.find({task_batch.input[0].begin, task_batch.input[1].begin});
            return (it != owner_map_.end()) ? it->second : -1;
        }


        // you need to define the exact argument(s) of operator() as tuple
        typedef std::tuple<const std::vector<Function<T, NDIM>>&,
                const std::vector<Function<T, NDIM>>&,
                const std::vector<Function<T, NDIM>>&> argtupleT;

        using resultT = std::vector<Function<T, NDIM>>;

        // you need to define an empty constructor for the result
        // resultT must implement operator+=(const resultT&)
        resultT allocator(World& world, const argtupleT& argtuple) const {
            std::size_t n = std::get<0>(argtuple).size();
            resultT result = zero_functions_compressed<T, NDIM>(world, n);
            return result;
        }

        /// Store the orbitals as owner-pinned batches, one record per batch.

        /// Called by the macrotask queue on the universe right after the argument tuple is
        /// stored. The batch boundaries and the record keys are derived exactly as the task
        /// side derives them, so no manifest has to be communicated.
        ///
        /// Every rank registers the routing for all records, which is local and needs no
        /// communication, and then each owner **pulls** the batches it owns into its own
        /// size-1 subworld and serializes them there. That is what spreads the ingest across
        /// the owners: serializing centrally instead funnels the whole orbital set through
        /// one rank's network interface.
        ///
        /// The symmetric case stores one set: bra == ket == vf, so the same batch serves as
        /// column, row and ket operand. The asymmetric case stores three, one per role, over two
        /// splits -- vf along the column dimension, bra and ket along the row dimension they share.
        void store_batches(World& world, World& subworld, Cloud& cloud, const argtupleT& argtuple,
                           const long nsubworld) {
            if (not owner_pinned) return;
            // Batch k is owned by rank (k mod nsubworld), so a batch index has to name a
            // universe rank, and that only holds with one subworld per rank. Anything else
            // would register the routing to the wrong ranks and read the wrong coefficients,
            // so say so rather than compute something wrong.
            MADNESS_CHECK_THROW(nsubworld == world.size(),
                                "owner-pinned exchange needs one subworld per rank");
            const vecfuncT& vf = std::get<0>(argtuple);
            const vecfuncT& mo_bra = std::get<1>(argtuple);
            const vecfuncT& mo_ket = std::get<2>(argtuple);
            // one fence up front, so store_batch does not need one per function
            world.gop.fence();

            // A cross-world copy picks its process map from the target world, but anything
            // on that path reading the process-wide default would route to universe ranks
            // that do not exist in a size-1 subworld. Point the default at the subworld for
            // the duration and restore it afterwards.
            auto saved_pmap = FunctionDefaults<NDIM>::get_pmap();
            FunctionDefaults<NDIM>::set_default_pmap(subworld);

            std::vector<std::pair<long, vecfuncT>> owned;
            // register the routing for one operand set, and pull the batches this rank owns. Batch k
            // goes to rank k, which is what puts a column batch on the rank running that column.
            auto stage = [&](const vecfuncT& v, const int dim, const std::vector<Batch_1D>& split) {
                for (long k = 0; k < long(split.size()); ++k) {
                    const Batch_1D& r = split[k];
                    const long record = exchange_batch_record_key(batch_salt_, dim, r);
                    cloud.register_batch_owner(record, ProcessID(k % nsubworld));
                    if (k % nsubworld == world.rank()) {
                        vecfuncT local(r.size());
                        for (long i = r.begin; i < r.end; ++i)
                            local[i - r.begin] = copy(subworld, v[i], /*fence=*/false);
                        owned.emplace_back(record, std::move(local));
                    }
                }
            };

            if (symmetric) {
                stage(mo_ket, EXCHANGE_BATCH_VF,
                      exchange_sym_owner_split(mo_ket.size(), nsubworld, granularity_level));
            } else {
                // bra and ket are indexed together by the operator's sum over pairs, so they share
                // the row boundaries; vf is independent and gets the column ones
                MADNESS_CHECK_THROW(mo_bra.size() == mo_ket.size(),
                                    "exchange: bra and ket are a paired set and must have equal length");
                const std::vector<Batch_1D> columns = exchange_row_owner_split(vf.size(), nsubworld);
                const std::vector<Batch_1D> rows = exchange_row_owner_split(mo_bra.size(), nsubworld);
                stage(vf, EXCHANGE_BATCH_VF, columns);
                stage(mo_bra, EXCHANGE_BATCH_BRA, rows);
                stage(mo_ket, EXCHANGE_BATCH_KET, rows);
            }
            // the source ranks' comm threads serve the pulls, so a fence on the size-1
            // subworld drains this owner's copies
            subworld.gop.fence();
            for (auto& [record, batch] : owned)
                cloud.store_batch(subworld, batch, world.rank(), record, /*fence=*/false);
            subworld.gop.fence();

            FunctionDefaults<NDIM>::set_pmap(saved_pmap);
            world.gop.fence();
        }

        /// Request the batch the next task will have to fetch, before computing this one.

        /// Called by the queue once per task, before the task body, with the batches of the
        /// next task this rank will run. Of that task's two batches one is normally owned here
        /// and reads locally, so at most one is worth requesting -- and requesting exactly one
        /// keeps the in-flight count within the bound PrefetchSlot documents.
        ///
        /// This is what makes the owner-pinned transport worth its machinery: without it every
        /// task pays the full latency of its remote batch with nothing to overlap it against.
        /// \param mo_ket  unused: the salt it used to be derived from is carried on the task.
        ///                It stays in the signature because the queue's detection trait matches
        ///                on it (has_sym_pipeline_advance_v).
        void sym_pipeline_advance(World& subworld, const vecfuncT&,
                                  const Batch_1D& next_col, const Batch_1D& next_row,
                                  const bool has_next) const {
            if (not owner_pinned) return;
            ensure_cache_world(subworld);
            // hand the previous task's request to this task, which is the one that reads it.
            // Retiring an unconsumed request is safe rather than merely tolerated: the reply
            // still lands, fills a result nobody reads, and the transport drops its pending
            // entry, so the cost is one wasted transfer and nothing is left dangling.
            prefetch_current_ = prefetch_next_;
            prefetch_next_ = PrefetchSlot();
            if (not has_next or cloud_ptr == nullptr) return;
            Cloud& cloud = *cloud_ptr;
            const long salt = batch_salt_;
            // Symmetric: either of the task's two ranges may be the remote one, so try both. In the
            // asymmetric case only the row rotates -- the column is held on this rank -- and of the
            // two records it carries, the bra is the one requested; the ket follows synchronously.
            const std::vector<std::pair<Batch_1D,int>> candidates =
                    symmetric ? std::vector<std::pair<Batch_1D,int>>{{next_col, EXCHANGE_BATCH_VF},
                                                                     {next_row, EXCHANGE_BATCH_VF}}
                              : std::vector<std::pair<Batch_1D,int>>{{next_row, EXCHANGE_BATCH_BRA}};
            for (const auto& [r, dim] : candidates) {
                const long record = exchange_batch_record_key(salt, dim, r);
                if (cloud.batch_owner(record) == universe_rank_) continue;   // reads locally
                if (batch_cache_.contains(record)) continue;                 // already resident
                if (prefetch_current_.valid and prefetch_current_.key == record) continue;
                prefetch_next_.key = record;
                prefetch_next_.fut = cloud.request_batch_bytes_async(record);
                prefetch_next_.valid = true;
                break;                                                      // one per task
            }
        }

        /// the owner-pinned path fetches its operand batches from the cloud itself
        bool handles_own_data_movement() const override { return owner_pinned; }

        /// Drop the cached batches while the subworld holding them is still alive.

        /// Leaving it to ensure_cache_world to notice a new subworld is too late: by then the
        /// cached functions refer to a destroyed world, and merely releasing them walks into
        /// it.
        void cleanup() override {
            clear_local_caches();
            cache_world_id_ = -1;
            // released here, while the subworld and node world still exist
            Kf_local_.clear();
            Kf_local_initialized_ = false;
            Kf_local_world_id_ = -1;
            Kf_node_.clear();
            Kf_node_initialized_ = false;
            Kf_node_world_id_ = -1;
            finalize_stage1_done_ = false;
            finalize_universe_done_ = false;
            // Both reducers go too, for the reason given at their declaration. Rebuilding one per
            // application costs a WorldObject construction, which is nothing beside an
            // application, and it also removes the shutdown hazard of a static outliving its
            // world.
            universe_reducer_.reset();
            universe_reducer_world_id_ = -1;
            node_reducer_.reset();
            node_reducer_world_id_ = -1;
            // cost_reference_ deliberately survives: it is what the next call places against.
            // cost_this_call_ survives too, until the caller has summed and committed it.
        }

        /// true if the task sums its own tile results and drains them in the finalize, rather
        /// than the queue moving every tile result into the universe result by itself
        bool accumulates_own_output() const override { return owner_pinned; }

        /// true if the drain goes subworld -> node -> universe rather than straight to the
        /// universe, so only one rank per node scatters across nodes
        /// not an override: the queue reaches this through its optional-hook detection, since
        /// the virtual it feeds lives on the internal task rather than on this base
        bool wants_node_local_reduction() const {
            return owner_pinned and accumulation_mode_ == 2;
        }

        /// The result entries this tile actually wrote.

        /// operator() scatters into a full-width Kf and leaves everything else zero, so summing all
        /// `nresult` entries would gaxpy mostly zeros -- a per-tile cost proportional to the whole
        /// result vector. Every tile writes its column range; a symmetric off-diagonal tile also
        /// writes its row range, reusing each intermediate for the transposed element, whereas an
        /// asymmetric tile contributes to its column alone. A full-size or absent range falls back
        /// to all of them.
        ///
        /// The result must be a *set*: the caller gaxpys one entry per index, so a repeated index is
        /// added twice. The two ranges do overlap in the asymmetric case, coming from separate splits
        /// of different lengths rather than from one split, where they are always equal or disjoint.
        std::vector<long> touched_result_indices() const {
            std::vector<long> idx;
            auto add_range = [&](const Batch_1D& b) {
                const long s = b.is_full_size() ? 0 : b.begin;
                const long e = b.is_full_size() ? long(nresult) : b.end;
                for (long i = s; i < e and i < long(nresult); ++i) idx.push_back(i);
            };
            if (batch.input.empty()) {
                for (long i = 0; i < long(nresult); ++i) idx.push_back(i);
                return idx;
            }
            add_range(batch.input[0]);
            if (symmetric and batch.input.size() > 1 and not (batch.input[1] == batch.input[0]))
                add_range(batch.input[1]);
            std::sort(idx.begin(), idx.end());
            idx.erase(std::unique(idx.begin(), idx.end()), idx.end());
            return idx;
        }

        /// sum one tile's result into this subworld's accumulator
        void accumulate_locally(World& subworld, const vecfuncT& result_subworld) const {
            const long wid = long(subworld.id());
            if (not Kf_local_initialized_ or Kf_local_world_id_ != wid) {
                Kf_local_ = zero_functions_compressed<T, NDIM>(subworld, nresult);
                Kf_local_initialized_ = true;
                Kf_local_world_id_ = wid;
            }
            // shallow handles for just the entries this tile wrote; the rest are structurally
            // zero, so skipping them is not an approximation
            const std::vector<long> touched = touched_result_indices();
            vecfuncT rs_sub, kf_sub;
            rs_sub.reserve(touched.size());
            kf_sub.reserve(touched.size());
            for (long i : touched) { rs_sub.push_back(result_subworld[i]); kf_sub.push_back(Kf_local_[i]); }
            if (rs_sub.empty()) return;
            const TreeState op_state =
                    rs_sub[0].get_impl()->get_tensor_type() == TT_FULL ? compressed : reconstructed;
            change_tree_state(rs_sub, op_state);
            gaxpy(1.0, kf_sub, 1.0, rs_sub, false);   // mutates kf_sub[k] == Kf_local_[touched[k]]
        }

        /// Collectively (re)build the node-shared accumulator in the node world.

        /// Must be reached by every rank of `nodeworld`, since constructing a Function is
        /// collective; the queue drives this uniformly across the replicated task list, so the
        /// initialized flag flips in lockstep. The process map is passed explicitly because the
        /// process-wide default is the subworld's during the finalize -- inheriting it would map
        /// keys to subworld rank indices for functions that live in the node world.
        void ensure_node_accumulator(World& nodeworld) const {
            const long wid = long(nodeworld.id());
            if (Kf_node_initialized_ and Kf_node_world_id_ == wid) return;
            auto node_pmap = FunctionDefaults<NDIM>::make_default_pmap(nodeworld);
            Kf_node_.resize(nresult);
            for (long i = 0; i < nresult; ++i)
                Kf_node_[i] = functionT(FunctionFactory<T, NDIM>(nodeworld)
                                                .pmap(node_pmap)
                                                .compressed(true)
                                                .fence(false));
            nodeworld.gop.fence();
            Kf_node_initialized_ = true;
            Kf_node_world_id_ = wid;
        }

        /// reduce this subworld's accumulator into the node-shared one
        void finalize_stage1(World& subworld, World* nodeworld) {
            if (not nodeworld) return;                 // one node: stage 2 drains directly
            if (finalize_stage1_done_) return;
            finalize_stage1_done_ = true;
            ensure_node_accumulator(*nodeworld);       // collective on the node
            if (Kf_local_initialized_) change_tree_state(Kf_local_, compressed);
            vecfuncT empty;
            vecfuncT& src = Kf_local_initialized_ ? Kf_local_ : empty;
            coalesced_gaxpy<T, NDIM>(*nodeworld, get_node_reducer(*nodeworld),
                                     Kf_node_, src, T(1.0), finalize_chunk_entries());
        }

        /// drain into the universe result, from the node accumulator if there is one
        void finalize_stage2(World& subworld, World* nodeworld, vecfuncT& universe_result) {
            if (finalize_universe_done_) return;
            finalize_universe_done_ = true;
            if (universe_result.empty()) return;
            World& universe = universe_result.front().world();
            vecfuncT empty;
            vecfuncT* src = &empty;
            if (nodeworld) {
                if (Kf_node_initialized_) src = &Kf_node_;
            } else {
                if (Kf_local_initialized_) {
                    change_tree_state(Kf_local_, compressed);
                    src = &Kf_local_;
                }
            }
            // every rank reaches this, including ones with nothing to send -- coalesced_gaxpy
            // is collective on the universe
            coalesced_gaxpy<T, NDIM>(universe, get_universe_reducer(universe),
                                     universe_result, *src, T(1.0), finalize_chunk_entries());
        }


        std::vector<Function<T, NDIM>>
        operator()(const std::vector<Function<T, NDIM>>& vf_batch,     // will be batched (column)
                   const std::vector<Function<T, NDIM>>& bra_batch,    // will be batched (row)
                   const std::vector<Function<T, NDIM>>& vket) {       // will not be batched

            // the operands are not necessarily this world's: on the owner-pinned path the
            // queue leaves them in the universe and the task fetches what it needs
            MADNESS_CHECK_THROW(subworld_ptr != nullptr, "MacroTaskExchangeSimple: subworld_ptr is null");
            World& world = *subworld_ptr;
            resultT Kf = zero_functions_compressed<T, NDIM>(world, nresult);

            bool diagonal_block = batch.input[0] == batch.input[1];
            auto& bra_range = batch.input[1];    // corresponds to vbra
            auto& vf_range = batch.input[0];       // corresponds to vf_batch

            if (vf_range.is_full_size()) vf_range.end = vf_batch.size();
            if (bra_range.is_full_size()) bra_range.end = bra_batch.size();

            MADNESS_CHECK(vf_range.end <= nresult);
            if (symmetric) MADNESS_CHECK(bra_range.end <= nresult);

            const double tile_wall_start = wall_time();
            const bool profiling = profile_active();
            if (profiling) {
                prof_.reset();
                prof_.task_id = prof_task_seq_++;
                prof_.universe_rank = universe_rank_;
                prof_.subworld_id = static_cast<unsigned long>(world.id());
                prof_.subworld_nrank = world.size();
                prof_.thresh = FunctionDefaults<NDIM>::get_thresh();
                prof_.k = FunctionDefaults<NDIM>::get_k();
                prof_.diagonal = diagonal_block;
                prof_.row_begin = bra_range.begin; prof_.row_end = bra_range.end;
                prof_.col_begin = vf_range.begin;  prof_.col_end = vf_range.end;
                prof_.wall_start = tile_wall_start;
            }

            // Owner-pinned: fetch the two batches this task needs from the cloud. Because
            // bra == ket == vf there, one batch per range serves every operand role, so the
            // ket over a range is just that range's batch.
            vecfuncT bra_owned, vf_owned, ket_owned;
            const vecfuncT* bra_work = &bra_batch;
            const vecfuncT* vf_work = &vf_batch;
            const vecfuncT* ket_work = nullptr;   // set only when the ket is fetched as its own role
            if (owner_pinned) {
                MADNESS_CHECK_THROW(cloud_ptr != nullptr, "owner-pinned exchange: cloud_ptr is null");
                Cloud& cloud = *cloud_ptr;
                const long salt = batch_salt_;
                // copy the handles out of the cache: the reference is only good until the
                // entry is evicted, and the second fetch may evict the first
                if (symmetric) {
                    bra_owned = fetch_batch(world, cloud,
                                            exchange_batch_record_key(salt, EXCHANGE_BATCH_VF, bra_range));
                    vf_owned = diagonal_block
                            ? bra_owned                     // one range, so one batch
                            : fetch_batch(world, cloud,
                                          exchange_batch_record_key(salt, EXCHANGE_BATCH_VF, vf_range));
                } else {
                    // The column batch is held: this task runs on the rank owning it, so the fetch
                    // is a cache hit and those coefficients never move. The row batch is the one
                    // that rotates, and it carries bra and ket over the same range.
                    vf_owned = fetch_batch(world, cloud,
                                           exchange_batch_record_key(salt, EXCHANGE_BATCH_VF, vf_range));
                    bra_owned = fetch_batch(world, cloud,
                                            exchange_batch_record_key(salt, EXCHANGE_BATCH_BRA, bra_range));
                    ket_owned = fetch_batch(world, cloud,
                                            exchange_batch_record_key(salt, EXCHANGE_BATCH_KET, bra_range));
                    ket_work = &ket_owned;
                }
                bra_work = &bra_owned;
                vf_work = &vf_owned;
            }
            // everything up to here was getting the operands in hand; everything after is compute
            const double compute_wall_start = wall_time();
            const double compute_cpu_start = process_cpu_time();
            if (profiling) prof_.wait_for_data_wall = compute_wall_start - tile_wall_start;

            if (symmetric and diagonal_block) {
                const vecfuncT ket_batch = owner_pinned ? *bra_work : bra_range.copy_batch(vket);
                vecfuncT resultcolumn = compute_diagonal_batch_in_symmetric_matrix(world, ket_batch, *bra_work,
                                                                                   *vf_work);

                for (int i = vf_range.begin; i < vf_range.end; ++i){
                    Kf[i] += resultcolumn[i - vf_range.begin];}

            } else if (symmetric and not diagonal_block) {
                const vecfuncT ket_rows = owner_pinned ? *bra_work : bra_range.copy_batch(vket);
                const vecfuncT ket_columns = owner_pinned ? *vf_work : vf_range.copy_batch(vket);
                auto[resultcolumn, resultrow]=compute_offdiagonal_batch_in_symmetric_matrix(world, ket_rows,
                                                                                            ket_columns,
                                                                                            *bra_work, *vf_work);

                for (int i = bra_range.begin; i < bra_range.end; ++i){
                    Kf[i] += resultcolumn[i - bra_range.begin];}
                for (int i = vf_range.begin; i < vf_range.end; ++i){
                    Kf[i] += resultrow[i - vf_range.begin];}
            } else {
                // the ket comes from its own record when pinned, and is sliced from the unbatched
                // argument otherwise -- which is the only form available without the batch store
                const vecfuncT ket_batch = ket_work ? *ket_work : bra_range.copy_batch(vket);
                vecfuncT resultcolumn = compute_batch_in_asymmetric_matrix(world, ket_batch,
                                                                           *bra_work, *vf_work);
                for (int i = vf_range.begin; i < vf_range.end; ++i)
                    Kf[i] += resultcolumn[i - vf_range.begin];
            }

            // this tile's cost, for the next call's placement. Keyed by its batch pair, so the
            // reference survives a different partition only if the pair count matches -- which
            // prepare_owner_assignment checks before using it.
            if (profiling) {
                prof_.compute_wall = wall_time() - compute_wall_start;
                prof_.compute_cpu = process_cpu_time() - compute_cpu_start;
                prof_.wall_end = wall_time();
                prof_.peak_rss_gb = get_rss_usage_in_GB();
                // one record per task, from the rank that ran it
                if (world.rank() == 0) exch_write_task_profile(prof_);
            }

            if (owner_pinned and not cost_this_call_.empty()) {
                const auto ic = batch_begin_to_index_.find(vf_range.begin);
                const auto jr = batch_begin_to_index_.find(bra_range.begin);
                if (ic != batch_begin_to_index_.end() and jr != batch_begin_to_index_.end()) {
                    const long t = exchange_sym_tri(ic->second, jr->second);
                    if (t < long(cost_this_call_.size()))
                        cost_this_call_[t] += wall_time() - tile_wall_start;
                }
            }
            return Kf;
        }

        /// compute a batch of the exchange matrix, with identical ranges, exploiting the matrix symmetry

        /// \param subworld     the world we're computing in
        /// \param cloud        where to store the results
        /// \param bra_batch    the bra batch of orbitals (including the nuclear correlation factor square)
        /// \param ket_batch    the ket batch of orbitals, i.e. the orbitals to premultiply with
        /// \param vf_batch     the argument of the exchange operator
        /// Streams the tile one row at a time, so only the intermediates of a single row
        /// are live at once where computing the tile in one go holds the whole triangle.
        /// Row i builds N_ij = P(bra[i] vf[j]) for j <= i, adds ket[i] N_ij to column j,
        /// and adds the mirrored ket[j] N_ij to column i.
        vecfuncT compute_diagonal_batch_in_symmetric_matrix(World& subworld,
                                                            const vecfuncT& ket_batch,      // is batched
                                                            const vecfuncT& bra_batch,      // is batched
                                                            const vecfuncT& vf_batch        // is batched
        ) const {
            MADNESS_CHECK_THROW(ket_batch.size() == bra_batch.size(),
                                "symmetric diagonal tile: ket/bra batch size mismatch");
            MADNESS_CHECK_THROW(vf_batch.size() == bra_batch.size(),
                                "symmetric diagonal tile: vf/bra batch size mismatch");

            const long n = long(vf_batch.size());
            vecfuncT resultcolumn = zero_functions_compressed<T, NDIM>(subworld, n);
            auto poisson = Exchange<double, 3>::ExchangeImpl::set_poisson(subworld, lo);

            // per-stage wall, for the profiler only; `tick` is a no-op when it is off
            const bool prof_on = profile_active();
            auto tick = [&](double& acc, const double t0) { if (prof_on) acc += wall_time() - t0; };

            for (long i = 0; i < n; ++i) {
                double cpu0 = cpu_time();
                double w0 = prof_on ? wall_time() : 0.0;
                const vecfuncT vf_subset(vf_batch.begin(), vf_batch.begin() + i + 1);
                vecfuncT psif = mul_sparse(subworld, bra_batch[i], vf_subset, mul_tol);
                tick(prof_.mul1_wall, w0);
                w0 = prof_on ? wall_time() : 0.0;
                truncate(subworld, psif);
                tick(prof_.truncate_wall, w0);
                double cpu1 = cpu_time();
                mul1_timer += long((cpu1 - cpu0) * 1000l);

                cpu0 = cpu_time();
                w0 = prof_on ? wall_time() : 0.0;
                psif = apply(subworld, *poisson.get(), psif);
                tick(prof_.apply_wall, w0);
                w0 = prof_on ? wall_time() : 0.0;
                truncate(subworld, psif);
                tick(prof_.truncate_wall, w0);
                cpu1 = cpu_time();
                apply_timer += long((cpu1 - cpu0) * 1000l);

                // rows overlap in the columns they touch, so this row's contribution is
                // assembled on its own and accumulated
                cpu0 = cpu_time();
                vecfuncT update = zero_functions_compressed<T, NDIM>(subworld, n);
                double w1 = prof_on ? wall_time() : 0.0;
                vecfuncT row_contrib = mul_sparse(subworld, ket_batch[i], psif, mul_tol);
                tick(prof_.mul2_wall, w1);
                compress(subworld, row_contrib);
                for (long j = 0; j <= i; ++j) update[j] += row_contrib[j];
                for (long j = 0; j < i; ++j) {
                    w1 = prof_on ? wall_time() : 0.0;
                    vecfuncT mirrored = mul_sparse(subworld, ket_batch[j], vecfuncT(1, psif[j]), mul_tol);
                    tick(prof_.mul2_wall, w1);
                    compress(subworld, mirrored);
                    update[i] += mirrored[0];
                }
                gaxpy(subworld, 1.0, resultcolumn, 1.0, update);
                cpu1 = cpu_time();
                mul2_timer += long((cpu1 - cpu0) * 1000l);
            }
            // !! NO TRUNCATION AT THIS POINT !!
            return resultcolumn;
        }

        /// compute a batch of the exchange matrix, with non-identical ranges

        /// \param subworld     the world we're computing in
        /// \param cloud        where to store the results
        /// \param bra_batch    the bra batch of orbitals (including the nuclear correlation factor square)
        /// \param ket_batch    the ket batch of orbitals, i.e. the orbitals to premultiply with
        /// \param vf_batch     the argument of the exchange operator
        vecfuncT compute_batch_in_asymmetric_matrix(World& subworld,
                                                    const vecfuncT& ket_batch,
                                                    const vecfuncT& bra_batch,
                                                    const vecfuncT& vf_batch) const {
            double symmetric = false;
            auto poisson = Exchange<double, 3>::ExchangeImpl::set_poisson(subworld, lo);
            return Exchange<T, NDIM>::ExchangeImpl::compute_K_tile(subworld, bra_batch, ket_batch, vf_batch, poisson, symmetric,
                                                     mul_tol);
        }

        /// compute a batch of the exchange matrix, with non-identical ranges

        /// The caller supplies the ket over each of the tile's two ranges: it is the one
        /// that knows the ranges, and where those orbitals come from depends on how the
        /// operands were supplied, which is not the kernel's concern.
        std::pair<vecfuncT, vecfuncT> compute_offdiagonal_batch_in_symmetric_matrix(World& subworld,
                                                                                    const vecfuncT& ket_rows,    // ket over the bra/row range
                                                                                    const vecfuncT& ket_columns, // ket over the vf/column range
                                                                                    const vecfuncT& bra_batch,   // batched
                                                                                    const vecfuncT& vf_batch) const; // batched

    };

    class MacroTaskExchangeRow : public MacroTaskOperationBase {

        long nresult;
        double lo = 1.e-4;
        double mul_tol = 1.e-7;
        bool symmetric = false;
        Algorithm algorithm_;

        /// custom partitioning for the exchange operator in exchangeoperator.h
        class MacroTaskPartitionerRow : public MacroTaskPartitioner {
        public:
            MacroTaskPartitionerRow() {
              max_batch_size=1;
            }       
        };

    public:
        MacroTaskExchangeRow(const long nresult, const double lo, const double mul_tol, const Algorithm algorithm)
                : nresult(nresult), lo(lo), mul_tol(mul_tol),  algorithm_(algorithm) {
            partitioner.reset(new MacroTaskPartitionerRow());
            name="MacroTaskExchangeRow";
        }

        // you need to define the exact argument(s) of operator() as tuple
        typedef std::tuple<const std::vector<Function<T, NDIM>>&,
                           const std::vector<Function<T, NDIM>>&,
                           const std::vector<Function<T, NDIM>>&> argtupleT;

        using resultT = std::vector<Function<T, NDIM>>;

        // you need to define an empty constructor for the result
        // resultT must implement operator+=(const resultT&)
        resultT allocator(World& world, const argtupleT& argtuple) const {
            std::size_t n = std::get<0>(argtuple).size();
            resultT result = zero_functions_compressed<T, NDIM>(world, n);
            return result;
        }

        /// compute exchange row-wise for a fixed orbital phi_i of vket

        /// create 2 worlds: one fetches the function coefficients from the universe, the other
        /// does the computation, then swap. The result is copied back to the universe
        std::vector<Function<T, NDIM>>
        operator()(const std::vector<Function<T, NDIM>>& vket,
                   const std::vector<Function<T, NDIM>>& mo_bra, 
                   const std::vector<Function<T, NDIM>>& mo_ket) {
            std::vector<Function<T,NDIM>> result;
            if (algorithm_==fetch_compute) {
                result=row_fetch_compute(vket,mo_bra,mo_ket);
            } else if (algorithm_==multiworld_efficient_row) {
                result=row(vket,mo_bra,mo_ket);
            } else {
                MADNESS_EXCEPTION("unknown algorithm in Exchange::MacroTaskExchangeRow::operator()",1);
            }
            return result;
        }

        std::vector<Function<T,NDIM>>
        row(const std::vector<Function<T, NDIM>>& vket,
            const std::vector<Function<T, NDIM>>& mo_bra,
            const std::vector<Function<T, NDIM>>& mo_ket) {

            double cpu0, cpu1;
            World& world = vket.front().world();

            resultT Kf = zero_functions_compressed<T, NDIM>(world, 1);
            vecfuncT psif = zero_functions_compressed<T,NDIM>(world, mo_bra.size());
            auto poisson = Exchange<double, 3>::ExchangeImpl::set_poisson(world, lo);

            // !! NO !! vket is batched, starts at batch.input[0].begin
            // auto& i = batch.input[0].begin;
            long i=0;
            MADNESS_CHECK_THROW(vket.size()==1,"out-of-bounds error in Exchange::MacroTaskExchangeRow::operator()");
            size_t min_tile = 10;
            size_t ntile = std::min(mo_bra.size(), min_tile);

            for (size_t ilo=0; ilo<mo_bra.size(); ilo+=ntile){
                cpu0 = cpu_time();
                size_t iend = std::min(ilo+ntile,mo_bra.size());
                vecfuncT tmp_mo_bra(mo_bra.begin()+ilo,mo_bra.begin()+iend);
                auto tmp_psif = mul_sparse(world, vket[i], tmp_mo_bra, mul_tol);
                truncate(world, tmp_psif);
                cpu1 = cpu_time();
                mul1_timer += long((cpu1 - cpu0) * 1000l);

                cpu0 = cpu_time();
                tmp_psif = apply(world, *poisson.get(), tmp_psif);
                truncate(world, tmp_psif);
                cpu1 = cpu_time();
                apply_timer += long((cpu1 - cpu0) * 1000l);

                cpu0 = cpu_time();
                vecfuncT tmp_mo_ket(mo_ket.begin()+ilo,mo_ket.begin()+iend);
                // screen the second multiplication too, at the same tolerance as the first
                auto tmp_Kf = dot(world, tmp_mo_ket, tmp_psif, true, true, mul_tol);
                cpu1 = cpu_time();
                mul2_timer += long((cpu1 - cpu0) * 1000l);

                Kf[0] += tmp_Kf;
                truncate(world, Kf);
            }

            return Kf;
        }

        std::vector<Function<T,NDIM>>
        row_fetch_compute(const std::vector<Function<T, NDIM>>& vket,
            const std::vector<Function<T, NDIM>>& mo_bra,
            const std::vector<Function<T, NDIM>>& mo_ket) {

            io_redirect_cout();
            double total_execution_time=0.0;
            double total_fetch_time=0.0;
            double total_fetch_spawn_time=0.0;

            resultT Kf = zero_functions_compressed<T, NDIM>(*subworld_ptr, 1);
            {
                // create the two worlds that will be used for fetching and computing
                // std::shared_ptr<World> executing_world(subworld_ptr);
                double cpu0=cpu_time();
                SafeMPI::Intracomm comm = subworld_ptr->mpi.comm();
                std::shared_ptr<World> fetching_world(new World(comm.Clone()));
                std::shared_ptr<World> executing_world(new World(comm.Clone()));
                double cpu1=cpu_time();
                print("time to create two worlds:",cpu1-cpu0,"seconds");
                print("executing_world.id()",executing_world->id(),"fetching_world.id()",fetching_world->id(),"in MacroTaskExchangeRow");

                {
                    auto poisson1 = Exchange<double, 3>::ExchangeImpl::set_poisson(*executing_world, lo);
                    auto poisson2 = Exchange<double, 3>::ExchangeImpl::set_poisson(*fetching_world, lo);

                    functionT phi1=copy(*executing_world,vket[0]);
                    functionT phi2=copy(*fetching_world,vket[0]);

                    // !! NO !! vket is batched, starts at batch.input[0].begin
                    // auto& i = batch.input[0].begin;
                    MADNESS_CHECK_THROW(vket.size()==1,"out-of-bounds error in Exchange::MacroTaskExchangeRow::operator()");
                    size_t min_tile = 10;
                    size_t ntile = std::min(mo_bra.size(), min_tile);

                    struct Tile {
                        size_t ilo;
                        size_t iend;
                    };


                    // copy the data from the universe bra and ket to subworld bra and ket
                    // returns a pair of vectors in the subworld which are still awaiting the function coefficients
                    auto fetch_data = [&](World& world, const Tile& tile) {
                        MADNESS_CHECK_THROW(mo_bra.size()==mo_ket.size(),
                                            "bra and ket size mismatch in Exchange::MacroTaskExchangeRow::execute()");

                        std::size_t sz=tile.iend-tile.ilo;
                        vecfuncT subworld_bra(sz);
                        vecfuncT subworld_ket;
                        for (size_t i=tile.ilo; i<tile.iend; ++i) {
                            auto f=copy(world,mo_bra[i],false);
                            subworld_bra[i-tile.ilo]=f;
                            subworld_ket.push_back(copy(world, mo_ket[i],false));
                        }
                        return std::make_pair(subworld_bra,subworld_ket);
                    };

                    // apply the exchange operator on phi for a a tile of mo_bra and mo_ket
                    auto execute = [&](World& world, auto poisson, const functionT& phi, const vecfuncT& mo_bra, const vecfuncT& mo_ket) {
                        MADNESS_CHECK_THROW(mo_bra.size()==mo_ket.size(),
                                            "bra and ket size mismatch in Exchange::MacroTaskExchangeRow::execute()");

                        auto world_id=world.id();
                        auto phi_id=phi.world().id();
                        auto bra_id=mo_bra.front().world().id();
                        auto ket_id=mo_ket.front().world().id();
                        std::string msg="world mismatch in Exchange::MacroTaskExchangeRow::execute(): ";
                        msg+="world.id()="+std::to_string(world_id)+", ";
                        msg+="phi.world().id()="+std::to_string(phi_id)+", ";
                        msg+="bra.world().id()="+std::to_string(bra_id)+", ";
                        msg+="ket.world().id()="+std::to_string(ket_id);
                        if (not (world_id==phi_id && world_id==bra_id && world_id==ket_id)) {
                            print(msg);
                        }
                        MADNESS_CHECK_THROW(world_id==phi_id && world_id==bra_id && world_id==ket_id,msg.c_str());

                        double cpu0 = cpu_time();
                        auto tmp_psif = mul_sparse(world, phi, mo_bra, mul_tol);
                        truncate(world, tmp_psif);
                        double cpu1 = cpu_time();
                        mul1_timer += long((cpu1 - cpu0) * 1000l);

                        cpu0 = cpu_time();
                        tmp_psif = apply(world, *poisson.get(), tmp_psif);
                        truncate(world, tmp_psif);
                        cpu1 = cpu_time();
                        apply_timer += long((cpu1 - cpu0) * 1000l);

                        cpu0 = cpu_time();
                        auto tmp_Kf = dot(world, mo_ket, tmp_psif);
                        cpu1 = cpu_time();
                        mul2_timer += long((cpu1 - cpu0) * 1000l);

                        return tmp_Kf.truncate();

                    };

                    std::vector<Tile> tiles;
                    for (size_t ilo=0; ilo<mo_bra.size(); ilo+=ntile) {
                        tiles.push_back(Tile{ilo,std::min(ilo+ntile,mo_bra.size())});
                    }

                    vecfuncT tmp_mo_bra1,tmp_mo_ket1;
                    vecfuncT tmp_mo_bra2,tmp_mo_ket2;

                    for (size_t itile=0; itile<tiles.size(); ++itile) {
                        Tile& tile = tiles[itile];

                        if (itile==0) {
                            double t0=cpu_time();
                            print("fetching tile",tile.ilo,"into world",executing_world->id());
                            std::tie(tmp_mo_bra1,tmp_mo_ket1)=fetch_data(*executing_world,tiles[itile]);
                            fetching_world->gop.set_forbid_fence(false);
                            double t2=cpu_time();
                            executing_world->gop.fence();
                            double t1=cpu_time();
                            total_fetch_time += (t1 - t0);
                            total_fetch_spawn_time += (t2 - t0);
                        }

                        double t0=cpu_time();
                        fetching_world->gop.set_forbid_fence(true);
                        if (itile<tiles.size()-1) {
                            // fetch data into fetching_world while computing in executing_world
                            print("fetching tile",tiles[itile+1].ilo,"into world",fetching_world->id()," at time ",wall_time());
                            std::tie(tmp_mo_bra2,tmp_mo_ket2)=fetch_data(*fetching_world,tiles[itile+1]);
                        }
                        fetching_world->gop.set_forbid_fence(false);
                        double t2=cpu_time();
                        // uncomment the next line to enforce that fetching is finished before executing
                        // fetching_world->gop.fence();
                        double t1=cpu_time();
                        total_fetch_time += (t1 - t0);
                        total_fetch_spawn_time += (t2 - t0);

                        print("executing tile",tile.ilo,"in world",executing_world->id());
                        double dpu0=cpu_time();
                        Kf[0]+=execute(*executing_world,poisson1,phi1,tmp_mo_bra1,tmp_mo_ket1);
                        double dpu1=cpu_time();
                        print("time to execute tile",tile.ilo,"in world",executing_world->id(),dpu1-dpu0,"seconds");
                        total_execution_time += dpu1-dpu0;

                        fetching_world->gop.fence();

                        // change roles of the two worlds
                        std::swap(poisson1,poisson2);
                        std::swap(phi1,phi2);
                        std::swap(tmp_mo_bra2,tmp_mo_bra1);
                        std::swap(tmp_mo_ket2,tmp_mo_ket1);
                        std::swap(executing_world,fetching_world);
                    }
                } // objects living in the two worlds must be destroyed before the worlds are freed

                // deferred destruction of WorldObjects happens here
                fetching_world->gop.fence();
                executing_world->gop.fence();
                double cpu2=cpu_time();
                print("overall time: ",cpu2-cpu0,"seconds");
                print("total execution time:",total_execution_time,"seconds");
                print("total fetch time:",total_fetch_time,"seconds");
                print("total fetch spawn time:",total_fetch_spawn_time,"seconds");
            } // worlds are destroyed here

            return Kf;
        }
    };
};

} /* namespace madness */

#endif /* SRC_APPS_CHEM_EXCHANGEOPERATOR_H_ */
