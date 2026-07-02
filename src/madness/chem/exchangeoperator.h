#ifndef SRC_APPS_CHEM_EXCHANGEOPERATOR_H_
#define SRC_APPS_CHEM_EXCHANGEOPERATOR_H_

#include<madness.h>
#include<madness/world/cloud.h>
#include<madness/mra/macrotaskq.h>
#include<madness/chem/SCFOperators.h>
#include<madness/world/ranks_and_hosts.h>
#include<unordered_map>
#include<queue>
#include<random>
#include<numeric>
#include<set>
#include<limits>
#include<fstream>

namespace madness {

// forward declaration
class SCF;
class Nemo;

/// Per-task performance record for the owner-pinned exchange (smallmem_sym_p2p_owner
/// primarily). Emitted as one JSON object per task to a per-subworld .jsonl file when
/// the env var MAD_EXCH_TASK_PROFILE is set; merge the files in post. Designed for a
/// Gantt view: [wall_start, wall_end] places the task; the interior splits into
/// wait-for-data, pure compute (further broken into mul1/apply/mul2/truncate), and
/// residual overhead. Wall times are honest per-component because the smallmem kernels
/// run every op with fence=true, so each op completes before the next (no added fence).
///   transient_class: worst fetch tier this task hit (lru < prefetch < sync); the owned
///   operand is always local, so the worst tier is the transient one — sync => the task
///   stalled waiting for remote bytes (waited=true), the key efficiency signal.
struct ExchTaskProfile {
    long task_id = -1;
    long universe_rank = 0;                     ///< keys the output file (one per process)
    unsigned long subworld_id = 0;
    int subworld_nrank = 0;
    double thresh = 0.0;                        ///< FunctionDefaults thresh at task time (protocol tier)
    long k = 0;                                 ///< FunctionDefaults wavelet order at task time
    bool diagonal = false;
    long row_begin = 0, row_end = 0, col_begin = 0, col_end = 0;
    long n_row = 0, n_col = 0;
    double wall_start = 0.0, wall_end = 0.0;   ///< absolute wall_time() clock
    double wait_for_data_wall = 0.0;            ///< bra+vf fetch wall (task-ready -> first compute)
    int transient_class = -1;                   ///< 0=lru,1=prefetch,2=sync,-1=n/a
    bool waited = false;                        ///< a synchronous (remote) fetch happened
    double compute_wall = 0.0, compute_cpu = 0.0;
    double mul1_wall = 0.0, apply_wall = 0.0, mul2_wall = 0.0, truncate_wall = 0.0;
    double peak_rss_gb = 0.0;                   ///< getrusage high-water-mark at task end

    void reset() { *this = ExchTaskProfile(); }
    void observe_fetch_tier(int tier) {         ///< tier: 0 lru, 1 prefetch, 2 sync
        if (tier > transient_class) transient_class = tier;
        if (tier == 2) waited = true;
    }
};

/// Is per-task exchange profiling enabled? Gated once on the env var MAD_EXCH_TASK_PROFILE.
inline bool exch_task_profile_enabled() {
    static const bool on = (std::getenv("MAD_EXCH_TASK_PROFILE") != nullptr);
    return on;
}

/// Append one task record as a JSON line to exch_taskprof.r<universe_rank>.jsonl.
/// Called only from subworld rank 0 when profiling is enabled. The output stream is
/// kept OPEN for the life of the process (one file per rank, appended across all
/// K-applications) — this bounds the file count at nproc regardless of iteration
/// count, and avoids a per-task open/close (costly metadata ops on a parallel FS).
/// subworld_id stays in each record so tasks can still be grouped by subworld/K-app.
inline void exch_write_task_profile(const ExchTaskProfile& p) {
    static std::ofstream os;   // one per process (this process has a single universe_rank)
    if (!os.is_open()) {
        char fname[64];
        std::snprintf(fname, sizeof(fname), "exch_taskprof.r%ld.jsonl", p.universe_rank);
        os.open(fname, std::ios::app);
        if (!os.is_open()) return;
    }
    auto tier_name = [](int t) {
        switch (t) { case 0: return "lru"; case 1: return "prefetch"; case 2: return "sync"; default: return "none"; }
    };
    const double residual = (p.wall_end - p.wall_start) - p.wait_for_data_wall - p.compute_wall;
    const double wall_total = p.wall_end - p.wall_start;
    const double compute_pct = (wall_total > 0.0) ? 100.0 * p.compute_wall / wall_total : 0.0;
    os.setf(std::ios::fixed); os.precision(6);
    os << "{"
       << "\"task_id\":" << p.task_id
       << ",\"universe_rank\":" << p.universe_rank
       << ",\"subworld_id\":" << p.subworld_id
       << ",\"subworld_nrank\":" << p.subworld_nrank
       << ",\"thresh\":" << p.thresh << ",\"k\":" << p.k
       << ",\"diagonal\":" << (p.diagonal ? "true" : "false")
       << ",\"row_range\":[" << p.row_begin << "," << p.row_end << "]"
       << ",\"col_range\":[" << p.col_begin << "," << p.col_end << "]"
       << ",\"n_row\":" << p.n_row << ",\"n_col\":" << p.n_col
       << ",\"wall_start\":" << p.wall_start << ",\"wall_end\":" << p.wall_end
       << ",\"wall_total\":" << wall_total
       << ",\"wait_for_data_wall\":" << p.wait_for_data_wall
       << ",\"transient_class\":\"" << tier_name(p.transient_class) << "\""
       << ",\"waited\":" << (p.waited ? "true" : "false")
       << ",\"compute_wall\":" << p.compute_wall
       << ",\"compute_cpu\":" << p.compute_cpu
       << ",\"accumulate_residual_wall\":" << residual
       << ",\"compute_components_wall\":{"
       <<   "\"mul1\":" << p.mul1_wall << ",\"apply\":" << p.apply_wall
       <<   ",\"mul2\":" << p.mul2_wall << ",\"truncate\":" << p.truncate_wall
       <<   ",\"other\":" << (p.compute_wall - p.mul1_wall - p.apply_wall - p.mul2_wall - p.truncate_wall)
       << "}"
       << ",\"compute_pct\":" << compute_pct
       << ",\"peak_rss_gb\":" << p.peak_rss_gb
       << "}\n";
    os.flush();   // keep the file fresh for mid-run inspection (no per-task open/close)
}

/// Split a vector of length n into exactly nbatch = min(nsubworld, n)
/// contiguous batches with sizes differing by at most 1. The first
/// (n mod nbatch) batches get size ceil(n/nbatch); the rest get
/// floor(n/nbatch). Eliminates the runt-of-1 case.
///
/// Single source of truth for the row-owner batch boundaries: used both by
/// MacroTaskPartitionerExchange (to build the task grid) and by the cloud
/// batch-storage side, so stored batches align with the partition by
/// construction. Depends only on (n, nsubworld).
inline std::vector<Batch_1D> exchange_row_owner_split(const std::size_t n, const long nsubworld) {
    std::vector<Batch_1D> out;
    if (n == 0) return out;
    const long nbatch = std::min<long>(std::max<long>(1, nsubworld), long(n));
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

/// Number of batches M for the symmetric round-robin (p2p) owner algorithm.
///
/// M = level * nsubworld: `level` (>= 1) is directly the number of batches per
/// rank, so EVERY level is selectable (level 1 -> 1 batch/rank, level 2 ->
/// 2/rank, ...) regardless of nsubworld parity. (The old formula carried a
/// parity device — base = 2 for even nsubworld — that made 1 batch/rank
/// UNREACHABLE for even nsubworld; removed so the granularity sweep can test
/// whether the even-tasks-per-rank constraint matters at all. The default knob
/// value is 2, preserving the previous even-nsubworld default of 2*nsubworld.)
/// The result is clamped to [1, n] (cannot have more batches than orbitals);
/// when n is below the formula value the caller is in the small-problem regime
/// and should fall back to small_memory_symmetric_mt_owner.
///
/// Each batch k is owned by rank (k mod nsubworld), so with the unclamped
/// formula every rank owns exactly `level` batches.
inline long exchange_sym_owner_nbatch(const std::size_t n, const long nsubworld,
                                      const long granularity_level) {
    if (n == 0) return 0;
    const long nsw = std::max<long>(1, nsubworld);
    const long level = std::max<long>(1, granularity_level);
    const long M = level * nsw;
    return std::min<long>(M, long(n));
}

/// Split a vector of length n into M = exchange_sym_owner_nbatch(...) contiguous
/// batches with sizes differing by at most 1 (same even-remainder spread as
/// exchange_row_owner_split). Single source of truth for the symmetric p2p-owner
/// batch boundaries: shared by the triangular partitioner, the cloud
/// batch-storage side, and the round-robin owner distribution, so the stored
/// batches, the task grid, and the ownership all align by construction.
/// Batch k is owned by rank (k mod nsubworld). Depends only on
/// (n, nsubworld, granularity_level).
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

/// Round-robin owner assignment for the lower-triangular batch-pair task matrix
/// of the symmetric p2p-owner algorithm. Pure function of (n, M):
///   - n = number of workers (= nsubworld = nproc); batch k owned by (k mod n).
///   - M = number of batches (= exchange_sym_owner_split(...).size()).
/// Task (i,j) with 0 <= j <= i < M is eligible for worker t iff
///   i==j: (i mod n)==t   (diagonal pinned to its single owner)
///   i!=j: (i mod n)==t or (j mod n)==t   (owns at least one operand)
/// Phase 1 round-robins eligible tasks across workers; Phase 2 rebalances
/// off-diagonal tasks (diagonals never move) until max-min load <= 1.
/// Returns a map (i,j) -> owner worker. Deterministic, identical on every rank.
inline std::map<std::pair<long,long>,long>
exchange_sym_round_robin_assign(const long n, const long M) {
    std::map<std::pair<long,long>,long> owner;
    if (n <= 0 or M <= 0) return owner;
    auto eligible = [n](long i, long j, long t) -> bool {
        if (i == j) return (i % n) == t;
        return ((i % n) == t) or ((j % n) == t);
    };
    // all tasks, ascending by (i,j)
    std::vector<std::pair<long,long>> remaining;
    remaining.reserve(std::size_t(M) * std::size_t(M + 1) / 2);
    for (long i = 0; i < M; ++i)
        for (long j = 0; j <= i; ++j)
            remaining.emplace_back(i, j);

    std::vector<std::vector<std::pair<long,long>>> T(n);

    // Phase 1 — round-robin placement
    long t = 0, misses = 0;
    while (not remaining.empty() and misses < n) {
        long picked = -1;
        for (long idx = 0; idx < long(remaining.size()); ++idx) {
            if (eligible(remaining[idx].first, remaining[idx].second, t)) { picked = idx; break; }
        }
        if (picked >= 0) {
            T[t].push_back(remaining[picked]);
            remaining.erase(remaining.begin() + picked);
            misses = 0;
        } else {
            ++misses;
        }
        t = (t + 1) % n;
    }
    // any leftovers (should be none — every task is eligible for someone): assign by i%n
    for (const auto& task : remaining) T[task.first % n].push_back(task);

    // Phase 2 — rebalance until max-min <= 1 (cap iterations defensively)
    for (long iter = 0; iter < 10000; ++iter) {
        long big = 0, small = 0;
        for (long p = 1; p < n; ++p) {
            if (long(T[p].size()) > long(T[big].size())) big = p;
            if (long(T[p].size()) < long(T[small].size())) small = p;
        }
        if (long(T[big].size()) - long(T[small].size()) <= 1) break;
        // move the first non-diagonal task in T[big] eligible for `small`
        bool moved = false;
        // re-derive sorted big/small pairs greedily: scan all (big,small) by load
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
        if (not moved) break;   // accept residual imbalance
    }

    for (long p = 0; p < n; ++p)
        for (const auto& tk : T[p])
            owner[tk] = p;
    return owner;
}


/// ---- Option C0: coalesced finalize transfer (exchange-local) ----
/// One coefficient node in transit during the coalesced finalize drain. Carries
/// the destination function index, the node key, and the full FunctionNode
/// (coeff + has_children + norms) — identical wire content to the per-node active
/// message in do_gaxpy_inplace, just batched. Serialized automatically as a
/// WorldObject::task argument.
template <typename T, std::size_t NDIM>
struct FinalizeNodeRec {
    std::size_t f = 0;                  ///< index into the destination function vector
    Key<NDIM> key;                      ///< tree-node key
    FunctionNode<T, NDIM> node;         ///< source node (full FunctionNode)
    template <typename Archive>
    void serialize(Archive& ar) { ar & f & key & node; }
};

/// Receiver endpoint for the coalesced finalize transfer. A WorldObject living in
/// the transport world (the universe for the subworld/node -> universe drains, the
/// node sub-World for stage1). Source nodes are pushed here in bulk chunks via
/// task() — so the deserialize + accumulate runs on a worker, never the comm
/// thread — and accumulated into the local slices of the destination functions,
/// reproducing FunctionNode::gaxpy_inplace exactly. Cached per transport world;
/// set_targets() refreshes the destination impls + beta before each finalize,
/// under the readiness barrier in coalesced_gaxpy.
template <typename T, std::size_t NDIM>
class FinalizeReducer : public WorldObject<FinalizeReducer<T, NDIM>> {
public:
    typedef FunctionImpl<T, NDIM> implT;
    typedef FinalizeNodeRec<T, NDIM> recT;

    explicit FinalizeReducer(World& world)
        : WorldObject<FinalizeReducer<T, NDIM>>(world) {
        this->process_pending();
    }

    /// (re)point at the destination function impls for the next drain. Local;
    /// every rank calls this before coalesced_gaxpy's readiness barrier.
    void set_targets(std::vector<implT*> dests, T beta) {
        dests_ = std::move(dests);
        beta_ = beta;
    }

    /// worker-side: accumulate a chunk of source nodes into the local slices of
    /// the destination functions. The accessor write-lock serializes concurrent
    /// senders to the same key; the math is identical to do_gaxpy_inplace
    /// (FunctionNode::gaxpy_inplace, funcimpl.h).
    void accumulate_chunk(const std::vector<recT>& recs) {
        for (const auto& r : recs) {
            MADNESS_ASSERT(r.f < dests_.size() and dests_[r.f] != nullptr);
            typename implT::dcT::accessor acc;
            dests_[r.f]->get_coeffs().insert(acc, r.key);   // get-or-create, locks key
            acc->second.template gaxpy_inplace<T, T>(T(1.0), r.node, beta_);
        }
    }

private:
    std::vector<implT*> dests_;
    T beta_ = T(1.0);
};

/// Coalesced replacement for the per-node gaxpy scatter in the exchange finalize.
/// Instead of one active message per source tree node, bucket the local source
/// nodes by their destination owner and push each bucket as bulk chunks
/// (task() -> worker-side accumulate). `dest_vec` and `src_vec` may live in
/// different worlds (subworld/node -> universe, or subworld -> node); routing uses
/// the DESTINATION pmap and the transport rides `transport_world` (== dest's
/// world). Completion is the caller's existing fence (node fence / universe
/// fence); the only fence issued here is the readiness barrier that guarantees
/// every rank has set its reducer targets before any chunk task can arrive.
///
/// MUST be called collectively on transport_world (every rank, once per finalize);
/// ranks with an empty `src_vec` still participate in the barrier (no buckets).
template <typename T, std::size_t NDIM>
void coalesced_gaxpy(World& transport_world,
                     FinalizeReducer<T, NDIM>& reducer,
                     std::vector<Function<T, NDIM>>& dest_vec,
                     std::vector<Function<T, NDIM>>& src_vec,
                     const T beta,
                     const std::size_t chunk_entries) {
    typedef FunctionImpl<T, NDIM> implT;
    typedef FinalizeNodeRec<T, NDIM> recT;

    // point the reducer at the destination impls (local), then barrier so no
    // chunk task can run on a rank that has not refreshed its targets yet.
    std::vector<implT*> dests(dest_vec.size(), nullptr);
    for (std::size_t f = 0; f < dest_vec.size(); ++f)
        dests[f] = dest_vec[f].get_impl().get();
    reducer.set_targets(dests, beta);
    transport_world.gop.fence();

    // bucket local source nodes by destination owner; flush as bulk chunks.
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
    static inline std::atomic<long> mul2_truncate_timer;
    static inline std::atomic<long> mul2_timer;
    static inline std::atomic<long> mul1_truncate_timer;
    static inline std::atomic<long> mul1_timer; ///< timing
    static inline std::atomic<long> owner_fetch_timer;
    static inline std::atomic<long> owner_compute_timer;       ///< CPU time inside compute kernels (process_cpu_time)
    static inline std::atomic<long> owner_compute_wall_timer;  ///< wall time wrapping the compute kernel call in operator() — measures "mere compute" including any in-compute communication/fences; subtract owner_compute_timer to attribute wait/communication overhead
    static inline std::atomic<long> sym_cache_hits_;           ///< sym_p2p: bounded-LRU hits on cloud-fetched batches (cache-aware sort payoff)
    static inline std::atomic<long> sym_cache_misses_;         ///< sym_p2p: bounded-LRU misses (a synchronous p2p fetch happened)
    static inline std::atomic<long> sym_prefetch_hits_;        ///< sym_p2p: Design-B prefetch consumes (transient requested one task ahead, overlapped with compute)
    static inline double elapsed_time;
    static inline double elapsed_process_cpu_time;

    static void reset_timer() {
        mul1_timer = 0l;
        mul1_truncate_timer = 0l;
        mul2_timer = 0l;
        mul2_truncate_timer = 0l;
        apply_timer = 0l;
        owner_fetch_timer = 0l;
        owner_compute_timer = 0l;
        owner_compute_wall_timer = 0l;
        sym_cache_hits_ = 0l;
        sym_cache_misses_ = 0l;
        sym_prefetch_hits_ = 0l;
        elapsed_time = 0.0;
        elapsed_process_cpu_time = 0.0;
    }

public:
    nlohmann::json gather_timings(World& world) const {
        double t1 = double(mul1_timer) * 0.001;
        double t2 = double(apply_timer) * 0.001;
        double t3 = double(mul2_timer) * 0.001;
        double t4 = double(mul1_truncate_timer) * 0.001;
        double t5 = double(mul2_truncate_timer) * 0.001;
        double t_fetch_owner = double(owner_fetch_timer) * 0.001;
        double t_compute_owner = double(owner_compute_timer) * 0.001;
        double t_compute_owner_wall = double(owner_compute_wall_timer) * 0.001;
        world.gop.sum(t1);
        world.gop.sum(t2);
        world.gop.sum(t3);
        world.gop.sum(t4);
        world.gop.sum(t5);
        world.gop.sum(t_fetch_owner);
        world.gop.sum(t_compute_owner);
        world.gop.sum(t_compute_owner_wall);
        nlohmann::json j;
        j["multiply1"] = t1;
        j["truncate1"] = t4;
        j["apply"] = t2;
        j["multiply2"] = t3;
        j["truncate2"] = t5;
        j["owner_fetch"] = t_fetch_owner;
        j["owner_compute"] = t_compute_owner;
        j["owner_compute_wall"] = t_compute_owner_wall;
        double sym_hits = double(sym_cache_hits_);
        double sym_misses = double(sym_cache_misses_);
        double sym_prefetch = double(sym_prefetch_hits_);
        world.gop.sum(sym_hits);
        world.gop.sum(sym_misses);
        world.gop.sum(sym_prefetch);
        j["sym_cache_hits"] = sym_hits;
        j["sym_cache_misses"] = sym_misses;
        j["sym_prefetch_hits"] = sym_prefetch;
        double total_cpu = elapsed_process_cpu_time;
        world.gop.sum(total_cpu);
        j["total_cpu"] = total_cpu;
        j["total"] = elapsed_time;
        return j;
    }

    void print_timer(World& world) const {
        auto timings= gather_timings(world);
        if (world.rank() == 0) {
            printf(" cpu time spent in multiply1   %8.2fs\n", timings["multiply1"].template get<double>());
            printf(" cpu time spent in truncate1   %8.2fs\n", timings["truncate1"].template get<double>());
            printf(" cpu time spent in apply       %8.2fs\n", timings["apply"].template get<double>());
            printf(" cpu time spent in multiply2   %8.2fs\n", timings["multiply2"].template get<double>());
            printf(" cpu time spent in truncate2   %8.2fs\n", timings["truncate2"].template get<double>());
            printf(" cpu time owner fetch          %8.2fs\n", timings["owner_fetch"].template get<double>());
            printf(" cpu time owner compute        %8.2fs\n", timings["owner_compute"].template get<double>());
            printf(" wall time owner compute       %8.2fs\n", timings["owner_compute_wall"].template get<double>());
            {
                const double h = timings["sym_cache_hits"].template get<double>();
                const double m = timings["sym_cache_misses"].template get<double>();
                const double pf = timings["sym_prefetch_hits"].template get<double>();
                if (h + m + pf > 0.0)
                    printf(" sym_p2p batch cache           %.0f LRU-hit / %.0f prefetch / %.0f sync (of %.0f lookups, %.1f%% non-sync)\n",
                           h, pf, m, h + m + pf, 100.0 * (h + pf) / (h + m + pf));
            }
            printf(" total process cpu time        %8.2fs\n", timings["total_cpu"].template get<double>());
            printf(" total wall time               %8.2fs\n", timings["total"].template get<double>());
        }
    }


    typedef Exchange<T,NDIM>::ExchangeAlgorithm Algorithm;
    Algorithm algorithm_ = multiworld_efficient_row;
    MacroTaskInfo macro_task_info = MacroTaskInfo::preset("default");
    bool replicate_for_debug_ = false;  ///< if true, use StoreFunction policy to pre-replicate all data (zero communication during tasks)
    int accumulation_mode_ = 2;         ///< owner-aware result-finalize mode: 0=per-task gaxpy into universe; 1=subworld-local accumulate + one bulk subworld->universe gaxpy; 2=node-local hierarchical reduction (default; intra-node reduce then one inter-node node-sum gaxpy, auto-degrades to mode 1 single-node)
    bool use_mflex_ = true;             ///< if true (and using owner-aware algorithm), run the m-flex peel search to load-balance the owner assignment
    long mflex_max_exhaustive_ = 5000;  ///< upper bound on C(R,m) for the exhaustive arm of the m-flex peel search
    bool use_cloud_batch_fetch_ = true; ///< if true (and algorithm==small_memory_mt_owner), fetch input batches from the cloud's owner-pinned batch container instead of copy()ing from the universe. Default off (A/B baseline).
    long batch_granularity_ = 2;        ///< granularity level (>=1) = batches per rank for small_memory_symmetric_p2p_owner: M = level*nproc batches. Default 2 (= 2 batches/rank, the previous default). Runtime tunable.
    bool bra_equals_ket_ = false;       ///< true iff set_bra_and_ket received the same vector for bra and ket (moldft / plain HF). Set in set_bra_and_ket. Gates the one-set small_memory_symmetric_p2p_owner path (nemo has bra=R^2*ket != ket).

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
        // Detect bra==ket (moldft passes the SAME vector for both; nemo passes
        // R^2*ket as bra). Compare impl handles BEFORE the deep copy. Cheap, no
        // collective. Gates the one-set small_memory_symmetric_p2p_owner path.
        // NOTE: this is a STRUCTURAL test (same impl handle), distinct from the
        // SEMANTIC value set in the ctors (real SCF -> true, complex/nemo -> false).
        // The moldft K(amo) path never calls this (bra/ket come from the ctor), so
        // the two rules don't collide there. But a caller passing semantically-equal
        // but structurally-distinct vectors (e.g. set_bra_and_ket(conj(amo), amo) for
        // real amo -> conj is a fresh-impl copy) would get false here and trip the
        // hard sym_p2p guard in K_macrotask_efficient. Fail-safe (rejects a valid
        // case, never silently wrong); revisit if a real symmetric caller hits it.
        bra_equals_ket_ = (bra.size() == ket.size());
        for (std::size_t i = 0; bra_equals_ket_ and i < bra.size(); ++i)
            bra_equals_ket_ = (bra[i].get_impl().get() == ket[i].get_impl().get());
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

    /// Enable fetching input batches from the cloud's owner-pinned batch container
    /// (only effective for small_memory_mt_owner). Default off keeps the existing
    /// copy()-from-universe path so the two backends can be A/B compared.
    ExchangeImpl& set_use_cloud_batch_fetch(const bool flag) {
        use_cloud_batch_fetch_ = flag;
        return *this;
    }

    /// granularity level (>=1) for small_memory_symmetric_p2p_owner; clamped to >=1.
    ExchangeImpl& set_batch_granularity(const long& level) {
        batch_granularity_ = std::max<long>(1, level);
        return *this;
    }

    ExchangeImpl& set_max_batch_size(const long& n) {
        max_batch_size_ = std::max<long>(1, n);
        return *this;
    }
    
    ExchangeImpl& set_min_batch_size(const long& n) {
        min_batch_size_ = std::max<long>(1, n);
        return *this;
    }

    ExchangeImpl& set_replicate_for_debug(const bool flag) {
        replicate_for_debug_ = flag;
        return *this;
    }

    ExchangeImpl& set_use_mflex(const bool flag) {
        use_mflex_ = flag;
        return *this;
    }

    ExchangeImpl& set_mflex_max_exhaustive(const long& n) {
        mflex_max_exhaustive_ = std::max<long>(0, n);
        return *this;
    }

    ExchangeImpl& set_accumulation_mode(const int mode) {
        accumulation_mode_ = mode;
        return *this;
    }

    /// TEMP debug knob. See declaration of smallmem_mul_tol_ above.
    ExchangeImpl& set_smallmem_mul_tol(const double tol) {
        smallmem_mul_tol_ = tol;
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

    /// computing the upper triangle and mirrored contributions for symmetric bra/ket/vket
    vecfuncT K_small_memory_symmetric(const vecfuncT& vket, const double mul_tol = 0.0) const;

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
    long min_batch_size_ = 5;
    long max_batch_size_ = 30;
    double mul_tol = FunctionDefaults<NDIM>::get_thresh()*0.1;
    /// TEMP debug knob: if >= 0, used directly as the screening tol at every mul_sparse/dot
    /// call site inside the smallmem_*_mt_owner kernels (no *0.1 factor applied).
    /// Negative => legacy behavior (mul_tol*0.1). Remove after the sensitivity test.
    double smallmem_mul_tol_ = -1.0;

    mutable nlohmann::json statistics;  ///< statistics of the Cloud (timings, memory)  and of the parameters of this run

    class MacroTaskExchangeSimple : public MacroTaskOperationBase {

        long nresult;
        double lo = 1.e-4;
        double mul_tol = 1.e-7;
        bool symmetric = false;
        Algorithm algorithm_ = multiworld_efficient;
    public:
        bool replicate_for_debug_ = false;
        // Result-finalize mode (active when use_owner_aware_fetch() is true).
        // 0: original per-task subworld->universe gaxpy (no buffering).
        // 1: accumulate each tile into the subworld-local buffer Kf_local_ and
        //    do a single subworld->universe gaxpy in finalize_stage2() (default).
        // 2: node-local hierarchical reduction -- stage1 reduces every node rank's
        //    Kf_local_ into a node-shared Kf_node_ (intra-node, cheap), stage2 does
        //    one inter-node gaxpy of the node-sums into the universe result, cutting
        //    simultaneous cross-node scatterers from nranks to nnodes.
        // Set from ExchangeImpl::accumulation_mode_ before submitting the macrotask.
        int accumulation_mode_ = 2;
        // When true, run the m-flex peel search before fold_and_assign so the
        // owner assignment minimizes per-rank load delta. Set from ExchangeImpl
        // before submitting the macrotask. Default true.
        bool use_mflex_ = true;
        // Cap on C(R, m) for the exhaustive arm of the m-flex peel search.
        long mflex_max_exhaustive_ = 5000;
        // TEMP debug knob; mirrors ExchangeImpl::smallmem_mul_tol_. Set from
        // ExchangeImpl before submitting the macrotask. Negative => legacy mul_tol*0.1.
        double smallmem_mul_tol_ = -1.0;
        // When true (and the algorithm is small_memory_mt_owner), fetch input
        // batches from the cloud's owner-pinned batch container (store_batch /
        // fetch_batch) instead of copy()ing from the universe functions. Set from
        // ExchangeImpl::use_cloud_batch_fetch_ before submitting. Default off so
        // the existing copy()-from-universe path stays the A/B baseline.
        bool use_cloud_batch_fetch_ = true;
        // Granularity level (>=1) = batches per rank for small_memory_symmetric_p2p_owner.
        // Drives M = level*nproc batches via exchange_sym_owner_split, shared
        // by the partitioner, store_batches, and the round-robin owner
        // assignment. Set from ExchangeImpl::batch_granularity_ before submitting.
        long batch_granularity_ = 2;
        // When true, dump owner_map_ and (for cloud-batch) the cloud's batch
        // record routing at universe rank 0, once per K-application, right
        // after the partition + owner assignment is built. Set from
        // ExchangeImpl::printdebug() in K_macrotask_efficient.
        bool log_diagnostics_ = false;
        // This process's universe rank, set from ExchangeImpl::world.rank() before
        // submitting. Used only to key the per-task profile file (MAD_EXCH_TASK_PROFILE)
        // to ONE file per process (appended across all K-applications), instead of one
        // file per subworld-id which would explode with iterations x nproc.
        long universe_rank_ = 0;
    private:
        // Per-task performance record (MAD_EXCH_TASK_PROFILE). NON-static: scoped to
        // the single operator() call on this task object; the kernels and
        // sym_batch_lru_fetch write into it via `this` during that call. Reset at
        // operator() entry, emitted at operator() exit on subworld rank 0.
        mutable ExchTaskProfile prof_;
        // Monotonic per-process task sequence for profiling identity (the task object
        // does not know its framework `element`); unique together with subworld_id.
        static inline std::atomic<long> prof_task_seq_{0};
        static inline std::unordered_map<long, functionT> bra_cache_;
        static inline std::unordered_map<long, functionT> ket_cache_;
        // Dedicated cache for the held i-batch (vf) in small_memory_mt_owner.
        // Keyed by global vf index. Separate from bra_cache_/ket_cache_ to avoid
        // conceptual collision with the symmetric algorithm, which already uses
        // those caches for the bra/ket dimensions.
        static inline std::unordered_map<long, functionT> held_vf_cache_;
        // Subworld-local accumulator for own-output accumulation.
        // Populated task-by-task via accumulate_locally(); drained once per
        // subworld via finalize_into() into the universe-resident result.
        static inline vecfuncT Kf_local_;
        static inline bool Kf_local_initialized_ = false;
        static inline long Kf_local_world_id_ = -1;
        // Node-shared accumulator for accumulation_mode_==2 (node-local reduction).
        // Lives in the node sub-World (Split_type(SHARED)); collectively constructed
        // by every node rank in finalize_stage1() via ensure_node_accumulator().
        // Each node rank's Kf_node_ handle references the same node-distributed
        // functions. Stage1 reduces Kf_local_ -> Kf_node_ (intra-node); stage2 does
        // one Kf_node_ -> universe gaxpy (only nnode inter-node scatterers per key).
        static inline vecfuncT Kf_node_;
        static inline bool Kf_node_initialized_ = false;
        static inline long Kf_node_world_id_ = -1;

        // ---- Option C0: coalesced finalize transfer ----
        // Cached receiver endpoints, one per transport world (universe for the
        // ->universe drains, the node sub-World for stage1). Rebuilt only when the
        // world id changes; no per-iteration WorldObject construction.
        static inline std::shared_ptr<FinalizeReducer<T, NDIM>> universe_reducer_;
        static inline long universe_reducer_world_id_ = -1;
        static inline std::shared_ptr<FinalizeReducer<T, NDIM>> node_reducer_;
        static inline long node_reducer_world_id_ = -1;
        // Once-per-finalize guards: the coalesced drain (with its collective
        // readiness barrier) must run exactly once per rank, uniformly, even on
        // ranks that owned no tasks (replicated taskq => first call lands on the
        // same loop iteration on every rank). Replaces the old data-less
        // early-return guards. Reset in cleanup().
        static inline bool finalize_stage1_done_ = false;
        static inline bool finalize_universe_done_ = false;   // shared by finalize_into + stage2 node->universe

        /// chunk size (entries) for coalesced_gaxpy, derived from k so a worst-case
        /// all-coeff chunk (~(2k)^NDIM*sizeof(T) per node) stays under the eager AM
        /// buffer (~1.5 MB default); we target ~1 MB.
        static std::size_t finalize_chunk_entries() {
            const long k = FunctionDefaults<NDIM>::get_k();
            const double node_bytes = std::pow(double(2 * k), double(NDIM)) * double(sizeof(T));
            const double target = 1024.0 * 1024.0;
            return std::max<std::size_t>(1, std::size_t(target / std::max(1.0, node_bytes)));
        }

        /// cached reducer on the universe transport world (collective on first use)
        static FinalizeReducer<T, NDIM>& get_universe_reducer(World& world) {
            const long wid = long(world.id());
            if (not universe_reducer_ or universe_reducer_world_id_ != wid) {
                universe_reducer_ = std::make_shared<FinalizeReducer<T, NDIM>>(world);
                universe_reducer_world_id_ = wid;
            }
            return *universe_reducer_;
        }
        /// cached reducer on the node sub-World transport (collective on first use)
        static FinalizeReducer<T, NDIM>& get_node_reducer(World& world) {
            const long wid = long(world.id());
            if (not node_reducer_ or node_reducer_world_id_ != wid) {
                node_reducer_ = std::make_shared<FinalizeReducer<T, NDIM>>(world);
                node_reducer_world_id_ = wid;
            }
            return *node_reducer_;
        }

        struct VfPrefetchState {
            Batch_1D current_range;
            vecfuncT current_data;
            bool has_current = false;

            Batch_1D next_range;
            vecfuncT next_data;
            bool has_next = false;

            Batch_1D next_hint;
            bool has_hint = false;
        };
        static inline VfPrefetchState vf_prefetch_;

        /// Self-contained prefetch state for the rotating k-batch (mo_bra + mo_ket
        /// paired by the inner-sum index) used by small_memory_mt_owner. Independent
        /// of vf_prefetch_ so the row-owner algorithm can rotate this dimension
        /// while the existing symmetric path keeps rotating its own (column / vf).
        struct KBatchPrefetchState {
            Batch_1D current_range;
            vecfuncT current_bra;
            vecfuncT current_ket;
            bool has_current = false;

            Batch_1D next_range;
            vecfuncT next_bra;
            vecfuncT next_ket;
            bool has_next = false;

            Batch_1D next_hint;
            bool has_hint = false;
        };
        static inline KBatchPrefetchState kbatch_prefetch_;

        /// Cloud-batch k-batch prefetch state (Step 5b).
        /// Single-slot: prefetch_next_bra_async (called at run() TAIL) issues
        /// async finds for the NEXT task's bra/ket records and stores the
        /// futures here; operator() at the next task's start consumes them
        /// via deserialize_batch. ensure_cache_world resets the slot on
        /// subworld change (same lifecycle as the copy-path caches).
        /// (A two-slot run-head/run-tail variant was tried and reverted: the
        /// synchronized cross-rank issuance burst it produced triggered
        /// 2-sided MPI deadlocks on large batch replies — see Fix B in the
        /// plan's deferred items.)
        // One in-flight (or completed) cloud prefetch of a bra+ket k-batch.
        struct CloudKBatchSlot {
            bool valid = false;
            Batch_1D range;
            typename Cloud::keyT bra_key = 0;
            typename Cloud::keyT ket_key = 0;
            Future<typename Cloud::batch_bytesT> bra_fut;
            Future<typename Cloud::batch_bytesT> ket_fut;
        };
        // Double buffer for the row-owner cloud pipeline (Design B): `next` is
        // issued one task ahead (overlapping the current task's compute) and
        // promoted to `current` at the head of the next task, where it is
        // consumed. See cloud_kbatch_pipeline_advance.
        struct CloudKBatchPrefetch {
            CloudKBatchSlot current;   // promoted from the previous task's `next`; consumed by this task
            CloudKBatchSlot next;      // issued ahead during this task's compute; for the next task
        };
        static inline CloudKBatchPrefetch cloud_kb_;

        /// Bounded LRU for the one-set sym_p2p algorithm's cloud-fetched ψ batches.
        /// Keyed by the cloud RECORD key (encodes salt+dim+range -> globally unique
        /// per K-application, so a stale entry from a previous application can never
        /// false-hit). The cache-aware task ordering (shuffle_partition_by_owner)
        /// keeps an owner's local batch and the current transient adjacent, so a
        /// small capacity captures the reuse while bounding the footprint (smallmem).
        /// Capacity is set in prepare_owner_assignment to (owned-per-rank + 1).
        /// Cleared subworld-scoped via clear_local_caches / ensure_cache_world,
        /// exactly like bra_cache_ — NOT the cloud-side cache_result path (that one
        /// segfaults across protocol transitions; see the plan's deferred items).
        struct SymBatchLRU {
            std::vector<std::pair<typename Cloud::keyT, vecfuncT>> slots;  // front = most-recently-used
            std::size_t capacity = 2;
        };
        static inline SymBatchLRU sym_lru_;

        /// Design-B prefetch (Option A) for sym_p2p: a STRICT single-slot double
        /// buffer mirroring the asymmetric cloud_kb_ discipline (NOT a flat map — an
        /// earlier flat cap-4 map with drop-oldest allowed too many concurrent /
        /// orphaned in-flight byte requests and corrupted the 2-sided transport on
        /// LARGE batch replies → MPI_ERR_TRUNCATE / wrong-sized deserialize; see the
        /// plan's "Fix B" precedent). At most ONE request is issued per task (the next
        /// owned task's transient), and every pre-compute promotes `next`→`current`
        /// (retiring the old current), so at most TWO are ever in flight — the same
        /// bound the proven asymmetric path uses. Cleared subworld-scoped in
        /// clear_local_caches (same lifecycle as cloud_kb_).
        struct SymPrefetchSlot {
            bool valid = false;
            typename Cloud::keyT key = 0;
            Future<typename Cloud::batch_bytesT> fut;
        };
        static inline SymPrefetchSlot sym_pf_current_;   // promoted from prev task's `next`; consumed this task
        static inline SymPrefetchSlot sym_pf_next_;       // issued ahead during this task's compute

        static inline long cache_world_id_ = -1;

        static void clear_kbatch_prefetch() {
            kbatch_prefetch_.current_bra.clear();
            kbatch_prefetch_.current_ket.clear();
            kbatch_prefetch_.next_bra.clear();
            kbatch_prefetch_.next_ket.clear();
            kbatch_prefetch_.has_current = false;
            kbatch_prefetch_.has_next = false;
            kbatch_prefetch_.has_hint = false;
            kbatch_prefetch_.current_range = Batch_1D();
            kbatch_prefetch_.next_range = Batch_1D();
            kbatch_prefetch_.next_hint = Batch_1D();
            // drop any in-flight cloud prefetch (both slots reset to unset)
            cloud_kb_.current = CloudKBatchSlot();
            cloud_kb_.next = CloudKBatchSlot();
        }

        /// pre-computed owner map: (col_begin, row_begin) -> owner rank
        /// populated by prepare_owner_assignment() using the fold algorithm
        std::map<std::pair<long,long>, long> owner_map_;

        /// peel indices chosen by m-flex search (empty if not used).
        /// Stored for diagnostic logging in prepare_owner_assignment().
        std::vector<long> chosen_peel_;

        /// if true, shuffle per-owner task order to reduce synchronized fetch contention
        /// disabled: shuffling destroys row-range locality needed for cache reuse and prefetch hits
        bool shuffle_task_order_ = false;

        bool use_owner_aware_fetch() const {
            return (algorithm_==small_memory_symmetric_mt_owner or algorithm_==small_memory_mt_owner
                    or algorithm_==small_memory_symmetric_p2p_owner)
                   and not replicate_for_debug_;
        }

        /// true iff the new symmetric round-robin p2p-owner algorithm is active.
        /// Triangular grid over the exchange_sym_owner_split batches, owner assigned
        /// by exchange_sym_round_robin_assign, one-set cloud fetch (bra==ket).
        bool use_sym_p2p_owner_algorithm() const {
            return algorithm_==small_memory_symmetric_p2p_owner;
        }

        /// true iff the new row-owner algorithm is active. The full-grid partition
        /// produced by row_owner_partition is INCOMPATIBLE with the upper-triangle
        /// symmetric kernels (compute_*_batch_in_symmetric_matrix), so this branch
        /// must fire unconditionally for the algorithm — including when
        /// replicate_for_debug_ is on. The owner-aware fetch path (cache + prefetch)
        /// is still gated separately by use_owner_aware_fetch().
        bool use_row_owner_algorithm() const {
            return algorithm_==small_memory_mt_owner;
        }

        /// true iff input batches should be fetched from the cloud's owner-pinned
        /// batch container instead of copy()ing from the universe functions.
        /// Only honored for the row-owner algorithm; default off (A/B baseline).
        bool use_cloud_batch_fetch() const {
            return use_cloud_batch_fetch_ and (use_row_owner_algorithm() or use_sym_p2p_owner_algorithm());
        }

        static bool same_range(const Batch_1D& a, const Batch_1D& b) {
            return (a.begin == b.begin) and (a.end == b.end);
        }

        static bool hint_matches_range(const Batch_1D& hint, const Batch_1D& range) {
            if (hint.begin != range.begin) return false;
            if (hint.end < 0) return true;
            return hint.end == range.end;
        }

        static Batch_1D normalize_range(const Batch_1D& range, const long full_size) {
            Batch_1D normalized = range;
            if (normalized.is_full_size()) normalized.end = full_size;
            return normalized;
        }

        static void clear_vf_prefetch() {
            vf_prefetch_.current_data.clear();
            vf_prefetch_.next_data.clear();
            vf_prefetch_.has_current = false;
            vf_prefetch_.has_next = false;
            vf_prefetch_.has_hint = false;
            vf_prefetch_.current_range = Batch_1D();
            vf_prefetch_.next_range = Batch_1D();
            vf_prefetch_.next_hint = Batch_1D();
        }

        static void clear_local_caches() {
            bra_cache_.clear();
            ket_cache_.clear();
            held_vf_cache_.clear();
            sym_lru_.slots.clear();   // capacity (set in prepare_owner_assignment) survives
            sym_pf_current_ = SymPrefetchSlot();   // drop any in-flight prefetch from the old subworld
            sym_pf_next_ = SymPrefetchSlot();
            clear_vf_prefetch();
            clear_kbatch_prefetch();
        }

        void add_owner_fetch_time(const double cpu0, const double cpu1) const {
            if (use_owner_aware_fetch()) owner_fetch_timer += long((cpu1 - cpu0) * 1000l);
        }

        void add_owner_compute_time(const double cpu0, const double cpu1) const {
            if (use_owner_aware_fetch()) owner_compute_timer += long((cpu1 - cpu0) * 1000l);
            if (exch_task_profile_enabled()) prof_.compute_cpu = cpu1 - cpu0;
        }

        /// wall-time companion to add_owner_compute_time. Called from operator() to
        /// accumulate wall time wrapping the compute kernel call. Combined with the
        /// CPU-time owner_compute_timer, lets the user attribute the gap to
        /// communication/fences that happened during the compute phase.
        void add_owner_compute_wall_time(const double wall0, const double wall1) const {
            if (use_owner_aware_fetch()) owner_compute_wall_timer += long((wall1 - wall0) * 1000l);
            if (exch_task_profile_enabled()) {
                // wall0 = compute start (post fetch), wall1 = compute end. Everything
                // from task entry to compute start is data-wait (bra+vf[+ket] fetch).
                prof_.compute_wall = wall1 - wall0;
                prof_.wait_for_data_wall = wall0 - prof_.wall_start;
            }
        }

        void ensure_cache_world(World& world) const {
            if (cache_world_id_ != world.id()) {
                clear_local_caches();
                cache_world_id_ = world.id();
            }
        }

        vecfuncT fetch_batch_with_cache(World& world, const vecfuncT& batch, const Batch_1D& range,
                                        std::unordered_map<long, functionT>& cache) const {
            const double cpu0 = process_cpu_time();
            MADNESS_CHECK_THROW(long(batch.size())==range.size(),
                                "batch/range size mismatch in fetch_batch_with_cache");
            ensure_cache_world(world);
            vecfuncT result;
            result.reserve(batch.size());
            for (long local_index = 0; local_index < long(batch.size()); ++local_index) {
                const long global_index = range.begin + local_index;
                auto it = cache.find(global_index);
                if (it==cache.end()) {
                    functionT local_copy = copy(world, batch[local_index], false);
                    it = cache.emplace(global_index, std::move(local_copy)).first;
                }
                result.push_back(it->second);
            }
            world.gop.fence();
            const double cpu1 = process_cpu_time();
            add_owner_fetch_time(cpu0, cpu1);
            return result;
        }

        vecfuncT fetch_batch_transient(World& world, const vecfuncT& batch) const {
            const double cpu0 = process_cpu_time();
            vecfuncT result;
            result.reserve(batch.size());
            for (const auto& f : batch) {
                result.push_back(copy(world, f, false));
            }
            world.gop.fence();
            const double cpu1 = process_cpu_time();
            add_owner_fetch_time(cpu0, cpu1);
            return result;
        }

        vecfuncT fetch_range_with_cache(World& world, const vecfuncT& source, const Batch_1D& range,
                                        std::unordered_map<long, functionT>& cache) const {
            const double cpu0 = process_cpu_time();
            ensure_cache_world(world);
            vecfuncT result;
            result.reserve(range.size());
            for (long global_index = range.begin; global_index < range.end; ++global_index) {
                auto it = cache.find(global_index);
                if (it==cache.end()) {
                    functionT local_copy = copy(world, source[global_index], false);
                    it = cache.emplace(global_index, std::move(local_copy)).first;
                }
                result.push_back(it->second);
            }
            world.gop.fence();
            const double cpu1 = process_cpu_time();
            add_owner_fetch_time(cpu0, cpu1);
            return result;
        }

        vecfuncT fetch_range_transient(World& world, const vecfuncT& source, const Batch_1D& range) const {
            const double cpu0 = process_cpu_time();
            vecfuncT result;
            result.reserve(range.size());
            for (long global_index = range.begin; global_index < range.end; ++global_index) {
                result.push_back(copy(world, source[global_index], false));
            }
            world.gop.fence();
            const double cpu1 = process_cpu_time();
            add_owner_fetch_time(cpu0, cpu1);
            return result;
        }

        /// True iff `key` is currently resident in the sym_p2p LRU.
        bool sym_lru_contains(const typename Cloud::keyT& key) const {
            for (const auto& s : sym_lru_.slots) if (s.first == key) return true;
            return false;
        }

        /// Insert a freshly fetched/deserialized batch at the LRU front, evicting the
        /// tail past capacity, and return a ref to it.
        const vecfuncT& sym_lru_insert(const typename Cloud::keyT& key, vecfuncT&& data) const {
            sym_lru_.slots.insert(sym_lru_.slots.begin(), std::make_pair(key, std::move(data)));
            if (sym_lru_.slots.size() > std::max<std::size_t>(1, sym_lru_.capacity))
                sym_lru_.slots.pop_back();
            return sym_lru_.slots.front().second;
        }

        /// Bounded-LRU fetch for the one-set sym_p2p path, keyed by cloud record key.
        /// Three tiers: (1) LRU hit — already-deserialized batch reused (no p2p);
        /// (2) Design-B prefetch consume — bytes were requested one task ahead by
        /// sym_pipeline_advance, so we deserialize the (likely already-landed) future
        /// instead of stalling; (3) synchronous miss — a cold fetch_batch_p2p. Every
        /// fetched batch is inserted at the LRU front (evicting the tail past capacity).
        /// The cache-aware task ordering keeps an owner's local batch and the rotating
        /// transient adjacent, so a small capacity (owned-per-rank + 1) captures the
        /// reuse while bounding the footprint. Returns a const ref into the cache; the
        /// caller must copy the (shallow) handles out before the next fetch (a later
        /// insert may relocate the slot). World-scoped via ensure_cache_world.
        const vecfuncT& sym_batch_lru_fetch(World& world, Cloud& cloud,
                                            const typename Cloud::keyT& key) const {
            const double cpu0 = process_cpu_time();
            ensure_cache_world(world);
            // (1) LRU hit
            for (std::size_t s = 0; s < sym_lru_.slots.size(); ++s) {
                if (sym_lru_.slots[s].first == key) {
                    if (s != 0)
                        std::rotate(sym_lru_.slots.begin(),
                                    sym_lru_.slots.begin() + s,
                                    sym_lru_.slots.begin() + s + 1);
                    ++sym_cache_hits_;
                    prof_.observe_fetch_tier(0);
                    add_owner_fetch_time(cpu0, process_cpu_time());
                    return sym_lru_.slots.front().second;
                }
            }
            // (2) Design-B prefetch consume: the bytes for this key were requested one
            //     task ahead; await + deserialize (overlapped with the prior compute).
            //     Check `current` (the promoted slot) first, then `next`.
            for (SymPrefetchSlot* slot : {&sym_pf_current_, &sym_pf_next_}) {
                if (slot->valid and slot->key == key) {
                    vecfuncT data = cloud.deserialize_batch_p2p<T, NDIM>(
                            world, slot->fut, key, /*cache_result=*/false);
                    *slot = SymPrefetchSlot();   // retire the consumed slot
                    ++sym_prefetch_hits_;
                    prof_.observe_fetch_tier(1);
                    add_owner_fetch_time(cpu0, process_cpu_time());
                    return sym_lru_insert(key, std::move(data));
                }
            }
            // (3) synchronous miss
            ++sym_cache_misses_;
            prof_.observe_fetch_tier(2);
            vecfuncT data = cloud.fetch_batch_p2p<T, NDIM>(world, key, /*cache_result=*/false);
            add_owner_fetch_time(cpu0, process_cpu_time());
            return sym_lru_insert(key, std::move(data));
        }

        /// custom partitioning for the exchange operator in exchangeoperator.h

        /// arguments are: result[i] += sum_k vket[k] \int 1/r vbra[k] f[i]
        /// with f and vbra being batched, result and vket being passed on as a whole
        class MacroTaskPartitionerExchange : public MacroTaskPartitioner {
        public:
            MacroTaskPartitionerExchange(const bool symmetric, const long min_batch_size_input,
                                         const long max_batch_size_input,
                                         const bool row_owner = false,
                                         const bool sym_p2p_owner = false,
                                         const long granularity_level = 1)
                    : symmetric(symmetric), row_owner_(row_owner),
                      sym_p2p_owner_(sym_p2p_owner),
                      granularity_level_(std::max<long>(1, granularity_level)) {
                const long min_bs = std::max<long>(1, min_batch_size_input);
                const long max_bs = std::max<long>(min_bs, std::max<long>(1, max_batch_size_input));
                min_batch_size=min_bs;
                max_batch_size=max_bs;
            }

            bool symmetric = false;
            /// when true, ignore min/max batch size and produce exactly nsubworld
            /// (or fewer, if n < nsubworld) batches per dimension via even-remainder
            /// spread. Used by small_memory_mt_owner.
            bool row_owner_ = false;
            /// when true, produce a lower-triangular grid over M=level*nsubworld
            /// batches (single split, both dims) via exchange_sym_owner_split. Used by
            /// small_memory_symmetric_p2p_owner.
            bool sym_p2p_owner_ = false;
            /// granularity level for the sym-p2p split (>=1).
            long granularity_level_ = 1;

            /// Row-owner batch boundaries for this partitioner's nsubworld.
            /// Thin wrapper over the shared free helper exchange_row_owner_split,
            /// which is the single source of truth shared with the cloud
            /// batch-storage side (see exchangeoperator.h).
            std::vector<Batch_1D> row_owner_split(std::size_t n) const {
                return exchange_row_owner_split(n, long(nsubworld));
            }

            /// Symmetric p2p-owner batch boundaries (single split, granularity-aware).
            /// Shared source of truth with store_batches and prepare_owner_assignment.
            std::vector<Batch_1D> sym_owner_split(std::size_t n) const {
                return exchange_sym_owner_split(n, long(nsubworld), granularity_level_);
            }

            partitionT do_partitioning(const std::size_t& vsize1, const std::size_t& vsize2,
                                       const std::string policy) const override {

                if (sym_p2p_owner_) {
                    // Lower-triangular grid over a single granularity-aware split.
                    // col = input[0] = batch i, row = input[1] = batch j, with j <= i.
                    // owner assigned later by prepare_owner_assignment (round-robin).
                    std::vector<Batch_1D> batches = sym_owner_split(vsize1);
                    partitionT result;
                    for (long i = 0; i < long(batches.size()); ++i) {
                        for (long j = 0; j <= i; ++j) {
                            Batch batch(batches[i], batches[j], _);
                            double priority = compute_priority(batch);
                            result.push_back(std::make_pair(batch, priority));
                        }
                    }
                    return result;
                }

                if (row_owner_) {
                    // Strict bs = ceil(n/nsubworld), even remainder spread, full grid.
                    // For asymmetric exchange only; no symmetry exploitation here.
                    std::vector<Batch_1D> col_batches = row_owner_split(vsize1); // input[0]: held / vf / i-batch
                    std::vector<Batch_1D> row_batches = row_owner_split(vsize2); // input[1]: rotating / bra / k-batch
                    partitionT result;
                    for (const auto& c : col_batches) {
                        for (const auto& r : row_batches) {
                            Batch batch(c, r, _);
                            double priority = compute_priority(batch);
                            result.push_back(std::make_pair(batch, priority));
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
        MacroTaskExchangeSimple(const long nresult, const double lo, const double mul_tol, const bool symmetric,
                                const long min_batch_size,
                                const long max_batch_size,
                                const Algorithm algorithm = multiworld_efficient,
                                const long batch_granularity = 1)
                : nresult(nresult), lo(lo), mul_tol(mul_tol), symmetric(symmetric), algorithm_(algorithm) {
            batch_granularity_ = std::max<long>(1, batch_granularity);
            const bool row_owner = (algorithm == small_memory_mt_owner);
            const bool sym_p2p = (algorithm == small_memory_symmetric_p2p_owner);
            partitioner.reset(new MacroTaskPartitionerExchange(symmetric, min_batch_size, max_batch_size,
                                                               row_owner, sym_p2p, batch_granularity_));
            name="MacroTaskExchangeSimple";
        }

        /// Per-task entry, indexed by (col_begin, row_begin) and weighted
        /// by partition priority (= nrow * ncol).
        struct FoldTaskEntry {
            std::pair<long,long> key;  // (col_begin, row_begin)
            double priority;
        };

        /// One row group of the symmetric upper triangle (all tasks sharing
        /// a row_begin). Cost is the sum of task priorities in this group.
        struct FoldRowGroup {
            long row_begin = 0;
            double cost = 0.0;
            std::vector<FoldTaskEntry> tasks;
        };

        /// Build row groups from a partition, ordered by row_begin.
        /// Shared by fold_and_assign() and the m-flex search.
        static std::vector<FoldRowGroup> build_row_groups(
                const MacroTaskPartitioner::partitionT& partition) {
            std::map<long, FoldRowGroup> row_group_map;
            for (const auto& [batch, prio] : partition) {
                const Batch_1D& col_range = batch.input[0];
                const Batch_1D& row_range = (batch.input.size() > 1) ? batch.input[1] : batch.input[0];
                auto& rg = row_group_map[row_range.begin];
                rg.row_begin = row_range.begin;
                rg.cost += prio;
                rg.tasks.push_back({{col_range.begin, row_range.begin}, prio});
            }
            std::vector<FoldRowGroup> row_groups;
            row_groups.reserve(row_group_map.size());
            for (auto& [key, rg] : row_group_map) row_groups.push_back(std::move(rg));
            return row_groups;
        }

        /// Pre-compute a load-balanced owner assignment for all tasks in the partition.
        ///
        /// Uses a fold algorithm inspired by the triangle-to-rectangle transformation:
        /// row groups from opposite ends of the triangular task matrix are paired so
        /// that each "rectangle row" has approximately equal task count. Rectangle rows
        /// are then assigned to ranks via greedy (least-loaded-first) scheduling.
        ///
        /// When use_mflex_ is true, an m-flex peel search runs first and selects a
        /// set of m row groups (m = R mod 2*nsubworld) to peel out of the fold and
        /// distribute round-robin across the remaining rectangle rows. This restores
        /// load balance when R is not a clean multiple of 2*nsubworld (e.g. when a
        /// runt batch dilutes one row group's cost).
        ///
        /// Called automatically by the MacroTask framework (via SFINAE hook) after
        /// partitioning and before per-task owner_hint queries.
        void prepare_owner_assignment(const MacroTaskPartitioner::partitionT& partition, long nsubworld) {
            if (!use_owner_aware_fetch() || nsubworld <= 0) return;

            // Row-owner (asymmetric, held i-batch): owner depends only on input[0]
            // (= col / vf / i-batch). All tasks sharing a col.begin go to the same
            // rank, in col-begin order. With nbatch = min(nsubworld, n_orb)
            // batches per dimension (see MacroTaskPartitionerExchange::row_owner_split),
            // this gives an injective rank assignment up to nsubworld ranks.
            if (algorithm_ == small_memory_mt_owner) {
                owner_map_.clear();
                chosen_peel_.clear();
                std::map<long, long> col_to_rank;
                long next_rank = 0;
                for (const auto& [batch, prio] : partition) {
                    MADNESS_CHECK_THROW(batch.input.size() >= 2,
                                        "small_memory_mt_owner expects 2D batches (col, row)");
                    const Batch_1D& col = batch.input[0];
                    const Batch_1D& row = batch.input[1];
                    auto it = col_to_rank.find(col.begin);
                    if (it == col_to_rank.end()) {
                        const long rank = next_rank % nsubworld;
                        ++next_rank;
                        it = col_to_rank.emplace(col.begin, rank).first;
                    }
                    owner_map_[{col.begin, row.begin}] = it->second;
                }
                return;
            }

            // Symmetric round-robin p2p-owner: triangular grid over the single
            // granularity-aware split. Owner per (i,j) batch-index pair comes from
            // exchange_sym_round_robin_assign (balanced + locality: owner holds one
            // operand). Map batch-begin offsets -> indices to key the rr assignment.
            if (algorithm_ == small_memory_symmetric_p2p_owner) {
                owner_map_.clear();
                chosen_peel_.clear();
                const std::vector<Batch_1D> split =
                        exchange_sym_owner_split(nresult, nsubworld, batch_granularity_);
                const long M = long(split.size());
                // Size the reuse LRU to hold an owner's local batch(es) plus one
                // rotating transient: owned-per-rank = ceil(M/nsubworld), +1 for the
                // transient. With the cache-aware ordering this captures the reuse
                // while keeping the footprint bounded (the smallmem invariant).
                const long owned_per_rank = (M + nsubworld - 1) / nsubworld;
                sym_lru_.capacity = std::size_t(std::max<long>(2, owned_per_rank + 1));
                std::map<long,long> begin_to_idx;
                for (long k = 0; k < M; ++k) begin_to_idx[split[k].begin] = k;
                const std::map<std::pair<long,long>,long> rr =
                        exchange_sym_round_robin_assign(nsubworld, M);
                for (const auto& [batch, prio] : partition) {
                    MADNESS_CHECK_THROW(batch.input.size() >= 2,
                                        "small_memory_symmetric_p2p_owner expects 2D batches (col, row)");
                    const Batch_1D& col = batch.input[0];
                    const Batch_1D& row = batch.input[1];
                    auto ic = begin_to_idx.find(col.begin);
                    auto jr = begin_to_idx.find(row.begin);
                    MADNESS_CHECK_THROW(ic != begin_to_idx.end() and jr != begin_to_idx.end(),
                                        "sym_p2p: task batch offset not found in split");
                    const long a = std::max(ic->second, jr->second);   // rr keyed by (i>=j)
                    const long b = std::min(ic->second, jr->second);
                    auto it = rr.find({a, b});
                    owner_map_[{col.begin, row.begin}] =
                            (it != rr.end()) ? it->second : (a % nsubworld);
                }
                return;
            }

            // Existing path (symmetric folded fold + LPT + m-flex peel) is
            // only valid for the symmetric algorithm.
            if (!symmetric) return;

            // Small-problem fallback: the fold needs R = #row_groups >= 2*nsubworld
            // to give every rank at least one rectangle row. When that floor cannot
            // be reached (typically n_MO < 2*nsubworld, where the partitioner already
            // produced bs = 1), bypass the fold and LPT-distribute individual tasks
            // directly across ranks. Result is identical-by-rank because LPT is
            // deterministic (sort by priority desc, tie-break by (col_begin,
            // row_begin)).
            std::vector<FoldRowGroup> row_groups = build_row_groups(partition);
            const long R = static_cast<long>(row_groups.size());
            if (R > 0 && R < 2 * nsubworld) {
                owner_map_.clear();
                chosen_peel_.clear();
                struct TaskRef {
                    std::pair<long,long> key;
                    double priority;
                };
                std::vector<TaskRef> all_tasks;
                all_tasks.reserve(partition.size());
                for (const auto& rg : row_groups) {
                    for (const auto& te : rg.tasks) {
                        all_tasks.push_back({te.key, te.priority});
                    }
                }
                std::sort(all_tasks.begin(), all_tasks.end(),
                          [](const TaskRef& a, const TaskRef& b) {
                              if (a.priority != b.priority) return a.priority > b.priority;
                              return a.key < b.key;
                          });
                using heap_entry = std::pair<double, long>;
                std::priority_queue<heap_entry, std::vector<heap_entry>,
                                    std::greater<heap_entry>> heap;
                for (long p = 0; p < nsubworld; ++p) heap.push({0.0, p});
                for (const auto& t : all_tasks) {
                    auto [load, rank] = heap.top();
                    heap.pop();
                    owner_map_[t.key] = rank;
                    heap.push({load + t.priority, rank});
                }
                return;
            }

            std::vector<long> peel;
            if (use_mflex_) {
                peel = find_best_peel(partition, nsubworld);
            }
            chosen_peel_ = peel;
            owner_map_ = fold_and_assign(partition, nsubworld, peel);
        }

        /// Fold the triangular task list and assign owners for load balance.
        ///
        /// 1. Group tasks by row range (input[1]), giving R row groups ordered by
        ///    row index. Cost is accumulated from task priority (batch area), not
        ///    raw task count, so runt batches are weighted correctly.
        /// 2. Determine the peel set:
        ///    - If peel_indices is non-empty, use it directly. Constraint:
        ///      (R - |peel|) must be even.
        ///    - Else, use defaults: {R/2} when R is odd, {} when R is even.
        /// 3. Fold the non-peeled row groups symmetrically by index:
        ///    pair non_peeled[k] with non_peeled[N-1-k] to form rectangle rows.
        /// 4. Distribute peeled row groups' tasks round-robin across the
        ///    rectangle rows.
        /// 5. Greedy (LPT) assignment: rectangle rows to ranks by cost.
        /// 6. Build a map from (col_begin, row_begin) -> owner rank for each task.
        static std::map<std::pair<long,long>, long> fold_and_assign(
                const MacroTaskPartitioner::partitionT& partition, long nsubworld,
                const std::vector<long>& peel_indices = {}) {

            std::vector<FoldRowGroup> row_groups = build_row_groups(partition);
            const long R = static_cast<long>(row_groups.size());
            if (R == 0) return {};

            // Resolve peel set. Empty input means "use default":
            //   even R -> no peel; odd R -> peel = {R/2}
            std::vector<long> peel = peel_indices;
            if (peel.empty() && (R % 2 != 0)) peel.push_back(R / 2);
            std::sort(peel.begin(), peel.end());
            peel.erase(std::unique(peel.begin(), peel.end()), peel.end());
            for (long p : peel) {
                MADNESS_CHECK_THROW(p >= 0 && p < R, "peel index out of range");
            }
            MADNESS_CHECK_THROW((R - static_cast<long>(peel.size())) % 2 == 0,
                                "R - |peel| must be even for clean fold");

            // Non-peeled indices, sorted (input is sorted, peel is sorted -> use set_difference).
            std::vector<long> non_peeled;
            non_peeled.reserve(R);
            {
                std::vector<long> all(R);
                std::iota(all.begin(), all.end(), 0L);
                std::set_difference(all.begin(), all.end(),
                                    peel.begin(), peel.end(),
                                    std::back_inserter(non_peeled));
            }
            const long N = static_cast<long>(non_peeled.size());
            const long half = N / 2;

            // -- Step 3: fold non-peeled row groups symmetrically by index --
            struct RectRow {
                double cost = 0.0;
                std::vector<std::pair<long,long>> task_keys;
            };
            std::vector<RectRow> rect_rows(half);
            for (long k = 0; k < half; ++k) {
                const long top_idx = non_peeled[k];
                const long bottom_idx = non_peeled[N - 1 - k];
                const auto& top = row_groups[top_idx];
                const auto& bottom = row_groups[bottom_idx];
                auto& rr = rect_rows[k];
                rr.cost = top.cost + bottom.cost;
                for (const auto& te : top.tasks) rr.task_keys.push_back(te.key);
                for (const auto& te : bottom.tasks) rr.task_keys.push_back(te.key);
            }

            // -- Step 4: distribute peeled row groups round-robin across rect rows --
            // Concatenate peeled tasks in (sorted) row-group order, then assign
            // task t to rect_rows[t % half].
            std::map<std::pair<long,long>, long> owner_map;
            if (half > 0) {
                long t_idx = 0;
                for (long p : peel) {
                    const auto& rg = row_groups[p];
                    for (const auto& te : rg.tasks) {
                        auto& rr = rect_rows[t_idx % half];
                        rr.task_keys.push_back(te.key);
                        rr.cost += te.priority;
                        ++t_idx;
                    }
                }
            } else {
                // half == 0: every row group is peeled (or R == 1). Assign all
                // tasks to rank 0. (Reachable only on extreme small inputs.)
                for (long p : peel) {
                    for (const auto& te : row_groups[p].tasks) owner_map[te.key] = 0;
                }
                return owner_map;
            }

            // -- Step 5: greedy (LPT) assignment of rectangle rows to ranks --
            using heap_entry = std::pair<double, long>;
            std::priority_queue<heap_entry, std::vector<heap_entry>, std::greater<heap_entry>> heap;
            for (long p = 0; p < nsubworld; ++p) heap.push({0.0, p});

            std::sort(rect_rows.begin(), rect_rows.end(),
                      [](const RectRow& a, const RectRow& b) { return a.cost > b.cost; });

            for (const auto& rr : rect_rows) {
                auto [load, rank] = heap.top();
                heap.pop();
                for (const auto& key : rr.task_keys) owner_map[key] = rank;
                heap.push({load + rr.cost, rank});
            }
            return owner_map;
        }

        /// Compute the max-min rank load delta resulting from running
        /// fold_and_assign with the given peel set. Used as the scoring
        /// function for the m-flex peel search. Recomputes the rect-row
        /// costs and runs LPT (without populating an owner_map -- cheaper
        /// than fold_and_assign for inner-loop scoring).
        static double load_delta_for_peel(
                const std::vector<FoldRowGroup>& row_groups, long nsubworld,
                const std::vector<long>& peel_sorted) {
            const long R = static_cast<long>(row_groups.size());
            const long M = static_cast<long>(peel_sorted.size());
            const long N = R - M;
            if (N % 2 != 0 || nsubworld <= 0) return std::numeric_limits<double>::infinity();
            const long half = N / 2;
            if (half == 0) return std::numeric_limits<double>::infinity();

            // Build non_peeled (set_difference)
            std::vector<long> non_peeled;
            non_peeled.reserve(N);
            std::vector<long> all(R);
            std::iota(all.begin(), all.end(), 0L);
            std::set_difference(all.begin(), all.end(),
                                peel_sorted.begin(), peel_sorted.end(),
                                std::back_inserter(non_peeled));

            // Rect-row costs from symmetric pair fold
            std::vector<double> rr_cost(half, 0.0);
            for (long k = 0; k < half; ++k) {
                rr_cost[k] = row_groups[non_peeled[k]].cost
                           + row_groups[non_peeled[N - 1 - k]].cost;
            }

            // Round-robin peeled task priorities across rect rows
            long t_idx = 0;
            for (long p : peel_sorted) {
                for (const auto& te : row_groups[p].tasks) {
                    rr_cost[t_idx % half] += te.priority;
                    ++t_idx;
                }
            }

            // LPT to nsubworld ranks
            using heap_entry = std::pair<double, long>;
            std::priority_queue<heap_entry, std::vector<heap_entry>, std::greater<heap_entry>> heap;
            for (long p = 0; p < nsubworld; ++p) heap.push({0.0, p});
            std::sort(rr_cost.begin(), rr_cost.end(), std::greater<double>());

            std::vector<double> rank_load(nsubworld, 0.0);
            for (double c : rr_cost) {
                auto [load, rank] = heap.top();
                heap.pop();
                rank_load[rank] = load + c;
                heap.push({load + c, rank});
            }

            const double mx = *std::max_element(rank_load.begin(), rank_load.end());
            const double mn = *std::min_element(rank_load.begin(), rank_load.end());
            return mx - mn;
        }

        /// Generate candidate peel sets of size m for R row groups.
        /// Heuristics target peels that keep the remaining row groups'
        /// fold pair-sums balanced; exhaustive enumeration is added when
        /// C(R, m) <= max_exhaustive.
        static std::vector<std::vector<long>> generate_peel_candidates(
                long R, long m, const std::vector<FoldRowGroup>& row_groups,
                long max_exhaustive) {
            std::set<std::vector<long>> cands;
            if (m == 0) { cands.insert({}); }
            else if (m >= R) {
                std::vector<long> all(R);
                std::iota(all.begin(), all.end(), 0L);
                cands.insert(all);
            } else {
                auto add_sorted = [&](std::vector<long> v) {
                    std::sort(v.begin(), v.end());
                    v.erase(std::unique(v.begin(), v.end()), v.end());
                    if (static_cast<long>(v.size()) == m) cands.insert(std::move(v));
                };

                // 1. Symmetric balanced: head [0, head_n) + tail [R-tail_n, R)
                {
                    const long head_n = m / 2;
                    const long tail_n = m - head_n;
                    std::vector<long> v;
                    for (long i = 0; i < head_n; ++i) v.push_back(i);
                    for (long i = R - tail_n; i < R; ++i) v.push_back(i);
                    add_sorted(std::move(v));
                }
                // 2. Reverse split (other parity)
                {
                    const long head_n = m - m / 2;
                    const long tail_n = m - head_n;
                    std::vector<long> v;
                    for (long i = 0; i < head_n; ++i) v.push_back(i);
                    for (long i = R - tail_n; i < R; ++i) v.push_back(i);
                    add_sorted(std::move(v));
                }
                // 3. Contiguous start
                {
                    std::vector<long> v(m);
                    std::iota(v.begin(), v.end(), 0L);
                    add_sorted(std::move(v));
                }
                // 4. Contiguous end
                {
                    std::vector<long> v(m);
                    std::iota(v.begin(), v.end(), R - m);
                    add_sorted(std::move(v));
                }
                // 5. Outlier-anchored: top-m residuals from least-squares linear fit
                if (R >= 2) {
                    double sx = 0, sy = 0, sxy = 0, sx2 = 0;
                    const double n = static_cast<double>(R);
                    for (long i = 0; i < R; ++i) {
                        const double x = static_cast<double>(i);
                        const double y = row_groups[i].cost;
                        sx += x; sy += y; sxy += x*y; sx2 += x*x;
                    }
                    const double denom = n * sx2 - sx * sx;
                    if (std::abs(denom) > 0.0) {
                        const double a = (n * sxy - sx * sy) / denom;
                        const double b = (sy - a * sx) / n;
                        std::vector<std::pair<long, double>> res(R);
                        for (long i = 0; i < R; ++i) {
                            res[i] = {i, std::abs(row_groups[i].cost - (a * i + b))};
                        }
                        std::sort(res.begin(), res.end(),
                                  [](const auto& p, const auto& q) { return p.second > q.second; });
                        std::vector<long> v;
                        for (long i = 0; i < m; ++i) v.push_back(res[i].first);
                        add_sorted(std::move(v));
                    }
                }
                // 6. Exhaustive if feasible: enumerate combinations of m from R
                {
                    // C(R, m) <= max_exhaustive: bounded multiplication that
                    // short-circuits if it exceeds the cap.
                    long count = 1;
                    bool fits = true;
                    for (long i = 0; i < m && fits; ++i) {
                        count *= (R - i);
                        count /= (i + 1);
                        if (count > max_exhaustive) fits = false;
                    }
                    if (fits) {
                        std::vector<long> v(m);
                        std::iota(v.begin(), v.end(), 0L);
                        // Generate all combinations via index advancement.
                        while (true) {
                            cands.insert(v);
                            // advance: find rightmost that can be incremented
                            long i = m - 1;
                            while (i >= 0 && v[i] == R - m + i) --i;
                            if (i < 0) break;
                            ++v[i];
                            for (long j = i + 1; j < m; ++j) v[j] = v[j-1] + 1;
                        }
                    }
                }
            }

            std::vector<std::vector<long>> out;
            out.reserve(cands.size());
            for (const auto& c : cands) out.push_back(c);
            return out;
        }

        /// Run the m-flex peel search: pick the peel set of size
        /// m_min = R mod (2*nsubworld) that minimizes per-rank load delta.
        /// Returns an empty peel when no improvement is possible (R < 2*nsubworld
        /// or partition empty).
        std::vector<long> find_best_peel(
                const MacroTaskPartitioner::partitionT& partition, long nsubworld) const {
            if (nsubworld <= 0) return {};
            std::vector<FoldRowGroup> row_groups = build_row_groups(partition);
            const long R = static_cast<long>(row_groups.size());
            const long two_np = 2 * nsubworld;
            if (R < two_np) return {};

            const long m_min = R % two_np;
            auto cands = generate_peel_candidates(R, m_min, row_groups, mflex_max_exhaustive_);
            if (cands.empty()) return {};

            std::vector<long> best;
            double best_delta = std::numeric_limits<double>::infinity();
            for (auto& c : cands) {
                const double d = load_delta_for_peel(row_groups, nsubworld, c);
                if (d < best_delta) {
                    best_delta = d;
                    best = c;
                }
            }
            return best;
        }

        /// Shuffle the partition list so that each owner's tasks appear in random
        /// order, reducing synchronized fetch contention across ranks.
        /// The shuffling is deterministic (seeded by nsubworld) for reproducibility.
        /// Tasks from different owners are interleaved round-robin so the queue
        /// keeps all ranks fed from the start.
        void shuffle_partition_by_owner(MacroTaskPartitioner::partitionT& partition, long nsubworld) const {
            if (!use_owner_aware_fetch() || owner_map_.empty() || nsubworld <= 0) return;

            // sym_p2p uses a cache-aware ORDER (locality-preserving, runs even with the
            // random de-sync shuffle disabled); the other owner-aware algorithms use the
            // optional random shuffle. Bail if neither applies.
            const bool cache_sort = use_sym_p2p_owner_algorithm();
            if (!cache_sort && !shuffle_task_order_) return;

            // Group partition entries by owner
            std::map<long, std::vector<std::pair<Batch,double>>> per_owner;
            for (auto& entry : partition) {
                const Batch_1D& col_range = entry.first.input[0];
                const Batch_1D& row_range = (entry.first.input.size() > 1)
                        ? entry.first.input[1] : entry.first.input[0];
                auto key = std::make_pair(col_range.begin, row_range.begin);
                auto it = owner_map_.find(key);
                long owner = (it != owner_map_.end()) ? it->second : -1;
                per_owner[owner].push_back(std::move(entry));
            }

            if (cache_sort) {
                // Cache-aware ordering (sort_for_caching). Per owner: no-fetch tasks
                // (diagonal i==j, or both operands owned) first, then group by the
                // transient (non-owned) batch index so the bounded LRU hits on
                // consecutive same-transient tasks and a future prefetch can stage
                // the next transient. begin->index from the owner_map_ keys (the split
                // is contiguous ascending, so sorted-unique begins map 1:1 to indices).
                std::set<long> begins;
                for (const auto& kv : owner_map_) {
                    begins.insert(kv.first.first);
                    begins.insert(kv.first.second);
                }
                std::map<long,long> begin_to_idx;
                { long idx = 0; for (long b : begins) begin_to_idx[b] = idx++; }
                const long n = nsubworld;
                auto idx_of = [&](const Batch& bt, int dim) -> long {
                    const Batch_1D& r = (dim == 1 && bt.input.size() > 1) ? bt.input[1] : bt.input[0];
                    auto it = begin_to_idx.find(r.begin);
                    return (it != begin_to_idx.end()) ? it->second : -1;
                };
                auto transient_of = [&](const Batch& bt, long owner) -> long {
                    const long i = idx_of(bt, 0);   // col
                    const long j = idx_of(bt, 1);   // row
                    if (i == j) return -1;                       // diagonal: zero fetch
                    const bool ci = (i % n) == owner;
                    const bool cj = (j % n) == owner;
                    if (ci && cj) return -1;                     // both-owned: zero fetch
                    return ci ? j : i;                           // the rotating transient batch
                };
                for (auto& [owner, tasks] : per_owner) {
                    std::stable_sort(tasks.begin(), tasks.end(),
                        [&](const std::pair<Batch,double>& a, const std::pair<Batch,double>& b) {
                            const long ta = transient_of(a.first, owner);
                            const long tb = transient_of(b.first, owner);
                            if (ta != tb) return ta < tb;        // -1 (no-fetch) first, then grouped
                            // deterministic tie-break within a transient group
                            return std::make_pair(idx_of(a.first,0), idx_of(a.first,1))
                                 < std::make_pair(idx_of(b.first,0), idx_of(b.first,1));
                        });
                }
            } else {
                // Shuffle each owner's task list independently (random de-sync)
                for (auto& [owner, tasks] : per_owner) {
                    std::mt19937 rng(static_cast<unsigned>(owner * 31 + nsubworld));
                    std::shuffle(tasks.begin(), tasks.end(), rng);
                }
            }

            // Rebuild partition by round-robin interleaving across owners
            partition.clear();
            bool any_remaining = true;
            std::vector<long> owner_ids;
            for (const auto& [owner, tasks] : per_owner) {
                owner_ids.push_back(owner);
            }
            std::vector<std::size_t> indices(owner_ids.size(), 0);

            while (any_remaining) {
                any_remaining = false;
                for (std::size_t g = 0; g < owner_ids.size(); ++g) {
                    const auto& tasks = per_owner[owner_ids[g]];
                    if (indices[g] < tasks.size()) {
                        partition.push_back(tasks[indices[g]]);
                        indices[g]++;
                        if (indices[g] < tasks.size()) any_remaining = true;
                    }
                }
            }
        }

        long owner_hint(const Batch& task_batch, const long nsubworld) const override {
            if (not use_owner_aware_fetch() or nsubworld<=0) return -1;
            MADNESS_CHECK_THROW(task_batch.input.size()>0, "empty task batch in owner_hint");

            // Use pre-computed fold-based assignment if available
            if (!owner_map_.empty()) {
                const Batch_1D& col_range = task_batch.input[0];
                const Batch_1D& row_range = (task_batch.input.size() > 1) ? task_batch.input[1] : task_batch.input[0];
                auto key = std::make_pair(col_range.begin, row_range.begin);
                auto it = owner_map_.find(key);
                if (it != owner_map_.end()) return it->second;
            }

            // Fallback: modulo assignment (should not be reached after prepare_owner_assignment)
            const Batch_1D& row_range = (task_batch.input.size()>1) ? task_batch.input[1] : task_batch.input[0];
            return std::max<long>(0,row_range.begin) % nsubworld;
        }

        /// the exchange task manages its own data movement via owner-aware fetch
        bool handles_own_data_movement() const override { return use_owner_aware_fetch(); }

        /// opt into own-output accumulation: each task folds its result into a
        /// subworld-local buffer; a single subworld->universe gaxpy happens in
        /// finalize_into(). Only active for the owner-aware algorithm.
        bool accumulates_own_output() const override {
            return use_owner_aware_fetch() and accumulation_mode_ >= 1;
        }

        /// opt into the framework's node-local hierarchical finalize (mode 2):
        /// stage1 reduces Kf_local_ -> node-shared Kf_node_ (intra-node), stage2
        /// does one node-sum -> universe gaxpy. Only meaningful when the task also
        /// accumulates_own_output(). The framework provides the node sub-World.
        bool wants_node_local_reduction() const {
            return use_owner_aware_fetch() and accumulation_mode_ == 2;
        }

        /// fold a task's subworld-local result vector into the subworld-local
        /// accumulator Kf_local_. Lazily initialized on first call per subworld.
        /// Typed as vecfuncT (== resultT) because resultT is declared later in
        /// the class body.
        void accumulate_locally(World& subworld, const vecfuncT& result_subworld) const {
            const long wid = long(subworld.id());
            if (not Kf_local_initialized_ or Kf_local_world_id_ != wid) {
                Kf_local_ = zero_functions_compressed<T, NDIM>(subworld, nresult);
                Kf_local_initialized_ = true;
                Kf_local_world_id_ = wid;
            }
            // match MacroTaskInternal::accumulate_into_final_result state handling
            vecfuncT& rs = const_cast<vecfuncT&>(result_subworld);
            TreeState op_state = rs[0].get_impl()->get_tensor_type()==TT_FULL
                                 ? compressed : reconstructed;
            change_tree_state(rs, op_state);
            gaxpy(1.0, Kf_local_, 1.0, rs, false);
        }

        /// single subworld->universe gaxpy, draining Kf_local_ into the
        /// universe-resident result. No-op if the subworld ran zero owned
        /// tasks (accumulator never initialized) or was already finalized.
        ///
        /// For small_memory_mt_owner the local accumulator is truncated before
        /// the universe gaxpy: each rank's Kf_local_ contains the full sum over
        /// its row contribution (sum over all k-batches), and truncation here
        /// is consistent with that being the per-row "final" result. Truncation
        /// is skipped for small_memory_symmetric_mt_owner because individual
        /// row contributions are spread across multiple ranks there and the
        /// final truncate happens once in K_macrotask_efficient on the universe Kf.
        void finalize_into(World& subworld, vecfuncT& universe_result) {
            // Once per rank (uniform across the replicated taskq => every rank's
            // first call lands on the same loop iteration). Replaces the old
            // `if (not Kf_local_initialized_) return;` data-less early-return so
            // that data-less ranks still participate in coalesced_gaxpy's
            // collective readiness barrier.
            if (finalize_universe_done_) return;
            finalize_universe_done_ = true;
            if (universe_result.empty()) return;   // nresult>0 in practice; uniform
            if (Kf_local_initialized_) {
                change_tree_state(Kf_local_, compressed);
                if (algorithm_ == small_memory_mt_owner) truncate(subworld, Kf_local_);
            }
            vecfuncT empty;
            vecfuncT& src = Kf_local_initialized_ ? Kf_local_ : empty;
            World& universe = universe_result.front().world();
            // Coalesced subworld->universe drain (option C0): bulk per-owner chunks
            // instead of one AM per tree node. Completes at the framework's universe
            // fence (so the lifetime contract below is unchanged).
            coalesced_gaxpy<T, NDIM>(universe, get_universe_reducer(universe),
                                     universe_result, src, T(1.0), finalize_chunk_entries());
            // Do NOT clear Kf_local_ here: keep the original lifetime (released at
            // cleanup() after the universe fence). (C0 copies node data into transit
            // records, so Kf_local_ is no longer referenced by in-flight tasks —
            // early release is a possible follow-up memory win; kept conservative.)
            Kf_local_initialized_ = false;
            Kf_local_world_id_ = -1;
        }

        /// Collectively (re)build the node-shared accumulator Kf_node_ in the node
        /// sub-World. MUST be called by every rank in `nodeworld` (collective
        /// Function construction) -- the framework guarantees this by driving
        /// finalize_stage1 uniformly across the identical replicated taskq, and the
        /// Kf_node_initialized_ guard flips in lockstep on all node ranks. Built with
        /// an explicit node pmap because the global FunctionDefaults pmap is the
        /// 1-rank subworld pmap during finalize -- letting zero_functions() inherit
        /// it would map keys to subworld rank indices while living in nodeworld.
        void ensure_node_accumulator(World& nodeworld) const {
            const long wid = long(nodeworld.id());
            if (Kf_node_initialized_ and Kf_node_world_id_ == wid) return;
            auto node_pmap = FunctionDefaults<NDIM>::make_default_pmap(nodeworld);
            Kf_node_.resize(nresult);
            for (long i = 0; i < nresult; ++i)
                Kf_node_[i] = functionT(FunctionFactory<T, NDIM>(nodeworld)
                                            .pmap(node_pmap)
                                            .compressed(true)
                                            .initial_level(1)
                                            .fence(false));
            // make Kf_node_ exist on every node rank before any stage1 send targets
            // it (intra-node, once per node sub-World).
            nodeworld.gop.fence();
            Kf_node_initialized_ = true;
            Kf_node_world_id_ = wid;
        }

        /// Mode-2 stage 1: reduce this rank's Kf_local_ into the node-shared Kf_node_
        /// (subworld(1-rank) -> nodeworld gaxpy; all sends intra-node). The gaxpy is
        /// fenceless and completes at the framework's node fence between the two
        /// stages, so Kf_local_ must outlive that fence (released at cleanup()).
        /// nodeworld==nullptr => mode 0/1, nothing to do here.
        void finalize_stage1(World& subworld, World* nodeworld) {
            if (nodeworld == nullptr) return;
            ensure_node_accumulator(*nodeworld);   // collective on all node ranks
            // Once per rank (uniform). Replaces the data-less early-return so
            // data-less ranks still hit coalesced_gaxpy's collective barrier.
            if (finalize_stage1_done_) return;
            finalize_stage1_done_ = true;
            if (Kf_local_initialized_) change_tree_state(Kf_local_, compressed);
            vecfuncT empty;
            vecfuncT& src = Kf_local_initialized_ ? Kf_local_ : empty;
            // Coalesced subworld->node drain (intra-node); completes at the node fence.
            coalesced_gaxpy<T, NDIM>(*nodeworld, get_node_reducer(*nodeworld),
                                     Kf_node_, src, T(1.0), finalize_chunk_entries());
            // Kf_local_ drained; do NOT clear (see finalize_into lifetime note).
            Kf_local_initialized_ = false;
            Kf_local_world_id_ = -1;
        }

        /// Mode-2 stage 2: drain the node-shared Kf_node_ into the universe result
        /// (nodeworld -> universe gaxpy; each result key has one inter-node sender
        /// per node). Falls back to the mode-0/1 direct Kf_local_ -> universe gaxpy
        /// when nodeworld==nullptr. The gaxpy is fenceless and completes at the
        /// framework's universe fence, so Kf_node_ must outlive it (released at
        /// cleanup()). Truncation matches finalize_into: applied for the asymmetric
        /// row owner (per-row-final node-sums), skipped for the symmetric algos
        /// (per-row partials that the universe-level truncate handles).
        void finalize_stage2(World& subworld, World* nodeworld, vecfuncT& universe_result) {
            if (nodeworld == nullptr) {
                finalize_into(subworld, universe_result);
                return;
            }
            // Once per rank (uniform). Kf_node_initialized_ is collectively true on
            // every node rank (ensure_node_accumulator), so this is not the data-less
            // case — but the once-guard still ensures the collective barrier fires once.
            if (finalize_universe_done_) return;
            finalize_universe_done_ = true;
            if (universe_result.empty()) return;
            if (Kf_node_initialized_) {
                change_tree_state(Kf_node_, compressed);
                if (algorithm_ == small_memory_mt_owner) truncate(*nodeworld, Kf_node_);
            }
            vecfuncT empty;
            vecfuncT& src = Kf_node_initialized_ ? Kf_node_ : empty;
            World& universe = universe_result.front().world();
            // Coalesced node->universe drain (option C0): bulk per-owner chunks
            // instead of one AM per tree node. Completes at the framework's universe fence.
            coalesced_gaxpy<T, NDIM>(universe, get_universe_reducer(universe),
                                     universe_result, src, T(1.0), finalize_chunk_entries());
            Kf_node_initialized_ = false;
            Kf_node_world_id_ = -1;
        }

        void set_next_vf_hint(const Batch_1D& next_hint, const bool has_hint) {
            if (not use_owner_aware_fetch()) return;
            // In small_memory_mt_owner the vf dimension is held, not rotated, so
            // skip the vf-side prefetch state entirely (k-batch state is rotated instead).
            // sym_p2p with cloud fetches vf from the cloud in operator(), so skip the
            // copy-from-universe vf prefetch too (it would be unused + wasteful).
            if (use_row_owner_algorithm() or use_cloud_batch_fetch()) return;
            vf_prefetch_.has_hint = has_hint;
            if (has_hint) {
                vf_prefetch_.next_hint = next_hint;
                if (vf_prefetch_.has_next and not hint_matches_range(next_hint, vf_prefetch_.next_range)) {
                    vf_prefetch_.next_data.clear();
                    vf_prefetch_.has_next = false;
                    vf_prefetch_.next_range = Batch_1D();
                }
            } else {
                vf_prefetch_.next_hint = Batch_1D();
                vf_prefetch_.next_data.clear();
                vf_prefetch_.has_next = false;
                vf_prefetch_.next_range = Batch_1D();
            }
        }

        void prefetch_next_vf_async(World& world, const vecfuncT& vf_full) const {
            if (not use_owner_aware_fetch()) return;
            if (use_row_owner_algorithm() or use_cloud_batch_fetch()) return;     // vf held (row-owner) or cloud-fetched (sym_p2p)
            const double cpu0 = process_cpu_time();
            ensure_cache_world(world);
            if (not vf_prefetch_.has_hint) return;

            const Batch_1D hint = normalize_range(vf_prefetch_.next_hint, long(vf_full.size()));
            MADNESS_CHECK_THROW(hint.begin >= 0 and hint.end >= hint.begin and hint.end <= long(vf_full.size()),
                                "prefetch_next_vf_async: invalid next vf hint range");

            if (vf_prefetch_.has_current and same_range(vf_prefetch_.current_range, hint)) return;
            if (vf_prefetch_.has_next and same_range(vf_prefetch_.next_range, hint)) return;

            vecfuncT prefetched;
            prefetched.reserve(hint.size());
            for (long global_index = hint.begin; global_index < hint.end; ++global_index) {
                prefetched.push_back(copy(world, vf_full[global_index], false));
            }
            // no fence here: overlap prefetch with current task compute
            vf_prefetch_.next_range = hint;
            vf_prefetch_.next_data = std::move(prefetched);
            vf_prefetch_.has_next = true;
            const double cpu1 = process_cpu_time();
            add_owner_fetch_time(cpu0, cpu1);
        }

        const vecfuncT& acquire_current_vf(World& world, const vecfuncT& vf_batch, const Batch_1D& vf_range) const {
            ensure_cache_world(world);
            const Batch_1D normalized = normalize_range(vf_range, long(vf_batch.size()));

            if (vf_prefetch_.has_current and same_range(vf_prefetch_.current_range, normalized)) {
                return vf_prefetch_.current_data;
            }

            if (vf_prefetch_.has_next and same_range(vf_prefetch_.next_range, normalized)) {
                // Do not fence here; keep overlap with currently outstanding prefetches.
                // MADNESS function kernels will synchronize as needed when data is touched.
                vf_prefetch_.current_range = vf_prefetch_.next_range;
                vf_prefetch_.current_data = std::move(vf_prefetch_.next_data);
                vf_prefetch_.has_current = true;
                vf_prefetch_.next_data.clear();
                vf_prefetch_.has_next = false;
                vf_prefetch_.next_range = Batch_1D();
                return vf_prefetch_.current_data;
            }

            vf_prefetch_.current_data = fetch_batch_transient(world, vf_batch);
            vf_prefetch_.current_range = normalized;
            vf_prefetch_.has_current = true;
            return vf_prefetch_.current_data;
        }

        void release_finished_vf(const Batch_1D& vf_range) const {
            if (not use_owner_aware_fetch()) return;
            if (vf_prefetch_.has_current and same_range(vf_prefetch_.current_range, vf_range)) {
                const bool keep_current_for_next = vf_prefetch_.has_hint and hint_matches_range(vf_prefetch_.next_hint, vf_range);
                if (not keep_current_for_next) {
                    vf_prefetch_.current_data.clear();
                    vf_prefetch_.has_current = false;
                    vf_prefetch_.current_range = Batch_1D();
                }
            }

            if (vf_prefetch_.has_next) {
                const bool keep_next = vf_prefetch_.has_hint and hint_matches_range(vf_prefetch_.next_hint, vf_prefetch_.next_range);
                if (not keep_next) {
                    vf_prefetch_.next_data.clear();
                    vf_prefetch_.has_next = false;
                    vf_prefetch_.next_range = Batch_1D();
                }
            }
        }

        // ---- small_memory_mt_owner: rotating k-batch (bra+ket) prefetch ----
        // Mirror of the vf-side hint/prefetch/acquire/release machinery, but for the
        // k-batch dimension (mo_bra + mo_ket paired by inner-sum index). The framework
        // calls set_next_bra_hint() with the next owned task's input[1] (= bra_range),
        // then prefetch_next_bra_async() with args<1>/<2> (= mo_bra / mo_ket) so the
        // next k-batch is staged while the current task computes.

        void set_next_bra_hint(const Batch_1D& next_hint, const bool has_hint) {
            if (not use_row_owner_algorithm()) return;
            kbatch_prefetch_.has_hint = has_hint;
            if (has_hint) {
                kbatch_prefetch_.next_hint = next_hint;
                if (kbatch_prefetch_.has_next
                    and not hint_matches_range(next_hint, kbatch_prefetch_.next_range)) {
                    kbatch_prefetch_.next_bra.clear();
                    kbatch_prefetch_.next_ket.clear();
                    kbatch_prefetch_.has_next = false;
                    kbatch_prefetch_.next_range = Batch_1D();
                }
            } else {
                kbatch_prefetch_.next_hint = Batch_1D();
                kbatch_prefetch_.next_bra.clear();
                kbatch_prefetch_.next_ket.clear();
                kbatch_prefetch_.has_next = false;
                kbatch_prefetch_.next_range = Batch_1D();
            }
        }

        /// Row-owner cloud pipeline (Design B), called PRE-compute by the framework
        /// (after the argtuple load, before operator()). Two responsibilities:
        ///   1. promote: the prefetch issued during the previous task (`next`)
        ///      becomes `current`, ready to be consumed by this task's operator().
        ///   2. issue ahead: launch the async byte fetch for the *next* owned
        ///      task's bra/ket records into `next`, so the round-trip overlaps
        ///      THIS task's compute. Consumed (promoted+deserialized) one task later.
        /// Only the trigger AMs go on the wire here; the heavy transfer lands on
        /// the owner's worker task (BatchTransport). No-op unless this is the
        /// row-owner cloud path.
        void cloud_kbatch_pipeline_advance(World& world,
                                           const vecfuncT& mo_bra_full,
                                           const vecfuncT& mo_ket_full,
                                           const Batch_1D& next_hint,
                                           const bool has_next_task) const {
            // Design-B prefetch is asymmetric-only for now; the sym_p2p path runs
            // with prefetch OFF (synchronous fetch in operator()) in the first cut.
            if (not use_cloud_batch_fetch() or not use_row_owner_algorithm()) return;
            // Establish the cache world FIRST. Each SCF iteration uses a fresh
            // subworld id, so the first task triggers clear_local_caches (which
            // clears cloud_kb_). Doing it here, before we issue, means the
            // operator()'s later ensure_cache_world is a no-op and won't wipe the
            // prefetch we are about to launch (otherwise the next task misses).
            ensure_cache_world(world);
            // promote previous `next` -> `current` (the slot this task consumes)
            cloud_kb_.current = std::move(cloud_kb_.next);
            cloud_kb_.next = CloudKBatchSlot();
            if (not has_next_task) return;        // last owned task: nothing to prefetch
            if (cloud_ptr == nullptr) return;
            MADNESS_CHECK_THROW(mo_bra_full.size() == mo_ket_full.size(),
                                "cloud_kbatch_pipeline_advance: bra/ket size mismatch");
            const Batch_1D hint = normalize_range(next_hint, long(mo_bra_full.size()));
            const typename Cloud::keyT salt = batch_salt(mo_ket_full);  // full ket
            CloudKBatchSlot& s = cloud_kb_.next;
            s.range   = hint;
            s.bra_key = batch_record_key(salt, DIM_BRA, hint);
            s.ket_key = batch_record_key(salt, DIM_KET, hint);
            s.bra_fut = cloud_ptr->request_batch_bytes_async(s.bra_key);
            s.ket_fut = cloud_ptr->request_batch_bytes_async(s.ket_key);
            s.valid   = true;
            if (cloud_ptr->is_batch_debug()) {
                const ProcessID bo = cloud_ptr->batch_owner(s.bra_key);
                const ProcessID ko = cloud_ptr->batch_owner(s.ket_key);
                print("BATCH_PREFETCH rank=", world.rank(), " sw=", world.id(),
                      " t=", wall_time(), " next_row=[", hint.begin, ",", hint.end, ")",
                      " bra_owner=", bo, " ket_owner=", ko);
            }
        }

        /// Design-B prefetch (Option A) for sym_p2p. Called PRE-compute by the
        /// framework (after the argtuple load, before operator()) with the NEXT owned
        /// task's col (input[0]) and row (input[1]) ranges. For each distinct record
        /// key NOT already resident in sym_lru_ and NOT already pending, issue an async
        /// byte request so the transfer overlaps THIS task's compute. Option A: the
        /// owned batch is resident in the LRU and thus skipped, so only the transient
        /// is prefetched — no subworld→universe rank math. Diagonal / both-owned next
        /// tasks prefetch nothing (both keys resident). No-op unless this is the sym_p2p
        /// cloud path. salt = ψ (full ket) matches store_batches / operator().
        /// MUST be public: the framework's external SFINAE dispatch helper calls it.
        void sym_pipeline_advance(World& world, const vecfuncT& psi_full,
                                  const Batch_1D& next_col, const Batch_1D& next_row,
                                  const bool has_next_task) const {
            if (not use_sym_p2p_owner_algorithm() or not use_cloud_batch_fetch()) return;
            // Establish the cache world FIRST so operator()'s later ensure_cache_world
            // is a no-op and won't wipe the prefetch we issue here (mirrors
            // cloud_kbatch_pipeline_advance).
            ensure_cache_world(world);
            // Promote the previous task's `next` to `current` (consumed by this task's
            // operator()); retire whatever was in `current` (consumed already, or a
            // misprediction whose in-flight request still completes harmlessly in
            // BatchTransport — dropping our Future handle does not orphan the MPI op).
            sym_pf_current_ = std::move(sym_pf_next_);
            sym_pf_next_ = SymPrefetchSlot();
            if (not has_next_task or cloud_ptr == nullptr) return;
            const typename Cloud::keyT salt = batch_salt(psi_full);
            const Batch_1D ncol = normalize_range(next_col, long(psi_full.size()));
            const Batch_1D nrow = normalize_range(next_row, long(psi_full.size()));
            // Issue AT MOST ONE request — the next task's transient (its first
            // NON-resident operand; Option A: the owned operand is already in the LRU
            // and is skipped). Strict single-issue keeps in-flight ≤ 2 (current+next),
            // the bound the proven asymmetric path uses.
            auto issue = [&](const Batch_1D& r) -> bool {
                const typename Cloud::keyT key = batch_record_key(salt, DIM_VF, r);
                if (sym_lru_contains(key)) return false;                 // owned/resident -> skip (Option A)
                if (sym_pf_current_.valid and sym_pf_current_.key == key) return false;  // already staged
                sym_pf_next_.valid = true;
                sym_pf_next_.key = key;
                sym_pf_next_.fut = cloud_ptr->request_batch_bytes_async(key);
                if (cloud_ptr->is_batch_debug())
                    print("SYM_PREFETCH rank=", world.rank(), " sw=", world.id(), " t=", wall_time(),
                          " row=[", r.begin, ",", r.end, ") key=", key);
                return true;
            };
            if (not issue(ncol) and not same_range(nrow, ncol)) issue(nrow);  // one transient only
        }

        void prefetch_next_bra_async(World& world,
                                     const vecfuncT& mo_bra_full,
                                     const vecfuncT& mo_ket_full) const {
            if (not use_row_owner_algorithm()) return;
            // The cloud path prefetches pre-compute via cloud_kbatch_pipeline_advance
            // (Design B double buffer); the tail call is a no-op for it.
            if (use_cloud_batch_fetch()) return;

            const double cpu0 = process_cpu_time();
            ensure_cache_world(world);
            if (not kbatch_prefetch_.has_hint) return;

            MADNESS_CHECK_THROW(mo_bra_full.size() == mo_ket_full.size(),
                                "prefetch_next_bra_async: bra/ket size mismatch");
            const Batch_1D hint = normalize_range(kbatch_prefetch_.next_hint, long(mo_bra_full.size()));
            MADNESS_CHECK_THROW(hint.begin >= 0 and hint.end >= hint.begin
                                and hint.end <= long(mo_bra_full.size()),
                                "prefetch_next_bra_async: invalid next k-batch hint range");

            if (kbatch_prefetch_.has_current and same_range(kbatch_prefetch_.current_range, hint)) return;
            if (kbatch_prefetch_.has_next    and same_range(kbatch_prefetch_.next_range,    hint)) return;

            vecfuncT pb, pk;
            pb.reserve(hint.size());
            pk.reserve(hint.size());
            for (long g = hint.begin; g < hint.end; ++g) {
                pb.push_back(copy(world, mo_bra_full[g], false));
                pk.push_back(copy(world, mo_ket_full[g], false));
            }
            // no fence: overlap with current task compute
            kbatch_prefetch_.next_range = hint;
            kbatch_prefetch_.next_bra   = std::move(pb);
            kbatch_prefetch_.next_ket   = std::move(pk);
            kbatch_prefetch_.has_next   = true;
            const double cpu1 = process_cpu_time();
            add_owner_fetch_time(cpu0, cpu1);
        }

        /// Return pointers to the current k-batch's bra and ket vectors.
        ///
        /// Tries prefetch hit / next-promote first; falls back to a synchronous
        /// cold-path fetch. Mirrors acquire_current_vf in shape but is asymmetric
        /// in its sync-path inputs: when called from operator(), `bra_batch_local`
        /// is the BATCHED mo_bra slice (size == k_range.size(), local indices)
        /// while `mo_ket_full` is the UNBATCHED universe-resident mo_ket vector
        /// (size == nresult, global indices). The prefetch path is fed by the
        /// SFINAE hook with both args as unbatched (see prefetch_next_bra_async),
        /// so a prefetch hit short-circuits before the cold path reads either
        /// argument.
        std::pair<const vecfuncT*, const vecfuncT*>
        acquire_current_kbatch(World& world,
                               const vecfuncT& bra_batch_local,
                               const vecfuncT& mo_ket_full,
                               const Batch_1D& k_range) const {
            ensure_cache_world(world);
            // For the prefetch lookup we need a normalized range; we don't yet know
            // whether the cold path will run, so don't dereference inputs to size them.
            // k_range is already a partition-produced range with end > 0.
            const Batch_1D r = normalize_range(k_range, long(mo_ket_full.size()));

            if (kbatch_prefetch_.has_current and same_range(kbatch_prefetch_.current_range, r)) {
                return {&kbatch_prefetch_.current_bra, &kbatch_prefetch_.current_ket};
            }
            if (kbatch_prefetch_.has_next and same_range(kbatch_prefetch_.next_range, r)) {
                // promote next -> current. Do not fence here; let downstream
                // function-kernel calls synchronize as needed.
                kbatch_prefetch_.current_range = kbatch_prefetch_.next_range;
                kbatch_prefetch_.current_bra   = std::move(kbatch_prefetch_.next_bra);
                kbatch_prefetch_.current_ket   = std::move(kbatch_prefetch_.next_ket);
                kbatch_prefetch_.has_current = true;
                kbatch_prefetch_.next_bra.clear();
                kbatch_prefetch_.next_ket.clear();
                kbatch_prefetch_.has_next = false;
                kbatch_prefetch_.next_range = Batch_1D();
                return {&kbatch_prefetch_.current_bra, &kbatch_prefetch_.current_ket};
            }

            // Cold path: synchronous fetch.
            // bra side: range-for over the batched local vector (size matches r.size()).
            // ket side: index into the full vector by global index.
            const double cpu0 = process_cpu_time();
            MADNESS_CHECK_THROW(long(bra_batch_local.size()) == r.size(),
                                "acquire_current_kbatch cold path: bra batch size mismatch with k-range");
            MADNESS_CHECK_THROW(r.begin >= 0 and r.end <= long(mo_ket_full.size()),
                                "acquire_current_kbatch cold path: k-range out of mo_ket bounds");
            vecfuncT pb, pk;
            pb.reserve(r.size());
            pk.reserve(r.size());
            for (long local = 0; local < r.size(); ++local) {
                pb.push_back(copy(world, bra_batch_local[local], false));
                pk.push_back(copy(world, mo_ket_full[r.begin + local], false));
            }
            world.gop.fence();
            const double cpu1 = process_cpu_time();
            add_owner_fetch_time(cpu0, cpu1);
            kbatch_prefetch_.current_range = r;
            kbatch_prefetch_.current_bra   = std::move(pb);
            kbatch_prefetch_.current_ket   = std::move(pk);
            kbatch_prefetch_.has_current   = true;
            return {&kbatch_prefetch_.current_bra, &kbatch_prefetch_.current_ket};
        }

        void release_finished_kbatch(const Batch_1D& k_range) const {
            if (not use_row_owner_algorithm()) return;
            if (kbatch_prefetch_.has_current and same_range(kbatch_prefetch_.current_range, k_range)) {
                const bool keep_current_for_next = kbatch_prefetch_.has_hint
                        and hint_matches_range(kbatch_prefetch_.next_hint, k_range);
                if (not keep_current_for_next) {
                    kbatch_prefetch_.current_bra.clear();
                    kbatch_prefetch_.current_ket.clear();
                    kbatch_prefetch_.has_current = false;
                    kbatch_prefetch_.current_range = Batch_1D();
                }
            }
            if (kbatch_prefetch_.has_next) {
                const bool keep_next = kbatch_prefetch_.has_hint
                        and hint_matches_range(kbatch_prefetch_.next_hint, kbatch_prefetch_.next_range);
                if (not keep_next) {
                    kbatch_prefetch_.next_bra.clear();
                    kbatch_prefetch_.next_ket.clear();
                    kbatch_prefetch_.has_next = false;
                    kbatch_prefetch_.next_range = Batch_1D();
                }
            }
        }

        void cleanup() override {
            clear_local_caches();
            Kf_local_.clear();
            Kf_local_initialized_ = false;
            Kf_local_world_id_ = -1;
            // Kf_node_ released here (after the protective universe fence), same
            // lifetime contract as Kf_local_: its stage2 gaxpy is fenceless.
            Kf_node_.clear();
            Kf_node_initialized_ = false;
            Kf_node_world_id_ = -1;
            cache_world_id_ = -1;
            // reset the coalesced-finalize once-guards for the next K() invocation.
            finalize_stage1_done_ = false;
            finalize_universe_done_ = false;
            // Destroy the cached reducers HERE, while their transport worlds are
            // still alive (cleanup() runs in run_all's cleanup_after_run, before the
            // taskq — and thus nodeworld_ptr — is destroyed). The node sub-World is
            // recreated every K() invocation, so a reducer cached across invocations
            // would outlive its world and segfault in ~WorldObject -> unregister_ptr
            // when later reassigned. Reconstructed fresh (collectively) on next use.
            node_reducer_.reset();
            node_reducer_world_id_ = -1;
            universe_reducer_.reset();
            universe_reducer_world_id_ = -1;
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

        // ---- StoreFunctionBatched: owner-pinned cloud batches (small_memory_mt_owner) ----
        // The input vectors are serialized into the cloud's batch container as
        // owner-pinned batches whose ranges match the row-owner partition
        // (exchange_row_owner_split), so each task fetches exactly the batch it
        // needs owner-to-owner. Records are keyed by (salt, dimension, range) so the
        // store side (store_batches) and the fetch side (operator(), Step 5) derive
        // identical record keys without a separately stored manifest; owner = split
        // position (the rank<->batch bijection), recorded in the cloud's
        // CloudOwnerPmap by store_batch.
        enum BatchDim { DIM_VF = 0, DIM_BRA = 1, DIM_KET = 2 };

        /// per-invocation salt from the ket identities; identical on every rank.
        /// Uses the ket vector because it is available in full both here (store
        /// side, mo_ket) and in operator() (the un-batched 3rd argument vket), so
        /// both sides derive the same record keys.
        static typename Cloud::keyT batch_salt(const vecfuncT& ket) {
            std::size_t k = 0x5a17ull;
            for (const auto& f : ket) hash_combine(k, hash_value(f.get_impl()->id()));
            return typename Cloud::keyT(k);
        }

        /// deterministic record key for a batch: (salt, dimension, range).
        /// range-keyed so store and fetch agree without communicating a manifest.
        static typename Cloud::keyT batch_record_key(const typename Cloud::keyT salt,
                                                     const int dim, const Batch_1D& r) {
            std::size_t k = std::size_t(salt);
            hash_combine(k, std::size_t(dim));
            hash_combine(k, std::size_t(r.begin));
            hash_combine(k, std::size_t(r.end));
            return typename Cloud::keyT(k);
        }

        /// store the inputs as owner-pinned batches aligned with the row-owner
        /// partition. Collective on `world` (the universe); no-op unless the
        /// cloud-batch fetch path is active. Called by the framework via the
        /// store_batches_or_noop hook right after the pointer argtuple is stored.
        void store_batches(World& world, World& subworld, Cloud& cloud, const argtupleT& argtuple,
                           const long nsubworld) {
            if (not use_cloud_batch_fetch()) return;
            // Phase 4a-i: store_batch now skips the per-function internal gop.fences,
            // so guarantee the input functions are complete with ONE fence up front
            // (replaces the ~480 per-function fences the gather used to do).
            world.gop.fence();
            const vecfuncT& vf = std::get<0>(argtuple);      // input[0]: held col / vf
            const vecfuncT& mo_bra = std::get<1>(argtuple);  // input[1]: rotating row / bra
            const vecfuncT& mo_ket = std::get<2>(argtuple);  // paired with bra by inner index
            const typename Cloud::keyT salt = batch_salt(mo_ket);  // ket: full in both store and operator()

            // Symmetric one-set store: bra==ket==vf==psi (guaranteed: nemo with
            // bra!=ket falls back to copy in K_macrotask_efficient). Store psi in
            // M=level*nproc batches over a single split; batch k pinned to
            // rank (k mod nsubworld). The symmetric kernels read this same set for
            // the column (vf), row (bra), and ket roles. One record/batch (DIM_VF).
            if (use_sym_p2p_owner_algorithm()) {
                // DISTRIBUTED per-owner store (B2-via-copy). The old path serialized
                // every batch via a collective ParallelOutputArchive gather to rank 0
                // (96 sequential MPI_Gathervs → rank 0 ingests the whole ~24 GB set
                // through one NIC). Instead: every rank registers routing for ALL
                // records (local, no comm), then each rank distributively pulls the
                // batches IT owns into its own size-1 subworld (copy() → owner pulls
                // p2p from every source rank, no rank-0 funnel) and serializes them
                // there. All owners run concurrently → the ingest spreads across owners
                // (~0.5 GB each) instead of funneling. Records are semantically
                // identical (the size-1-subworld store emits the same parallel-format
                // bytes; pair order differs but deserialize is order-independent), so
                // the fetch path is untouched.
                static const bool phase_dbg = (std::getenv("MAD_BATCH_STORE_PHASES") != nullptr);
                const std::vector<Batch_1D> split =
                        exchange_sym_owner_split(mo_ket.size(), nsubworld, batch_granularity_);
                const double t0 = wall_time();
                // OPTION-2-STEP-1 PROBE: replicate run_all's context for the cross-world
                // copy by making the GLOBAL FunctionDefaults pmap the size-1 subworld
                // (in setup it is the universe / load-balanced pmap; if any code on the
                // copy path reads the global pmap it would route to universe ranks that
                // do not exist in the size-1 subworld -> MPI_ERR_RANK). Restored below.
                auto b2_saved_pmap = FunctionDefaults<NDIM>::get_pmap();
                FunctionDefaults<NDIM>::set_default_pmap(subworld);
                // phase 0+1: register routing on every rank; owner pulls its batches.
                std::vector<std::pair<typename Cloud::keyT, vecfuncT>> pending;
                for (long k = 0; k < long(split.size()); ++k) {
                    const Batch_1D& r = split[k];
                    const ProcessID owner = ProcessID(k % nsubworld);
                    const typename Cloud::keyT record = batch_record_key(salt, DIM_VF, r);
                    cloud.register_batch_owner(record, owner);   // all ranks, local
                    if (owner == world.rank()) {
                        vecfuncT local(r.end - r.begin);
                        for (long i = r.begin; i < r.end; ++i)
                            local[i - r.begin] = copy(subworld, mo_ket[i], /*fence=*/false);
                        pending.emplace_back(record, std::move(local));
                    }
                }
                const double t_issue = wall_time();
                // Complete the cross-world copies on the SUBWORLD (size-1), matching the
                // copy-fetch idiom: the source ranks' comm threads serve the
                // serialize_remote_coeffs requests, so a subworld fence drains this
                // owner's pending inserts. (The global-pmap wrap above is what keeps the
                // cross-world copy routing inside the size-1 subworld.)
                subworld.gop.fence();
                const double t_copy = wall_time();
                // phase 2: each owner serializes its now-local batches into the record.
                // store_batch over the size-1 subworld → trivial-local gather → same
                // parallel-format bytes; its internal set_owner re-registers harmlessly.
                for (auto& pr : pending)
                    cloud.store_batch(subworld, pr.second, world.rank(), pr.first, /*fence=*/false);
                subworld.gop.fence();   // complete the (local) record writes
                FunctionDefaults<NDIM>::set_pmap(b2_saved_pmap);   // restore the universe pmap
                world.gop.fence();      // collective alignment; no cross-world ops pending now
                if (phase_dbg && world.rank() == 0) {
                    printf("BATCH_STORE_PHASES(B2) nbatch=%ld owned=%ld issue=%.2f copy=%.2f store=%.2f\n",
                           long(split.size()), long(pending.size()),
                           t_issue - t0, t_copy - t_issue, wall_time() - t_copy);
                }
                return;
            }

            const std::vector<Batch_1D> col_split = exchange_row_owner_split(vf.size(), nsubworld);
            const std::vector<Batch_1D> row_split = exchange_row_owner_split(mo_bra.size(), nsubworld);

            // vf col-batches: owner == split index (held locally by that rank).
            // Suppress per-call fences inside store_batch; emit one collective
            // fence after all three loops to amortize the ~3*nsubworld fences
            // that would otherwise dominate setup time.
            for (long i = 0; i < long(col_split.size()); ++i) {
                const Batch_1D& r = col_split[i];
                vecfuncT slice(vf.begin() + r.begin, vf.begin() + r.end);
                cloud.store_batch(world, slice, ProcessID(i), batch_record_key(salt, DIM_VF, r),
                                  /*fence=*/false);
            }
            // bra/ket row-batches: separate records, owner == split index
            for (long j = 0; j < long(row_split.size()); ++j) {
                const Batch_1D& r = row_split[j];
                vecfuncT bra_slice(mo_bra.begin() + r.begin, mo_bra.begin() + r.end);
                vecfuncT ket_slice(mo_ket.begin() + r.begin, mo_ket.begin() + r.end);
                cloud.store_batch(world, bra_slice, ProcessID(j), batch_record_key(salt, DIM_BRA, r),
                                  /*fence=*/false);
                cloud.store_batch(world, ket_slice, ProcessID(j), batch_record_key(salt, DIM_KET, r),
                                  /*fence=*/false);
            }
            world.gop.fence();  // single collective fence covering all batch writes above
        }

        /// Dump owner_map_ and (for cloud-batch) the CloudOwnerPmap batch
        /// routing at universe rank 0. Called by the framework via the
        /// log_owner_layout_or_noop hook right after store_batches; no-op
        /// unless log_diagnostics_ is set (by K_macrotask_efficient when
        /// printdebug() is on).
        void log_owner_layout(World& universe, Cloud& cloud, const long nsubworld) const {
            if (not log_diagnostics_) return;
            if (universe.rank() != 0) return;
            print("===== owner layout for", name, "=====");
            print("  algorithm=", int(algorithm_), " nsubworld=", nsubworld,
                  " use_cloud_batch_fetch=", int(use_cloud_batch_fetch()),
                  " accumulation_mode=", accumulation_mode_,
                  " ntasks=", owner_map_.size());
            print("  ---- task -> rank (owner_map_) ----");
            for (const auto& kv : owner_map_) {
                print("    (col_begin=", kv.first.first,
                      ", row_begin=", kv.first.second, ") -> rank", kv.second);
            }
            if (use_cloud_batch_fetch()) {
                print("  ---- batch record routing (CloudOwnerPmap) ----");
                cloud.print_batch_owner_map(universe, name);
            }
            print("===== end owner layout =====");
        }

        std::vector<Function<T, NDIM>>
        operator()(const std::vector<Function<T, NDIM>>& vf_batch,     // will be batched (column)
                   const std::vector<Function<T, NDIM>>& bra_batch,    // will be batched (row)
                   const std::vector<Function<T, NDIM>>& vket) {       // will not be batched

            MADNESS_CHECK_THROW(subworld_ptr!=0, "MacroTaskExchangeSimple: subworld_ptr is null");
            World& world = *subworld_ptr;
            resultT Kf = zero_functions_compressed<T, NDIM>(world, nresult);

            bool diagonal_block = batch.input[0] == batch.input[1];
            auto& bra_range = batch.input[1];    // corresponds to vbra
            auto& vf_range = batch.input[0];       // corresponds to vf_batch

            if (vf_range.is_full_size()) vf_range.end = vf_batch.size();
            if (bra_range.is_full_size()) bra_range.end = bra_batch.size();

            MADNESS_CHECK(vf_range.end <= nresult);
            if (symmetric) MADNESS_CHECK(bra_range.end <= nresult);

            const double t_op_start = wall_time();

            // per-task profiling (MAD_EXCH_TASK_PROFILE): reset + stamp identity/geometry.
            // Interior timers (wait/compute/components/cpu) are filled by the accumulator
            // helpers and the smallmem kernels; wall_end + peak RSS + emit happen at exit.
            const bool prof_on = exch_task_profile_enabled();
            if (prof_on) {
                prof_.reset();
                prof_.task_id = prof_task_seq_++;
                prof_.universe_rank = universe_rank_;
                prof_.subworld_id = static_cast<unsigned long>(world.id());
                prof_.subworld_nrank = world.size();
                prof_.thresh = FunctionDefaults<NDIM>::get_thresh();   // current protocol tier
                prof_.k = FunctionDefaults<NDIM>::get_k();
                prof_.diagonal = diagonal_block;
                prof_.row_begin = bra_range.begin; prof_.row_end = bra_range.end;
                prof_.col_begin = vf_range.begin;  prof_.col_end = vf_range.end;
                prof_.n_row = bra_range.end - bra_range.begin;
                prof_.n_col = vf_range.end - vf_range.begin;
                prof_.wall_start = t_op_start;
            }

            // -------- small_memory_mt_owner: row-owner asymmetric branch --------
            // For this algorithm each rank exclusively owns one i-batch (= vf_range).
            // The held i-batch is fetched once per subworld into held_vf_cache_ and
            // reused across all tasks. The k-batch (= mo_bra + mo_ket paired by inner
            // index = bra_range) rotates task-by-task and is prefetched by the
            // framework's set_next_bra_hint / prefetch_next_bra_async hooks.
            if (use_row_owner_algorithm()) {

                // ---- cloud-batch fetch path ----
                // Fetch the input batches from the cloud's owner-pinned batch
                // container instead of copy()ing from the universe. vf is owned by
                // this rank (local read); the bra/ket k-batch is owned by another
                // rank (remote read). Records are keyed by (salt, dim, range),
                // matching store_batches; owner routing is in the CloudOwnerPmap.
                // The k-batch is consumed from the prefetch futures issued by the
                // previous task (Step 5b) when available; else fetched synchronously.
                if (use_cloud_batch_fetch()) {
                    MADNESS_CHECK_THROW(cloud_ptr != nullptr, "cloud_ptr is null in cloud-batch fetch");
                    // reset stale static prefetch state if the subworld changed
                    // (clears cloud_kb_ futures from a previous application/cloud)
                    ensure_cache_world(world);
                    Cloud& cloud = *cloud_ptr;
                    const bool op_dbg = cloud.is_batch_debug();
                    if (op_dbg) {
                        fprintf(stderr, "FIXB_DBG op() ENTER sw_id=%lu col=[%ld,%ld) row=[%ld,%ld) cur_valid=%d t=%.6f\n",
                                static_cast<unsigned long>(world.id()),
                                vf_range.begin, vf_range.end, bra_range.begin, bra_range.end,
                                int(cloud_kb_.current.valid), wall_time());
                        fflush(stderr);
                    }
                    const typename Cloud::keyT salt = batch_salt(vket);  // full ket
                    const double cf0 = process_cpu_time();

                    // vf is held: dedup across this rank's whole task sequence (one
                    // fetch + many reuses), so we keep it in the cloud cache.
                    if (op_dbg) { fprintf(stderr, "FIXB_DBG op() vf_fetch BEGIN t=%.6f\n", wall_time()); fflush(stderr); }
                    vecfuncT vf_local = cloud.fetch_batch_p2p<T, NDIM>(
                            world, batch_record_key(salt, DIM_VF, vf_range),
                            /*cache_result=*/true);                              // local (owner==self)
                    const double t_vf_done_c = wall_time();
                    if (op_dbg) { fprintf(stderr, "FIXB_DBG op() vf_fetch DONE t=%.6f\n", wall_time()); fflush(stderr); }

                    // The k-batch rotates: each row-batch is fetched exactly once
                    // per rank in pure rotation, so caching has no dedup benefit
                    // and only grows memory. cache_result=false releases it once
                    // operator() returns and the local vecfuncTs go out of scope,
                    // matching the copy() path's release_finished_kbatch behavior.
                    const typename Cloud::keyT bra_key = batch_record_key(salt, DIM_BRA, bra_range);
                    const typename Cloud::keyT ket_key = batch_record_key(salt, DIM_KET, bra_range);
                    // `current` was filled by the previous task's PRE-compute
                    // cloud_kbatch_pipeline_advance (Design B): the trigger AMs went
                    // on the wire one task earlier and the transfer overlapped the
                    // previous task's compute, so it should be ready here.
                    const bool cloud_prefetch_hit = cloud_kb_.current.valid
                        and same_range(cloud_kb_.current.range, normalize_range(bra_range, long(bra_batch.size())))
                        and cloud_kb_.current.bra_key == bra_key
                        and cloud_kb_.current.ket_key == ket_key;
                    vecfuncT bra_k, ket_k;
                    if (cloud_prefetch_hit) {
                        // consume the prefetch issued (and overlapped) one task ago
                        cloud.tally_batch_prefetch_hit();
                        if (op_dbg) { fprintf(stderr, "FIXB_DBG op() consume_hit BEFORE_BRA t=%.6f\n", wall_time()); fflush(stderr); }
                        bra_k = cloud.deserialize_batch_p2p<T, NDIM>(world, cloud_kb_.current.bra_fut, bra_key,
                                                                 /*cache_result=*/false);
                        if (op_dbg) { fprintf(stderr, "FIXB_DBG op() consume_hit AFTER_BRA  t=%.6f\n", wall_time()); fflush(stderr); }
                        ket_k = cloud.deserialize_batch_p2p<T, NDIM>(world, cloud_kb_.current.ket_fut, ket_key,
                                                                 /*cache_result=*/false);
                        if (op_dbg) { fprintf(stderr, "FIXB_DBG op() consume_hit AFTER_KET  t=%.6f\n", wall_time()); fflush(stderr); }
                        cloud_kb_.current.valid = false;
                    } else {
                        cloud.tally_batch_prefetch_miss();
                        if (op_dbg) { fprintf(stderr, "FIXB_DBG op() sync_miss BEFORE_BRA t=%.6f\n", wall_time()); fflush(stderr); }
                        bra_k = cloud.fetch_batch_p2p<T, NDIM>(world, bra_key,
                                                           /*cache_result=*/false);   // remote, synchronous (p2p)
                        if (op_dbg) { fprintf(stderr, "FIXB_DBG op() sync_miss AFTER_BRA  t=%.6f\n", wall_time()); fflush(stderr); }
                        ket_k = cloud.fetch_batch_p2p<T, NDIM>(world, ket_key,
                                                           /*cache_result=*/false);   // remote, synchronous (p2p)
                        if (op_dbg) { fprintf(stderr, "FIXB_DBG op() sync_miss AFTER_KET  t=%.6f\n", wall_time()); fflush(stderr); }
                    }
                    const double t_k_done_c = wall_time();
                    add_owner_fetch_time(cf0, process_cpu_time());

                    if (op_dbg) { fprintf(stderr, "FIXB_DBG op() compute BEGIN t=%.6f\n", wall_time()); fflush(stderr); }
                    vecfuncT resultcolumn = compute_batch_in_asymmetric_matrix_smallmem(
                            world, ket_k, bra_k, vf_local);
                    const double t_compute_end_c = wall_time();
                    if (op_dbg) { fprintf(stderr, "FIXB_DBG op() compute DONE  t=%.6f\n", wall_time()); fflush(stderr); }
                    add_owner_compute_wall_time(t_k_done_c, t_compute_end_c);

                    for (int i = vf_range.begin; i < vf_range.end; ++i) {
                        Kf[i] += resultcolumn[i - vf_range.begin];
                    }
                    if (op_dbg) { fprintf(stderr, "FIXB_DBG op() RETURN t=%.6f\n", wall_time()); fflush(stderr); }
                    if (world.rank() == 0) {
                        print("OVERLAP_OP row_owner cloud_batch prefetch_hit=", int(cloud_prefetch_hit),
                              " col=[", vf_range.begin, ",", vf_range.end,
                              ") row=[", bra_range.begin, ",", bra_range.end,
                              ") vf_fetch=", t_vf_done_c - t_op_start,
                              " k_fetch=", t_k_done_c - t_vf_done_c,
                              " compute=", t_compute_end_c - t_k_done_c);
                    }
                    return Kf;
                }
                // ---- end cloud-batch fetch path ----

                const bool kbatch_prefetch_hit = kbatch_prefetch_.has_next
                    and same_range(kbatch_prefetch_.next_range, normalize_range(bra_range, long(bra_batch.size())));

                // Held i-batch (vf): one-shot fetch into held_vf_cache_; cached by global vf index.
                vecfuncT vf_local = fetch_batch_with_cache(world, vf_batch, vf_range, held_vf_cache_);
                const double t_vf_done = wall_time();

                // Rotating k-batch: try prefetch state, else synchronous fetch from universe.
                auto [bra_k_ptr, ket_k_ptr] = acquire_current_kbatch(world, bra_batch, vket, bra_range);
                const double t_k_done = wall_time();

                vecfuncT resultcolumn = compute_batch_in_asymmetric_matrix_smallmem(
                        world, *ket_k_ptr, *bra_k_ptr, vf_local);
                const double t_compute_end = wall_time();
                add_owner_compute_wall_time(t_k_done, t_compute_end);

                for (int i = vf_range.begin; i < vf_range.end; ++i) {
                    Kf[i] += resultcolumn[i - vf_range.begin];
                }

                release_finished_kbatch(bra_range);

                if (world.rank() == 0) {
                    print("OVERLAP_OP row_owner kbatch_prefetch_hit=", int(kbatch_prefetch_hit),
                          " col=[", vf_range.begin, ",", vf_range.end,
                          ") row=[", bra_range.begin, ",", bra_range.end,
                          ") vf_fetch=", t_vf_done - t_op_start,
                          " k_fetch=", t_k_done - t_vf_done,
                          " compute=", t_compute_end - t_k_done);
                }
                return Kf;
            }
            // -------- end row-owner branch --------

            // detect prefetch hit: vf_prefetch_.has_next matching this vf_range means
            // the previous task's run() pre-fetched our column batch
            const bool prefetch_hit = use_owner_aware_fetch()
                and vf_prefetch_.has_next
                and same_range(vf_prefetch_.next_range, normalize_range(vf_range, long(vf_batch.size())));

            vecfuncT bra_local;
            vecfuncT vf_local_sym;        // sym_p2p: cloud-fetched vf (column); aliased by vf_work
            const vecfuncT* bra_work = &bra_batch;
            const vecfuncT* vf_work = &vf_batch;
            if (use_sym_p2p_owner_algorithm() and use_cloud_batch_fetch()) {
                // One-set cloud fetch: psi over the row (bra) and column (vf) ranges.
                // One range is owned by this rank (local read), the other is transient
                // (remote p2p). Both come from the single DIM_VF record set. Because
                // bra==ket==vf, these batches also serve the kernels' ket role
                // (the offdiagonal kernel uses bra_work/vf_work as ket_rows/ket_columns).
                MADNESS_CHECK_THROW(cloud_ptr != nullptr, "cloud_ptr is null in sym_p2p fetch");
                Cloud& cloud = *cloud_ptr;
                const typename Cloud::keyT salt = batch_salt(vket);
                // Bounded-LRU fetch (sym_batch_lru_fetch): the cache-aware task order
                // (shuffle_partition_by_owner) keeps an owner's local batch and the
                // rotating transient adjacent, so consecutive tasks reuse the resident
                // batch instead of re-fetching. The LRU is keyed by the cloud record
                // key and is subworld-scoped (clear_local_caches) — NOT the cloud-side
                // cache_result path, which segfaults across protocol transitions.
                // Copy the handles out of the returned ref before the next fetch, since
                // an insert there may relocate the slot.
                bra_local = sym_batch_lru_fetch(world, cloud, batch_record_key(salt, DIM_VF, bra_range));
                if (diagonal_block) {
                    vf_local_sym = bra_local;     // same range -> same data
                } else {
                    vf_local_sym = sym_batch_lru_fetch(world, cloud, batch_record_key(salt, DIM_VF, vf_range));
                }
                bra_work = &bra_local;
                vf_work = &vf_local_sym;
            } else if (use_owner_aware_fetch()) {
                bra_local = fetch_batch_with_cache(world, bra_batch, bra_range, bra_cache_);
                bra_work = &bra_local;
                vf_work = &acquire_current_vf(world, vf_batch, vf_range);
            }
            const double t_bravf_done = wall_time();

            if (symmetric and diagonal_block) {
                vecfuncT ket_batch = use_sym_p2p_owner_algorithm()
                        ? *bra_work        // one-set: ket == bra over the (diagonal) range
                        : use_owner_aware_fetch()
                        ? fetch_range_with_cache(world, vket, bra_range, ket_cache_)
                        : bra_range.copy_batch(vket);
                const double t_ket_done = wall_time();
                vecfuncT resultcolumn;
                if (algorithm_==small_memory_symmetric_mt or algorithm_==small_memory_symmetric_mt_owner
                    or algorithm_==small_memory_symmetric_p2p_owner) {
                    resultcolumn = compute_diagonal_batch_in_symmetric_matrix_smallmem_symmetric(world, ket_batch,
                                                                                                  *bra_work, *vf_work);
                } else {
                    resultcolumn = compute_diagonal_batch_in_symmetric_matrix(world, ket_batch, *bra_work, *vf_work);
                }
                const double t_compute_end = wall_time();
                add_owner_compute_wall_time(t_ket_done, t_compute_end);

                for (int i = vf_range.begin; i < vf_range.end; ++i){
                    Kf[i] += resultcolumn[i - vf_range.begin];}

                if (world.rank()==0) {
                    print("OVERLAP_OP task diag prefetch_hit=",int(prefetch_hit),
                          " col=[",vf_range.begin,",",vf_range.end,
                          ") row=[",bra_range.begin,",",bra_range.end,
                          ") bravf_fetch=",t_bravf_done-t_op_start,
                          " ket_fetch=",t_ket_done-t_bravf_done,
                          " compute=",t_compute_end-t_ket_done);
                }

            } else if (symmetric and not diagonal_block) {
                std::pair<vecfuncT, vecfuncT> resultpair;
                const double t_ket_done = wall_time(); // ket fetch is inside compute for offdiag
                if (algorithm_==small_memory_symmetric_mt or algorithm_==small_memory_symmetric_mt_owner
                    or algorithm_==small_memory_symmetric_p2p_owner) {
                    resultpair = compute_offdiagonal_batch_in_symmetric_matrix_smallmem_symmetric(world, vket,
                                                                                                    *bra_work, *vf_work);
                } else {
                    resultpair = compute_offdiagonal_batch_in_symmetric_matrix(world, vket, *bra_work, *vf_work);
                }
                const double t_compute_end = wall_time();
                add_owner_compute_wall_time(t_ket_done, t_compute_end);
                auto& resultcolumn = resultpair.first;
                auto& resultrow = resultpair.second;

                for (int i = bra_range.begin; i < bra_range.end; ++i){
                    Kf[i] += resultcolumn[i - bra_range.begin];}
                for (int i = vf_range.begin; i < vf_range.end; ++i){
                    Kf[i] += resultrow[i - vf_range.begin];}

                if (world.rank()==0) {
                    print("OVERLAP_OP task offdiag prefetch_hit=",int(prefetch_hit),
                          " col=[",vf_range.begin,",",vf_range.end,
                          ") row=[",bra_range.begin,",",bra_range.end,
                          ") bravf_fetch=",t_bravf_done-t_op_start,
                          " ket_fetch=",t_ket_done-t_bravf_done,
                          " compute=",t_compute_end-t_ket_done);
                }

            } else {
                vecfuncT ket_batch = use_owner_aware_fetch()
                        ? fetch_range_with_cache(world, vket, bra_range, ket_cache_)
                        : bra_range.copy_batch(vket);
                const double t_ket_done = wall_time();
                vecfuncT resultcolumn = compute_batch_in_asymmetric_matrix(world, ket_batch, *bra_work, *vf_work);
                const double t_compute_end = wall_time();
                add_owner_compute_wall_time(t_ket_done, t_compute_end);
                for (int i = vf_range.begin; i < vf_range.end; ++i)
                    Kf[i] += resultcolumn[i - vf_range.begin];

                if (world.rank()==0) {
                    print("OVERLAP_OP task asym prefetch_hit=",int(prefetch_hit),
                          " col=[",vf_range.begin,",",vf_range.end,
                          ") row=[",bra_range.begin,",",bra_range.end,
                          ") bravf_fetch=",t_bravf_done-t_op_start,
                          " ket_fetch=",t_ket_done-t_bravf_done,
                          " compute=",t_compute_end-t_ket_done);
                }
            }
            if (use_owner_aware_fetch()) release_finished_vf(vf_range);

            // per-task profiling: stamp end + peak RSS and emit one JSON line (rank 0).
            // Wall end is taken here (after the compute branches) so it bounds the whole
            // operator() span; the accumulate/finalize overhead is downstream in run().
            if (prof_on) {
                prof_.wall_end = wall_time();
                prof_.peak_rss_gb = madness::get_rss_usage_in_GB();
                if (world.rank() == 0) exch_write_task_profile(prof_);
            }
            return Kf;
        }

        /// compute a batch of the exchange matrix, with identical ranges, exploiting the matrix symmetry

        /// \param subworld     the world we're computing in
        /// \param cloud        where to store the results
        /// \param bra_batch    the bra batch of orbitals (including the nuclear correlation factor square)
        /// \param ket_batch    the ket batch of orbitals, i.e. the orbitals to premultiply with
        /// \param vf_batch     the argument of the exchange operator
        vecfuncT compute_diagonal_batch_in_symmetric_matrix(World& subworld,
                                                            const vecfuncT& ket_batch,      // is batched
                                                            const vecfuncT& bra_batch,      // is batched
                                                            const vecfuncT& vf_batch        // is batched
        ) const {
            double mul_tol = 0.0;
            double symmetric = true;
            auto poisson = Exchange<double, 3>::ExchangeImpl::set_poisson(subworld, lo);
            return Exchange<T, NDIM>::ExchangeImpl::compute_K_tile(subworld, bra_batch, ket_batch, vf_batch, poisson, symmetric,
                                                     mul_tol);
        }

        /// scaffold for small-memory symmetric diagonal tiles
        vecfuncT compute_diagonal_batch_in_symmetric_matrix_smallmem_symmetric(World& subworld,
                                                                                const vecfuncT& ket_batch,
                                                                                const vecfuncT& bra_batch,
                                                                                const vecfuncT& vf_batch) const {
            const double cpu0 = process_cpu_time();
            MADNESS_CHECK_THROW(ket_batch.size()==bra_batch.size(),
                                "smallmem_sym_mt diagonal: ket/bra batch size mismatch");
            MADNESS_CHECK_THROW(vf_batch.size()==bra_batch.size(),
                                "smallmem_sym_mt diagonal: vf/bra batch size mismatch");

            const long n = vf_batch.size();
            vecfuncT resultcolumn = zero_functions_compressed<T, NDIM>(subworld, n);
            auto poisson = Exchange<double, 3>::ExchangeImpl::set_poisson(subworld, lo);

            // TEMP debug knob.
            const double eff_tol = (smallmem_mul_tol_ >= 0.0) ? smallmem_mul_tol_ : mul_tol*0.1;

            // per-component wall timing (MAD_EXCH_TASK_PROFILE). Honest without added
            // fences: every op below runs with fence=true, so each completes before the
            // next. `tick` folds the elapsed into the per-task accumulator.
            const bool prof_on = exch_task_profile_enabled();
            auto tick = [&](double& acc, double t0) { if (prof_on) acc += wall_time() - t0; };

            for (long i = 0; i < n; ++i) {
                // Build N_ij for the upper-triangle of this diagonal tile row.
                vecfuncT vf_subset(vf_batch.begin(), vf_batch.begin() + i + 1);
                double _t = prof_on ? wall_time() : 0.0;
                vecfuncT psif = mul_sparse(subworld, bra_batch[i], vf_subset, eff_tol, true, false);
                tick(prof_.mul1_wall, _t);
                if (log_diagnostics_) {
                    print("smallmem_sym_mt diagonal i=", i, " psif.size()=", psif.size());
                }
                _t = prof_on ? wall_time() : 0.0;
                truncate(subworld, psif);
                tick(prof_.truncate_wall, _t);
                _t = prof_on ? wall_time() : 0.0;
                psif = apply(subworld, *poisson.get(), psif);
                tick(prof_.apply_wall, _t);
                _t = prof_on ? wall_time() : 0.0;
                truncate(subworld, psif);
                tick(prof_.truncate_wall, _t);
                make_redundant(subworld, psif, true);

                // Assemble full row update within the tile and accumulate by vector gaxpy.
                vecfuncT update_i = zero_functions_compressed<T, NDIM>(subworld, n);
                compress(subworld, update_i);

                _t = prof_on ? wall_time() : 0.0;
                vecfuncT row_contrib = mul_sparse(subworld, ket_batch[i], psif, eff_tol, true, false);
                tick(prof_.mul2_wall, _t);
                if (log_diagnostics_) {
                    print("smallmem_sym_mt diagonal i=", i, " vpsi.size()=", row_contrib.size());
                }
                compress(subworld, row_contrib);
                for (long j = 0; j <= i; ++j) {
                    update_i[j] += row_contrib[j];
                }

                for (long j = 0; j < i; ++j) {
                    vecfuncT psif_single(1, psif[j]);
                    _t = prof_on ? wall_time() : 0.0;
                    vecfuncT mirrored = mul_sparse(subworld, ket_batch[j], psif_single, eff_tol, true, false);
                    tick(prof_.mul2_wall, _t);
                    compress(subworld, mirrored);
                    update_i[i] += mirrored[0];
                }

                gaxpy(subworld, 1.0, resultcolumn, 1.0, update_i);
            }

            const double cpu1 = process_cpu_time();
            add_owner_compute_time(cpu0, cpu1);
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
            const double cpu0 = process_cpu_time();
            double mul_tol = 0.0;
            double symmetric = false;
            auto poisson = Exchange<double, 3>::ExchangeImpl::set_poisson(subworld, lo);
            auto result = Exchange<T, NDIM>::ExchangeImpl::compute_K_tile(subworld, bra_batch, ket_batch, vf_batch, poisson, symmetric,
                                                                           mul_tol);
            const double cpu1 = process_cpu_time();
            add_owner_compute_time(cpu0, cpu1);
            return result;
        }

        /// Small-memory asymmetric tile kernel for small_memory_mt_owner.
        ///
        /// Computes the contribution of one k-batch (size n, the rotating inner-sum
        /// dimension) into the held i-batch (size m, the result/row dimension):
        ///
        ///     resultcolumn[i] = \sum_{k in k-batch} mo_ket[k] \cdot P( mo_bra[k] * vf[i] )
        ///
        /// Memory: holds m psif functions at a time (vs m*n for compute_K_tile).
        /// Iteration is over the k dimension so each psif vector spans only the held m.
        ///
        /// State contract: bra_batch / ket_batch / vf_batch are pre-staged in redundant
        /// state with tree norms by ExchangeImpl::operator() (see make_redundant in the
        /// redundant branch there). We therefore pass do_make_redundant=false to both
        /// mul_sparse calls and skip the per-iteration make_redundant + fence pair that
        /// the default would otherwise issue. psif is reset to redundant after apply +
        /// truncate so it's valid input for the second mul_sparse. mul_tol*0.1 matches
        /// the tighter screening used by K_small_memory; important for larger systems.
        vecfuncT compute_batch_in_asymmetric_matrix_smallmem(
                World& subworld,
                const vecfuncT& ket_batch,       // size n (paired with bra by inner-sum index)
                const vecfuncT& bra_batch,       // size n
                const vecfuncT& vf_batch) const  // size m (held i-batch)
        {
            const double cpu0 = process_cpu_time();
            const long n = static_cast<long>(bra_batch.size());
            const long m = static_cast<long>(vf_batch.size());
            MADNESS_CHECK_THROW(ket_batch.size() == bra_batch.size(),
                                "smallmem_mt_owner: ket/bra batch size mismatch");

            vecfuncT resultcolumn = zero_functions_compressed<T, NDIM>(subworld, m);
            auto poisson = Exchange<double, 3>::ExchangeImpl::set_poisson(subworld, lo);

            // TEMP debug knob.
            const double eff_tol = (smallmem_mul_tol_ >= 0.0) ? smallmem_mul_tol_ : mul_tol*0.1;

            for (long k = 0; k < n; ++k) {
                // psif[i] = bra[k] * vf[i]. Inputs are already redundant (universe-side).
                double cpu_phase0 = process_cpu_time();
                vecfuncT psif = mul_sparse(subworld, bra_batch[k], vf_batch, eff_tol, true, false);
                if (subworld.rank() == 0) {
                    print("smallmem_mt_owner asym k=", k, " psif.size()=", psif.size());
                }
                truncate(subworld, psif);
                double cpu_phase1 = process_cpu_time();
                mul1_timer += long((cpu_phase1 - cpu_phase0) * 1000l);

                cpu_phase0 = process_cpu_time();
                psif = apply(subworld, *poisson.get(), psif);
                truncate(subworld, psif);
                cpu_phase1 = process_cpu_time();
                apply_timer += long((cpu_phase1 - cpu_phase0) * 1000l);

                // apply + truncate leave psif in compressed state; convert to redundant
                // (compress(redundant, ...) computes the tree norms as part of the
                // conversion) so it's a valid input for the second mul_sparse.
                cpu_phase0 = process_cpu_time();
                make_redundant(subworld, psif, true);

                // update[i] = ket[k] * psif[i]; accumulate into resultcolumn.
                vecfuncT update = mul_sparse(subworld, ket_batch[k], psif, eff_tol, true, false);
                if (subworld.rank() == 0) {
                    print("smallmem_mt_owner asym k=", k, " vpsi.size()=", update.size());
                }
                compress(subworld, update);
                gaxpy(subworld, 1.0, resultcolumn, 1.0, update);
                cpu_phase1 = process_cpu_time();
                mul2_timer += long((cpu_phase1 - cpu_phase0) * 1000l);
            }

            const double cpu1 = process_cpu_time();
            add_owner_compute_time(cpu0, cpu1);
            return resultcolumn;
        }

        /// compute a batch of the exchange matrix, with non-identical ranges

        /// \param subworld     the world we're computing in
        /// \param cloud        where to store the results
        /// \param bra_batch    the bra batch of orbitals (including the nuclear correlation factor square)
        /// \param ket_batch    the ket batch of orbitals, i.e. the orbitals to premultiply with
        /// \param vf_batch     the argument of the exchange operator
        std::pair<vecfuncT, vecfuncT> compute_offdiagonal_batch_in_symmetric_matrix(World& subworld,
                                                                                    const vecfuncT& ket, // not batched
                                                                                    const vecfuncT& bra_batch, // batched
                                                                                    const vecfuncT& vf_batch) const; // batched

        /// scaffold for small-memory symmetric offdiagonal tiles
        std::pair<vecfuncT, vecfuncT> compute_offdiagonal_batch_in_symmetric_matrix_smallmem_symmetric(
                World& subworld, const vecfuncT& ket, const vecfuncT& bra_batch, const vecfuncT& vf_batch) const {
            MADNESS_CHECK_THROW(bra_batch.size() > 0 and vf_batch.size() > 0,
                                "smallmem_sym_mt offdiagonal: empty tile batch");

            // Row range corresponds to bra_batch, column range corresponds to vf_batch.
            auto row_range = batch.input[1];
            auto column_range = batch.input[0];
            MADNESS_CHECK_THROW(row_range.size() == long(bra_batch.size()),
                                "smallmem_sym_mt offdiagonal: row range mismatch");
            MADNESS_CHECK_THROW(column_range.size() == long(vf_batch.size()),
                                "smallmem_sym_mt offdiagonal: column range mismatch");

            // ket vectors corresponding to tile row and tile column ranges.
            // One-set (sym_p2p): bra==ket==vf, so ket_rows==bra_batch and
            // ket_columns==vf_batch — reuse them directly instead of re-fetching
            // (the full `ket` is unused on this path).
            vecfuncT ket_rows = use_sym_p2p_owner_algorithm()
                    ? bra_batch
                    : use_owner_aware_fetch()
                    ? fetch_range_with_cache(subworld, ket, row_range, ket_cache_)
                    : row_range.copy_batch(ket);
            vecfuncT ket_columns = use_sym_p2p_owner_algorithm()
                    ? vf_batch
                    : use_owner_aware_fetch()
                    ? fetch_range_transient(subworld, ket, column_range)
                    : column_range.copy_batch(ket);
            MADNESS_CHECK_THROW(ket_rows.size() == bra_batch.size(),
                                "smallmem_sym_mt offdiagonal: ket_rows size mismatch");
            MADNESS_CHECK_THROW(ket_columns.size() == vf_batch.size(),
                                "smallmem_sym_mt offdiagonal: ket_columns size mismatch");

            const double cpu0 = process_cpu_time();
            const long nrow = bra_batch.size();
            const long ncolumn = vf_batch.size();
            vecfuncT resultcolumn = zero_functions_compressed<T, NDIM>(subworld, nrow);   // maps to bra_range
            vecfuncT resultrow = zero_functions_compressed<T, NDIM>(subworld, ncolumn);    // maps to vf_range
            auto poisson = Exchange<double, 3>::ExchangeImpl::set_poisson(subworld, lo);

            // TEMP debug knob.
            const double eff_tol = (smallmem_mul_tol_ >= 0.0) ? smallmem_mul_tol_ : mul_tol*0.1;

            // per-component wall timing (MAD_EXCH_TASK_PROFILE); honest without added
            // fences (every op runs with fence=true). See the diagonal kernel.
            const bool prof_on = exch_task_profile_enabled();
            auto tick = [&](double& acc, double t0) { if (prof_on) acc += wall_time() - t0; };

            for (long irow = 0; irow < nrow; ++irow) {
                // Build N_ij for this offdiagonal tile row, all tile columns.
                double _t = prof_on ? wall_time() : 0.0;
                vecfuncT psif = mul_sparse(subworld, bra_batch[irow], vf_batch, eff_tol, true, false);
                tick(prof_.mul1_wall, _t);
                if (log_diagnostics_) {
                    print("smallmem_sym_mt offdiag irow=", irow, " psif.size()=", psif.size());
                }
                _t = prof_on ? wall_time() : 0.0;
                truncate(subworld, psif);
                tick(prof_.truncate_wall, _t);
                _t = prof_on ? wall_time() : 0.0;
                psif = apply(subworld, *poisson.get(), psif);
                tick(prof_.apply_wall, _t);
                _t = prof_on ? wall_time() : 0.0;
                truncate(subworld, psif);
                tick(prof_.truncate_wall, _t);
                make_redundant(subworld, psif, true);

                // Row accumulation for vf-range: resultrow[j] += ket_row[irow] * N_ij.
                _t = prof_on ? wall_time() : 0.0;
                vecfuncT row_update = mul_sparse(subworld, ket_rows[irow], psif, eff_tol, true, false);
                tick(prof_.mul2_wall, _t);
                if (log_diagnostics_) {
                    print("smallmem_sym_mt offdiag irow=", irow, " vpsi.size()=", row_update.size());
                }
                compress(subworld, row_update);
                gaxpy(subworld, 1.0, resultrow, 1.0, row_update);

                // Mirror accumulation for bra-range: resultcolumn[irow] += sum_j ket_column[j] * N_ij.
                _t = prof_on ? wall_time() : 0.0;
                auto column_update = dot(subworld, ket_columns, psif, true, false, eff_tol);
                tick(prof_.mul2_wall, _t);
                vecfuncT single_column_update = zero_functions_compressed<T, NDIM>(subworld, nrow);
                single_column_update[irow] = copy(column_update);
                gaxpy(subworld, 1.0, resultcolumn, 1.0, single_column_update);
            }

            const double cpu1 = process_cpu_time();
            add_owner_compute_time(cpu0, cpu1);
            return std::make_pair(resultcolumn, resultrow);
        }

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
                size_t iend = std::min(ilo+ntile,mo_bra.size());
                vecfuncT tmp_mo_bra(mo_bra.begin()+ilo,mo_bra.begin()+iend);

                // mul_sparse legacy for reference
                //
                //cpu0 = process_cpu_time();
                //auto tmp_psif = mul_sparse(world, vket[i], tmp_mo_bra, mul_tol*0.1, true, false);
                //cpu1 = process_cpu_time();
                //for (unsigned int i=0; i<tmp_psif.size(); ++i){
                //    print(ilo+i, "mul_sparse output ", tmp_psif[i].tree_size());
                //}
                //mul1_timer += long((cpu1 - cpu0) * 1000l);
                //truncate(world, tmp_psif);
                //for (unsigned int i=0; i<tmp_psif.size(); ++i){
                //    print(ilo+i, "mul_sparse truncated ", tmp_psif[i].tree_size());
                //}

                // vector size 1 function call (debug)
                //
                //vecfuncT tmp_psif2;
                //for (unsigned int i=0; i<tmp_mo_bra.size(); ++i){
                //    //vecfuncT v;
                //    //v.push_back(copy(tmp_mo_bra[i]));
                //    //auto res = mul_sparse(world, vket[0], v, mul_tol*0.1, true, false);
                //    auto res = mul_sparse_debug(vket[0], tmp_mo_bra[i], mul_tol*0.1, true, false, false, true);
                //    print(ilo+i, "mw_mul output size ", res.tree_size());
                //    res.truncate();
                //    print(ilo+i, "mw_mul truncated size ", res.tree_size());
                //    tmp_psif2.push_back(res);
                //}

        //        for (unsigned int i=0; i<tmp_mo_bra_redundant.size(); ++i){
        //            print("START ", ilo+i);
        //            tmp_mo_bra_redundant[i].print_tree();
        //            print("END ", ilo+i);
        //        }
        //        print("start vket");
        //        vket_redundant.print_tree();
        //        print("end vket");

                cpu0 = process_cpu_time();
                auto tmp_psif2 = mul_sparse(world, vket[0], tmp_mo_bra, mul_tol*0.1, true, false);
                cpu1 = process_cpu_time();
                mul1_timer += long((cpu1 - cpu0) * 1000l);
                
                // mul_sparse (new) for mw screening
                //
              //  cpu0 = process_cpu_time();
              //  //vket_redundant.compress(true);
              //  //vket_redundant.make_redundant(true);
              //  //compress(tmp_mo_bra_redundant);
              //  //make_redundant(tmp_mo_bra, true);
              //  vecfuncT tmp_psif2;
              //  for (unsigned int i=0; i<tmp_mo_bra.size(); ++i){
              //      //tmp_mo_bra_redundant[i].compress(true);
              //      //tmp_mo_bra[i].make_redundant(true);
              //      auto res = mul_sparse_debug(vket[0], tmp_mo_bra[i], mul_tol*0.1, true, false, false, true);
              //      tmp_psif2.push_back(res);
              //  }
              //  cpu1 = process_cpu_time();
              //  mul1_timer += long((cpu1 - cpu0) * 1000l);
                for (unsigned int i=0; i<tmp_psif2.size(); ++i){
                    print(ilo+i, "mw_mul output ", tmp_psif2[i].tree_size());
                }
                cpu0 = process_cpu_time();
                truncate(world, tmp_psif2);
                cpu1 = process_cpu_time();
                mul1_truncate_timer += long((cpu1 - cpu0) * 1000l);
                for (unsigned int i=0; i<tmp_psif2.size(); ++i){
                    print(ilo+i, "mw_mul truncated ", tmp_psif2[i].tree_size());
                }

                cpu0 = process_cpu_time();
                //tmp_psif = apply(world, *poisson.get(), tmp_psif);
                tmp_psif2 = apply(world, *poisson.get(), tmp_psif2);
                cpu1 = process_cpu_time();
                print("finished apply");
                apply_timer += long((cpu1 - cpu0) * 1000l);
                //truncate(world, tmp_psif);
                truncate(world, tmp_psif2);
                print("finished truncate");

                cpu0 = process_cpu_time();
                vecfuncT tmp_mo_ket(mo_ket.begin()+ilo,mo_ket.begin()+iend);
                //auto tmp_Kf = dot(world, tmp_mo_ket, tmp_psif, true, false, false, mul_tol*0.01);
                //for (unsigned int i=0; i<tmp_mo_ket.size(); ++i){
                //    tmp_psif2[i].compress(true);
                //    tmp_psif2[i].make_redundant(true);
                //}
                make_redundant(world, tmp_psif2, true);
                auto tmp_Kf = dot(world, tmp_mo_ket, tmp_psif2, true, false, mul_tol*0.1);
                print("finished dot");
                cpu1 = process_cpu_time();
                mul2_timer += long((cpu1 - cpu0) * 1000l);

                Kf[0] += tmp_Kf;
                cpu0 = process_cpu_time();
                truncate(world, Kf);
                cpu1 = process_cpu_time();
                mul2_truncate_timer += long((cpu1 - cpu0) * 1000l);
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
                double cpu0=process_cpu_time();
                SafeMPI::Intracomm comm = subworld_ptr->mpi.comm();
                std::shared_ptr<World> fetching_world(new World(comm.Clone()));
                std::shared_ptr<World> executing_world(new World(comm.Clone()));
                double cpu1=process_cpu_time();
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

                        double cpu0 = process_cpu_time();
                        auto tmp_psif = mul_sparse(world, phi, mo_bra, mul_tol);
                        truncate(world, tmp_psif);
                        double cpu1 = process_cpu_time();
                        mul1_timer += long((cpu1 - cpu0) * 1000l);

                        cpu0 = process_cpu_time();
                        tmp_psif = apply(world, *poisson.get(), tmp_psif);
                        truncate(world, tmp_psif);
                        cpu1 = process_cpu_time();
                        apply_timer += long((cpu1 - cpu0) * 1000l);

                        cpu0 = process_cpu_time();
                        auto tmp_Kf = dot(world, mo_ket, tmp_psif);
                        cpu1 = process_cpu_time();
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
                            double t0=process_cpu_time();
                            print("fetching tile",tile.ilo,"into world",executing_world->id());
                            std::tie(tmp_mo_bra1,tmp_mo_ket1)=fetch_data(*executing_world,tiles[itile]);
                            fetching_world->gop.set_forbid_fence(false);
                            double t2=process_cpu_time();
                            executing_world->gop.fence();
                            double t1=process_cpu_time();
                            total_fetch_time += (t1 - t0);
                            total_fetch_spawn_time += (t2 - t0);
                        }

                        if (itile>=0) {
                            double t0=process_cpu_time();
                            fetching_world->gop.set_forbid_fence(true);
                            if (itile<tiles.size()-1) {
                                // fetch data into fetching_world while computing in executing_world
                                print("fetching tile",tiles[itile+1].ilo,"into world",fetching_world->id()," at time ",wall_time());
                                std::tie(tmp_mo_bra2,tmp_mo_ket2)=fetch_data(*fetching_world,tiles[itile+1]);
                            }
                            fetching_world->gop.set_forbid_fence(false);
                            double t2=process_cpu_time();
                            // uncomment the next line to enforce that fetching is finished before executing
                            // fetching_world->gop.fence();
                            double t1=process_cpu_time();
                            total_fetch_time += (t1 - t0);
                            total_fetch_spawn_time += (t2 - t0);

                            print("executing tile",tile.ilo,"in world",executing_world->id());
                            double dpu0=process_cpu_time();
                            Kf[0]+=execute(*executing_world,poisson1,phi1,tmp_mo_bra1,tmp_mo_ket1);
                            double dpu1=process_cpu_time();
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
                    }
                } // objects living in the two worlds must be destroyed before the worlds are freed

                // deferred destruction of WorldObjects happens here
                fetching_world->gop.fence();
                executing_world->gop.fence();
                double cpu2=process_cpu_time();
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
