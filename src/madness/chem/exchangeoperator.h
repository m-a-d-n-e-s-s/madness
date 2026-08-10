#ifndef SRC_APPS_CHEM_EXCHANGEOPERATOR_H_
#define SRC_APPS_CHEM_EXCHANGEOPERATOR_H_

#include<madness.h>
#include<madness/world/cloud.h>
#include<madness/mra/macrotaskq.h>
#include<madness/chem/SCFOperators.h>

#include <algorithm>
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
    }

public:
    nlohmann::json gather_timings(World& world) const {
        double t1 = double(mul1_timer) * 0.001;
        double t2 = double(apply_timer) * 0.001;
        double t3 = double(mul2_timer) * 0.001;
        world.gop.sum(t1);
        world.gop.sum(t2);
        world.gop.sum(t3);
        nlohmann::json j;
        j["multiply1"] = t1;
        j["apply"] = t2;
        j["multiply2"] = t3;
        j["total"] = elapsed_time;
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
            const long resident = MacroTaskExchangeSimple::batch_cache_hits();
            const long fetched = MacroTaskExchangeSimple::batch_cache_misses();
            if (resident + fetched > 0)
                printf(" operand batches resident/fetched %6ld /%6ld\n", resident, fetched);
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
        // how often a task found its operand batch already resident instead of fetching it
        // from its owner; zero fetches at all means the owner-pinned path did not run
        j["batch_cache_hits"] = MacroTaskExchangeSimple::batch_cache_hits();
        j["batch_cache_misses"] = MacroTaskExchangeSimple::batch_cache_misses();
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

    mutable nlohmann::json statistics;  ///< statistics of the Cloud (timings, memory)  and of the parameters of this run

    class MacroTaskExchangeSimple : public MacroTaskOperationBase {

        long nresult;
        double lo = 1.e-4;
        double mul_tol = 1.e-7;
        bool symmetric = false;
        /// pin each task to a rank that owns one of its batches, over the owner-pinned split
        bool owner_pinned = false;
        long granularity_level = 1;
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

        static void clear_local_caches() { batch_cache_.clear(); }

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
                return *resident;
            }
            ++batch_cache_misses_;
            const bool owned = (cloud.batch_owner(record) == universe_rank_);
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
                                const long granularity_level = 1, const long universe_rank = 0)
                : nresult(nresult), lo(lo), mul_tol(mul_tol), symmetric(symmetric),
                  owner_pinned(owner_pinned), granularity_level(granularity_level),
                  universe_rank_(universe_rank) {
            partitioner.reset(new MacroTaskPartitionerExchange(symmetric, owner_pinned, granularity_level));
        }

        /// how often a task's operand batch was already resident, and how often it had to be
        /// fetched from the rank owning it. No fetches at all means the path never ran.
        static long batch_cache_hits() { return batch_cache_hits_; }
        static long batch_cache_misses() { return batch_cache_misses_; }

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

            const std::vector<Batch_1D> split =
                    exchange_sym_owner_split(nresult, nsubworld, granularity_level);
            const long M = long(split.size());
            std::map<long,long> begin_to_index;
            for (long k = 0; k < M; ++k) begin_to_index[split[k].begin] = k;

            const std::map<std::pair<long,long>,long> assignment =
                    exchange_sym_round_robin_assign(nsubworld, M);

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
        /// column, row and ket operand.
        void store_batches(World& world, World& subworld, Cloud& cloud, const argtupleT& argtuple,
                           const long nsubworld) {
            if (not owner_pinned) return;
            // Batch k is owned by rank (k mod nsubworld), so a batch index has to name a
            // universe rank, and that only holds with one subworld per rank. Anything else
            // would register the routing to the wrong ranks and read the wrong coefficients,
            // so say so rather than compute something wrong.
            MADNESS_CHECK_THROW(nsubworld == world.size(),
                                "owner-pinned exchange needs one subworld per rank");
            const vecfuncT& mo_ket = std::get<2>(argtuple);
            // one fence up front, so store_batch does not need one per function
            world.gop.fence();
            const long salt = exchange_batch_salt(mo_ket);
            const std::vector<Batch_1D> split =
                    exchange_sym_owner_split(mo_ket.size(), nsubworld, granularity_level);

            // A cross-world copy picks its process map from the target world, but anything
            // on that path reading the process-wide default would route to universe ranks
            // that do not exist in a size-1 subworld. Point the default at the subworld for
            // the duration and restore it afterwards.
            auto saved_pmap = FunctionDefaults<NDIM>::get_pmap();
            FunctionDefaults<NDIM>::set_default_pmap(subworld);

            std::vector<std::pair<long, vecfuncT>> owned;
            for (long k = 0; k < long(split.size()); ++k) {
                const Batch_1D& r = split[k];
                const long record = exchange_batch_record_key(salt, EXCHANGE_BATCH_VF, r);
                cloud.register_batch_owner(record, ProcessID(k % nsubworld));
                if (k % nsubworld == world.rank()) {
                    vecfuncT local(r.size());
                    for (long i = r.begin; i < r.end; ++i)
                        local[i - r.begin] = copy(subworld, mo_ket[i], /*fence=*/false);
                    owned.emplace_back(record, std::move(local));
                }
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

        /// the owner-pinned path fetches its operand batches from the cloud itself
        bool handles_own_data_movement() const override { return owner_pinned; }

        /// Drop the cached batches while the subworld holding them is still alive.

        /// Leaving it to ensure_cache_world to notice a new subworld is too late: by then the
        /// cached functions refer to a destroyed world, and merely releasing them walks into
        /// it.
        void cleanup() override {
            clear_local_caches();
            cache_world_id_ = -1;
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

            // Owner-pinned: fetch the two batches this task needs from the cloud. Because
            // bra == ket == vf there, one batch per range serves every operand role, so the
            // ket over a range is just that range's batch.
            vecfuncT bra_owned, vf_owned;
            const vecfuncT* bra_work = &bra_batch;
            const vecfuncT* vf_work = &vf_batch;
            if (owner_pinned) {
                MADNESS_CHECK_THROW(cloud_ptr != nullptr, "owner-pinned exchange: cloud_ptr is null");
                Cloud& cloud = *cloud_ptr;
                const long salt = exchange_batch_salt(vket);
                // copy the handles out of the cache: the reference is only good until the
                // entry is evicted, and the second fetch may evict the first
                bra_owned = fetch_batch(world, cloud,
                                        exchange_batch_record_key(salt, EXCHANGE_BATCH_VF, bra_range));
                vf_owned = diagonal_block
                        ? bra_owned                     // one range, so one batch
                        : fetch_batch(world, cloud,
                                      exchange_batch_record_key(salt, EXCHANGE_BATCH_VF, vf_range));
                bra_work = &bra_owned;
                vf_work = &vf_owned;
            }

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
                auto ket_batch = bra_range.copy_batch(vket);
                vecfuncT resultcolumn = compute_batch_in_asymmetric_matrix(world, ket_batch, bra_batch, vf_batch);
                for (int i = vf_range.begin; i < vf_range.end; ++i)
                    Kf[i] += resultcolumn[i - vf_range.begin];
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

            for (long i = 0; i < n; ++i) {
                double cpu0 = cpu_time();
                const vecfuncT vf_subset(vf_batch.begin(), vf_batch.begin() + i + 1);
                vecfuncT psif = mul_sparse(subworld, bra_batch[i], vf_subset, mul_tol);
                truncate(subworld, psif);
                double cpu1 = cpu_time();
                mul1_timer += long((cpu1 - cpu0) * 1000l);

                cpu0 = cpu_time();
                psif = apply(subworld, *poisson.get(), psif);
                truncate(subworld, psif);
                cpu1 = cpu_time();
                apply_timer += long((cpu1 - cpu0) * 1000l);

                // rows overlap in the columns they touch, so this row's contribution is
                // assembled on its own and accumulated
                cpu0 = cpu_time();
                vecfuncT update = zero_functions_compressed<T, NDIM>(subworld, n);
                vecfuncT row_contrib = mul_sparse(subworld, ket_batch[i], psif, mul_tol);
                compress(subworld, row_contrib);
                for (long j = 0; j <= i; ++j) update[j] += row_contrib[j];
                for (long j = 0; j < i; ++j) {
                    vecfuncT mirrored = mul_sparse(subworld, ket_batch[j], vecfuncT(1, psif[j]), mul_tol);
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
