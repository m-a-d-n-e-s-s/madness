#ifndef SRC_APPS_CHEM_EXCHANGEOPERATOR_H_
#define SRC_APPS_CHEM_EXCHANGEOPERATOR_H_

#include<madness.h>
#include<madness/world/cloud.h>
#include<madness/mra/macrotaskq.h>
#include<madness/chem/SCFOperators.h>
#include<unordered_map>
#include<queue>
#include<random>
#include<numeric>
#include<set>
#include<limits>

namespace madness {

// forward declaration
class SCF;
class Nemo;


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
            printf(" total process cpu time        %8.2fs\n", timings["total_cpu"].template get<double>());
            printf(" total wall time               %8.2fs\n", timings["total"].template get<double>());
        }
    }


    typedef Exchange<T,NDIM>::ExchangeAlgorithm Algorithm;
    Algorithm algorithm_ = multiworld_efficient_row;
    MacroTaskInfo macro_task_info = MacroTaskInfo::preset("default");
    bool replicate_for_debug_ = false;  ///< if true, use StoreFunction policy to pre-replicate all data (zero communication during tasks)
    bool local_accumulation_ = true;    ///< if true (and using owner-aware algorithm), accumulate task results subworld-locally and do one final subworld->universe gaxpy
    bool use_mflex_ = true;             ///< if true (and using owner-aware algorithm), run the m-flex peel search to load-balance the owner assignment
    long mflex_max_exhaustive_ = 5000;  ///< upper bound on C(R,m) for the exhaustive arm of the m-flex peel search

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

    ExchangeImpl& set_local_accumulation(const bool flag) {
        local_accumulation_ = flag;
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

    mutable nlohmann::json statistics;  ///< statistics of the Cloud (timings, memory)  and of the parameters of this run

    class MacroTaskExchangeSimple : public MacroTaskOperationBase {

        long nresult;
        double lo = 1.e-4;
        double mul_tol = 1.e-7;
        bool symmetric = false;
        Algorithm algorithm_ = multiworld_efficient;
    public:
        bool replicate_for_debug_ = false;
        // When true (and use_owner_aware_fetch() is true), the task accumulates
        // each tile's contribution into a subworld-local buffer and performs a
        // single subworld->universe gaxpy in finalize_into() after all tasks on
        // this subworld have run. Defaults to true; flip off to contrast against
        // the original per-task universe accumulation path.
        bool local_accumulation_ = true;
        // When true, run the m-flex peel search before fold_and_assign so the
        // owner assignment minimizes per-rank load delta. Set from ExchangeImpl
        // before submitting the macrotask. Default true.
        bool use_mflex_ = true;
        // Cap on C(R, m) for the exhaustive arm of the m-flex peel search.
        long mflex_max_exhaustive_ = 5000;
    private:
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
            return (algorithm_==small_memory_symmetric_mt_owner or algorithm_==small_memory_mt_owner)
                   and not replicate_for_debug_;
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
            clear_vf_prefetch();
            clear_kbatch_prefetch();
        }

        void add_owner_fetch_time(const double cpu0, const double cpu1) const {
            if (use_owner_aware_fetch()) owner_fetch_timer += long((cpu1 - cpu0) * 1000l);
        }

        void add_owner_compute_time(const double cpu0, const double cpu1) const {
            if (use_owner_aware_fetch()) owner_compute_timer += long((cpu1 - cpu0) * 1000l);
        }

        /// wall-time companion to add_owner_compute_time. Called from operator() to
        /// accumulate wall time wrapping the compute kernel call. Combined with the
        /// CPU-time owner_compute_timer, lets the user attribute the gap to
        /// communication/fences that happened during the compute phase.
        void add_owner_compute_wall_time(const double wall0, const double wall1) const {
            if (use_owner_aware_fetch()) owner_compute_wall_timer += long((wall1 - wall0) * 1000l);
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

        /// custom partitioning for the exchange operator in exchangeoperator.h

        /// arguments are: result[i] += sum_k vket[k] \int 1/r vbra[k] f[i]
        /// with f and vbra being batched, result and vket being passed on as a whole
        class MacroTaskPartitionerExchange : public MacroTaskPartitioner {
        public:
            MacroTaskPartitionerExchange(const bool symmetric, const long min_batch_size_input,
                                         const long max_batch_size_input,
                                         const bool row_owner = false)
                    : symmetric(symmetric), row_owner_(row_owner) {
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

            /// Split a vector of length n into exactly nbatch = min(nsubworld, n)
            /// contiguous batches with sizes differing by at most 1. The first
            /// (n mod nbatch) batches get size ceil(n/nbatch); the rest get
            /// floor(n/nbatch). Eliminates the runt-of-1 case.
            std::vector<Batch_1D> row_owner_split(std::size_t n) const {
                std::vector<Batch_1D> out;
                if (n == 0) return out;
                const long nbatch = std::min<long>(std::max<long>(1, long(nsubworld)), long(n));
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

            partitionT do_partitioning(const std::size_t& vsize1, const std::size_t& vsize2,
                                       const std::string policy) const override {

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
                                const Algorithm algorithm = multiworld_efficient)
                : nresult(nresult), lo(lo), mul_tol(mul_tol), symmetric(symmetric), algorithm_(algorithm) {
            const bool row_owner = (algorithm == small_memory_mt_owner);
            partitioner.reset(new MacroTaskPartitionerExchange(symmetric, min_batch_size, max_batch_size, row_owner));
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
            if (!use_owner_aware_fetch() || !shuffle_task_order_ || owner_map_.empty() || nsubworld <= 0) return;

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

            // Shuffle each owner's task list independently
            for (auto& [owner, tasks] : per_owner) {
                std::mt19937 rng(static_cast<unsigned>(owner * 31 + nsubworld));
                std::shuffle(tasks.begin(), tasks.end(), rng);
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
            return use_owner_aware_fetch() and local_accumulation_;
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
            if (not Kf_local_initialized_) return;
            change_tree_state(Kf_local_, compressed);
            if (algorithm_ == small_memory_mt_owner) {
                truncate(subworld, Kf_local_);
            }
            gaxpy(1.0, universe_result, 1.0, Kf_local_, false);
            Kf_local_.clear();
            Kf_local_initialized_ = false;
            Kf_local_world_id_ = -1;
        }

        void set_next_vf_hint(const Batch_1D& next_hint, const bool has_hint) {
            if (not use_owner_aware_fetch()) return;
            // In small_memory_mt_owner the vf dimension is held, not rotated, so
            // skip the vf-side prefetch state entirely (k-batch state is rotated instead).
            if (use_row_owner_algorithm()) return;
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
            if (use_row_owner_algorithm()) return;     // vf is held in this algorithm
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

        void prefetch_next_bra_async(World& world,
                                     const vecfuncT& mo_bra_full,
                                     const vecfuncT& mo_ket_full) const {
            if (not use_row_owner_algorithm()) return;
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
            cache_world_id_ = -1;
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

            // -------- small_memory_mt_owner: row-owner asymmetric branch --------
            // For this algorithm each rank exclusively owns one i-batch (= vf_range).
            // The held i-batch is fetched once per subworld into held_vf_cache_ and
            // reused across all tasks. The k-batch (= mo_bra + mo_ket paired by inner
            // index = bra_range) rotates task-by-task and is prefetched by the
            // framework's set_next_bra_hint / prefetch_next_bra_async hooks.
            if (use_row_owner_algorithm()) {
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
            const vecfuncT* bra_work = &bra_batch;
            const vecfuncT* vf_work = &vf_batch;
            if (use_owner_aware_fetch()) {
                bra_local = fetch_batch_with_cache(world, bra_batch, bra_range, bra_cache_);
                bra_work = &bra_local;
                vf_work = &acquire_current_vf(world, vf_batch, vf_range);
            }
            const double t_bravf_done = wall_time();

            if (symmetric and diagonal_block) {
                vecfuncT ket_batch = use_owner_aware_fetch()
                        ? fetch_range_with_cache(world, vket, bra_range, ket_cache_)
                        : bra_range.copy_batch(vket);
                const double t_ket_done = wall_time();
                vecfuncT resultcolumn;
                if (algorithm_==small_memory_symmetric_mt or algorithm_==small_memory_symmetric_mt_owner) {
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
                if (algorithm_==small_memory_symmetric_mt or algorithm_==small_memory_symmetric_mt_owner) {
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

            for (long i = 0; i < n; ++i) {
                // Build N_ij for the upper-triangle of this diagonal tile row.
                vecfuncT vf_subset(vf_batch.begin(), vf_batch.begin() + i + 1);
                vecfuncT psif = mul_sparse(subworld, bra_batch[i], vf_subset, mul_tol*0.1, true, false);
                if (subworld.rank() == 0) {
                    print("smallmem_sym_mt diagonal i=", i, " psif.size()=", psif.size());
                }
                truncate(subworld, psif);
                psif = apply(subworld, *poisson.get(), psif);
                truncate(subworld, psif);
                make_redundant(subworld, psif, true);

                // Assemble full row update within the tile and accumulate by vector gaxpy.
                vecfuncT update_i = zero_functions_compressed<T, NDIM>(subworld, n);
                compress(subworld, update_i);

                vecfuncT row_contrib = mul_sparse(subworld, ket_batch[i], psif, mul_tol*0.1, true, false);
                if (subworld.rank() == 0) {
                    print("smallmem_sym_mt diagonal i=", i, " vpsi.size()=", row_contrib.size());
                }
                compress(subworld, row_contrib);
                for (long j = 0; j <= i; ++j) {
                    update_i[j] += row_contrib[j];
                }

                for (long j = 0; j < i; ++j) {
                    vecfuncT psif_single(1, psif[j]);
                    vecfuncT mirrored = mul_sparse(subworld, ket_batch[j], psif_single, mul_tol*0.1, true, false);
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

            for (long k = 0; k < n; ++k) {
                // psif[i] = bra[k] * vf[i]. Inputs are already redundant (universe-side).
                double cpu_phase0 = process_cpu_time();
                vecfuncT psif = mul_sparse(subworld, bra_batch[k], vf_batch, mul_tol*0.1, true, false);
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
                vecfuncT update = mul_sparse(subworld, ket_batch[k], psif, mul_tol*0.1, true, false);
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
            vecfuncT ket_rows = use_owner_aware_fetch()
                    ? fetch_range_with_cache(subworld, ket, row_range, ket_cache_)
                    : row_range.copy_batch(ket);
            vecfuncT ket_columns = use_owner_aware_fetch()
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

            for (long irow = 0; irow < nrow; ++irow) {
                // Build N_ij for this offdiagonal tile row, all tile columns.
                vecfuncT psif = mul_sparse(subworld, bra_batch[irow], vf_batch, mul_tol*0.1, true, false);
                if (subworld.rank() == 0) {
                    print("smallmem_sym_mt offdiag irow=", irow, " psif.size()=", psif.size());
                }
                truncate(subworld, psif);
                psif = apply(subworld, *poisson.get(), psif);
                truncate(subworld, psif);
                make_redundant(subworld, psif, true);

                // Row accumulation for vf-range: resultrow[j] += ket_row[irow] * N_ij.
                vecfuncT row_update = mul_sparse(subworld, ket_rows[irow], psif, mul_tol*0.1, true, false);
                if (subworld.rank() == 0) {
                    print("smallmem_sym_mt offdiag irow=", irow, " vpsi.size()=", row_update.size());
                }
                compress(subworld, row_update);
                gaxpy(subworld, 1.0, resultrow, 1.0, row_update);

                // Mirror accumulation for bra-range: resultcolumn[irow] += sum_j ket_column[j] * N_ij.
                auto column_update = dot(subworld, ket_columns, psif, true, false, mul_tol*0.1);
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
