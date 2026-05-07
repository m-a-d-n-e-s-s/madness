#ifndef SRC_APPS_CHEM_EXCHANGEOPERATOR_H_
#define SRC_APPS_CHEM_EXCHANGEOPERATOR_H_

#include<madness.h>
#include<madness/world/cloud.h>
#include<madness/mra/macrotaskq.h>
#include<madness/chem/SCFOperators.h>
#include<unordered_map>
#include<queue>
#include<random>

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
    static inline std::atomic<long> owner_compute_timer;
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
        world.gop.sum(t1);
        world.gop.sum(t2);
        world.gop.sum(t3);
        world.gop.sum(t4);
        world.gop.sum(t5);
        world.gop.sum(t_fetch_owner);
        world.gop.sum(t_compute_owner);
        nlohmann::json j;
        j["multiply1"] = t1;
        j["truncate1"] = t4;
        j["apply"] = t2;
        j["multiply2"] = t3;
        j["truncate2"] = t5;
        j["owner_fetch"] = t_fetch_owner;
        j["owner_compute"] = t_compute_owner;
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
            printf(" total process cpu time        %8.2fs\n", timings["total_cpu"].template get<double>());
            printf(" total wall time               %8.2fs\n", timings["total"].template get<double>());
        }
    }


    typedef Exchange<T,NDIM>::ExchangeAlgorithm Algorithm;
    Algorithm algorithm_ = multiworld_efficient_row;
    MacroTaskInfo macro_task_info = MacroTaskInfo::preset("default");
    bool replicate_for_debug_ = false;  ///< if true, use StoreFunction policy to pre-replicate all data (zero communication during tasks)
    bool local_accumulation_ = true;    ///< if true (and using owner-aware algorithm), accumulate task results subworld-locally and do one final subworld->universe gaxpy

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
    private:
        static inline std::unordered_map<long, functionT> bra_cache_;
        static inline std::unordered_map<long, functionT> ket_cache_;
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
        static inline long cache_world_id_ = -1;

        /// pre-computed owner map: (col_begin, row_begin) -> owner rank
        /// populated by prepare_owner_assignment() using the fold algorithm
        std::map<std::pair<long,long>, long> owner_map_;

        /// if true, shuffle per-owner task order to reduce synchronized fetch contention
        /// disabled: shuffling destroys row-range locality needed for cache reuse and prefetch hits
        bool shuffle_task_order_ = false;

        bool use_owner_aware_fetch() const { return algorithm_==small_memory_symmetric_mt_owner and not replicate_for_debug_; }

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
            clear_vf_prefetch();
        }

        void add_owner_fetch_time(const double cpu0, const double cpu1) const {
            if (use_owner_aware_fetch()) owner_fetch_timer += long((cpu1 - cpu0) * 1000l);
        }

        void add_owner_compute_time(const double cpu0, const double cpu1) const {
            if (use_owner_aware_fetch()) owner_compute_timer += long((cpu1 - cpu0) * 1000l);
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
                                         const long max_batch_size_input)
                    : symmetric(symmetric) {
                const long min_bs = std::max<long>(1, min_batch_size_input);
                const long max_bs = std::max<long>(min_bs, std::max<long>(1, max_batch_size_input));
                min_batch_size=min_bs;
                max_batch_size=max_bs;
            }

            bool symmetric = false;

            partitionT do_partitioning(const std::size_t& vsize1, const std::size_t& vsize2,
                                       const std::string policy) const override {

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
            partitioner.reset(new MacroTaskPartitionerExchange(symmetric, min_batch_size, max_batch_size));
            name="MacroTaskExchangeSimple";
        }

        /// Pre-compute a load-balanced owner assignment for all tasks in the partition.
        ///
        /// Uses a fold algorithm inspired by the triangle-to-rectangle transformation:
        /// row groups from opposite ends of the triangular task matrix are paired so
        /// that each "rectangle row" has approximately equal task count. Rectangle rows
        /// are then assigned to ranks via greedy (least-loaded-first) scheduling.
        ///
        /// Called automatically by the MacroTask framework (via SFINAE hook) after
        /// partitioning and before per-task owner_hint queries.
        void prepare_owner_assignment(const MacroTaskPartitioner::partitionT& partition, long nsubworld) {
            if (!use_owner_aware_fetch() || nsubworld <= 0 || !symmetric) return;
            owner_map_ = fold_and_assign(partition, nsubworld);
        }

        /// Fold the triangular task list and assign owners for load balance.
        ///
        /// 1. Group tasks by row range (input[1]), giving R row groups ordered by
        ///    row index. Cost is accumulated from task priority (batch area), not
        ///    raw task count, so runt batches are weighted correctly.
        /// 2. Fold: pair row group (half-1-k) with row group (half+k) to form
        ///    rectangle row k. Each rectangle row has ~uniform cost.
        ///    For odd R, the middle row group's tasks are distributed round-robin
        ///    across rectangle rows.
        /// 3. Greedy (LPT) assignment: rectangle rows to ranks by cost.
        /// 4. Build a map from (col_begin, row_begin) -> owner rank for each task.
        static std::map<std::pair<long,long>, long> fold_and_assign(
                const MacroTaskPartitioner::partitionT& partition, long nsubworld) {

            // -- Step 1: group tasks by row range --
            // Use row_begin as the key; accumulate priority-weighted cost.
            struct TaskEntry {
                std::pair<long,long> key;  // (col_begin, row_begin)
                double priority;
            };
            struct RowGroup {
                long row_begin = 0;
                double cost = 0.0;
                std::vector<TaskEntry> tasks;
            };

            // Discover row groups from partition (order by row_begin)
            std::map<long, RowGroup> row_group_map;
            for (const auto& [batch, prio] : partition) {
                const Batch_1D& col_range = batch.input[0];
                const Batch_1D& row_range = (batch.input.size() > 1) ? batch.input[1] : batch.input[0];
                auto& rg = row_group_map[row_range.begin];
                rg.row_begin = row_range.begin;
                rg.cost += prio;
                rg.tasks.push_back({{col_range.begin, row_range.begin}, prio});
            }

            // Flatten into ordered vector
            std::vector<RowGroup> row_groups;
            row_groups.reserve(row_group_map.size());
            for (auto& [key, rg] : row_group_map) {
                row_groups.push_back(std::move(rg));
            }
            const long R = static_cast<long>(row_groups.size());

            // Degenerate case
            if (R == 0) return {};

            // -- Step 2: fold --
            // Pair row group (half-1-k) with row group (half+k) to form rectangle rows.
            // For even R: half = R/2, R/2 rectangle rows.
            // For odd R:  half = R/2, R/2 rectangle rows; the middle row group's
            //             tasks are distributed round-robin across the rectangle rows
            //             to avoid dumping ~50% extra load onto a single rank.
            const long half = R / 2;
            const bool odd = (R % 2 != 0);

            struct RectRow {
                double cost = 0.0;
                std::vector<std::pair<long,long>> task_keys;
            };
            std::vector<RectRow> rect_rows(half);

            for (long k = 0; k < half; ++k) {
                auto& rr = rect_rows[k];
                const long bottom_idx = half + (odd ? 1 : 0) + k;
                const long top_idx = half - 1 - k;

                const auto& bottom = row_groups[bottom_idx];
                const auto& top = row_groups[top_idx];
                rr.cost = bottom.cost + top.cost;
                for (const auto& te : bottom.tasks) rr.task_keys.push_back(te.key);
                for (const auto& te : top.tasks) rr.task_keys.push_back(te.key);
            }

            // Distribute middle row group's tasks round-robin across rectangle rows
            if (odd && half > 0) {
                const auto& mid = row_groups[half];
                for (long t = 0; t < static_cast<long>(mid.tasks.size()); ++t) {
                    auto& rr = rect_rows[t % half];
                    rr.task_keys.push_back(mid.tasks[t].key);
                    rr.cost += mid.tasks[t].priority;
                }
            }

            // -- Step 3: greedy assignment of rectangle rows to ranks --
            // Min-heap: (load, rank_id)
            using heap_entry = std::pair<double, long>;
            std::priority_queue<heap_entry, std::vector<heap_entry>, std::greater<heap_entry>> heap;
            for (long p = 0; p < nsubworld; ++p) {
                heap.push({0.0, p});
            }

            // Assign rectangle rows (sorted by cost descending for LPT scheduling)
            std::sort(rect_rows.begin(), rect_rows.end(),
                      [](const RectRow& a, const RectRow& b) { return a.cost > b.cost; });

            std::map<std::pair<long,long>, long> owner_map;
            for (const auto& rr : rect_rows) {
                auto [load, rank] = heap.top();
                heap.pop();
                for (const auto& key : rr.task_keys) {
                    owner_map[key] = rank;
                }
                heap.push({load + rr.cost, rank});
            }

            // Edge case: odd R with half==0 means R==1, single row group, assign directly
            if (odd && half == 0) {
                const auto& mid = row_groups[0];
                for (const auto& te : mid.tasks) {
                    owner_map[te.key] = 0;
                }
            }

            return owner_map;
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
        void finalize_into(World& /*subworld*/, vecfuncT& universe_result) {
            if (not Kf_local_initialized_) return;
            change_tree_state(Kf_local_, compressed);
            gaxpy(1.0, universe_result, 1.0, Kf_local_, false);
            Kf_local_.clear();
            Kf_local_initialized_ = false;
            Kf_local_world_id_ = -1;
        }

        void set_next_vf_hint(const Batch_1D& next_hint, const bool has_hint) {
            if (not use_owner_aware_fetch()) return;
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
