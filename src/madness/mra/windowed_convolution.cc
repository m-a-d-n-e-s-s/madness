/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
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

/// \file windowed_convolution.cc
/// \brief Out-of-line implementation and explicit instantiation of
///        `WindowedSeparatedConvolution<double, 3>` and opt-in `SeparatedConvolution` bridge.

#include <madness/mra/mra.h>
#include <madness/mra/windowed_quadrature.h>
#include <madness/mra/windowed_ns_stencil.h>
#include <madness/mra/windowed_ns_halo.h>
#include <madness/mra/windowed_ns_sumfact.h>
#include <madness/mra/windowed_convolution.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <mutex>
#include <utility>
#include <vector>

namespace madness {

namespace {

double check_and_get_cubic_cell_width_3d() {
    const Tensor<double>& cw = FunctionDefaults<3>::get_cell_width();
    const double L0 = cw(0L);
    MADNESS_CHECK_THROW(L0 > 0.0,
                        "WindowedSeparatedConvolution: simulation cell width must be positive");
    const double tol = 1e-13 * L0;
    MADNESS_CHECK_THROW(std::abs(cw(1L) - L0) <= tol && std::abs(cw(2L) - L0) <= tol,
                        "WindowedSeparatedConvolution: isotropic 3D kernels require a cubic simulation cell");
    return L0;
}

struct PlanCacheEntry {
    WindowedKernelKind kind = WindowedKernelKind::Coulomb;
    int k = 0;
    double eps = 0.0;
    double mu_user = 0.0;
    double cell_width_L = 0.0;
    int max_tree_level = 0;
    int oversampling_M = 0;
    int gc_rank_Q = 0;
    std::vector<double> custom_coeffs_user;
    std::vector<double> custom_expnts_user;
    std::shared_ptr<const WindowedQuadraturePlan> plan;
    std::shared_ptr<const WindowedStencilCache1D> stencil;
};

std::pair<std::shared_ptr<const WindowedQuadraturePlan>,
          std::shared_ptr<const WindowedStencilCache1D>>
get_or_build_plan_and_stencil(const WindowedQuadratureConfig& cfg) {
    static std::mutex cache_mtx;
    static std::vector<PlanCacheEntry> cache;

    // Hold lock across plan + stencil construction (< 0.5 ms per unique config)
    // so concurrent multiworld MacroTaskQ threads cannot race on
    // GaussianConvolution1DCache<double>::get duplicate-key insertion.
    std::lock_guard<std::mutex> lock(cache_mtx);
    for (const auto& entry : cache) {
        if (entry.kind == cfg.kind &&
            entry.k == cfg.k &&
            entry.eps == cfg.eps &&
            entry.mu_user == cfg.mu_user &&
            entry.cell_width_L == cfg.cell_width_L &&
            entry.max_tree_level == cfg.max_tree_level &&
            entry.oversampling_M == cfg.oversampling_M &&
            entry.gc_rank_Q == cfg.gc_rank_Q &&
            entry.custom_coeffs_user == cfg.custom_coeffs_user &&
            entry.custom_expnts_user == cfg.custom_expnts_user) {
            return {entry.plan, entry.stencil};
        }
    }

    auto plan = std::make_shared<const WindowedQuadraturePlan>(cfg);
    auto stencil = std::make_shared<const WindowedStencilCache1D>(*plan);

    if (cache.size() >= 64) {
        cache.erase(cache.begin());
    }
    cache.push_back(PlanCacheEntry{
        cfg.kind,
        cfg.k,
        cfg.eps,
        cfg.mu_user,
        cfg.cell_width_L,
        cfg.max_tree_level,
        cfg.oversampling_M,
        cfg.gc_rank_Q,
        cfg.custom_coeffs_user,
        cfg.custom_expnts_user,
        plan,
        stencil});
    return {plan, stencil};
}

struct ResolvedOperatorPlan {
    const WindowedQuadraturePlan* plan = nullptr;
    const WindowedStencilCache1D* stencil = nullptr;
    std::shared_ptr<const WindowedQuadraturePlan> ext_plan;
    std::shared_ptr<const WindowedStencilCache1D> ext_stencil;
};

ResolvedOperatorPlan resolve_plan_and_stencil_for_depth(
    const WindowedSeparatedConvolution<double, 3>& op,
    double cw_L,
    int max_d) {
    ResolvedOperatorPlan res;
    res.plan = &op.plan();
    res.stencil = &op.stencil_cache();

    const auto& opts = op.options();
    const bool cell_changed =
        std::abs(cw_L - op.cell_width()) > 1e-13 * std::max(1.0, op.cell_width());
    if (cell_changed || max_d > res.plan->n_max()) {
        const int max_allowed = static_cast<int>(MAXLEVEL) - 1;
        const int req_n_max = std::min(max_allowed, std::max(res.plan->n_max(), max_d + 2));
        WindowedQuadratureConfig cfg;
        cfg.k = op.k();
        cfg.eps = windowed_sumfact_detail::tighten_loose_windowed_eps(op.eps(), opts.eps_loc);
        cfg.mu_user = op.mu_user();
        cfg.cell_width_L = cw_L;
        cfg.max_tree_level = req_n_max;
        cfg.oversampling_M = opts.oversampling_M;
        cfg.gc_rank_Q = opts.gc_rank_Q;
        cfg.kind = op.kind();
        cfg.custom_coeffs_user = op.custom_coeffs_user();
        cfg.custom_expnts_user = op.custom_expnts_user();
        auto built = get_or_build_plan_and_stencil(cfg);
        res.ext_plan = std::move(built.first);
        res.ext_stencil = std::move(built.second);
        res.plan = res.ext_plan.get();
        res.stencil = res.ext_stencil.get();
    }
    return res;
}

template <typename Fn>
class WindowedBatchWorkerTask final : public PoolTaskInterface {
public:
    WindowedBatchWorkerTask(Fn* fn, std::atomic<int>* remaining)
        : PoolTaskInterface(TaskAttributes::hipri())
        , fn_(fn)
        , remaining_(remaining) {}

    ~WindowedBatchWorkerTask() override {
        if (remaining_ != nullptr) {
            remaining_->fetch_sub(1, std::memory_order_release);
        }
    }

    void run(const TaskThreadEnv&) override {
        (*fn_)();
    }

private:
    void get_id(std::pair<void*, unsigned short>& id) const override {
        PoolTaskInterface::make_id(id, &WindowedBatchWorkerTask::run);
    }

    Fn* fn_;
    std::atomic<int>* remaining_;
};

template <typename R>
Function<R, 3> apply_pure_spiky_same_world(World& target_world,
                                           const Function<R, 3>& f_same,
                                           const WindowedQuadraturePlan& plan,
                                           double thresh,
                                           bool fence,
                                           bool destructive_same_world) {
    mTxmq_init();
    Function<R, 3>& ff = const_cast<Function<R, 3>&>(f_same);
    const TreeState orig_st = ff.get_impl()->get_tree_state();
    const bool need_reconstruct =
        (orig_st != reconstructed && orig_st != nonstandard_with_leaves);
    if (need_reconstruct) {
        ff.reconstruct(true);
    }

    const int k = ff.k();
    const long kk = static_cast<long>(k);
    const double eff_thresh = (thresh > 0.0) ? thresh : ff.thresh();

    Function<R, 3> result =
        FunctionFactory<R, 3>(target_world).k(k).thresh(eff_thresh).pmap(ff.get_pmap()).empty();
    auto* res_impl = result.get_impl().get();
    auto& res_coeffs = res_impl->get_coeffs();
    res_coeffs.process_pending();

    struct LeafJob {
        Key<3> key;
        int level = 0;
        Tensor<R> s_in;
    };

    std::vector<Key<3>> internal_keys;
    std::vector<LeafJob> leaf_jobs;
    const auto& src_coeffs = ff.get_impl()->get_coeffs();
    internal_keys.reserve(src_coeffs.size());
    leaf_jobs.reserve(src_coeffs.size());
    for (auto it = src_coeffs.begin(); it != src_coeffs.end(); ++it) {
        const Key<3>& key = it->first;
        const FunctionNode<R, 3>& node = it->second;
        if (node.has_children()) {
            internal_keys.push_back(key);
        } else if (node.has_coeff()) {
            leaf_jobs.push_back(
                LeafJob{key, static_cast<int>(key.level()), node.coeff().full_tensor_copy()});
        } else {
            Tensor<R> zero_s(kk, kk, kk);
            leaf_jobs.push_back(LeafJob{key, static_cast<int>(key.level()), std::move(zero_s)});
        }
    }

    std::vector<Tensor<R>> leaf_outs(leaf_jobs.size());
    windowed_parallel_detail::parallel_for_chunks(
        leaf_jobs.size(),
        windowed_parallel_detail::kDefaultParallelGrain,
        [&](std::size_t b, std::size_t e_hi, std::size_t) {
            Tensor<R> tmp;
            Tensor<R> work1;
            Tensor<R> work2;
            for (std::size_t idx = b; idx < e_hi; ++idx) {
                Tensor<R> out(kk, kk, kk);
                plan.accumulate_gc_leaf_3d(
                    leaf_jobs[idx].level, leaf_jobs[idx].s_in, out, tmp, work1, work2);
                leaf_outs[idx] = std::move(out);
            }
        });

    const TensorArgs full_targs(0.0, TT_FULL);
    for (const Key<3>& key : internal_keys) {
        res_coeffs.replace(key, FunctionNode<R, 3>(GenTensor<R>(), true));
    }
    for (std::size_t idx = 0; idx < leaf_jobs.size(); ++idx) {
        GenTensor<R> coeff(leaf_outs[idx], full_targs);
        res_coeffs.replace(leaf_jobs[idx].key, FunctionNode<R, 3>(coeff, false));
    }
    res_impl->set_tree_state(reconstructed);

    if (destructive_same_world) {
        ff.clear(false);
    } else if (need_reconstruct && orig_st == compressed) {
        ff.compress(false);
    }
    if (fence) {
        target_world.gop.fence();
    }
    return result;
}

template <typename R>
Function<R, 3> apply_windowed_impl(World& target_world,
                                   const WindowedSeparatedConvolution<double, 3>& op,
                                   const Function<R, 3>& f,
                                   bool fence) {
    MADNESS_CHECK_THROW(f.is_initialized(),
                        "WindowedSeparatedConvolution: input Function must be initialized");
    MADNESS_CHECK_THROW(!f.is_on_demand(),
                        "WindowedSeparatedConvolution: on-demand input Function is not supported");
    MADNESS_CHECK_THROW(!FunctionDefaults<3>::get_bc().is_periodic_any(),
                        "WindowedSeparatedConvolution: periodic boundary conditions are not supported");
    MADNESS_CHECK_THROW(f.k() == op.k(),
                        "WindowedSeparatedConvolution: mismatch between Function::k() and operator k");

    const double cw_L = check_and_get_cubic_cell_width_3d();
    const auto& opts = op.options();

    auto run_on_same_world = [&](const Function<R, 3>& f_same,
                                 bool inner_fence,
                                 bool destructive_same_world) -> Function<R, 3> {
        int max_d = static_cast<int>(f_same.max_local_depth());
        if (target_world.size() > 1) {
            target_world.gop.max(max_d);
        }

        const ResolvedOperatorPlan resolved =
            resolve_plan_and_stencil_for_depth(op, cw_L, max_d);

        const double thresh = (opts.eps_loc > 0.0)
                                  ? opts.eps_loc
                                  : ((f_same.thresh() > 0.0) ? f_same.thresh() : op.eps());

        if (resolved.plan->is_purely_spiky_for_depth(max_d)) {
            return apply_pure_spiky_same_world<R>(
                target_world, f_same, *resolved.plan, thresh, inner_fence, destructive_same_world);
        }

        return apply_windowed_ns_3d<R>(target_world,
                                       f_same,
                                       *resolved.plan,
                                       *resolved.stencil,
                                       thresh,
                                       opts.enable_dynamic_refinement,
                                       opts.dynamic_refine_factor,
                                       inner_fence,
                                       destructive_same_world,
                                       opts.eps_loc);
    };

    if (&f.world() == &target_world) {
        return run_on_same_world(f, fence, op.destructive());
    }

    Function<R, 3> result;
    {
        Function<R, 3> f_sub = copy(target_world, f, true);
        result = run_on_same_world(f_sub, true, true);
    }
    if (op.destructive()) {
        const_cast<Function<R, 3>&>(f).clear(fence);
    } else if (fence) {
        target_world.gop.fence();
    }
    return result;
}

template <typename R>
std::vector<Function<R, 3>> apply_windowed_batch_impl(
    World& target_world,
    const std::vector<const WindowedSeparatedConvolution<double, 3>*>& ops,
    const std::vector<Function<R, 3>>& f,
    bool fence,
    bool ignore_op_destructive,
    bool force_destructive) {
    MADNESS_CHECK_THROW(
        ops.size() == f.size(),
        "WindowedSeparatedConvolution::apply_batch_on_world: operator/function vector size mismatch");
    const std::size_t N = f.size();
    if (N == 0) {
        if (fence) {
            target_world.gop.fence();
        }
        return {};
    }

    MADNESS_CHECK_THROW(!FunctionDefaults<3>::get_bc().is_periodic_any(),
                        "WindowedSeparatedConvolution: periodic boundary conditions are not supported");
    const double cw_L = check_and_get_cubic_cell_width_3d();

    for (std::size_t i = 0; i < N; ++i) {
        MADNESS_CHECK_THROW(ops[i] != nullptr,
                            "WindowedSeparatedConvolution: null operator pointer in batch apply");
        MADNESS_CHECK_THROW(f[i].is_initialized(),
                            "WindowedSeparatedConvolution: input Function must be initialized");
        MADNESS_CHECK_THROW(!f[i].is_on_demand(),
                            "WindowedSeparatedConvolution: on-demand input Function is not supported");
        MADNESS_CHECK_THROW(f[i].k() == ops[i]->k(),
                            "WindowedSeparatedConvolution: mismatch between Function::k() and operator k");
    }

    mTxmq_init();

    // Prepare working handles on `target_world` (copying across worlds only when needed).
    std::vector<Function<R, 3>> f_work(N);
    bool any_subworld_copy = false;
    for (std::size_t i = 0; i < N; ++i) {
        if (&f[i].world() == &target_world) {
            f_work[i] = f[i];
        } else {
            f_work[i] = copy(target_world, f[i], false);
            any_subworld_copy = true;
        }
    }
    if (any_subworld_copy) {
        target_world.gop.fence();
    }

    // Deduplicate `f_work` by `FunctionImpl*` so duplicate handles in `f` never race
    // or double-mutate tree state during batched Stage 1a / Stage 5b.
    struct UniqueInputEntry {
        FunctionImpl<R, 3>* impl = nullptr;
        Function<R, 3>* func = nullptr;
        bool was_compressed = false;
        bool destructive = false;
    };
    std::vector<UniqueInputEntry> unique_inputs;
    unique_inputs.reserve(N);

    struct SubworldClearEntry {
        FunctionImpl<R, 3>* impl = nullptr;
        Function<R, 3>* func = nullptr;
    };
    std::vector<SubworldClearEntry> subworld_clears;

    for (std::size_t i = 0; i < N; ++i) {
        const bool item_destructive =
            (!ignore_op_destructive) && (force_destructive || ops[i]->destructive());
        const bool is_subworld = (&f[i].world() != &target_world);
        const bool work_destructive = is_subworld ? true : item_destructive;

        FunctionImpl<R, 3>* impl_ptr = f_work[i].get_impl().get();
        auto it = std::find_if(unique_inputs.begin(),
                               unique_inputs.end(),
                               [impl_ptr](const UniqueInputEntry& e) { return e.impl == impl_ptr; });
        if (it == unique_inputs.end()) {
            unique_inputs.push_back(UniqueInputEntry{impl_ptr, &f_work[i], false, work_destructive});
        } else {
            it->destructive = it->destructive || work_destructive;
        }

        if (is_subworld && item_destructive) {
            FunctionImpl<R, 3>* orig_impl = f[i].get_impl().get();
            auto oit = std::find_if(subworld_clears.begin(),
                                    subworld_clears.end(),
                                    [orig_impl](const SubworldClearEntry& e) {
                                        return e.impl == orig_impl;
                                    });
            if (oit == subworld_clears.end()) {
                subworld_clears.push_back(
                    SubworldClearEntry{orig_impl, &const_cast<Function<R, 3>&>(f[i])});
            }
        }
    }

    // Stage 1a (batched across all unique inputs):
    // - On `target_world.size() == 1`, `reconstructed`, `compressed`, and `nonstandard_with_leaves`
    //   are ingested read-only with zero `WorldContainer` mutations.
    // - On `target_world.size() > 1`, prepare `nonstandard_with_leaves` across ranks.
    bool any_need_reconstruct = false;
    for (auto& ui : unique_inputs) {
        const TreeState st = ui.impl->get_tree_state();
        ui.was_compressed = (st == compressed);
        if (target_world.size() == 1) {
            if (st != reconstructed && st != compressed && st != nonstandard_with_leaves) {
                ui.func->reconstruct(false);
                any_need_reconstruct = true;
            }
        } else {
            if (st != reconstructed && st != nonstandard_with_leaves) {
                ui.func->reconstruct(false);
                any_need_reconstruct = true;
            }
        }
    }
    if (any_need_reconstruct) {
        target_world.gop.fence();
    }

    if (target_world.size() > 1) {
        bool any_need_ns = false;
        for (auto& ui : unique_inputs) {
            if (ui.impl->get_tree_state() != nonstandard_with_leaves) {
                ui.func->make_nonstandard(true, false);
                any_need_ns = true;
            }
        }
        if (any_need_ns) {
            target_world.gop.fence();
        }
    }

    // Pre-compute max tree depths across `target_world` in a single collective,
    // pre-warm `FunctionCommonData`, resolve plans/stencils, and allocate output Functions
    // deterministically in `0..N-1` order on the calling thread.
    std::vector<int> max_depths(N, 0);
    for (std::size_t i = 0; i < N; ++i) {
        max_depths[i] = static_cast<int>(f_work[i].max_local_depth());
    }
    if (target_world.size() > 1) {
        target_world.gop.max(max_depths.data(), N);
    }

    std::vector<ResolvedOperatorPlan> resolved_plans(N);
    std::vector<double> eff_thresh(N, 0.0);
    std::vector<double> eff_refine_fac(N, 1.0);
    std::vector<Function<R, 3>> result(N);

    for (std::size_t i = 0; i < N; ++i) {
        const auto& op_i = *ops[i];
        const auto& opts_i = op_i.options();
        (void)FunctionCommonData<R, 3>::get(op_i.k());
        resolved_plans[i] = resolve_plan_and_stencil_for_depth(op_i, cw_L, max_depths[i]);

        const double raw_thresh =
            (opts_i.eps_loc > 0.0)
                ? opts_i.eps_loc
                : ((f_work[i].thresh() > 0.0) ? f_work[i].thresh() : op_i.eps());
        eff_thresh[i] = (raw_thresh > 0.0) ? raw_thresh : f_work[i].thresh();
        eff_refine_fac[i] =
            (opts_i.dynamic_refine_factor > 0.0) ? opts_i.dynamic_refine_factor : 1.0;

        result[i] = FunctionFactory<R, 3>(target_world)
                        .k(op_i.k())
                        .thresh(eff_thresh[i])
                        .pmap(f_work[i].get_pmap())
                        .empty();
        result[i].get_impl()->get_coeffs().process_pending();
    }
    if (target_world.size() > 1) {
        target_world.gop.fence();
    }

    // Stages 1b..5a:
    // - On single-rank (`target_world.size() == 1`) with `N >= 2` and available `ThreadPool`
    //   workers outside nested batch tasks, dispatch `(f_work[idx], result[idx])` across
    //   pool threads + calling thread via atomic work-stealing under `WindowedBatchTaskScope`.
    // - On multi-rank MPI (`target_world.size() > 1`) or `N == 1`, execute sequentially in
    //   deterministic `0..N-1` order so `WindowedNSHaloExchange` collectives match across ranks.
    bool use_function_parallel = false;
#if !defined(HAVE_PARSEC)
    if (target_world.size() == 1 && N >= 2 && !g_in_windowed_batch_task && ThreadPool::size() > 0) {
        use_function_parallel = true;
    }
#endif

    if (use_function_parallel) {
#if !defined(HAVE_PARSEC)
        // LPT (Longest-Processing-Time-First) dispatch order: sort indices `0..N-1` by
        // descending local node count so deep nuclear-cusp trees (e.g., `O 1s`) start first
        // at `t = 0` while keeping `result[0..N-1]` in exact caller order.
        std::vector<std::size_t> node_counts(N, 0);
        std::vector<std::size_t> exec_order(N, 0);
        for (std::size_t i = 0; i < N; ++i) {
            node_counts[i] = f_work[i].get_impl()->get_coeffs().size();
            exec_order[i] = i;
        }
        std::stable_sort(exec_order.begin(),
                         exec_order.end(),
                         [&node_counts](std::size_t a, std::size_t b) {
                             return node_counts[a] > node_counts[b];
                         });

        std::atomic<std::size_t> next_slot{0};
        std::mutex exc_mtx;
        std::exception_ptr captured_exc;

#ifdef HAVE_INTEL_TBB
        const std::size_t total_threads = ThreadPool::size() + 1;
        const std::size_t n_workers = std::min<std::size_t>(N, total_threads);
        auto run_tbb_worker = [&]() noexcept {
            WindowedActiveBatchWorkerGuard active_guard;
            WindowedBatchTaskScope batch_scope;
            for (;;) {
                const std::size_t slot = next_slot.fetch_add(1, std::memory_order_relaxed);
                if (slot >= N) break;
                const std::size_t idx = exec_order[slot];
                try {
                    const auto& op_i = *ops[idx];
                    const auto& opts_i = op_i.options();
                    apply_windowed_ns_3d_prepared<R>(target_world,
                                                     f_work[idx],
                                                     result[idx],
                                                     *resolved_plans[idx].plan,
                                                     *resolved_plans[idx].stencil,
                                                     eff_thresh[idx],
                                                     opts_i.enable_dynamic_refinement,
                                                     eff_refine_fac[idx],
                                                     opts_i.eps_loc);
                } catch (...) {
                    std::lock_guard<std::mutex> lock(exc_mtx);
                    if (!captured_exc) {
                        captured_exc = std::current_exception();
                    }
                    next_slot.store(N, std::memory_order_relaxed);
                    break;
                }
            }
        };

        if (ThreadPool::tbb_arena) {
            ThreadPool::tbb_arena->execute([&]() {
                tbb::parallel_for(
                    tbb::blocked_range<std::size_t>(0, n_workers, 1),
                    [&](const tbb::blocked_range<std::size_t>& r) {
                        for (std::size_t w = r.begin(); w != r.end(); ++w) {
                            run_tbb_worker();
                        }
                    },
                    tbb::simple_partitioner());
            });
        } else {
            run_tbb_worker();
        }
#else
        std::atomic<int> remaining{0};

        struct PreReservedBatchWorkerRelease {
            ~PreReservedBatchWorkerRelease() {
                g_active_windowed_batch_workers.fetch_sub(1, std::memory_order_acq_rel);
            }
        };

        auto run_worker = [&]() noexcept {
            PreReservedBatchWorkerRelease worker_release;
            WindowedBatchTaskScope batch_scope;
            for (;;) {
                const std::size_t slot = next_slot.fetch_add(1, std::memory_order_relaxed);
                if (slot >= N) break;
                const std::size_t idx = exec_order[slot];
                try {
                    const auto& op_i = *ops[idx];
                    const auto& opts_i = op_i.options();
                    apply_windowed_ns_3d_prepared<R>(target_world,
                                                     f_work[idx],
                                                     result[idx],
                                                     *resolved_plans[idx].plan,
                                                     *resolved_plans[idx].stencil,
                                                     eff_thresh[idx],
                                                     opts_i.enable_dynamic_refinement,
                                                     eff_refine_fac[idx],
                                                     opts_i.eps_loc);
                } catch (...) {
                    std::lock_guard<std::mutex> lock(exc_mtx);
                    if (!captured_exc) {
                        captured_exc = std::current_exception();
                    }
                    next_slot.store(N, std::memory_order_relaxed);
                    break;
                }
            }
        };

        const std::size_t n_helpers = std::min<std::size_t>(N - 1, ThreadPool::size());
        using WorkerTask = WindowedBatchWorkerTask<decltype(run_worker)>;
        std::vector<WorkerTask*> helper_tasks(n_helpers, nullptr);
        for (std::size_t w = 0; w < n_helpers; ++w) {
            helper_tasks[w] = new WorkerTask(&run_worker, &remaining);
        }

        const int total_batch_workers = static_cast<int>(n_helpers + 1);
        g_active_windowed_batch_workers.fetch_add(total_batch_workers, std::memory_order_acq_rel);
        remaining.store(static_cast<int>(n_helpers), std::memory_order_relaxed);
        for (std::size_t w = 0; w < n_helpers; ++w) {
            ThreadPool::add(helper_tasks[w]);
        }
        thread_purge();

        run_worker();

        ThreadPool::await(
            [&remaining]() { return remaining.load(std::memory_order_acquire) == 0; }, true);
#endif

        if (captured_exc) {
            std::rethrow_exception(captured_exc);
        }
#endif
    } else {
        for (std::size_t i = 0; i < N; ++i) {
            const auto& op_i = *ops[i];
            const auto& opts_i = op_i.options();
            apply_windowed_ns_3d_prepared<R>(target_world,
                                             f_work[i],
                                             result[i],
                                             *resolved_plans[i].plan,
                                             *resolved_plans[i].stencil,
                                             eff_thresh[i],
                                             opts_i.enable_dynamic_refinement,
                                             eff_refine_fac[i],
                                             opts_i.eps_loc);
        }
    }

    // Stage 5b (batched): restore unique input Functions and reconstruct any output Functions
    // that used `nonstandard_after_apply` (multi-rank or dynamic leaf/virtual promotion).
    for (auto& ui : unique_inputs) {
        restore_windowed_ns_input_3d(*ui.func, ui.was_compressed, ui.destructive, false);
    }
    for (std::size_t i = 0; i < N; ++i) {
        if (result[i].get_impl()->get_tree_state() != reconstructed) {
            result[i].reconstruct(false);
        }
    }
    if (fence || any_subworld_copy) {
        target_world.gop.fence();
    }

    for (auto& sc : subworld_clears) {
        sc.func->clear(fence);
    }

    return result;
}

} // namespace

template <typename Q, std::size_t NDIM>
WindowedSeparatedConvolution<Q, NDIM>::WindowedSeparatedConvolution(
    World& world,
    WindowedKernelKind kind,
    double mu_user,
    double eps,
    int k,
    const WindowedConvolutionOptions& opts)
    : world_(&world)
    , kind_(kind)
    , mu_user_((kind == WindowedKernelKind::Coulomb) ? 0.0 : mu_user)
    , eps_(eps)
    , k_(k)
    , cell_width_L_(check_and_get_cubic_cell_width_3d())
    , opts_(opts)
    , destructive_(false) {
    MADNESS_CHECK_THROW(k_ >= 1,
                        "WindowedSeparatedConvolution: wavelet order k must be >= 1");
    MADNESS_CHECK_THROW(eps_ > 0.0 && eps_ < 1.0,
                        "WindowedSeparatedConvolution: accuracy eps must lie in (0, 1)");
    MADNESS_CHECK_THROW(mu_user >= 0.0,
                        "WindowedSeparatedConvolution: screening parameter mu must be >= 0");
    MADNESS_CHECK_THROW(kind_ != WindowedKernelKind::Coulomb || mu_user == 0.0,
                        "WindowedSeparatedConvolution: Coulomb kernel requires mu == 0");

    const int max_allowed = static_cast<int>(MAXLEVEL) - 1;
    const int n_max = std::clamp(std::max(1, opts_.max_level), 1, max_allowed);

    WindowedQuadratureConfig cfg;
    cfg.k = k_;
    cfg.eps = windowed_sumfact_detail::tighten_loose_windowed_eps(eps_, opts_.eps_loc);
    cfg.mu_user = mu_user_;
    cfg.cell_width_L = cell_width_L_;
    cfg.max_tree_level = n_max;
    cfg.oversampling_M = opts_.oversampling_M;
    cfg.gc_rank_Q = opts_.gc_rank_Q;
    cfg.kind = kind_;

    auto built = get_or_build_plan_and_stencil(cfg);
    plan_ = std::move(built.first);
    stencil_cache_ = std::move(built.second);
}

template <typename Q, std::size_t NDIM>
WindowedSeparatedConvolution<Q, NDIM>::WindowedSeparatedConvolution(
    World& world,
    const Tensor<double>& coeff,
    const Tensor<double>& expnt,
    double eps,
    int k,
    const WindowedConvolutionOptions& opts)
    : world_(&world)
    , kind_(WindowedKernelKind::CustomGaussian)
    , mu_user_(0.0)
    , eps_(eps)
    , k_(k)
    , cell_width_L_(check_and_get_cubic_cell_width_3d())
    , opts_(opts)
    , destructive_(false) {
    MADNESS_CHECK_THROW(k_ >= 1,
                        "WindowedSeparatedConvolution: wavelet order k must be >= 1");
    MADNESS_CHECK_THROW(eps_ > 0.0 && eps_ < 1.0,
                        "WindowedSeparatedConvolution: accuracy eps must lie in (0, 1)");
    MADNESS_CHECK_THROW(coeff.ndim() == 1 && expnt.ndim() == 1 &&
                            coeff.size() > 0 && coeff.size() == expnt.size(),
                        "WindowedSeparatedConvolution: custom Gaussian coeff/expnt must be matching non-empty 1D tensors");

    const long n_terms = coeff.size();
    custom_coeffs_user_.resize(static_cast<std::size_t>(n_terms));
    custom_expnts_user_.resize(static_cast<std::size_t>(n_terms));
    for (long idx = 0; idx < n_terms; ++idx) {
        const double c_val = coeff(idx);
        const double e_val = expnt(idx);
        MADNESS_CHECK_THROW(e_val > 0.0 && std::isfinite(e_val) && std::isfinite(c_val),
                            "WindowedSeparatedConvolution: custom Gaussian exponents must be finite and > 0");
        custom_coeffs_user_[static_cast<std::size_t>(idx)] = c_val;
        custom_expnts_user_[static_cast<std::size_t>(idx)] = e_val;
    }

    const int max_allowed = static_cast<int>(MAXLEVEL) - 1;
    const int n_max = std::clamp(std::max(1, opts_.max_level), 1, max_allowed);

    WindowedQuadratureConfig cfg;
    cfg.k = k_;
    cfg.eps = windowed_sumfact_detail::tighten_loose_windowed_eps(eps_, opts_.eps_loc);
    cfg.mu_user = 0.0;
    cfg.cell_width_L = cell_width_L_;
    cfg.max_tree_level = n_max;
    cfg.oversampling_M = opts_.oversampling_M;
    cfg.gc_rank_Q = opts_.gc_rank_Q;
    cfg.kind = WindowedKernelKind::CustomGaussian;
    cfg.custom_coeffs_user = custom_coeffs_user_;
    cfg.custom_expnts_user = custom_expnts_user_;

    auto built = get_or_build_plan_and_stencil(cfg);
    plan_ = std::move(built.first);
    stencil_cache_ = std::move(built.second);
}

template <typename Q, std::size_t NDIM>
Function<double, 3>
WindowedSeparatedConvolution<Q, NDIM>::apply_on_world(
    World& world, const Function<double, 3>& f, bool fence) const {
    return apply_windowed_impl<double>(world, *this, f, fence);
}

template <typename Q, std::size_t NDIM>
Function<double_complex, 3>
WindowedSeparatedConvolution<Q, NDIM>::apply_on_world(
    World& world, const Function<double_complex, 3>& f, bool fence) const {
    return apply_windowed_impl<double_complex>(world, *this, f, fence);
}

template <typename Q, std::size_t NDIM>
Function<double, 3>
WindowedSeparatedConvolution<Q, NDIM>::apply_func(
    const Function<double, 3>& f, bool fence) const {
    World& target_world = f.is_initialized() ? f.world() : *world_;
    return apply_windowed_impl<double>(target_world, *this, f, fence);
}

template <typename Q, std::size_t NDIM>
Function<double_complex, 3>
WindowedSeparatedConvolution<Q, NDIM>::apply_func(
    const Function<double_complex, 3>& f, bool fence) const {
    World& target_world = f.is_initialized() ? f.world() : *world_;
    return apply_windowed_impl<double_complex>(target_world, *this, f, fence);
}

template <typename Q, std::size_t NDIM>
std::vector<Function<double, 3>>
WindowedSeparatedConvolution<Q, NDIM>::apply_batch_on_world(
    World& world,
    const std::vector<const WindowedSeparatedConvolution<Q, NDIM>*>& ops,
    const std::vector<Function<double, 3>>& f,
    bool fence,
    bool ignore_op_destructive,
    bool force_destructive) {
    static_assert(std::is_same_v<Q, double> && NDIM == 3,
                  "WindowedSeparatedConvolution::apply_batch_on_world requires Q = double and NDIM = 3.");
    return apply_windowed_batch_impl<double>(
        world, ops, f, fence, ignore_op_destructive, force_destructive);
}

template <typename Q, std::size_t NDIM>
std::vector<Function<double_complex, 3>>
WindowedSeparatedConvolution<Q, NDIM>::apply_batch_on_world(
    World& world,
    const std::vector<const WindowedSeparatedConvolution<Q, NDIM>*>& ops,
    const std::vector<Function<double_complex, 3>>& f,
    bool fence,
    bool ignore_op_destructive,
    bool force_destructive) {
    static_assert(std::is_same_v<Q, double> && NDIM == 3,
                  "WindowedSeparatedConvolution::apply_batch_on_world requires Q = double and NDIM = 3.");
    return apply_windowed_batch_impl<double_complex>(
        world, ops, f, fence, ignore_op_destructive, force_destructive);
}

template class WindowedSeparatedConvolution<double, 3>;

namespace windowed_ns_bridge {

bool env_opt_in_enabled() noexcept {
    const char* val = std::getenv("MADNESS_WINDOWED_NS");
    if (val == nullptr || val[0] == '\0') {
        return false;
    }
    if (std::strcmp(val, "1") == 0 ||
        std::strcmp(val, "true") == 0 ||
        std::strcmp(val, "TRUE") == 0 ||
        std::strcmp(val, "on") == 0 ||
        std::strcmp(val, "ON") == 0 ||
        std::strcmp(val, "yes") == 0 ||
        std::strcmp(val, "YES") == 0) {
        return true;
    }
    return false;
}

std::shared_ptr<const WindowedSeparatedConvolution<double, 3>>
get_or_make(World& world,
            std::shared_ptr<const WindowedSeparatedConvolution<double, 3>>& slot,
            int op_type,
            double mu,
            double thresh,
            int k,
            const Tensor<double>& stored_coeffs,
            const Tensor<double>& stored_exponents) {
    if (!(thresh > 0.0 && thresh < 1.0) || k < 1) {
        return nullptr;
    }
    const Tensor<double>& cw = FunctionDefaults<3>::get_cell_width();
    const double cw_L = cw(0L);
    if (!(cw_L > 0.0)) {
        return nullptr;
    }
    if (std::abs(cw(1L) - cw_L) > 1e-13 * cw_L ||
        std::abs(cw(2L) - cw_L) > 1e-13 * cw_L) {
        return nullptr;
    }

    const bool is_coulomb_or_bsh =
        (op_type == static_cast<int>(OT_G12) || op_type == static_cast<int>(OT_BSH));
    if (is_coulomb_or_bsh) {
        if (mu < 0.0) {
            return nullptr;
        }
        const WindowedKernelKind kind =
            (op_type == static_cast<int>(OT_G12)) ? WindowedKernelKind::Coulomb
                                                  : WindowedKernelKind::BSH;
        const double eff_mu = (kind == WindowedKernelKind::Coulomb) ? 0.0 : mu;

        static std::mutex bridge_mtx;
        std::lock_guard<std::mutex> lock(bridge_mtx);
        if (slot &&
            slot->kind() == kind &&
            slot->k() == k &&
            slot->eps() == thresh &&
            slot->mu_user() == eff_mu &&
            std::abs(slot->cell_width() - cw_L) <= 1e-13 * std::max(1.0, cw_L)) {
            return slot;
        }

        slot = std::make_shared<const WindowedSeparatedConvolution<double, 3>>(
            world, kind, eff_mu, thresh, k);
        return slot;
    }

    if (stored_coeffs.ndim() != 1 || stored_exponents.ndim() != 1 ||
        stored_coeffs.size() <= 0 || stored_coeffs.size() != stored_exponents.size()) {
        return nullptr;
    }
    const long n_terms = stored_coeffs.size();
    for (long idx = 0; idx < n_terms; ++idx) {
        if (!(stored_exponents(idx) > 0.0) || !std::isfinite(stored_exponents(idx)) ||
            !std::isfinite(stored_coeffs(idx))) {
            return nullptr;
        }
    }

    static std::mutex custom_bridge_mtx;
    std::lock_guard<std::mutex> lock(custom_bridge_mtx);
    if (slot &&
        slot->kind() == WindowedKernelKind::CustomGaussian &&
        slot->k() == k &&
        slot->eps() == thresh &&
        std::abs(slot->cell_width() - cw_L) <= 1e-13 * std::max(1.0, cw_L) &&
        static_cast<long>(slot->custom_coeffs_user().size()) == n_terms) {
        bool same = true;
        for (long idx = 0; idx < n_terms; ++idx) {
            if (slot->custom_coeffs_user()[static_cast<std::size_t>(idx)] != stored_coeffs(idx) ||
                slot->custom_expnts_user()[static_cast<std::size_t>(idx)] != stored_exponents(idx)) {
                same = false;
                break;
            }
        }
        if (same) {
            return slot;
        }
    }

    slot = std::make_shared<const WindowedSeparatedConvolution<double, 3>>(
        world, stored_coeffs, stored_exponents, thresh, k);
    return slot;
}

Function<double, 3>
apply_on_world(World& world,
               const WindowedSeparatedConvolution<double, 3>& wop,
               const Function<double, 3>& f,
               bool fence) {
    return wop.apply_on_world(world, f, fence);
}

Function<double_complex, 3>
apply_on_world(World& world,
               const WindowedSeparatedConvolution<double, 3>& wop,
               const Function<double_complex, 3>& f,
               bool fence) {
    return wop.apply_on_world(world, f, fence);
}

std::vector<Function<double, 3>>
apply_batch_single_op(World& world,
                      const WindowedSeparatedConvolution<double, 3>& wop,
                      const std::vector<Function<double, 3>>& f,
                      bool fence,
                      bool destructive_input) {
    std::vector<const WindowedSeparatedConvolution<double, 3>*> ops(f.size(), &wop);
    return WindowedSeparatedConvolution<double, 3>::apply_batch_on_world(
        world, ops, f, fence, false, destructive_input);
}

std::vector<Function<double_complex, 3>>
apply_batch_single_op(World& world,
                      const WindowedSeparatedConvolution<double, 3>& wop,
                      const std::vector<Function<double_complex, 3>>& f,
                      bool fence,
                      bool destructive_input) {
    std::vector<const WindowedSeparatedConvolution<double, 3>*> ops(f.size(), &wop);
    return WindowedSeparatedConvolution<double, 3>::apply_batch_on_world(
        world, ops, f, fence, false, destructive_input);
}

std::vector<Function<double, 3>>
apply_batch_multi_op(World& world,
                     const std::vector<const WindowedSeparatedConvolution<double, 3>*>& wops,
                     const std::vector<Function<double, 3>>& f,
                     bool fence) {
    return WindowedSeparatedConvolution<double, 3>::apply_batch_on_world(
        world, wops, f, fence, true, false);
}

std::vector<Function<double_complex, 3>>
apply_batch_multi_op(World& world,
                     const std::vector<const WindowedSeparatedConvolution<double, 3>*>& wops,
                     const std::vector<Function<double_complex, 3>>& f,
                     bool fence) {
    return WindowedSeparatedConvolution<double, 3>::apply_batch_on_world(
        world, wops, f, fence, true, false);
}

} // namespace windowed_ns_bridge

} // namespace madness
