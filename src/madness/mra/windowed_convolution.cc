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

inline double select_operator_gamma_spiky_override(int k) noexcept {
    // At low wavelet order `k < 8` (e.g. loose SCF protocol `1e-4, k = 6`), floor Stage 2's
    // box-local moment threshold `gamma_spiky` at the `k = 8` stiffness (`172.8`, keeping `W = 4`,
    // `M = 2`) so `erfc` boundary leakage across leaf faces stays < 1e-11 without `W = 5` bloat.
    return (k < 8) ? 172.8 : 0.0;
}

std::size_t get_multi_rank_batch_tile() noexcept {
    static const std::size_t tile = []() -> std::size_t {
        constexpr std::size_t kDefaultTile = 16;
        const char* s = std::getenv("MADNESS_WINDOWED_MPI_BATCH_TILE");
        if (s == nullptr || *s == '\0') return kDefaultTile;
        const int v = std::atoi(s);
        if (v <= 0) return kDefaultTile;
        return static_cast<std::size_t>(std::clamp(v, 1, 256));
    }();
    return tile;
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
    double gamma_spiky_override = 0.0;
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
            entry.gamma_spiky_override == cfg.gamma_spiky_override &&
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
        cfg.gamma_spiky_override,
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
        cfg.gamma_spiky_override = select_operator_gamma_spiky_override(op.k());
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

    Function<R, 3> result;
    {
        windowed_parallel_detail::SingleRankWorldObjLock obj_lock(target_world);
        result = FunctionFactory<R, 3>(target_world)
                     .k(k)
                     .thresh(eff_thresh)
                     .truncate_mode(ff.get_impl()->get_truncate_mode())
                     .pmap(ff.get_pmap())
                     .empty();
    }
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
        target_world,
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

double compute_effective_apply_thresh(
    const WindowedSeparatedConvolution<double, 3>& op,
    double f_thresh) noexcept {
    const auto& opts = op.options();
    if (opts.eps_loc > 0.0) {
        return opts.eps_loc;
    }
    if (f_thresh > 0.0) {
        return f_thresh;
    }
    return op.eps();
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
        int max_d = (target_world.size() == 1)
                        ? compute_serial_max_tree_level(f_same)
                        : static_cast<int>(f_same.max_local_depth());
        if (target_world.size() > 1) {
            target_world.gop.max(max_d);
        }

        const ResolvedOperatorPlan resolved =
            resolve_plan_and_stencil_for_depth(op, cw_L, max_d);

        const double thresh = compute_effective_apply_thresh(op, f_same.thresh());

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

    // On multi-rank MPI (`target_world.size() > 1`), bound simultaneous nonstandard input/output
    // tree memory (`nonstandard_with_leaves` + `nonstandard_after_apply`) to tiles of `<= 16`
    // functions by default (overridable via `MADNESS_WINDOWED_MPI_BATCH_TILE` in `[1, 256]`)
    // whenever all `N` input handles refer to distinct `FunctionImpl` instances, amortizing
    // lockstep halo & restriction collectives across orbitals in the tile.
    const std::size_t mpi_batch_tile = get_multi_rank_batch_tile();
    if (target_world.size() > 1 && N > mpi_batch_tile) {
        bool all_distinct_inputs = true;
        for (std::size_t i = 0; i < N && all_distinct_inputs; ++i) {
            const auto* pi = f[i].get_impl().get();
            for (std::size_t j = 0; j < i; ++j) {
                if (f[j].get_impl().get() == pi) {
                    all_distinct_inputs = false;
                    break;
                }
            }
        }
        if (all_distinct_inputs) {
            std::vector<Function<R, 3>> tiled_res;
            tiled_res.reserve(N);
            for (std::size_t ilo = 0; ilo < N; ilo += mpi_batch_tile) {
                const std::size_t ihi = std::min(ilo + mpi_batch_tile, N);
                const bool is_last_tile = (ihi == N);
                std::vector<const WindowedSeparatedConvolution<double, 3>*> sub_ops(
                    ops.begin() + static_cast<std::ptrdiff_t>(ilo),
                    ops.begin() + static_cast<std::ptrdiff_t>(ihi));
                std::vector<Function<R, 3>> sub_f(
                    f.begin() + static_cast<std::ptrdiff_t>(ilo),
                    f.begin() + static_cast<std::ptrdiff_t>(ihi));
                auto sub_out = apply_windowed_batch_impl<R>(target_world,
                                                            sub_ops,
                                                            sub_f,
                                                            is_last_tile ? fence : true,
                                                            ignore_op_destructive,
                                                            force_destructive);
                for (auto& fn_out : sub_out) {
                    tiled_res.push_back(std::move(fn_out));
                }
            }
            return tiled_res;
        }
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
    //   are ingested read-only with zero `WorldContainer` mutations (and compressed union-skeleton
    //   bloat is pruned on the fly inside `ingest_local_real_nodes_serial_compressed`).
    // - On `target_world.size() > 1`, pre-prune union-skeleton bloat on `TreeState::compressed`
    //   inputs (`W_rem(q) * ||d^Q||_F <= tau_keep(q)`) into slim compressed working copies in ONE
    //   coalesced ancestor-key exchange before `reconstruct()` + `make_nonstandard()`.
    if (target_world.size() > 1) {
        struct CompPruneCand {
            std::size_t ui_idx = 0;
            int k = 0;
            double eff_tol = 0.0;
            const WindowedQuadraturePlan* base_plan = nullptr;
            std::vector<Key<3>> local_seeds;
            long local_orig_int = 0;
        };
        std::vector<CompPruneCand> cands;
        for (std::size_t u = 0; u < unique_inputs.size(); ++u) {
            auto& ui = unique_inputs[u];
            if (ui.impl->get_tree_state() != compressed) continue;

            bool all_dyn = true;
            double min_eff_tol = 1e300;
            const WindowedQuadraturePlan* bplan = nullptr;
            int k_u = ui.func->k();
            for (std::size_t i = 0; i < N; ++i) {
                if (f_work[i].get_impl().get() != ui.impl) continue;
                const auto& op_i = *ops[i];
                const auto& opts_i = op_i.options();
                if (!opts_i.enable_dynamic_refinement) {
                    all_dyn = false;
                    break;
                }
                const double eth = compute_effective_apply_thresh(op_i, f_work[i].thresh());
                const double min_pt = (eth > 0.0 && opts_i.eps_loc <= 0.0)
                                          ? std::max(op_i.eps(), 0.1 * eth)
                                          : op_i.eps();
                const double rtol = (eth > 0.0) ? std::min(min_pt, eth) : op_i.eps();
                const double etol =
                    windowed_sumfact_detail::tighten_loose_windowed_eps(rtol, opts_i.eps_loc);
                min_eff_tol = std::min(min_eff_tol, etol);
                if (bplan == nullptr) bplan = &op_i.plan();
            }
            if (!all_dyn || bplan == nullptr || !(min_eff_tol > 0.0 && min_eff_tol < 1.0)) {
                continue;
            }

            CompPruneCand cp;
            cp.ui_idx = u;
            cp.k = k_u;
            cp.eff_tol = min_eff_tol;
            cp.base_plan = bplan;

            const long twok = 2L * static_cast<long>(k_u);
            const auto& coeffs = ui.impl->get_coeffs();
            for (auto it = coeffs.begin(); it != coeffs.end(); ++it) {
                const Key<3>& key = it->first;
                const int j = static_cast<int>(key.level());
                const auto& node = it->second;
                bool has_2k = false;
                double dn = 0.0;
                if (node.has_coeff()) {
                    const auto& gen_c = node.coeff();
                    if (gen_c.has_data() && gen_c.ndim() == 3 && gen_c.dim(0) == twok) {
                        has_2k = true;
                        const bool is_contig =
                            (gen_c.tensor_type() == TT_FULL && gen_c.get_tensor().iscontiguous());
                        const Tensor<R> full_c =
                            is_contig ? gen_c.get_tensor() : gen_c.full_tensor_copy();
                        dn = windowed_halo_detail::wavelet_octants_normf(full_c, k_u);
                    }
                }
                const bool is_exp = node.has_children() || (has_2k && (j > 0 || dn > 0.0));
                if (!is_exp) continue;
                ++cp.local_orig_int;

                if (j <= 1) {
                    cp.local_seeds.push_back(key);
                    continue;
                }
                const double w_bound = bplan->subtree_prune_op_bound(j);
                const double loc_t =
                    ui.impl->truncate_tol(min_eff_tol, Key<3>(static_cast<Level>(j)));
                const double dtol = (w_bound > 0.0) ? ((0.02 * loc_t) / w_bound) : 1e300;
                if (dn > dtol) {
                    cp.local_seeds.push_back(key);
                }
            }
            cands.push_back(std::move(cp));
        }

        if (!cands.empty()) {
            const std::size_t M_c = cands.size();
            std::vector<long> counts_gop(2 * M_c, 0L);
            for (std::size_t c = 0; c < M_c; ++c) {
                counts_gop[2 * c + 0] = cands[c].local_orig_int;
                counts_gop[2 * c + 1] = static_cast<long>(cands[c].local_seeds.size());
            }
            target_world.gop.sum(counts_gop.data(), counts_gop.size());

            std::vector<std::size_t> active_cands;
            for (std::size_t c = 0; c < M_c; ++c) {
                const long orig_int = counts_gop[2 * c + 0];
                const long seed_int = counts_gop[2 * c + 1];
                if (orig_int >= 32 && seed_int * 10L <= orig_int * 8L) {
                    active_cands.push_back(c);
                }
            }

            if (!active_cands.empty()) {
                const std::size_t A_c = active_cands.size();
                const ProcessID me = target_world.rank();
                const int nproc = target_world.size();

                std::vector<std::vector<Key<3>>> kept_local_keys(A_c);
                std::vector<Function<R, 3>> slim_funcs(A_c);

                {
                    windowed_halo_detail::WindowedNSHaloExchange<R> prune_ex(target_world, A_c);
                    std::vector<std::vector<windowed_halo_detail::HaloCandidateMsg>> anc_outbox(
                        static_cast<std::size_t>(nproc));
                    std::vector<std::vector<Key<3>>> per_rank_keys(static_cast<std::size_t>(nproc));

                    for (std::size_t ai = 0; ai < A_c; ++ai) {
                        const CompPruneCand& cp = cands[active_cands[ai]];
                        const auto& ui = unique_inputs[cp.ui_idx];
                        const auto pmap = ui.func->get_pmap();
                        for (int p = 0; p < nproc; ++p) per_rank_keys[static_cast<std::size_t>(p)].clear();

                        for (const Key<3>& skey : cp.local_seeds) {
                            Key<3> cur = skey;
                            for (;;) {
                                const ProcessID owner = pmap->owner(cur);
                                per_rank_keys[static_cast<std::size_t>(owner)].push_back(cur);
                                if (cur.level() == 0) break;
                                cur = cur.parent();
                            }
                        }

                        const std::uint16_t slot = static_cast<std::uint16_t>(ai);
                        for (int p = 0; p < nproc; ++p) {
                            auto& kvec = per_rank_keys[static_cast<std::size_t>(p)];
                            if (kvec.empty()) continue;
                            windowed_halo_detail::sort_unique_keys(kvec);
                            if (p == me) {
                                kept_local_keys[ai] = std::move(kvec);
                            } else {
                                auto& dst = anc_outbox[static_cast<std::size_t>(p)];
                                for (const Key<3>& k_anc : kvec) {
                                    dst.push_back(
                                        windowed_halo_detail::HaloCandidateMsg{k_anc, 1u, slot});
                                }
                            }
                        }

                        {
                            windowed_parallel_detail::SingleRankWorldObjLock obj_lock(target_world);
                            slim_funcs[ai] = FunctionFactory<R, 3>(target_world)
                                                 .k(cp.k)
                                                 .thresh(ui.func->thresh())
                                                 .truncate_mode(ui.impl->get_truncate_mode())
                                                 .pmap(pmap)
                                                 .empty();
                        }
                        slim_funcs[ai].get_impl()->get_coeffs().process_pending();
                    }

                    for (int p = 0; p < nproc; ++p) {
                        auto& batch = anc_outbox[static_cast<std::size_t>(p)];
                        if (batch.empty()) continue;
                        prune_ex.send(
                            static_cast<ProcessID>(p),
                            &windowed_halo_detail::WindowedNSHaloExchange<R>::recv_halo_candidates,
                            batch);
                    }
                    target_world.gop.fence();

                    for (std::size_t ai = 0; ai < A_c; ++ai) {
                        auto& kvec = kept_local_keys[ai];
                        for (const auto& msg : prune_ex.take_halo_candidates(ai)) {
                            kvec.push_back(msg.key);
                        }
                        windowed_halo_detail::sort_unique_keys(kvec);

                        const CompPruneCand& cp = cands[active_cands[ai]];
                        auto& ui = unique_inputs[cp.ui_idx];
                        const int k_u = cp.k;
                        const long kk = static_cast<long>(k_u);
                        const long kk3 = kk * kk * kk;
                        const long twok = 2L * kk;
                        const long twok3 = twok * twok * twok;
                        const Slice s0(0, kk - 1);
                        const auto& cdata = FunctionCommonData<R, 3>::get(k_u);
                        const Key<3>& key0 = cdata.key0;
                        const auto pmap = ui.func->get_pmap();
                        auto& src_coeffs = ui.impl->get_coeffs();
                        auto* slim_impl = slim_funcs[ai].get_impl().get();
                        auto& slim_coeffs = slim_impl->get_coeffs();
                        const TensorArgs full_targs(-1.0, TT_FULL);

                        if (pmap->owner(key0) == me) {
                            const bool root_kept_int = std::binary_search(
                                kvec.begin(),
                                kvec.end(),
                                key0,
                                windowed_halo_detail::windowed_key_less);
                            auto rit = src_coeffs.find(key0).get();
                            Tensor<R> c0;
                            if (rit != src_coeffs.end() && rit->second.has_coeff()) {
                                const auto& gc0 = rit->second.coeff();
                                if (gc0.has_data()) {
                                    c0 = (gc0.tensor_type() == TT_FULL &&
                                          gc0.get_tensor().iscontiguous())
                                             ? gc0.get_tensor()
                                             : gc0.full_tensor_copy();
                                }
                            }
                            if (root_kept_int) {
                                if (rit != src_coeffs.end() && c0.size() == twok3) {
                                    slim_coeffs.replace(
                                        key0, FunctionNode<R, 3>(rit->second.coeff(), true));
                                } else {
                                    Tensor<R> ns0(twok, twok, twok);
                                    if (c0.size() > 0) {
                                        ns0(s0, s0, s0) =
                                            windowed_halo_detail::extract_scaling_cube(c0, k_u);
                                    }
                                    slim_coeffs.replace(
                                        key0, FunctionNode<R, 3>(GenTensor<R>(ns0, full_targs), true));
                                }
                            } else {
                                Tensor<R> s0_cube =
                                    (c0.size() == kk3 && c0.iscontiguous())
                                        ? c0
                                        : ((c0.size() > 0)
                                               ? windowed_halo_detail::extract_scaling_cube(c0, k_u)
                                               : Tensor<R>(kk, kk, kk));
                                slim_coeffs.replace(
                                    key0,
                                    FunctionNode<R, 3>(GenTensor<R>(s0_cube, full_targs), false));
                            }
                        }

                        for (const Key<3>& k_int : kvec) {
                            if (k_int.level() == 0) continue;
                            auto it = src_coeffs.find(k_int).get();
                            if (it != src_coeffs.end() && it->second.has_coeff() &&
                                it->second.coeff().has_data() &&
                                it->second.coeff().dim(0) == twok) {
                                slim_coeffs.replace(
                                    k_int, FunctionNode<R, 3>(it->second.coeff(), true));
                            } else {
                                slim_coeffs.replace(k_int, FunctionNode<R, 3>(GenTensor<R>(), true));
                            }
                        }
                        slim_impl->set_tree_state(compressed);

                        FunctionImpl<R, 3>* old_impl = ui.impl;
                        if (ui.destructive) {
                            ui.func->clear(false);
                        }
                        std::size_t first_match = SIZE_MAX;
                        for (std::size_t i = 0; i < N; ++i) {
                            if (f_work[i].get_impl().get() == old_impl) {
                                f_work[i] = slim_funcs[ai];
                                if (first_match == SIZE_MAX) first_match = i;
                            }
                        }
                        MADNESS_CHECK(first_match != SIZE_MAX);
                        ui.func = &f_work[first_match];
                        ui.impl = slim_impl;
                        ui.was_compressed = false;
                        ui.destructive = true;
                    }

                    target_world.gop.fence();
                }
            }
        }
    }

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
        max_depths[i] = (target_world.size() == 1)
                            ? compute_serial_max_tree_level(f_work[i])
                            : static_cast<int>(f_work[i].max_local_depth());
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

        const double raw_thresh = compute_effective_apply_thresh(op_i, f_work[i].thresh());
        eff_thresh[i] = (raw_thresh > 0.0) ? raw_thresh : f_work[i].thresh();
        eff_refine_fac[i] =
            (opts_i.dynamic_refine_factor > 0.0) ? opts_i.dynamic_refine_factor : 1.0;

        {
            windowed_parallel_detail::SingleRankWorldObjLock obj_lock(target_world);
            result[i] = FunctionFactory<R, 3>(target_world)
                            .k(op_i.k())
                            .thresh(eff_thresh[i])
                            .truncate_mode(f_work[i].get_impl()->get_truncate_mode())
                            .pmap(f_work[i].get_pmap())
                            .empty();
        }
        result[i].get_impl()->get_coeffs().process_pending();
    }
    if (target_world.size() > 1) {
        target_world.gop.fence();
    }

    // Stages 1b..5a:
    // - Route single-rank vector batches (`target_world.size() == 1`) directly onto MADNESS's
    //   `WorldTaskQueue` / `tbb_arena` (`windowed_for_each_outer_tasks`) in LPT order:
    //     * Small-tree batches (`max_nodes < 512`, e.g., `HF Exchange` pair batches): submit all
    //       `N` functions as outer MADNESS tasks so fast P-cores and slower E-cores drain both
    //       outer pair tasks and inner `madness::Range` sub-tasks from the shared task queue.
    //     * Medium-tree batches (`512 <= max_nodes < 2048`): execute in memory-bounded 2-function
    //       tiles to overlap stage barriers while bounding `WindowedNSHaloWorkspace` RSS.
    //     * Large-tree batches (`max_nodes >= 2048`, e.g., `Vpsi` in `Apply BSH`): execute one
    //       orbital at a time (`1x` Peak RSS) while `8x` inner `madness::Range` over-decomposition
    //       keeps all heterogeneous/dynamically clocked cores saturated.
    std::size_t outer_batch_tile = 1u;
    std::vector<std::size_t> exec_order;
    if (target_world.size() == 1 && N >= 2 && !g_in_windowed_batch_task && ThreadPool::size() > 0) {
        static const int env_batch_workers = []() {
            const char* s = std::getenv("MADNESS_WINDOWED_BATCH_WORKERS");
            if (s == nullptr || *s == '\0') return 0;
            const int v = std::atoi(s);
            return (v > 0) ? v : 0;
        }();

        std::vector<std::size_t> node_counts(N, 0);
        exec_order.resize(N, 0);
        std::size_t max_nodes = 0;
        for (std::size_t i = 0; i < N; ++i) {
            node_counts[i] = f_work[i].get_impl()->get_coeffs().size();
            max_nodes = std::max(max_nodes, node_counts[i]);
            exec_order[i] = i;
        }

        if (env_batch_workers > 0) {
            outer_batch_tile =
                std::min<std::size_t>(N, static_cast<std::size_t>(env_batch_workers));
        } else if (max_nodes < 1024u) {
            outer_batch_tile = N;
        } else if (max_nodes < 8192u) {
            outer_batch_tile = std::min<std::size_t>(N, ThreadPool::size() + 1u);
        } else if (max_nodes < 16384u) {
            outer_batch_tile = std::min<std::size_t>(N, 4u);
        } else if (max_nodes < 32768u) {
            outer_batch_tile = std::min<std::size_t>(N, 2u);
        } else {
            outer_batch_tile = 1u;
        }
#ifdef HAVE_PARSEC
        // Under PaRSEC, outer tasks execute inner sweeps inline without nested `ForEachTask`
        // blocking; keep `outer_batch_tile = N` for small trees (`max_nodes < 1024`, e.g., HF
        // exchange pairs and loose-protocol orbitals) and use `outer_batch_tile = 1` on larger
        // trees (`max_nodes >= 1024`) so all PaRSEC threads collaborate on each orbital.
        if (!windowed_parallel_detail::can_spawn_and_wait_parsec() ||
            (env_batch_workers == 0 && max_nodes >= 1024u)) {
            outer_batch_tile = 1u;
        }
#endif

        if (outer_batch_tile >= 2u) {
            std::stable_sort(exec_order.begin(),
                             exec_order.end(),
                             [&node_counts](std::size_t a, std::size_t b) {
                                 return node_counts[a] > node_counts[b];
                             });
        }
    }

    if (outer_batch_tile >= 2u) {
        for (std::size_t s0 = 0; s0 < N; s0 += outer_batch_tile) {
            const std::size_t s1 = std::min<std::size_t>(N, s0 + outer_batch_tile);
            windowed_parallel_detail::windowed_for_each_outer_tasks_if_heavy(
                target_world,
                s1 - s0,
                [&](std::size_t k_off) {
                    WindowedBatchTaskScope batch_scope;
                    const std::size_t idx = exec_order[s0 + k_off];
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
                },
                [&](std::size_t k_off) noexcept {
                    const auto* impl = f_work[exec_order[s0 + k_off]].get_impl().get();
                    const std::size_t thresh_nodes =
                        (impl->get_tree_state() == compressed) ? 32u : 256u;
                    return impl->get_coeffs().size() >= thresh_nodes;
                });
        }
    } else if (target_world.size() > 1) {
        std::vector<const WindowedQuadraturePlan*> plan_ptrs(N, nullptr);
        std::vector<const WindowedStencilCache1D*> stencil_ptrs(N, nullptr);
        std::vector<std::uint8_t> dyn_flags(N, 0u);
        std::vector<double> eps_locs(N, 0.0);
        for (std::size_t i = 0; i < N; ++i) {
            const auto& opts_i = ops[i]->options();
            plan_ptrs[i] = resolved_plans[i].plan;
            stencil_ptrs[i] = resolved_plans[i].stencil;
            dyn_flags[i] = static_cast<std::uint8_t>(opts_i.enable_dynamic_refinement ? 1u : 0u);
            eps_locs[i] = opts_i.eps_loc;
        }
        apply_windowed_ns_3d_batch_distributed<R>(target_world,
                                                  f_work,
                                                  result,
                                                  max_depths,
                                                  plan_ptrs,
                                                  stencil_ptrs,
                                                  eff_thresh,
                                                  dyn_flags,
                                                  eff_refine_fac,
                                                  eps_locs);
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
    cfg.gamma_spiky_override = select_operator_gamma_spiky_override(k_);
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
    cfg.gamma_spiky_override = select_operator_gamma_spiky_override(k_);
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
