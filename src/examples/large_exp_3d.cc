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

/// \file large_exp_3d.cc
/// \brief Extracted from `test6` in `large_exp.cc`: 3D convolution of normalized
///        Gaussians ranging from diffuse (`b = 0.5`) to ultra-spiky (`b = 5.50e11`)
///        against a displaced Gaussian of fixed exponent `a = 16` in `[-11, +11]^3`,
///        with full octree / box-count diagnostics comparing Legacy `SeparatedConvolution`
///        ("Before") against Tail-Completed Two-Threshold `WindowedSeparatedConvolution` ("After").

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/windowed_convolution.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

using namespace madness;

namespace {

class NormalizedGaussian : public FunctionFunctorInterface<double, 3> {
private:
    double expnt_;
    double coeff_;

public:
    explicit NormalizedGaussian(double expnt)
        : expnt_(expnt)
        , coeff_(std::pow(expnt / constants::pi, 1.5)) {}

    double operator()(const coord_3d& r) const override {
        // Displaced off the origin so adaptive refinement straddles the x = 0 dyadic divide.
        const double x = r[0] - 0.1;
        return coeff_ * std::exp(-expnt_ * (x * x + r[1] * r[1] + r[2] * r[2]));
    }
};

struct TreeInspection {
    long total_boxes = 0;
    long leaf_boxes = 0;
    long internal_boxes = 0;
    int min_leaf_level = 0;
    int max_leaf_level = 0;
    std::vector<long> leaves_per_level;
    std::vector<long> internals_per_level;
};

TreeInspection examine_tree(const real_function_3d& fn) {
    World& world = fn.world();
    const auto& coeffs = fn.get_impl()->get_coeffs();

    int local_max_lvl = 0;
    int local_min_leaf = 1000;
    long local_total = 0;
    long local_leaves = 0;
    long local_internals = 0;

    for (auto it = coeffs.begin(); it != coeffs.end(); ++it) {
        const int lvl = static_cast<int>(it->first.level());
        local_max_lvl = std::max(local_max_lvl, lvl);
        ++local_total;
        if (it->second.has_children()) {
            ++local_internals;
        } else {
            ++local_leaves;
            local_min_leaf = std::min(local_min_leaf, lvl);
        }
    }

    int global_max_lvl = local_max_lvl;
    int global_min_leaf = local_min_leaf;
    if (world.size() > 1) {
        world.gop.max(global_max_lvl);
        world.gop.min(global_min_leaf);
    }
    if (global_min_leaf > global_max_lvl) {
        global_min_leaf = 0;
    }

    std::vector<long> leaves_hist(static_cast<std::size_t>(global_max_lvl + 1), 0L);
    std::vector<long> internals_hist(static_cast<std::size_t>(global_max_lvl + 1), 0L);
    for (auto it = coeffs.begin(); it != coeffs.end(); ++it) {
        const std::size_t lvl = static_cast<std::size_t>(it->first.level());
        if (it->second.has_children()) {
            ++internals_hist[lvl];
        } else {
            ++leaves_hist[lvl];
        }
    }

    if (world.size() > 1) {
        world.gop.sum(local_total);
        world.gop.sum(local_leaves);
        world.gop.sum(local_internals);
        world.gop.sum(leaves_hist.data(), leaves_hist.size());
        world.gop.sum(internals_hist.data(), internals_hist.size());
    }

    TreeInspection out;
    out.total_boxes = local_total;
    out.leaf_boxes = local_leaves;
    out.internal_boxes = local_internals;
    out.min_leaf_level = global_min_leaf;
    out.max_leaf_level = global_max_lvl;
    out.leaves_per_level = std::move(leaves_hist);
    out.internals_per_level = std::move(internals_hist);
    return out;
}

std::string format_level_histogram(const std::vector<long>& hist) {
    std::ostringstream oss;
    bool first = true;
    for (std::size_t lvl = 0; lvl < hist.size(); ++lvl) {
        if (hist[lvl] == 0) continue;
        if (!first) oss << " ";
        oss << "L" << lvl << ":" << hist[lvl];
        first = false;
    }
    return first ? std::string("L0:0") : oss.str();
}

std::string classify_regime(const WindowedQuadraturePlan& plan, int max_tree_depth) {
    if (plan.terms().empty()) return "empty";
    const auto& t = plan.terms().front();
    if (t.n_diff < 0) {
        return "root-diff";
    }
    if (t.n_spiky >= max_tree_depth) {
        return "pure-GC";
    }
    std::ostringstream oss;
    const int j_hi = std::min(max_tree_depth, t.n_diff);
    oss << "win[" << t.j0 << ".." << j_hi << "]";
    return oss.str();
}

struct ScaleDetailRecord {
    int iscale = 0;
    double b = 0.0;
    std::string regime;
    TreeInspection exact_tree;
    TreeInspection before_raw;
    TreeInspection before_trunc;
    TreeInspection after_raw;
    TreeInspection after_trunc;
};

} // namespace

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    {
    startup(world, argc, argv, true);

    const int k = 10;
    const double thresh = 1e-8;
    const double L = 11.0;
    const double lo = 1e-6;

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_truncate_mode(0);
    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_cubic_cell(-L, L);

    const double a = 16.0;
    const double sigma_a = 1.0 / std::sqrt(2.0 * a);
    real_function_3d f = real_factory_3d(world).functor(NormalizedGaussian(a));
    f.truncate();
    f.reconstruct();

    const TreeInspection f_tree = examine_tree(f);
    const Level depth = f.max_depth();
    const double userh = (2.0 * L) * std::pow(2.0, -depth);
    const double input_err = f.err(NormalizedGaussian(a));
    const double f_norm = f.norm2();

    if (world.rank() == 0) {
        std::printf("\n====================================================================================================\n");
        std::printf(" 3D Gaussian Convolution Benchmark & Tree Inspection (extracted from large_exp.cc::test6)\n");
        std::printf("====================================================================================================\n");
        std::printf(" Cubic domain   : [-%.1f, +%.1f]^3 (width = %.1f), k = %d, thresh = %.1e, truncate_mode = 0\n",
                    L, L, 2.0 * L, k, thresh);
        std::printf(" Input Gaussian : a = %.6f, sigma_a = %.6f, center = (0.1, 0, 0), ||f||_2 = %.6f\n",
                    a, sigma_a, f_norm);
        std::printf(" Input MRA tree : total_boxes = %ld, leaf_boxes = %ld, internal_boxes = %ld, levels = [%d..%d]\n",
                    f_tree.total_boxes,
                    f_tree.leaf_boxes,
                    f_tree.internal_boxes,
                    f_tree.min_leaf_level,
                    f_tree.max_leaf_level);
        std::printf(" Input leaves   : %s\n", format_level_histogram(f_tree.leaves_per_level).c_str());
        std::printf(" Finest h_user  : %.6f, input projection L2 error = %.2e\n\n", userh, input_err);

        std::printf(" isc  userh/sig     expnt_b     regime     ex_leaf |  bef_raw  bef_trn  bef_dep   bef_err   bef_ms |  aft_raw  aft_trn  aft_dep   aft_err   aft_ms | leaf_rat  speedup\n");
        std::printf(" ---  ----------  ----------  -----------  ------- | -------- -------- -------- --------- -------- | -------- -------- -------- --------- -------- | -------- --------\n");
        std::fflush(stdout);
    }

    // Warm up libxsmm / thread pool / caches once before timed sweep.
    {
        Tensor<double> w_exp(1L), w_coeff(1L);
        w_exp[0] = 1.0;
        w_coeff[0] = std::pow(1.0 / constants::pi, 1.5);
        SeparatedConvolution<double, 3> Gw0(world, w_coeff, w_exp, lo, thresh);
        Gw0.enable_windowed_ns(false);
        real_function_3d tmp0 = Gw0(f);
        SeparatedConvolution<double, 3> Gw1(world, w_coeff, w_exp, lo, thresh);
        Gw1.enable_windowed_ns(true);
        real_function_3d tmp1 = Gw1(f);
    }

    std::vector<ScaleDetailRecord> detail_records;
    double sum_before_ms = 0.0;
    double sum_after_ms = 0.0;
    long sum_before_raw_leaves = 0;
    long sum_after_raw_leaves = 0;
    long sum_before_trunc_leaves = 0;
    long sum_after_trunc_leaves = 0;
    double max_before_err = 0.0;
    double max_after_err = 0.0;

    for (int iscale = 0; iscale <= 40; ++iscale) {
        const double sigma = std::pow(2.0, -0.5 * iscale);
        const double b = 1.0 / (2.0 * sigma * sigma);
        const double ab_eff = (a * b) / (a + b);
        NormalizedGaussian exact_fn(ab_eff);

        real_function_3d exact_proj = real_factory_3d(world).functor(exact_fn);
        exact_proj.truncate();
        exact_proj.reconstruct();
        const TreeInspection exact_tree = examine_tree(exact_proj);

        Tensor<double> expnt(1L), coeff(1L);
        expnt[0] = b;
        coeff[0] = std::pow(b / constants::pi, 1.5);

        // 1. BEFORE: Legacy SeparatedConvolution
        SeparatedConvolution<double, 3> G_before(world, coeff, expnt, lo, thresh);
        G_before.enable_windowed_ns(false);
        world.gop.fence();
        const double t0_bef = wall_time();
        real_function_3d ff_before = G_before(f);
        world.gop.fence();
        const double bef_ms = 1000.0 * (wall_time() - t0_bef);

        const TreeInspection bef_raw = examine_tree(ff_before);
        const double bef_err = ff_before.err(exact_fn);
        real_function_3d ff_before_trn = copy(ff_before);
        ff_before_trn.truncate();
        ff_before_trn.reconstruct();
        const TreeInspection bef_trn = examine_tree(ff_before_trn);

        // 2. AFTER: WindowedSeparatedConvolution (via SeparatedConvolution opt-in bridge)
        SeparatedConvolution<double, 3> G_after(world, coeff, expnt, lo, thresh);
        G_after.enable_windowed_ns(true);
        const auto wop_ptr = G_after.get_windowed_ns_op();
        const std::string regime =
            wop_ptr ? classify_regime(wop_ptr->plan(), static_cast<int>(depth)) : std::string("fallback");

        world.gop.fence();
        const double t0_aft = wall_time();
        real_function_3d ff_after = G_after(f);
        world.gop.fence();
        const double aft_ms = 1000.0 * (wall_time() - t0_aft);

        const TreeInspection aft_raw = examine_tree(ff_after);
        const double aft_err = ff_after.err(exact_fn);
        real_function_3d ff_after_trn = copy(ff_after);
        ff_after_trn.truncate();
        ff_after_trn.reconstruct();
        const TreeInspection aft_trn = examine_tree(ff_after_trn);

        sum_before_ms += bef_ms;
        sum_after_ms += aft_ms;
        sum_before_raw_leaves += bef_raw.leaf_boxes;
        sum_after_raw_leaves += aft_raw.leaf_boxes;
        sum_before_trunc_leaves += bef_trn.leaf_boxes;
        sum_after_trunc_leaves += aft_trn.leaf_boxes;
        max_before_err = std::max(max_before_err, bef_err);
        max_after_err = std::max(max_after_err, aft_err);

        const double leaf_ratio =
            (aft_raw.leaf_boxes > 0)
                ? (static_cast<double>(bef_raw.leaf_boxes) / static_cast<double>(aft_raw.leaf_boxes))
                : 1.0;
        const double speedup = (aft_ms > 0.0) ? (bef_ms / aft_ms) : 1.0;

        if (world.rank() == 0) {
            std::printf(" %3d  %10.3e  %10.3e  %-11s  %7ld | %8ld %8ld %8d %9.2e %8.2f | %8ld %8ld %8d %9.2e %8.2f | %7.2fx %7.2fx\n",
                        iscale,
                        userh / sigma,
                        b,
                        regime.c_str(),
                        exact_tree.leaf_boxes,
                        bef_raw.leaf_boxes,
                        bef_trn.leaf_boxes,
                        bef_raw.max_leaf_level,
                        bef_err,
                        bef_ms,
                        aft_raw.leaf_boxes,
                        aft_trn.leaf_boxes,
                        aft_raw.max_leaf_level,
                        aft_err,
                        aft_ms,
                        leaf_ratio,
                        speedup);
            std::fflush(stdout);
        }

        if (iscale == 0 || iscale == 4 || iscale == 8 || iscale == 10 ||
            iscale == 12 || iscale == 16 || iscale == 24 || iscale == 40) {
            detail_records.push_back(ScaleDetailRecord{
                iscale,
                b,
                regime,
                exact_tree,
                bef_raw,
                bef_trn,
                aft_raw,
                aft_trn});
        }
    }

    if (world.rank() == 0) {
        std::printf(" ---  ----------  ----------  -----------  ------- | -------- -------- -------- --------- -------- | -------- -------- -------- --------- -------- | -------- --------\n");
        std::printf(" SUM / WORST OVER 41 SCALES (iscale = 0..40)       | %8ld %8ld        - %9.2e %8.2f | %8ld %8ld        - %9.2e %8.2f | %7.2fx %7.2fx\n",
                    sum_before_raw_leaves,
                    sum_before_trunc_leaves,
                    max_before_err,
                    sum_before_ms,
                    sum_after_raw_leaves,
                    sum_after_trunc_leaves,
                    max_after_err,
                    sum_after_ms,
                    static_cast<double>(sum_before_raw_leaves) /
                        std::max(1.0, static_cast<double>(sum_after_raw_leaves)),
                    sum_before_ms / std::max(1e-9, sum_after_ms));

        std::printf("\n====================================================================================================\n");
        std::printf(" Detailed Octree Per-Level Leaf-Box Distributions at Representative Scales\n");
        std::printf("====================================================================================================\n");
        std::printf(" Input f (a=16)  : total=%5ld leaves=%5ld [%d..%d] | %s\n",
                    f_tree.total_boxes,
                    f_tree.leaf_boxes,
                    f_tree.min_leaf_level,
                    f_tree.max_leaf_level,
                    format_level_histogram(f_tree.leaves_per_level).c_str());
        for (const auto& rec : detail_records) {
            std::printf("\n [iscale = %2d, b = %.3e, regime = %s]\n",
                        rec.iscale, rec.b, rec.regime.c_str());
            std::printf("   Exact target  : total=%5ld leaves=%5ld [%d..%d] | %s\n",
                        rec.exact_tree.total_boxes,
                        rec.exact_tree.leaf_boxes,
                        rec.exact_tree.min_leaf_level,
                        rec.exact_tree.max_leaf_level,
                        format_level_histogram(rec.exact_tree.leaves_per_level).c_str());
            std::printf("   Before (raw)  : total=%5ld leaves=%5ld [%d..%d] | %s\n",
                        rec.before_raw.total_boxes,
                        rec.before_raw.leaf_boxes,
                        rec.before_raw.min_leaf_level,
                        rec.before_raw.max_leaf_level,
                        format_level_histogram(rec.before_raw.leaves_per_level).c_str());
            std::printf("   Before (trunc): total=%5ld leaves=%5ld [%d..%d] | %s\n",
                        rec.before_trunc.total_boxes,
                        rec.before_trunc.leaf_boxes,
                        rec.before_trunc.min_leaf_level,
                        rec.before_trunc.max_leaf_level,
                        format_level_histogram(rec.before_trunc.leaves_per_level).c_str());
            std::printf("   After  (raw)  : total=%5ld leaves=%5ld [%d..%d] | %s\n",
                        rec.after_raw.total_boxes,
                        rec.after_raw.leaf_boxes,
                        rec.after_raw.min_leaf_level,
                        rec.after_raw.max_leaf_level,
                        format_level_histogram(rec.after_raw.leaves_per_level).c_str());
            std::printf("   After  (trunc): total=%5ld leaves=%5ld [%d..%d] | %s\n",
                        rec.after_trunc.total_boxes,
                        rec.after_trunc.leaf_boxes,
                        rec.after_trunc.min_leaf_level,
                        rec.after_trunc.max_leaf_level,
                        format_level_histogram(rec.after_trunc.leaves_per_level).c_str());
        }
        std::printf("====================================================================================================\n");
    }
    world.gop.fence();
    }

    finalize();
    return 0;
}
