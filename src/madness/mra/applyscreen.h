#ifndef MADNESS_MRA_APPLYSCREEN_H__INCLUDED
#define MADNESS_MRA_APPLYSCREEN_H__INCLUDED

/// \file applyscreen.h
/// \brief Destination screening of the apply of an integral operator.
///
/// The result of an apply is often multiplied with another function afterwards (the exchange operator:
/// K f = sum_i ket_i * P(bra_i f)). A contribution to a box in which that multiplier is negligible is then
/// never used, and need not be computed (Yanai, Fann, Gan, Harrison, Beylkin, JCP 121, 6680 (2004), eq. 17).
/// A screen answers, per destination box, a bound on the multiplier's largest magnitude there (supbound.h).
/// FunctionImpl::do_apply skips a block whose contribution times that bound, capped at one, is negligible,
/// computes a kept block with its per-term tolerance loosened by the same factor, and accumulates it only
/// if its norm times the factor is significant; a skipped block still counts as used by the shell loop.
/// A skipped block leaves at most the block's share of the tolerance in the product function, since
/// ||s du||_2 <= sup|s| ||du||_2. The tasks of one apply carry the screen (ApplyScreenHandle; apply_only,
/// and the vector apply(world, op, f, screens)); without one, do_apply is unchanged.

#include <madness/mra/key.h>
#include <madness/mra/supbound.h>
#include <madness/world/madness_exception.h>
#include <madness/world/world.h>
#include <cstdio>
#include <algorithm>
#include <array>
#include <atomic>
#include <functional>
#include <memory>
#include <vector>

namespace madness {

    /// per destination box, a bound on the size of what the result of the apply will be multiplied with
    template <std::size_t NDIM>
    using apply_screen_t = std::function<double(const Key<NDIM>&)>;

    /// The screen as the tasks of an apply carry it: shared ownership, so the screen lives as long as
    /// any task of the apply does, and nothing is stored on the function. Empty means no screen. The
    /// tasks stay on the rank that creates them (apply_only refuses apply_randomize), so the handle
    /// does not travel.
    template <std::size_t NDIM>
    class ApplyScreenHandle {
        std::shared_ptr<const apply_screen_t<NDIM>> screen_;
    public:
        ApplyScreenHandle() = default;
        explicit ApplyScreenHandle(const apply_screen_t<NDIM>& screen)
            : screen_(screen ? std::make_shared<const apply_screen_t<NDIM>>(screen) : nullptr) {}

        explicit operator bool() const { return bool(screen_); }
        double operator()(const Key<NDIM>& key) const { return (*screen_)(key); }

        template <typename Archive>
        void serialize(Archive&) {
            MADNESS_EXCEPTION("no serialization of ApplyScreenHandle: the apply screen stays on its rank", 1);
        }
    };

    /// Counters of the apply screen, process-wide: the blocks that passed the source test while a screen
    /// was installed and the ones it skipped, and the sup bounds built for screens with the time it took.
    /// Reset and reported by whoever installs the screens (the exchange operator, summed over the ranks);
    /// a caller that runs the exchange tasks in its own queue reports after the queue ran.
    struct ApplyScreenStats {
        static inline std::atomic<long> blocks{0}, skipped{0}, sup_maps{0}, sup_build_ns{0};
        static void reset() { blocks = 0; skipped = 0; sup_maps = 0; sup_build_ns = 0; }

        /// one line on rank 0 with the counters summed over the ranks. Collective.
        static void report(World& world) {
            std::array<long, 4> n{long(blocks), long(skipped), long(sup_maps), long(sup_build_ns)};
            world.gop.sum(n.data(), n.size());
            if (world.rank() != 0) return;
            printf(" apply screen: blocks %ld, skipped %ld (%.1f%%)", n[0], n[1],
                   n[0] > 0 ? 100.0 * double(n[1]) / double(n[0]) : 0.0);
            if (n[2] > 0)
                printf("; sup bounds of %ld multipliers built in %.2fs (%.1f ms each)", n[2], 1e-9 * double(n[3]),
                       1e-6 * double(n[3]) / double(n[2]));
            printf("\n");
        }
    };

    /// The screen built from the multipliers' sup bounds (SupBoundMap, FunctionImpl::sup_bound_map): the
    /// largest of their bounds in the box. A box none of the maps holds (not local to this rank, or without
    /// coefficients) answers unbounded, which screens nothing.
    template <std::size_t NDIM>
    class SupBoundScreen {
        std::vector<std::shared_ptr<const SupBoundMap<NDIM>>> sup_;
    public:
        static constexpr double unbounded() { return SupBoundMap<NDIM>::unbounded(); }

        explicit SupBoundScreen(std::vector<std::shared_ptr<const SupBoundMap<NDIM>>> sup) : sup_(std::move(sup)) {
            MADNESS_CHECK_THROW(not sup_.empty(), "SupBoundScreen: at least one multiplier");
            for (const auto& s : sup_) MADNESS_CHECK_THROW(bool(s), "SupBoundScreen: a missing sup bound");
        }

        double operator()(const Key<NDIM>& key) const {
            double s = 0.0;
            for (const auto& m : sup_) s = std::max(s, (*m)(key));
            return s;
        }
    };

}   // namespace madness

#endif   // MADNESS_MRA_APPLYSCREEN_H__INCLUDED
