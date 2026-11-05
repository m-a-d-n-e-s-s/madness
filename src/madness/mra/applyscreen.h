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
/// ||s du||_2 <= sup|s| ||du||_2. The screen rides on the result of one apply (apply_only, and the vector
/// apply(world, op, f, screens)); without one, do_apply is unchanged.
///
/// MAD_APPLY_SCREEN=0 in the environment switches every destination screen in the process off: the applies
/// then run as without a screen, whatever a caller installs. It cannot switch a screen on.

#include <madness/mra/key.h>
#include <madness/mra/supbound.h>
#include <madness/world/madness_exception.h>
#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <memory>
#include <vector>

namespace madness {

    /// per destination box, a bound on the size of what the result of the apply will be multiplied with
    template <std::size_t NDIM>
    using apply_screen_t = std::function<double(const Key<NDIM>&)>;

    /// Is the destination screen switched off for this process (MAD_APPLY_SCREEN=0)? Read once.
    inline bool apply_screen_disabled() {
        static const bool off = [] {
            const char* v = std::getenv("MAD_APPLY_SCREEN");
            return v != nullptr and (std::strcmp(v, "0") == 0 or std::strcmp(v, "off") == 0 or std::strcmp(v, "false") == 0);
        }();
        return off;
    }

    /// Counters of the apply screen, process-wide: the blocks that passed the source test while a screen
    /// was installed and the ones it skipped, and the sup bounds built for screens with the time it took.
    /// Reset and reported by whoever installs the screens (the exchange operator, summed over the ranks).
    struct ApplyScreenStats {
        static inline std::atomic<long> blocks{0}, skipped{0}, sup_maps{0}, sup_build_ns{0};
        static void reset() { blocks = 0; skipped = 0; sup_maps = 0; sup_build_ns = 0; }
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
