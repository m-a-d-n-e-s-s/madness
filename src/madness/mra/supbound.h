#ifndef MADNESS_MRA_SUPBOUND_H__INCLUDED
#define MADNESS_MRA_SUPBOUND_H__INCLUDED

/// \file supbound.h
/// \brief Bounds on the largest magnitude of a function inside a box, from its scaling coefficients.
///
/// Inside the box (n, l) a function with scaling coefficients c is, in user coordinates,
///
///     f(x) = 2^{NDIM n/2} V^{-1/2} sum_i c_i prod_d phi_{i_d}(2^n x_d - l_d)
///
/// (x in simulation coordinates, V the cell volume; FunctionCommonData::coeffs2values), and the normalized
/// Legendre scaling functions phi_i(t) = sqrt(2i+1) P_i(2t-1) peak at sqrt(2i+1) on the box faces. Hence
///
///     sup_box |f| <= 2^{NDIM n/2} V^{-1/2} sum_i |c_i| prod_d sqrt(2 i_d + 1),
///
/// rigorously, at k^NDIM multiply-adds per box. The bound lies above the box's rms magnitude
/// ||c|| / sqrt(V_box) by a factor between 1 (a constant) and k^NDIM; for a function that is smooth across the
/// box it is close to 1, and it is the larger the more f varies inside the box.
///
/// SupBoundMap holds that bound for the boxes of a tree: a leaf's from its coefficients, an interior node's as
/// the max over its children (the sup over a union is the max of the sups), and a box below a leaf takes the
/// leaf's (the sup over a sub-box is at most the sup over the box). An interior node's own coefficients (a
/// redundant tree) are a projection of f and bound nothing, so only the leaves' coefficients are read.
/// FunctionImpl::sup_bound_map() builds the map for a function.

#include <madness/mra/key.h>
#include <madness/mra/funcdefaults.h>
#include <madness/tensor/tensor.h>
#include <madness/world/worldhash.h>

#include <cmath>
#include <unordered_map>

namespace madness {

    /// Bounds on sup|f| per box of a function's tree (supbound.h)
    template <std::size_t NDIM>
    class SupBoundMap {
    public:
        typedef std::unordered_map<Key<NDIM>, double, Hash<Key<NDIM>>> mapT;

        /// the answer for a box the map holds no bound for
        static constexpr double unbounded() { return 1e300; }

        /// the weights of the leaf bound, prod_d sqrt(2 i_d + 1), for k coefficients per dimension
        static Tensor<double> leaf_weights(const int k) {
            Tensor<double> w1(k);
            for (int i = 0; i < k; ++i) w1(i) = std::sqrt(2.0 * i + 1.0);
            Tensor<double> w = w1;
            for (std::size_t d = 1; d < NDIM; ++d) w = outer(w, w1);
            return w;
        }

        /// A bound on sup|f| over the box `key`, from f's scaling coefficients `c` there and the weights of
        /// leaf_weights(k). No coefficients means f vanishes in the box.
        template <typename T>
        static double leaf_bound(const Key<NDIM>& key, const Tensor<T>& c, const Tensor<double>& weights) {
            if (c.size() == 0) return 0.0;
            MADNESS_CHECK_THROW(c.size() == weights.size(), "SupBoundMap: coefficients and weights differ in size");
            const Tensor<T> cc = c.iscontiguous() ? c : copy(c);
            const T* p = cc.ptr();
            const double* w = weights.ptr();
            double s = 0.0;
            for (long i = 0; i < cc.size(); ++i) s += std::abs(p[i]) * w[i];
            return s * std::pow(2.0, 0.5 * NDIM * key.level()) / std::sqrt(FunctionDefaults<NDIM>::get_cell_volume());
        }

        SupBoundMap() = default;
        explicit SupBoundMap(mapT bound) : bound_(std::move(bound)) {}

        /// The bound for the box `key`: its own, or the leaf's it lies below; unbounded() for a box the map
        /// holds neither for. `below_leaf`, if given, is set when the answer is a leaf's above `key`.
        double operator()(Key<NDIM> key, bool* below_leaf = nullptr) const {
            for (bool below = false;; below = true) {
                const auto it = bound_.find(key);
                if (it != bound_.end()) {
                    if (below_leaf) *below_leaf = below;
                    return it->second;
                }
                if (key.level() == 0) return unbounded();
                key = key.parent();
            }
        }

        /// the number of boxes held
        std::size_t size() const { return bound_.size(); }

    private:
        mapT bound_;
    };

}   // namespace madness

#endif   // MADNESS_MRA_SUPBOUND_H__INCLUDED
