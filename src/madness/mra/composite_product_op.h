/*
  This file is part of MADNESS.

  Unified adaptive NS-form product/inner operator.

  Generalises three pre-existing operators in funcimpl.h:
    * multiply_op<LDIM>           : h(1,2) = f(1,2) * g(p)
    * Vphi_op_NS<opT,LDIM>        : result = sum_t product_t       (sum of single-multiplications)
    * Vphi_ij_u_op_NS<opT,LDIM>   : u * i(1) * j(2)  and  <eri | u*i*j>

  A composite_term carries a ket (or a hartree pair p1*p2), a chain of LDIM factors per
  particle, and an optional NDIM "extra" factor (typically an on-demand operator such as
  an ERI).  The result is sum over terms of the products.  An inner-product mode contracts
  each term against an NDIM "bra" function and accumulates a scalar at key0.

  Error control follows Vphi_ij_u_op_NS::compute_error_from_inaccurate_refinement:
    (s,d)-norm cross-terms over {ket, factors_1, factors_2, extra} summed across the
    multiplication chain, plus tnorm/oversampling error from each pointwise_multiplier.
*/
#ifndef MADNESS_MRA_COMPOSITE_PRODUCT_OP_H__
#define MADNESS_MRA_COMPOSITE_PRODUCT_OP_H__

#include <madness/mra/funcimpl.h>

namespace madness {

/// One term of a composite product: ket(1,2) * factor1(1) * factor2(2) * extra(1,2)
///
/// The ket is either an explicit NDIM Function (`ket`) or a low-rank sum of hartree
/// products  ket(1,2) = Σ_k p1[k](1) ⊗ p2[k](2).  Provide `ket` *or* (p1, p2); not both.
template <typename T, std::size_t NDIM>
struct composite_term {
    static constexpr std::size_t LDIM = NDIM / 2;
    static_assert(NDIM == 2 * LDIM, "composite_term requires NDIM = 2*LDIM");

    using FN = Function<T, NDIM>;
    using FL = Function<T, LDIM>;

    /// NDIM ket factor. If uninitialized, the term ket is Σ_k p1[k] ⊗ p2[k].
    FN ket;
    /// parallel vectors forming the ket as a sum of hartree products (p1.size() == p2.size())
    std::vector<FL> p1, p2;
    /// single LDIM factor per particle (uninitialized = skip)
    FL factor1, factor2;
    /// optional single NDIM extra factor (e.g. an on-demand ERI).  Multiplied into the
    /// chain in function mode; ignored in inner mode (use the global `bra` instead).
    FN extra;
};

/// Adaptive sum-of-products operator over a vector of composite_terms.
template <typename T, std::size_t NDIM, typename opT, std::size_t LDIM>
struct composite_product_op_NS {
    static_assert(NDIM == 2 * LDIM, "composite_product_op_NS requires NDIM = 2*LDIM");

    enum class Mode { function, inner };

    using this_type = composite_product_op_NS<T, NDIM, opT, LDIM>;
    using implT     = FunctionImpl<T, NDIM>;
    using implL     = FunctionImpl<T, LDIM>;
    using keyT      = Key<NDIM>;
    using keyL      = Key<LDIM>;
    using nodeT     = typename implT::nodeT;
    using coeffT    = typename implT::coeffT;
    using tensorT   = Tensor<T>;
    using ctT       = CoeffTracker<T, NDIM>;
    using ctL       = CoeffTracker<T, LDIM>;

    /// internal serialisable view of a term (CoeffTracker-based)
    struct TermImpl {
        ctT iaket;                       ///< NDIM ket (optional, exclusive with iap*)
        std::vector<ctL> iap1, iap2;     ///< sum-of-hartree-products ket if iaket is null
        ctL ifactor1, ifactor2;          ///< single per-particle LDIM factor (each optional)
        const implT* extra = nullptr;    ///< NDIM extra factor (function mode only)

        bool have_ket()     const { return iaket.get_impl(); }
        bool have_factor1() const { return ifactor1.get_impl(); }
        bool have_factor2() const { return ifactor2.get_impl(); }
        bool have_extra()   const { return extra; }

        template <typename Archive> void serialize(Archive& ar) {
            ar & iaket & iap1 & iap2 & ifactor1 & ifactor2 & extra;
        }
    };

    bool randomize() const { return true; }

    implT*               result;            ///< output tree (function mode) or scalar sink (inner)
    opT                  leaf_op;
    std::vector<TermImpl> terms;
    Mode                 mode = Mode::function;
    const implT*         bra  = nullptr;    ///< inner-mode bra (e.g., on-demand ERI)
    double               target_precision = 0.0;
    int                  oversampling = 1;

    composite_product_op_NS() : result(nullptr) {}

    composite_product_op_NS(implT* result, const opT& leaf_op,
                            const std::vector<TermImpl>& terms,
                            Mode mode,
                            const implT* bra,
                            double target_precision,
                            int oversampling)
        : result(result), leaf_op(leaf_op), terms(terms),
          mode(mode), bra(bra),
          target_precision(target_precision), oversampling(oversampling) {
        for (const auto& t : terms)
            if (t.extra) MADNESS_ASSERT(t.extra->is_on_demand() || not t.extra->is_on_demand());
        if (mode == Mode::inner) MADNESS_ASSERT(bra && bra->is_on_demand());
    }

    // ---------------------------------------------------------------------
    // error model: snorm/dnorm cross-terms over the factor chain of one term
    // ---------------------------------------------------------------------

    /// retrieve scaling- and detail-norm at the given key from a CoeffTracker
    template <typename CT, typename KEY>
    static std::pair<double, double> sd_norm(const CT& ct, const KEY& key) {
        if (!ct.get_impl()) return {1.0, 0.0};
        return {ct.coeff(key).normf(), ct.dnorm(key)};
    }

    /// scaling- and detail-norm for an on-demand NDIM extra factor.  Returns
    /// (snorm, dnorm) ~ (||P_n eri||, ||Q_n eri||), with dnorm computed from the
    /// NS s0 / non-s0 split of the eri coeff tensor.
    static std::pair<double, double> sd_norm_extra(const implT* extra, const keyT& key) {
        if (!extra) return {1.0, 0.0};
        tensorT c = eri_coeffs(extra, key);
        const double sn = c(extra->get_cdata().s0).normf();
        tensorT d = copy(c);
        d(extra->get_cdata().s0) = 0.0;
        return {sn, d.normf()};
    }

    /// Σ over all non-trivial cross-products of (s_i + d_i) − ∏ s_i.
    /// Each factor contributes (s,d); error = total − product_of_s.
    static double cross_term_error(const std::vector<std::pair<double, double>>& sd) {
        double total = 1.0, sprod = 1.0;
        for (const auto& [s, d] : sd) { total *= (s + d); sprod *= s; }
        return total - sprod;
    }

    /// (s,d) of a sum-of-hartree-products ket: snorm = Σ_k s1_k s2_k,
    /// dnorm = Σ_k (s1_k d2_k + d1_k s2_k + d1_k d2_k).  Upper-bound on the true
    /// (s,d) of the rank-r sum.
    std::pair<double, double> sd_norm_pair_sum(const std::vector<ctL>& iap1,
                                                const std::vector<ctL>& iap2,
                                                const keyL& k1, const keyL& k2) const {
        MADNESS_ASSERT(iap1.size() == iap2.size());
        double s = 0.0, d = 0.0;
        for (std::size_t k = 0; k < iap1.size(); ++k) {
            const auto [s1, d1] = sd_norm(iap1[k], k1);
            const auto [s2, d2] = sd_norm(iap2[k], k2);
            s += s1 * s2;
            d += s1 * d2 + d1 * s2 + d1 * d2;
        }
        return {s, d};
    }

    /// full inaccurate-refinement error for one term
    double term_refinement_error(const TermImpl& t, const keyT& key,
                                 const keyL& k1, const keyL& k2) const {
        std::vector<std::pair<double, double>> sd;
        if (t.have_ket()) sd.push_back(sd_norm(t.iaket, key));
        else              sd.push_back(sd_norm_pair_sum(t.iap1, t.iap2, k1, k2));
        if (t.have_factor1()) sd.push_back(sd_norm(t.ifactor1, k1));
        if (t.have_factor2()) sd.push_back(sd_norm(t.ifactor2, k2));
        if (t.have_extra())   sd.push_back(sd_norm_extra(t.extra, key));
        return cross_term_error(sd);
    }

    // ---------------------------------------------------------------------
    // multiplication chain for one term: current = ket → ket·f1 → … → ket·…·extra
    // ---------------------------------------------------------------------

    /// fetch coeffs of an on-demand NDIM function at `key`
    static tensorT eri_coeffs(const implT* eri, const keyT& key) {
        MADNESS_ASSERT(eri);
        if (eri->get_functor()->provides_coeff()) {
            return eri->get_functor()->coeff(key).full_tensor();
        } else {
            tensorT val(eri->get_cdata().vk);
            eri->fcube(key, *(eri->get_functor()), eri->get_cdata().quad_x, val);
            return eri->values2coeffs(key, val);
        }
    }

    /// build the NDIM ket coefficients at `key` for a term: either the explicit NDIM ket
    /// or the sum Σ_k outer(p1[k](k1), p2[k](k2)).
    coeffT make_term_ket(const TermImpl& t, const keyT& key,
                         const keyL& k1, const keyL& k2) const {
        if (t.have_ket()) return t.iaket.coeff(key);
        MADNESS_ASSERT(!t.iap1.empty() && t.iap1.size() == t.iap2.size());
        coeffT acc = outer(t.iap1[0].coeff(k1), t.iap2[0].coeff(k2),
                           result->get_tensor_args());
        for (std::size_t k = 1; k < t.iap1.size(); ++k) {
            acc += outer(t.iap1[k].coeff(k1), t.iap2[k].coeff(k2),
                         result->get_tensor_args());
        }
        acc.reduce_rank(result->get_tensor_args().thresh);
        return acc;
    }

    /// multiply one term: ket * factor1(1) * factor2(2) * extra(1,2)
    /// Returns (coeffs at key, accumulated tnorm/oversampling error).
    std::pair<coeffT, double> make_term_coeffs(const TermImpl& t, const keyT& key) const {
        keyL k1, k2; key.break_apart(k1, k2);
        coeffT current = make_term_ket(t, key, k1, k2);
        double err = 0.0;

        if (t.have_factor1()) {
            typename implT::template pointwise_multiplier<LDIM> pm(key, current);
            pm.oversampling = oversampling;
            current = pm(key, t.ifactor1.coeff(k1).full_tensor(), 1);
            err += pm.error;
        }
        if (t.have_factor2()) {
            typename implT::template pointwise_multiplier<LDIM> pm(key, current);
            pm.oversampling = oversampling;
            current = pm(key, t.ifactor2.coeff(k2).full_tensor(), 2);
            err += pm.error;
        }
        if (t.have_extra()) {
            tensorT ceri = eri_coeffs(t.extra, key);
            typename implT::template pointwise_multiplier<LDIM> pm(key, current);
            pm.oversampling = oversampling;
            tensorT v = pm(key, ceri);
            current = coeffT(v, result->get_tensor_args());
            err += pm.error;
        }
        return {current, err};
    }

    /// Σ over terms.  Function mode: returns the sum of NDIM coeffs.
    ///                Inner mode:    returns a degenerate coeff (the scalar lives elsewhere);
    ///                               here we return the sum of contractions packed as the
    ///                               (0,…) element of a vk-shaped tensor.
    std::pair<coeffT, double> make_sum_coeffs(const keyT& key) const {
        keyL k1, k2; key.break_apart(k1, k2);

        double error = 0.0;
        if (mode == Mode::function) {
            coeffT cresult(result->get_cdata().vk, result->get_tensor_args());
            for (const auto& t : terms) {
                error += term_refinement_error(t, key, k1, k2);
                auto [c, e] = make_term_coeffs(t, key);
                cresult += c;
                error += e;
            }
            if (terms.size() > 1) cresult.reduce_rank(result->get_tensor_args().thresh);
            return {cresult, error};
        } else {  // inner mode
            const tensorT cbra = eri_coeffs(bra, key);
            T sum = T(0);
            for (const auto& t : terms) {
                error += term_refinement_error(t, key, k1, k2);
                auto [c, e] = make_term_coeffs(t, key);
                error += e;
                sum += cbra.trace(c.full_tensor_copy());
            }
            // wrap scalar in a vk tensor at key0 convention
            tensorT scalar(result->get_cdata().vk);
            scalar.flat()(0L) = sum;
            return {coeffT(scalar, TensorArgs(TT_FULL, -1.0)), error};
        }
    }

    // ---------------------------------------------------------------------
    // leaf action: insert into result tree (function), or accumulate scalar at key0 (inner)
    // ---------------------------------------------------------------------

    void accumulate_into_result(const keyT& key, const coeffT& coeff) const {
        result->get_coeffs().task(key, &nodeT::accumulate, coeff,
                                  result->get_coeffs(), key, result->get_tensor_args());
    }

    void leaf_action(const keyT& key, const coeffT& coeff) const {
        if (mode == Mode::function) {
            accumulate_into_result(key, coeff);
        } else {
            accumulate_into_result(result->get_cdata().key0, coeff);
        }
    }

    // ---------------------------------------------------------------------
    // tree traversal — mirrors Vphi_ij_u_op_NS::operator()
    // ---------------------------------------------------------------------

    std::pair<bool, coeffT> operator()(const keyT& key) const {
        MADNESS_ASSERT(result->get_coeffs().is_local(key));

        if (leaf_op.do_pre_screening()) {
            if (leaf_op.pre_screening(key)) {
                auto [c, err] = make_sum_coeffs(key);
                leaf_action(key, c);
                return {true, coeffT()};
            }
            return continue_recursion(std::vector<bool>(1 << NDIM, false), tensorT(), key);
        }

        const size_t il = result->get_initial_level();
        const int min_lvl = FunctionDefaults<NDIM>::get_refine() ? int(il) + 1 : int(il);
        if (key.level() < min_lvl)
            return continue_recursion(std::vector<bool>(1 << NDIM, false), tensorT(), key);

        if (key.level() < result->get_special_level() &&
            leaf_op.special_refinement_needed(key))
            return continue_recursion(std::vector<bool>(1 << NDIM, false), tensorT(), key);

        auto [c, err] = make_sum_coeffs(key);

        const double tol = result->truncate_tol(target_precision, key);
        if (leaf_op.post_screening(key, c) || err < tol) {
            leaf_action(key, c);
            return {true, coeffT()};
        }
        return continue_recursion(std::vector<bool>(1 << NDIM, false), tensorT(), key);
    }

    std::pair<bool, coeffT> continue_recursion(const std::vector<bool> child_is_leaf,
                                                const tensorT& coeffs,
                                                const keyT& key) const {
        std::size_t i = 0;
        for (KeyChildIterator<NDIM> kit(key); kit; ++kit, ++i) {
            const keyT child = kit.key();
            if (child_is_leaf[i]) {
                insert_op<T, NDIM> iop(result);
                iop(child, coeffT(copy(coeffs(result->child_patch(child))),
                                  result->get_tensor_args()), true);
            } else {
                this_type child_op = this->make_child(child);
                noop<T, NDIM> no;
                const ProcessID p = result->get_coeffs().owner(child);
                void (implT::*ft)(const this_type&, const noop<T, NDIM>&, const keyT&) const
                    = &implT::template forward_traverse<this_type, noop<T, NDIM>>;
                result->task(p, ft, child_op, no, child);
            }
        }
        return {true, coeffT()};
    }

    this_type make_child(const keyT& child) const {
        keyL k1, k2; child.break_apart(k1, k2);
        std::vector<TermImpl> child_terms;
        child_terms.reserve(terms.size());
        for (const auto& t : terms) {
            TermImpl ct;
            ct.iaket = t.have_ket() ? t.iaket.make_child(child) : ctT();
            if (!t.have_ket()) {
                ct.iap1.reserve(t.iap1.size());
                ct.iap2.reserve(t.iap2.size());
                for (const auto& p : t.iap1) ct.iap1.push_back(p.make_child(k1));
                for (const auto& p : t.iap2) ct.iap2.push_back(p.make_child(k2));
            }
            ct.ifactor1 = t.have_factor1() ? t.ifactor1.make_child(k1) : ctL();
            ct.ifactor2 = t.have_factor2() ? t.ifactor2.make_child(k2) : ctL();
            ct.extra = t.extra;
            child_terms.push_back(std::move(ct));
        }
        return this_type(result, leaf_op, child_terms, mode, bra,
                         target_precision, oversampling);
    }

    Future<this_type> activate() const {
        // Activate every CoeffTracker in every Term.  Existing ops bind a fixed
        // number of `Future<ct*>` to the forward_ctor task; we instead synchronously
        // materialise the activated trackers via Future::get() since the per-term
        // tracker count is variable.  CoeffTracker::activate() returns a ready future
        // for null / on-demand impls, so the .get() is free in the common case.
        std::vector<TermImpl> activated;
        activated.reserve(terms.size());
        for (const auto& t : terms) {
            TermImpl at;
            if (t.have_ket()) at.iaket = t.iaket.activate().get();
            else {
                at.iap1.reserve(t.iap1.size());
                at.iap2.reserve(t.iap2.size());
                for (const auto& p : t.iap1) at.iap1.push_back(p.activate().get());
                for (const auto& p : t.iap2) at.iap2.push_back(p.activate().get());
            }
            if (t.have_factor1()) at.ifactor1 = t.ifactor1.activate().get();
            if (t.have_factor2()) at.ifactor2 = t.ifactor2.activate().get();
            at.extra = t.extra;
            activated.push_back(std::move(at));
        }
        return result->world.taskq.add(
            detail::wrap_mem_fn(*const_cast<this_type*>(this), &this_type::forward_ctor),
            result, leaf_op, activated, mode, bra, target_precision, oversampling);
    }

    this_type forward_ctor(implT* result1, const opT& leaf_op1,
                           const std::vector<TermImpl>& terms1, Mode mode1,
                           const implT* bra1, double target_precision1, int oversampling1) {
        return this_type(result1, leaf_op1, terms1, mode1, bra1,
                         target_precision1, oversampling1);
    }

    template <typename Archive> void serialize(Archive& ar) {
        ar & result & leaf_op & terms & mode & bra & target_precision & oversampling;
    }
};

// =====================================================================================
//  ORCHESTRATORS — FunctionImpl-level free helpers used by mra.h wrappers.
// =====================================================================================

namespace detail {

/// Bring every input of a composite_term into redundant state.
template <typename T, std::size_t NDIM>
void make_term_inputs_redundant(composite_term<T, NDIM>& term) {
    if (term.ket.is_initialized())     term.ket.change_tree_state(redundant, false);
    for (auto& p : term.p1)            p.change_tree_state(redundant, false);
    for (auto& p : term.p2)            p.change_tree_state(redundant, false);
    if (term.factor1.is_initialized()) term.factor1.change_tree_state(redundant, false);
    if (term.factor2.is_initialized()) term.factor2.change_tree_state(redundant, false);
    // extra is on-demand → no state change
}

/// Restore every input of a composite_term to reconstructed state.
template <typename T, std::size_t NDIM>
void make_term_inputs_reconstructed(composite_term<T, NDIM>& term) {
    if (term.ket.is_initialized())     term.ket.change_tree_state(reconstructed, false);
    for (auto& p : term.p1)            p.change_tree_state(reconstructed, false);
    for (auto& p : term.p2)            p.change_tree_state(reconstructed, false);
    if (term.factor1.is_initialized()) term.factor1.change_tree_state(reconstructed, false);
    if (term.factor2.is_initialized()) term.factor2.change_tree_state(reconstructed, false);
}

/// Translate a composite_term<T,NDIM> into the op's CoeffTracker-based TermImpl.
template <typename T, std::size_t NDIM, typename OP>
typename OP::TermImpl make_term_impl(const composite_term<T, NDIM>& term) {
    constexpr std::size_t LDIM = NDIM / 2;
    typename OP::TermImpl x;
    if (term.ket.is_initialized()) {
        x.iaket = CoeffTracker<T, NDIM>(term.ket.get_impl().get());
    } else {
        MADNESS_ASSERT(!term.p1.empty() && term.p1.size() == term.p2.size());
        x.iap1.reserve(term.p1.size());
        x.iap2.reserve(term.p2.size());
        for (const auto& p : term.p1) x.iap1.emplace_back(p.get_impl().get());
        for (const auto& p : term.p2) x.iap2.emplace_back(p.get_impl().get());
    }
    if (term.factor1.is_initialized())
        x.ifactor1 = CoeffTracker<T, LDIM>(term.factor1.get_impl().get());
    if (term.factor2.is_initialized())
        x.ifactor2 = CoeffTracker<T, LDIM>(term.factor2.get_impl().get());
    if (term.extra.is_initialized()) x.extra = term.extra.get_impl().get();
    return x;
}

}  // namespace detail

/// build result = Σ_t product_t adaptively
///
/// All input functions are brought into redundant state (cheap from reconstructed).
/// The result is left in reconstructed state on exit.
template <typename T, std::size_t NDIM, typename leaf_opT>
void composite_product_apply(FunctionImpl<T, NDIM>* result,
                             const leaf_opT& leaf_op,
                             std::vector<composite_term<T, NDIM>>& terms,
                             double target_precision,
                             int oversampling,
                             bool fence) {
    constexpr std::size_t LDIM = NDIM / 2;
    using OP = composite_product_op_NS<T, NDIM, leaf_opT, LDIM>;

    World& world = result->world;
    for (auto& t : terms) detail::make_term_inputs_redundant(t);
    world.gop.fence();

    if (target_precision == 0.0) target_precision = result->get_thresh();

    std::vector<typename OP::TermImpl> impl_terms;
    impl_terms.reserve(terms.size());
    for (const auto& t : terms) {
        impl_terms.push_back(detail::make_term_impl<T, NDIM, OP>(t));
    }

    OP coeff_op(result, leaf_op, impl_terms, OP::Mode::function,
                /*bra*/ nullptr, target_precision, oversampling);
    noop<T, NDIM> apply_op;

    if (world.rank() == result->get_coeffs().owner(result->get_cdata().key0)) {
        result->task(world.rank(),
                     &FunctionImpl<T, NDIM>::template forward_traverse<OP, noop<T, NDIM>>,
                     coeff_op, apply_op, result->get_cdata().key0);
    }

    result->set_tree_state(redundant_after_merge);
    if (fence) {
        world.gop.fence();
        result->flo_unary_op_node_inplace(
            typename FunctionImpl<T, NDIM>::do_consolidate_buffer(result->get_tensor_args()),
            true);
        result->sum_down(true);
        result->set_tree_state(reconstructed);

        // Restore input tree state so callers can keep using the same Function handles
        // (e.g., to compare against the old multiply / make_Vphi paths that expect
        // reconstructed inputs).
        for (auto& t : terms) detail::make_term_inputs_reconstructed(t);
        world.gop.fence();
    }
}

/// compute <bra | Σ_t product_t> adaptively.
///
/// `bra` must be on-demand (e.g. an on-demand ERI).  Returns the scalar.
template <typename T, std::size_t NDIM, typename leaf_opT>
T composite_inner_apply(World& world,
                        const FunctionImpl<T, NDIM>* like,
                        const leaf_opT& leaf_op,
                        std::vector<composite_term<T, NDIM>>& terms,
                        const FunctionImpl<T, NDIM>* bra,
                        double target_precision,
                        int oversampling) {
    constexpr std::size_t LDIM = NDIM / 2;
    using OP = composite_product_op_NS<T, NDIM, leaf_opT, LDIM>;

    MADNESS_ASSERT(bra && bra->is_on_demand());

    for (auto& t : terms) {
        MADNESS_ASSERT(!t.extra.is_initialized() &&
                       "inner mode cannot mix per-term extra factor with global bra");
        detail::make_term_inputs_redundant(t);
    }
    world.gop.fence();

    // Build a same-shape accumulator tree (empty initially); the scalar accumulates at key0.
    FunctionImpl<T, NDIM> acc(*like, like->get_pmap(), false);
    if (target_precision == 0.0) target_precision = acc.get_thresh() * 0.01;

    std::vector<typename OP::TermImpl> impl_terms;
    impl_terms.reserve(terms.size());
    for (const auto& t : terms) {
        impl_terms.push_back(detail::make_term_impl<T, NDIM, OP>(t));
    }

    OP coeff_op(&acc, leaf_op, impl_terms, OP::Mode::inner,
                bra, target_precision, oversampling);
    noop<T, NDIM> apply_op;

    if (world.rank() == acc.get_coeffs().owner(acc.get_cdata().key0)) {
        acc.task(world.rank(),
                 &FunctionImpl<T, NDIM>::template forward_traverse<OP, noop<T, NDIM>>,
                 coeff_op, apply_op, acc.get_cdata().key0);
    }
    world.gop.fence();

    // consolidate the buffer at key0 and read the scalar
    acc.flo_unary_op_node_inplace(
        typename FunctionImpl<T, NDIM>::do_consolidate_buffer(acc.get_tensor_args()), true);

    T scalar = T(0);
    if (world.rank() == acc.get_coeffs().owner(acc.get_cdata().key0)) {
        typename FunctionImpl<T, NDIM>::dcT::const_accessor a;
        if (acc.get_coeffs().find(a, acc.get_cdata().key0))
            scalar = a->second.coeff().full_tensor_copy().flat()(0L);
    }
    world.gop.sum(scalar);

    // Restore input tree state for the caller.
    for (auto& t : terms) detail::make_term_inputs_reconstructed(t);
    world.gop.fence();

    return scalar;
}

}  // namespace madness

#endif  // MADNESS_MRA_COMPOSITE_PRODUCT_OP_H__
