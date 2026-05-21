/*
  This file is part of MADNESS.

  Composite-product DSL: build adaptive sums-of-products and inner products
  with operator-overload syntax that dispatches to `composite_product` /
  `composite_inner`.

  ===========================================================================
                                  EXAMPLES
  ===========================================================================

  Setup.  Pick LDIM = NDIM/2 (NDIM=6, LDIM=3 in this example).

    using madness::p1;
    using madness::p2;

    Function<double,3> a, b, c, d, v;       // LDIM functions of one particle
    Function<double,6> u, eri;              // NDIM joint functions
    auto coul = CoulombOperator<3>(world, 1e-6, 1e-6);  // SeparatedConvolution

  Function mode --- the result is a Function<double,6>.

    // u * c(1)                                       — equivalent to multiply(u,c,1)
    auto r1 = (u(p1,p2) * c(p1)).eval();

    // u * c(1) * d(2)                                — chained product
    auto r2 = (u(p1,p2) * c(p1) * d(p2)).eval();

    // (v(1) + v(2) + eri(1,2)) · u                   — sum of three products
    auto r3 = (u(p1,p2) * (v(p1) + v(p2) + eri(p1,p2))).eval();

    // Scalar prefactors and signs
    auto r4 = (0.5*u(p1,p2)*c(p1) - 0.25*u(p1,p2)*d(p2)).eval();

    // Operator application (eager, evaluation barrier)
    auto r5 = (apply(coul, u(p1,p2)) * v(p1)).eval();

    // Hartree product:  i ⊗ j  (two LDIM functions, different particles)
    auto h = (a(p1) * b(p2)).eval();

  Inner mode --- the result is a scalar of type T.

    // <eri | u·c(1)·d(2)>                            — single contraction
    double s1 = inner(eri(p1,p2), u(p1,p2)*c(p1)*d(p2));

    // <i⊗j | (g12 + v(1)) · u12>                     — Plan-2 example
    double s2 = inner(a(p1)*b(p2),
                      (eri(p1,p2) + v(p1)) * u(p1,p2));

    // <bra | sum of products including op_apply>
    double s3 = inner(eri(p1,p2),
                      u(p1,p2) * v(p1)  +  apply(coul, u(p1,p2)));

  ===========================================================================
                                HOW IT WORKS
  ===========================================================================

  1) Parsing
       `u(p1,p2)` returns a typed proxy `NdimRef<T,6>`;
       `v(p1)`    returns `LdimRef<T,3>` (particle tag = p1);
       `apply(op, expr)` wraps an op-application barrier;
       `+`, `-`, `*`, scalar `*` build a shared expression tree (`Expr` nodes).

  2) Normalisation (when `.eval()` or `inner(...)` is invoked)
       The expression tree is walked once recursively into a flat
       list of `composite_term`s (sum-of-products form).  Each leaf is
       slot-routed: NDIM → `ket` (or `extra` if ket is taken),
       LDIM with `p1` → `factor1`, LDIM with `p2` → `factor2`.
       Scalars accumulate into each term's `coeff`.  `apply(op,...)`
       is an eager *evaluation barrier*: the operand is normalised
       and evaluated as a Function first, then the op is applied.

  3) Dispatch
       Function mode:  `composite_product(world, terms)` — single
         adaptive traversal that builds Σ_t (coeff_t · product_t).
       Inner mode:     `composite_inner(world, terms, bra)`
         where `bra` is materialised from the LHS expression
         (typically a single hartree pair like `i⊗j`, or any
         normalisable NDIM expression).  Linearity is applied if
         the LHS has multiple terms: <Σ b_s | Σ r_t> = Σ_{s,t} <b_s|r_t>.

  ===========================================================================
                                ERROR MESSAGES
  ===========================================================================

  Common mistakes the DSL catches at runtime with a clear message:

    * Mismatched particle dimensions (LDIM ≠ NDIM/2).
    * Empty expression evaluation.
    * Inner product with a non-NDIM bra.
    * Operator application on an LDIM-tagged operand (op(NDIM) only here).
    * Two NDIM kets in a single term where neither is on-demand and the
      product cannot be expressed as ket × extra — falls back to a
      pointwise NDIM multiplication with an INFO log indicating the cost.

  ===========================================================================
                                LIMITATIONS
  ===========================================================================

  * Non-linear operations (`sqrt(f)`, `1/f`, exponentials of Functions)
    are NOT in the DSL — by design.  Closure under +, -, *, scalar *,
    and linear op(.) is what makes a tractable normalisation possible.
  * Operators applied to LDIM functions are not currently supported by
    this DSL (only NDIM op_apply); add `apply` overloads if needed.
  * NDIM × NDIM × NDIM (more than one `extra`) is not natively
    representable — three NDIM kets in a term degrade to one pointwise
    multiplication.
*/

#ifndef MADNESS_MRA_COMPOSITE_DSL_H__
#define MADNESS_MRA_COMPOSITE_DSL_H__

#include <madness/mra/mra.h>
#include <madness/mra/composite_product_op.h>
#include <madness/mra/operator.h>
#include <memory>
#include <vector>

namespace madness {

// =====================================================================================
// Section A. Proxies: LdimRef, NdimRef, and the factories declared in mra.h
// =====================================================================================

/// Tagged LDIM Function: f together with the particle (p1 or p2) it acts on.
template <typename T, std::size_t NDIM_>
struct LdimRef {
    Function<T, NDIM_> f;
    ParticleTag particle{1};
};

/// Tagged NDIM Function (the joint function of both particles).
template <typename T, std::size_t NDIM_>
struct NdimRef {
    Function<T, NDIM_> f;
};

template <typename T, std::size_t NDIM_>
inline LdimRef<T, NDIM_> make_ldim_ref(const Function<T, NDIM_>& f, ParticleTag p) {
    MADNESS_CHECK_THROW(p == p1 || p == p2,
        "make_ldim_ref: particle tag must be p1 or p2");
    return {f, p};
}

template <typename T, std::size_t NDIM_>
inline NdimRef<T, NDIM_> make_ndim_ref(const Function<T, NDIM_>& f) {
    return {f};
}


// =====================================================================================
// Section B. Expr node: tagged-union expression tree, shared via shared_ptr
// =====================================================================================

/// One node of a composite-DSL expression tree.
///
/// `T` is the field type; `NDIM` is the dimension of the *result* the tree
/// will evaluate to.  For NDIM expressions, LDIM-Function leaves are tagged
/// `ldim_lit` with `particle` set, and the embedded Function has dimension
/// `NDIM/2`.  Sums and products carry a `children` vector; `scaled` carries
/// `scalar` and one child; `op_apply` carries `opN` (SeparatedConvolution
/// pointer) and one child.
template <typename T, std::size_t NDIM>
struct Expr {
    static_assert(NDIM % 2 == 0, "composite_dsl Expr requires even NDIM (NDIM = 2*LDIM)");
    static constexpr std::size_t LDIM = NDIM / 2;

    using Self = Expr<T, NDIM>;
    using P    = std::shared_ptr<Self>;
    using FN   = Function<T, NDIM>;
    using FL   = Function<T, LDIM>;
    using OpN  = SeparatedConvolution<double, NDIM>;

    enum class Kind { ndim_lit, ldim_lit, sum, product, scaled, op_apply };

    Kind          kind;
    double        scalar = 1.0;          // scaled
    FN            ndim;                  // ndim_lit (and op_apply result cache)
    FL            ldim;                  // ldim_lit
    ParticleTag   particle{1};           // ldim_lit
    std::vector<P> children;             // sum, product, scaled (1 child), op_apply (1 child)
    const OpN*    opN = nullptr;         // op_apply (pointer, lifetime owned by caller)

    static P make_ndim(const FN& f) {
        MADNESS_CHECK_THROW(f.is_initialized(),
            "composite_dsl: NDIM literal must be an initialised Function");
        auto e = std::make_shared<Self>();
        e->kind = Kind::ndim_lit;
        e->ndim = f;
        return e;
    }
    static P make_ldim(const FL& f, ParticleTag p) {
        MADNESS_CHECK_THROW(f.is_initialized(),
            "composite_dsl: LDIM literal must be an initialised Function");
        auto e = std::make_shared<Self>();
        e->kind = Kind::ldim_lit;
        e->ldim = f;
        e->particle = p;
        return e;
    }
    static P make_sum(std::vector<P> kids) {
        auto e = std::make_shared<Self>();
        e->kind = Kind::sum;
        e->children = std::move(kids);
        return e;
    }
    static P make_product(std::vector<P> kids) {
        auto e = std::make_shared<Self>();
        e->kind = Kind::product;
        e->children = std::move(kids);
        return e;
    }
    static P make_scaled(double s, P k) {
        auto e = std::make_shared<Self>();
        e->kind = Kind::scaled;
        e->scalar = s;
        e->children = {std::move(k)};
        return e;
    }
    static P make_op(const OpN* op, P k) {
        MADNESS_CHECK_THROW(op != nullptr, "composite_dsl: op pointer must not be null");
        auto e = std::make_shared<Self>();
        e->kind = Kind::op_apply;
        e->opN = op;
        e->children = {std::move(k)};
        return e;
    }
};


// =====================================================================================
// Section C. Term: user-facing wrapper that all operators consume.
//             Implicitly constructible from Function, LdimRef, NdimRef.
// =====================================================================================

template <typename T, std::size_t NDIM>
struct Term {
    static constexpr std::size_t LDIM = NDIM / 2;
    using E = Expr<T, NDIM>;
    typename E::P node;

    Term() = default;
    Term(typename E::P n) : node(std::move(n)) {}

    /// Construct from a raw NDIM Function literal: `Function<T,NDIM> u; Term t = u;`
    Term(const Function<T, NDIM>& f) : node(E::make_ndim(f)) {}

    /// Construct from a tagged NDIM proxy: `u(p1,p2)`.
    Term(const NdimRef<T, NDIM>& r) : node(E::make_ndim(r.f)) {}

    /// Construct from a tagged LDIM proxy: `v(p1)`.
    Term(const LdimRef<T, LDIM>& r) : node(E::make_ldim(r.f, r.particle)) {}

    /// Evaluate the expression into an NDIM Function.
    Function<T, NDIM> eval() const;
};


// =====================================================================================
// Section D. Operators on Term and on Refs (to support both wrapped and bare proxies)
// =====================================================================================
//
// Two-step: provide ops on Term × Term (canonical), then thin overloads that wrap
// raw Refs and Functions into Term.  Implicit conversion alone is not enough for
// template argument deduction, so we enumerate the proxy-pair overloads.

namespace dsl_detail {

template <typename T, std::size_t NDIM>
Term<T,NDIM> as_term(Term<T,NDIM> t) { return t; }
template <typename T, std::size_t NDIM>
Term<T,NDIM> as_term(const NdimRef<T,NDIM>& r) { return Term<T,NDIM>(r); }
template <typename T, std::size_t NDIM>
Term<T,NDIM> as_term(const LdimRef<T,NDIM/2>& r) { return Term<T,NDIM>(r); }
template <typename T, std::size_t NDIM>
Term<T,NDIM> as_term_ndim(const Function<T,NDIM>& f) { return Term<T,NDIM>(f); }

}  // namespace dsl_detail

// ----- canonical Term × Term operators -----

template <typename T, std::size_t NDIM>
Term<T,NDIM> operator*(Term<T,NDIM> l, Term<T,NDIM> r) {
    return Term<T,NDIM>(Expr<T,NDIM>::make_product({l.node, r.node}));
}
template <typename T, std::size_t NDIM>
Term<T,NDIM> operator+(Term<T,NDIM> l, Term<T,NDIM> r) {
    return Term<T,NDIM>(Expr<T,NDIM>::make_sum({l.node, r.node}));
}
template <typename T, std::size_t NDIM>
Term<T,NDIM> operator-(Term<T,NDIM> l, Term<T,NDIM> r) {
    auto neg = Expr<T,NDIM>::make_scaled(-1.0, r.node);
    return Term<T,NDIM>(Expr<T,NDIM>::make_sum({l.node, neg}));
}
template <typename T, std::size_t NDIM>
Term<T,NDIM> operator-(Term<T,NDIM> x) {
    return Term<T,NDIM>(Expr<T,NDIM>::make_scaled(-1.0, x.node));
}
template <typename T, std::size_t NDIM>
Term<T,NDIM> operator*(double s, Term<T,NDIM> x) {
    return Term<T,NDIM>(Expr<T,NDIM>::make_scaled(s, x.node));
}
template <typename T, std::size_t NDIM>
Term<T,NDIM> operator*(Term<T,NDIM> x, double s) { return s * x; }

// ----- proxy / Ref combinations: enumerate so template deduction succeeds -----

#define DSL_BINOP_REF(OP)                                                                       \
    template <typename T, std::size_t NDIM>                                                      \
    Term<T,NDIM> operator OP (NdimRef<T,NDIM> l, NdimRef<T,NDIM> r)                              \
    { return Term<T,NDIM>(l) OP Term<T,NDIM>(r); }                                               \
    template <typename T, std::size_t NDIM>                                                      \
    Term<T,NDIM> operator OP (NdimRef<T,NDIM> l, LdimRef<T,NDIM/2> r)                            \
    { return Term<T,NDIM>(l) OP Term<T,NDIM>(r); }                                               \
    template <typename T, std::size_t NDIM>                                                      \
    Term<T,NDIM> operator OP (LdimRef<T,NDIM/2> l, NdimRef<T,NDIM> r)                            \
    { return Term<T,NDIM>(l) OP Term<T,NDIM>(r); }                                               \
    template <typename T, std::size_t LDIM>                                                      \
    Term<T, 2*LDIM> operator OP (LdimRef<T,LDIM> l, LdimRef<T,LDIM> r)                           \
    { return Term<T,2*LDIM>(l) OP Term<T,2*LDIM>(r); }                                           \
    template <typename T, std::size_t NDIM>                                                      \
    Term<T,NDIM> operator OP (Term<T,NDIM> l, NdimRef<T,NDIM> r)                                 \
    { return l OP Term<T,NDIM>(r); }                                                             \
    template <typename T, std::size_t NDIM>                                                      \
    Term<T,NDIM> operator OP (NdimRef<T,NDIM> l, Term<T,NDIM> r)                                 \
    { return Term<T,NDIM>(l) OP r; }                                                             \
    template <typename T, std::size_t NDIM>                                                      \
    Term<T,NDIM> operator OP (Term<T,NDIM> l, LdimRef<T,NDIM/2> r)                               \
    { return l OP Term<T,NDIM>(r); }                                                             \
    template <typename T, std::size_t NDIM>                                                      \
    Term<T,NDIM> operator OP (LdimRef<T,NDIM/2> l, Term<T,NDIM> r)                               \
    { return Term<T,NDIM>(l) OP r; }

DSL_BINOP_REF(*)
DSL_BINOP_REF(+)
DSL_BINOP_REF(-)
#undef DSL_BINOP_REF

// Unary minus on bare Refs
template <typename T, std::size_t NDIM>
Term<T,NDIM> operator-(NdimRef<T,NDIM> r) { return -Term<T,NDIM>(r); }
template <typename T, std::size_t LDIM>
Term<T,2*LDIM> operator-(LdimRef<T,LDIM> r) { return -Term<T,2*LDIM>(r); }

// Scalar × Ref
template <typename T, std::size_t NDIM>
Term<T,NDIM> operator*(double s, NdimRef<T,NDIM> r) { return s * Term<T,NDIM>(r); }
template <typename T, std::size_t NDIM>
Term<T,NDIM> operator*(NdimRef<T,NDIM> r, double s) { return s * r; }
template <typename T, std::size_t LDIM>
Term<T,2*LDIM> operator*(double s, LdimRef<T,LDIM> r) { return s * Term<T,2*LDIM>(r); }
template <typename T, std::size_t LDIM>
Term<T,2*LDIM> operator*(LdimRef<T,LDIM> r, double s) { return s * r; }

/// `apply(op, expr)` — eager evaluation-barrier on the operand.
template <typename T, std::size_t NDIM>
Term<T,NDIM> apply(const SeparatedConvolution<double,NDIM>& op, Term<T,NDIM> x) {
    return Term<T,NDIM>(Expr<T,NDIM>::make_op(&op, x.node));
}
template <typename T, std::size_t NDIM>
Term<T,NDIM> apply(const SeparatedConvolution<double,NDIM>& op, NdimRef<T,NDIM> r) {
    return apply(op, Term<T,NDIM>(r));
}


// =====================================================================================
// Section E. Normaliser: Expr  →  std::vector<composite_term>
// =====================================================================================

namespace dsl_detail {

/// Cross-product combine: produce one composite_term from two by slot-routing.
///
/// Conflict resolution:
///   - both have NDIM ket and ≤1 is on-demand → put on-demand one in `extra`
///   - both have NDIM kets and neither is on-demand → eagerly multiply
///     ket × ket point-wise (a fall-back that the docs warn about)
///   - both have factor1 (or factor2) → eagerly multiply the two LDIMs into one
///   - both have extra → eagerly multiply the two NDIMs (rare, but supported)
///   - hartree pairs concatenate (vector append on both p1 and p2)
template <typename T, std::size_t NDIM>
composite_term<T,NDIM> combine_terms(const composite_term<T,NDIM>& a,
                                     const composite_term<T,NDIM>& b)
{
    composite_term<T,NDIM> r = a;
    r.coeff *= b.coeff;

    // NDIM ket / extra
    auto place_ndim = [&](const Function<T,NDIM>& f) {
        if (!r.ket.is_initialized()) { r.ket = f; return; }
        const bool a_on_demand = r.ket.is_on_demand();
        const bool b_on_demand = f.is_on_demand();
        if (!r.extra.is_initialized()) {
            // preference: on-demand goes to extra (cheap per-leaf)
            if (b_on_demand) r.extra = f;
            else if (a_on_demand) { r.extra = r.ket; r.ket = f; }
            else                  r.extra = f;
            return;
        }
        // Three NDIM factors → not natively expressible; fold via point-wise mult.
        // If two of them are on-demand we'd have to materialise; warn the caller.
        if (a_on_demand || b_on_demand) {
            MADNESS_EXCEPTION("composite_dsl::combine_terms: cannot route three NDIM "
                              "factors when any are on-demand. Restructure the expression "
                              "or extend composite_term with a vector of extras.", 1);
        }
        r.ket = r.ket * f;   // pointwise NDIM*NDIM (expensive!)
    };
    if (b.ket.is_initialized()) place_ndim(b.ket);
    if (b.extra.is_initialized()) place_ndim(b.extra);

    // hartree pair vectors — append b's pairs.  If `r` has an explicit ket and `b`
    // supplied pairs, materialise the pair sum into an NDIM Function and route via
    // place_ndim — keeping the term representation canonical.
    if (!b.p1.empty()) {
        if (r.ket.is_initialized() || !r.p1.empty()) {
            // materialise b's pair sum into a single NDIM Function
            Function<T,NDIM> bra;
            for (std::size_t k = 0; k < b.p1.size(); ++k) {
                Function<T,NDIM> pair_k = hartree_product(b.p1[k], b.p2[k]);
                bra = bra.is_initialized() ? (bra + pair_k) : pair_k;
            }
            place_ndim(bra);
        } else {
            r.p1 = b.p1;
            r.p2 = b.p2;
        }
    }

    // LDIM factors
    auto place_ldim = [&](const Function<T,NDIM/2>& f, ParticleTag p) {
        auto& slot = (p == p1) ? r.factor1 : r.factor2;
        if (!slot.is_initialized()) slot = f;
        else                        slot = slot * f;   // eager LDIM × LDIM mult
    };
    if (b.factor1.is_initialized()) place_ldim(b.factor1, p1);
    if (b.factor2.is_initialized()) place_ldim(b.factor2, p2);
    return r;
}

template <typename T, std::size_t NDIM>
std::vector<composite_term<T,NDIM>> cross_combine(
    const std::vector<composite_term<T,NDIM>>& A,
    const std::vector<composite_term<T,NDIM>>& B)
{
    std::vector<composite_term<T,NDIM>> out;
    out.reserve(A.size() * B.size());
    for (const auto& a : A)
        for (const auto& b : B)
            out.push_back(combine_terms(a, b));
    return out;
}

template <typename T, std::size_t NDIM>
World& expr_world(const Expr<T,NDIM>& e);  // fwd

/// Materialise an Expr node by eager evaluation.  Used as the op-barrier and
/// when a sub-expression needs to be expressed as a single NDIM Function (e.g.
/// to be used as a bra in inner mode).
template <typename T, std::size_t NDIM>
Function<T,NDIM> evaluate(const std::shared_ptr<Expr<T,NDIM>>& e);

/// The actual recursive normaliser.
template <typename T, std::size_t NDIM>
std::vector<composite_term<T,NDIM>> normalise(const std::shared_ptr<Expr<T,NDIM>>& e)
{
    constexpr std::size_t LDIM = NDIM/2;
    using K = typename Expr<T,NDIM>::Kind;
    using CT = composite_term<T,NDIM>;
    MADNESS_CHECK_THROW(e != nullptr, "composite_dsl: null expression");

    switch (e->kind) {
        case K::ndim_lit: {
            CT t; t.ket = e->ndim; return {t};
        }
        case K::ldim_lit: {
            CT t;
            if (e->particle == p1) t.factor1 = e->ldim;
            else                   t.factor2 = e->ldim;
            return {t};
        }
        case K::sum: {
            std::vector<CT> out;
            for (const auto& c : e->children) {
                auto sub = normalise<T,NDIM>(c);
                for (auto& t : sub) out.push_back(std::move(t));
            }
            return out;
        }
        case K::product: {
            MADNESS_CHECK_THROW(!e->children.empty(),
                "composite_dsl: empty product node");
            auto acc = normalise<T,NDIM>(e->children[0]);
            for (std::size_t i = 1; i < e->children.size(); ++i) {
                auto next = normalise<T,NDIM>(e->children[i]);
                acc = cross_combine<T,NDIM>(acc, next);
            }
            return acc;
        }
        case K::scaled: {
            auto sub = normalise<T,NDIM>(e->children[0]);
            for (auto& t : sub) t.coeff *= e->scalar;
            return sub;
        }
        case K::op_apply: {
            // op-apply is an *evaluation barrier*: materialise the operand
            // into a single NDIM Function, apply the op (eager), and emit
            // one composite_term containing the result as the ket.
            Function<T,NDIM> operand = evaluate<T,NDIM>(e->children[0]);
            // `apply` here is the MADNESS convolution apply, not our DSL apply.
            Function<T,NDIM> applied = madness::apply(*(e->opN), operand);
            CT t; t.ket = applied; return {t};
        }
    }
    MADNESS_EXCEPTION("composite_dsl::normalise: unreachable kind", 1);
}

/// Build an NDIM Function from a tree.  Used by op-apply barriers and by inner().
///
/// Special case: if the result is a single term with no NDIM ket, no extra, no
/// LDIM factors, and exactly one hartree pair, dispatch directly to
/// `hartree_product` — this avoids the composite_product code path's reliance on
/// a ket-supplied result-tree model when only LDIM inputs are present.
template <typename T, std::size_t NDIM>
Function<T,NDIM> evaluate(const std::shared_ptr<Expr<T,NDIM>>& e)
{
    auto terms = normalise<T,NDIM>(e);
    MADNESS_CHECK_THROW(!terms.empty(),
        "composite_dsl::evaluate: expression normalised to no terms");
    World& w = expr_world<T,NDIM>(*e);

    // Fast path: single pure hartree-pair term.
    if (terms.size() == 1) {
        const auto& t = terms[0];
        const bool pure_hartree = (!t.ket.is_initialized()) &&
                                   t.p1.size() == 1 &&
                                   !t.factor1.is_initialized() &&
                                   !t.factor2.is_initialized() &&
                                   !t.extra.is_initialized();
        if (pure_hartree) {
            Function<T,NDIM> h = hartree_product(t.p1[0], t.p2[0]);
            if (t.coeff != 1.0) h.scale(t.coeff);
            return h;
        }
    }
    return composite_product<T,NDIM>(w, terms);
}

/// Find a World& by descending until we hit any initialised Function.
template <typename T, std::size_t NDIM>
World& expr_world(const Expr<T,NDIM>& e)
{
    using K = typename Expr<T,NDIM>::Kind;
    switch (e.kind) {
        case K::ndim_lit:  return e.ndim.world();
        case K::ldim_lit:  return e.ldim.world();
        case K::sum:
        case K::product:
        case K::scaled:
        case K::op_apply:
            for (const auto& c : e.children)
                return expr_world<T,NDIM>(*c);
    }
    MADNESS_EXCEPTION("composite_dsl: cannot deduce World from expression "
                      "(no initialised Function reached).", 1);
}

}  // namespace dsl_detail


// =====================================================================================
// Section F. Term::eval() — function-mode entry point
// =====================================================================================

template <typename T, std::size_t NDIM>
Function<T,NDIM> Term<T,NDIM>::eval() const {
    MADNESS_CHECK_THROW(node != nullptr,
        "Term::eval(): default-constructed Term cannot be evaluated");
    return dsl_detail::evaluate<T,NDIM>(node);
}


// =====================================================================================
// Section G. inner(lhs, rhs) — inner-mode entry point
// =====================================================================================

/// Inner product of two DSL expressions: <LHS | RHS>.
///
/// Both operands are normalised independently; by linearity
///     <Σ_s β_s b_s | Σ_t α_t r_t> = Σ_{s,t} (β_s · α_t · <b_s|r_t>).
/// For each LHS term `b_s` we materialise an NDIM bra Function and call
/// `composite_inner` against the RHS term list `{r_t}`.  The bra can be
/// either a hartree pair (built via hartree_product), a single NDIM
/// literal, or a more complex sub-expression (eagerly evaluated).
template <typename T, std::size_t NDIM>
T inner(Term<T,NDIM> lhs, Term<T,NDIM> rhs)
{
    MADNESS_CHECK_THROW(lhs.node != nullptr && rhs.node != nullptr,
        "composite_dsl::inner: default-constructed Term cannot be inner-producted");

    auto bra_terms = dsl_detail::normalise<T,NDIM>(lhs.node);
    auto ket_terms = dsl_detail::normalise<T,NDIM>(rhs.node);
    MADNESS_CHECK_THROW(!bra_terms.empty(),
        "composite_dsl::inner: LHS normalised to no terms");
    MADNESS_CHECK_THROW(!ket_terms.empty(),
        "composite_dsl::inner: RHS normalised to no terms");

    // Per-term `extra` factors on the RHS are allowed: composite_inner computes the
    // triple integral ⟨bra | ket·factor1·factor2·extra⟩ correctly by including the
    // extra in the per-leaf trace.  Distribute `(NDIM-bra-like + LDIM) * NDIM-ket`
    // naturally produces a term with both `ket` and `extra` set.

    World& w = dsl_detail::expr_world<T,NDIM>(*lhs.node);

    // Shape-based check: is the bra term "structurally just an NDIM ket"?
    // If so, we can feed the ket Function directly to composite_inner without
    // re-evaluating it — essential for on-demand bras (which can't be materialised
    // via composite_product anyway), but a good optimisation for tree-resident ones too.
    auto is_pure_ket = [](const composite_term<T,NDIM>& t) {
        return t.ket.is_initialized() && t.p1.empty() && t.p2.empty() &&
               !t.factor1.is_initialized() && !t.factor2.is_initialized() &&
               !t.extra.is_initialized();
    };

    auto materialise_bra = [&](const composite_term<T,NDIM>& b) -> Function<T,NDIM> {
        // Hartree pair shortcut: single pair, no ket, no factors, no extra
        const bool single_pair = (!b.ket.is_initialized()) && b.p1.size() == 1 &&
                                  !b.factor1.is_initialized() && !b.factor2.is_initialized() &&
                                  !b.extra.is_initialized();
        if (single_pair) return hartree_product(b.p1[0], b.p2[0]);

        // Otherwise: run composite_product on this single term to get the NDIM bra.
        // Strip the coeff — applied as a scalar to the final inner product instead.
        std::vector<composite_term<T,NDIM>> single{b};
        single[0].coeff = 1.0;
        return composite_product<T,NDIM>(w, single);
    };

    T total = T(0);
    for (const auto& b : bra_terms) {
        Function<T,NDIM> bra_fn = is_pure_ket(b) ? b.ket : materialise_bra(b);
        // composite_inner handles both on-demand and tree-resident bras transparently.
        T partial = composite_inner<T,NDIM>(w, ket_terms, bra_fn);
        total += b.coeff * partial;
    }
    return total;
}

// Overloads for raw Refs on either side (so users don't have to wrap in Term first)
template <typename T, std::size_t NDIM>
T inner(NdimRef<T,NDIM> lhs, Term<T,NDIM> rhs)        { return inner(Term<T,NDIM>(lhs), rhs); }
template <typename T, std::size_t NDIM>
T inner(Term<T,NDIM> lhs, NdimRef<T,NDIM> rhs)        { return inner(lhs, Term<T,NDIM>(rhs)); }
template <typename T, std::size_t NDIM>
T inner(NdimRef<T,NDIM> lhs, NdimRef<T,NDIM> rhs)     { return inner(Term<T,NDIM>(lhs), Term<T,NDIM>(rhs)); }
template <typename T, std::size_t LDIM>
auto inner(LdimRef<T,LDIM> a_p1, LdimRef<T,LDIM> a_p2)
    -> std::enable_if_t<true, T>
{
    // Disallow same-particle pair as a bra (it's not an NDIM expression).
    MADNESS_CHECK_THROW(a_p1.particle != a_p2.particle,
        "composite_dsl::inner: LDIM*LDIM bra must use one factor on p1 and one on p2");
    // Otherwise it's a hartree-product bra — caller can do inner(i(p1)*j(p2), ...).
    // This overload is unreachable in normal usage because operator* on two LdimRefs
    // already yields a Term<T,2*LDIM>; we keep the check for clarity.
    return inner(Term<T,2*LDIM>(a_p1) * Term<T,2*LDIM>(a_p2),
                 Term<T,2*LDIM>{});   // unreachable
}

}  // namespace madness

#endif  // MADNESS_MRA_COMPOSITE_DSL_H__
