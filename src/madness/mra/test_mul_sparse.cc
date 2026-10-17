/// Unit tests for the screened multiplication mul_sparse.
///
/// The operand pairs have closed-form products, so the reference is an
/// independent projection of the analytic product rather than another call to
/// the kernel under test.

#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/world/test_utilities.h>
#include <madness/world/vector_archive.h>
#include <madness/mra/legendre.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace madness;

static const std::size_t D = 3;

/// projection threshold of the reference
static const double THRESH_REF = 1.e-9;

/// Error floor of the exact path (tol=0): the product of two order-k polynomials
/// is not representable at order k, and multiplication does not autorefine.
/// Measured max 2.3e-8 over the pairs below.
static const double EXACT_FLOOR = 1.e-7;

/// Bounds on the screening error over the tol ladder, all measured, see the
/// header of test_screening_accuracy.
static const double C_MAX = 100.0;   ///< err <= C_MAX*tol; measured max 34
static const double P_MIN = 0.4;     ///< err ~ tol^p; measured p = 0.58 .. 0.75
/// T19: ||q(tol) - q(0)|| <= C_TOL0 tol at tol 1e-10 .. 1e-14; measured max 147 (today's criterion, gauss-offset
/// at 1e-14; the Hoelder form 22): the guard is against losing orders of magnitude, as July 2026 did
static const double C_TOL0 = 1000.0;

/// exp(-a |r-R|^2)
class Gauss : public FunctionFunctorInterface<double, D> {
    double a;
    Vector<double, D> R;
public:
    Gauss(double a, const Vector<double, D>& R) : a(a), R(R) {}
    double operator()(const Vector<double, D>& r) const override {
        double r2 = 0.0;
        for (std::size_t i = 0; i < D; ++i) r2 += (r[i]-R[i]) * (r[i]-R[i]);
        return std::exp(-a * r2);
    }
    std::vector<Vector<double,D>> special_points() const override { return {R}; }
};

/// c * exp(-a |r|), cusped at the origin: this is what drives the tree deep
class Slater : public FunctionFunctorInterface<double, D> {
    double a, c;
public:
    Slater(double a, double c=1.0) : a(a), c(c) {}
    double operator()(const Vector<double, D>& r) const override {
        return c * std::exp(-a * r.normf());
    }
    std::vector<Vector<double,D>> special_points() const override {
        return {Vector<double,D>(0.0)};
    }
};

/// the criterion of mulXXveca under test (MAD_MUL_SCREEN, set by main for T1-T6), tagging titles and lines
static std::string screen_tag() { return mul_screen() == MulScreen::hoelder ? " [hoelder]" : " [cs]"; }

static std::string fmt_tol(double tol) {
    char buf[32];
    snprintf(buf, sizeof(buf), "%.0e", tol);
    return std::string(buf);
}

/// an operand pair plus the independently projected analytic product
struct operand_pair {
    std::string name;
    Function<double,D> f, g, fg;
};

static Function<double,D> project(World& world,
                                  const std::shared_ptr<FunctionFunctorInterface<double,D>>& functor) {
    Function<double,D> f = FunctionFactory<double,D>(world).functor(functor);
    f.truncate();
    return f;
}

/// coincident Gaussians, separated Gaussians (where the screening actually
/// fires, as it does for distant localized orbitals), and a Slater cusp
static std::vector<operand_pair> make_pairs(World& world) {
    typedef std::shared_ptr<FunctionFunctorInterface<double,D>> functorT;
    const Vector<double,D> origin(0.0);
    Vector<double,D> Rl(0.0), Rr(0.0);
    Rl[0] = -1.5;
    Rr[0] =  1.5;

    // exp(-a|r-Rl|^2) exp(-a|r-Rr|^2) = exp(-a|Rl-Rr|^2/2) exp(-2a|r-(Rl+Rr)/2|^2)
    const double a = 1.0;
    const double prefac = std::exp(-0.5 * a * (Rl[0]-Rr[0]) * (Rl[0]-Rr[0]));

    std::vector<operand_pair> pairs;
    pairs.push_back({"gauss",
                     project(world, functorT(new Gauss(a, origin))),
                     project(world, functorT(new Gauss(a, origin))),
                     project(world, functorT(new Gauss(2*a, origin)))});
    pairs.push_back({"gauss-offset",
                     project(world, functorT(new Gauss(a, Rl))),
                     project(world, functorT(new Gauss(a, Rr))),
                     project(world, functorT(new Gauss(2*a, origin)))});
    pairs.back().fg.scale(prefac);
    pairs.push_back({"slater",
                     project(world, functorT(new Slater(3.0))),
                     project(world, functorT(new Slater(3.0))),
                     project(world, functorT(new Slater(6.0)))});
    return pairs;
}

/// number of nodes violating the reconstructed layout: coefficients live on the
/// leaves only
static long interior_nodes_with_coeffs(World& world, const Function<double,D>& f) {
    long bad = 0;
    for (const auto& datum : f.get_impl()->get_coeffs()) {
        const auto& node = datum.second;
        if (node.has_children() and node.has_coeff()) ++bad;
    }
    world.gop.sum(bad);
    return bad;
}

/// T1: the exact path reproduces the analytic product
int test_exact(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T1" + screen_tag() + ": tol=0 vs an independent projection");
    t.set_do_print(world.rank() == 0);
    for (auto& p : pairs) {
        Function<double,D> q = mul_sparse(p.f, p.g, 0.0);
        const double err = (q - p.fg).norm2();
        if (world.rank() == 0)
            printf("  T1%s %-13s err %.3e   bound %.1e\n", screen_tag().c_str(), p.name.c_str(), err, EXACT_FLOOR);
        t.checkpoint(err, EXACT_FLOOR, "T1" + screen_tag() + " " + p.name);
    }
    return t.end();
}

/// T2: the result is a reconstructed tree, both scalar and vector entry points
int test_tree_state(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T2" + screen_tag() + ": output tree state");
    t.set_do_print(world.rank() == 0);
    operand_pair& p = pairs.front();

    for (double tol : {0.0, 1.e-5}) {
        Function<double,D> q = mul_sparse(p.f, p.g, tol);
        world.gop.fence();
        const long bad = interior_nodes_with_coeffs(world, q);
        t.checkpoint(get_tree_state(q) == reconstructed,
                     "T2" + screen_tag() + " scalar state is reconstructed, tol=" + fmt_tol(tol));
        t.checkpoint(bad == 0L,
                     "T2" + screen_tag() + " scalar coefficients on leaves only, tol=" + fmt_tol(tol));
        t.checkpoint(q.get_impl()->verify_tree_state_local(),
                     "T2" + screen_tag() + " scalar tree state is self-consistent, tol=" + fmt_tol(tol));
    }

    // the vector entry point is a separate code path (vmulXX) and mislabeled its
    // results at one point during development
    std::vector<Function<double,D>> v = {pairs[0].g, pairs[1].g};
    std::vector<Function<double,D>> qv = mul_sparse(world, p.f, v, 1.e-5);
    world.gop.fence();
    bool state_ok = true, layout_ok = true;
    for (auto& q : qv) {
        state_ok = state_ok and (get_tree_state(q) == reconstructed)
                            and q.get_impl()->verify_tree_state_local();
        layout_ok = layout_ok and (interior_nodes_with_coeffs(world, q) == 0L);
    }
    t.checkpoint(state_ok, "T2" + screen_tag() + " vector state is reconstructed");
    t.checkpoint(layout_ok, "T2" + screen_tag() + " vector coefficients on leaves only");
    return t.end();
}

/// T3: the screening error over a tol ladder, against the analytic product.
///
/// The bounds are measured, not derived: err ~ tol^p with p short of 1, and tree
/// quantization (adjacent tol yielding the same tree) leaves err flat and inflates
/// err/tol by up to a decade. Hence a loose C_MAX plus two tighter statements --
/// err falls monotonically with tol, and its decay exponent does not collapse.
int test_screening_accuracy(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T3" + screen_tag() + ": screening error vs tol");
    t.set_do_print(world.rank() == 0);
    const std::vector<double> tols = {1.e-4, 1.e-5, 1.e-6, 1.e-7, 1.e-8};

    for (auto& p : pairs) {
        std::vector<double> errs;
        bool bounded = true;
        for (double tol : tols) {
            Function<double,D> q = mul_sparse(p.f, p.g, tol);
            const double err = (q - p.fg).norm2();
            errs.push_back(err);
            bounded = bounded and (err < C_MAX * tol);
            if (world.rank() == 0)
                printf("  T3%s %-13s tol %.1e  err %.3e  C %6.2f\n",
                       screen_tag().c_str(), p.name.c_str(), tol, err, err / tol);
        }
        // 5% slack absorbs the quantization plateaus
        bool monotone = true;
        for (std::size_t i = 1; i < errs.size(); ++i)
            monotone = monotone and (errs[i] <= 1.05 * errs[i-1]);

        const double decades = std::log10(tols.front() / tols.back());
        const double pfit = std::log10(errs.front() / errs.back()) / decades;
        if (world.rank() == 0)
            printf("  T3%s %-13s decay exponent p %.2f\n", screen_tag().c_str(), p.name.c_str(), pfit);

        t.checkpoint(bounded,  "T3" + screen_tag() + " " + p.name + " err < C*tol");
        t.checkpoint(monotone, "T3" + screen_tag() + " " + p.name + " err falls with tol");
        t.checkpoint(pfit > P_MIN, "T3" + screen_tag() + " " + p.name + " decay exponent");
    }
    return t.end();
}

/// T4: screened and exact multiplication agree to within tol, no analytic
/// reference needed
int test_self_consistency(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T4" + screen_tag() + ": screened vs exact");
    t.set_do_print(world.rank() == 0);

    for (auto& p : pairs) {
        Function<double,D> exact = mul_sparse(p.f, p.g, 0.0);
        for (double tol : {1.e-4, 1.e-6, 1.e-8}) {
            Function<double,D> q = mul_sparse(p.f, p.g, tol);
            const double err = (q - exact).norm2();
            const double C = err / tol;
            if (world.rank() == 0)
                printf("  T4%s %-13s tol %.1e  err %.3e  C %.2f\n",
                       screen_tag().c_str(), p.name.c_str(), tol, err, C);
            t.checkpoint(C < C_MAX, "T4" + screen_tag() + " " + p.name + " tol=" + fmt_tol(tol));
        }
    }
    return t.end();
}

/// T5: dot_sparse routes through the screened kernel; check it against the analytic sum
int test_dot_sparse(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T5" + screen_tag() + ": dot vs the analytic sum of products");
    t.set_do_print(world.rank() == 0);

    std::vector<Function<double,D>> a, b;
    Function<double,D> ref = copy(pairs[0].fg);
    a.push_back(pairs[0].f);  b.push_back(pairs[0].g);
    for (std::size_t i = 1; i < pairs.size(); ++i) {
        a.push_back(pairs[i].f);
        b.push_back(pairs[i].g);
        ref += pairs[i].fg;
    }

    for (double tol : {0.0, 1.e-6}) {
        Function<double,D> q = dot_sparse(world, a, b, tol, true, true);
        const double err = (q - ref).norm2();
        const double bound = (tol == 0.0) ? EXACT_FLOOR : C_MAX * tol;
        if (world.rank() == 0)
            printf("  T5%s dot_sparse tol %.1e  err %.3e  bound %.1e\n", screen_tag().c_str(), tol, err, bound);
        t.checkpoint(err, bound, "T5" + screen_tag() + " dot_sparse tol=" + fmt_tol(tol));
    }
    return t.end();
}

/// T6: mul_sparse leaves its operands redundant. The global norm2()/trace()
/// entry points must stay correct across that state.
int test_operand_state(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T6" + screen_tag() + ": operand state and the norms read off it");
    t.set_do_print(world.rank() == 0);
    // this test changes the state of both operands, and the
    // pairs are shared with every test that runs after it
    operand_pair p{pairs.front().name, copy(pairs.front().f),
                   copy(pairs.front().g), copy(pairs.front().fg)};

    // reference values, taken while the operands are still reconstructed
    p.f.reconstruct();
    p.g.reconstruct();
    const double norm_before    = p.f.norm2();
    const double trace_before   = p.f.trace();
    const double g_trace_before = p.g.trace();

    Function<double,D> q = mul_sparse(p.f, p.g, 1.e-6);
    world.gop.fence();

    t.checkpoint(get_tree_state(p.f) == redundant,
                 "T6" + screen_tag() + " left operand is left redundant");
    t.checkpoint(get_tree_state(p.g) == redundant,
                 "T6" + screen_tag() + " right operand is left redundant");

    // norm2()/trace() read the redundant tree instead of converting it, so both
    // can be measured on the same operand, and the operand a caller hands to
    // the next mul_sparse() is the one it had
    const double norm_after    = p.f.norm2();
    const double trace_after   = p.f.trace();
    const double g_trace_after = p.g.trace();
    if (world.rank() == 0)
        printf("  T6%s norm2  %.10e -> %.10e\n  T6 trace  %.10e -> %.10e\n",
               screen_tag().c_str(), norm_before, norm_after, g_trace_before, g_trace_after);

    t.checkpoint(norm_after,    norm_before,    1.e-10, "T6" + screen_tag() + " norm2 on a redundant function");
    t.checkpoint(trace_after,   trace_before,   1.e-10, "T6" + screen_tag() + " trace on a redundant function");
    t.checkpoint(g_trace_after, g_trace_before, 1.e-10, "T6" + screen_tag() + " trace on the other redundant operand");
    t.checkpoint(get_tree_state(p.f) == redundant and
                 get_tree_state(p.g) == redundant,
                 "T6" + screen_tag() + " norm2/trace leave the operand redundant");

    // the vector overload takes the same path and must agree
    t.checkpoint(get_tree_state(p.f) == redundant,
                 "T6" + screen_tag() + " norm2s operand is redundant on entry");
    std::vector<Function<double,D>> v = {p.f};
    t.checkpoint(norm2s(world, v)[0], norm_before, 1.e-10,
                 "T6" + screen_tag() + " vector norm2s agrees");
    t.checkpoint(get_tree_state(p.f) == redundant,
                 "T6" + screen_tag() + " norm2s leaves the operand redundant");

    // p.f has been through reconstructed -> redundant and back; its trace must
    // still be the value measured at the top
    p.f.reconstruct();
    t.checkpoint(p.f.trace(), trace_before, 1.e-10,
                 "T6" + screen_tag() + " trace round-trips through redundant form");
    return t.end();
}

/// T7: norm_tree means ||f|| on the box, in every tree state. compress()
/// clears leaf coefficients, and the leaf norm has to be read before that
/// happens or a zero propagates to the root.
int test_norm_tree_states(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T7: norm_tree across tree states");
    t.set_do_print(world.rank() == 0);

    Function<double,D> f = copy(pairs.front().f);
    f.reconstruct();
    const double exact = f.norm2();

    auto root_norm_tree = [&world](const Function<double,D>& g) {
        const Key<D> key0 = g.get_impl()->get_cdata().key0;
        double nt = 0.0;
        if (g.get_impl()->get_coeffs().probe(key0))
            nt = g.get_impl()->get_coeffs().find(key0).get()->second.get_norm_tree();
        world.gop.max(nt);
        return nt;
    };

    auto root_dnorm_tree = [&world](const Function<double,D>& g) {
        const Key<D> key0 = g.get_impl()->get_cdata().key0;
        double dnt = 0.0;
        if (g.get_impl()->get_coeffs().probe(key0))
            dnt = g.get_impl()->get_coeffs().find(key0).get()->second.get_dnorm_tree();
        world.gop.max(dnt);
        return dnt;
    };

    auto root_snorm = [&world](const Function<double,D>& g) {
        const Key<D> key0 = g.get_impl()->get_cdata().key0;
        double sn = 0.0;
        if (g.get_impl()->get_coeffs().probe(key0))
            sn = g.get_impl()->get_coeffs().find(key0).get()->second.get_snorm();
        world.gop.max(sn);
        return sn;
    };

    f.make_redundant(true);
    const double nt_redundant = root_norm_tree(f);
    f.reconstruct();
    f.compress();
    const double nt_compressed = root_norm_tree(f);

    if (world.rank() == 0)
        printf("  T7 ||f|| %.10e  redundant %.10e  compressed %.10e\n",
               exact, nt_redundant, nt_compressed);

    t.checkpoint(std::abs(nt_redundant  - exact) < 1.e-10 * exact,
                 "T7 root norm_tree on a redundant tree equals ||f||");
    t.checkpoint(std::abs(nt_compressed - exact) < 1.e-10 * exact,
                 "T7 root norm_tree on a compressed tree equals ||f||");

    // ||f||^2 = ||s||^2 + ||d||^2, so dnorm_tree at the root must be strictly
    // below ||f||: equality means the scaling block was counted as detail
    const double dnt_compressed = root_dnorm_tree(f);
    if (world.rank() == 0)
        printf("  T7 root dnorm_tree compressed %.10e  (||f|| %.10e)\n",
               dnt_compressed, exact);
    t.checkpoint(dnt_compressed < 0.999 * exact,
                 "T7 root dnorm_tree excludes the scaling block");

    // The multiwavelet basis is orthonormal, so the root scaling block and every
    // difference coefficient below it are mutually orthogonal:
    //   ||f||^2 == ||s_root||^2 + sum_nodes ||d||^2 == snorm^2 + dnorm_tree^2
    // This pins dnorm_tree rather than merely bounding it, and so also catches a
    // dnorm that comes out too small.
    const double sn_compressed = root_snorm(f);
    const double recomposed = std::sqrt(sn_compressed * sn_compressed
                                      + dnt_compressed * dnt_compressed);
    if (world.rank() == 0)
        printf("  T7 snorm %.14e  dnorm_tree %.14e\n"
               "  T7 sqrt(s^2+d^2) %.14e  ||f|| %.14e\n",
               sn_compressed, dnt_compressed, recomposed, exact);
    t.checkpoint(std::abs(recomposed - exact) < 1.e-12 * exact,
                 "T7 snorm^2 + dnorm_tree^2 == ||f||^2 at the root");
    return t.end();
}

/// T8: broaden() resets the per-node norms through zero_norm_tree(). The
/// screening criterion reads norm_tree and dnorm_tree as a pair, so a reset
/// that clears one and leaves the other in place is worse than no reset at
/// all: norm_tree = 0 with a small stale dnorm_tree passes the test in
/// mulXXveca() and multiplies at a level that is too coarse. sup_tree, the
/// Hoelder screen's bound, goes with them.
int test_broaden_resets_norms(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T8: broaden resets norm_tree and dnorm_tree");
    t.set_do_print(world.rank() == 0);

    Function<double,D> f = copy(pairs.front().f);
    f.reconstruct();

    // redundant form writes both norms on every node, and undoing it only
    // drops the internal coefficients -- so the tree carries computed norms
    // into broaden()
    f.make_redundant(true);
    f.reconstruct();

    // counts local nodes, summed over ranks: computed dnorm_tree, nonzero
    // norm_tree, computed sup_tree
    auto count = [&world](const Function<double,D>& g) {
        long computed = 0, nonzero = 0, sup = 0;
        const auto& c = g.get_impl()->get_coeffs();
        for (auto it = c.begin(); it != c.end(); ++it) {
            if (it->second.get_dnorm_tree() != NORM_TREE_UNCOMPUTED) ++computed;
            if (it->second.get_norm_tree() != 0.0) ++nonzero;
            if (it->second.get_sup_tree() != NORM_TREE_UNCOMPUTED) ++sup;
        }
        world.gop.sum(computed);
        world.gop.sum(nonzero);
        world.gop.sum(sup);
        return std::array<long, 3>{computed, nonzero, sup};
    };

    const auto before = count(f);
    if (world.rank() == 0)
        printf("  T8 before broaden: %ld nodes with a computed dnorm_tree, "
               "%ld with a nonzero norm_tree, %ld with a computed sup_tree\n", before[0], before[1], before[2]);
    // without this the test would pass on an empty premise
    t.checkpoint(before[0] > 0 and before[1] > 0 and before[2] > 0,
                 "T8 the tree carries computed norms into broaden");

    f.broaden();
    world.gop.fence();

    const auto after = count(f);
    if (world.rank() == 0)
        printf("  T8 after broaden:  %ld nodes with a computed dnorm_tree, "
               "%ld with a nonzero norm_tree, %ld with a computed sup_tree\n", after[0], after[1], after[2]);
    t.checkpoint(after[1] == 0, "T8 broaden clears norm_tree");
    t.checkpoint(after[0] == 0, "T8 broaden clears dnorm_tree");
    t.checkpoint(after[2] == 0, "T8 broaden clears sup_tree");
    return t.end();
}

/// T9: the vector norms accept an empty vector. norm2s()/norm2s_T() index
/// norms[0] to get the buffer address, which is undefined behaviour on an
/// empty container even with a length of zero.
int test_empty_vector_norms(World& world) {
    test_output t("mul_sparse T9: the vector norms on an empty vector");
    t.set_do_print(world.rank() == 0);

    const std::vector<Function<double,D>> empty;

    const std::vector<double> ns = norm2s(world, empty);
    t.checkpoint(ns.empty(), "T9 norm2s returns an empty vector");

    const Tensor<double> nt = norm2s_T(world, empty);
    t.checkpoint(nt.size() == 0, "T9 norm2s_T returns an empty tensor");

    t.checkpoint(norm2(world, empty), 0.0, 1.e-15, "T9 norm2 returns zero");
    return t.end();
}

/// T10: norm2() and trace() answer from a redundant or nonstandard-with-leaves
/// tree without converting it. Adding up every node overcounts there -- the
/// scaling coefficients sit on every level -- but the leaves alone are exactly
/// the reconstructed tree, so the answer costs no mutation and no extra fence.
/// mul_sparse() leaves its operands redundant and callers take norms of them in
/// a loop, so converting here would be paid back by the next multiplication.
int test_norms_preserve_tree_state(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T10: norm2/trace do not convert the tree");
    t.set_do_print(world.rank() == 0);

    Function<double,D> f = copy(pairs.front().f);
    f.reconstruct();
    const double exact_norm = f.norm2();
    const double exact_trace = f.trace();

    for (const TreeState state : {redundant, nonstandard_with_leaves}) {
        f.reconstruct();
        f.change_tree_state(state, true);
        t.checkpoint(f.get_impl()->get_tree_state() == state,
                     "T10 operand is in the state under test");
        // without this the internal nodes would not be there to overcount
        t.checkpoint(f.get_impl()->has_coefficients_on_leaves_only(),
                     "T10 the state keeps coefficients on the internal nodes too");

        const double n = f.norm2();
        const double tr = f.trace();
        if (world.rank() == 0)
            printf("  T10 %s: norm %.14e (exact %.14e)  trace %.14e (exact %.14e)\n",
                   (state == redundant ? "redundant" : "nonstandard_with_leaves"),
                   n, exact_norm, tr, exact_trace);
        t.checkpoint(std::abs(n - exact_norm) < 1.e-12 * exact_norm,
                     "T10 norm2 equals ||f||");
        t.checkpoint(std::abs(tr - exact_trace) < 1.e-12 * std::abs(exact_trace),
                     "T10 trace equals int(f)");
        t.checkpoint(f.get_impl()->get_tree_state() == state,
                     "T10 norm2/trace leave the tree state alone");

        // the vector forms take the same path
        const std::vector<Function<double,D>> v{f};
        const std::vector<double> ns = norm2s(world, v);
        t.checkpoint(std::abs(ns[0] - exact_norm) < 1.e-12 * exact_norm,
                     "T10 norm2s equals ||f||");
        t.checkpoint(std::abs(norm2(world, v) - exact_norm) < 1.e-12 * exact_norm,
                     "T10 vector norm2 equals ||f||");
        t.checkpoint(f.get_impl()->get_tree_state() == state,
                     "T10 the vector norms leave the tree state alone");
    }
    return t.end();
}

/// T11: compress() clears the leaf coefficients, and snorm has to go with them.
/// recur_down_for_contraction_map() reads snorm > 0 as "this node holds s
/// coefficients" and then skips building them, so a stale snorm on an emptied
/// node sends the recursion into unfilter() with an empty tensor. The norm the
/// parent filters with is norm_tree, which does keep the leaf norm -- reading it
/// after the coefficients are gone is what propagates a zero to the root.
int test_compress_clears_snorm(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T11: compress clears snorm with the coefficients");
    t.set_do_print(world.rank() == 0);

    Function<double,D> f = copy(pairs.front().f);
    f.reconstruct();
    const double exact = f.norm2();
    f.compress();

    long empty_nodes = 0, empty_with_snorm = 0;
    const auto& c = f.get_impl()->get_coeffs();
    for (auto it = c.begin(); it != c.end(); ++it) {
        if (it->second.has_coeff()) continue;
        ++empty_nodes;
        if (it->second.get_snorm() > 0.0) ++empty_with_snorm;
    }
    world.gop.sum(empty_nodes);
    world.gop.sum(empty_with_snorm);

    if (world.rank() == 0)
        printf("  T11 %ld nodes without coefficients, %ld of them with a positive snorm\n",
               empty_nodes, empty_with_snorm);
    // without this the test would pass on an empty premise
    t.checkpoint(empty_nodes > 0, "T11 compress leaves nodes without coefficients");
    t.checkpoint(empty_with_snorm == 0,
                 "T11 no coefficient-less node carries a positive snorm");

    // ... and the leaf norm still reached the root
    t.checkpoint(f.norm2(), exact, 1.e-12 * exact, "T11 the compressed norm is still ||f||");
    return t.end();
}

/// number of nodes of scaled whose norms (sup_tree included) are not |q| times
/// those of orig, or whose dnorm_tree or sup_tree is uncomputed; orig and scaled
/// share a tree
static long nodes_with_unscaled_norms(World& world, const Function<double,D>& orig,
                                      const Function<double,D>& scaled, double q) {
    auto agrees = [q](double s, double o) {
        return std::abs(s - std::abs(q)*o) <= 1.e-14 * std::abs(q)*o;
    };
    long bad = 0;
    const auto& oc = orig.get_impl()->get_coeffs();
    const auto& sc = scaled.get_impl()->get_coeffs();
    for (auto it = sc.begin(); it != sc.end(); ++it) {
        const auto& onode = oc.find(it->first).get()->second;
        const auto& snode = it->second;
        if (snode.get_dnorm_tree() >= NORM_TREE_UNCOMPUTED
            or snode.get_sup_tree() >= NORM_TREE_UNCOMPUTED
            or not agrees(snode.get_sup_tree(),   onode.get_sup_tree())
            or not agrees(snode.get_norm_tree(),  onode.get_norm_tree())
            or not agrees(snode.get_dnorm_tree(), onode.get_dnorm_tree())
            or not agrees(snode.get_snorm(),      onode.get_snorm())
            or not agrees(snode.get_dnorm(),      onode.get_dnorm())) ++bad;
    }
    world.gop.sum(bad);
    return bad;
}

/// number of nodes where exactly one of orig and scaled carries the
/// norm_tree == -1.0 broadened marker; orig and scaled share a tree
static long nodes_with_moved_marker(World& world, const Function<double,D>& orig,
                                    const Function<double,D>& scaled) {
    long moved = 0;
    const auto& oc = orig.get_impl()->get_coeffs();
    const auto& sc = scaled.get_impl()->get_coeffs();
    for (auto it = sc.begin(); it != sc.end(); ++it) {
        const auto& onode = oc.find(it->first).get()->second;
        const auto& snode = it->second;
        if ((snode.get_norm_tree() == -1.0) != (onode.get_norm_tree() == -1.0)) ++moved;
    }
    world.gop.sum(moved);
    return moved;
}

/// number of nodes carrying the norm_tree == -1.0 broadened marker
static long nodes_with_marker(World& world, const Function<double,D>& f) {
    long marked = 0;
    for (const auto& datum : f.get_impl()->get_coeffs()) {
        if (datum.second.get_norm_tree() == -1.0) ++marked;
    }
    world.gop.sum(marked);
    return marked;
}

/// ||mul_sparse(qf, g, tol) - ref||, or 1e300 if the kernel throws
static double sparse_product_error(World& world, const Function<double,D>& qf,
                                   const Function<double,D>& g,
                                   const Function<double,D>& ref, double tol,
                                   const std::string& tag) {
    double err = 1.e300;
    try {
        Function<double,D> prod = mul_sparse(qf, g, tol);
        err = (prod - ref).norm2();
    } catch (const MadnessException& e) {
        if (world.rank() == 0) printf("  %s mul_sparse threw: %s\n", tag.c_str(), e.what());
    }
    return err;
}

/// T12: q*f copies f's tree state, so a scaled copy of a redundant function is
/// redundant too, and mul_sparse() takes that label at its word and skips
/// make_redundant(). The copy therefore has to carry the norms the screen reads,
/// scaled by |q|; with them left uncomputed the kernel's redundant-form check
/// fires on the scaled operand.
int test_scaled_redundant_operand(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T12: a scaled redundant function keeps its norms");
    t.set_do_print(world.rank() == 0);

    const operand_pair& p = pairs.front();
    Function<double,D> f = copy(p.f);
    Function<double,D> g = copy(p.g);
    f.make_redundant(true);

    for (const double q : {-1.0, 2.5}) {
        Function<double,D> qf = q*f;
        t.checkpoint(qf.get_impl()->get_tree_state() == redundant,
                     "T12 q*f inherits the redundant state, q=" + std::to_string(q));

        const long bad = nodes_with_unscaled_norms(world, f, qf, q);
        if (world.rank() == 0)
            printf("  T12 q=%4.1f  %ld nodes with norms that are not |q| times f's\n", q, bad);
        t.checkpoint(bad == 0, "T12 q*f carries |q| times f's norms, q=" + std::to_string(q));

        // the scaled copy goes straight into the kernel
        const double err = sparse_product_error(world, qf, g, q*p.fg, 1.e-6, "T12");
        if (world.rank() == 0) printf("  T12 q=%4.1f  err %.3e\n", q, err);
        t.checkpoint(err, C_MAX * 1.e-6, "T12 mul_sparse(q*f, g), q=" + std::to_string(q));
    }
    return t.end();
}

/// T13: f.scale(q) scales in place and keeps the tree state, so the node norms
/// have to follow the coefficients. Left alone they stay computed but wrong by
/// |q|: the kernel's check passes, and the screen silently judges the product
/// against norms |q| times too small or too large.
int test_scale_inplace_norms(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T13: in-place scaling scales the norms");
    t.set_do_print(world.rank() == 0);

    const operand_pair& p = pairs.front();
    Function<double,D> f = copy(p.f);
    Function<double,D> g = copy(p.g);
    f.make_redundant(true);

    for (const double q : {-1.0, 2.5}) {
        // copy() keeps both the tree state and the norms of f
        Function<double,D> qf = copy(f);
        qf.scale(q);
        t.checkpoint(qf.get_impl()->get_tree_state() == redundant,
                     "T13 f.scale(q) keeps the redundant state, q=" + std::to_string(q));

        const long bad = nodes_with_unscaled_norms(world, f, qf, q);
        if (world.rank() == 0)
            printf("  T13 q=%4.1f  %ld nodes with norms that are not |q| times f's\n", q, bad);
        t.checkpoint(bad == 0, "T13 f.scale(q) scales the norms by |q|, q=" + std::to_string(q));

        const double err = sparse_product_error(world, qf, g, q*p.fg, 1.e-6, "T13");
        if (world.rank() == 0) printf("  T13 q=%4.1f  err %.3e\n", q, err);
        t.checkpoint(err, C_MAX * 1.e-6, "T13 mul_sparse(f.scale(q), g), q=" + std::to_string(q));
    }
    return t.end();
}

/// T14: broaden and refine_op mark handled nodes with norm_tree == -1.0,
/// and a later broaden skips those nodes by exact comparison. Scaling must
/// preserve the marker or broadening can repeat before norms are recomputed.
int test_scale_preserves_broaden_marker(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T14: scaling preserves the broadened marker");
    t.set_do_print(world.rank() == 0);

    Function<double,D> f = copy(pairs.front().f);
    f.broaden(FunctionDefaults<D>::get_bc(), false);
    world.gop.fence();

    const long marked = nodes_with_marker(world, f);
    if (world.rank() == 0) printf("  T14 %ld nodes carry the broadened marker\n", marked);
    t.checkpoint(marked > 0, "T14 broaden(false) leaves markers to preserve");

    for (const double q : {0.0, 2.5}) {
        Function<double,D> qf = q*f;
        long moved = nodes_with_moved_marker(world, f, qf);
        if (world.rank() == 0)
            printf("  T14 q=%4.1f  q*f moved %ld broadened markers\n", q, moved);
        t.checkpoint(moved == 0,
                     "T14 q*f preserves broadened markers, q=" + std::to_string(q));

        qf = copy(f);
        qf.scale(q);
        moved = nodes_with_moved_marker(world, f, qf);
        if (world.rank() == 0)
            printf("  T14 q=%4.1f  f.scale(q) moved %ld broadened markers\n", q, moved);
        t.checkpoint(moved == 0,
                     "T14 f.scale(q) preserves broadened markers, q=" + std::to_string(q));
    }
    return t.end();
}

/// T15: sup_tree, the bound on sup|f| per box that the Hoelder screen reads, comes out of the passes that compute
/// the norms (make_redundant, compress): computed on every node, and wherever a box's subtree is local equal, bit
/// for bit, to the bound FunctionImpl::sup_bound_map() builds from the leaves
int test_sup_tree(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T15: sup_tree from the passes against sup_bound_map");
    t.set_do_print(world.rank() == 0);
    // for hunting the 2-rank failure: MSP_T15_REPEAT repeats the whole check, MSP_T15_FENCE adds fences (bit 1: after the
    // copy, 2: after the reconstruct, 8: before the copy)
    const int reps = std::getenv("MSP_T15_REPEAT") ? std::atoi(std::getenv("MSP_T15_REPEAT")) : 1;
    const int fences = std::getenv("MSP_T15_FENCE") ? std::atoi(std::getenv("MSP_T15_FENCE")) : 0;
    for (int rep = 0; rep < reps; ++rep)
    for (auto& p : pairs) {
        for (const TreeState state : {redundant, compressed}) {
            if (fences & 8) world.gop.fence();
            Function<double,D> f = copy(p.f);
            if (fences & 1) world.gop.fence();
            f.reconstruct();
            if (fences & 2) world.gop.fence();
            const SupBoundMap<D> map = f.get_impl()->sup_bound_map();
            // what each local node looked like when the map was built, for the record of a difference below:
            // has_children, has_coeff, the sup_tree it carried, its norm_tree
            std::map<Key<D>, std::array<double, 4>> before;
            for (auto it = f.get_impl()->get_coeffs().begin(); it != f.get_impl()->get_coeffs().end(); ++it)
                before[it->first] = {double(it->second.has_children()), double(it->second.has_coeff()),
                                     it->second.get_sup_tree(), it->second.get_norm_tree()};
            // the map and the snapshot read each rank's own nodes, while the state change below is collective: once
            // one rank starts it, its tasks run on the other rank's nodes (compress clears the leaves' coefficients).
            // Without this fence a rank still reading saw leaves already cleared, about one check in thirty on 2 ranks.
            world.gop.fence();
            f.change_tree_state(state, true);
            // nodes, nodes without a computed sup_tree, nodes the map holds, and of those the ones that differ
            std::array<long, 4> n{0, 0, 0, 0};
            const auto& c = f.get_impl()->get_coeffs();
            for (auto it = c.begin(); it != c.end(); ++it) {
                const double st = it->second.get_sup_tree();
                const double m = map(it->first);
                ++n[0];
                if (not (st < NORM_TREE_UNCOMPUTED)) ++n[1];
                if (m == SupBoundMap<D>::unbounded()) continue;
                ++n[2];
                if (st != m) {
                    // the first few, for the record: which boxes, and by how much
                    if (++n[3] <= 4) {
                        std::ostringstream ks;
                        ks << it->first;
                        const auto b = before.find(it->first);
                        const bool had = (b != before.end());
                        printf("    T15 differs on rank %d: key %s %s, sup_tree %.17g, sup_bound_map %.17g; when the map "
                               "was built: %s, children %d, coeff %d, sup_tree %.17g, norm_tree %.6g; now norm_tree "
                               "%.6g\n", world.rank(), ks.str().c_str(), it->second.has_children() ? "interior" : "leaf",
                               st, m, had ? "present" : "absent", had ? int(b->second[0]) : -1,
                               had ? int(b->second[1]) : -1, had ? b->second[2] : -1.0, had ? b->second[3] : -1.0,
                               it->second.get_norm_tree());
                    }
                }
            }
            world.gop.sum(n.data(), n.size());
            const std::string tag = "T15 " + p.name + (state == redundant ? " redundant" : " compressed");
            if (world.rank() == 0)
                printf("  %s: %ld nodes, %ld without sup_tree, %ld held by sup_bound_map, %ld differ\n", tag.c_str(),
                       n[0], n[1], n[2], n[3]);
            if (n[3] > 0) {
                // after the fact, so it cannot hide the failure: are the operand's and the copy's trees consistent
                // (every node's parent present and marked as having children)?
                world.gop.fence();
                const bool pf = p.f.get_impl()->verify_parents_and_children();
                const bool fc = f.get_impl()->verify_parents_and_children();
                printf("    T15 after the failure on rank %d: operand tree %s (state %d), copy tree %s\n", world.rank(),
                       pf ? "consistent" : "INCONSISTENT", int(p.f.get_impl()->get_tree_state()),
                       fc ? "consistent" : "INCONSISTENT");
                world.gop.fence();
            }
            t.checkpoint(n[0] > 0 and n[2] > 0, tag + " has nodes to compare");
            t.checkpoint(n[1] == 0, tag + " sup_tree computed on every node");
            t.checkpoint(n[3] == 0, tag + " sup_tree equals sup_bound_map");
        }
    }
    return t.end();
}

/// T16: every way a node is made or copied keeps what it is given. The short constructors, and the explicit-norm one
/// without a bound, leave sup_tree uncomputed; the explicit-norm one with a bound and every copy (copy constructor,
/// assignment, convert, serialization) carry dnorm_tree and sup_tree. July's shuffle bug was a copy path that lost
/// dnorm_tree, and the screen then stopped early on the operands it reached.
int test_node_paths(World& world) {
    test_output t("mul_sparse T16: node constructors and copies keep dnorm_tree and sup_tree");
    t.set_do_print(world.rank() == 0);
    typedef FunctionNode<double,D> nodeT;
    const double U = NORM_TREE_UNCOMPUTED;
    Tensor<double> tc(8L, 8L, 8L);
    tc.fillrandom();
    const nodeT::coeffT c(tc);
    auto check = [&t](const auto& node, const double dnt, const double st, const std::string& what) {
        t.checkpoint(node.get_dnorm_tree() == dnt and node.get_sup_tree() == st, "T16 " + what);
    };
    check(nodeT(), U, U, "default constructor: uncomputed");
    check(nodeT(c, true), U, U, "constructor from coefficients: uncomputed");
    check(nodeT(c, 1.0, true), U, U, "constructor with norm_tree: uncomputed");
    check(nodeT(c, 1.0, 0.5, 0.9, 0.4, true), 0.5, U, "explicit norms without a bound: sup_tree uncomputed");
    const nodeT full(c, 1.0, 0.5, 0.9, 0.4, true, 2.0);
    check(full, 0.5, 2.0, "explicit norms with a bound");
    check(nodeT(full), 0.5, 2.0, "copy constructor");
    nodeT assigned;
    assigned = full;
    check(assigned, 0.5, 2.0, "assignment");
    check(full.convert<double_complex>(), 0.5, 2.0, "convert to complex");
    std::vector<unsigned char> buf;
    {
        archive::VectorOutputArchive oar(buf);
        nodeT out(full);
        oar & out;
    }
    nodeT in;
    archive::VectorInputArchive iar(buf);
    iar & in;
    check(in, 0.5, 2.0, "serialization round trip");
    return t.end();
}

/// T17: the archive carries sup_tree, and an archive written before sup_tree was part of the node format is refused
/// by its header magic instead of being misread
int test_archive(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T17: save/load keeps sup_tree, the previous format is refused");
    t.set_do_print(world.rank() == 0);
    const char* name = "test_mul_sparse_T17";
    Function<double,D> f = copy(pairs.back().f);
    f.make_redundant(true);
    {
        archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> out(world, name, 1);
        out & f;
    }
    Function<double,D> g;
    {
        archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> in(world, name, 1);
        in & g;
        in.close();
        in.remove();
    }
    // g has f's distribution (both the default process map): every node of g has f's sup_tree at the same box
    long nodes = 0, differ = 0;
    const auto& fc = f.get_impl()->get_coeffs();
    for (auto it = g.get_impl()->get_coeffs().begin(); it != g.get_impl()->get_coeffs().end(); ++it) {
        ++nodes;
        if (not fc.is_local(it->first)) {
            ++differ;
            continue;
        }
        const auto fit = fc.find(it->first).get();
        if (fit == fc.end() or fit->second.get_sup_tree() != it->second.get_sup_tree()) ++differ;
    }
    world.gop.sum(nodes);
    world.gop.sum(differ);
    if (world.rank() == 0) printf("  T17 %ld nodes loaded, %ld with a sup_tree that is not f's\n", nodes, differ);
    t.checkpoint(get_tree_state(g) == redundant, "T17 the loaded function is redundant");
    t.checkpoint(nodes > 0 and differ == 0, "T17 the loaded function carries f's sup_tree on every node");

    // the header of the previous format, magic 7776770 (FunctionNode without _sup_tree)
    {
        archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> out(world, name, 1);
        out & long(7776770) & long(TensorTypeData<double>::id) & long(D) & long(FunctionDefaults<D>::get_k())
            & FunctionDefaults<D>::get_cell();
    }
    bool refused = false;
    {
        archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> in(world, name, 1);
        try {
            Function<double,D> h;
            in & h;
        } catch (const MadnessException&) {
            refused = true;
        }
        in.close();
        in.remove();
    }
    t.checkpoint(refused, "T17 an archive with the previous magic is refused");
    return t.end();
}

/// largest |value| of the polynomial with scaling coefficients c on a box at `level`, on a 24^3 grid with the faces
/// and corners (the leaf bracket of T18; test_supbound_quantify refines it)
static double grid_max(const Tensor<double>& c, const int level) {
    const long k = c.dim(0);
    Tensor<double> phi(k, 24L);
    std::vector<double> v(k);
    for (long j = 0; j < 24; ++j) {
        legendre_scaling_functions(double(j) / 23.0, k, v.data());
        for (long i = 0; i < k; ++i) phi(i, j) = v[i];
    }
    return transform(c, phi).absmax() * std::pow(2.0, 0.5 * D * level) / std::sqrt(FunctionDefaults<D>::get_cell_volume());
}

/// median of v (0 if empty)
static double median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

/// the leaf bracket of f, per level: M (largest value at the quadrature points), S (the sum bound) and T (largest
/// value on the 24^3 grid), as medians of M/T and S/T; printed, not asserted
static void print_leaf_bracket(World& world, const std::string& name, const Function<double,D>& f) {
    const auto* impl = f.get_impl().get();
    std::vector<double> lev, mt, st;
    for (auto it = impl->get_coeffs().begin(); it != impl->get_coeffs().end(); ++it) {
        if (it->second.has_children() or not it->second.has_coeff()) continue;
        const Tensor<double> c = it->second.coeff().full_tensor_copy();
        const int n = it->first.level();
        const double T = std::max(grid_max(c, n), impl->coeffs2values(it->first, c).absmax());
        if (T <= FunctionDefaults<D>::get_thresh()) continue;
        lev.push_back(n);
        mt.push_back(impl->coeffs2values(it->first, c).absmax() / T);
        st.push_back(impl->coeff_sup_bound(it->first, c) / T);
    }
    lev = world.gop.concat0(lev);
    mt = world.gop.concat0(mt);
    st = world.gop.concat0(st);
    if (world.rank() != 0) return;
    printf("  T18 leaf bracket of %s (leaves with T > thresh): level, leaves, M/T median, S/T median\n", name.c_str());
    for (int l = 0; l < 30; ++l) {
        std::vector<double> m, s;
        for (std::size_t i = 0; i < lev.size(); ++i)
            if (int(lev[i]) == l) {
                m.push_back(mt[i]);
                s.push_back(st[i]);
            }
        if (not m.empty()) printf("    %2d %5zu  %.3f  %.3f\n", l, m.size(), median(m), median(s));
    }
}

/// ||s_f s_g - q_B||_B: the collocation error of do_mul at B, the product of the operands' polynomials there against
/// the product the kernel formed from them, by Gauss quadrature with 2k points per dimension (exact: the integrand is
/// a polynomial of degree 4k - 4 per dimension)
static double collocation_error(const Key<D>& B, const Tensor<double>& cf, const Tensor<double>& cg,
                                const Tensor<double>& cq) {
    const long k = cf.dim(0);
    const int npt = 2 * k;
    std::vector<double> x(npt), w(npt), v(k);
    gauss_legendre(npt, 0.0, 1.0, x.data(), w.data());
    Tensor<double> phi(k, long(npt));
    for (int j = 0; j < npt; ++j) {
        legendre_scaling_functions(x[j], k, v.data());
        for (long i = 0; i < k; ++i) phi(i, j) = v[i];
    }
    const double scale = std::pow(2.0, 0.5 * D * B.level()) / std::sqrt(FunctionDefaults<D>::get_cell_volume());
    Tensor<double> sfg = transform(cf, phi);
    sfg.emul(transform(cg, phi));
    const Tensor<double> q = transform(cq, phi);
    double c2 = 0.0;
    for (int a = 0; a < npt; ++a)
        for (int b = 0; b < npt; ++b)
            for (int e = 0; e < npt; ++e) {
                const double d = scale * scale * sfg(a, b, e) - scale * q(a, b, e);
                c2 += w[a] * w[b] * w[e] * d * d;
            }
    return std::sqrt(c2 * FunctionDefaults<D>::get_cell_volume() * std::pow(0.5, double(D) * B.level()));
}

/// T18 (the hand-off's T7): where the Hoelder screen forms the product above the union leaves, at a box B, the error it
/// commits by stopping there is bounded by its estimate S(c_f(B)) ||d_g|| + ||d_f|| sup_tree(g)(B). Against fg, the
/// independent projection of the analytic product, the product's error on B also holds the collocation error of
/// do_mul at B, C(B) (outside every estimate, measured here), and the floor of the representation (EXACT_FLOOR):
/// ||q - fg||_B <= estimate + C(B) + EXACT_FLOOR at every such box. Printed per level: the boxes, the median and the
/// minimum of estimate / (||q - fg||_B - C(B)) where that is positive, and the median share of C(B) in the error;
/// the operands' leaf bracket once. Boxes whose operand coefficients live on another rank are skipped (counted).
int test_hoelder_per_box(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T18: the Hoelder estimate bounds the error where the product is formed");
    t.set_do_print(world.rank() == 0);
    for (auto& p : pairs) print_leaf_bracket(world, p.name + " f", p.f);
    const MulScreen mode = mul_screen();
    mul_screen() = MulScreen::hoelder;
    for (auto& p : pairs) {
        for (const double tol : {1.e-4, 1.e-6}) {
            Function<double,D> q = mul_sparse(p.f, p.g, tol);
            Function<double,D> q0 = mul_sparse(p.f, p.g, 0.0);
            // the screened product's error on every box: norm_tree of the redundant difference
            Function<double,D> d = q - p.fg;
            d.make_redundant(true);
            q.reconstruct();    // the subtraction compressed it: the formed products back on its leaves
            const auto* fi = p.f.get_impl().get();
            const auto* gi = p.g.get_impl().get();
            // the operand's coefficients, dnorm_tree and sup_tree at B: its node's, or below its leaf the leaf's
            // polynomial restricted to B (no detail); false if the coefficients are on another rank
            auto at = [](const FunctionImpl<double,D>* impl, const Key<D>& B, Tensor<double>& c, double& dn, double& sup) {
                const auto& coeffs = impl->get_coeffs();
                for (Key<D> a = B;; a = a.parent()) {
                    if (not coeffs.is_local(a)) return false;
                    const auto it = coeffs.find(a).get();
                    if (it != coeffs.end()) {
                        if (a == B) {
                            c = it->second.coeff().full_tensor_copy();
                            dn = it->second.get_dnorm_tree();
                            sup = it->second.get_sup_tree();
                        } else {
                            c = impl->parent_to_child(it->second.coeff(), a, B).full_tensor_copy();
                            dn = 0.0;
                            sup = impl->coeff_sup_bound(B, c);
                        }
                        return true;
                    }
                    if (a.level() == 0) return false;
                }
            };
            // per box above the union leaves: level, estimate / error; boxes checked, beyond the floor, skipped
            std::vector<double> lev, ratio, share;
            long checked = 0, over = 0, skipped = 0;
            const auto& qc = q.get_impl()->get_coeffs();
            const auto& q0c = q0.get_impl()->get_coeffs();
            const auto& dc = d.get_impl()->get_coeffs();
            for (auto it = qc.begin(); it != qc.end(); ++it) {
                const Key<D>& B = it->first;
                if (it->second.has_children()) continue;
                const auto u = q0c.find(B).get();
                if (u != q0c.end() and not u->second.has_children()) continue;     // a union leaf: formed as at tol 0
                Tensor<double> cf, cg;
                double dnf = 0, supf = 0, dng = 0, supg = 0;
                if (not at(fi, B, cf, dnf, supf) or not at(gi, B, cg, dng, supg)) {
                    ++skipped;
                    continue;
                }
                const double est = fi->coeff_sup_bound(B, cf) * dng + dnf * supg;
                const auto dit = dc.find(B).get();
                MADNESS_CHECK(dit != dc.end());
                const double err = dit->second.get_norm_tree();
                const double coll = collocation_error(B, cf, cg, it->second.coeff().full_tensor_copy());
                ++checked;
                if (err > est + coll + EXACT_FLOOR) ++over;
                if (err > 0.0) {
                    lev.push_back(B.level());
                    ratio.push_back(err > coll ? est / (err - coll) : -1.0);     // -1: all of it collocation
                    share.push_back(std::min(1.0, coll / err));
                }
            }
            world.gop.sum(checked);
            world.gop.sum(over);
            world.gop.sum(skipped);
            lev = world.gop.concat0(lev);
            ratio = world.gop.concat0(ratio);
            share = world.gop.concat0(share);
            const std::string tag = "T18 " + p.name + " tol=" + fmt_tol(tol);
            if (world.rank() == 0) {
                printf("  %s: %ld boxes above the union leaves checked, %ld beyond estimate + C + EXACT_FLOOR, %ld on "
                       "another rank; per level: boxes, those whose error is all collocation, C/error median, and over the "
                       "rest estimate/(error - C) median and min:", tag.c_str(), checked, over, skipped);
                for (int l = 0; l < 30; ++l) {
                    std::vector<double> r, c;
                    long n = 0;
                    for (std::size_t i = 0; i < lev.size(); ++i)
                        if (int(lev[i]) == l) {
                            ++n;
                            c.push_back(share[i]);
                            if (ratio[i] >= 0.0) r.push_back(ratio[i]);
                        }
                    if (n == 0) continue;
                    printf(" %d: %ld %ld %.2f", l, n, n - long(r.size()), median(c));
                    if (not r.empty()) printf(" %.2g %.2g", median(r), *std::min_element(r.begin(), r.end()));
                    printf(";");
                }
                printf("\n");
            }
            t.checkpoint(over == 0, tag + " the error is within the estimate, the collocation and the floor");
        }
    }
    mul_screen() = mode;
    return t.end();
}

/// T19: tol -> 0 loses nothing: under either criterion the product at a tiny tolerance is the tol = 0 product to the
/// tolerance (July 2026 lost 3e-2 Eh in one iteration at an effective 1e-11 on a development branch, never diagnosed)
int test_tol_to_zero(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T19: a tolerance going to zero gives the exact product");
    t.set_do_print(world.rank() == 0);
    const MulScreen mode = mul_screen();
    for (const MulScreen m : {MulScreen::cs, MulScreen::hoelder}) {
        mul_screen() = m;
        for (auto& p : pairs) {
            const Function<double,D> q0 = mul_sparse(p.f, p.g, 0.0);
            for (const double tol : {1.e-10, 1.e-12, 1.e-14}) {
                const double err = (mul_sparse(p.f, p.g, tol) - q0).norm2();
                if (world.rank() == 0)
                    printf("  T19%s %-13s tol %.0e  ||q(tol) - q(0)|| %.3e\n", screen_tag().c_str(), p.name.c_str(), tol, err);
                t.checkpoint(err <= C_TOL0 * tol, "T19" + screen_tag() + " " + p.name + " tol=" + fmt_tol(tol));
            }
        }
    }
    mul_screen() = mode;
    return t.end();
}

/// T20: the identity pair. Without MAD_MUL_SCREEN the criterion is today's; and at tol = 0 the two criteria form the
/// same product, node for node and bit for bit (both stop only where neither operand has detail left)
int test_identity(World& world, std::vector<operand_pair>& pairs, const MulScreen initial) {
    test_output t("mul_sparse T20: the default criterion, and tol = 0 under both");
    t.set_do_print(world.rank() == 0);
    if (std::getenv("MAD_MUL_SCREEN") == nullptr)
        t.checkpoint(initial == MulScreen::cs, "T20 without MAD_MUL_SCREEN the criterion is today's");
    const MulScreen mode = mul_screen();
    for (auto& p : pairs) {
        mul_screen() = MulScreen::cs;
        const Function<double,D> a = mul_sparse(p.f, p.g, 0.0);
        mul_screen() = MulScreen::hoelder;
        const Function<double,D> b = mul_sparse(p.f, p.g, 0.0);
        long nodes = 0, differ = 0;
        const auto& bc = b.get_impl()->get_coeffs();
        for (auto it = a.get_impl()->get_coeffs().begin(); it != a.get_impl()->get_coeffs().end(); ++it) {
            ++nodes;
            const auto jt = bc.find(it->first).get();
            if (jt == bc.end() or jt->second.has_children() != it->second.has_children()
                or jt->second.has_coeff() != it->second.has_coeff()
                or (it->second.has_coeff()
                    and (jt->second.coeff().full_tensor_copy() - it->second.coeff().full_tensor_copy()).absmax() != 0.0))
                ++differ;
        }
        long bnodes = bc.size();
        world.gop.sum(nodes);
        world.gop.sum(differ);
        world.gop.sum(bnodes);
        if (world.rank() == 0)
            printf("  T20 %-13s tol=0: %ld nodes under cs, %ld under hoelder, %ld differ\n", p.name.c_str(), nodes,
                   bnodes, differ);
        t.checkpoint(differ == 0 and nodes == bnodes, "T20 " + p.name + " tol=0 is the same product under both criteria");
    }
    mul_screen() = mode;
    return t.end();
}

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    startup(world, argc, argv, true);

    FunctionDefaults<D>::set_k(8);
    FunctionDefaults<D>::set_thresh(THRESH_REF);
    FunctionDefaults<D>::set_cubic_cell(-8.0, 8.0);
    FunctionDefaults<D>::set_refine(true);
    FunctionDefaults<D>::set_initial_level(2);
    FunctionDefaults<D>::set_autorefine(false);

    // the criterion the environment asks for (MAD_MUL_SCREEN), restored after the tests that set it
    const MulScreen initial = mul_screen();
    // for hunting order-dependent failures: MSP_TESTS selects tests by number ("1-6,15"; all by default), MSP_FENCE=1
    // fences before each test
    std::vector<char> run(21, 1);
    if (const char* sel = std::getenv("MSP_TESTS")) {
        std::fill(run.begin(), run.end(), 0);
        std::string list(sel);
        for (std::size_t q = 0; q < list.size();) {
            const std::size_t e = std::min(list.find(',', q), list.size());
            const std::string item = list.substr(q, e - q);
            const std::size_t dash = item.find('-');
            const int lo = std::stoi(item.substr(0, dash));
            const int hi = (dash == std::string::npos) ? lo : std::stoi(item.substr(dash + 1));
            for (int n = lo; n <= hi and n <= 20; ++n) run[n] = 1;
            q = e + 1;
        }
    }
    const bool fence_each = (std::getenv("MSP_FENCE") != nullptr);
    auto want = [&](const int n) {
        if (fence_each) world.gop.fence();
        return bool(run[n]);
    };
    int success = 0;
    {
        std::vector<operand_pair> pairs = make_pairs(world);
        for (const MulScreen mode : {MulScreen::cs, MulScreen::hoelder}) {
            mul_screen() = mode;
            if (want(1)) success += test_exact(world, pairs);
            if (want(2)) success += test_tree_state(world, pairs);
            if (want(3)) success += test_screening_accuracy(world, pairs);
            if (want(4)) success += test_self_consistency(world, pairs);
            if (want(5)) success += test_dot_sparse(world, pairs);
            if (want(6)) success += test_operand_state(world, pairs);
        }
        mul_screen() = initial;
        if (want(7)) success += test_norm_tree_states(world, pairs);
        if (want(8)) success += test_broaden_resets_norms(world, pairs);
        if (want(9)) success += test_empty_vector_norms(world);
        if (want(10)) success += test_norms_preserve_tree_state(world, pairs);
        if (want(11)) success += test_compress_clears_snorm(world, pairs);
        if (want(12)) success += test_scaled_redundant_operand(world, pairs);
        if (want(13)) success += test_scale_inplace_norms(world, pairs);
        if (want(14)) success += test_scale_preserves_broaden_marker(world, pairs);
        if (want(15)) success += test_sup_tree(world, pairs);
        if (want(16)) success += test_node_paths(world);
        if (want(17)) success += test_archive(world, pairs);
        if (want(18)) success += test_hoelder_per_box(world, pairs);
        if (want(19)) success += test_tol_to_zero(world, pairs);
        if (want(20)) success += test_identity(world, pairs, initial);
    }

    world.gop.fence();
    finalize();
    return success;
}
