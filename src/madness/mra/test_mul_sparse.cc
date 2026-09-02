/// Unit tests for the screened multiplication mul_sparse.
///
/// The operand pairs have closed-form products, so the reference is an
/// independent projection of the analytic product rather than another call to
/// the kernel under test.

#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/world/test_utilities.h>

#include <cmath>
#include <cstdio>
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
    test_output t("mul_sparse T1: tol=0 vs an independent projection");
    t.set_do_print(world.rank() == 0);
    for (auto& p : pairs) {
        Function<double,D> q = mul_sparse(p.f, p.g, 0.0);
        const double err = (q - p.fg).norm2();
        if (world.rank() == 0)
            printf("  T1 %-13s err %.3e   bound %.1e\n", p.name.c_str(), err, EXACT_FLOOR);
        t.checkpoint(err, EXACT_FLOOR, "T1 " + p.name);
    }
    return t.end();
}

/// T2: the result is a reconstructed tree, both scalar and vector entry points
int test_tree_state(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T2: output tree state");
    t.set_do_print(world.rank() == 0);
    operand_pair& p = pairs.front();

    for (double tol : {0.0, 1.e-5}) {
        Function<double,D> q = mul_sparse(p.f, p.g, tol);
        world.gop.fence();
        const long bad = interior_nodes_with_coeffs(world, q);
        t.checkpoint(get_tree_state(q) == reconstructed,
                     "T2 scalar state is reconstructed, tol=" + fmt_tol(tol));
        t.checkpoint(bad == 0L,
                     "T2 scalar coefficients on leaves only, tol=" + fmt_tol(tol));
        t.checkpoint(q.get_impl()->verify_tree_state_local(),
                     "T2 scalar tree state is self-consistent, tol=" + fmt_tol(tol));
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
    t.checkpoint(state_ok, "T2 vector state is reconstructed");
    t.checkpoint(layout_ok, "T2 vector coefficients on leaves only");
    return t.end();
}

/// T3: the screening error over a tol ladder, against the analytic product.
///
/// The bounds are measured, not derived: err ~ tol^p with p short of 1, and tree
/// quantization (adjacent tol yielding the same tree) leaves err flat and inflates
/// err/tol by up to a decade. Hence a loose C_MAX plus two tighter statements --
/// err falls monotonically with tol, and its decay exponent does not collapse.
int test_screening_accuracy(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T3: screening error vs tol");
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
                printf("  T3 %-13s tol %.1e  err %.3e  C %6.2f\n",
                       p.name.c_str(), tol, err, err / tol);
        }
        // 5% slack absorbs the quantization plateaus
        bool monotone = true;
        for (std::size_t i = 1; i < errs.size(); ++i)
            monotone = monotone and (errs[i] <= 1.05 * errs[i-1]);

        const double decades = std::log10(tols.front() / tols.back());
        const double pfit = std::log10(errs.front() / errs.back()) / decades;
        if (world.rank() == 0)
            printf("  T3 %-13s decay exponent p %.2f\n", p.name.c_str(), pfit);

        t.checkpoint(bounded,  "T3 " + p.name + " err < C*tol");
        t.checkpoint(monotone, "T3 " + p.name + " err falls with tol");
        t.checkpoint(pfit > P_MIN, "T3 " + p.name + " decay exponent");
    }
    return t.end();
}

/// T4: screened and exact multiplication agree to within tol, no analytic
/// reference needed
int test_self_consistency(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T4: screened vs exact");
    t.set_do_print(world.rank() == 0);

    for (auto& p : pairs) {
        Function<double,D> exact = mul_sparse(p.f, p.g, 0.0);
        for (double tol : {1.e-4, 1.e-6, 1.e-8}) {
            Function<double,D> q = mul_sparse(p.f, p.g, tol);
            const double err = (q - exact).norm2();
            const double C = err / tol;
            if (world.rank() == 0)
                printf("  T4 %-13s tol %.1e  err %.3e  C %.2f\n",
                       p.name.c_str(), tol, err, C);
            t.checkpoint(C < C_MAX, "T4 " + p.name + " tol=" + fmt_tol(tol));
        }
    }
    return t.end();
}

/// T5: dot routes through the screened kernel; check it against the analytic sum
int test_dot(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T5: dot vs the analytic sum of products");
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
        Function<double,D> q = dot(world, a, b, true, true, tol);
        const double err = (q - ref).norm2();
        const double bound = (tol == 0.0) ? EXACT_FLOOR : C_MAX * tol;
        if (world.rank() == 0)
            printf("  T5 dot tol %.1e  err %.3e  bound %.1e\n", tol, err, bound);
        t.checkpoint(err, bound, "T5 dot tol=" + fmt_tol(tol));
    }
    return t.end();
}

/// T6: mul_sparse leaves its operands redundant. The global norm2()/trace()
/// entry points must stay correct across that state.
int test_operand_state(World& world, std::vector<operand_pair>& pairs) {
    test_output t("mul_sparse T6: operand state and the norms read off it");
    t.set_do_print(world.rank() == 0);
    operand_pair& p = pairs.front();

    // reference values, taken while the operands are still reconstructed
    p.f.reconstruct();
    p.g.reconstruct();
    const double norm_before    = p.f.norm2();
    const double trace_before   = p.f.trace();
    const double g_trace_before = p.g.trace();

    Function<double,D> q = mul_sparse(p.f, p.g, 1.e-6);
    world.gop.fence();

    t.checkpoint(get_tree_state(p.f) == redundant,
                 "T6 left operand is left redundant");
    t.checkpoint(get_tree_state(p.g) == redundant,
                 "T6 right operand is left redundant");

    // each entry point has to be called on a tree that is still redundant at
    // the point of the call, and norm2()/trace() reconstruct as a side effect
    // -- so they cannot be measured one after the other on the same operand
    const double norm_after  = p.f.norm2();     // p.f: redundant -> reconstructed
    const double trace_after = p.g.trace();     // p.g: redundant -> reconstructed
    if (world.rank() == 0)
        printf("  T6 norm2  %.10e -> %.10e\n  T6 trace  %.10e -> %.10e\n",
               norm_before, norm_after, g_trace_before, trace_after);

    t.checkpoint(norm_after,  norm_before,    1.e-10, "T6 norm2 on a redundant function");
    t.checkpoint(trace_after, g_trace_before, 1.e-10, "T6 trace on a redundant function");
    t.checkpoint(get_tree_state(p.f) == reconstructed and
                 get_tree_state(p.g) == reconstructed,
                 "T6 norm2/trace leave the operand reconstructed");

    // the vector overload already reconstructs; it must still agree.  p.f is
    // reconstructed by now, so put it back into redundant form first
    p.f.make_redundant(true);
    t.checkpoint(get_tree_state(p.f) == redundant,
                 "T6 norm2s operand is redundant on entry");
    std::vector<Function<double,D>> v = {p.f};
    t.checkpoint(norm2s(world, v)[0], norm_before, 1.e-10,
                 "T6 vector norm2s agrees");

    // p.f has now been through reconstructed -> redundant -> reconstructed;
    // its trace must still be the value measured at the top
    t.checkpoint(p.f.trace(), trace_before, 1.e-10,
                 "T6 trace round-trips through redundant form");
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

    t.checkpoint(std::abs(nt_redundant  - exact) < 1.e-10,
                 "T7 root norm_tree on a redundant tree equals ||f||");
    t.checkpoint(std::abs(nt_compressed - exact) < 1.e-10,
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
/// mulXXveca() and multiplies at a level that is too coarse.
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

    // counts local nodes, summed over ranks: (computed dnorm_tree, nonzero
    // norm_tree)
    auto count = [&world](const Function<double,D>& g) {
        long computed = 0, nonzero = 0;
        const auto& c = g.get_impl()->get_coeffs();
        for (auto it = c.begin(); it != c.end(); ++it) {
            if (it->second.get_dnorm_tree() != NORM_TREE_UNCOMPUTED) ++computed;
            if (it->second.get_norm_tree() != 0.0) ++nonzero;
        }
        world.gop.sum(computed);
        world.gop.sum(nonzero);
        return std::make_pair(computed, nonzero);
    };

    const auto before = count(f);
    if (world.rank() == 0)
        printf("  T8 before broaden: %ld nodes with a computed dnorm_tree, "
               "%ld with a nonzero norm_tree\n", before.first, before.second);
    // without this the test would pass on an empty premise
    t.checkpoint(before.first > 0 and before.second > 0,
                 "T8 the tree carries computed norms into broaden");

    f.broaden();
    world.gop.fence();

    const auto after = count(f);
    if (world.rank() == 0)
        printf("  T8 after broaden:  %ld nodes with a computed dnorm_tree, "
               "%ld with a nonzero norm_tree\n", after.first, after.second);
    t.checkpoint(after.second == 0, "T8 broaden clears norm_tree");
    t.checkpoint(after.first == 0,  "T8 broaden clears dnorm_tree");
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

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    startup(world, argc, argv, true);

    FunctionDefaults<D>::set_k(8);
    FunctionDefaults<D>::set_thresh(THRESH_REF);
    FunctionDefaults<D>::set_cubic_cell(-8.0, 8.0);
    FunctionDefaults<D>::set_refine(true);
    FunctionDefaults<D>::set_initial_level(2);
    FunctionDefaults<D>::set_autorefine(false);

    int success = 0;
    {
        std::vector<operand_pair> pairs = make_pairs(world);
        success += test_exact(world, pairs);
        success += test_tree_state(world, pairs);
        success += test_screening_accuracy(world, pairs);
        success += test_self_consistency(world, pairs);
        success += test_dot(world, pairs);
        success += test_operand_state(world, pairs);
        success += test_norm_tree_states(world, pairs);
        success += test_broaden_resets_norms(world, pairs);
        success += test_empty_vector_norms(world);
    }

    world.gop.fence();
    finalize();
    return success;
}
