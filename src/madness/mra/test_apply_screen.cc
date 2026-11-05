/// Unit tests for the destination screen of the nonstandard apply (applyscreen.h).
///
/// The source is a Gaussian whose Coulomb potential is known in closed form, so the reference is an
/// independent projection of the analytic potential rather than another call to the kernel under test.
/// The screen is installed through the vector apply, apply(world, op, f, screens), so the skip is live;
/// with MAD_APPLY_SCREEN=0 set nothing is screened and the test refuses to run.
///
/// T1  a screen that never screens (an unbounded multiplier) reproduces the unscreened apply; its count of
///     the blocks that pass the source test is the unscreened shell loops' reach, which T2-T4 compare with.
/// T2  the error weighted by the multiplier stays within twice the unscreened one or C_T2 thresh, while
///     blocks are skipped and the result tree shrinks; the same for redundant multipliers. Measured (k 8,
///     thresh 1e-6, one rank): weighted errors unscreened m1 1.8e-7, m2 1.5e-7; screened, as C = error /
///     thresh, m1 0.18 and m2 0.15 with 1744 of 2417 blocks skipped.
/// T3  the shell loop (k 6, thresh 1e-4): a multiplier that lives only in a box B four bohr away from the
///     source, negligible everywhere between them. A block skipped by the screen must still count as used
///     for the shell loop, or a source box whose near shells lie outside B ends its loop before it reaches
///     B. Every block that reaches B is computed as without the screen, so the potential in B is the
///     unscreened one to roundoff (measured relative 1e-15) and the loops examine as many blocks as the
///     unscreened ones (2129, also checked in T2 and T4). With the skipped blocks not counted (the nused++
///     of the skip removed) the loops stop early (1445 blocks) and the potential in B is off by 5e-8
///     relative; T2 and T4 then fail the reach check too.
/// T3b the same box with a small but not negligible multiplier outside it (1e-3, 1e-2): kept blocks there
///     are computed but not accumulated, and must still count for the shell loop by their own size.
///     Counting them only when accumulated (their norm times the multiplier's size) stops the loops early
///     (1625 and 1921 blocks) and puts the potential in B off by 5e-8 relative at 1e-3.
/// T4  a complex function under a real screen.
/// T5  the sup bound per box (supbound.h, FunctionImpl::sup_bound_map) of T2's multipliers and of a cusp
///     exp(-35 r) refined to level 15 at the origin: on every leaf it lies above the function's largest
///     magnitude at the quadrature points and at the box corners and 50 random points (evaluated from the
///     leaf's coefficients, to roundoff), and above the analytic function less SLACK_T5 thresh (the
///     projection's pointwise error, which reaches 1.04 thresh at two points of m2); an interior box holds
///     the max over its children, a box below a leaf the leaf's bound. Prints how far the bound lies above
///     the leaf's rms magnitude and above its largest value at the quadrature points. Measured on the leaves
///     with max|f| > thresh: S / max|values| median 1.2-1.4, max 2.1-4.3; S / rms median 3-10, 15-150 at
///     the levels 3-5 of the Gaussians, 1.2-3 at the cusp's levels 10-13.
/// T6  the exchange's geometry: the source is the multiplier squared (the diagonal product bra_i vf_i), so
///     the multiplier is O(1) on top of the source at every level the source has. The multiplier-weighted
///     error against the unscreened potential stays below C_T6 thresh. Measured, C and the share of blocks
///     skipped: a Gaussian exp(-1.5 r^2), whose tree stops at level 5, 0.02 (75%); a cusp exp(-35 r),
///     refined to level 13, 3e-14 (58%): the deep levels, where a box is far smaller than one bohr^3, are
///     never skipped.

#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/mra/operator.h>
#include <madness/mra/legendre.h>
#include <madness/world/test_utilities.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace madness;

static const std::size_t D = 3;

/// projection threshold of the analytic reference potential
static const double THRESH_REF = 1.e-9;

/// T2, T4: the weighted error of the screened potential stays below C_T2 * thresh (or twice the unscreened
/// error); measured C 0.18 and 0.15, see the header
static const double C_T2 = 1.0;

/// T3, T3b: inside B the screened potential is the unscreened one to roundoff, relative to |chi u|;
/// measured at most 1.4e-15 (shell loops stopped too early: 5e-8), see the header
static const double REL_T3 = 1.e-10;

/// T6: the weighted error of the screened potential against the unscreened one stays below C_T6 * thresh;
/// measured 0.02 (Gaussian) and 3e-14 (cusp), see the header
static const double C_T6 = 1.0;

/// T5: the bound must lie above the analytic function less SLACK_T5 thresh, the projection's pointwise error;
/// measured worst 1.04 thresh (m2, two of its 20000 points), every other point within the bound
static const double SLACK_T5 = 2.0;

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

/// the Coulomb potential of exp(-a r^2): (pi/a)^{3/2} erf(sqrt(a) r) / r, 2 pi/a at the origin
class GaussPotential : public FunctionFunctorInterface<double, D> {
    double a;
public:
    explicit GaussPotential(double a) : a(a) {}
    double operator()(const Vector<double, D>& r) const override {
        const double rr = r.normf();
        if (rr < 1.e-10) return 2.0 * constants::pi / a;
        return std::pow(constants::pi / a, 1.5) * std::erf(std::sqrt(a) * rr) / rr;
    }
    std::vector<Vector<double,D>> special_points() const override { return {Vector<double,D>(0.0)}; }
};

static real_function_3d project(World& world, const std::shared_ptr<FunctionFunctorInterface<double,D>>& functor,
                                const double thresh) {
    real_function_3d f = real_factory_3d(world).functor(functor).thresh(thresh);
    f.truncate(thresh);
    return f;
}

/// blocks that passed the source test, and blocks the screen skipped, since the last reset, all ranks
static std::array<long, 2> screen_counts(World& world) {
    std::array<long, 2> n{long(ApplyScreenStats::blocks), long(ApplyScreenStats::skipped)};
    world.gop.sum(n.data(), n.size());
    return n;
}

/// one apply through the vector overload, with `screen` (none if empty); the counters cover this apply
template <typename T>
static Function<T, D> apply_with(World& world, const real_convolution_3d& P, const Function<T, D>& f,
                                 const apply_screen_t<D>& screen, std::array<long, 2>* counts = nullptr) {
    world.gop.fence();
    ApplyScreenStats::reset();
    world.gop.fence();
    std::vector<apply_screen_t<D>> screens;
    if (screen) screens.push_back(screen);
    Function<T, D> u = apply(world, P, std::vector<Function<T, D>>{f}, screens)[0];
    if (counts) *counts = screen_counts(world);
    return u;
}

/// a screen that never screens: the multiplier is unbounded everywhere, so the apply runs as without a screen
/// while the counters still see every block that passes the source test
static apply_screen_t<D> never_screen() {
    return [](const Key<D>&) { return SupBoundScreen<D>::unbounded(); };
}

/// the screen built from the multipliers' sup bounds
static apply_screen_t<D> sup_screen(const std::vector<real_function_3d>& m) {
    std::vector<std::shared_ptr<const SupBoundMap<D>>> sup;
    for (const auto& mi : m) sup.push_back(std::make_shared<const SupBoundMap<D>>(mi.get_impl()->sup_bound_map()));
    return SupBoundScreen<D>(std::move(sup));
}

/// the screen of T3: `inside` in the boxes that intersect B = [4,8] x [-2,2] x [-2,2], `outside` elsewhere
static apply_screen_t<D> box_b_screen(const double inside, const double outside) {
    return [inside, outside](const Key<D>& key) -> double {
        const Tensor<double>& cell = FunctionDefaults<D>::get_cell();
        const Tensor<double>& width = FunctionDefaults<D>::get_cell_width();
        const double blo[D] = {4.0, -2.0, -2.0}, bhi[D] = {8.0, 2.0, 2.0};
        const double h = std::pow(0.5, double(key.level()));
        bool intersects = true;
        for (std::size_t d = 0; d < D; ++d) {
            const double lo = cell(long(d), 0L) + width(long(d)) * h * double(key.translation()[d]);
            const double hi = lo + width(long(d)) * h;
            intersects = intersects and (hi > blo[d]) and (lo < bhi[d]);
        }
        return intersects ? inside : outside;
    };
}

static void set_defaults(const int k, const double thresh) {
    FunctionDefaults<D>::set_k(k);
    FunctionDefaults<D>::set_thresh(thresh);
    FunctionDefaults<D>::set_cubic_cell(-20.0, 20.0);
    FunctionDefaults<D>::set_refine(true);
    FunctionDefaults<D>::set_initial_level(2);
    FunctionDefaults<D>::set_truncate_mode(1);
}

/// the operands of T1, T2, T4 at k 8, thresh 1e-6
struct operands {
    real_function_3d g, uexact, ufull;
    std::vector<real_function_3d> m;     ///< the multipliers, reconstructed
    std::vector<std::string> mname;
    long blocks_ref = 0;                 ///< blocks the unscreened shell loops examine (T1)
};

static operands make_operands(World& world, const real_convolution_3d& P, const double thresh) {
    typedef std::shared_ptr<FunctionFunctorInterface<double,D>> functorT;
    const double a = 1.5;
    operands o;
    o.g = project(world, functorT(new Gauss(a, Vector<double,D>(0.0))), thresh);
    o.uexact = project(world, functorT(new GaussPotential(a)), THRESH_REF);
    o.ufull = apply_with(world, P, o.g, apply_screen_t<D>());
    Vector<double,D> R1(0.0), R2(0.0);
    R1[0] = 3.0;
    R2[1] = -2.5;  R2[2] = 0.5;
    o.m.push_back(project(world, functorT(new Gauss(1.0, R1)), thresh));
    o.m.push_back(project(world, functorT(new Gauss(2.0, R2)), thresh));
    o.mname = {"m1", "m2"};
    return o;
}

/// T1: a screen that never screens reproduces the unscreened apply
int test_unbounded_screen(World& world, const real_convolution_3d& P, operands& o) {
    test_output t("apply screen T1: an unbounded screen reproduces the unscreened apply");
    t.set_do_print(world.rank() == 0);
    std::array<long, 2> n{};
    real_function_3d u = apply_with(world, P, o.g, never_screen(), &n);
    o.blocks_ref = n[0];
    real_function_3d ufull = copy(o.ufull);
    // the same truncate on both; the accumulation order of the blocks is not fixed, so the two applies may
    // differ in the last bits
    u.truncate();
    ufull.truncate();
    const double diff = (u - ufull).norm2();
    const long nodes = long(u.tree_size()), nodes_full = long(ufull.tree_size());
    if (world.rank() == 0)
        printf("  T1 blocks %ld skipped %ld  nodes %ld / %ld  |u - u_full| %.3e\n", n[0], n[1], nodes, nodes_full, diff);
    t.checkpoint(n[0] > 0, "T1 the screen was consulted");
    t.checkpoint(n[1] == 0, "T1 nothing skipped");
    t.checkpoint(nodes == nodes_full, "T1 same number of nodes");
    t.checkpoint(diff < 1.e-12 * o.ufull.norm2(), "T1 same function");
    return t.end();
}

/// T2: the weighted error of the screened potential, for reconstructed and redundant multipliers
int test_accuracy(World& world, const real_convolution_3d& P, operands& o, const double thresh) {
    test_output t("apply screen T2: accuracy where the multipliers live");
    t.set_do_print(world.rank() == 0);

    std::vector<double> err_full;
    for (std::size_t i = 0; i < o.m.size(); ++i) err_full.push_back((o.m[i] * (o.ufull - o.uexact)).norm2());
    const long nodes_full = long(o.ufull.tree_size());

    std::vector<real_function_3d> mred;
    for (const auto& mi : o.m) {
        real_function_3d r = copy(mi);
        r.make_redundant();
        mred.push_back(r);
    }

    for (const std::string form : {"reconstructed", "redundant"}) {
        const std::vector<real_function_3d>& mult = (form == "redundant") ? mred : o.m;
        std::array<long, 2> n{};
        real_function_3d u = apply_with(world, P, o.g, sup_screen(mult), &n);
        const long nodes = long(u.tree_size());
        const std::string tag = "T2 " + form;
        if (world.rank() == 0)
            printf("  %s: blocks %ld skipped %ld (%.1f%%)  nodes %ld (unscreened %ld)\n", tag.c_str(),
                   n[0], n[1], 100.0 * double(n[1]) / double(std::max(n[0], 1L)), nodes, nodes_full);
        t.checkpoint(n[1] > 0, tag + " blocks skipped");
        t.checkpoint(n[0] >= o.blocks_ref, tag + " the shell loops reach as far as unscreened");
        t.checkpoint(nodes < nodes_full, tag + " the result tree is smaller");
        for (std::size_t i = 0; i < o.m.size(); ++i) {
            const double err = (o.m[i] * (u - o.uexact)).norm2();
            const double bound = std::max(2.0 * err_full[i], C_T2 * thresh);
            if (world.rank() == 0)
                printf("    %s: weighted error %.3e  unscreened %.3e  C %.2f  bound %.2e\n", o.mname[i].c_str(),
                       err, err_full[i], err / thresh, bound);
            t.checkpoint(err < bound, tag + " " + o.mname[i] + " weighted error");
        }
    }
    return t.end();
}

/// T3, T3b: the shell loop must reach the boxes where the multiplier lives
int test_shell_loop(World& world) {
    test_output t("apply screen T3: the shell loop reaches a distant multiplier");
    t.set_do_print(world.rank() == 0);

    const int k = 6;
    const double thresh = 1.e-4;
    set_defaults(k, thresh);
    typedef std::shared_ptr<FunctionFunctorInterface<double,D>> functorT;
    real_convolution_3d P = CoulombOperator(world, 1.e-4, thresh);
    real_function_3d g = project(world, functorT(new Gauss(1.5, Vector<double,D>(0.0))), thresh);
    Vector<double,D> centre(0.0);
    centre[0] = 6.0;
    // weighs B: one at its centre, 6e-6 at its faces, truncated to zero not far beyond
    real_function_3d chi = project(world, functorT(new Gauss(3.0, centre)), thresh);

    real_function_3d ufull = apply_with(world, P, g, apply_screen_t<D>());
    std::array<long, 2> nref{};
    apply_with(world, P, g, never_screen(), &nref);
    const double ref = (chi * ufull).norm2();
    if (world.rank() == 0) printf("  T3 unscreened: blocks %ld, |chi u| %.3e\n", nref[0], ref);
    for (const double outside : {1.e-12, 1.e-3, 1.e-2}) {
        std::array<long, 2> n{};
        real_function_3d u = apply_with(world, P, g, box_b_screen(1.0, outside), &n);
        const double err = (chi * (u - ufull)).norm2();
        char tag[64];
        snprintf(tag, sizeof(tag), "%s outside B %.0e", outside < 1.e-6 ? "T3 " : "T3b", outside);
        if (world.rank() == 0)
            printf("  %s: blocks %ld skipped %ld (%.1f%%)  error in B %.3e (relative %.1e)\n", tag, n[0], n[1],
                   100.0 * double(n[1]) / double(std::max(n[0], 1L)), err, err / ref);
        t.checkpoint(n[1] > 0, std::string(tag) + " blocks skipped");
        t.checkpoint(n[0] >= nref[0], std::string(tag) + " the shell loops reach as far as unscreened");
        t.checkpoint(err < REL_T3 * ref, std::string(tag) + " the potential in B is unchanged");
    }
    return t.end();
}

/// T4: a complex function under a real screen
int test_complex(World& world, const real_convolution_3d& P, operands& o, const double thresh) {
    test_output t("apply screen T4: complex function, real screen");
    t.set_do_print(world.rank() == 0);
    complex_function_3d gc = convert<double, double_complex, D>(o.g);
    std::array<long, 2> n{};
    complex_function_3d uc = apply_with(world, P, gc, sup_screen({o.m[0]}), &n);
    const double err = (o.m[0] * (real(uc) - o.uexact)).norm2();
    const double err_full = (o.m[0] * (o.ufull - o.uexact)).norm2();
    const double imag_norm = imag(uc).norm2();
    if (world.rank() == 0)
        printf("  T4 blocks %ld skipped %ld  weighted error of the real part %.3e (unscreened %.3e) C %.2f, |imag| %.1e\n",
               n[0], n[1], err, err_full, err / thresh, imag_norm);
    t.checkpoint(n[1] > 0, "T4 blocks skipped");
    t.checkpoint(n[0] >= o.blocks_ref, "T4 the shell loops reach as far as unscreened");
    t.checkpoint(err < std::max(2.0 * err_full, C_T2 * thresh), "T4 weighted error");
    t.checkpoint(imag_norm < 1.e-12, "T4 the imaginary part stays zero");
    return t.end();
}

/// exp(-Z r), refined to level 15 at the origin: a nuclear cusp
class Cusp : public FunctionFunctorInterface<double, D> {
    double Z;
public:
    explicit Cusp(double Z) : Z(Z) {}
    double operator()(const Vector<double, D>& r) const override { return std::exp(-Z * r.normf()); }
    std::vector<Vector<double,D>> special_points() const override { return {Vector<double,D>(0.0)}; }
    Level special_level() const override { return 15; }
};

/// median, 90th percentile and max
static std::array<double, 3> distribution(std::vector<double> v) {
    if (v.empty()) return {0.0, 0.0, 0.0};
    std::sort(v.begin(), v.end());
    auto at = [&v](const double q) { return v[std::min(v.size() - 1, std::size_t(q * double(v.size())))]; };
    return {at(0.5), at(0.9), v.back()};
}

/// the function in the box `key` at the point x in [0,1]^D of the box, from the box's scaling coefficients
static double eval_in_box(const Key<D>& key, const Tensor<double>& c, const std::array<double, D>& x) {
    const long k = c.dim(0);
    std::vector<double> p[D];
    for (std::size_t d = 0; d < D; ++d) {
        p[d].resize(k);
        legendre_scaling_functions(x[d], k, p[d].data());
    }
    double s = 0.0;
    for (long i = 0; i < k; ++i)
        for (long j = 0; j < k; ++j)
            for (long l = 0; l < k; ++l) s += c(i, j, l) * p[0][i] * p[1][j] * p[2][l];
    return s * std::pow(2.0, 0.5 * D * key.level()) / std::sqrt(FunctionDefaults<D>::get_cell_volume());
}

/// T5: the sup bound per box bounds the function
int test_sup_bound(World& world, const double thresh) {
    test_output t("apply screen T5: the sup bound per box bounds the function");
    t.set_do_print(world.rank() == 0);
    typedef std::shared_ptr<FunctionFunctorInterface<double,D>> functorT;
    Vector<double,D> R1(0.0), R2(0.0);
    R1[0] = 3.0;
    R2[1] = -2.5;  R2[2] = 0.5;
    const std::vector<std::pair<std::string, functorT>> functions{
            {"m1", functorT(new Gauss(1.0, R1))}, {"m2", functorT(new Gauss(2.0, R2))},
            {"cusp", functorT(new Cusp(35.0))}};
    const Tensor<double>& cell = FunctionDefaults<D>::get_cell();
    const Tensor<double>& width = FunctionDefaults<D>::get_cell_width();
    const Tensor<double> weights = SupBoundMap<D>::leaf_weights(FunctionDefaults<D>::get_k());
    std::mt19937 rng(20261014 + world.rank());
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    const int npoint = 50;

    for (const auto& [name, functor] : functions) {
        const real_function_3d f = project(world, functor, thresh);
        const FunctionImpl<double, D>& impl = *f.get_impl();
        double tbuild = wall_time();
        const SupBoundMap<D> S = impl.sup_bound_map();
        tbuild = wall_time() - tbuild;

        // leaves; failures of the bound at the quadrature points, at the corners and random points, against
        // the analytic function less SLACK_T5 thresh; boxes whose map entry is not what the construction
        // promises; boxes in the map
        std::array<double, 6> n{0, 0, 0, 0, 0, double(S.size())};
        double worst_exact = -1e300;    // max over the points of (|f analytic| - S) / thresh
        // per leaf with coefficients: its level, S / rms, S / max|values|, and whether max|values| > thresh
        std::vector<double> lev, s_rms, s_max, big;
        for (auto it = impl.get_coeffs().begin(); it != impl.get_coeffs().end(); ++it) {
            const Key<D>& key = it->first;
            const auto& node = it->second;
            const double s = S(key);
            if (node.has_children()) {
                // an interior box: the max over its children, or unbounded when one of them is not here
                double cmax = 0.0;
                for (KeyChildIterator<D> kit(key); kit; ++kit) cmax = std::max(cmax, S(kit.key()));
                if (s != cmax) n[4] += 1;
                continue;
            }
            n[0] += 1;
            if (not node.has_coeff()) {
                if (s != 0.0) n[4] += 1;
                continue;
            }
            const Tensor<double> c = node.coeff().full_tensor();
            bool below = false;
            if (s != SupBoundMap<D>::leaf_bound(key, c, weights) or S(KeyChildIterator<D>(key).key(), &below) != s
                or not below)
                n[4] += 1;
            const double vmax = impl.coeffs2values(key, c).absmax();
            const double vbox = FunctionDefaults<D>::get_cell_volume() * std::pow(2.0, -double(D) * key.level());
            const double rms = c.normf() / std::sqrt(vbox);
            if (vmax > s * (1.0 + 1.e-12)) n[1] += 1;
            if (rms > 0.0 and vmax > 0.0) {
                lev.push_back(key.level());
                s_rms.push_back(s / rms);
                s_max.push_back(s / vmax);
                big.push_back(vmax > thresh ? 1.0 : 0.0);
            }
            // the corners of the box, where the scaling functions peak, and random points
            for (int ip = 0; ip < 8 + npoint; ++ip) {
                std::array<double, D> x{};
                for (std::size_t d = 0; d < D; ++d) x[d] = (ip < 8) ? double((ip >> d) & 1) : unif(rng);
                if (std::abs(eval_in_box(key, c, x)) > s * (1.0 + 1.e-12)) n[2] += 1;
                Vector<double, D> r;
                for (std::size_t d = 0; d < D; ++d)
                    r[d] = cell(long(d), 0L) + width(long(d)) * std::pow(0.5, double(key.level()))
                                               * (double(key.translation()[d]) + x[d]);
                const double exact = std::abs((*functor)(r));
                worst_exact = std::max(worst_exact, (exact - s) / thresh);
                if (exact - SLACK_T5 * thresh > s) n[3] += 1;
            }
        }
        world.gop.sum(n.data(), n.size());
        world.gop.max(worst_exact);
        world.gop.max(tbuild);
        lev = world.gop.concat0(lev);
        s_rms = world.gop.concat0(s_rms);
        s_max = world.gop.concat0(s_max);
        big = world.gop.concat0(big);
        if (world.rank() == 0) {
            // the ratios of the leaves at `level` (all levels if negative), of all of them or of those with max|f| > thresh
            auto pick = [&](const std::vector<double>& v, const int level, const bool big_only) {
                std::vector<double> out;
                for (std::size_t i = 0; i < v.size(); ++i)
                    if ((level < 0 or lev[i] == level) and (not big_only or big[i] > 0.0)) out.push_back(v[i]);
                return out;
            };
            printf("  T5 %-4s: %.0f leaves, %.0f boxes in the map, built in %.1f ms (slowest rank)\n", name.c_str(),
                   n[0], n[5], 1e3 * tbuild);
            printf("    level  leaves (max|f| > thresh)   S/rms median  p90    max  (of those)   S/max|values| median  p90   max\n");
            int lmin = 100, lmax = -1;
            for (double l : lev) { lmin = std::min(lmin, int(l)); lmax = std::max(lmax, int(l)); }
            for (int l = lmin - 1; l <= lmax; ++l) {
                // the first line is all levels together
                const int level = (l < lmin) ? -1 : l;
                const auto all = pick(s_rms, level, false), sb = pick(s_rms, level, true), mb = pick(s_max, level, true);
                if (all.empty()) continue;
                const auto dr = distribution(sb), dm = distribution(mb);
                printf("    %5s  %6zu (%6zu)              %7.2f %6.2f %7.2f              %7.2f %6.2f %6.2f\n",
                       level < 0 ? "all" : std::to_string(level).c_str(), all.size(), sb.size(), dr[0], dr[1], dr[2],
                       dm[0], dm[1], dm[2]);
            }
            const auto dr = distribution(s_rms), dm = distribution(s_max);
            printf("    every leaf: S/rms median %.2f p90 %.2f max %.2f, S/max|values| median %.2f p90 %.2f max %.2f\n",
                   dr[0], dr[1], dr[2], dm[0], dm[1], dm[2]);
            printf("    S below: the quadrature max %.0f, the corners or random points %.0f, |f analytic| - %.0f thresh %.0f;"
                   " worst (|f analytic| - S) / thresh %.2f; map entries off %.0f\n", n[1], n[2], SLACK_T5, n[3],
                   worst_exact, n[4]);
        }
        t.checkpoint(n[0] > 0, "T5 " + name + " has leaves");
        t.checkpoint(n[1] == 0, "T5 " + name + " the bound is above the values at the quadrature points");
        t.checkpoint(n[2] == 0, "T5 " + name + " the bound is above the function at the corners and random points");
        t.checkpoint(n[3] == 0, "T5 " + name + " the bound is above the analytic function less thresh");
        t.checkpoint(n[4] == 0, "T5 " + name + " interior boxes hold the max over their children, boxes below a leaf the leaf's");
    }
    return t.end();
}

/// T6: the exchange's own geometry. The source is the product of the multiplier with itself, as in the
/// diagonal product bra_i vf_i of the exchange (vf = bra = ket), so the multiplier is O(1) on top of the
/// source at every level of the source. Two multipliers: a Gaussian, whose tree stops at level 5, and a
/// cusp, refined to level 13 at the nucleus, where a box's volume is far below one bohr^3.
int test_exchange_geometry(World& world, const real_convolution_3d& P, const double thresh) {
    test_output t("apply screen T6: the exchange's geometry, the multiplier on top of the source");
    t.set_do_print(world.rank() == 0);
    typedef std::shared_ptr<FunctionFunctorInterface<double,D>> functorT;
    const std::vector<std::pair<std::string, functorT>> multipliers{
            {"gauss", functorT(new Gauss(1.5, Vector<double,D>(0.0)))}, {"cusp", functorT(new Cusp(35.0))}};
    for (const auto& [name, functor] : multipliers) {
        real_function_3d s = project(world, functor, thresh);
        real_function_3d g = s * s;
        g.truncate();
        s.reconstruct();

        const real_function_3d ufull = apply_with(world, P, g, apply_screen_t<D>());
        std::array<long, 2> nref{};
        apply_with(world, P, g, never_screen(), &nref);
        std::array<long, 2> n{};
        const real_function_3d u = apply_with(world, P, g, sup_screen({s}), &n);
        const double err = (s * (u - ufull)).norm2();
        const std::string tag = "T6 " + name;
        if (world.rank() == 0)
            printf("  %s: unscreened blocks %ld; screened blocks %ld skipped %ld (%.1f%%)  weighted error against the"
                   " unscreened potential %.3e  C %.2f\n", tag.c_str(), nref[0], n[0], n[1],
                   100.0 * double(n[1]) / double(std::max(n[0], 1L)), err, err / thresh);
        t.checkpoint(n[0] >= nref[0], tag + " the shell loops reach as far as unscreened");
        t.checkpoint(n[1] > 0, tag + " blocks skipped");
        t.checkpoint(err < C_T6 * thresh, tag + " weighted error");
    }
    return t.end();
}

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    startup(world, argc, argv, true);
    if (apply_screen_disabled()) {
        if (world.rank() == 0) print("test_apply_screen: unset MAD_APPLY_SCREEN, the screen is switched off");
        finalize();
        return 1;
    }

    int success = 0;
    {
        const double thresh = 1.e-6;
        set_defaults(8, thresh);
        real_convolution_3d P = CoulombOperator(world, 1.e-4, thresh);
        operands o = make_operands(world, P, thresh);
        success += test_unbounded_screen(world, P, o);
        success += test_accuracy(world, P, o, thresh);
        success += test_complex(world, P, o, thresh);
        success += test_sup_bound(world, thresh);
        success += test_exchange_geometry(world, P, thresh);
    }
    success += test_shell_loop(world);

    world.gop.fence();
    finalize();
    return success;
}
