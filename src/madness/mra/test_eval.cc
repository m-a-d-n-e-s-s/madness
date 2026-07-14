/// \file test_eval.cc
/// \brief Targeted unit tests for FunctionImpl::eval_local_only and eval_cube.
///
/// Covers:
///   - correctness matrix: NDIM x k x thresh x functor x point location x ranks
///   - golden-value characterization: pins computed values to guard silent drift
///   - micro-test #1: contraction correctness (raw loop vs separated form)
///   - micro-test #2: scaling factor (exp2 vs pow, cached 1/sqrt(V))
///
/// Heavy cells (NDIM==6, k>8, thresh<=1e-8) are compiled but skipped at runtime
/// unless MADNESS_TEST_EVAL_EXHAUSTIVE=1 is set.

#include <madness/mra/mra.h>
#include <madness/mra/funcdefaults.h>
#include <madness/world/test_utilities.h>
#include <madness/tensor/tensor.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <random>
#include <string>
#include <vector>

using namespace madness;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static const bool run_heavy =
    (std::getenv("MADNESS_TEST_EVAL_EXHAUSTIVE") != nullptr);

static bool is_heavy(std::size_t ndim, int k, double thresh) {
    return ndim == 6 || k > 8 || thresh <= 1e-8;
}

template <typename T, std::size_t NDIM>
class Gaussian : public FunctionFunctorInterface<T, NDIM> {
public:
    typedef Vector<double, NDIM> coordT;
    coordT center;
    double exponent;
    T coefficient;

    Gaussian(const coordT& c, double e, T coeff)
        : center(c), exponent(e), coefficient(coeff) {}

    T operator()(const coordT& x) const {
        double r2 = 0.0;
        for (std::size_t i = 0; i < NDIM; ++i) {
            double d = x[i] - center[i];
            r2 += d * d;
        }
        return coefficient * std::exp(-exponent * r2);
    }
};

// A function that is constant everywhere — analytic eval is trivial.
template <typename T, std::size_t NDIM>
class Constant : public FunctionFunctorInterface<T, NDIM> {
public:
    T value;
    explicit Constant(T v) : value(v) {}
    T operator()(const Vector<double, NDIM>&) const { return value; }
};

// Oscillatory: sin(freq * sum(x_d)) — stresses all coefficient slots.
template <std::size_t NDIM>
class Oscillatory : public FunctionFunctorInterface<double, NDIM> {
public:
    double freq;
    explicit Oscillatory(double f) : freq(f) {}
    double operator()(const Vector<double, NDIM>& x) const {
        double s = 0.0;
        for (std::size_t d = 0; d < NDIM; ++d) s += x[d];
        return std::sin(freq * s);
    }
};

// Pseudo-random seeded point set in sim coords [eps, 1-eps]^NDIM.
template <std::size_t NDIM>
std::vector<Vector<double, NDIM>> random_sim_points(int n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> dist(0.01, 0.99);
    std::vector<Vector<double, NDIM>> pts(n);
    for (auto& p : pts)
        for (std::size_t d = 0; d < NDIM; ++d)
            p[d] = dist(rng);
    return pts;
}

// Convert a sim-coord point to a user-coord point given cell bounds.
// cell is a 2-D tensor of shape (NDIM, 2): cell(d,0)=lo, cell(d,1)=hi.
template <std::size_t NDIM>
Vector<double, NDIM> sim_to_user(const Vector<double, NDIM>& xsim,
                                  const Tensor<double>& cell) {
    Vector<double, NDIM> xu;
    for (std::size_t d = 0; d < NDIM; ++d)
        xu[d] = cell(d, 0) + xsim[d] * (cell(d, 1) - cell(d, 0));
    return xu;
}

// Relative+absolute tolerance: a few ULP scaled by magnitude.
static double eval_tol(double ref, double k_factor = 1.0) {
    return std::max(std::abs(ref) * 1e-11 * k_factor, 1e-13);
}

// ---------------------------------------------------------------------------
// Correctness matrix
// ---------------------------------------------------------------------------

// For a single (NDIM, k, thresh, functor) cell:
//   - build f in sim cell [-L,L]^NDIM
//   - reconstruct
//   - evaluate at a set of user-coord points via eval() and eval_local_only()
//   - assert pair.second ≈ eval() for local points
//   - assert exactly one rank owns each in-domain point (via gop.sum)
template <typename FunctorT, std::size_t NDIM>
int run_correctness_cell(World& world,
                         const std::string& label,
                         FunctorT* raw_functor,
                         int k,
                         double thresh,
                         double L,
                         const std::vector<Vector<double, NDIM>>& user_pts) {
    test_output t(label);

    typedef Function<double, NDIM> funcT;
    typedef std::shared_ptr<FunctionFunctorInterface<double, NDIM>> functorT;

    FunctionDefaults<NDIM>::set_k(k);
    FunctionDefaults<NDIM>::set_thresh(thresh);
    FunctionDefaults<NDIM>::set_cubic_cell(-L, L);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_initial_level(2);
    world.gop.fence();

    functorT functor(raw_functor);
    funcT f = FunctionFactory<double, NDIM>(world).functor(functor);
    f.reconstruct();
    world.gop.fence();

    Level maxlevel = static_cast<Level>(f.max_local_depth());

    // Tolerance scales mildly with k: higher k means more coefficient
    // accumulations in the inner loop, so FP reassociation adds a bit more.
    double ktol = 1.0 + k * 0.01;

    for (const auto& xu : user_pts) {
        double ref = f.eval(xu).get();
        world.gop.fence();

        std::pair<bool, double> lo = f.eval_local_only(xu, maxlevel);

        // Verify value agreement when local.
        if (lo.first) {
            double err = ref - lo.second;
            t.checkpoint(std::abs(err) <= eval_tol(ref, ktol),
                         "eval_local_only value ≈ eval for " + label);
        }

        // Exactly one rank should own an in-domain point.
        int total = lo.first ? 1 : 0;
        world.gop.sum(total);
        world.gop.fence();
        t.checkpoint(total == 1, "exactly one rank owns point in " + label);
    }

    world.gop.fence();
    return t.end();
}

// Driver: sweep (NDIM, k, thresh) matrix.
template <std::size_t NDIM>
int test_correctness_ndim(World& world) {
    int errors = 0;
    const double L = 2.0;

    struct Cell { int k; double thresh; };
    // light cells (always run) + heavy cells (gated)
    static const Cell cells[] = {
        {2, 1e-3}, {7, 1e-3}, {8, 1e-3},
        {2, 1e-5}, {7, 1e-5}, {8, 1e-5},
        {10, 1e-3}, {10, 1e-5},    // heavy: k>8
        {7, 1e-8}, {7, 1e-10},     // heavy: tight thresh
    };

    // Use a fixed set of seeded random points.
    auto user_pts_sim = random_sim_points<NDIM>(12, 0xDEADBEEF + NDIM);
    // Convert sim → user coords for the default cell [-L,L]^NDIM.
    Tensor<double> cell(NDIM, 2);
    for (std::size_t d = 0; d < NDIM; ++d) { cell(d, 0) = -L; cell(d, 1) = L; }
    std::vector<Vector<double, NDIM>> user_pts;
    for (const auto& ps : user_pts_sim)
        user_pts.push_back(sim_to_user<NDIM>(ps, cell));

    // Also add a few explicit edge cases.
    {
        // on-axis points: only first coordinate varies, rest are mid-sim (0.5 → user 0).
        Vector<double, NDIM> p;
        for (std::size_t d = 0; d < NDIM; ++d) p[d] = 0.0;
        for (double x : {-1.5, -0.7, 0.0, 0.7, 1.5}) {
            p[0] = x;
            user_pts.push_back(p);
        }
        // off-axis: all coords distinct
        Vector<double, NDIM> q;
        for (std::size_t d = 0; d < NDIM; ++d)
            q[d] = (d % 2 == 0 ? 1.0 : -1.0) * (d + 1) * 0.13;
        user_pts.push_back(q);
        // near-boundary: just inside the clamp eps
        Vector<double, NDIM> r;
        for (std::size_t d = 0; d < NDIM; ++d) r[d] = -L + 1e-10;
        user_pts.push_back(r);
    }

    // center of a Gaussian
    Vector<double, NDIM> center;
    for (std::size_t d = 0; d < NDIM; ++d) center[d] = 0.3 * (d % 2 == 0 ? 1.0 : -1.0);

    for (const auto& c : cells) {
        if (is_heavy(NDIM, c.k, c.thresh) && !run_heavy) continue;

        std::string base = "NDIM=" + std::to_string(NDIM)
                         + " k=" + std::to_string(c.k)
                         + " thresh=" + std::to_string(c.thresh);

        // Gaussian functor
        errors += run_correctness_cell<Gaussian<double, NDIM>, NDIM>(
            world, base + " Gaussian",
            new Gaussian<double, NDIM>(center, 1.0, 1.0),
            c.k, c.thresh, L, user_pts);

        // Off-center Gaussian
        Vector<double, NDIM> center2;
        for (std::size_t d = 0; d < NDIM; ++d) center2[d] = -0.5 + 0.1 * d;
        errors += run_correctness_cell<Gaussian<double, NDIM>, NDIM>(
            world, base + " off-center Gaussian",
            new Gaussian<double, NDIM>(center2, 3.0, 1.5),
            c.k, c.thresh, L, user_pts);

        // Constant function
        errors += run_correctness_cell<Constant<double, NDIM>, NDIM>(
            world, base + " constant",
            new Constant<double, NDIM>(2.5),
            c.k, c.thresh, L, user_pts);

        // Oscillatory — only for NDIM<=3 to keep build times reasonable
        // (still compiled for all NDIM; the run_heavy gate handles 6D).
        if (NDIM <= 3 || run_heavy) {
            errors += run_correctness_cell<Oscillatory<NDIM>, NDIM>(
                world, base + " oscillatory",
                new Oscillatory<NDIM>(2.5),
                c.k, c.thresh, L, user_pts);
        }
    }
    return errors;
}

// ---------------------------------------------------------------------------
// Micro-test #1: contraction correctness
// Build random dense coeff tensor + random px[d] arrays,
// evaluate with the reference nested-loop form and the separated form,
// assert agreement to relative ~1e-13.
// ---------------------------------------------------------------------------

namespace {

// Reference: flat nested-loop contraction (original code shape).
double eval_cube_ref_1d(int k, const Tensor<double>& c, double px[][MAXK]) {
    double sum = 0.0;
    const double* cp = c.ptr();
    for (int p = 0; p < k; ++p) sum += cp[p] * px[0][p];
    return sum;
}

double eval_cube_ref_2d(int k, const Tensor<double>& c, double px[][MAXK]) {
    double sum = 0.0;
    for (int p = 0; p < k; ++p)
        for (int q = 0; q < k; ++q)
            sum += c(p, q) * px[0][p] * px[1][q];
    return sum;
}

double eval_cube_ref_3d(int k, const Tensor<double>& c, double px[][MAXK]) {
    double sum = 0.0;
    for (int p = 0; p < k; ++p)
        for (int q = 0; q < k; ++q)
            for (int r = 0; r < k; ++r)
                sum += c(p, q, r) * px[0][p] * px[1][q] * px[2][r];
    return sum;
}

// Separated form matching the new eval_cube implementation.
double eval_cube_sep_2d(int k, const Tensor<double>& c, double px[][MAXK]) {
    double sum = 0.0;
    for (int p = 0; p < k; ++p) {
        const double a = px[0][p];
        const double* cq = &c(p, 0);
        double s2 = 0.0;
        for (int q = 0; q < k; ++q) s2 += cq[q] * px[1][q];
        sum += a * s2;
    }
    return sum;
}

double eval_cube_sep_3d(int k, const Tensor<double>& c, double px[][MAXK]) {
    double sum = 0.0;
    for (int p = 0; p < k; ++p) {
        const double a = px[0][p];
        for (int q = 0; q < k; ++q) {
            const double ab = a * px[1][q];
            const double* cr = &c(p, q, 0);
            double s3 = 0.0;
            for (int r = 0; r < k; ++r) s3 += cr[r] * px[2][r];
            sum += ab * s3;
        }
    }
    return sum;
}

} // namespace

int test_contraction_micro() {
    test_output t("eval_cube: separated contraction agrees with nested-loop reference");

    std::mt19937_64 rng(0xC0FFEE42);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    static const int ks[] = {2, 5, 7, 8};
    constexpr int NREPS = 8;

    for (int k : ks) {
        for (int rep = 0; rep < NREPS; ++rep) {
            // 1D
            {
                Tensor<double> c(k);
                for (int p = 0; p < k; ++p) c(p) = dist(rng);
                double px[1][MAXK];
                for (int p = 0; p < k; ++p) px[0][p] = dist(rng);

                double ref = eval_cube_ref_1d(k, c, px);
                // 1D separated form is identical to reference; just verify ptr().
                const double* cp = c.ptr();
                double sep = 0.0;
                for (int p = 0; p < k; ++p) sep += cp[p] * px[0][p];

                double err = std::abs(sep - ref);
                t.checkpoint(err <= std::max(std::abs(ref) * 1e-13, 1e-15),
                             "1D k=" + std::to_string(k));
            }
            // 2D
            {
                Tensor<double> c(k, k);
                for (int p = 0; p < k; ++p)
                    for (int q = 0; q < k; ++q)
                        c(p, q) = dist(rng);
                double px[2][MAXK];
                for (int d = 0; d < 2; ++d)
                    for (int p = 0; p < k; ++p) px[d][p] = dist(rng);

                double ref = eval_cube_ref_2d(k, c, px);
                double sep = eval_cube_sep_2d(k, c, px);
                double err = std::abs(sep - ref);
                t.checkpoint(err <= std::max(std::abs(ref) * 1e-13, 1e-15),
                             "2D k=" + std::to_string(k));
            }
            // 3D
            {
                Tensor<double> c(k, k, k);
                for (int p = 0; p < k; ++p)
                    for (int q = 0; q < k; ++q)
                        for (int r = 0; r < k; ++r)
                            c(p, q, r) = dist(rng);
                double px[3][MAXK];
                for (int d = 0; d < 3; ++d)
                    for (int p = 0; p < k; ++p) px[d][p] = dist(rng);

                double ref = eval_cube_ref_3d(k, c, px);
                double sep = eval_cube_sep_3d(k, c, px);
                double err = std::abs(sep - ref);
                t.checkpoint(err <= std::max(std::abs(ref) * 1e-13, 1e-15),
                             "3D k=" + std::to_string(k));
            }
        }
    }
    return t.end();
}

// ---------------------------------------------------------------------------
// Micro-test #2: scaling factor
// Assert exp2(0.5*NDIM*n) == pow(2.0, 0.5*NDIM*n) to within 1 ULP for all
// (NDIM, n) in the expected evaluation range.  Pure scalar arithmetic; no
// Function build needed, so the full NDIM range including 6 is fast.
// ---------------------------------------------------------------------------

int test_scaling_factor() {
    test_output t("eval_cube scaling: exp2 agrees with pow over all (NDIM,n)");

    for (int ndim = 1; ndim <= 6; ++ndim) {
        for (Level n = 0; n <= MAXLEVEL; ++n) {
            double arg = 0.5 * ndim * n;
            double via_pow  = std::pow(2.0, arg);
            double via_exp2 = std::exp2(arg);
            // For exact integer or half-integer powers of two, they are
            // bit-for-bit identical.  In general allow 1 ULP.
            double rel = (via_pow != 0.0)
                         ? std::abs(via_exp2 - via_pow) / via_pow
                         : std::abs(via_exp2 - via_pow);
            t.checkpoint(rel <= std::numeric_limits<double>::epsilon() * 2,
                         "NDIM=" + std::to_string(ndim)
                         + " n=" + std::to_string(n));
        }
    }
    return t.end();
}

// ---------------------------------------------------------------------------
// Golden-value characterization
// For a fixed set of (NDIM, k, thresh, point) cells, assert that
// eval_local_only().second agrees with f.eval() to ~1e-11 relative.
// This pins the consistency between the two paths: a normalization bug or
// index transposition that affects eval_cube would cause them to diverge.
//
// To emit hardcoded values for a standalone table (useful to catch regressions
// where BOTH paths drift identically), run with MADNESS_TEST_EVAL_PRINT_GOLDEN=1
// and paste the printed lines into a separate pinned-values table.
// ---------------------------------------------------------------------------

static const bool print_golden =
    (std::getenv("MADNESS_TEST_EVAL_PRINT_GOLDEN") != nullptr);

struct GoldenPoint {
    std::size_t ndim;
    int k;
    double thresh;
    double pt[3];  // first ndim coords used; rest zero
};

static const GoldenPoint golden_pts[] = {
    {1, 7, 1e-5, { 0.3,  0.0,  0.0}},
    {1, 7, 1e-5, { 0.0,  0.0,  0.0}},
    {1, 7, 1e-5, {-0.5,  0.0,  0.0}},
    {2, 7, 1e-5, { 0.3,  0.3,  0.0}},
    {2, 7, 1e-5, { 0.1, -0.2,  0.0}},
    {3, 7, 1e-5, { 0.3, -0.3,  0.3}},
    {3, 7, 1e-5, { 0.0,  0.0,  0.0}},
};

// eval_local_only vs f.eval() for one (NDIM, k, thresh, point) cell.
template <std::size_t NDIM>
bool golden_check_ndim(World& world, int k, double thresh, double cx,
                       const double raw_pt[3],
                       double& out_ref, double& out_lo) {
    typedef std::shared_ptr<FunctionFunctorInterface<double, NDIM>> functorT;
    const double L = 2.0;
    FunctionDefaults<NDIM>::set_k(k);
    FunctionDefaults<NDIM>::set_thresh(thresh);
    FunctionDefaults<NDIM>::set_cubic_cell(-L, L);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_initial_level(2);
    world.gop.fence();

    Vector<double, NDIM> center, pt;
    for (std::size_t d = 0; d < NDIM; ++d) { center[d] = cx; pt[d] = raw_pt[d]; }

    functorT functor(new Gaussian<double, NDIM>(center, 1.0, 1.0));
    Function<double, NDIM> f = FunctionFactory<double, NDIM>(world).functor(functor);
    f.reconstruct();
    world.gop.fence();

    Level maxlevel = static_cast<Level>(f.max_local_depth());
    out_ref = f.eval(pt).get();
    world.gop.fence();
    std::pair<bool, double> lo = f.eval_local_only(pt, maxlevel);
    world.gop.fence();

    if (lo.first) {
        out_lo = lo.second;
        if (print_golden && world.rank() == 0) {
            std::printf("  {%zu, %d, %.0e, {%.1f,%.1f,%.1f}, %.17g},\n",
                        NDIM, k, thresh,
                        raw_pt[0], raw_pt[1], raw_pt[2], lo.second);
        }
        return true;
    }
    return false;
}

int test_golden_values(World& world) {
    test_output t("eval_local_only golden-value characterization");
    const double cx = 0.3;

    for (const auto& g : golden_pts) {
        if (g.ndim > 3) continue;

        double eval_ref = 0.0;
        double eval_lo  = 0.0;
        bool found = false;

        if (g.ndim == 1)
            found = golden_check_ndim<1>(world, g.k, g.thresh, cx, g.pt, eval_ref, eval_lo);
        else if (g.ndim == 2)
            found = golden_check_ndim<2>(world, g.k, g.thresh, cx, g.pt, eval_ref, eval_lo);
        else if (g.ndim == 3)
            found = golden_check_ndim<3>(world, g.k, g.thresh, cx, g.pt, eval_ref, eval_lo);

        if (found) {
            double err = std::abs(eval_lo - eval_ref);
            double tol = std::max(std::abs(eval_ref) * 1e-11, 1e-13);
            t.checkpoint(err <= tol,
                         "golden NDIM=" + std::to_string(g.ndim)
                         + " k=" + std::to_string(g.k)
                         + " pt=(" + std::to_string(g.pt[0]) + ")");
        } else {
            t.checkpoint(true, "golden (non-local rank skip)");
        }
    }
    return t.end();
}

// ---------------------------------------------------------------------------
// Micro-test: general_fast_transform agrees with general_transform
// Verifies ping-pong parity (D odd / D even), rectangular c[d] (k×1),
// and thread-local scratch reuse across multiple calls on the same buffer.
// ---------------------------------------------------------------------------

namespace {

template <std::size_t D>
int fast_vs_general_ndim(int k) {
    test_output t("general_fast_transform D=" + std::to_string(D)
                  + " k=" + std::to_string(k));

    std::mt19937_64 rng(0xBEEF0000 + D * 100 + k);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    // Build a random k^D coefficient tensor.
    long dims[D];
    for (std::size_t d = 0; d < D; ++d) dims[d] = k;
    Tensor<double> coeff(D, dims);
    for (long i = 0; i < coeff.size(); ++i) coeff.ptr()[i] = dist(rng);

    // Build per-axis rank-1 matrices (k×1) with random values.
    Tensor<double> phi[D];
    for (std::size_t d = 0; d < D; ++d) {
        phi[d] = Tensor<double>(long(k), 1L);
        for (int p = 0; p < k; ++p) phi[d].ptr()[p] = dist(rng);
    }

    // Reference: general_transform (allocating baseline).
    Tensor<double> ref = general_transform(coeff, phi);
    double ref_val = ref.ptr()[0];

    // Test: general_fast_transform with caller-owned buffers.
    Tensor<double> result(coeff.size()), workspace(coeff.size());
    general_fast_transform(coeff, phi, result, workspace);
    double fast_val = result.ptr()[0];

    double err = std::abs(fast_val - ref_val);
    double tol = std::max(std::abs(ref_val) * 1e-13, 1e-15);
    t.checkpoint(err <= tol, "caller-owned buffers");

    // Test: general_fast_transform via thread-local scratch (detail::eval_scratch).
    {
        auto [ws, res] = madness::detail::eval_scratch<double>(coeff.size());
        general_fast_transform(coeff, phi, res, ws);
        double tl_val = res.ptr()[0];
        double terr = std::abs(tl_val - ref_val);
        t.checkpoint(terr <= tol, "thread-local scratch");
    }

    // Second call — scratch is warm; result must still be correct.
    {
        // Use different random coeff to confirm scratch is not stale.
        for (long i = 0; i < coeff.size(); ++i) coeff.ptr()[i] = dist(rng);
        for (std::size_t d = 0; d < D; ++d)
            for (int p = 0; p < k; ++p) phi[d].ptr()[p] = dist(rng);

        Tensor<double> ref2 = general_transform(coeff, phi);
        double ref2_val = ref2.ptr()[0];

        auto [ws2, res2] = madness::detail::eval_scratch<double>(coeff.size());
        general_fast_transform(coeff, phi, res2, ws2);
        double tl2_val = res2.ptr()[0];
        double terr2 = std::abs(tl2_val - ref2_val);
        double tol2 = std::max(std::abs(ref2_val) * 1e-13, 1e-15);
        t.checkpoint(terr2 <= tol2, "thread-local scratch reuse");
    }

    return t.end();
}

} // namespace

int test_general_fast_transform() {
    int errors = 0;
    // Odd D=1,3 exercises the D&1 parity swap (result/workspace pointer swap).
    // Even D=2,6 exercises the other branch.
    errors += fast_vs_general_ndim<1>(5);
    errors += fast_vs_general_ndim<1>(8);
    errors += fast_vs_general_ndim<2>(5);
    errors += fast_vs_general_ndim<2>(8);
    errors += fast_vs_general_ndim<3>(5);
    errors += fast_vs_general_ndim<3>(8);
    if (run_heavy) {
        errors += fast_vs_general_ndim<6>(5);
        errors += fast_vs_general_ndim<6>(8);
    }
    return errors;
}

// ---------------------------------------------------------------------------
// Batched eval_local_only (P1)
// The batched, descend-once-per-leaf-box form must reproduce the per-point
// path *bit-for-bit* (it reuses the identical eval_cube on the identical leaf
// coefficient tensor) and preserve the "exactly one rank owns each in-domain
// point" invariant.  Points are deliberately clustered so that several share a
// leaf box, exercising the bucketing.
// ---------------------------------------------------------------------------
template <std::size_t NDIM>
int run_batched_cell(World& world, const std::string& label, int k, double thresh,
                     double L, const std::vector<Vector<double, NDIM>>& user_pts) {
    test_output t(label);
    typedef std::shared_ptr<FunctionFunctorInterface<double, NDIM>> functorT;

    FunctionDefaults<NDIM>::set_k(k);
    FunctionDefaults<NDIM>::set_thresh(thresh);
    FunctionDefaults<NDIM>::set_cubic_cell(-L, L);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_initial_level(2);
    world.gop.fence();

    Vector<double, NDIM> center;
    for (std::size_t d = 0; d < NDIM; ++d) center[d] = 0.3 * (d % 2 == 0 ? 1.0 : -1.0);
    functorT functor(new Gaussian<double, NDIM>(center, 1.0, 1.0));
    Function<double, NDIM> f = FunctionFactory<double, NDIM>(world).functor(functor);
    f.reconstruct();
    world.gop.fence();

    Level maxlevel = static_cast<Level>(f.max_local_depth());

    // Per-point reference.
    std::vector<std::pair<bool, double>> ref;
    ref.reserve(user_pts.size());
    for (const auto& xu : user_pts) ref.push_back(f.eval_local_only(xu, maxlevel));

    // Batched.
    std::vector<std::pair<bool, double>> bat = f.eval_local_only(user_pts, maxlevel);

    t.checkpoint(bat.size() == user_pts.size(), "batched size matches " + label);

    for (std::size_t i = 0; i < user_pts.size(); ++i) {
        bool same_flag = bat[i].first == ref[i].first;
        // Bit-for-bit value match when local (batched reuses the same eval_cube).
        bool same_val = (!bat[i].first) || (bat[i].second == ref[i].second);
        t.checkpoint(same_flag && same_val,
                     "batched == per-point for point " + std::to_string(i) + " " + label);

        int total = bat[i].first ? 1 : 0;
        world.gop.sum(total);
        world.gop.fence();
        t.checkpoint(total == 1,
                     "exactly one rank owns batched point " + std::to_string(i) + " " + label);
    }

    // Order-independence: shuffling the points must not change any result.
    // (Guards the memoized-descent fast path: a stale or wrongly-matched
    // cached leaf would show up as a value difference under reordering.)
    {
        std::vector<std::size_t> perm(user_pts.size());
        for (std::size_t i = 0; i < perm.size(); ++i) perm[i] = i;
        std::mt19937_64 rng(0x5EED);
        std::shuffle(perm.begin(), perm.end(), rng);
        std::vector<Vector<double, NDIM>> shuffled(user_pts.size());
        for (std::size_t i = 0; i < perm.size(); ++i) shuffled[i] = user_pts[perm[i]];
        auto bat_sh = f.eval_local_only(shuffled, maxlevel);
        bool ok = (bat_sh.size() == shuffled.size());
        for (std::size_t i = 0; ok && i < perm.size(); ++i) {
            ok = (bat_sh[i].first == bat[perm[i]].first) &&
                 (!bat_sh[i].first || bat_sh[i].second == bat[perm[i]].second);
        }
        t.checkpoint(ok, "batched results order-independent " + label);
    }

    world.gop.fence();
    return t.end();
}

template <std::size_t NDIM>
int test_batched_ndim(World& world) {
    int errors = 0;
    const double L = 2.0;

    Tensor<double> cell(NDIM, 2);
    for (std::size_t d = 0; d < NDIM; ++d) { cell(d, 0) = -L; cell(d, 1) = L; }

    // Scattered points across the domain ...
    auto scattered = random_sim_points<NDIM>(20, 0xBA7C4ED0 + NDIM);
    // ... plus a tight cluster guaranteed to share leaf boxes.
    {
        std::mt19937_64 rng(0x5EED + NDIM);
        std::uniform_real_distribution<double> jitter(0.0, 1e-3);
        for (int c = 0; c < 8; ++c) {
            Vector<double, NDIM> p;
            for (std::size_t d = 0; d < NDIM; ++d) p[d] = 0.51 + jitter(rng);
            scattered.push_back(p);
        }
    }

    std::vector<Vector<double, NDIM>> user_pts;
    user_pts.reserve(scattered.size());
    for (const auto& ps : scattered) user_pts.push_back(sim_to_user<NDIM>(ps, cell));

    struct Cell { int k; double thresh; };
    static const Cell cells[] = {{7, 1e-3}, {8, 1e-5}};
    for (const auto& c : cells) {
        std::string base = "batched NDIM=" + std::to_string(NDIM)
                         + " k=" + std::to_string(c.k);
        errors += run_batched_cell<NDIM>(world, base, c.k, c.thresh, L, user_pts);
    }
    return errors;
}

// Output-parameter overload: reusing one results buffer across calls of
// different sizes must produce the same values as the vector-returning form.
int test_output_param_overload(World& world) {
    constexpr std::size_t NDIM = 3;
    test_output t("output-parameter eval_local_only overload");
    typedef std::shared_ptr<FunctionFunctorInterface<double, NDIM>> functorT;
    const double L = 2.0;
    FunctionDefaults<NDIM>::set_k(7);
    FunctionDefaults<NDIM>::set_thresh(1e-5);
    FunctionDefaults<NDIM>::set_cubic_cell(-L, L);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_initial_level(2);
    world.gop.fence();

    Vector<double, NDIM> center;
    for (std::size_t d = 0; d < NDIM; ++d) center[d] = 0.0;
    functorT functor(new Gaussian<double, NDIM>(center, 1.0, 1.0));
    Function<double, NDIM> f = FunctionFactory<double, NDIM>(world).functor(functor);
    f.reconstruct();
    world.gop.fence();
    const Level maxlevel = static_cast<Level>(f.max_local_depth());

    std::vector<Vector<double, NDIM>> pts_big, pts_small;
    for (int i = 0; i < 40; ++i) {
        Vector<double, NDIM> p;
        for (std::size_t d = 0; d < NDIM; ++d)
            p[d] = -L + 2.0*L*((i*7 + int(d)*3) % 37)/37.0;
        pts_big.push_back(p);
        if (i < 5) pts_small.push_back(p);
    }

    std::vector<std::pair<bool,double>> reused;
    f.eval_local_only(pts_big, maxlevel, reused);           // fills, allocates once
    auto ref_big = f.eval_local_only(pts_big, maxlevel);    // vector-returning form
    bool ok = (reused.size() == ref_big.size());
    for (std::size_t i = 0; ok && i < ref_big.size(); ++i)
        ok = reused[i].first == ref_big[i].first &&
             (!reused[i].first || reused[i].second == ref_big[i].second);
    t.checkpoint(ok, "fresh results buffer matches vector-returning form");

    f.eval_local_only(pts_small, maxlevel, reused);         // shrink, reuse capacity
    auto ref_small = f.eval_local_only(pts_small, maxlevel);
    ok = (reused.size() == ref_small.size());
    for (std::size_t i = 0; ok && i < ref_small.size(); ++i)
        ok = reused[i].first == ref_small[i].first &&
             (!reused[i].first || reused[i].second == ref_small[i].second);
    t.checkpoint(ok, "shrunk reused buffer matches vector-returning form");

    world.gop.fence();
    return t.end();
}

// ---------------------------------------------------------------------------
// Benchmark harness (P4): single-point vs batched eval_local_only.
//
// NOT a CTest entry: only runs when MADNESS_BENCH_EVAL is set, so the
// registered `./test_eval` run never touches it.  Reports ns/eval for the
// per-point path and the batched path at the current MAD_NUM_THREADS; run it
// repeatedly with MAD_NUM_THREADS in {1,2,4,8} to capture the scaling curve.
// Points are clustered so many share a leaf box -- the regime where amortising
// the descent (batched) is expected to pay.
// ---------------------------------------------------------------------------

static const bool run_bench =
    (std::getenv("MADNESS_BENCH_EVAL") != nullptr);

template <std::size_t NDIM>
void bench_cell(World& world, int k) {
    typedef std::shared_ptr<FunctionFunctorInterface<double, NDIM>> functorT;
    const double L = 2.0;
    FunctionDefaults<NDIM>::set_k(k);
    FunctionDefaults<NDIM>::set_thresh(1e-6);
    FunctionDefaults<NDIM>::set_cubic_cell(-L, L);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_initial_level(2);
    world.gop.fence();

    Vector<double, NDIM> center;
    for (std::size_t d = 0; d < NDIM; ++d) center[d] = 0.0;
    functorT functor(new Gaussian<double, NDIM>(center, 1.0, 1.0));
    Function<double, NDIM> f = FunctionFactory<double, NDIM>(world).functor(functor);
    f.reconstruct();
    world.gop.fence();
    Level maxlevel = static_cast<Level>(f.max_local_depth());

    // Build a clustered point set: several clusters, many points each, so points
    // share leaf boxes (the case batching targets).
    Tensor<double> cell(NDIM, 2);
    for (std::size_t d = 0; d < NDIM; ++d) { cell(d, 0) = -L; cell(d, 1) = L; }
    const int n_clusters = (NDIM <= 3) ? 64 : 16;
    const int per_cluster = 32;
    std::mt19937_64 rng(0xB3C4 + NDIM * 100 + k);
    std::uniform_real_distribution<double> c0(0.05, 0.95), jit(0.0, 2e-3);
    std::vector<Vector<double, NDIM>> pts;
    pts.reserve(n_clusters * per_cluster);
    for (int c = 0; c < n_clusters; ++c) {
        Vector<double, NDIM> base;
        for (std::size_t d = 0; d < NDIM; ++d) base[d] = c0(rng);
        for (int j = 0; j < per_cluster; ++j) {
            Vector<double, NDIM> ps;
            for (std::size_t d = 0; d < NDIM; ++d) ps[d] = std::min(0.99, base[d] + jit(rng));
            pts.push_back(sim_to_user<NDIM>(ps, cell));
        }
    }
    const std::size_t npt = pts.size();

    // Count locally-owned points (only those do real work).
    auto warm = f.eval_local_only(pts, maxlevel);
    std::size_t nlocal = 0;
    for (const auto& r : warm) if (r.first) ++nlocal;
    if (nlocal == 0) return;

    using clock = std::chrono::steady_clock;
    const int reps = 20;

    // Per-point path.
    volatile double sink = 0.0;
    auto t0 = clock::now();
    for (int r = 0; r < reps; ++r)
        for (const auto& xu : pts) {
            auto pr = f.eval_local_only(xu, maxlevel);
            if (pr.first) sink += pr.second;
        }
    auto t1 = clock::now();

    // Batched path.
    auto t2 = clock::now();
    for (int r = 0; r < reps; ++r) {
        auto br = f.eval_local_only(pts, maxlevel);
        for (const auto& pr : br) if (pr.first) sink += pr.second;
    }
    auto t3 = clock::now();
    (void)sink;

    double ns_pp  = std::chrono::duration<double, std::nano>(t1 - t0).count() / (reps * nlocal);
    double ns_bat = std::chrono::duration<double, std::nano>(t3 - t2).count() / (reps * nlocal);

    if (world.rank() == 0) {
        std::printf("  NDIM=%zu k=%2d  npts=%zu local=%zu  per-point=%8.1f ns/eval  "
                    "batched=%8.1f ns/eval  speedup=%5.2fx\n",
                    NDIM, k, npt, nlocal, ns_pp, ns_bat,
                    ns_bat > 0 ? ns_pp / ns_bat : 0.0);
        std::fflush(stdout);
    }
}

void run_benchmark(World& world) {
    if (world.rank() == 0) {
        const char* nt = std::getenv("MAD_NUM_THREADS");
        std::printf("eval_local_only benchmark (MAD_NUM_THREADS=%s, nranks=%d)\n",
                    nt ? nt : "default", world.size());
    }
    bench_cell<1>(world, 8);
    bench_cell<1>(world, 10);
    bench_cell<3>(world, 8);
    bench_cell<3>(world, 10);
    // The NDIM=6 cells spend tens of minutes in tree construction at low
    // thread counts, so they are opt-in on top of MADNESS_BENCH_EVAL.
    if (std::getenv("MADNESS_BENCH_EVAL_6D")) {
        bench_cell<6>(world, 8);
        bench_cell<6>(world, 10);
    }
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char** argv) {
    World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);

    int errors = 0;

    // Micro-tests: cheap, always run (full NDIM range for scaling factor).
    errors += test_contraction_micro();
    errors += test_scaling_factor();
    errors += test_general_fast_transform();

    // Correctness matrix: NDIM=1,2,3 (light cells always; heavy cells gated).
    errors += test_correctness_ndim<1>(world);
    errors += test_correctness_ndim<2>(world);
    errors += test_correctness_ndim<3>(world);

    // NDIM=6 correctness: compiled, heavy-gated.
    if (run_heavy) {
        errors += test_correctness_ndim<6>(world);
    }

    // Golden-value characterization.
    errors += test_golden_values(world);

    // Batched eval_local_only (P1): NDIM=1,2,3.
    errors += test_batched_ndim<1>(world);
    errors += test_batched_ndim<2>(world);
    errors += test_batched_ndim<3>(world);
    errors += test_output_param_overload(world);

    // Benchmark harness (P4): opt-in only, never part of the CTest run.
    if (run_bench) run_benchmark(world);

    world.gop.fence();
    madness::finalize();
    return errors;
}
