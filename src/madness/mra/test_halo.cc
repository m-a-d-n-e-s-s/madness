/// Unit tests for the derivative operator's neighbor-halo pre-staging and pool submission.
///
/// Neither may change results, so the reference is the derivative taken by the plain path and
/// agreement with it is required to round-off, not to a tolerance.
///
/// The operands are cusped so the tree is deep and unevenly refined, which is what makes
/// neighbors resolve by all three routes the halo must reproduce: a same-level leaf, an interior
/// node to descend into, and a coarser ancestor it does not cover and which falls back to the fetch.

#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/world/test_utilities.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace madness;

static const std::size_t D = 3;

/// The paths perform identical arithmetic, so any difference is a defect. Round-off rather than
/// zero only to avoid depending on the order the difference itself is accumulated in.
static const double EXACT = 1.e-13;

/// c * exp(-a |r-R|), cusped at R
class Slater : public FunctionFunctorInterface<double, D> {
    double a, c;
    Vector<double, D> R;
public:
    Slater(double a, double c, const Vector<double, D>& R) : a(a), c(c), R(R) {}
    double operator()(const Vector<double, D>& r) const override {
        double r2 = 0.0;
        for (std::size_t i = 0; i < D; ++i) r2 += (r[i]-R[i]) * (r[i]-R[i]);
        return c * std::exp(-a * std::sqrt(r2));
    }
    std::vector<Vector<double,D>> special_points() const override { return {R}; }
};

/// Cusps at different places and of different widths, so the ranks own unequal trees.
static std::vector<Function<double,D>> make_operands(World& world) {
    std::vector<Vector<double,D>> centres = {
        {{ 0.0,  0.0,  0.0}},
        {{ 1.3, -0.7,  0.4}},
        {{-0.9,  1.1, -1.2}},
        {{ 2.1,  0.3, -0.6}},
    };
    std::vector<double> alpha = {1.0, 1.7, 2.5, 0.8};
    std::vector<Function<double,D>> v;
    v.reserve(centres.size());
    for (std::size_t i = 0; i < centres.size(); ++i)
        v.emplace_back(FunctionFactory<double,D>(world)
                       .functor(std::shared_ptr<FunctionFunctorInterface<double,D>>(
                                    new Slater(alpha[i], 1.0 + 0.1*double(i), centres[i]))));
    reconstruct(world, v);
    return v;
}

/// Differentiate v along every axis, optionally staging the halo and/or submitting in parallel.
static std::vector<Function<double,D>>
differentiate(World& world,
              std::vector<std::shared_ptr<Derivative<double,D>>>& grad,
              const std::vector<Function<double,D>>& v,
              bool halo, bool parallel_submit,
              std::size_t& staged) {
    for (std::size_t axis = 0; axis < D; ++axis) grad[axis]->parallel_submit_ = parallel_submit;

    if (halo) {
        for (const auto& f : v) f.get_impl()->halo_enable();
        world.gop.fence();
        for (std::size_t axis = 0; axis < D; ++axis)
            for (const auto& f : v) grad[axis]->stage_halo(f.get_impl().get());
        world.gop.fence();
        staged = 0;
        for (const auto& f : v) staged += f.get_impl()->halo_size();
    }

    std::vector<Function<double,D>> dv;
    for (std::size_t axis = 0; axis < D; ++axis) {
        std::vector<Function<double,D>> d = apply(world, *grad[axis], v, false);
        dv.insert(dv.end(), d.begin(), d.end());
    }
    world.gop.fence();

    if (halo) for (const auto& f : v) f.get_impl()->halo_clear();
    return dv;
}

/// The halo and the parallel submission, separately and together, against the plain path.
static int test_equivalence(World& world, const std::vector<Function<double,D>>& v) {
    test_output t("neighbor halo and parallel submission reproduce the plain derivative");
    t.set_cout_to_logger();
    std::vector<std::shared_ptr<Derivative<double,D>>> grad = gradient_operator<double,D>(world);

    std::size_t staged = 0;
    std::vector<Function<double,D>> ref = differentiate(world, grad, v, false, false, staged);

    struct cell { bool halo, submit; const char* name; };
    const std::vector<cell> cells = {
        {true,  false, "halo"},
        {false, true,  "parallel submit"},
        {true,  true,  "halo + parallel submit"},
    };

    for (const cell& c : cells) {
        staged = 0;
        std::vector<Function<double,D>> dv = differentiate(world, grad, v, c.halo, c.submit, staged);
        double err = 0.0;
        for (std::size_t i = 0; i < ref.size(); ++i)
            err = std::max(err, (dv[i] - ref[i]).norm2());
        if (world.rank() == 0)
            printf("  %-24s max |d - d_ref| = %.3e   staged nodes %zu\n", c.name, err, staged);
        t.checkpoint(err, EXACT, std::string(c.name) + " matches the plain derivative");

        // a halo that moved nothing would pass the comparison above without exercising anything.
        // On one rank every consumer is local and nothing is pushed, by design.
        if (c.halo && world.size() > 1)
            t.checkpoint(staged > 0, std::string(c.name) + " staged a non-empty halo");
    }
    return t.end();
}

/// A staged halo must be inert for a lone derivative too, not just for the vector form the
/// kinetic-energy matrix uses.
static int test_single_function(World& world, const std::vector<Function<double,D>>& v) {
    test_output t("a staged halo does not disturb a single-function derivative");
    t.set_cout_to_logger();
    std::vector<std::shared_ptr<Derivative<double,D>>> grad = gradient_operator<double,D>(world);

    const Function<double,D>& f = v[0];
    Function<double,D> ref = (*grad[0])(f, true);

    f.get_impl()->halo_enable();
    world.gop.fence();
    grad[0]->stage_halo(f.get_impl().get());
    world.gop.fence();
    Function<double,D> df = (*grad[0])(f, true);
    const std::size_t staged = f.get_impl()->halo_size();
    f.get_impl()->halo_clear();

    const double err = (df - ref).norm2();
    if (world.rank() == 0)
        printf("  single function          max |d - d_ref| = %.3e   staged nodes %zu\n", err, staged);
    t.checkpoint(err, EXACT, "single-function derivative matches");
    return t.end();
}

/// halo_clear must restore the plain path, so that a stale halo cannot leak into a later
/// derivative of the same function.
static int test_clear(World& world, const std::vector<Function<double,D>>& v) {
    test_output t("halo_clear restores the plain path");
    t.set_cout_to_logger();
    std::vector<std::shared_ptr<Derivative<double,D>>> grad = gradient_operator<double,D>(world);

    const Function<double,D>& f = v[0];
    Function<double,D> ref = (*grad[0])(f, true);

    f.get_impl()->halo_enable();
    world.gop.fence();
    grad[0]->stage_halo(f.get_impl().get());
    world.gop.fence();
    f.get_impl()->halo_clear();
    t.checkpoint(f.get_impl()->halo_enabled() == false, "halo is disabled after clear");

    Function<double,D> df = (*grad[0])(f, true);
    const double err = (df - ref).norm2();
    if (world.rank() == 0) printf("  after clear              max |d - d_ref| = %.3e\n", err);
    t.checkpoint(err, EXACT, "derivative after clear matches");
    return t.end();
}

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    startup(world, argc, argv, true);

    FunctionDefaults<D>::set_k(6);
    FunctionDefaults<D>::set_thresh(1.e-6);
    FunctionDefaults<D>::set_cubic_cell(-8.0, 8.0);
    FunctionDefaults<D>::set_refine(true);
    FunctionDefaults<D>::set_initial_level(2);

    int success = 0;
    {
        std::vector<Function<double,D>> v = make_operands(world);
        success += test_equivalence(world, v);
        success += test_single_function(world, v);
        success += test_clear(world, v);
    }

    world.gop.fence();
    finalize();
    return success;
}
