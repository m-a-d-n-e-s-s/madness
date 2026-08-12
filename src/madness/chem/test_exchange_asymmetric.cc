/// \file test_exchange_asymmetric.cc
/// \brief tests the exchange operator with bra, ket and vf all distinct

/// The generic form of the operator is K f_j = sum_i ket_i G[bra_i f_j], with bra and ket a paired
/// set of one length and f a separate set of another. **No consumer has ever exercised that**: HF
/// exchange has bra == ket == f, and nemo, oep and znemo have bra != ket but still apply K to the
/// ket itself, so f == ket and the two lengths agree. molresponse is the first caller with f
/// genuinely distinct from ket.
///
/// moldft cannot test this at all -- it always applies K to its own orbitals -- so a unit test is the
/// only route to the case, which is why this file exists rather than a calculation.
///
/// The reference is `smallmem`, whose kernel is the plain double loop over distinct nocc and nf with
/// no symmetry anywhere. Comparing one algorithm against another makes the test self-checking: there
/// are no reference numbers to maintain, and both sides see the identical operands.

#include <madness/chem/exchangeoperator.h>
#include <madness/chem/SCFOperators.h>
#include <madness/world/test_utilities.h>

#include <vector>

using namespace madness;

namespace {

/// a gaussian at an offset origin, so bra, ket and vf are genuinely different functions rather
/// than scaled copies of one
struct shifted_gaussian {
    double exponent = 1.0;
    coord_3d origin{};
    shifted_gaussian() = default;
    shifted_gaussian(const double e, const coord_3d& o) : exponent(e), origin(o) {}
    double operator()(const coord_3d& r) const {
        const double x = r[0] - origin[0], y = r[1] - origin[1], z = r[2] - origin[2];
        return exp(-exponent * (x * x + y * y + z * z));
    }
};

std::vector<real_function_3d> make_set(World& world, const long n, const double e0,
                                       const double shift) {
    std::vector<real_function_3d> v;
    for (long i = 0; i < n; ++i) {
        const coord_3d o{shift * double(i), 0.3 * shift, -0.2 * shift};
        v.push_back(real_factory_3d(world)
                            .functor(shifted_gaussian(e0 + 0.4 * double(i), o))
                            .truncate_on_project());
    }
    truncate(world, v);
    return v;
}

}   // namespace

/// K with bra, ket and vf distinct, and nf != nocc, against the smallmem reference
int test_generic_exchange(World& world) {
    test_output t1("testing exchange with bra, ket and vf all distinct", world.rank() == 0);

    // nf != nocc on purpose: an implementation that used one length for both dimensions would pass
    // every existing test, since every consumer so far has had them equal
    const long nocc = 5, nf = 3;
    const auto bra = make_set(world, nocc, 1.0, 0.00);
    const auto ket = make_set(world, nocc, 1.3, 0.25);
    const auto vf = make_set(world, nf, 0.8, -0.35);

    t1.checkpoint(bra.size() == std::size_t(nocc) and ket.size() == std::size_t(nocc)
                          and vf.size() == std::size_t(nf),
                  "the operands have the intended, differing lengths");

    const double lo = 1.e-4;
    auto apply_with = [&](const Exchange<double, 3>::ExchangeAlgorithm alg) {
        Exchange<double, 3> K(world, lo, FunctionDefaults<3>::get_thresh());
        K.set_bra_and_ket(bra, ket);
        K.set_symmetric(false);
        K.set_algorithm(alg);
        return K(vf);
    };

    const auto reference = apply_with(Exchange<double, 3>::small_memory);
    t1.checkpoint(reference.size() == std::size_t(nf), "the reference returns one function per vf");

    const auto got = apply_with(Exchange<double, 3>::multiworld_efficient);
    t1.checkpoint(got.size() == std::size_t(nf), "the macrotask path returns one function per vf");

    // compare the functions themselves, not their norms: two different functions can share a norm
    double maxdiff = 0.0;
    for (std::size_t i = 0; i < reference.size() and i < got.size(); ++i)
        maxdiff = std::max(maxdiff, (reference[i] - got[i]).norm2());
    if (world.rank() == 0) print("largest |K_multiworld - K_smallmem| over vf:", maxdiff);
    t1.checkpoint(maxdiff < 10.0 * FunctionDefaults<3>::get_thresh(),
                  "the macrotask path agrees with the reference on every function");

    // a non-trivial result, so agreement cannot be two zeroes matching
    double smallest = std::numeric_limits<double>::max();
    for (const auto& f : reference) smallest = std::min(smallest, f.norm2());
    t1.checkpoint(smallest > 1.e-6, "the reference result is non-trivial");

    return t1.end();
}

int main(int argc, char** argv) {
    madness::initialize(argc, argv);
    int result = 0;
    {
        World& world = World::get_default();
        startup(world, argc, argv, true);
        FunctionDefaults<3>::set_thresh(1.e-5);
        FunctionDefaults<3>::set_k(7);
        FunctionDefaults<3>::set_cubic_cell(-16.0, 16.0);
        result += test_generic_exchange(world);
        world.gop.fence();
    }
    if (result == 0) print("\n --> all tests \033[32m passed \033[0m \n");
    else print("\n --> all tests \033[31m failed \033[0m \n");
    madness::finalize();
    return result;
}
