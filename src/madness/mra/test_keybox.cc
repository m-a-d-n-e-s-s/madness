/// \file test_keybox.cc
/// \brief Tests for Key::is_neighbor_of(simpt, bperiodic) — the geometric
///        replacement for the legacy
///        `simpt2key(simpt, n).is_neighbor_of(*this, bperiodic)` idiom.

#include <madness/mra/mra.h>
#include <madness/mra/funcdefaults.h>
#include <madness/mra/key.h>
#include <madness/world/test_utilities.h>

using namespace madness;

namespace {

template <std::size_t NDIM>
Key<NDIM> make_key(Level n, std::initializer_list<Translation> ls) {
    Vector<Translation, NDIM> v;
    std::size_t i = 0;
    for (auto x : ls) { v[i++] = x; }
    return Key<NDIM>(n, v);
}

// 1D basic adjacency: at level 2 the 4 boxes are
// l=0:[0,.25] l=1:[.25,.5] l=2:[.5,.75] l=3:[.75,1].
// is_neighbor_of(simpt) is true on this box and the two boxes adjacent to
// the box(es) containing simpt.
int test_1d_basic(World& world) {
    test_output t("Key<1>::is_neighbor_of(simpt) — 1D");

    array_of_bools<1> nonper{false};

    auto key0 = make_key<1>(2, {0});
    auto key1 = make_key<1>(2, {1});
    auto key2 = make_key<1>(2, {2});
    auto key3 = make_key<1>(2, {3});

    // simpt at sim=0.30 (interior of l=1)
    Vector<double,1> p_interior{0.30};
    t.checkpoint( key0.is_neighbor_of(p_interior, nonper),
                  "interior pt: l=0 is neighbor of containing l=1");
    t.checkpoint( key1.is_neighbor_of(p_interior, nonper),
                  "interior pt: l=1 contains pt");
    t.checkpoint( key2.is_neighbor_of(p_interior, nonper),
                  "interior pt: l=2 is neighbor of containing l=1");
    t.checkpoint(!key3.is_neighbor_of(p_interior, nonper),
                  "interior pt: l=3 is NOT neighbor of containing l=1");

    // simpt at sim=0.5 lies on wall between l=1 and l=2; both adjacent boxes
    // and their neighbors should be flagged.
    Vector<double,1> p_wall{0.5};
    t.checkpoint( key0.is_neighbor_of(p_wall, nonper),
                  "wall pt: l=0 (neighbor of l=1) flagged");
    t.checkpoint( key1.is_neighbor_of(p_wall, nonper),
                  "wall pt: l=1 (contains pt) flagged");
    t.checkpoint( key2.is_neighbor_of(p_wall, nonper),
                  "wall pt: l=2 (contains pt) flagged");
    t.checkpoint( key3.is_neighbor_of(p_wall, nonper),
                  "wall pt: l=3 (neighbor of l=2) flagged");

    return t.end();
}

// Symmetry across a wall: the set of flagged boxes should be the same for
// simpt = 0.5, 0.5 + ε, and 0.5 − ε. The legacy code (simpt2key + is_neighbor_of)
// flips which box `specialkey` lands in and shifts the 3-box window.
int test_symmetry_at_boundary(World& world) {
    test_output t("Key::is_neighbor_of(simpt) — boundary symmetry");

    array_of_bools<1> nonper{false};
    constexpr double eps = 1e-15;
    constexpr Level n = 4;
    const Translation NB = Translation(1) << n;  // 16 boxes

    for (Translation li = 0; li < NB; ++li) {
        auto key = make_key<1>(n, {li});
        bool hit_below = key.is_neighbor_of(Vector<double,1>{0.5 - eps}, nonper);
        bool hit_above = key.is_neighbor_of(Vector<double,1>{0.5 + eps}, nonper);
        bool hit_exact = key.is_neighbor_of(Vector<double,1>{0.5}, nonper);
        // wall sits between l=7 and l=8; expected neighborhood is {6,7,8,9}
        bool expected = (li >= 6 && li <= 9);
        t.checkpoint(hit_exact == expected,
                     "exact wall point: l=" + std::to_string(li));
        t.checkpoint(hit_below == expected,
                     "wall - eps: l=" + std::to_string(li));
        t.checkpoint(hit_above == expected,
                     "wall + eps: l=" + std::to_string(li));
    }
    return t.end();
}

// 3D corner-of-cell: point at sim coord (0.5, 0.5, 0.5) touches 8 boxes at
// every level n>=1. The 1-box-width neighborhood of those 8 contains 4^3=64
// boxes: {l[i] ∈ {MID-2, MID-1, MID, MID+1}}.
int test_3d_corner(World& world) {
    test_output t("Key<3>::is_neighbor_of(simpt) — 3D corner");

    array_of_bools<3> nonper{false, false, false};
    constexpr Level n = 3;
    const Translation MID = Translation(1) << (n - 1);  // = 4
    const Vector<double,3> origin{0.5, 0.5, 0.5};

    int hits = 0;
    const Translation NB = Translation(1) << n;
    for (Translation i = 0; i < NB; ++i)
      for (Translation j = 0; j < NB; ++j)
        for (Translation k = 0; k < NB; ++k) {
            auto key = make_key<3>(n, {i, j, k});
            if (key.is_neighbor_of(origin, nonper)) ++hits;
        }
    t.checkpoint(hits == 64,
                 "64 boxes within one box-width of the corner (4^NDIM)");

    // explicit symmetry across the corner
    auto key_lo = make_key<3>(n, {MID - 2, MID - 2, MID - 2});
    auto key_hi = make_key<3>(n, {MID + 1, MID + 1, MID + 1});
    t.checkpoint(key_lo.is_neighbor_of(origin, nonper) ==
                 key_hi.is_neighbor_of(origin, nonper),
                 "neighborhood is symmetric across the corner");

    return t.end();
}

// Periodic seam: point at sim coord 0.0 (== 1.0 mod 1) — boxes at both ends
// of the domain should be flagged when bperiodic[i] is true.
int test_periodic_seam(World& world) {
    test_output t("Key<1>::is_neighbor_of(simpt) — periodic seam");

    array_of_bools<1> per{true};
    constexpr Level n = 3;  // 8 boxes
    auto key_lo = make_key<1>(n, {0});  // [0/8, 1/8]
    auto key_hi = make_key<1>(n, {7});  // [7/8, 8/8] — wraps to 0/8

    t.checkpoint( key_lo.is_neighbor_of(Vector<double,1>{0.0}, per),
                  "sim=0 contained in low box (periodic)");
    t.checkpoint( key_hi.is_neighbor_of(Vector<double,1>{0.0}, per),
                  "sim=0 wraps into high box (periodic)");
    t.checkpoint( key_lo.is_neighbor_of(Vector<double,1>{1.0}, per),
                  "sim=1 wraps into low box (periodic)");
    t.checkpoint( key_hi.is_neighbor_of(Vector<double,1>{1.0}, per),
                  "sim=1 contained in high box (periodic)");

    // Non-periodic: no wrap; sim=1 is contained only in (and adjacent to) high boxes
    array_of_bools<1> nonper{false};
    auto key_mid = make_key<1>(n, {3});  // [3/8, 4/8] — far from sim=1
    t.checkpoint(!key_mid.is_neighbor_of(Vector<double,1>{1.0}, nonper),
                  "sim=1 does NOT reach interior box (non-periodic)");
    t.checkpoint( key_hi.is_neighbor_of(Vector<double,1>{1.0}, nonper),
                  "sim=1 hits high box (non-periodic)");

    return t.end();
}

// Match interior behavior of the legacy code: for points strictly inside a
// box, key.is_neighbor_of(simpt) and the legacy
// simpt2key(simpt).is_neighbor_of(key) idiom must agree.
int test_legacy_interior_equivalence(World& world) {
    test_output t("Key::is_neighbor_of(simpt) ≡ legacy on interior");

    array_of_bools<2> nonper{false, false};
    constexpr Level n = 4;
    std::vector<Vector<double,2>> pts = {
        {0.123, 0.456},
        {0.781, 0.219},
        {0.9999, 0.0001}};

    const Translation NB = Translation(1) << n;
    for (auto& pt : pts) {
        Vector<double,2> simpt = pt;
        Key<2> specialkey = simpt2key(simpt, n);
        for (Translation i = 0; i < NB; ++i)
          for (Translation j = 0; j < NB; ++j) {
              auto key = make_key<2>(n, {i, j});
              bool legacy = specialkey.is_neighbor_of(key, nonper);
              bool now = key.is_neighbor_of(simpt, nonper);
              t.checkpoint(legacy == now,
                           "legacy ≡ is_neighbor_of(simpt) at interior pt");
          }
    }
    return t.end();
}

}  // namespace

int main(int argc, char** argv) {
    World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);

    int errors = 0;
    errors += test_1d_basic(world);
    errors += test_symmetry_at_boundary(world);
    errors += test_3d_corner(world);
    errors += test_periodic_seam(world);
    errors += test_legacy_interior_equivalence(world);

    world.gop.fence();
    madness::finalize();
    return errors;
}
