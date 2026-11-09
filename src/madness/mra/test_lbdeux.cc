/// Unit tests for the LBDeuxPmap process map.
///
/// Each map is written by hand, so that every case reaches one lookup branch or one constructor
/// check. The ranks are only labels here, checked against nproc, not against the size of the World.

#include <madness/mra/mra.h>
#include <madness/mra/lbdeux.h>
#include <madness/world/test_utilities.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace madness;

template <std::size_t NDIM>
using Assignment = std::vector<std::pair<Key<NDIM>, ProcessID>>;

/// Rank bound for the hand-written maps, which use the ranks 0, 1 and 2
static constexpr ProcessID nproc = 3;

static Key<1> key1(Level n, Translation l) {
    return {n, vec(l)};
}

static Key<3> key3(Level n, Translation x, Translation y, Translation z) {
    return {n, vec(x, y, z)};
}

/// True only for MadnessException. An old-style `throw "..."` or a clean return is a failure.
template <std::size_t NDIM>
static bool constructor_throws(const Assignment<NDIM>& assignment) {
    try {
        LBDeuxPmap<NDIM> pmap(assignment, nproc);
    } catch (const MadnessException&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

/// Mapped tree, with a distinct rank for each entry:
///
///     (0,0) -> 0
///       (1,1) -> 1
///         (2,2) -> 2
static int test_nearest_ancestor(World& world) {
    test_output t("owner() returns the nearest mapped ancestor", world.rank() == 0);
    const LBDeuxPmap<1> pmap(Assignment<1>{{Key<1>(0), 0}, {key1(1, 1), 1}, {key1(2, 2), 2}}, nproc);

    t.checkpoint(pmap.owner(Key<1>(0)), 0, "exact key at the root");
    t.checkpoint(pmap.owner(key1(1, 1)), 1, "exact key at an interior level");
    t.checkpoint(pmap.owner(key1(2, 2)), 2, "exact key at the deepest level");
    t.checkpoint(pmap.owner(key1(3, 4)), 2, "child of the deepest entry");
    t.checkpoint(pmap.owner(key1(4, 9)), 2, "grandchild of the deepest entry");
    t.checkpoint(pmap.owner(key1(2, 3)), 1, "mapped sibling is not an ancestor");
    t.checkpoint(pmap.owner(key1(3, 1)), 0, "grandchild of an unmapped key goes to the root");
    return t.end();
}

/// Each level-20 key below has the translation of its level-3 ancestor in the high 3 bits. The
/// two level-3 entries are siblings, so a lookup that stops one level too high finds neither.
static int test_key_below_deepest_level(World& world) {
    test_output t("owner() resolves a key below the deepest entry", world.rank() == 0);
    const LBDeuxPmap<3> pmap(Assignment<3>{{Key<3>(0), 0}, {key3(3, 5, 2, 7), 1}, {key3(3, 5, 2, 6), 2}}, nproc);

    t.checkpoint(pmap.owner(key3(20, 667705, 262145, 1048575)), 1, "level 20 inside (3,[5,2,7])");
    t.checkpoint(pmap.owner(key3(20, 667705, 262145, 917503)), 2, "level 20 inside (3,[5,2,6])");
    t.checkpoint(pmap.owner(key3(20, 0, 0, 0)), 0, "level 20 outside both subtrees");
    return t.end();
}

/// A level-0 key with a nonzero translation is not in the box, and only such a key can reach this
/// branch, because the constructor requires the root.
static int test_unmatched_level_zero_key_throws(World& world) {
    test_output t("owner() throws if no ancestor matches", world.rank() == 0);
    const LBDeuxPmap<1> pmap(Assignment<1>{{Key<1>(0), 0}}, nproc);
    bool threw = false;
    try {
        pmap.owner(key1(0, 1));
    } catch (const MadnessException&) {
        threw = true;
    } catch (...) {
    }
    t.checkpoint(threw, "level-0 key with a nonzero translation");
    return t.end();
}

static int test_missing_root_throws(World& world) {
    test_output t("constructor rejects a map without the root", world.rank() == 0);
    t.checkpoint(constructor_throws<1>({{key1(1, 0), 0}, {key1(1, 1), 1}}), "only level-1 keys");
    t.checkpoint(constructor_throws<1>({}), "empty assignment");
    return t.end();
}

static int test_duplicate_key_throws(World& world) {
    test_output t("constructor rejects a duplicate key", world.rank() == 0);
    t.checkpoint(constructor_throws<1>({{Key<1>(0), 0}, {key1(1, 1), 1}, {key1(1, 1), 1}}), "same rank");
    t.checkpoint(constructor_throws<1>({{Key<1>(0), 0}, {key1(1, 1), 1}, {key1(1, 1), 2}}), "different rank");
    return t.end();
}

static int test_rank_out_of_range_throws(World& world) {
    test_output t("constructor rejects a rank out of range", world.rank() == 0);
    t.checkpoint(constructor_throws<1>({{Key<1>(0), 0}, {key1(1, 1), -1}}), "negative rank");
    t.checkpoint(constructor_throws<1>({{Key<1>(0), 0}, {key1(1, 1), nproc}}), "rank equal to nproc");
    t.checkpoint(!constructor_throws<1>({{Key<1>(0), 0}, {key1(1, 1), nproc - 1}}), "rank nproc-1 is accepted");
    return t.end();
}

static double gaussian(const coord_3d& r) {
    return exp(-(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]));
}

/// A negative cost at the root stands in for roundoff in its residual. The real residual can only
/// go below zero after all children split off, which needs a large fac in load_balance().
struct NegativeRootCost {
    static constexpr double root_cost = -1.e-3;
    double operator()(const Key<3>& key, const FunctionNode<double, 3>&) const {
        return (key.level() == 0) ? root_cost : 1.0;
    }
};

/// A guard, not a regression test: a zero root cost passed the old `cost >= 0` filter too.
struct ZeroCost {
    double operator()(const Key<3>&, const FunctionNode<double, 3>&) const {
        return 0.0;
    }
};

struct NaNCost {
    double operator()(const Key<3>& key, const FunctionNode<double, 3>&) const {
        return (key.level() == 2) ? std::nan("") : 1.0;
    }
};

/// A non-finite cost must stop load_balance() on every rank, before it reaches the sort.
static int test_driver_rejects_nonfinite_cost(World& world) {
    test_output t("load_balance() rejects a non-finite cost", world.rank() == 0);
    const real_function_3d f = real_factory_3d(world).f(gaussian);
    LoadBalanceDeux<3> lb(world);
    lb.add_tree(f, NaNCost(), true);
    bool threw = false;
    try {
        lb.load_balance();
    } catch (const MadnessException&) {
        threw = true;
    }
    t.checkpoint(threw, "NaN cost below the root");
    return t.end();
}

/// True if load_balance(fac) splits every child off the root, so that only the root's own cost
/// stays as its residual. partition() tests the remaining cost against 1.1 times the average before
/// it subtracts each child, in descending order, so the last test sees the smallest child.
template <typename costT>
static bool all_children_split_off(World& world, const real_function_3d& f, const costT& cost,
                                   const double root_cost, const double fac) {
    std::vector<double> child_cost(8, 0.0);
    const auto& coeffs = f.get_impl()->get_coeffs();
    for (auto it = coeffs.begin(); it != coeffs.end(); ++it) {
        const Key<3>& key = it->first;
        if (key.level() == 0) continue;
        const Key<3> child = key.parent(key.level() - 1);
        child_cost[4 * child.translation()[0] + 2 * child.translation()[1] + child.translation()[2]] +=
            cost(key, it->second);
    }
    world.gop.sum(child_cost.data(), child_cost.size());
    double total = root_cost;
    for (const double c : child_cost) total += c;
    const double threshold = 1.1 * total / (world.size() * fac);
    const double smallest = *std::min_element(child_cost.begin(), child_cost.end());
    return smallest + root_cost > threshold;
}

/// Runs the real driver, then asks the returned map for the root and for every local function key.
template <typename costT>
static int test_driver_keeps_root(World& world, const std::string& name, const costT& cost, const double fac) {
    test_output t(name, world.rank() == 0);
    const real_function_3d f = real_factory_3d(world).f(gaussian);

    if constexpr (std::is_same_v<costT, NegativeRootCost>) {
        // Without this, the root residual stays positive and the test passes without the fix.
        t.checkpoint(all_children_split_off(world, f, cost, NegativeRootCost::root_cost, fac),
                     "precondition: the root residual goes negative");
    }

    LoadBalanceDeux<3> lb(world);
    lb.add_tree(f, cost, true);

    std::shared_ptr<WorldDCPmapInterface<Key<3>>> pmap;
    try {
        pmap = lb.load_balance(fac);
    } catch (const MadnessException&) {
        // The constructor sees the same broadcast map on every rank, so every rank comes here.
        t.checkpoint(false, "load_balance() returns a map");
        return t.end();
    }
    t.checkpoint(true, "load_balance() returns a map");

    const auto valid_rank = [&world](const ProcessID p) { return p >= 0 && p < world.size(); };
    t.checkpoint(valid_rank(pmap->owner(Key<3>(0))), "the map owns the root");

    bool all_keys_owned = true;
    const auto& coeffs = f.get_impl()->get_coeffs();
    for (auto it = coeffs.begin(); it != coeffs.end(); ++it) {
        all_keys_owned = all_keys_owned && valid_rank(pmap->owner(it->first));
    }
    t.checkpoint(all_keys_owned, "the map owns every local function key");
    return t.end();
}

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    startup(world, argc, argv, true);

    FunctionDefaults<3>::set_k(4);
    FunctionDefaults<3>::set_thresh(1.e-3);
    FunctionDefaults<3>::set_cubic_cell(-4.0, 4.0);

    int failures = 0;
    failures += test_nearest_ancestor(world);
    failures += test_key_below_deepest_level(world);
    failures += test_missing_root_throws(world);
    failures += test_duplicate_key_throws(world);
    failures += test_rank_out_of_range_throws(world);
    failures += test_unmatched_level_zero_key_throws(world);
    failures += test_driver_keeps_root(world, "load_balance() keeps a root with negative residual", NegativeRootCost(), 20.0);
    failures += test_driver_keeps_root(world, "load_balance() keeps a root with zero cost", ZeroCost(), 1.0);
    failures += test_driver_rejects_nonfinite_cost(world);

    world.gop.fence();
    finalize();
    return failures;
}
