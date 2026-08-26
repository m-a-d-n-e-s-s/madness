#include <madness/mra/mra.h>
#include <madness/mra/displacements.h>
#include <madness/world/test_utilities.h>

using namespace madness;

namespace {

template <std::size_t NDIM>
using Radii = std::array<std::optional<std::int64_t>, NDIM>;

/// the range boundary, in boxes, for a kernel range of \p N half-simulation-cells at level \p n
/// mirrors the conversion in the BoxSurfaceDisplacementRange constructor
Translation radius_in_boxes(std::int64_t N, Level n) {
  return (n == 0) ? (N + 1) / 2 : (N * Translation(1) << (n - 1));
}

/// the center box of the level-\p n grid
template <std::size_t NDIM>
Key<NDIM> centered_key(Level n) {
  return Key<NDIM>(n, Vector<Translation, NDIM>(n == 0 ? 0 : (Translation(1) << (n - 1))));
}

/// keeps only destinations inside the simulation cell
template <std::size_t NDIM>
typename BoxSurfaceDisplacementRange<NDIM>::Validator in_domain_only() {
  return [](const Level level, const typename BoxSurfaceDisplacementRange<NDIM>::PointPattern& dest,
            std::optional<Key<NDIM>>&) -> bool {
    const auto twon = (Translation(1) << level);
    for (std::size_t d = 0; d != NDIM; ++d)
      if (dest[d].has_value() && (*dest[d] < 0 || *dest[d] >= twon)) return false;
    return true;
  };
}

/// No displacement may be enumerated twice.
/// We test any and cases where the logic splits, including lattice summed or not,
/// and odd vs even N
int test_no_duplicates(World& world) {
  test_output t("BoxSurfaceDisplacementRange: no duplicate displacements", world.rank() == 0);

  const auto check_unique = [&](auto ndim_tag, Level n, std::array<std::int64_t, decltype(ndim_tag)::value> N,
                                bool lattice_summed, bool filter_to_domain, const std::string& what) {
    constexpr std::size_t NDIM = decltype(ndim_tag)::value;
    Radii<NDIM> box_radius, surface_thickness;
    for (std::size_t d = 0; d != NDIM; ++d) {
      box_radius[d] = N[d];
      surface_thickness[d] = 1;
    }
    // N.B. the domain filter is only meaningful when the range boundary can land
    // inside the cell at all.
    BoxSurfaceDisplacementRange<NDIM> range(
        centered_key<NDIM>(n), box_radius, surface_thickness, array_of_bools<NDIM>{lattice_summed},
        filter_to_domain ? in_domain_only<NDIM>() : typename BoxSurfaceDisplacementRange<NDIM>::Validator{});
    std::vector<Key<NDIM>> disps;
    for (auto&& disp : range) disps.push_back(disp);
    std::sort(disps.begin(), disps.end());
    const auto last = std::unique(disps.begin(), disps.end());
    t.checkpoint(!disps.empty(), what + ": surface is non-empty");
    t.checkpoint(last == disps.end(), what + ": all displacements distinct");
  };

  constexpr auto d2 = std::integral_constant<std::size_t, 2>{};
  constexpr auto d3 = std::integral_constant<std::size_t, 3>{};

  // odd N and its even counterpart
  check_unique(d2, 4, {1, 1}, false, true, "2D n=4 N={1,1} plain, in-domain");
  check_unique(d2, 4, {1, 1}, false, false, "2D n=4 N={1,1} plain");
  check_unique(d2, 4, {2, 2}, false, false, "2D n=4 N={2,2} plain");
  // ... and both parities again with lattice summation on
  check_unique(d2, 4, {1, 1}, true, false, "2D n=4 N={1,1} lattice-summed");
  check_unique(d2, 4, {2, 2}, true, false, "2D n=4 N={2,2} lattice-summed");
  check_unique(d2, 4, {2, 3}, true, false, "2D n=4 N={2,3} lattice-summed (mixed parity)");
  // ... and in 3D, where each face has two free axes rather than one
  check_unique(d3, 3, {1, 1, 1}, true, false, "3D n=3 N={1,1,1} lattice-summed");
  check_unique(d3, 3, {2, 2, 2}, true, false, "3D n=3 N={2,2,2} lattice-summed");
  check_unique(d3, 3, {1, 2, 3}, true, false, "3D n=3 N={1,2,3} lattice-summed (mixed parity)");

  return t.end();
}

/// Test the cases where there is an odd box radius.
int test_odds(World& world) {
  test_output t("BoxSurfaceDisplacementRange: odd dimensions", world.rank() == 0);

  Level level = 4;

  const auto check_probe = [&](Level n, std::array<std::int64_t, 3> N, std::array<std::int64_t, 3> expected,
                               const std::string& what) {
    constexpr std::size_t NDIM = 3;
    Radii<NDIM> box_radius, surface_thickness;
    for (std::size_t d = 0; d != NDIM; ++d) {
      box_radius[d] = N[d];
      surface_thickness[d] = 1;
    }
    BoxSurfaceDisplacementRange<NDIM> range(centered_key<NDIM>(n), box_radius, surface_thickness,
                                            array_of_bools<NDIM>{true} /* lattice_summed */);
    const auto probe = range.probing_displacement().translation();
    for (size_t d = 0; d < NDIM; d++) {
      t.checkpoint(probe[d] == expected[d], "3D n=4 N=" + what + " component " + std::to_string(d));
    }
  };

  check_probe(level, {3, 2, 1}, {0, 0, radius_in_boxes(1, level)}, "{3,2,1}");
  check_probe(level, {1, 2, 3}, {radius_in_boxes(1, level), 0, 0}, "{1,2,3}");
  check_probe(level, {1, 2, 1}, {radius_in_boxes(1, level), 0, 0}, "{1,2,1}");
  check_probe(level, {3, 2, 3}, {radius_in_boxes(3, level), 0, 0}, "{3,2,3}");
  return t.end();
}

/// Test the cases where there is an odd box radius.
int test_singular(World& world) {
  test_output t("BoxSurfaceDisplacementRange: singular cases", world.rank() == 0);

  const auto check_probe = [&]<std::size_t NDIM>(Level n, std::array<std::int64_t, NDIM> N, std::array<std::int64_t, NDIM> expected,
                               const std::string& what) {
    Radii<NDIM> box_radius, surface_thickness;
    for (std::size_t d = 0; d != NDIM; ++d) {
      box_radius[d] = N[d];
      surface_thickness[d] = 1;
    }
    BoxSurfaceDisplacementRange<NDIM> range(centered_key<NDIM>(n), box_radius, surface_thickness,
                                            array_of_bools<NDIM>{true} /* lattice_summed */);
    const auto probe = range.probing_displacement().translation();
    for (size_t d = 0; d < NDIM; d++) {
      t.checkpoint(probe[d] == expected[d], std::to_string(NDIM) + "D n=" + std::to_string(n) + " N=" + what + " component " + std::to_string(d));
    }
  };

  Level level = 4;
  check_probe(level, std::array<std::int64_t, 1>{4}, {radius_in_boxes(4, level)}, "{2}");
  check_probe(0, std::array<std::int64_t, 3>{3, 2, 1}, {0, 0, 1}, "{3,2,1}");
  return t.end();
}

/// Test mixed boundary conditions
int test_mixed(World& world) {
  test_output t("BoxSurfaceDisplacementRange: mixed cases", world.rank() == 0);

  Level level = 4;

  const auto check_probe = [&]<std::size_t NDIM>(Level n, Radii<NDIM> box_radius, std::array<std::int64_t, NDIM> expected,
                               const std::string& what) {
    Radii<NDIM> surface_thickness;
    for (std::size_t d = 0; d != NDIM; ++d) {
      if (box_radius[d]) surface_thickness[d] = 1;
    }
    BoxSurfaceDisplacementRange<NDIM> range(centered_key<NDIM>(n), box_radius, surface_thickness,
                                            array_of_bools<NDIM>{true} /* lattice_summed */);
    const auto probe = range.probing_displacement().translation();
    for (size_t d = 0; d < NDIM; d++) {
      t.checkpoint(probe[d] == expected[d], what + " component " + std::to_string(d));
    }
  };


  Radii<2> radii2d;
  radii2d[1] = 2;
  check_probe(level, radii2d, {-radius_in_boxes(1, level), radius_in_boxes(2, level)}, "2D n=4 N={*,2}");
  radii2d[1] = 3;
  check_probe(level, radii2d, {0, radius_in_boxes(3, level)}, "2D n=4 N={*,3}");
  Radii<3> radii3d;
  radii3d[1] = 2;
  radii3d[2] = 3;
  check_probe(level, radii3d, {0, 0, radius_in_boxes(3, level)}, "3D n=4 N={*,2,3}");
  return t.end();
}

/// Test the cases where there is an odd box radius.
int test_all_evens(World& world) {
  test_output t("BoxSurfaceDisplacementRange: pure evens", world.rank() == 0);

  Level level = 4;

  const auto check_probe = [&](Level n, std::array<std::int64_t, 3> N, std::array<std::int64_t, 3> expected,
                               const std::string& what) {
    constexpr std::size_t NDIM = 3;
    Radii<NDIM> box_radius, surface_thickness;
    for (std::size_t d = 0; d != NDIM; ++d) {
      box_radius[d] = N[d];
      surface_thickness[d] = 1;
    }
    BoxSurfaceDisplacementRange<NDIM> range(centered_key<NDIM>(n), box_radius, surface_thickness,
                                            array_of_bools<NDIM>{true} /* lattice_summed */);
    const auto probe = range.probing_displacement().translation();
    for (size_t d = 0; d < NDIM; d++) {
      t.checkpoint(probe[d] == expected[d], std::to_string(NDIM) + "D n=" + std::to_string(n) + " N=" + what + " component " + std::to_string(d));
    }
  };

  check_probe(level, {4, 6, 2}, {radius_in_boxes(1, level), 0, radius_in_boxes(2, level)}, "{4,6,2}");
  check_probe(level, {6, 2, 2}, {0, radius_in_boxes(2, level), radius_in_boxes(1, level)}, "{6,2,2}");
  check_probe(level, {8, 2, 2}, {0, radius_in_boxes(2, level), radius_in_boxes(1, level)}, "{8,2,2}");
  check_probe(level, {6, 4, 2}, {0, radius_in_boxes(1, level), radius_in_boxes(2, level)}, "{6,4,2}");
  check_probe(level, {4, 6, 4}, {radius_in_boxes(4, level), 0, radius_in_boxes(1, level)}, "{4,6,2}");

  return t.end();
}

}

int main(int argc, char** argv) {
  World& world = madness::initialize(argc, argv);
  startup(world, argc, argv, true);

  int errors = 0;
  errors += test_no_duplicates(world);
  errors += test_odds(world);
  errors += test_singular(world);
  errors += test_mixed(world);
  errors += test_all_evens(world);

  world.gop.fence();
  madness::finalize();
  return errors;
}