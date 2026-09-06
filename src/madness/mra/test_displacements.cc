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
    if (lattice_summed) {  // displacements that differ by a period along a lattice-summed axis are the same displacement
      const auto period = Translation(1) << n;
      std::vector<Key<NDIM>> canonical;
      for (const auto& disp : disps) {
        auto l = disp.translation();
        for (std::size_t d = 0; d != NDIM; ++d) l[d] = ((l[d] % period) + period) % period;
        canonical.emplace_back(n, l);
      }
      std::sort(canonical.begin(), canonical.end());
      const auto ndup = canonical.end() - std::unique(canonical.begin(), canonical.end());
      t.checkpoint(ndup == 0, what + ": all displacements distinct modulo the lattice (" + std::to_string(ndup) + " duplicates)");
    }
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

/// The iteration must not depend on whether a validator is given: a validator that accepts
/// everything must yield the same displacements as no validator at all.
int test_validator_agnostic(World& world) {
  test_output t("BoxSurfaceDisplacementRange: iteration does not depend on the presence of a validator", world.rank() == 0);

  constexpr std::size_t NDIM = 2;
  const auto disps_of = [](typename BoxSurfaceDisplacementRange<NDIM>::Validator validator) {
    Radii<NDIM> box_radius, surface_thickness;
    for (std::size_t d = 0; d != NDIM; ++d) {
      box_radius[d] = 1;
      surface_thickness[d] = 1;
    }
    BoxSurfaceDisplacementRange<NDIM> range(centered_key<NDIM>(4), box_radius, surface_thickness,
                                            array_of_bools<NDIM>{false}, std::move(validator));
    std::vector<Key<NDIM>> disps;
    for (auto&& disp : range) disps.push_back(disp);
    std::sort(disps.begin(), disps.end());
    return disps;
  };

  const auto without = disps_of({});
  const auto with = disps_of([](const Level, const typename BoxSurfaceDisplacementRange<NDIM>::PointPattern&,
                                std::optional<Key<NDIM>>&) { return true; });
  t.checkpoint(!with.empty(), "2D n=4 N={1,1}: surface is non-empty");
  t.checkpoint(without.size() == with.size(), "2D n=4 N={1,1}: same number of displacements with and without a validator (" +
                                                  std::to_string(without.size()) + " vs " + std::to_string(with.size()) + ")");
  t.checkpoint(without == with, "2D n=4 N={1,1}: same displacements with and without a validator");
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
      t.checkpoint(probe[d] == expected[d], "3D n=" + std::to_string(n) + " N=" + what + " component " + std::to_string(d));
    }
  };

  check_probe(level, {3, 2, 1}, {0, 0, radius_in_boxes(1, level)}, "{3,2,1}");
  check_probe(level, {1, 2, 3}, {radius_in_boxes(1, level), 0, 0}, "{1,2,3}");
  check_probe(level, {1, 2, 1}, {radius_in_boxes(1, level), 0, 0}, "{1,2,1}");
  check_probe(level, {3, 2, 3}, {radius_in_boxes(3, level), 0, 0}, "{3,2,3}");
  return t.end();
}

/// Test the degenerate cases in which criterion 2 is out of reach: a single dimension, and level 0.
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
  check_probe(level, std::array<std::int64_t, 1>{4}, {radius_in_boxes(4, level)}, "{4}");
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

/// Test the case in which every box radius is even, so that the face dimension alone cannot satisfy
/// criterion 2 and a half-simulation-cell offset along a second dimension is needed.
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
  check_probe(level, {4, 6, 4}, {radius_in_boxes(4, level), 0, radius_in_boxes(1, level)}, "{4,6,4}");

  return t.end();
}

/// Test differences between lattice-summed and non-lattice summed axes. In particular, offsets
/// are only needed for lattice-summed axes.
int test_lattice_summation_awareness(World& world) {
  test_output t("BoxSurfaceDisplacementRange: offset only where lattice summed", world.rank() == 0);

  const auto check_probe = [&]<std::size_t NDIM>(Level n, Radii<NDIM> box_radius,
                                                 array_of_bools<NDIM> is_lattice_summed,
                                                 std::array<std::int64_t, NDIM> expected,
                                                 const std::string& what) {
    Radii<NDIM> surface_thickness;
    for (std::size_t d = 0; d != NDIM; ++d) {
      if (box_radius[d]) surface_thickness[d] = 1;
    }
    BoxSurfaceDisplacementRange<NDIM> range(centered_key<NDIM>(n), box_radius, surface_thickness,
                                            is_lattice_summed);
    const auto probe = range.probing_displacement().translation();
    for (size_t d = 0; d < NDIM; d++) {
      t.checkpoint(probe[d] == expected[d], what + " component " + std::to_string(d));
    }
  };

  Level level = 4;
  const auto r = [&](std::int64_t N) { return radius_in_boxes(N, level); };

  // all even, nothing lattice summed => no offset, and the probe is the nearest surface point
  check_probe(level, Radii<3>{2, 2, 2}, array_of_bools<3>{false}, {r(2), 0, 0}, "3D N={2,2,2} not summed");
  check_probe(level, Radii<3>{4, 2, 6}, array_of_bools<3>{false}, {0, r(2), 0}, "3D N={4,2,6} not summed");
  // ... the same radii with lattice summation on do need the offset
  check_probe(level, Radii<3>{2, 2, 2}, array_of_bools<3>{true}, {r(2), r(1), 0}, "3D N={2,2,2} summed");

  // with only the first dimension lattice summed, that dimension is the best face even though it has to
  // pay for an offset: it folds to a half cell (8 boxes), whereas y or z would leave the probe a full
  // 2 half-cells out with nothing to fold it back
  check_probe(level, Radii<3>{2, 2, 2}, array_of_bools<3>{true, false, false}, {r(2), r(1), 0},
                 "3D N={2,2,2} summed along x only");
  // ... and all the more so when the unsummed dimensions are wider
  check_probe(level, Radii<3>{2, 4, 4}, array_of_bools<3>{true, false, false}, {r(2), r(1), 0},
                 "3D N={2,4,4} summed along x only");

  // an unrestricted dimension absorbs the offset, but again only when one is needed
  Radii<2> radii2d;
  radii2d[1] = 2;
  check_probe(level, radii2d, array_of_bools<2>{false, true}, {-r(1), r(2)}, "2D N={*,2} summed along y");
  check_probe(level, radii2d, array_of_bools<2>{false, false}, {0, r(2)}, "2D N={*,2} not summed");

  return t.end();
}

/// The offset for all-even radii is set by the caller.
int test_probe_offset_radius(World& world) {
  test_output t("BoxSurfaceDisplacementRange: offset magnitude follows the filter shortrange", world.rank() == 0);

  constexpr std::size_t NDIM = 3;
  const Level n = 6;                       // half-cell = 32 boxes, so the clamp is far away
  const auto check = [&](std::optional<Translation> cutoff, std::array<std::int64_t, NDIM> expected,
                         const std::string& what) {
    Radii<NDIM> box_radius, surface_thickness;
    for (std::size_t d = 0; d != NDIM; ++d) { box_radius[d] = 2; surface_thickness[d] = 1; }  // all even
    BoxSurfaceDisplacementRange<NDIM> range(centered_key<NDIM>(n), box_radius, surface_thickness,
                                            array_of_bools<NDIM>{true} /* lattice_summed */,
                                            typename BoxSurfaceDisplacementRange<NDIM>::Validator{}, cutoff);
    const auto probe = range.probing_displacement().translation();
    for (size_t d = 0; d < NDIM; d++)
      t.checkpoint(probe[d] == expected[d], what + " component " + std::to_string(d));
  };

  const auto face = radius_in_boxes(2, n);         // 64
  const auto half_cell = radius_in_boxes(1, n);    // 32

  // unknown shortrange cutoff => the maximal (half-cell) offset, i.e. the behavior when nothing better is known
  check(std::nullopt, {face, half_cell, 0}, "cutoff unknown");
  // a cutoff inside the half cell is used verbatim -- this is the case that matters in practice, since
  // the cutoff is capped by bmax_default() (4 in 3D) while the half cell grows as 2^{n-1}
  check(Translation(5), {face, 5, 0}, "cutoff=5");
  check(Translation(2), {face, 2, 0}, "cutoff=2");
  // a cutoff beyond the half cell is clamped: a step of 2^{n-1}+k folds back down to 2^{n-1}-k
  check(Translation(1000), {face, half_cell, 0}, "cutoff beyond half cell");
  // a cutoff of zero says nothing is filtered, so the surface reaches the source and no offset helps;
  // fall back to the bare face probe, whose norm is the on-site norm and screens nothing
  check(Translation(0), {face, 0, 0}, "cutoff=0 (nothing filtered)");

  // the same clamping applies when the offset lands on an unrestricted dimension
  {
    Radii<2> box_radius, surface_thickness;
    box_radius[1] = 2; surface_thickness[1] = 1;
    BoxSurfaceDisplacementRange<2> range(centered_key<2>(n), box_radius, surface_thickness,
                                         array_of_bools<2>{true},
                                         typename BoxSurfaceDisplacementRange<2>::Validator{}, Translation(3));
    const auto probe = range.probing_displacement().translation();
    t.checkpoint(probe[0] == -3, "2D N={*,2} cutoff=3 component 0");
    t.checkpoint(probe[1] == radius_in_boxes(2, n), "2D N={*,2} cutoff=3 component 1");
  }

  return t.end();
}

}

int main(int argc, char** argv) {
  World& world = madness::initialize(argc, argv);
  startup(world, argc, argv, true);

  int errors = 0;
  errors += test_no_duplicates(world);
  errors += test_validator_agnostic(world);
  errors += test_odds(world);
  errors += test_singular(world);
  errors += test_mixed(world);
  errors += test_all_evens(world);
  errors += test_lattice_summation_awareness(world);
  errors += test_probe_offset_radius(world);

  world.gop.fence();
  madness::finalize();
  return errors;
}