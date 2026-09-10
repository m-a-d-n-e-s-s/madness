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

/// where the probe of a face sits: on its innermost layer, i.e. one box (the surface thickness used throughout) inside the boundary
Translation probe_radius(std::int64_t N, Level n) {
  return radius_in_boxes(N, n) - 1;
}

/// the center box of the level-\p n grid
template <std::size_t NDIM>
Key<NDIM> centered_key(Level n) {
  return Key<NDIM>(n, Vector<Translation, NDIM>(n == 0 ? 0 : (Translation(1) << (n - 1))));
}

template <std::size_t NDIM>
using Reach = StandardDisplacementsReach<NDIM>;
template <std::size_t NDIM>
using Validator = BoxSurfaceDisplacementValidator<NDIM>;

/// standard displacements of unit-width cells reached out to real distance squared \p max_distsq
template <std::size_t NDIM>
Reach<NDIM> unit_reach(double max_distsq) {
  Reach<NDIM> reach{max_distsq, {}};
  reach.cell_width.fill(1.);
  return reach;
}

/// standard displacements that reached (far) beyond bmax boxes, so that the least unfiltered
/// offset is bmax+1 boxes, capped at a half cell
template <std::size_t NDIM>
Reach<NDIM> far_reach() {
  return unit_reach<NDIM>(1e6);
}
template <std::size_t NDIM>
Translation far_offset(Level n) {
  return std::min<Translation>(Displacements<NDIM>::bmax_default() + 1, Translation(1) << (n - 1));
}

/// keeps only destinations inside the simulation cell (finite domain), nothing else is filtered
template <std::size_t NDIM>
Validator<NDIM> in_domain_only(const array_of_bools<NDIM>& is_lattice_summed) {
  return Validator<NDIM>(array_of_bools<NDIM>{false}, is_lattice_summed);
}

/// keeps everything (infinite domain, no standard displacements to deduplicate against)
template <std::size_t NDIM>
Validator<NDIM> keep_all(const array_of_bools<NDIM>& is_lattice_summed) {
  return Validator<NDIM>(array_of_bools<NDIM>{true}, is_lattice_summed);
}

/// surface with unit thickness; if only `reach` is given the validator keeps everything but the standard displacements
template <std::size_t NDIM>
BoxSurfaceDisplacementRange<NDIM> make_range(Level n, const Radii<NDIM>& box_radius,
                                             const array_of_bools<NDIM>& is_lattice_summed,
                                             std::optional<Reach<NDIM>> reach = {},
                                             std::optional<Validator<NDIM>> validator = {},
                                             std::optional<Key<NDIM>> center = {}) {
  Radii<NDIM> surface_thickness;
  for (std::size_t d = 0; d != NDIM; ++d) {
    if (box_radius[d]) surface_thickness[d] = 1;
  }
  if (reach && !validator) validator.emplace(array_of_bools<NDIM>{true}, is_lattice_summed, std::move(reach));
  return BoxSurfaceDisplacementRange<NDIM>(center.value_or(centered_key<NDIM>(n)), box_radius, surface_thickness,
                                           is_lattice_summed, std::move(validator));
}

/// the probing displacement of every face, as a translation; null for dimensions of unlimited size (no face)
template <std::size_t NDIM>
using Probes = std::array<std::optional<std::array<std::int64_t, NDIM>>, NDIM>;

template <std::size_t NDIM>
Probes<NDIM> probes_of(const BoxSurfaceDisplacementRange<NDIM>& range) {
  Probes<NDIM> result;
  for (std::size_t d = 0; d != NDIM; ++d) {
    if (!range.probing_displacements()[d]) continue;
    std::array<std::int64_t, NDIM> v;
    for (std::size_t e = 0; e != NDIM; ++e) v[e] = range.probing_displacement(d).translation()[e];
    result[d] = v;
  }
  return result;
}

template <std::size_t NDIM>
Probes<NDIM> probes_of(Level n, const Radii<NDIM>& box_radius, const array_of_bools<NDIM>& is_lattice_summed,
                       std::optional<Reach<NDIM>> reach = {}) {
  return probes_of(make_range<NDIM>(n, box_radius, is_lattice_summed, std::move(reach)));
}

template <std::size_t NDIM>
void check_probes(test_output& t, const Probes<NDIM>& actual, const Probes<NDIM>& expected, const std::string& what) {
  for (std::size_t d = 0; d != NDIM; ++d) {
    const auto face = what + " face " + std::to_string(d);
    t.checkpoint(actual[d].has_value() == expected[d].has_value(), face + " exists iff radius is finite");
    if (actual[d] && expected[d]) {
      for (std::size_t e = 0; e != NDIM; ++e)
        t.checkpoint((*actual[d])[e] == (*expected[d])[e], face + " component " + std::to_string(e));
    }
  }
}

/// shorthand for all-finite radii
template <std::size_t NDIM>
Radii<NDIM> finite(std::array<std::int64_t, NDIM> N) {
  Radii<NDIM> result;
  for (std::size_t d = 0; d != NDIM; ++d) result[d] = N[d];
  return result;
}

/// No displacement may be enumerated twice.
/// We test any and cases where the logic splits, including lattice summed or not,
/// and odd vs even N
int test_no_duplicates(World& world) {
  test_output t("BoxSurfaceDisplacementRange: no duplicate displacements", world.rank() == 0);

  const auto check_unique = [&](auto ndim_tag, Level n, std::array<std::int64_t, decltype(ndim_tag)::value> N,
                                std::array<bool, decltype(ndim_tag)::value> lattice_summed, bool filter_to_domain,
                                const std::string& what) {
    constexpr std::size_t NDIM = decltype(ndim_tag)::value;
    array_of_bools<NDIM> summed{false};
    for (std::size_t d = 0; d != NDIM; ++d) summed[d] = lattice_summed[d];
    // N.B. the domain filter is only meaningful when the range boundary can land
    // inside the cell at all.
    const auto range = make_range<NDIM>(n, finite<NDIM>(N), summed, {},
                                        filter_to_domain ? std::optional{in_domain_only<NDIM>(summed)} : std::nullopt);
    std::vector<Key<NDIM>> disps;
    for (auto&& disp : range) disps.push_back(disp);
    std::sort(disps.begin(), disps.end());
    const auto last = std::unique(disps.begin(), disps.end());
    t.checkpoint(!disps.empty(), what + ": surface is non-empty");
    t.checkpoint(last == disps.end(), what + ": all displacements distinct");
    if (summed.any()) {  // displacements that differ by a period along a lattice-summed axis are the same displacement
      const auto period = Translation(1) << n;
      std::vector<Key<NDIM>> canonical;
      for (const auto& disp : disps) {
        auto l = disp.translation();
        for (std::size_t d = 0; d != NDIM; ++d)
          if (summed[d]) l[d] = ((l[d] % period) + period) % period;
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
  check_unique(d2, 4, {1, 1}, {false, false}, true, "2D n=4 N={1,1} plain, in-domain");
  check_unique(d2, 4, {1, 1}, {false, false}, false, "2D n=4 N={1,1} plain");
  check_unique(d2, 4, {2, 2}, {false, false}, false, "2D n=4 N={2,2} plain");
  // ... and both parities again with lattice summation on
  check_unique(d2, 4, {1, 1}, {true, true}, false, "2D n=4 N={1,1} lattice-summed");
  check_unique(d2, 4, {2, 2}, {true, true}, false, "2D n=4 N={2,2} lattice-summed");
  check_unique(d2, 4, {2, 3}, {true, true}, false, "2D n=4 N={2,3} lattice-summed (mixed parity)");
  // ... and in 3D, where each face has two free axes rather than one
  check_unique(d3, 3, {1, 1, 1}, {true, true, true}, false, "3D n=3 N={1,1,1} lattice-summed");
  check_unique(d3, 3, {2, 2, 2}, {true, true, true}, false, "3D n=3 N={2,2,2} lattice-summed");
  check_unique(d3, 3, {1, 2, 3}, {true, true, true}, false, "3D n=3 N={1,2,3} lattice-summed (mixed parity)");
  // ... and with lattice summation along some axes only, so that a summed axis (one period, one face)
  // is processed next to an unsummed one (box plus thickness, two faces) in either order
  check_unique(d2, 4, {2, 2}, {true, false}, false, "2D n=4 N={2,2} summed along x only");
  check_unique(d2, 4, {2, 2}, {false, true}, false, "2D n=4 N={2,2} summed along y only");
  check_unique(d3, 3, {2, 2, 2}, {true, false, true}, false, "3D n=3 N={2,2,2} summed along x and z");
  check_unique(d3, 3, {1, 2, 3}, {false, true, true}, false, "3D n=3 N={1,2,3} summed along y and z (mixed parity)");

  return t.end();
}

/// The iteration must not depend on whether a validator is given: a validator that accepts
/// everything must yield the same displacements as no validator at all.
int test_validator_agnostic(World& world) {
  test_output t("BoxSurfaceDisplacementRange: iteration does not depend on the presence of a validator", world.rank() == 0);

  constexpr std::size_t NDIM = 2;
  const auto disps_of = [](std::optional<Validator<NDIM>> validator) {
    const auto range = make_range<NDIM>(4, finite<NDIM>({1, 1}), array_of_bools<NDIM>{false}, {}, std::move(validator));
    std::vector<Key<NDIM>> disps;
    for (auto&& disp : range) disps.push_back(disp);
    std::sort(disps.begin(), disps.end());
    return disps;
  };

  const auto without = disps_of(std::nullopt);
  const auto with = disps_of(keep_all<NDIM>(array_of_bools<NDIM>{false}));
  t.checkpoint(!with.empty(), "2D n=4 N={1,1}: surface is non-empty");
  t.checkpoint(without.size() == with.size(), "2D n=4 N={1,1}: same number of displacements with and without a validator (" +
                                                  std::to_string(without.size()) + " vs " + std::to_string(with.size()) + ")");
  t.checkpoint(without == with, "2D n=4 N={1,1}: same displacements with and without a validator");
  return t.end();
}

/// Every face gets its own probe: odd (lattice-summed) faces sit a half cell out and need no offset,
/// even ones fold onto the source and need an offset along another axis.
int test_odds(World& world) {
  test_output t("BoxSurfaceDisplacementRange: mixed-parity radii", world.rank() == 0);

  const Level level = 4;
  const auto r = [&](std::int64_t N) { return probe_radius(N, level); };
  const auto off = far_offset<3>(level);
  const auto probes = [&](std::array<std::int64_t, 3> N) {
    return probes_of<3>(level, finite<3>(N), array_of_bools<3>{true}, far_reach<3>());
  };

  // the offset of an even face goes along the finite axis with the smallest radius
  check_probes<3>(t, probes({3, 2, 1}), {{{{r(3), 0, 0}}, {{0, r(2), off}}, {{0, 0, r(1)}}}}, "3D n=4 N={3,2,1}");
  check_probes<3>(t, probes({1, 2, 3}), {{{{r(1), 0, 0}}, {{off, r(2), 0}}, {{0, 0, r(3)}}}}, "3D n=4 N={1,2,3}");
  // ... ties are broken by index
  check_probes<3>(t, probes({1, 2, 1}), {{{{r(1), 0, 0}}, {{off, r(2), 0}}, {{0, 0, r(1)}}}}, "3D n=4 N={1,2,1}");
  return t.end();
}

/// Test the degenerate cases in which criterion 2 is out of reach: a single dimension, and level 0.
int test_singular(World& world) {
  test_output t("BoxSurfaceDisplacementRange: singular cases", world.rank() == 0);

  const Level level = 4;
  check_probes<1>(t, probes_of<1>(level, finite<1>({4}), array_of_bools<1>{true}, far_reach<1>()),
                  {{{{probe_radius(4, level)}}}}, "1D n=4 N={4}");
  // at level 0 a half cell is 1 box: N=3 -> 2 boxes, N=2 -> 1 box, N=1 -> 1 box; no offsets, and with unit
  // thickness the innermost layer of the N=2 and N=1 faces is the source box itself
  check_probes<3>(t, probes_of<3>(0, finite<3>({3, 2, 1}), array_of_bools<3>{true}, far_reach<3>()),
                  {{{{1, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}, "3D n=0 N={3,2,1}");
  return t.end();
}

/// Test mixed boundary conditions
int test_mixed(World& world) {
  test_output t("BoxSurfaceDisplacementRange: mixed cases", world.rank() == 0);

  const Level level = 4;
  const auto r = [&](std::int64_t N) { return probe_radius(N, level); };

  Radii<2> radii2d;
  radii2d[1] = 2;
  // no face normal to the unrestricted dimension; the even face's offset goes along it, toward the side with more room
  check_probes<2>(t, probes_of<2>(level, radii2d, array_of_bools<2>{true}, far_reach<2>()),
                  {{std::nullopt, {{-far_offset<2>(level), r(2)}}}}, "2D n=4 N={*,2}");
  radii2d[1] = 3;
  check_probes<2>(t, probes_of<2>(level, radii2d, array_of_bools<2>{true}, far_reach<2>()),
                  {{std::nullopt, {{0, r(3)}}}}, "2D n=4 N={*,3}");
  Radii<3> radii3d;
  radii3d[1] = 2;
  radii3d[2] = 3;
  // a finite axis is preferred over an unrestricted one for the offset (the probe is guaranteed to stay on the face)
  check_probes<3>(t, probes_of<3>(level, radii3d, array_of_bools<3>{true}, far_reach<3>()),
                  {{std::nullopt, {{0, r(2), far_offset<3>(level)}}, {{0, 0, r(3)}}}}, "3D n=4 N={*,2,3}");
  return t.end();
}

/// Test the case in which every box radius is even, so that every face folds onto the source
/// and needs an offset along a second dimension.
int test_all_evens(World& world) {
  test_output t("BoxSurfaceDisplacementRange: pure evens", world.rank() == 0);

  const Level level = 4;
  const auto r = [&](std::int64_t N) { return probe_radius(N, level); };
  const auto off = far_offset<3>(level);
  const auto probes = [&](std::array<std::int64_t, 3> N) {
    return probes_of<3>(level, finite<3>(N), array_of_bools<3>{true}, far_reach<3>());
  };

  check_probes<3>(t, probes({2, 2, 2}), {{{{r(2), off, 0}}, {{off, r(2), 0}}, {{off, 0, r(2)}}}}, "3D n=4 N={2,2,2}");
  check_probes<3>(t, probes({4, 6, 2}), {{{{r(4), 0, off}}, {{0, r(6), off}}, {{off, 0, r(2)}}}}, "3D n=4 N={4,6,2}");
  check_probes<3>(t, probes({6, 2, 2}), {{{{r(6), off, 0}}, {{0, r(2), off}}, {{0, off, r(2)}}}}, "3D n=4 N={6,2,2}");

  return t.end();
}

/// Test differences between lattice-summed and non-lattice summed axes. In particular, offsets
/// are only needed for lattice-summed axes.
int test_lattice_summation_awareness(World& world) {
  test_output t("BoxSurfaceDisplacementRange: offset only where lattice summed", world.rank() == 0);

  const Level level = 4;
  const auto r = [&](std::int64_t N) { return probe_radius(N, level); };
  const auto off3 = far_offset<3>(level);
  const auto probes = [&](std::array<std::int64_t, 3> N, array_of_bools<3> summed) {
    return probes_of<3>(level, finite<3>(N), summed, far_reach<3>());
  };

  // all even, nothing lattice summed => no offsets, each probe is the nearest point of its face
  check_probes<3>(t, probes({2, 2, 2}, array_of_bools<3>{false}), {{{{r(2), 0, 0}}, {{0, r(2), 0}}, {{0, 0, r(2)}}}},
                  "3D N={2,2,2} not summed");
  check_probes<3>(t, probes({4, 2, 6}, array_of_bools<3>{false}), {{{{r(4), 0, 0}}, {{0, r(2), 0}}, {{0, 0, r(6)}}}},
                  "3D N={4,2,6} not summed");
  // ... the same radii with lattice summation on do need the offsets
  check_probes<3>(t, probes({2, 2, 2}, array_of_bools<3>{true}), {{{{r(2), off3, 0}}, {{off3, r(2), 0}}, {{off3, 0, r(2)}}}},
                  "3D N={2,2,2} summed");
  // ... with only the first dimension lattice summed, only its face needs the offset
  check_probes<3>(t, probes({2, 2, 2}, array_of_bools<3>{true, false, false}),
                  {{{{r(2), off3, 0}}, {{0, r(2), 0}}, {{0, 0, r(2)}}}}, "3D N={2,2,2} summed along x only");

  // an unrestricted dimension absorbs the offset, but again only when one is needed
  Radii<2> radii2d;
  radii2d[1] = 2;
  check_probes<2>(t, probes_of<2>(level, radii2d, array_of_bools<2>{false, true}, far_reach<2>()),
                  {{std::nullopt, {{-far_offset<2>(level), r(2)}}}}, "2D N={*,2} summed along y");
  check_probes<2>(t, probes_of<2>(level, radii2d, array_of_bools<2>{false, false}, far_reach<2>()),
                  {{std::nullopt, {{0, r(2)}}}}, "2D N={*,2} not summed");

  return t.end();
}

/// The offset for an even face is derived from the real-space reach of the standard displacements,
/// as supplied by the caller: it is the least number of boxes that the validator does not filter out.
int test_standard_reach(World& world) {
  test_output t("BoxSurfaceDisplacementRange: offset follows the reach of the standard displacements", world.rank() == 0);

  constexpr std::size_t NDIM = 3;
  const Translation bmax = Displacements<NDIM>::bmax_default();  // 4
  // all radii even and lattice summed; check the probe of face `face`
  const auto check = [&](Level n, std::optional<Reach<NDIM>> reach, std::size_t face, std::array<std::int64_t, NDIM> expected,
                         const std::string& what) {
    const auto probe = probes_of<NDIM>(n, finite<NDIM>({2, 2, 2}), array_of_bools<NDIM>{true}, std::move(reach))[face];
    t.checkpoint(probe.has_value(), what + " face " + std::to_string(face) + " exists");
    if (probe) {
      for (size_t d = 0; d < NDIM; d++)
        t.checkpoint((*probe)[d] == expected[d], what + " face " + std::to_string(face) + " component " + std::to_string(d));
    }
  };
  const auto anisotropic_reach = [](double max_distsq, std::array<double, NDIM> cell_width) {
    return Reach<NDIM>{max_distsq, cell_width};
  };

  {
    const Level n = 6;                     // half-cell = 32 boxes, so the cap is far away
    const auto face = probe_radius(2, n);  // 63: the innermost layer of the face at 64

    // nothing is known to be filtered => the surface reaches the source and no offset helps;
    // fall back to the bare face probe, whose norm is the on-site norm and screens nothing
    check(n, std::nullopt, 0, {face, 0, 0}, "nothing filtered");

    // unit cell widths: the real distance of a displacement l along an axis is |l|-1 (see Key::real_distsq_bc),
    // so the least offset beyond real distance sqrt(max_distsq) is floor(sqrt(max_distsq))+2 boxes ...
    check(n, unit_reach<NDIM>(0.), 0, {face, 2, 0}, "unit widths, only the nearest neighbors reached");
    check(n, unit_reach<NDIM>(6.25), 0, {face, 4, 0}, "unit widths, sqrt(max_distsq)=2.5");
    check(n, unit_reach<NDIM>(9.), 0, {face, 5, 0}, "unit widths, sqrt(max_distsq)=3");
    // ... but never more than bmax+1, beyond which nothing is a standard displacement
    check(n, unit_reach<NDIM>(1e6), 0, {face, bmax + 1, 0}, "unit widths, reached beyond bmax");

    // anisotropic cell: the offset goes along the axis with the least *real-space* offset.
    // With sqrt(max_distsq)=5 and width 1 the offset is capped at bmax+1=5 boxes, i.e. 4 real units ...
    check(n, anisotropic_reach(25., {10., 1., 10.}), 0, {face, 5, 0}, "widths {10,1,10}");
    check(n, anisotropic_reach(25., {10., 1., 10.}), 2, {0, 5, face}, "widths {10,1,10}");
    check(n, anisotropic_reach(25., {10., 10., 1.}), 0, {face, 0, 5}, "widths {10,10,1}");
    // ... whereas a wider axis can need fewer boxes yet be farther in real space (width 2: 4 boxes = 6 units) ...
    check(n, anisotropic_reach(25., {1., 2., 100.}), 0, {face, 4, 0}, "widths {1,2,100}");
    // ... or fewer boxes *and* nearer in real space (width 2.6: 3 boxes = 5.2 units vs. width 2: 4 boxes = 6 units)
    check(n, anisotropic_reach(25., {1., 2., 2.6}), 0, {face, 0, 3}, "widths {1,2,2.6}");
  }

  // the offset is capped at a half cell: a step of 2^{n-1}+k folds back down to 2^{n-1}-k
  {
    const Level n = 2;  // half cell = 2 boxes < bmax+1
    check(n, unit_reach<NDIM>(1e6), 0, {probe_radius(2, n), 2, 0}, "n=2, reached beyond bmax");
  }

  // the same cap applies when the offset lands on an unrestricted dimension
  {
    const Level n = 6;
    Radii<2> box_radius;
    box_radius[1] = 2;
    check_probes<2>(t, probes_of<2>(n, box_radius, array_of_bools<2>{true}, unit_reach<2>(4.)),
                    {{std::nullopt, {{-4, probe_radius(2, n)}}}}, "2D N={*,2} sqrt(max_distsq)=2");
  }

  return t.end();
}

/// Skipping a face drops exactly the boxes that lie on no other face; the edge boxes it shares
/// with the remaining faces are still visited through them.
int test_skip_face(World& world) {
  test_output t("BoxSurfaceDisplacementRange: skipping faces", world.rank() == 0);

  // displacements that differ by a period along a lattice-summed axis are equivalent, and which representative
  // is produced depends on the order in which faces are visited; so compare them mapped to [0, 2^n)
  const auto sorted = [](const auto& range, const auto& lattice_summed) {
    using key_type = std::decay_t<decltype(*range.begin())>;
    constexpr std::size_t NDIM = key_type::static_size;
    std::vector<key_type> disps;
    for (auto&& disp : range) {
      auto l = disp.translation();
      const auto period = Translation(1) << disp.level();
      for (std::size_t d = 0; d != NDIM; ++d)
        if (lattice_summed[d]) l[d] = ((l[d] % period) + period) % period;
      disps.emplace_back(disp.level(), l);
    }
    std::sort(disps.begin(), disps.end());
    return disps;
  };

  const auto check = [&](auto ndim_tag, Level n, std::array<std::int64_t, decltype(ndim_tag)::value> N,
                         std::array<bool, decltype(ndim_tag)::value> lattice_summed, std::size_t skipped,
                         const std::string& what) {
    constexpr std::size_t NDIM = decltype(ndim_tag)::value;
    const auto radii = finite<NDIM>(N);
    array_of_bools<NDIM> summed{false};
    for (std::size_t d = 0; d != NDIM; ++d) summed[d] = lattice_summed[d];

    const auto all = sorted(make_range<NDIM>(n, radii, summed), summed);

    auto range = make_range<NDIM>(n, radii, summed);
    range.skip_face(skipped);
    t.checkpoint(range.face_skipped(skipped), what + ": face " + std::to_string(skipped) + " marked skipped");
    const auto rest = sorted(range, summed);

    // is `disp` on the layers of the faces normal to `d`? surface thickness is 1 box on each side of the boundary;
    // along a lattice-summed dimension only the + side is iterated over
    const auto on_face = [&](const Key<NDIM>& disp, std::size_t d) {
      const auto r = radius_in_boxes(N[d], n);
      auto l = disp.translation()[d];
      if (summed[d]) {  // `disp` is canonical, so map the face position into [0, 2^n) as well
        const auto period = Translation(1) << n;
        for (Translation layer = r - 1; layer <= r + 1; ++layer)
          if (((layer % period) + period) % period == l) return true;
        return false;
      }
      return (l >= r - 1 && l <= r + 1) || (l >= -r - 1 && l <= -r + 1);
    };
    const auto on_another_face = [&](const Key<NDIM>& disp) {
      for (std::size_t d = 0; d != NDIM; ++d)
        if (d != skipped && on_face(disp, d)) return true;
      return false;
    };
    std::vector<Key<NDIM>> expected;
    std::copy_if(all.begin(), all.end(), std::back_inserter(expected), on_another_face);

    t.checkpoint(!rest.empty(), what + ": other faces remain");
    t.checkpoint(rest.size() < all.size(), what + ": fewer displacements than the full surface");
    t.checkpoint(rest == expected, what + ": exactly the boxes on the other faces remain");
  };

  constexpr auto d2 = std::integral_constant<std::size_t, 2>{};
  constexpr auto d3 = std::integral_constant<std::size_t, 3>{};
  // skipping the first face (whose edges would otherwise be excluded from the later faces) and a later one
  check(d2, 4, {2, 2}, {false, false}, 0, "2D n=4 N={2,2} plain, skip x");
  check(d2, 4, {2, 2}, {false, false}, 1, "2D n=4 N={2,2} plain, skip y");
  check(d2, 4, {1, 2}, {true, true}, 0, "2D n=4 N={1,2} lattice-summed, skip x");
  check(d3, 3, {1, 2, 1}, {true, true, true}, 0, "3D n=3 N={1,2,1} lattice-summed, skip x");
  check(d3, 3, {1, 2, 1}, {true, true, true}, 1, "3D n=3 N={1,2,1} lattice-summed, skip y");
  check(d3, 3, {1, 2, 1}, {true, true, true}, 2, "3D n=3 N={1,2,1} lattice-summed, skip z");
  // ... and with lattice summation along some axes only
  check(d2, 4, {2, 2}, {true, false}, 0, "2D n=4 N={2,2} summed along x only, skip x");
  check(d2, 4, {2, 2}, {true, false}, 1, "2D n=4 N={2,2} summed along x only, skip y");
  check(d3, 3, {2, 1, 2}, {true, false, true}, 1, "3D n=3 N={2,1,2} summed along x and z, skip y");
  check(d3, 3, {2, 1, 2}, {true, false, true}, 2, "3D n=3 N={2,1,2} summed along x and z, skip z");

  // skipping every face leaves nothing
  {
    auto range = make_range<3>(3, finite<3>({1, 2, 1}), array_of_bools<3>{true});
    for (std::size_t d = 0; d != 3; ++d) range.skip_face(d);
    t.checkpoint(range.begin() == range.end(), "3D n=3 N={1,2,1}: skipping all faces leaves nothing");
  }
  // a box that is not hollow along some dimension has every box on that face; skipping it changes nothing
  // since every box is on the other face as well
  {
    const auto radii = finite<2>({1, 1});  // at n=1 the radius is 1 box = the thickness
    const auto all = sorted(make_range<2>(1, radii, array_of_bools<2>{false}), array_of_bools<2>{false});
    auto range = make_range<2>(1, radii, array_of_bools<2>{false});
    range.skip_face(0);
    t.checkpoint(!all.empty() && sorted(range, array_of_bools<2>{false}) == all, "2D n=1 N={1,1}: skipping the face of a non-hollow dimension changes nothing");
  }
  // an unrestricted dimension has no face and the others are unaffected
  {
    Radii<2> radii;
    radii[1] = 2;
    const auto all = sorted(make_range<2>(4, radii, array_of_bools<2>{false}), array_of_bools<2>{false});
    auto range = make_range<2>(4, radii, array_of_bools<2>{false});
    range.skip_face(1);
    t.checkpoint(range.begin() == range.end(), "2D n=4 N={*,2}: skipping the only face leaves nothing");
    t.checkpoint(!all.empty(), "2D n=4 N={*,2}: the only face is non-empty");
  }

  return t.end();
}

/// Faces (or whole surfaces) lying entirely outside of a finite domain: the iterator must skip them
/// and produce exactly the in-domain part of the surface.
int test_faces_outside_domain(World& world) {
  test_output t("BoxSurfaceDisplacementRange: faces outside of the domain", world.rank() == 0);

  constexpr std::size_t NDIM = 2;
  const Level n = 4;
  const auto twon = Translation(1) << n;
  const array_of_bools<NDIM> not_summed{false};
  const auto in_domain = [&](const Key<NDIM>& disp, const Key<NDIM>& center) {
    for (std::size_t d = 0; d != NDIM; ++d) {
      const auto x = center.translation()[d] + disp.translation()[d];
      if (x < 0 || x >= twon) return false;
    }
    return true;
  };
  const auto check = [&](std::array<std::int64_t, NDIM> N, Key<NDIM> center, const std::string& what) {
    // without a validator: the whole surface, in and out of the domain
    std::vector<Key<NDIM>> expected;
    for (auto&& disp : make_range<NDIM>(n, finite<NDIM>(N), not_summed, {}, {}, center))
      if (in_domain(disp, center)) expected.push_back(disp);
    std::sort(expected.begin(), expected.end());
    // with the domain filter
    std::vector<Key<NDIM>> actual;
    for (auto&& disp : make_range<NDIM>(n, finite<NDIM>(N), not_summed, {}, in_domain_only<NDIM>(not_summed), center))
      actual.push_back(disp);
    std::sort(actual.begin(), actual.end());
    t.checkpoint(actual == expected, what + ": exactly the in-domain part of the surface (" + std::to_string(actual.size()) +
                                         " vs " + std::to_string(expected.size()) + ")");
    return expected.size();
  };

  // both faces normal to x are out (at -8 and 24), those normal to y are partly in
  t.checkpoint(check({2, 1}, Key<NDIM>(n, {8, 8}), "N={2,1} centered") > 0, "N={2,1} centered: something is in the domain");
  // the + face normal to x is out (at 20), the - face is in
  t.checkpoint(check({1, 1}, Key<NDIM>(n, {12, 8}), "N={1,1} off-center") > 0, "N={1,1} off-center: something is in the domain");
  // everything is out
  t.checkpoint(check({2, 2}, Key<NDIM>(n, {8, 8}), "N={2,2} centered") == 0, "N={2,2} centered: nothing is in the domain");

  return t.end();
}

/// The validator's notion of "standard displacement" must match Displacements: with lattice summation along
/// any axis the standard displacements are clipped to 2^n-1 boxes along *every* axis, summed or not.
int test_validator_bmax(World& world) {
  test_output t("BoxSurfaceDisplacementValidator: bmax of the standard displacements", world.rank() == 0);

  constexpr std::size_t NDIM = 2;
  const Level n = 2;  // 4 boxes per axis, fewer than bmax+1 in 2D
  const Translation bmax = Displacements<NDIM>::bmax_default();
  t.checkpoint(bmax > 3, "2D bmax exceeds 2^n-1 at n=2");
  // lattice summed along x only, infinite domain, and every standard displacement was reached
  const Validator<NDIM> v(array_of_bools<NDIM>{true}, array_of_bools<NDIM>{true, false}, unit_reach<NDIM>(1e6));
  const auto keeps = [&](std::int64_t lx, std::int64_t ly) {
    typename BoxSurfaceDisplacementRange<NDIM>::PointPattern dest;
    dest[0] = 2 + lx;
    dest[1] = 2 + ly;
    std::optional<Key<NDIM>> disp(Key<NDIM>(n, {lx, ly}));
    return v(n, dest, disp);
  };
  t.checkpoint(!keeps(0, 3), "(0,3): within 2^n-1 along the unsummed axis => standard => filtered");
  t.checkpoint(keeps(0, 4), "(0,4): beyond 2^n-1 along the unsummed axis => not standard => kept");
  t.checkpoint(!keeps(4, 0), "(4,0): folds onto 0 along the summed axis => standard => filtered");
  t.checkpoint(!keeps(3, 0), "(3,0): folds onto -1 along the summed axis => standard => filtered");

  // level 0 must not trip the bit tricks: everything folds onto the single box
  {
    typename BoxSurfaceDisplacementRange<NDIM>::PointPattern dest;
    dest[0] = 0;
    dest[1] = 0;
    std::optional<Key<NDIM>> disp(Key<NDIM>(0, {0, 0}));
    t.checkpoint(!v(0, dest, disp), "level 0: the on-site displacement is standard => filtered");
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
  errors += test_standard_reach(world);
  errors += test_skip_face(world);
  errors += test_faces_outside_domain(world);
  errors += test_validator_bmax(world);

  world.gop.fence();
  madness::finalize();
  return errors;
}
