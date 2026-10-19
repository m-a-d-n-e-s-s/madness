#include <madness/mra/mra.h>
#include <madness/mra/displacements.h>
#include <madness/world/test_utilities.h>

#include <limits>

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

/// ... and, if the face's dimension is lattice summed, folded into the cell to the representative nearest to the source:
/// -1 for even N (the face folds onto the source), 2^{n-1}-1 for odd N (a half cell away)
Translation probe_radius_summed(std::int64_t N, Level n) {
  const auto period = Translation(1) << n;
  auto l = ((probe_radius(N, n) % period) + period) % period;
  if (l > period / 2) l -= period;
  return l;
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
  const auto r = [&](std::int64_t N) { return probe_radius_summed(N, level); };
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
                  {{{{probe_radius_summed(4, level)}}}}, "1D n=4 N={4}");
  // at level 0 the cell is a single box, onto which every lattice-summed face folds; no offsets
  check_probes<3>(t, probes_of<3>(0, finite<3>({3, 2, 1}), array_of_bools<3>{true}, far_reach<3>()),
                  {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}}, "3D n=0 N={3,2,1}");
  return t.end();
}

/// Test mixed boundary conditions
int test_mixed(World& world) {
  test_output t("BoxSurfaceDisplacementRange: mixed cases", world.rank() == 0);

  const Level level = 4;
  const auto r = [&](std::int64_t N) { return probe_radius_summed(N, level); };

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
  const auto r = [&](std::int64_t N) { return probe_radius_summed(N, level); };
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
  const auto rs = [&](std::int64_t N) { return probe_radius_summed(N, level); };
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
  check_probes<3>(t, probes({2, 2, 2}, array_of_bools<3>{true}), {{{{rs(2), off3, 0}}, {{off3, rs(2), 0}}, {{off3, 0, rs(2)}}}},
                  "3D N={2,2,2} summed");
  // ... with only the first dimension lattice summed, only its face needs the offset
  check_probes<3>(t, probes({2, 2, 2}, array_of_bools<3>{true, false, false}),
                  {{{{rs(2), off3, 0}}, {{0, r(2), 0}}, {{0, 0, r(2)}}}}, "3D N={2,2,2} summed along x only");

  // an unrestricted dimension absorbs the offset, but again only when one is needed
  Radii<2> radii2d;
  radii2d[1] = 2;
  check_probes<2>(t, probes_of<2>(level, radii2d, array_of_bools<2>{false, true}, far_reach<2>()),
                  {{std::nullopt, {{-far_offset<2>(level), rs(2)}}}}, "2D N={*,2} summed along y");
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
    const auto face = probe_radius_summed(2, n);  // -1: the innermost layer of the face at 64, folded onto the source

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
    check(n, unit_reach<NDIM>(1e6), 0, {probe_radius_summed(2, n), 2, 0}, "n=2, reached beyond bmax");
  }

  // the same cap applies when the offset lands on an unrestricted dimension
  {
    const Level n = 6;
    Radii<2> box_radius;
    box_radius[1] = 2;
    check_probes<2>(t, probes_of<2>(n, box_radius, array_of_bools<2>{true}, unit_reach<2>(4.)),
                    {{std::nullopt, {{-4, probe_radius_summed(2, n)}}}}, "2D N={*,2} sqrt(max_distsq)=2");
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

/// The standard displacements must be ordered by the real-space distance measured with the FunctionDefaults cell,
/// since FunctionImpl::do_apply screens them shell by shell in that metric; the order must follow the cell when it changes.
int test_standard_displacements_order(World& world) {
  test_output t("Displacements: ordered by real-space distance in the FunctionDefaults cell", world.rank() == 0);

  constexpr std::size_t NDIM = 3;
  const Level n = 4;
  const auto is_sorted_by = [&](const std::vector<Key<NDIM>>& disps, auto distsq) {
    for (std::size_t i = 1; i < disps.size(); ++i)
      if (distsq(disps[i]) < distsq(disps[i - 1])) return false;
    return true;
  };
  const auto check_order = [&](const std::string& what) {
    const auto& width = FunctionDefaults<NDIM>::get_cell_width();
    Displacements<NDIM> displacements;  // builds the lists if this is the first use
    const auto& plain = displacements.get_disp(n, array_of_bools<NDIM>{false});
    const auto& summed = displacements.get_disp(n, array_of_bools<NDIM>{true});
    t.checkpoint(!plain.empty() && !summed.empty(), what + ": lists are populated");
    t.checkpoint(is_sorted_by(plain, [&](const Key<NDIM>& k) { return k.real_distsq(width); }),
                 what + ": plain displacements ordered by real distance");
    t.checkpoint(is_sorted_by(summed, [&](const Key<NDIM>& k) { return k.real_distsq_bc(array_of_bools<NDIM>{true}, width); }),
                 what + ": lattice-summed displacements ordered by real distance");
  };

  const Tensor<double> cell0 = copy(FunctionDefaults<NDIM>::get_cell());
  check_order("default (cubic) cell");

  // an anisotropic cell: the order in boxes and in real space differ, e.g. (0,2,0) is 1 box but 10 units out
  // while (3,0,0) is 2 boxes but 2 units out
  Tensor<double> cell(NDIM, 2);
  cell(0, 0) = 0.; cell(0, 1) = 1.;
  cell(1, 0) = 0.; cell(1, 1) = 10.;
  cell(2, 0) = 0.; cell(2, 1) = 10.;
  FunctionDefaults<NDIM>::set_cell(cell);
  check_order("anisotropic cell (1,10,10) set after the lists were built");

  FunctionDefaults<NDIM>::set_cell(cell0);
  check_order("cubic cell restored");

  return t.end();
}

/// The lists for a set of lattice-summed axes must be ordered by the metric of those axes, whether or not the
/// default boundary conditions are periodic along more of them: a list sorted for a superset would place a
/// displacement wrapped along an axis the kernel does not sum next to the source, where for that kernel it is
/// a cell away. The rest-of-crystal lists (get_disp_images) must be the same displacements ordered by the
/// distance to the nearest image other than the home cell.
int test_images_and_subset_displacements_order(World& world) {
  test_output t("Displacements: lists of every axis set ordered by their own metric, rest-of-crystal lists by image distance", world.rank() == 0);

  constexpr std::size_t NDIM = 3;
  const Level n = 4;
  const Translation twon = Translation(1) << n;
  const int bmax = Displacements<NDIM>::bmax_default();
  const auto is_sorted_by = [&](const std::vector<Key<NDIM>>& disps, auto distsq) {
    for (std::size_t i = 1; i < disps.size(); ++i)
      if (distsq(disps[i]) < distsq(disps[i - 1])) return false;
    return true;
  };
  const auto wraps_along = [&](const Key<NDIM>& k, std::size_t d) { return std::abs(k.translation()[d]) > bmax; };

  const auto check = [&](const array_of_bools<NDIM>& axes, const std::string& what) {
    const auto& width = FunctionDefaults<NDIM>::get_cell_width();
    Displacements<NDIM> displacements;
    const auto& summed = displacements.get_disp(n, axes);
    const auto& images = displacements.get_disp_images(n, axes);
    t.checkpoint(!summed.empty() && summed.size() == images.size() &&
                     std::is_permutation(summed.begin(), summed.end(), images.begin()),
                 what + ": the rest-of-crystal list is a reordering of the lattice-summed list");
    t.checkpoint(is_sorted_by(summed, [&](const Key<NDIM>& k) { return k.real_distsq_bc(axes, width); }),
                 what + ": lattice-summed list ordered by the distance modulo the lattice of these axes");
    t.checkpoint(is_sorted_by(images, [&](const Key<NDIM>& k) { return k.real_distsq_images(axes, width); }),
                 what + ": rest-of-crystal list ordered by the distance to the nearest non-home image of these axes");
    // a displacement wraps (|l| > bmax) only along a lattice-summed axis
    bool wrap_ok = true, any_wrap = false;
    for (const auto& k : summed)
      for (std::size_t d = 0; d < NDIM; ++d)
        if (wraps_along(k, d)) { any_wrap = true; if (!axes[d]) wrap_ok = false; }
    t.checkpoint(wrap_ok && any_wrap, what + ": displacements wrap along the lattice-summed axes only");
    // the first rest-of-crystal displacement is adjacent to an image: wrapped along a summed axis, a cell from the source
    const Key<NDIM>& first = images.front();
    t.checkpoint(first.distsq_images(axes) == 1 && first.real_distsq_images(axes, width) == 0.0 &&
                     std::count_if(first.translation().begin(), first.translation().end(),
                                   [&](Translation l) { return std::abs(l) == twon - 1; }) == 1,
                 what + ": the first rest-of-crystal displacement is one box from the nearest image");
    // ... along the narrowest summed axis: every image-adjacent displacement is at least distance 0 from its
    // image, so the tie is broken by the distance between the box centers, which is the axis width
    std::size_t narrowest = NDIM;
    for (std::size_t d = 0; d < NDIM; ++d)
      if (axes[d] && (narrowest == NDIM || width(static_cast<long>(d)) < width(static_cast<long>(narrowest)))) narrowest = d;
    t.checkpoint(wraps_along(first, narrowest), what + ": ... along the narrowest summed axis");
    // the image-adjacent displacements (real distance 0, one box) lead the list, ordered by center distance
    std::size_t nadjacent = 0;
    bool centers_ordered = true;
    for (std::size_t i = 0; i < images.size(); ++i) {
      if (images[i].real_distsq_images(axes, width) != 0.0 || images[i].distsq_images(axes) != 1) break;
      ++nadjacent;
      if (i > 0 && images[i].real_distsq_images_centers(axes, width) < images[i - 1].real_distsq_images_centers(axes, width))
        centers_ordered = false;
    }
    t.checkpoint(nadjacent == 2 * static_cast<std::size_t>(std::count(axes.begin(), axes.end(), true)) && centers_ordered,
                 what + ": the image-adjacent displacements lead the list, ordered by the distance between box centers");
    // ... and the home displacement is a cell away from its nearest image, along the narrowest summed axis,
    // so it follows every displacement adjacent to an image
    const Key<NDIM> home(n, Vector<Translation, NDIM>(0));
    double home_distsq = std::numeric_limits<double>::max();
    for (std::size_t d = 0; d < NDIM; ++d)
      if (axes[d]) home_distsq = std::min(home_distsq, std::pow(width(static_cast<long>(d)) * static_cast<double>(twon - 1), 2));
    const auto home_it = std::find(images.begin(), images.end(), home);
    const auto last_adjacent = std::find_if(images.rbegin(), images.rend(),
                                            [&](const Key<NDIM>& k) { return k.real_distsq_images(axes, width) == 0.0; });
    t.checkpoint(home.real_distsq_images(axes, width) == home_distsq && home_it != images.end() &&
                     last_adjacent != images.rend() && (last_adjacent.base() - 1) < home_it,
                 what + ": the home displacement is a cell from its image and follows every image-adjacent displacement");
  };

  check(array_of_bools<NDIM>{true}, "summed along every axis");
  check(array_of_bools<NDIM>(false, false, true), "summed along z only");
  check(array_of_bools<NDIM>(true, true, false), "summed along x and y");

  // the cell changes after the lists exist: every list, in both orders, follows
  const Tensor<double> cell0 = copy(FunctionDefaults<NDIM>::get_cell());
  Tensor<double> cell(NDIM, 2);
  cell(0, 0) = 0.; cell(0, 1) = 1.;
  cell(1, 0) = 0.; cell(1, 1) = 10.;
  cell(2, 0) = 0.; cell(2, 1) = 2.5;
  FunctionDefaults<NDIM>::set_cell(cell);
  check(array_of_bools<NDIM>{true}, "anisotropic cell, summed along every axis");
  check(array_of_bools<NDIM>(false, false, true), "anisotropic cell, summed along z only");
  // the narrowest summed axis is not x, so the lexicographic tie-break would put the wrong displacement first
  cell(0, 1) = 10.; cell(1, 1) = 1.;
  FunctionDefaults<NDIM>::set_cell(cell);
  check(array_of_bools<NDIM>{true}, "anisotropic cell, narrowest along y, summed along every axis");
  check(array_of_bools<NDIM>(true, true, false), "anisotropic cell, narrowest along y, summed along x and y");
  FunctionDefaults<NDIM>::set_cell(cell0);

  return t.end();
}

/// Key::distsq_images / real_distsq_images take the nearest image on every axis unless that is the home cell,
/// in which case one axis, the cheapest, moves to its next image. Check that shortcut against the minimum over
/// the lattice vectors R != 0 taken directly, for every axis set and an anisotropic cell.
int test_images_distance(World& world) {
  test_output t("Key::distsq_images: minimum over the lattice images other than the home cell", world.rank() == 0);

  constexpr std::size_t NDIM = 3;
  Tensor<double> width(3L);
  width(0L) = 1.0; width(1L) = 10.0; width(2L) = 2.5;

  // min over R != 0 (R_d = 0 along nonperiodic axes) of sum_d axis_distsq(d, l_d + R_d 2^n), R_d in [-3, 3]
  const auto brute_force = [&](const Key<NDIM>& k, const array_of_bools<NDIM>& per, auto axis_distsq) {
    const Translation twon = Translation(1) << k.level();
    std::optional<double> best;
    for (int rx = -3; rx <= 3; ++rx)
      for (int ry = -3; ry <= 3; ++ry)
        for (int rz = -3; rz <= 3; ++rz) {
          const int R[3] = {rx, ry, rz};
          bool allowed = true, nonzero = false;
          for (std::size_t d = 0; d < NDIM; ++d) {
            if (R[d] != 0 && !per[d]) allowed = false;
            if (R[d] != 0) nonzero = true;
          }
          if (!allowed || !nonzero) continue;
          double s = 0;
          for (std::size_t d = 0; d < NDIM; ++d) s += axis_distsq(d, k.translation()[d] + R[d] * twon);
          if (!best || s < *best) best = s;
        }
    return *best;
  };
  const auto box_axis = [](std::size_t, Translation l) { return double(l * l); };
  const auto real_axis = [&](std::size_t d, Translation l) {
    const double a = width(static_cast<long>(d)) * static_cast<double>(std::max<Translation>(std::abs(l) - 1, 0));
    return a * a;
  };

  std::size_t nchecked = 0, nbad = 0;
  // |l| up to 2^n + 2 exercises translations beyond the ones Displacements produces (|l| < 2^n)
  for (Level n = 0; n <= 4; ++n) {
    const Translation lmax = (Translation(1) << n) + 2;
    for (int mask = 1; mask < 8; ++mask) {
      array_of_bools<NDIM> per{false};
      for (std::size_t d = 0; d < NDIM; ++d) per[d] = (mask >> d) & 1;
      for (Translation x = -lmax; x <= lmax; ++x)
        for (Translation y = -lmax; y <= lmax; ++y)
          for (Translation z = -lmax; z <= lmax; ++z) {
            const Key<NDIM> k(n, Vector<Translation, NDIM>{x, y, z});
            ++nchecked;
            const double box = double(k.distsq_images(per)), box_bf = brute_force(k, per, box_axis);
            const double real = k.real_distsq_images(per, width), real_bf = brute_force(k, per, real_axis);
            if (box != box_bf || std::abs(real - real_bf) > 1e-12 * std::max(1.0, real_bf)) {
              if (++nbad <= 5 && world.rank() == 0)
                print("  mismatch: level", n, "axes", per, "l =", k.translation(), " boxes", box, "vs", box_bf, " real", real, "vs", real_bf);
            }
          }
    }
  }
  t.checkpoint(nchecked > 0 && nbad == 0, "matches the direct minimum for levels 0..4, every axis set, |l| <= 2^n + 2");

  // a few by hand, level 4 summed along z in a cell of width 2.5 along z: the image of the source is 16 boxes away
  {
    const array_of_bools<NDIM> per(false, false, true);
    const auto key = [](Translation x, Translation y, Translation z) { return Key<NDIM>(4, Vector<Translation, NDIM>{x, y, z}); };
    t.checkpoint(key(0, 0, 0).distsq_images(per) == 256 && key(0, 0, 0).real_distsq_images(per, width) == 2.5 * 15 * 2.5 * 15,
                 "the home displacement is a cell away from the nearest image");
    t.checkpoint(key(0, 0, 15).distsq_images(per) == 1 && key(0, 0, 15).real_distsq_images(per, width) == 0.0,
                 "l = 2^n - 1 is one box from the image");
    t.checkpoint(key(0, 0, 8).distsq_images(per) == 64 && key(0, 0, -8).distsq_images(per) == 64,
                 "l = 2^(n-1) is as far from the image as from the source");
    t.checkpoint(key(3, 0, 0).distsq_images(per) == 9 + 256,
                 "an unsummed axis adds its own distance");
  }
  // deep levels, where the offset to the next-nearest image (~2^n boxes) squared does not fit in 64 bits:
  // Displacements sorts every level up to 61 with this metric. The box metric saturates per axis, the
  // real metric (double) does not; levels <= 60 keep the brute force's l + 3 * 2^n within Translation
  {
    const uint64_t cap = Key<NDIM>::distsq_images_axis_max();
    const auto sat_axis = [cap](std::size_t, Translation l) {
      const double a = std::abs(double(l));
      return std::min(a * a, double(cap));
    };
    const array_of_bools<NDIM> per(false, false, true);
    std::size_t nbad_deep = 0;
    bool monotone = true;
    for (Level n : {31, 32, 40, 60}) {
      const Translation twon = Translation(1) << n;
      const auto key = [n](Translation x, Translation z) { return Key<NDIM>(n, Vector<Translation, NDIM>{x, 0, z}); };
      // l_z from the source (0) to one box from its image (2^n - 1): the box metric never increases
      uint64_t prev = std::numeric_limits<uint64_t>::max();
      for (Translation z : {Translation(0), Translation(1), twon >> 2, twon >> 1, twon - (Translation(1) << 20), twon - 2, twon - 1}) {
        for (Translation x : {Translation(0), Translation(3)}) {
          const Key<NDIM> k = key(x, z);
          const double box = double(k.distsq_images(per)), box_bf = brute_force(k, per, sat_axis);
          const double real = k.real_distsq_images(per, width), real_bf = brute_force(k, per, real_axis);
          if (box != box_bf || std::abs(real - real_bf) > 1e-12 * std::max(1.0, real_bf)) {
            if (++nbad_deep <= 5 && world.rank() == 0)
              print("  mismatch: level", n, "l =", k.translation(), " boxes", box, "vs", box_bf, " real", real, "vs", real_bf);
          }
        }
        const uint64_t d = key(0, z).distsq_images(per);
        monotone = monotone && d <= prev;
        prev = d;
      }
      monotone = monotone && prev == 1;
    }
    t.checkpoint(nbad_deep == 0, "levels 31..60: box metric saturates per axis, real metric matches the direct minimum");
    t.checkpoint(monotone, "levels 31..60: box metric non-increasing from the source to one box from the image");
    t.checkpoint(Key<NDIM>(40, Vector<Translation, NDIM>(0)).distsq_images(per) == cap, "level 40: the home displacement saturates");
  }
  // without a periodic axis there is no image
  {
    bool threw = false;
    try { Key<NDIM>(4, Vector<Translation, NDIM>(0)).distsq_images(array_of_bools<NDIM>{false}); }
    catch (const MadnessException&) { threw = true; }
    t.checkpoint(threw, "no periodic axis => throws");
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
  errors += test_standard_displacements_order(world);
  errors += test_images_and_subset_displacements_order(world);
  errors += test_images_distance(world);

  world.gop.fence();
  madness::finalize();
  return errors;
}
