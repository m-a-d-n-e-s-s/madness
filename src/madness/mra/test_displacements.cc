/// \file test_displacements.cc
/// \brief Tests for BoxSurfaceDisplacementRange — the surface-of-the-box
///        displacement enumerator used for range-restricted operators.
///
/// This file is very concerned with the case where all periodic dimensions are even.
/// Why is that special? `box_radius[d] == N` is the kernel range in
/// half-simulation-cells, so at level `n` the range boundary sits `r = N*2^{n-1}`
/// boxes from the center while the lattice period is `2^n` boxes. For odd `N`,
/// `r` is a half-period offset and the boundary hyperface maps onto a cell face.
/// For even `N`, `r` is an exact multiple of the period, so under lattice
/// summation the hyperface is equivalent to a hyperplane through the origin —
/// and the natural probing displacement (the hyperface origin, `r` along one
/// axis) canonicalizes to the *zero* displacement. This is expected to be large
/// due to the interaction between the center and itself, so it isn't usable.
///
/// The tests below cover both the iteration bounds and each branch of
/// `compute_probe()`.

#include <madness/mra/mra.h>
#include <madness/mra/displacements.h>
#include <madness/world/test_utilities.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

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

/// least nonnegative residue of \p x modulo \p m
Translation mod_positive(Translation x, Translation m) {
  const auto r = x % m;
  return (r < 0) ? r + m : r;
}

/// a validator that rejects everything it is asked about
/// forces compute_probe()'s shrink loop to stop on its first trial, exposing the
/// displacement vector exactly as the fallback block initialized it
template <std::size_t NDIM>
typename BoxSurfaceDisplacementRange<NDIM>::Validator reject_all() {
  return [](Level, const typename BoxSurfaceDisplacementRange<NDIM>::PointPattern&,
            std::optional<Key<NDIM>>&) -> bool { return false; };
}

/// a validator that accepts everything it is asked about
/// drives compute_probe()'s greedy stage all the way down to |disp| == 1, so the
/// shortest still-valid trial is the one that has to come back
template <std::size_t NDIM>
typename BoxSurfaceDisplacementRange<NDIM>::Validator accept_all() {
  return [](Level, const typename BoxSurfaceDisplacementRange<NDIM>::PointPattern&,
            std::optional<Key<NDIM>>&) -> bool { return true; };
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

}  // namespace

/// The surface must span the full extent of *every* dimension.
///
/// Regression guard: the per-dimension bounds are computed in one loop in the
/// constructor. If that loop ever returns early again (it used to, once the
/// probe computation was folded into it), the trailing dimensions keep
/// default-constructed bounds and the surface silently collapses along them.
/// Deliberately uses a different radius per axis so a collapsed axis cannot
/// hide behind a neighbour's bounds.
int test_surface_extent(World& world) {
  test_output t("BoxSurfaceDisplacementRange: surface spans every dimension", world.rank() == 0);

  const auto check_extent = [&](auto ndim_tag, Level n, std::array<std::int64_t, decltype(ndim_tag)::value> N,
                                std::int64_t thickness, const std::string& what) {
    constexpr std::size_t NDIM = decltype(ndim_tag)::value;
    Radii<NDIM> box_radius, surface_thickness;
    for (std::size_t d = 0; d != NDIM; ++d) {
      box_radius[d] = N[d];
      surface_thickness[d] = thickness;
    }
    BoxSurfaceDisplacementRange<NDIM> range(centered_key<NDIM>(n), box_radius, surface_thickness,
                                            array_of_bools<NDIM>{false});

    Vector<Translation, NDIM> lo(std::numeric_limits<Translation>::max());
    Vector<Translation, NDIM> hi(std::numeric_limits<Translation>::min());
    std::size_t count = 0;
    for (auto&& disp : range) {
      for (std::size_t d = 0; d != NDIM; ++d) {
        lo[d] = std::min(lo[d], disp.translation()[d]);
        hi[d] = std::max(hi[d], disp.translation()[d]);
      }
      ++count;
    }
    t.checkpoint(count > 0, what + ": surface is non-empty");
    for (std::size_t d = 0; d != NDIM; ++d) {
      const auto expected = radius_in_boxes(N[d], n) + thickness;
      t.checkpoint(lo[d] == -expected, what + ": dim " + std::to_string(d) + " reaches -" + std::to_string(expected));
      t.checkpoint(hi[d] == expected, what + ": dim " + std::to_string(d) + " reaches +" + std::to_string(expected));
    }
  };

  check_extent(std::integral_constant<std::size_t, 2>{}, 4, {1, 1}, 1, "2D n=4 N={1,1}");
  check_extent(std::integral_constant<std::size_t, 2>{}, 4, {1, 2}, 1, "2D n=4 N={1,2}");
  check_extent(std::integral_constant<std::size_t, 2>{}, 3, {2, 2}, 0, "2D n=3 N={2,2} hard");
  check_extent(std::integral_constant<std::size_t, 3>{}, 3, {1, 2, 3}, 1, "3D n=3 N={1,2,3}");

  return t.end();
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

/// Whenever the hyperface origin is a usable probe, it is the one chosen.
///
/// That covers every case except "lattice-summed with even N on all finite
/// axes": no lattice summation at all, or an odd N, or NDIM==1, or n==0.
int test_probe_hyperface_origin(World& world) {
  test_output t("BoxSurfaceDisplacementRange: probe on the hyperface origin", world.rank() == 0);

  const auto check_probe = [&](Level n, std::array<std::int64_t, 3> N, bool lattice_summed,
                               const std::string& what) {
    constexpr std::size_t NDIM = 3;
    Radii<NDIM> box_radius, surface_thickness;
    for (std::size_t d = 0; d != NDIM; ++d) {
      box_radius[d] = N[d];
      surface_thickness[d] = 1;
    }
    BoxSurfaceDisplacementRange<NDIM> range(centered_key<NDIM>(n), box_radius, surface_thickness,
                                            array_of_bools<NDIM>{lattice_summed}, in_domain_only<NDIM>());
    const auto probe = range.probing_displacement().translation();
    t.checkpoint(probe[0] == radius_in_boxes(N[0], n), what + ": probe[0] == range boundary");
    t.checkpoint(probe[1] == 0 && probe[2] == 0, what + ": probe is axial");
  };

  check_probe(4, {1, 1, 1}, false, "3D n=4 N={1,1,1} plain");
  check_probe(4, {2, 2, 2}, false, "3D n=4 N={2,2,2} plain");
  check_probe(4, {1, 1, 1}, true, "3D n=4 N={1,1,1} lattice-summed (odd)");
  check_probe(4, {3, 2, 2}, true, "3D n=4 N={3,2,2} lattice-summed (leading odd)");
  check_probe(0, {2, 2, 2}, true, "3D n=0 N={2,2,2} lattice-summed");

  return t.end();
}

/// The point of the whole exercise: for an even, lattice-summed range the probe
/// must not be lattice-equivalent to the zero displacement.
///
/// The hyperface origin is `r = N*2^{n-1}` along one axis, an exact multiple of
/// the period `2^n` when N is even, so it canonicalizes to zero. Screening the
/// entire range boundary on the norm of a zero displacement is meaningless —
/// `op->norm(level, probe, source)` would report the on-site norm and the
/// surface would never be screened correctly.
int test_probe_even_lattice_sum(World& world) {
  test_output t("BoxSurfaceDisplacementRange: even lattice-summed probe avoids the origin", world.rank() == 0);

  constexpr std::size_t NDIM = 3;
  const auto lattice_summed = array_of_bools<NDIM>{true};

  const auto check_probe = [&](Level n, std::array<std::int64_t, NDIM> N, const std::string& what) {
    Radii<NDIM> box_radius, surface_thickness;
    std::array<KernelRange, NDIM> kernel_range;
    for (std::size_t d = 0; d != NDIM; ++d) {
      box_radius[d] = N[d];
      surface_thickness[d] = 1;
      kernel_range[d] = KernelRange(static_cast<unsigned int>(N[d]));
    }
    const auto distsq = [&](const Key<NDIM>& disp) -> double {
      return disp.real_distsq_bc(lattice_summed, FunctionDefaults<NDIM>::get_cell_width());
    };
    // a standard-displacement loop that already reached everything nearby, so the
    // validator rejects any displacement it recognizes as already-computed
    BoxSurfaceDisplacementValidator<NDIM> validator(/* is_infinite_domain= */ array_of_bools<NDIM>{true},
                                                    lattice_summed, kernel_range, distsq, 1.e9);
    BoxSurfaceDisplacementRange<NDIM> range(centered_key<NDIM>(n), box_radius, surface_thickness,
                                            lattice_summed, validator);

    const auto probe = range.probing_displacement().translation();
    const auto period = Translation(1) << n;

    bool equivalent_to_origin = true;
    for (std::size_t d = 0; d != NDIM; ++d)
      if (mod_positive(probe[d], period) != 0) equivalent_to_origin = false;
    t.checkpoint(!equivalent_to_origin, what + ": probe is not lattice-equivalent to the origin");

    for (std::size_t d = 0; d != NDIM; ++d) {
      const auto bound = radius_in_boxes(N[d], n) + 1;  // surface thickness 1
      t.checkpoint(std::abs(probe[d]) <= bound, what + ": probe[" + std::to_string(d) + "] inside the box");
    }
  };

  check_probe(4, {2, 2, 2}, "3D n=4 N={2,2,2}");
  check_probe(5, {2, 2, 2}, "3D n=5 N={2,2,2}");
  check_probe(4, {4, 4, 4}, "3D n=4 N={4,4,4}");
  check_probe(4, {2, 4, 6}, "3D n=4 N={2,4,6}");

  // Same requirement when the validator never rejects anything: the greedy stage
  // then runs its trial displacement all the way down to |disp| == 1 without ever
  // being told to stop. That shortest trial is still valid, so it is the answer —
  // abandoning it and falling back to the bare face displacement would put the
  // probe right back on the origin.
  {
    const Level n = 4;
    Radii<NDIM> box_radius, surface_thickness;
    for (std::size_t d = 0; d != NDIM; ++d) {
      box_radius[d] = 2;
      surface_thickness[d] = 1;
    }
    BoxSurfaceDisplacementRange<NDIM> range(centered_key<NDIM>(n), box_radius, surface_thickness,
                                            lattice_summed, accept_all<NDIM>());
    const auto probe = range.probing_displacement().translation();
    const auto period = Translation(1) << n;
    bool equivalent_to_origin = true;
    for (std::size_t d = 0; d != NDIM; ++d)
      if (mod_positive(probe[d], period) != 0) equivalent_to_origin = false;
    t.checkpoint(!equivalent_to_origin, "3D n=4 N={2,2,2} accept-all: probe is not lattice-equivalent to the origin");
  }

  return t.end();
}

/// The fallback block builds its displacement in *boxes*, not in units of the
/// kernel range.
///
/// `box_radius` is in half-simulation-cells; a displacement is in boxes, and the
/// two differ by 2^{n-1}. Using the raw radius yields a probe that is 2^{n-1}
/// times too short — an interior displacement rather than a boundary one, and
/// inconsistent with the backup-face component, which is scaled.
///
/// The configuration below (small level, small even radii) is the one that
/// reaches the fallback: no axis can be displaced past bmax_default while
/// staying inside the box, so the greedy stage declines every dimension.
int test_probe_fallback_is_in_box_units(World& world) {
  test_output t("BoxSurfaceDisplacementRange: fallback probe is in box units", world.rank() == 0);

  constexpr std::size_t NDIM = 3;
  const Level n = 2;
  const std::array<std::int64_t, NDIM> N = {2, 2, 2};

  Radii<NDIM> box_radius, surface_thickness;
  for (std::size_t d = 0; d != NDIM; ++d) {
    box_radius[d] = N[d];
    surface_thickness[d] = 0;
  }
  BoxSurfaceDisplacementRange<NDIM> range(centered_key<NDIM>(n), box_radius, surface_thickness,
                                          array_of_bools<NDIM>{true}, reject_all<NDIM>());
  const auto probe = range.probing_displacement().translation();

  // r = N * 2^{n-1} = 2*2 = 4 boxes; the raw radius would be 2
  for (std::size_t d = 0; d != NDIM; ++d) {
    const auto r = radius_in_boxes(N[d], n);
    t.checkpoint(probe[d] == r, "probe[" + std::to_string(d) + "] == " + std::to_string(r) + " boxes");
  }

  return t.end();
}

/// Every component of the fallback probe must respect *its own* axis' radius.
///
/// The running maximum across axes is the loop bound for the shrink stage, not
/// the value to assign: handing axis d the largest radius seen so far puts the
/// probe outside the box whenever a later axis is shorter than an earlier one,
/// and makes the result depend on iteration order.
///
/// N={2,4,2} at n=1 keeps r == N, so this isolates the per-axis question from
/// the units question tested above.
int test_probe_fallback_respects_each_radius(World& world) {
  test_output t("BoxSurfaceDisplacementRange: fallback probe respects each axis' radius", world.rank() == 0);

  constexpr std::size_t NDIM = 3;
  const Level n = 1;
  const std::array<std::int64_t, NDIM> N = {2, 4, 2};  // shorter axis *after* the longer one

  Radii<NDIM> box_radius, surface_thickness;
  for (std::size_t d = 0; d != NDIM; ++d) {
    box_radius[d] = N[d];
    surface_thickness[d] = 0;
  }
  BoxSurfaceDisplacementRange<NDIM> range(centered_key<NDIM>(n), box_radius, surface_thickness,
                                          array_of_bools<NDIM>{true}, reject_all<NDIM>());
  const auto probe = range.probing_displacement().translation();

  for (std::size_t d = 0; d != NDIM; ++d) {
    const auto r = radius_in_boxes(N[d], n);
    t.checkpoint(std::abs(probe[d]) <= r,
                 "probe[" + std::to_string(d) + "] = " + std::to_string(probe[d]) + " within its own radius " +
                     std::to_string(r));
  }

  return t.end();
}

/// An axis of unlimited extent has no surface, so it contributes nothing to the
/// probe — and must not be asked for a radius it does not have.
int test_probe_partial_range(World& world) {
  test_output t("BoxSurfaceDisplacementRange: unrestricted axis contributes no probe component",
                world.rank() == 0);

  constexpr std::size_t NDIM = 3;
  const Level n = 2;

  Radii<NDIM> box_radius, surface_thickness;
  box_radius[0] = 2;               surface_thickness[0] = 0;
  box_radius[1] = std::nullopt;    surface_thickness[1] = std::nullopt;  // unlimited along axis 1
  box_radius[2] = 2;               surface_thickness[2] = 0;

  BoxSurfaceDisplacementRange<NDIM> range(centered_key<NDIM>(n), box_radius, surface_thickness,
                                          array_of_bools<NDIM>{true}, reject_all<NDIM>());
  const auto probe = range.probing_displacement().translation();

  t.checkpoint(probe[1] == 0, "unrestricted axis 1 has a zero probe component");
  t.checkpoint(probe[0] == radius_in_boxes(2, n), "axis 0 probe is the range boundary");
  t.checkpoint(probe[2] == radius_in_boxes(2, n), "axis 2 probe is the range boundary");

  return t.end();
}

/// Constructing without a validator must work.
///
/// `Validator` is a std::function and the constructor defaults it to empty;
/// funcimpl.h leaves it empty whenever the standard-displacement loop screened
/// everything out. compute_probe() must not call through it.
int test_probe_without_validator(World& world) {
  test_output t("BoxSurfaceDisplacementRange: empty validator is not called", world.rank() == 0);

  constexpr std::size_t NDIM = 3;
  const Level n = 4;

  Radii<NDIM> box_radius, surface_thickness;
  for (std::size_t d = 0; d != NDIM; ++d) {
    box_radius[d] = 2;  // even + lattice-summed => the branch that wants a validator
    surface_thickness[d] = 1;
  }

  bool threw = false;
  Translation probe0 = 0;
  try {
    BoxSurfaceDisplacementRange<NDIM> range(centered_key<NDIM>(n), box_radius, surface_thickness,
                                            array_of_bools<NDIM>{true});
    probe0 = range.probing_displacement().translation()[0];
  } catch (...) {
    threw = true;
  }
  t.checkpoint(!threw, "construction with a default-constructed validator does not throw");
  t.checkpoint(probe0 == radius_in_boxes(2, n), "probe falls back to the plain face displacement");

  return t.end();
}

int main(int argc, char** argv) {
  World& world = madness::initialize(argc, argv);
  startup(world, argc, argv, true);

  int errors = 0;
  errors += test_surface_extent(world);
  errors += test_no_duplicates(world);
  errors += test_probe_hyperface_origin(world);
  errors += test_probe_even_lattice_sum(world);
  errors += test_probe_fallback_is_in_box_units(world);
  errors += test_probe_fallback_respects_each_radius(world);
  errors += test_probe_partial_range(world);
  errors += test_probe_without_validator(world);

  world.gop.fence();
  madness::finalize();
  return errors;
}
