/// \file test_periodic_truncation.cc
/// \brief the truncated fit of an infinite lattice sum must stay accurate along the finite axes

/// A Gaussian of the kernel fit with exponent below 0.25/L^2 has a lattice sum along a
/// periodic axis of width L that is flat over the cell -- a gauge constant -- but the
/// same Gaussian is not flat along an axis that is not summed. Dropping such terms
/// outright is exact for a cell periodic in every direction and wrong otherwise, by an
/// amount that grows with the transverse cell size and does not shrink with the
/// threshold (a few 1e-3 in the potential below for the 100x100x18 cell, at eps 1e-6..1e-8).
/// GFit::truncate_mixed_expansion keeps those terms as far as the finite axes need them.
///
/// The test: spherical, neutral charge blobs on a chain periodic along z. Their images
/// produce no potential at all (a neutral sphere has none outside itself), so the
/// lattice-summed potential must equal the free-space one up to the gauge constant.

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

using namespace madness;

namespace {

constexpr int natoms = 5;
constexpr double spacing = 3.2, a_sharp = 300.0, a_diffuse = 6.0;
inline double center_z(int i) { return (i - 0.5*(natoms - 1))*spacing; }

double rho(const coord_3d& r) {
    double v = 0.0;
    for (int i = 0; i < natoms; ++i) {
        const double z = r[2] - center_z(i);
        const double rr = r[0]*r[0] + r[1]*r[1] + z*z;
        v += std::pow(a_sharp/constants::pi, 1.5)*std::exp(-a_sharp*rr) - std::pow(a_diffuse/constants::pi, 1.5)*std::exp(-a_diffuse*rr);
    }
    return v;
}

inline double erf_over_r(double a, double r) {   // free-space potential of the unit Gaussian (a/pi)^{3/2} exp(-a r^2)
    return (r < 1e-8) ? 2.0*std::sqrt(a/constants::pi) : std::erf(std::sqrt(a)*r)/r;
}

double v_free(const coord_3d& r) {
    double v = 0.0;
    for (int i = 0; i < natoms; ++i) {
        const double z = r[2] - center_z(i);
        const double rr = std::sqrt(r[0]*r[0] + r[1]*r[1] + z*z);
        v += erf_over_r(a_sharp, rr) - erf_over_r(a_diffuse, rr);
    }
    return v;
}

int check(World& world, double Lxy, double Lz, double bound_factor, const std::string& label) {
    Tensor<double> cell(3,2);
    cell(0,0) = -Lxy/2; cell(0,1) = Lxy/2; cell(1,0) = -Lxy/2; cell(1,1) = Lxy/2; cell(2,0) = -Lz/2; cell(2,1) = Lz/2;
    FunctionDefaults<3>::set_cell(cell);
    // keep the displacement ordering in step with a non-cubic cell; a no-op where set_cell
    // already does this (master since 46333af6b), and without it the apply, not the
    // truncation, is what fails in the 100x100x18 cell
    Displacements<3>::set_width(FunctionDefaults<3>::get_cell_width());
    BoundaryConditions<3> bc(BC_FREE); bc(2,0) = bc(2,1) = BC_PERIODIC;
    FunctionDefaults<3>::set_bc(bc);
    Displacements<3>().reset_periodic_axes(array_of_bools<3>(false, false, true));
    const double thresh = FunctionDefaults<3>::get_thresh();

    std::vector<coord_3d> centers;
    for (int i = 0; i < natoms; ++i) centers.push_back(coord_3d{0.0, 0.0, center_z(i)});
    real_function_3d f = real_factory_3d(world).f(rho).special_points(centers).special_level(6);
    f.truncate();
    real_function_3d vex = real_factory_3d(world).f(v_free).special_points(centers).special_level(6);

    std::array<LatticeRange,3> lr{LatticeRange(0), LatticeRange(0), LatticeRange(true)};
    OperatorInfo info(0.0, 1.e-4, thresh, OT_G12);   // infinite lattice sum along z: the fit is truncated
    real_convolution_3d coulomb(world, info, lr);
    real_function_3d v = coulomb(f);
    real_function_3d d = v - vex;
    const double avg = d.trace()/FunctionDefaults<3>::get_cell_volume();
    const double err = (d - avg).norm2();
    int errors = 0;
    if (world.rank() == 0)
        print(" ", label, ": rank", coulomb.get_rank(), " f size", f.size(), " |V_inf - V_free - avg| =", err, " (avg", avg, ")");
    if (err > bound_factor * thresh) { print("FAIL: the truncated lattice sum is off along the finite axes"); ++errors; }
    return errors;
}

}   // namespace

int main(int argc, char** argv) {
    World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);

    FunctionDefaults<3>::set_k(8);
    FunctionDefaults<3>::set_thresh(1.e-6);
    FunctionDefaults<3>::set_truncate_mode(0);

    // The bounds leave room for the apply's own floor, which in the 5.6:1 cell is set by the
    // isotropic displacement reach (~6e-5 here; ~1e-5 with the per-axis reach of PR #819) --
    // the old truncation is two orders of magnitude above either.
    int errors = 0;
    errors += check(world, 18,  18, 30.0,  "18^3, periodic z        ");
    errors += check(world, 100, 18, 100.0, "100x100x18, periodic z  ");

    if (world.rank() == 0) {
        if (errors == 0) print("\ntest_periodic_truncation passed\n");
        else print("\ntest_periodic_truncation FAILED with", errors, "errors\n");
    }
    world.gop.fence();
    madness::finalize();
    return errors;
}
