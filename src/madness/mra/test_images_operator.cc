/// \file test_images_operator.cc
/// \brief the rest-of-crystal operator: a lattice sum that excludes the home cell

/// A lattice-summed separable kernel sums, along each periodic axis, over the
/// images R = -maxR..maxR of every 1D Gaussian factor. The rest-of-crystal
/// operator is the same sum with the single lattice vector L = 0 removed, so
/// that   full = home + images   exactly, term by term, with the same fit.
///
/// Removing L = 0 is not "drop R = 0 on every axis": that would also drop
/// every L with any zero component, and for a chain periodic along z the only
/// images at all are (0, 0, n_z). The operator is instead a sum over the 2^p - 1
/// nonzero patterns of the p periodic axes, each a product of 1D factors that
/// sum either R != 0 or R = 0 only, so its rank is (2^p - 1) times the fit's.
///
/// Why it exists: MPQC forms the periodic energy correction as V - V^0, two full
/// applies on the cuspy total charge density whose difference is ~1e-11 for an
/// isolated system, so the input cannot be loosened (measured: 10x looser on
/// one side moves the energy 3e-6). Built directly, the images operator is
/// smooth and long-range, cancels nothing, and lets the home-cell potential be
/// split by density and cached.

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/funcdefaults.h>

using namespace madness;

namespace {

double f_func(const coord_3d& r) {
    const double a = (r[0]-0.5)*(r[0]-0.5) + (r[1]+0.3)*(r[1]+0.3) + (r[2]-1.2)*(r[2]-1.2);
    const double b = (r[0]+1.0)*(r[0]+1.0) + (r[1]-0.7)*(r[1]-0.7) + (r[2]+2.0)*(r[2]+2.0);
    return std::exp(-1.5*a) + 0.7*std::exp(-0.8*b);
}

/// full == home + images for the periodic axes given by `periodic`, lattice range `N`.
/// `home` is built from the full fit; an infinite `N` drops the low-exponent tail of
/// the fit from `full`, and the images operator carries that tail as negated
/// home-only terms, so the identity holds term by term in that case too.
int check_identity(World& world, const std::array<bool,3>& periodic, const LatticeRange& N, const std::string& label) {
    BoundaryConditions<3> bc(BC_FREE);
    std::array<LatticeRange,3> lr_full, lr_home;
    int nper = 0;
    for (int d = 0; d < 3; ++d) {
        lr_home[d] = LatticeRange(0);
        if (periodic[d]) { bc(d,0) = bc(d,1) = BC_PERIODIC; lr_full[d] = N; ++nper; }
        else lr_full[d] = LatticeRange(0);
    }
    FunctionDefaults<3>::set_bc(bc);
    Displacements<3>().reset_periodic_axes(array_of_bools<3>(periodic[0], periodic[1], periodic[2]));

    const double thresh = FunctionDefaults<3>::get_thresh();
    real_function_3d f = real_factory_3d(world).f(f_func);
    f.truncate();

    OperatorInfo info(0.0, 1.e-4, thresh, OT_G12);
    real_convolution_3d full(world, info, lr_full);
    real_convolution_3d home(world, info, lr_home);
    info.images_only = true;
    real_convolution_3d images(world, info, lr_full);

    real_function_3d vfull = full(f), vhome = home(f), vimg = images(f);
    const double err = (vfull - vhome - vimg).norm2();
    // the images operator carries the truncated fit's patterns plus, for an infinite sum, one
    // home-only term per Gaussian by which the truncated fit differs from the full one: at least
    // the dropped terms, at most every term if the truncation also rescales kept coefficients
    const int npat = (1 << nper) - 1;
    const int min_rank = npat * full.get_rank() + (home.get_rank() - full.get_rank());
    const int max_rank = npat * full.get_rank() + home.get_rank();

    int errors = 0;
    if (world.rank() == 0)
        print(" ", label, ": |full - home - images| =", err, "  ranks: full", full.get_rank(),
              " home", home.get_rank(), " images", images.get_rank(), "(expected", min_rank, "..", max_rank, ")  |images f| =", vimg.norm2());
    if (err > 20.0 * thresh) { print("FAIL: identity violated"); ++errors; }
    if (images.get_rank() < min_rank || images.get_rank() > max_rank) { print("FAIL: rank"); ++errors; }
    // the images potential must be non-trivial: it is what the identity is for
    if (vimg.norm2() < 1.e-3) { print("FAIL: images operator is (near) zero"); ++errors; }
    return errors;
}

/// A compact source near a periodic face of a non-cubic cell. The apply stops sweeping
/// displacements at the first shell that contributes nothing, which assumes a kernel that
/// decays away from the source; the images kernel instead vanishes near the source and turns
/// on near the lattice images, and in a 100x100x18 cell the empty shells along the short axis
/// ended the sweep before the wrapped displacement carrying the nearest-image interaction.
/// Term 42 of the fit lost 4% of its images potential, the total 2.5e-3, independent of the
/// threshold. SeparatedConvolution::screen_by_shell_decay() opts the images operator out.
double g_z0 = 7.2;
double near_face_f(const coord_3d& r) {
    const double a = 100.0, z = r[2] - g_z0;
    return std::pow(a/constants::pi, 1.5)*std::exp(-a*(r[0]*r[0] + r[1]*r[1] + z*z));
}
int check_near_face(World& world) {
    Tensor<double> cell(3,2);
    cell(0,0) = -50; cell(0,1) = 50; cell(1,0) = -50; cell(1,1) = 50; cell(2,0) = -9; cell(2,1) = 9;
    FunctionDefaults<3>::set_cell(cell);
    BoundaryConditions<3> bc(BC_FREE); bc(2,0) = bc(2,1) = BC_PERIODIC;
    FunctionDefaults<3>::set_bc(bc);
    Displacements<3>().reset_periodic_axes(array_of_bools<3>(false, false, true));
    const double thresh = FunctionDefaults<3>::get_thresh();

    std::array<LatticeRange,3> lr_full{LatticeRange(0), LatticeRange(0), LatticeRange(10)};
    std::array<LatticeRange,3> lr_home{LatticeRange(0), LatticeRange(0), LatticeRange(0)};
    OperatorInfo info(0.0, 1.e-4, thresh, OT_G12, false);
    real_convolution_3d full(world, info, lr_full), home(world, info, lr_home);
    info.images_only = true;
    real_convolution_3d images(world, info, lr_full);

    const std::vector<coord_3d> sp{coord_3d{0.0, 0.0, g_z0}};
    real_function_3d f = real_factory_3d(world).f(near_face_f).special_points(sp).special_level(8);
    f.truncate();
    real_function_3d vfull = full(f), vhome = home(f), vimg = images(f);
    const double err = (vfull - vhome - vimg).norm2();
    int errors = 0;
    if (world.rank() == 0)
        print("  100x100x18, source 1.8 bohr from the periodic face: |full - home - images| =", err,
              "  |images f| =", vimg.norm2());
    // the residual left is the displacement-reach floor shared by full and home in this cell
    // (each ~1e-4 from the exact potential at thresh 1e-6); before the fix it was 2.5e-3
    if (err > 100.0 * thresh) { print("FAIL: identity violated near the periodic face"); ++errors; }
    return errors;
}

}   // namespace

int main(int argc, char** argv) {
    World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);

    FunctionDefaults<3>::set_cubic_cell(-10.0, 10.0);
    FunctionDefaults<3>::set_k(8);
    FunctionDefaults<3>::set_thresh(1.e-5);
    FunctionDefaults<3>::set_truncate_mode(0);

    int errors = 0;
    errors += check_identity(world, {false, false, true}, LatticeRange(2),    "periodic z, N=2   (1 pattern)");
    errors += check_identity(world, {true,  true,  true}, LatticeRange(1),    "periodic xyz, N=1 (7 patterns)");
    errors += check_identity(world, {false, false, true}, LatticeRange(true), "periodic z, N=inf (1 pattern + dropped tail)");
    errors += check_identity(world, {true,  true,  true}, LatticeRange(true), "periodic xyz, N=inf (7 patterns + dropped tail)");
    FunctionDefaults<3>::set_thresh(1.e-6);
    errors += check_near_face(world);

    if (world.rank() == 0) {
        if (errors == 0) print("\ntest_images_operator passed\n");
        else print("\ntest_images_operator FAILED with", errors, "errors\n");
    }
    world.gop.fence();
    madness::finalize();
    return errors;
}
