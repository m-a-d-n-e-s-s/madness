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

/// full == home + images for the periodic axes given by `periodic`, lattice range N
int check_identity(World& world, const std::array<bool,3>& periodic, int N, const std::string& label) {
    BoundaryConditions<3> bc(BC_FREE);
    std::array<LatticeRange,3> lr_full, lr_home;
    int nper = 0;
    for (int d = 0; d < 3; ++d) {
        lr_home[d] = LatticeRange(0);
        if (periodic[d]) { bc(d,0) = bc(d,1) = BC_PERIODIC; lr_full[d] = LatticeRange(N); ++nper; }
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
    const int expected_rank = ((1 << nper) - 1) * full.get_rank();

    int errors = 0;
    if (world.rank() == 0)
        print(" ", label, ": |full - home - images| =", err, "  ranks: full", full.get_rank(),
              " images", images.get_rank(), "(expected", expected_rank, ")  |images f| =", vimg.norm2());
    if (err > 20.0 * thresh) { print("FAIL: identity violated"); ++errors; }
    if (images.get_rank() != expected_rank) { print("FAIL: rank"); ++errors; }
    // the images potential must be non-trivial: it is what the identity is for
    if (vimg.norm2() < 1.e-3) { print("FAIL: images operator is (near) zero"); ++errors; }
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
    errors += check_identity(world, {false, false, true}, 2, "periodic z, N=2   (1 pattern)");
    errors += check_identity(world, {true,  true,  true}, 1, "periodic xyz, N=1 (7 patterns)");

    if (world.rank() == 0) {
        if (errors == 0) print("\ntest_images_operator passed\n");
        else print("\ntest_images_operator FAILED with", errors, "errors\n");
    }
    world.gop.fence();
    madness::finalize();
    return errors;
}
