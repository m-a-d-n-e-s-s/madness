/// \file test_anisotropic_reach.cc
/// \brief the free-space Coulomb potential in a non-cubic cell against the analytic result

/// FunctionImpl::do_apply visits, at level n, the displacements with |l_d| <= bmax on every
/// axis and applies through them the blocks of the non-standard form that involve a level-n
/// wavelet. For a Gaussian fit term those blocks are O(1) along an axis where the term is
/// unresolved at level n (exponent a with a h^2 >~ 1, h the box width) and fall off along a
/// resolved axis only as the (s|s) matrix element, ~ exp(-a (l h)^2). In a non-cubic cell a
/// term unresolved on the long axis is resolved on a short one, where h is smaller, and its
/// block reaches |l| ~ (long box)/(short box) there -- beyond an isotropic bound of 4 boxes.
/// No coarser level holds that block (a level-n wavelet is orthogonal to coarser scaling
/// functions), so the bound must be isotropic in real space (Displacements::bmax_axes).
/// Measured here at eps 1e-6: 5.8e-5 in the 100x100x18 cell with an isotropic bound of
/// 4 boxes, 1.2e-5 with (4, 4, 23), 5e-6 in the 18^3 cell.

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

using namespace madness;

namespace {

// three neutral, spherical charge blobs along z: a sharp positive and a diffuse negative Gaussian each
static double a_sharp = 100.0, a_diffuse = 2.0, spacing = 1.6; static int natoms = 3;
inline double center_z(int i) { return (i - 0.5*(natoms-1))*spacing; }

double rho(const coord_3d& r) {
    double v = 0.0;
    for (int i = 0; i < natoms; ++i) {
        const double z = r[2] - center_z(i);
        const double rr = r[0]*r[0] + r[1]*r[1] + z*z;
        v += std::pow(a_sharp/constants::pi, 1.5)*std::exp(-a_sharp*rr) - std::pow(a_diffuse/constants::pi, 1.5)*std::exp(-a_diffuse*rr);
    }
    return v;
}

inline double erf_over_r(double a, double r) {   // Coulomb potential of the unit Gaussian (a/pi)^{3/2} exp(-a r^2)
    return (r < 1e-8) ? 2.0*std::sqrt(a/constants::pi) : std::erf(std::sqrt(a)*r)/r;
}

double v_exact(const coord_3d& r) {
    double v = 0.0;
    for (int i = 0; i < natoms; ++i) {
        const double z = r[2] - center_z(i);
        const double rr = std::sqrt(r[0]*r[0] + r[1]*r[1] + z*z);
        v += erf_over_r(a_sharp, rr) - erf_over_r(a_diffuse, rr);
    }
    return v;
}

int check(World& world, double Lxy, double Lz, int nat, double sp, double ash, double adf, const std::string& label) {
    natoms = nat; spacing = sp; a_sharp = ash; a_diffuse = adf;
    Tensor<double> cell(3,2);
    cell(0,0) = -Lxy/2; cell(0,1) = Lxy/2; cell(1,0) = -Lxy/2; cell(1,1) = Lxy/2; cell(2,0) = -Lz/2; cell(2,1) = Lz/2;
    FunctionDefaults<3>::set_cell(cell);
    const double thresh = FunctionDefaults<3>::get_thresh();

    std::vector<coord_3d> centers;
    for (int i = 0; i < natoms; ++i) centers.push_back(coord_3d{0.0, 0.0, center_z(i)});
    real_function_3d f = real_factory_3d(world).f(rho).special_points(centers).special_level(6);
    f.truncate();
    real_function_3d vex = real_factory_3d(world).f(v_exact).special_points(centers).special_level(6);

    real_convolution_3d coulomb = CoulombOperator(world, 1.e-4, thresh);
    const double t0 = wall_time();
    real_function_3d v = coulomb(f);
    const double t_apply = wall_time() - t0;
    const double err = (v - vex).norm2();
    const auto b = Displacements<3>::bmax_current();
    if (world.rank() == 0)
        print(" ", label, ": reach", b[0], b[1], b[2], " f size", f.size(), " |V - exact| =", err, " |exact| =", vex.norm2(), " apply", t_apply, "s");
    int errors = 0;
    if (err > 25.0 * thresh) { print("FAIL: Coulomb potential off the analytic result"); ++errors; }
    return errors;
}

}   // namespace

int main(int argc, char** argv) {
    World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);

    FunctionDefaults<3>::set_k(8);
    FunctionDefaults<3>::set_thresh(1.e-6);
    FunctionDefaults<3>::set_truncate_mode(0);
    FunctionDefaults<3>::set_bc(BoundaryConditions<3>(BC_FREE));

    // five neutral blobs 3.2 bohr apart along z, spanning 13 of the 18 bohr of the short axis;
    // the diffuse parts (width 0.4) are 2.6 bohr from the faces, so the density is not cut by the cell
    int errors = 0;
    errors += check(world, 18,  18, 5, 3.2, 300.0, 6.0, "cubic 18^3       ");
    errors += check(world, 100, 18, 5, 3.2, 300.0, 6.0, "100x100x18 (5.6:1)");

    if (world.rank() == 0) {
        if (errors == 0) print("\ntest_anisotropic_reach passed\n");
        else print("\ntest_anisotropic_reach FAILED with", errors, "errors\n");
    }
    world.gop.fence();
    madness::finalize();
    return errors;
}
