/// \file test_truncate_tol.cc
/// \brief the cell-width scaling of truncate_tol() in modes 1, 2 and 3

/// Modes 1-3 scale the truncation tolerance by the physical width of the box,
/// tol*min(1, 2^-level * L). For an anisotropic cell L is the geometric mean
/// width (volume^(1/NDIM)): the smallest dimension would tie the tolerance to
/// the cell's aspect ratio even though the resolution is set by the widest
/// dimension. Cubic cells are unaffected -- there the two agree.

#include <madness/mra/mra.h>
#include <madness/mra/funcdefaults.h>

using namespace madness;

namespace {

int check(const std::string& msg, double got, double expected, double tol = 1e-12) {
    const double err = std::abs(got - expected) / std::max(1.0, std::abs(expected));
    if (err > tol) {
        print("FAIL:", msg, " got", got, " expected", expected);
        return 1;
    }
    return 0;
}

/// the tolerance of an anisotropic cell scales with its geometric mean width
int test_anisotropic_cell(World& world) {
    int errors = 0;
    Tensor<double> cell(3, 2);
    cell(0, 0) = -200.0; cell(0, 1) = 200.0;     // 400
    cell(1, 0) = -200.0; cell(1, 1) = 200.0;     // 400
    cell(2, 0) = -20.0;  cell(2, 1) = 20.0;      //  40
    FunctionDefaults<3>::set_cell(cell);

    const double L = std::pow(400.0 * 400.0 * 40.0, 1.0 / 3.0);   // 185.66...
    errors += check("cell_min_width", FunctionDefaults<3>::get_cell_min_width(), 40.0);
    errors += check("cell_geometric_mean_width",
                    FunctionDefaults<3>::get_cell_geometric_mean_width(), L);

    real_function_3d f = real_factory_3d(world);     // empty; only its impl is needed
    auto impl = f.get_impl();
    const double tol = 1.e-6;
    const double fac3 = 1.0 / std::pow(2.0, 3 * 0.5);   // mode 3's sibling factor

    impl->set_truncate_mode(0);
    errors += check("mode 0 is unscaled",
                    impl->truncate_tol(tol, Key<3>(10, Vector<Translation,3>(0))), tol);

    impl->set_truncate_mode(1);
    for (int n : {0, 8, 10}) {
        const double expected = tol * std::min(1.0, std::pow(0.5, double(n)) * L);
        errors += check("mode 1, level " + std::to_string(n),
                        impl->truncate_tol(tol, Key<3>(n, Vector<Translation,3>(0))), expected);
    }

    impl->set_truncate_mode(2);
    for (int n : {0, 8}) {
        const double expected = tol * std::min(1.0, std::pow(0.25, double(n)) * L * L);
        errors += check("mode 2, level " + std::to_string(n),
                        impl->truncate_tol(tol, Key<3>(n, Vector<Translation,3>(0))), expected);
    }

    impl->set_truncate_mode(3);
    for (int n : {0, 8}) {
        const double expected = tol * fac3 * std::min(1.0, std::pow(0.5, double(n)) * L);
        errors += check("mode 3, level " + std::to_string(n),
                        impl->truncate_tol(tol, Key<3>(n, Vector<Translation,3>(0))), expected);
    }
    return errors;
}

/// a cubic cell is untouched: there the geometric mean is the minimum width
int test_cubic_cell_unchanged(World& world) {
    int errors = 0;
    FunctionDefaults<3>::set_cubic_cell(-50.0, 50.0);
    errors += check("cubic: min == geometric mean",
                    FunctionDefaults<3>::get_cell_geometric_mean_width(),
                    FunctionDefaults<3>::get_cell_min_width());

    real_function_3d f = real_factory_3d(world);
    auto impl = f.get_impl();
    impl->set_truncate_mode(1);
    const double tol = 1.e-6;
    errors += check("cubic cell, mode 1, level 8",
                    impl->truncate_tol(tol, Key<3>(8, Vector<Translation,3>(0))),
                    tol * std::min(1.0, std::pow(0.5, 8.0) * 100.0));
    return errors;
}

}   // namespace

int main(int argc, char** argv) {
    World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);

    int errors = 0;
    errors += test_anisotropic_cell(world);
    errors += test_cubic_cell_unchanged(world);

    if (world.rank() == 0) {
        if (errors == 0) print("\ntest_truncate_tol passed\n");
        else print("\ntest_truncate_tol FAILED with", errors, "errors\n");
    }
    world.gop.fence();
    madness::finalize();
    return errors;
}
