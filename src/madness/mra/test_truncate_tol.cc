/// \file test_truncate_tol.cc
/// \brief the cell-width scaling of truncate_tol() in modes 1, 2 and 3

/// Modes 1-3 scale the truncation tolerance by the physical width of the box,
/// tol*min(1, 2^-level * L). For an anisotropic cell L must be the geometric
/// mean width (volume^(1/NDIM)), not the smallest dimension: keying off the
/// smallest dimension tightens the tolerance by the cell's aspect ratio even
/// though the resolution is set by the widest dimension. Measured consequence
/// (MPQC, H10 chain in a 100x100x18 bohr cell): mode 1 cost 10-13x the time and
/// 8x the memory of mode 0, for LDA, PBE, TPSS and HF alike.

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

/// tolerance scales with the geometric mean width, and the old behaviour is still reachable
int test_anisotropic_cell(World& world) {
    int errors = 0;
    Tensor<double> cell(3, 2);
    cell(0, 0) = -200.0; cell(0, 1) = 200.0;     // 400
    cell(1, 0) = -200.0; cell(1, 1) = 200.0;     // 400
    cell(2, 0) = -20.0;  cell(2, 1) = 20.0;      //  40
    FunctionDefaults<3>::set_cell(cell);

    const double Lmin = 40.0;
    const double Lmean = std::pow(400.0 * 400.0 * 40.0, 1.0 / 3.0);   // 185.66...

    errors += check("cell_min_width", FunctionDefaults<3>::get_cell_min_width(), Lmin);
    errors += check("cell_geometric_mean_width",
                    FunctionDefaults<3>::get_cell_geometric_mean_width(), Lmean);
    errors += check("truncate width defaults to the geometric mean",
                    FunctionDefaults<3>::get_truncate_cell_width(), Lmean);

    real_function_3d f = real_factory_3d(world);     // empty; only its impl is needed
    auto impl = f.get_impl();
    const double tol = 1.e-6;
    const double fac3 = 1.0 / std::pow(2.0, 3 * 0.5);   // mode 3's sibling factor

    for (const bool by_min : {false, true}) {
        FunctionDefaults<3>::set_truncate_scale_by_min_width(by_min);
        const double L = by_min ? Lmin : Lmean;
        const std::string tag = by_min ? " [min]" : " [mean]";

        impl->set_truncate_mode(0);
        errors += check("mode 0 is unscaled" + tag, impl->truncate_tol(tol, Key<3>(10, Vector<Translation,3>(0))), tol);

        impl->set_truncate_mode(1);
        for (int n : {0, 8, 10}) {
            const double expected = tol * std::min(1.0, std::pow(0.5, double(n)) * L);
            errors += check("mode 1, level " + std::to_string(n) + tag,
                            impl->truncate_tol(tol, Key<3>(n, Vector<Translation,3>(0))), expected);
        }

        impl->set_truncate_mode(2);
        for (int n : {0, 8}) {
            const double expected = tol * std::min(1.0, std::pow(0.25, double(n)) * L * L);
            errors += check("mode 2, level " + std::to_string(n) + tag,
                            impl->truncate_tol(tol, Key<3>(n, Vector<Translation,3>(0))), expected);
        }

        impl->set_truncate_mode(3);
        for (int n : {0, 8}) {
            const double expected = tol * fac3 * std::min(1.0, std::pow(0.5, double(n)) * L);
            errors += check("mode 3, level " + std::to_string(n) + tag,
                            impl->truncate_tol(tol, Key<3>(n, Vector<Translation,3>(0))), expected);
        }
    }
    FunctionDefaults<3>::set_truncate_scale_by_min_width(false);
    return errors;
}

/// a cubic cell must be untouched by the change: there min == geometric mean
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
    for (const bool by_min : {false, true}) {
        FunctionDefaults<3>::set_truncate_scale_by_min_width(by_min);
        errors += check("cubic cell, mode 1, level 8",
                        impl->truncate_tol(tol, Key<3>(8, Vector<Translation,3>(0))),
                        tol * std::min(1.0, std::pow(0.5, 8.0) * 100.0));
    }
    FunctionDefaults<3>::set_truncate_scale_by_min_width(false);
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
