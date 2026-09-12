#include "DFConvergence.h"
#include "DFParameters.h"
#include <madness/world/MADworld.h>
#include <gtest/gtest.h>

using madness::DFConvergenceCriterion;
using madness::DFConvergenceMetrics;
using madness::df_iteration_converged;

TEST(DFConvergence, BshResidualUsesOnlyStrictResidual) {
    const DFConvergenceMetrics m{-75.0, -74.0, 1.e-6,
                                 1.e-2, 1.e-6, 9.e-7, 1.e-6};
    EXPECT_TRUE(df_iteration_converged(DFConvergenceCriterion::bsh_residual, m));

    auto at_threshold = m;
    at_threshold.max_bsh_residual = 1.e-6;
    EXPECT_TRUE(df_iteration_converged(
        DFConvergenceCriterion::bsh_residual, at_threshold));

    auto relaxed_only = m;
    relaxed_only.max_bsh_residual = 5.e-5;
    EXPECT_FALSE(df_iteration_converged(
        DFConvergenceCriterion::bsh_residual, relaxed_only));
}

TEST(DFConvergence, CombinedCriterionRejectsUnconvergedEnergy) {
    const DFConvergenceMetrics m{-75.0, -74.0, 1.e-6,
                                 1.e-8, 1.e-6, 5.e-5, 1.e-6};
    EXPECT_FALSE(df_iteration_converged(
        DFConvergenceCriterion::energy_density_residual, m));
}

TEST(DFConvergence, CombinedCriterionRequiresAllThreeMeasurements) {
    const auto combined = DFConvergenceCriterion::energy_density_residual;
    const DFConvergenceMetrics converged{-75.0, -75.0 + 1.e-7, 1.e-6,
                                         1.e-8, 1.e-6, 5.e-5, 1.e-6};
    EXPECT_TRUE(df_iteration_converged(combined, converged));

    auto bad_density = converged;
    bad_density.density_residual = 2.e-6;
    EXPECT_FALSE(df_iteration_converged(combined, bad_density));

    auto bad_bsh = converged;
    bad_bsh.max_bsh_residual = 2.e-4;
    EXPECT_FALSE(df_iteration_converged(combined, bad_bsh));
}

TEST(DFConvergence, CombinedCriterionUsesAbsoluteEnergyChangeAtZeroEnergy) {
    const auto combined = DFConvergenceCriterion::energy_density_residual;
    const DFConvergenceMetrics converged{0.0, 5.e-7, 1.e-6,
                                         1.e-8, 1.e-6, 5.e-5, 1.e-6};
    EXPECT_TRUE(df_iteration_converged(combined, converged));

    auto bad_energy = converged;
    bad_energy.previous_total_energy = 2.e-6;
    EXPECT_FALSE(df_iteration_converged(combined, bad_energy));
}

int main(int argc, char** argv) {
    madness::initialize(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    const int status = RUN_ALL_TESTS();
    madness::finalize();
    return status;
}
