#include "DFConvergence.h"
#include "DFParameters.h"
#include "DFRestart.h"
#include "InitParameters.h"
#include <madness/world/MADworld.h>
#include <madness/world/binary_fstream_archive.h>
#include <madness/world/parallel_archive.h>
#include <gtest/gtest.h>
#include <filesystem>
#include <sstream>
#include <string>
#include <unistd.h>

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

TEST(DFConvergence, StopMessagesAreUnambiguous) {
    using madness::DFStopReason;
    EXPECT_EQ(madness::df_stop_message(DFStopReason::converged_bsh),
              "Converged due to residuals");
    EXPECT_EQ(madness::df_stop_message(DFStopReason::converged_combined),
              "Converged due to energy, density, and residuals");
    EXPECT_EQ(madness::df_stop_message(DFStopReason::max_iterations),
              "WARNING: maximum iterations reached without convergence");
}

namespace {
madness::DFParameters read_parameters(const std::string& body) {
    std::istringstream input("DiracFock\n" + body + "\nend\n");
    madness::DFParameters parameters;
    parameters.read(input);
    return parameters;
}
}  // namespace

TEST(DFParameters, ExplicitThresholdsAreOrderIndependent) {
    const auto before = read_parameters(
        "dconv 2e-5\nthresh_mul 0.0\nthresh 1e-6");
    const auto after = read_parameters(
        "thresh 1e-6\ndconv 2e-5\nthresh_mul 0.0");
    EXPECT_DOUBLE_EQ(before.thresh, after.thresh);
    EXPECT_DOUBLE_EQ(before.dconv, 2.e-5);
    EXPECT_DOUBLE_EQ(after.dconv, 2.e-5);
    EXPECT_DOUBLE_EQ(before.thresh_mul, 0.0);
    EXPECT_DOUBLE_EQ(after.thresh_mul, 0.0);
}

TEST(DFParameters, UnspecifiedThresholdsFollowFinalThresh) {
    const auto parameters = read_parameters("thresh 3e-7");
    EXPECT_DOUBLE_EQ(parameters.dconv, 3.e-7);
    EXPECT_DOUBLE_EQ(parameters.thresh_mul, 3.e-7);
}

namespace {
madness::World* test_world = nullptr;

// The pid goes in front of the stem, not after it. clean_archive_filename()
// strips a trailing 5-digit ".NNNNN" chunk suffix. A pid at the end looks
// the same as that suffix.
std::string temp_archive_path(const char* stem) {
    const std::string unique = std::to_string(::getpid()) + "_" + stem;
    return (std::filesystem::temp_directory_path() / unique).string();
}

std::string collective_temp_archive_path(madness::World& world, const char* stem) {
    std::string path;
    if (world.rank() == 0) path = temp_archive_path(stem);
    world.gop.broadcast_serializable(path, 0);
    return path;
}

void write_legacy_df_header(const std::string& path, const double energy,
                            const bool krestricted, const bool closed_shell,
                            const unsigned int norbitals) {
    madness::archive::BinaryFstreamOutputArchive ar(path.c_str());
    madness::Tensor<double> energies(norbitals);
    ar & energy & krestricted & closed_shell & norbitals & energies;
}

void write_legacy_parallel_df_archive(madness::World& world,
                                      const std::string& path) {
    madness::archive::ParallelOutputArchive<> ar(world, path, 1);
    const double energy = -14.5;
    const bool krestricted = false;
    const bool closed_shell = true;
    const unsigned int norbitals = 0;
    const madness::Tensor<double> energies(norbitals);
    const double box_size = 20.0;
    const int wavelet_order = 8;
    const madness::Molecule molecule;
    ar & energy & krestricted & closed_shell & norbitals & energies
       & box_size & wavelet_order & molecule;
    ar.flush();
}
}  // namespace

TEST(DFRestart, VersionRoundTrips) {
    const std::string path = temp_archive_path("df_restart_version.bin");
    {
        madness::archive::BinaryFstreamOutputArchive ar(path.c_str());
        madness::write_df_restart_version(ar);
    }
    unsigned int version = 0;
    {
        madness::archive::BinaryFstreamInputArchive ar(path.c_str());
        version = madness::read_df_restart_version(ar);
    }
    EXPECT_EQ(version, madness::DF_RESTART_VERSION);
    EXPECT_NO_THROW(madness::require_supported_df_restart(version));
    std::filesystem::remove(path);
}

TEST(DFRestart, LegacyArchiveReadsAsUnversioned) {
    const std::string path = temp_archive_path("df_restart_legacy.bin");
    write_legacy_df_header(path, -14.5, false, true, 2);
    unsigned int version = madness::DF_RESTART_VERSION;
    {
        madness::archive::BinaryFstreamInputArchive ar(path.c_str());
        version = madness::read_df_restart_version(ar);
    }
    EXPECT_EQ(version, 0u);
    std::filesystem::remove(path);
}

TEST(DFRestart, UnsupportedVersionsExplainRecovery) {
    for (const unsigned int bad : {0u, madness::DF_RESTART_VERSION + 1u}) {
        try {
            madness::require_supported_df_restart(bad);
            FAIL() << "accepted DF restart version " << bad;
        } catch (const madness::MadnessException& error) {
            const std::string message = error.what();
            EXPECT_NE(message.find("DF restart archive"), std::string::npos);
            EXPECT_NE(message.find("start from a moldft archive"),
                      std::string::npos);
        }
    }
}

TEST(DFRestart, LegacyArchiveIsRejectedCollectively) {
    ASSERT_NE(test_world, nullptr);
    madness::World& world = *test_world;
    const std::string path =
        collective_temp_archive_path(world, "df_restart_collective_legacy");
    write_legacy_parallel_df_archive(world, path);
    world.gop.fence();

    madness::InitParameters parameters;
    EXPECT_THROW(parameters.read(world, path, 137.03599917697017,
                                 true, false),
                 madness::MadnessException);

    world.gop.fence();
    madness::archive::ParallelInputArchive<>::remove(world, path.c_str());
    world.gop.fence();
}

int main(int argc, char** argv) {
    madness::World& world = madness::initialize(argc, argv);
    test_world = &world;
    ::testing::InitGoogleTest(&argc, argv);
    const int status = RUN_ALL_TESTS();
    test_world = nullptr;
    madness::finalize();
    return status;
}
