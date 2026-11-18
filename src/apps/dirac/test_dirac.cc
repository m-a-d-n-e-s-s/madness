#include "DFConvergence.h"
#include "DFKramers.h"
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

TEST(DFKramers, QuartetDegeneracyPreservesKramersOrthogonality) {
    using cplx = std::complex<double>;
    // 3 Kramers pairs (m = 6): a p_{1/2} doublet (-0.854) and a 4-fold degenerate
    // p_{3/2} quartet (-0.848), matching Neon's 2p shell structure.
    const long n = 3;
    const long m = 2 * n;

    madness::Tensor<cplx> S11(n, n);
    S11(0, 0) = 1.0;               S11(0, 1) = cplx(0.03, -0.02); S11(0, 2) = cplx(-0.01, 0.04);
    S11(1, 0) = std::conj(S11(0, 1)); S11(1, 1) = 1.0;            S11(1, 2) = cplx(0.05, 0.01);
    S11(2, 0) = std::conj(S11(0, 2)); S11(2, 1) = std::conj(S11(1, 2)); S11(2, 2) = 1.0;

    madness::Tensor<cplx> S21_noisy(n, n);
    S21_noisy(0, 0) = cplx(1e-6, -2e-6); // numerical noise on diagonal
    S21_noisy(0, 1) = cplx(0.04, 0.03);
    S21_noisy(1, 0) = -S21_noisy(0, 1) + cplx(2e-6, 1e-6);
    S21_noisy(1, 2) = cplx(-0.06, 0.05); // non-zero <K*psi_1 | psi_2> overlap!
    S21_noisy(2, 1) = -S21_noisy(1, 2) + cplx(-1e-6, 3e-6);

    madness::Tensor<cplx> S_unperm =
        madness::assemble_kramers_hermitian_matrix(S11, S21_noisy, true);

    // Check Lowdin S^{-1/2} preserves Kramers block structure and gives Q^H S Q = I
    madness::Tensor<cplx> Q =
        madness::symplectify_unpermuted_hermitian(madness::hermitian_inv_sqrt(S_unperm));
    madness::Tensor<cplx> prod = madness::inner(madness::conj_transpose(Q), madness::inner(S_unperm, Q));
    for (long i = 0; i < m; ++i) {
        for (long j = 0; j < m; ++j) {
            const double target = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(std::abs(prod(i, j) - target), 0.0, 1.0e-13);
        }
    }

    // Build interleaved overlap and a scrambled p_{1/2} + p_{3/2} quaternionic Fock matrix
    madness::Tensor<double> P(m, m);
    for (long j = 0; j < n; ++j) {
        P(j, 2 * j) = 1.0;
        P(n + j, 2 * j + 1) = 1.0;
    }
    madness::Tensor<cplx> S_inter = madness::inner(madness::transpose(P), madness::inner(S_unperm, P));

    madness::Tensor<cplx> D_diag(m, m);
    D_diag(0, 0) = -0.854; D_diag(1, 1) = -0.854;
    D_diag(2, 2) = -0.848; D_diag(3, 3) = -0.848;
    D_diag(4, 4) = -0.848; D_diag(5, 5) = -0.848;

    // Scramble in an arbitrary symplectic gauge mixed with S_inter
    madness::Tensor<cplx> R_raw(m, m);
    R_raw(0, 0) = cplx(0.95, 0.1);  R_raw(1, 0) = cplx(0.05, -0.08);
    R_raw(2, 2) = cplx(0.6, 0.3);   R_raw(3, 2) = cplx(-0.2, 0.4);
    R_raw(4, 2) = cplx(0.3, -0.4);  R_raw(5, 2) = cplx(0.25, 0.15);
    R_raw(2, 4) = cplx(-0.35, 0.2); R_raw(3, 4) = cplx(0.3, 0.25);
    R_raw(4, 4) = cplx(0.7, -0.1);  R_raw(5, 4) = cplx(-0.15, 0.35);
    madness::Tensor<cplx> R_symp = madness::symplectify_interleaved_matrix(R_raw);
    madness::Tensor<cplx> VH(m, m), W(m, m);
    madness::Tensor<double> sig(m);
    madness::svd(R_symp, W, sig, VH);
    madness::Tensor<cplx> G = madness::symplectify_interleaved_matrix(madness::inner(W, VH));

    // V_exact = S^{-1/2} * G satisfies V_exact^H * S * V_exact = I, so
    // F = S * V_exact * D_diag * V_exact^H * S has exact generalized eigenvalues D_diag wrt S
    madness::Tensor<cplx> Q_inter = madness::inner(madness::transpose(P), madness::inner(Q, P));
    madness::Tensor<cplx> V_exact = madness::inner(Q_inter, G);
    madness::Tensor<cplx> F_inter = madness::inner(
        S_inter, madness::inner(V_exact, madness::inner(D_diag, madness::inner(madness::conj_transpose(V_exact), S_inter))));
    F_inter = 0.5 * (F_inter + madness::conj_transpose(F_inter));

    madness::Tensor<cplx> U;
    madness::Tensor<double> evals;
    madness::sygv(F_inter, S_inter, 1, U, evals);
    madness::kramers_pair_and_align_eigenvectors(F_inter, S_inter, U, evals, 1.0e-5, false);

    // Check exact partner degeneracy and eigenvalue recovery
    EXPECT_NEAR(evals(0L), -0.854, 1.0e-13);
    EXPECT_DOUBLE_EQ(evals(0L), evals(1L));
    for (long c = 2; c < 6; ++c) {
        EXPECT_NEAR(evals(c), -0.848, 1.0e-13);
    }
    EXPECT_DOUBLE_EQ(evals(2L), evals(3L));
    EXPECT_DOUBLE_EQ(evals(4L), evals(5L));

    // Check U^H * S_inter * U == I and U(:, 2k+1) == K_vec(U(:, 2k))
    madness::Tensor<cplx> US = madness::inner(madness::conj_transpose(U), madness::inner(S_inter, U));
    for (long i = 0; i < m; ++i) {
        for (long j = 0; j < m; ++j) {
            const double target = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(std::abs(US(i, j) - target), 0.0, 1.0e-13);
        }
    }
    for (long k = 0; k < n; ++k) {
        madness::Tensor<cplx> uk = madness::copy(U(madness::_, 2 * k));
        madness::Tensor<cplx> uk_bar = madness::copy(U(madness::_, 2 * k + 1));
        madness::Tensor<cplx> expected_bar = madness::kramers_interleaved_partner(uk);
        for (long r = 0; r < m; ++r) {
            EXPECT_NEAR(std::abs(uk_bar(r) - expected_bar(r)), 0.0, 1.0e-14);
        }
    }
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
