/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

/// \file madlcao.cc
/// \brief driver for Hartree-Fock (RHF or UHF) in a Gaussian basis with separated-kernel integrals

/// Reads the same input as moldft (geometry, `dft` group), so the molecule, its
/// orientation and the spin (nopen) are identical. The basis set and the other
/// LCAO settings live in the `lcao` group, as for moldft's `guess lcao`.
///
///   madlcao --geometry=water --lcao="basis 6-31g; guess core; check_mra true"

#include <madness/madness_config.h>
#include <madness/chem/CalculationParameters.h>
#include <madness/chem/Restart.h>
#include <madness/chem/lcao_cholesky.h>
#include <madness/chem/lcao_scf.h>
#include <madness/chem/molecular_functors.h>
#include <madness/chem/potentialmanager.h>
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/vmra.h>
#include <madness/tensor/cblas.h>
#include <madness/world/parallel_archive.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>

using namespace madness;

namespace {

/// compare the separated-kernel integrals with MRA quadrature of the same basis functions

/// An independent check of the 1D machinery, the component order and the
/// normalization convention. The nuclear attraction differs by MADNESS's
/// smoothing of the nuclear potential, which the LCAO calculation does not use.
void check_against_mra(World& world, const Molecule& molecule, const AtomicBasisSet& aobasis,
                       const lcao::LCAOSCF& scf, const double L) {
    const double thresh = 1.e-6;
    FunctionDefaults<3>::set_cubic_cell(-L, L);
    FunctionDefaults<3>::set_k(8);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_initial_level(2);
    FunctionDefaults<3>::set_truncate_mode(1);

    const double t0 = wall_time();
    const int nbf = aobasis.nbf(molecule);
    std::vector<real_function_3d> ao(nbf);
    for (int i = 0; i < nbf; ++i) {
        // the base-class pointer selects the right FunctionFactory::functor overload
        const std::shared_ptr<FunctionFunctorInterface<double,3>> f =
            std::make_shared<madchem::AtomicBasisFunctor>(aobasis.get_atomic_basis_function(molecule, i));
        ao[i] = real_factory_3d(world).functor(f).truncate_on_project().nofence();
    }
    world.gop.fence();

    const Tensor<double> S = matrix_inner(world, ao, ao, true);

    Tensor<double> T(nbf, nbf);
    const auto grad = gradient_operator<double,3>(world);
    for (int axis = 0; axis < 3; ++axis) {
        const std::vector<real_function_3d> d = apply(world, *grad[axis], ao);
        T += matrix_inner(world, d, d, true);
    }
    T.scale(0.5);

    PotentialManager pm(molecule, molecule.parameters.core_type());
    pm.make_nuclear_potential(world);
    const real_function_3d vnuc = pm.vnuclear();
    const std::vector<real_function_3d> vao = mul(world, vnuc, ao);
    Tensor<double> V = matrix_inner(world, ao, vao);
    V = 0.5 * (V + transpose(V));

    const std::vector<real_function_3d> rho = square(world, ao);
    const real_convolution_3d coulop = CoulombOperator(world, 1.e-6, 0.1 * thresh);
    const std::vector<real_function_3d> vrho = apply(world, coulop, rho);
    const Tensor<double> G = matrix_inner(world, rho, vrho, true);
    Tensor<double> Glcao(nbf, nbf);
    for (int i = 0; i < nbf; ++i)
        for (int j = 0; j < nbf; ++j) Glcao(i, j) = scf.eri()(i, i, j, j);

    if (world.rank() == 0) {
        printf("\nseparated-kernel integrals against MRA quadrature of the same functions (k=8, thresh %.0e, %.1fs)\n",
               thresh, wall_time() - t0);
        printf("   overlap             max |diff| %.2e   (max |S| %.2e)\n", (S - scf.overlap()).absmax(), S.absmax());
        printf("   kinetic             max |diff| %.2e   (max |T| %.2e)\n", (T - scf.kinetic()).absmax(), T.absmax());
        printf("   nuclear attraction  max |diff| %.2e   (max |V| %.2e; MRA smooths the nuclei, LCAO does not)\n",
               (V - scf.nuclear_attraction()).absmax(), V.absmax());
        printf("   (ii|jj)             max |diff| %.2e   (max (ii|jj) %.2e)\n", (G - Glcao).absmax(), G.absmax());
    }
}

/// the screened nuclear attraction of the SCF against every kernel term, nucleus and primitive pair computed
void check_onee(World& world, const Molecule& molecule, const lcao::LCAOSCF& scf, const LCAOParameters& lparam) {
    const lcao::GaussianKernel kernel =
        lcao::GaussianKernel::coulomb(lparam.kernel_lo(), lparam.kernel_hi(), lparam.kernel_eps());
    const lcao::SeparatedGaussianIntegrals ints(scf.shells(), kernel);
    const double t0 = wall_time();
    const Tensor<double> V = ints.nuclear_attraction(world, molecule, 0.0);
    if (world.rank() == 0)
        printf("\nnuclear attraction, screened (as in the SCF) against every term (%.2fs): largest |diff| %.2e   "
               "(largest |V| %.2e)\n", wall_time() - t0, (scf.nuclear_attraction() - V).absmax(), V.absmax());
}

/// compare the two-electron integrals with the unoptimized reference implementation (v1's loop)
void check_eri(World& world, const lcao::LCAOSCF& scf, const LCAOParameters& lparam) {
    const lcao::GaussianKernel kernel =
        lcao::GaussianKernel::coulomb(lparam.kernel_lo(), lparam.kernel_hi(), lparam.kernel_eps());
    const lcao::SeparatedGaussianIntegrals ints(scf.shells(), kernel);
    const double t0 = wall_time();
    const Tensor<double> G = ints.eri_reference();
    const lcao::PackedERI& P = scf.eri();
    double maxdiff = 0.0;
    for (long i = 0; i < P.nbf(); ++i)
        for (long j = 0; j <= i; ++j)
            for (long k = 0; k <= i; ++k)
                for (long l = 0; l <= k; ++l)
                    if (lcao::PackedERI::pair(k, l) <= lcao::PackedERI::pair(i, j))
                        maxdiff = std::max(maxdiff, std::abs(G(i, j, k, l) - P(i, j, k, l)));
    if (world.rank() == 0) {
        printf("\ntwo-electron integrals against the reference implementation (%.2fs)\n", wall_time() - t0);
        printf("   max |diff| %.2e   (max |(ij|kl)| %.2e)\n", maxdiff, G.absmax());
    }

    // The column engine of the Cholesky decomposition against the stored integrals. It keeps
    // every kernel term and skips nothing, so the stored integrals must be unscreened too.
    if (lparam.kernel_screen() != 0.0 or lparam.schwarz() != 0.0) {
        if (world.rank() == 0) print("\ncolumn engine: not checked (needs kernel_screen 0; schwarz 0)");
        return;
    }
    const double t1 = wall_time();
    const lcao::FunctionPairs fp = ints.function_pairs();
    const long n = P.nbf();
    std::vector<int> seen(n * n, 0);
    for (const auto& f : fp.functions) ++seen[std::max(f[0], f[1]) * n + std::min(f[0], f[1])];
    bool once = true;
    for (long i = 0; i < n; ++i)
        for (long j = 0; j <= i; ++j) once = once and seen[i * n + j] == 1;

    const Tensor<double> D = ints.pair_diagonal(world, fp);
    double ddiff = 0.0;
    for (std::size_t r = 0; r < fp.size(); ++r) {
        const auto& f = fp.functions[r];
        ddiff = std::max(ddiff, std::abs(D(long(r)) - P(f[0], f[1], f[0], f[1])));
    }

    // eri() computed (ab|cd) with the larger group pair as the bra, so the columns of a
    // smaller ket group pair must agree bitwise; the others only to rounding
    std::vector<std::size_t> all(fp.ngroup_pairs());
    std::iota(all.begin(), all.end(), 0);
    double cdiff = 0.0, cdiff_same = 0.0;
    for (std::size_t g = 0; g < fp.ngroup_pairs(); ++g) {
        const Tensor<double> W = ints.eri_columns(world, fp, g, all);
        for (std::size_t j = 0; j < fp.count[g]; ++j) {
            const auto& c = fp.functions[fp.first[g] + j];
            for (std::size_t r = 0; r < fp.size(); ++r) {
                const auto& f = fp.functions[r];
                const double d = std::abs(W(long(j), long(r)) - P(f[0], f[1], c[0], c[1]));
                cdiff = std::max(cdiff, d);
                if (fp.group_pair[r] > g) cdiff_same = std::max(cdiff_same, d);
            }
        }
    }
    if (world.rank() == 0) {
        printf("\ncolumn engine against the stored integrals (%.2fs)\n", wall_time() - t1);
        printf("   function pairs %zu in %zu group pairs, each unordered pair once: %s\n", fp.size(),
               fp.ngroup_pairs(), once ? "yes" : "NO");
        printf("   diagonal: largest |diff| %.2e   columns: largest |diff| %.2e (bra group pair > ket, as eri() "
               "computed them: %.2e)\n", ddiff, cdiff, cdiff_same);
    }
}

/// the Cholesky decomposition of the two-electron integrals at cholesky_tol against the stored integrals
///
/// The bound |V - L L^T| <= cholesky_tol holds for V without kernel screening and Schwarz
/// skips, which the decomposition does not use, so the stored integrals must be unscreened too.
void check_cholesky(World& world, const lcao::LCAOSCF& scf, const LCAOParameters& lparam) {
    if (lparam.kernel_screen() != 0.0 or lparam.schwarz() != 0.0) {
        if (world.rank() == 0) print("\nCholesky decomposition: not checked (needs kernel_screen 0; schwarz 0)");
        return;
    }
    const lcao::GaussianKernel kernel =
        lcao::GaussianKernel::coulomb(lparam.kernel_lo(), lparam.kernel_hi(), lparam.kernel_eps());
    const lcao::SeparatedGaussianIntegrals ints(scf.shells(), kernel);
    const double t0 = wall_time();
    const auto chol = std::make_shared<const lcao::CholeskyERIDecomposition>(world, ints, lparam.cholesky_tol());
    const double t1 = wall_time();

    // L L^T on the kept rows, one dgemm; the screened rows have no vector components
    const lcao::FunctionPairs& fp = chol->pairs();
    const long nk = long(chol->nkept()), m = chol->nvec();
    Tensor<double> LLT(std::max(nk, 1L), std::max(nk, 1L));
    if (m > 0)
        cblas::gemm(cblas::NoTrans, cblas::Trans, nk, nk, m, 1.0, chol->vectors(), nk, chol->vectors(), nk, 0.0,
                    LLT.ptr(), nk);
    std::vector<long> kept(fp.size(), -1);
    for (std::size_t i = 0; i < chol->rows().size(); ++i) kept[chol->rows()[i]] = long(i);
    const lcao::PackedERI& P = scf.eri();
    double maxerr = 0.0;
    for (std::size_t r = 0; r < fp.size(); ++r) {
        const auto& a = fp.functions[r];
        for (std::size_t c = 0; c <= r; ++c) {
            const auto& b = fp.functions[c];
            const double llt = (kept[r] >= 0 and kept[c] >= 0) ? LLT(kept[r], kept[c]) : 0.0;
            maxerr = std::max(maxerr, std::abs(P(a[0], a[1], b[0], b[1]) - llt));
        }
    }

    // J and K of the converged densities from the vectors against those from the stored integrals
    const lcao::InCoreERI incore(world, std::shared_ptr<const lcao::PackedERI>(&P, [](const lcao::PackedERI*) {}));
    const lcao::CholeskyERI fromvectors(world, chol, P.nbf());
    const bool open = not scf.restricted();
    const Tensor<double> Pa = scf.density(0);
    const Tensor<double> Pb = open ? scf.density(1) : Tensor<double>();
    Tensor<double> J0, Ka0, Kb0, J1, Ka1, Kb1;
    incore.jk(Pa, Pb, J0, Ka0, Kb0);
    const double t2 = wall_time();
    fromvectors.jk(Pa, Pb, J1, Ka1, Kb1);
    const double t3 = wall_time();

    if (world.rank() == 0) {
        const lcao::CholeskyERIDecomposition::Stats& s = chol->stats();
        printf("\nCholesky decomposition of the two-electron integrals (cholesky_tol %.0e, span %.0e, %.2fs)\n",
               chol->tol(), chol->span(), t1 - t0);
        printf("   %zu of %zu function pairs kept, %ld vectors = %.2f N, %zu integral columns in %zu batches\n",
               s.nkept, s.npairs, m, double(m) / double(P.nbf()), s.ncolumns, s.nbatches);
        printf("   time: diagonal %.2fs, integrals %.2fs, updates %.2fs\n", s.t_diagonal, s.t_integrals, s.t_updates);
        printf("   largest |V - L L^T| %.2e: %s\n", maxerr,
               maxerr <= chol->tol() ? "within cholesky_tol" : "EXCEEDS cholesky_tol");
        printf("   J and K of the converged density, vectors against stored integrals (%.3fs): largest |dJ| %.2e, "
               "|dKa| %.2e, |dKb| %.2e\n", t3 - t2, (J1 - J0).absmax(), (Ka1 - Ka0).absmax(),
               open ? (Kb1 - Kb0).absmax() : 0.0);
    }
}

/// the polynomial order moldft uses for a threshold, unless k is given (SCF::set_protocol)
int k_for_thresh(const double thresh) {
    if (thresh >= 0.9e-2) return 4;
    if (thresh >= 0.9e-4) return 6;
    if (thresh >= 0.9e-6) return 8;
    if (thresh >= 0.9e-8) return 10;
    return 12;
}


/// write the occupied LCAO orbitals as <prefix>.restartdata, for moldft to start from

/// The orbitals are projected at the thresh and k of rung seed_rung of the dft
/// group's protocol and Loewdin-orthonormalized. The header is moldft's own
/// (RestartMetadata), so `restart auto` reads the archive like any other. It
/// claims convergence at the previous rung, so moldft starts at seed_rung, and
/// leaves the density convergence unknown, so moldft always iterates.
void write_seed(World& world, const Molecule& molecule, const AtomicBasisSet& aobasis, const lcao::LCAOSCF& scf,
                const CalculationParameters& param, const LCAOParameters& lparam) {
    const std::vector<double> protocol = param.protocol();
    const int rung = lparam.seed_rung();
    MADNESS_CHECK_THROW(rung >= 0 and rung < int(protocol.size()), "seed_rung is outside the dft group's protocol");
    const double thresh = protocol[rung];
    const int k = (param.k() > 0) ? param.k() : k_for_thresh(thresh);

    // the function defaults moldft uses for this rung (SCF ctor and SCF::set_protocol)
    FunctionDefaults<3>::set_cubic_cell(-param.L(), param.L());
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_initial_level(2);
    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_apply_randomize(false);
    FunctionDefaults<3>::set_project_randomize(false);
    FunctionDefaults<3>::set_truncate_mode(1);

    const double t0 = wall_time();
    // the occupied orbitals of both spins in one projection, so the basis functions are
    // projected once; then each spin is Loewdin-orthonormalized on its own. maxdev is the
    // largest deviation of the projections from orthonormality.
    const bool unrestricted = not scf.restricted();
    const long na = scf.nalpha(), nb = unrestricted ? scf.nbeta() : 0;
    Tensor<double> C(scf.nbf(), na + nb);
    C(_, Slice(0, na - 1)) = scf.coefficients(0)(_, Slice(0, na - 1));
    if (nb > 0) C(_, Slice(na, na + nb - 1)) = scf.coefficients(1)(_, Slice(0, nb - 1));
    const std::vector<real_function_3d> all = lcao::project_orbitals(world, molecule, aobasis, C, na + nb);
    double maxdev = 0.0;
    const auto orthonormalize = [&](const std::vector<real_function_3d>& mo) {
        const Tensor<double> S = matrix_inner(world, mo, mo, true);
        for (long i = 0; i < long(mo.size()); ++i)
            for (long j = 0; j < long(mo.size()); ++j) maxdev = std::max(maxdev, std::abs(S(i, j) - (i == j ? 1.0 : 0.0)));
        return orthonormalize_symmetric(mo, S);
    };
    const std::vector<real_function_3d> amo = orthonormalize({all.begin(), all.begin() + na});
    const std::vector<real_function_3d> bmo =
        nb > 0 ? orthonormalize({all.begin() + na, all.end()}) : std::vector<real_function_3d>();

    // moldft's format (SCF::save_mos): per spin the number of orbitals, their energies,
    // occupations (one electron per spin orbital) and localization sets, then the functions
    const auto write_block = [&](auto& ar, const int spin, const std::vector<real_function_3d>& mo) {
        const long nocc = mo.size();
        const Tensor<double> eps = copy(scf.orbital_energies(spin)(Slice(0, nocc - 1)));
        Tensor<double> occ(nocc);
        occ.fill(1.0);
        const std::vector<int> set(nocc, 0);
        ar & static_cast<unsigned int>(nocc);
        ar & eps & occ & set;
        for (const real_function_3d& f : mo) ar & f;
    };

    RestartMetadata meta;
    meta.current_energy = scf.energies().total;
    meta.spin_restricted = not unrestricted;
    meta.L = param.L();
    meta.k = k;
    meta.molecule = molecule;
    meta.xc = param.xc();
    meta.localize = param.localize_method();
    meta.converged_for_thresh = (rung > 0) ? protocol[rung - 1] : 1.e10;
    meta.converged_for_dconv = 1.e10;
    meta.representation = Representation::mo;
    meta.eprec = molecule.parameters.eprec();
    meta.madness_version = MADNESS_PACKAGE_VERSION;

    const std::string name = param.prefix() + ".restartdata";
    {
        archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(world, name.c_str(),
                                                                                param.get<int>("nio"));
        meta.write(ar);
        write_block(ar, 0, amo);
        if (not bmo.empty()) write_block(ar, 1, bmo);
    }
    world.gop.fence();

    if (world.rank() == 0) {
        printf("\nseed: wrote %zu alpha and %zu beta occupied orbitals to %s (rung %d: thresh %.0e, k %d, %.1fs)\n",
               amo.size(), unrestricted ? bmo.size() : amo.size(), name.c_str(), rung, thresh, k, wall_time() - t0);
        printf("      max |S_ij - delta_ij| of the projections before Loewdin: %.2e\n", maxdev);
    }
}

} // namespace


int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    int status = 0;
    {
        startup(world, argc, argv);
        commandlineparser parser(argc, argv);

        if (parser.key_exists("help")) {
            if (world.rank() == 0) {
                print("madlcao: Hartree-Fock (RHF or UHF) in a Gaussian basis, with all integrals from");
                print("MADNESS's Gaussian fit of 1/r (prototype v1)\n");
                print("usage: madlcao --geometry=water --lcao=\"basis 6-31g; guess core\"\n");
                print("the lcao group holds:");
                LCAOParameters().print("lcao", "end");
            }
        } else {
            try {
                const Molecule molecule(world, parser);
                CalculationParameters param(world, parser);
                param.set_derived_values(molecule);
                const LCAOParameters lparam(world, parser);
                MADNESS_CHECK_THROW(lparam.eri() == "incore" or
                                    not (lparam.check_eri() or lparam.check_mra() or lparam.check_cholesky()),
                                    "check_eri, check_mra and check_cholesky compare with the stored integrals: "
                                    "they need eri incore");

                AtomicBasisSet aobasis;
                aobasis.read_file(lparam.basis());
                if (world.rank() == 0) {
                    print("\n madlcao: Hartree-Fock (RHF or UHF) in a Gaussian basis (separated-kernel integrals)\n");
                    molecule.print();
                    lparam.print("lcao", "end");
                    print("\nbasis set", lparam.basis(), "with", aobasis.nbf(molecule), "functions,",
                          param.nalpha(), "alpha and", param.nbeta(), "beta electrons\n");
                }

                lcao::LCAOSCF scf(world, molecule, aobasis, param.nalpha(), param.nbeta(), lparam);
                const double t0 = wall_time();
                const double energy = scf.solve();
                if (world.rank() == 0) printf("final energy=%16.8f  (%.2fs)\n", energy, wall_time() - t0);
                if (not scf.converged()) status = 1;

                if (lparam.check_onee()) check_onee(world, molecule, scf, lparam);
                if (lparam.check_eri()) check_eri(world, scf, lparam);
                if (lparam.check_cholesky()) check_cholesky(world, scf, lparam);
                if (lparam.check_mra()) check_against_mra(world, molecule, aobasis, scf, param.L());
                if (lparam.seed()) write_seed(world, molecule, aobasis, scf, param, lparam);

            } catch (const madness::MadnessException& e) {
                print(e);
                status = 1;
            } catch (const std::exception& e) {
                print("madlcao failed:", e.what());
                status = 1;
            }
        }
        world.gop.fence();
    }
    finalize();
    return status;
}
