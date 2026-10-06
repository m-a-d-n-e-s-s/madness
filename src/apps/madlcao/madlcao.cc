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
/// \brief driver for closed-shell Hartree-Fock in a Gaussian basis with separated-kernel integrals

/// Reads the same input as moldft (geometry, `dft` group), so the molecule and
/// its orientation are identical; the basis set is the `aobasis` of the `dft`
/// group. LCAO-specific settings live in the `lcao` group.
///
///   madlcao --geometry=water --dft="aobasis 6-31g" [--lcao="guess core; check_mra true"]

#include <madness/madness_config.h>
#include <madness/chem/CalculationParameters.h>
#include <madness/chem/Restart.h>
#include <madness/chem/lcao_scf.h>
#include <madness/chem/molecular_functors.h>
#include <madness/chem/potentialmanager.h>
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/vmra.h>
#include <madness/world/parallel_archive.h>

#include <algorithm>
#include <cmath>
#include <cstdio>

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
    const long nocc = param.nalpha();
    std::vector<real_function_3d> mo = lcao::project_orbitals(world, molecule, aobasis, scf.coefficients(), nocc);
    const Tensor<double> S = matrix_inner(world, mo, mo, true);
    double maxdev = 0.0;
    for (long i = 0; i < nocc; ++i)
        for (long j = 0; j < nocc; ++j) maxdev = std::max(maxdev, std::abs(S(i, j) - (i == j ? 1.0 : 0.0)));
    mo = orthonormalize_symmetric(mo, S);

    const Tensor<double> eps = copy(scf.orbital_energies()(Slice(0, nocc - 1)));
    Tensor<double> occ(nocc);
    occ.fill(1.0);                              // one electron per spin orbital; moldft doubles for closed shells
    const std::vector<int> set(nocc, 0);

    RestartMetadata meta;
    meta.current_energy = scf.energies().total;
    meta.spin_restricted = true;
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
        ar & static_cast<unsigned int>(mo.size());
        ar & eps & occ & set;
        for (const real_function_3d& f : mo) ar & f;
    }
    world.gop.fence();

    if (world.rank() == 0) {
        printf("\nseed: wrote %ld occupied orbitals to %s (rung %d: thresh %.0e, k %d, %.1fs)\n",
               nocc, name.c_str(), rung, thresh, k, wall_time() - t0);
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
                print("madlcao: closed-shell Hartree-Fock in a Gaussian basis, with all integrals from");
                print("MADNESS's Gaussian fit of 1/r (prototype v1)\n");
                print("usage: madlcao --geometry=water --dft=\"aobasis 6-31g\" [--lcao=\"guess core\"]\n");
                print("the basis set is the aobasis keyword of the dft group; the lcao group holds:");
                LCAOParameters().print("lcao", "end");
            }
        } else {
            try {
                const Molecule molecule(world, parser);
                CalculationParameters param(world, parser);
                param.set_derived_values(molecule);
                const LCAOParameters lparam(world, parser);
                MADNESS_CHECK_THROW(param.nalpha() == param.nbeta(), "madlcao v1: closed-shell RHF only");

                AtomicBasisSet aobasis;
                aobasis.read_file(param.aobasis());
                if (world.rank() == 0) {
                    print("\n madlcao: closed-shell Hartree-Fock in a Gaussian basis (separated-kernel integrals)\n");
                    molecule.print();
                    lparam.print("lcao", "end");
                    print("\nbasis set", param.aobasis(), "with", aobasis.nbf(molecule), "functions,",
                          param.nalpha(), "doubly occupied orbitals\n");
                }

                lcao::LCAOSCF scf(world, molecule, aobasis, param.nalpha(), lparam);
                const double t0 = wall_time();
                const double energy = scf.solve();
                if (world.rank() == 0) printf("final energy=%16.8f  (%.2fs)\n", energy, wall_time() - t0);
                if (not scf.converged()) status = 1;

                if (lparam.check_eri()) check_eri(world, scf, lparam);
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
