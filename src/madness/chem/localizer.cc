//
// Created by Florian Bischoff on 11/1/21.
//

#include<madness/chem/localizer.h>
#include <madness/mra/mra.h>
#include <madness/mra/function_factory.h>
#include <madness/tensor/jacobi.h>

using namespace madness;

namespace madness {

Localizer::Localizer(World& world, const AtomicBasisSet& aobasis, const Molecule& molecule,
                     const std::vector<Function<double, 3>>& ao) : aobasis(aobasis), molecule(molecule), ao(ao) {
    aobasis.atoms_to_bfn(molecule, at_to_bf, at_nbf);
    thresh_degenerate = FunctionDefaults<3>::get_thresh();
}

template<typename T, std::size_t NDIM>
MolecularOrbitals<T, NDIM> Localizer::localize(const MolecularOrbitals<T, NDIM>& mo_in, bool randomize) const {

    Tensor<T> Fock ;
    return localize(mo_in, Fock, randomize);
}

template<typename T, std::size_t NDIM>
MolecularOrbitals<T, NDIM> Localizer::localize(const MolecularOrbitals<T, NDIM>& mo_in, const Tensor<T>& Fock,
                                               bool randomize) const {

    World& world = mo_in.get_mos()[0].world();
    long nmo = mo_in.get_mos().size();
    Tensor<T> U(nmo, nmo);
    Tensor<T> fnew(nmo, nmo);
    if (enforce_core_valence_separation) {
        MADNESS_CHECK(size_t(Fock.dim(0)) == mo_in.get_mos().size());

        MolecularOrbitals<T, NDIM> mo_cv_separated = mo_in;
        Tensor<double> U1 = compute_core_valence_separation_transformation_matrix(world, mo_in, Fock);
        mo_cv_separated.update_mos(transform(world, mo_in.get_mos(), U1));

        auto slices = MolecularOrbitals<T, NDIM>::convert_set_to_slice(mo_in.get_localize_sets());

        // localize all orbital sets separately
        Tensor<T> UT2(nmo, nmo);
        for (int set = 0; set <= mo_cv_separated.get_localize_sets().back(); ++set) {
            Slice s = slices[set];
            MolecularOrbitals<T, NDIM> mo1 = mo_cv_separated.get_subset(set);
            Tensor<T> block_UT = compute_localization_matrix(world, mo1, randomize);
            UT2(s, s) = block_UT;
        }
        U = inner(U1, transpose(UT2));
        fnew = inner(U, inner(Fock, U, 1, 0), 0, 0);
    } else {
        U = transpose(compute_localization_matrix(world, mo_in, randomize));
    }

    // transform orbitals
    std::vector<Function<T, NDIM>> result = truncate(transform(world, mo_in.get_mos(), U));

    // transform orbital energies if given
    Tensor<double> eps(mo_in.get_mos().size());
    if (fnew.has_data()) for (int i = 0; i < fnew.dim(0); ++i) eps(i) = fnew(i, i);

    MolecularOrbitals<T, NDIM> lmo(result, eps, {}, Tensor<double>(), mo_in.get_localize_sets());
    lmo.set_all_orbitals_occupied();
    return lmo;

}


template<typename T, std::size_t NDIM>
Tensor<T> Localizer::compute_localization_matrix(World& world, const MolecularOrbitals<T, NDIM>& mo_in,
                                                 bool randomize) const {
    // localize using the reconstructed orbitals
    std::vector<Function<T, NDIM>> psi = metric.is_initialized() ? metric * mo_in.get_mos() : mo_in.get_mos();

    DistributedMatrix<T> dUT;
    if (method == "pm") {
        dUT = localize_PM(world, psi, mo_in.get_localize_sets(), tolloc, randomize, false);
    } else if (method == "boys") {
        dUT = localize_boys(world, psi, mo_in.get_localize_sets(), tolloc, randomize);
    } else if (method == "new") {
        dUT = localize_new(world, psi, mo_in.get_localize_sets(), tolloc, randomize, false);
    } else {
        print("unknown localization method", method);
        MADNESS_EXCEPTION("unknown localization method", 1);
    }
    long nmo = mo_in.get_mos().size();
    Tensor<T> UT(nmo, nmo);
    dUT.copy_to_replicated(UT);
    return UT;
}

template<typename T, std::size_t NDIM>
MolecularOrbitals<T, NDIM> Localizer::separate_core_valence(const MolecularOrbitals<T, NDIM>& mo_in,
                                                            const Tensor<T>& Fock) const {
    World& world = mo_in.get_mos()[0].world();

    Tensor<double> UT = compute_core_valence_separation_transformation_matrix(world, mo_in, Fock);

    std::vector<Function<T, NDIM>> result = transform(world, mo_in.get_mos(), UT);
    Tensor<T> fnew = inner(UT, inner(Fock, UT, 1, 0), 0, 0);
    Tensor<double> eps(mo_in.get_mos().size());
    for (int i = 0; i < eps.size(); ++i) eps(i) = std::real(fnew(i, i));
    MolecularOrbitals<T, NDIM> mo_out(result, eps, {}, Tensor<double>(), mo_in.get_localize_sets());
    mo_out.set_all_orbitals_occupied();
    return mo_out;

}

template<typename T, std::size_t NDIM>
Tensor<T> Localizer::compute_core_valence_separation_transformation_matrix(World& world,
                         const MolecularOrbitals<T, NDIM>& mo_in, const Tensor<T>& Fock) const {

    // block diagonalize the fock matrix
    Tensor<T> U(Fock.dim(0),Fock.dim(1));
    Tensor<T> f=copy(Fock);
    jacobi(f,U,mo_in.get_localize_sets());
    U=transpose(U);

    world.gop.broadcast(U.ptr(), U.size(), 0);

    // compute new fock matrix
    Tensor<T> fnew = inner(U, inner(Fock, U, 1, 0), 0, 0);
    bool success = check_core_valence_separation(fnew, mo_in.get_localize_sets());
    if (not success) {
        print("fnew in compute_core_valence_separation_transformation_matrix");
        print(fnew);
    }
    MADNESS_CHECK(success);
    return U;
}

template<typename T>
bool Localizer::check_core_valence_separation(const Tensor<T>& Fock, const std::vector<int>& localized_set,
                                              const bool silent) {
    std::vector<Slice> slices = MolecularOrbitals<T, 3>::convert_set_to_slice(localized_set);
    Tensor<T> F = copy(Fock);
    for (auto s: slices) F(s, s) = 0.0;
    double error = F.absmax();
    bool success = (error < FunctionDefaults<3>::get_thresh() * 5.0);
    if (silent) return success;
    if (not success) {
        print("faulty localization: core-valence separation requested but Fock matrix not block diagonal");
        print("error norm", error);
        print("max error", F.absmax());
        print(Fock);
    }
    return success;
}


template<typename T, std::size_t NDIM>
DistributedMatrix<T> Localizer::localize_boys(World& world, const std::vector<Function<T, NDIM>>& mo,
                                              const std::vector<int>& set, const double thresh1,
                                              const bool randomize,
                                              const bool doprint) const {
    typedef Tensor<T> tensorT;
    typedef Function<T, NDIM> functionT;
    typedef FunctionFactory<T, NDIM> factoryT;

    double vtol = 1.e-7;
    double thresh = thresh1;
    if (thresh < 1e-6) thresh = 1e-6; //<<<<<<<<<<<<<<<<<<<<< need to implement new line search like in pm routine

    long nmo = mo.size();
    tensorT dip(nmo, nmo, 3);
    for (int axis = 0; axis < 3; ++axis) {
        auto dipole = [&axis](const Vector<double, NDIM>& x) { return x[axis]; };
        functionT fdip = factoryT(world).functor(dipole).initial_level(4);
        dip(_, _, axis) = matrix_inner(world, mo, mul_sparse(world, fdip, mo, vtol), true);
    }

    tensorT U(nmo, nmo);
    default_random_generator.setstate(
            182041 + world.rank() * 10101); // To help with reproducibility for debugging, etc.
    if (world.rank() == 0) {
        for (long i = 0; i < nmo; ++i)
            U(i, i) = 1.0;

        tensorT xprev; // previous search direction
        tensorT gprev; // previous gradient
        bool rprev = true; // if true previous iteration restricted step or did incomplete search (so don't do conjugate)
        const int N = (nmo * (nmo - 1)) / 2;
        for (long iter = 0; iter < 1200; ++iter) {
            tensorT g(nmo, nmo);
            double W = 0.0;
            // cannot restrict size of individual gradients if want to do line search --- should instead modify line search direction
            for (long i = 0; i < nmo; ++i) {
                W += DIP(dip, i, i, i, i);
                for (long j = 0; j < i; ++j) {
                    g(j, i) = (DIP(dip, i, i, i, j) - DIP(dip, j, j, j, i));
                    if (randomize && iter == 0) g(j, i) += 0.1 * (RandomValue<double>() - 0.5);
                    g(i, j) = -g(j, i);
                }
            }
            double maxg = g.absmax();
            if (doprint)
                printf("iteration %ld W=%.8f maxg=%.2e\n", iter, W, maxg);
            if (maxg < thresh) break;

            // construct search direction using conjugate gradient approach
            tensorT x = copy(g);
            if (!rprev) { // Only apply conjugacy if did LS with real gradient
                double gamma = g.trace(g - gprev) / gprev.trace(gprev);
                if (doprint) print("gamma", gamma);
                x.gaxpy(1.0, xprev, gamma);
            }

            // Perform the line search.
            rprev = false;
            double dxgrad = x.trace(g) * 2.0;
            if (dxgrad < 0 || ((iter + 1) % N) == 0) {
                if (doprint) print("resetting since dxgrad -ve or due to dimension", dxgrad, iter, N);
                x = copy(g);
                dxgrad = x.trace(g) * 2.0; // 2*2 = 4 which should be prefactor on integrals in gradient
            }
            xprev = x; // Save for next iteration, noting shallow copy
            gprev = g;

            double mu = 0.01 / std::max(0.1, maxg); // Restrict intial step mu by size of max gradient
            tensorT dU = matrix_exponential(x * mu);
            tensorT newdip = inner(dU, dip, 0, 1); // can optimize this since only want (ii|ii)
            newdip = inner(dU, newdip, 0, 1);
            double newW = 0.0;
            for (long i = 0; i < nmo; ++i) {
                newW += DIP(newdip, i, i, i, i);
            }

            if (randomize && iter == 0) {
                rprev = true; // since did not use real gradient
            } else { // perform quadratic fit using f(0), df(0)/dx=dxgrad, f(mu)
                double f0 = W;
                double f1 = newW;
                double hess = 2.0 * (f1 - f0 - mu * dxgrad) / (mu * mu);
                if (hess >= 0) {
                    if (doprint) print("+ve hessian", hess);
                    hess = -2.0 * dxgrad; // force a bigish step to get out of bad region
                    rprev = true; // since did not do line search
                }
                double mu2 = -dxgrad / hess;
                if (mu2 * maxg > 0.5) {
                    mu2 = 0.5 / maxg; // pi/6 = 0.524, pi/4=0.785
                    rprev = true; // since did not do line search
                }
                double f2p = f0 + dxgrad * mu2 + 0.5 * hess * mu2 * mu2;
                if (doprint) print(f0, f1, f2p, "dxg", dxgrad, "hess", hess, "mu", mu, "mu2", mu2);
                mu = mu2;
            }

            dU = matrix_exponential(x * mu);
            U = inner(U, dU, 1, 0);
            dip = inner(dU, dip, 0, 1);
            dip = inner(dU, dip, 0, 1);
        }

        bool switched = true;
        while (switched) {
            switched = false;
            for (int i = 0; i < nmo; i++) {
                for (int j = i + 1; j < nmo; j++) {
                    if (set[i] == set[j]) {
                        double sold = U(i, i) * U(i, i) + U(j, j) * U(j, j);
                        double snew = U(i, j) * U(i, j) + U(j, i) * U(j, i);
                        if (snew > sold) {
                            tensorT tmp = copy(U(_, i));
                            U(_, i) = U(_, j);
                            U(_, j) = tmp;
                            switched = true;
                        }
                    }
                }
            }
        }

        // Fix phases.
        for (long i = 0; i < nmo; ++i) {
            if (U(i, i) < 0.0)
                U(_, i).scale(-1.0);
        }

    }

    world.gop.broadcast(U.ptr(), U.size(), 0);

    DistributedMatrix<double> dUT = column_distributed_matrix<double>(world, nmo, nmo);
    dUT.copy_from_replicated(transpose(U));

    return dUT;
}


template<typename T, std::size_t NDIM>
DistributedMatrix<T> Localizer::localize_PM(World& world, const std::vector<Function<T, NDIM>>& mo,
                                            const std::vector<int>& set, const double thresh,
                                            const bool randomize, const bool doprint) const {

    DistributedMatrix<T> dUT = distributed_localize_PM(world, mo, ao, set, at_to_bf, at_nbf,
                                                       thresh, thetamax, randomize, doprint);
    return dUT;
}


template<typename T, std::size_t NDIM>
DistributedMatrix<T> Localizer::localize_new(World& world, const std::vector<Function<T, NDIM>>& mo,
                                             const std::vector<int>& set, double thresh,
                                             const bool randomize, const bool doprint) const {
    // PROFILE_MEMBER_FUNC(SCF);
    typedef Tensor<T> tensorT;

    int nmo = mo.size();
    int nao = ao.size();

    tensorT C = matrix_inner(world, mo, ao);
    std::vector<int> at_to_bf, at_nbf; // OVERRIDE DATA IN CLASS OBJ TO USE ATOMS OR SHELLS FOR TESTING

    bool use_atomic_evecs = true;
    if (use_atomic_evecs) {
        // Transform from AOs to orthonormal atomic eigenfunctions
        int ilo = 0;
        for (size_t iat = 0; iat < molecule.natom(); ++iat) {
            const tensorT& avec = aobasis.get_avec(molecule, iat);
            int ihi = ilo + avec.dim(1);
            Slice s(ilo, ihi - 1);
            C(_, s) = inner(C(_, s), avec);

            // This commented out code would put each shell (s, p, d, etc.) into a separate set, however this
            // led to delocalized valence orbitals due to the inability to localize hybrid orbitals (e.g., sp3)
            // into a single set.  The initial workaround was to put the 1s core orbital into one set, and
            // everything else into a second.  However, this was not optimal for heavier atoms.  So now for
            // atoms with Z>=12 (Mg and beyond) we separate out the 2s+p sets.

            // // generate shell dimensions for atomic eigenfunctions
            // // ... this relies upon spherical symmetry being enforced
            // // when making atomic states
            // const tensorT& aeps = aobasis.get_aeps(molecule, iat);
            // //print(aeps);
            // double prev = aeps(0L);
            // int start = 0;
            // int i; // used after loop
            // for (i = 0; i < aeps.dim(0); ++i) {
            //     //print(" ... ", i, prev, aeps(i), (std::abs(aeps(i)-prev) > 1e-2*std::abs(prev)));
            //     if (std::abs(aeps(i) - prev) > 1e-2 * std::abs(prev)) {
            //         at_to_bf.push_back(ilo + start);
            //         at_nbf.push_back(i - start);
            //         //print("    ", start, i-start);
            //         start = i;
            //     }
            //     prev = aeps(i);
            // }
            // at_to_bf.push_back(ilo + start);
            // at_nbf.push_back(i - start);
            // //print("    ", start, i-start);

            at_to_bf.push_back(ilo);
            at_nbf.push_back(1); // 1s core orbital on atom
            if (molecule.get_atomic_number(iat) < 12) {
                // For lighter atoms everything else put into one set
                if (avec.dim(1) > 1) {
                    at_to_bf.push_back(ilo + 1);
                    at_nbf.push_back(avec.dim(1) - 1);
                }
            }
            else {            
                // For heavier atoms separate out 2s2p
                at_to_bf.push_back(ilo+1);
                at_nbf.push_back(4);
                if (avec.dim(1) > 5) {
                    at_to_bf.push_back(ilo + 5);
                    at_nbf.push_back(avec.dim(1) - 5);
                }
            }

            ilo = ihi;
        }
        MADNESS_ASSERT(ilo == nao);
        MADNESS_ASSERT(std::accumulate(at_nbf.begin(), at_nbf.end(), 0) == nao);
//        MADNESS_ASSERT(at_to_bf.back() + at_nbf.back() == nao);
        //print("newloc", at_to_bf, at_nbf);
    } else {
        aobasis.shells_to_bfn(molecule, at_to_bf, at_nbf);
        //aobasis.atoms_to_bfn(molecule, at_to_bf, at_nbf);
    }

    // Below here atoms may be shells or atoms --- by default shells

    int natom = at_to_bf.size();

    tensorT U(nmo, nmo);
    for (int i = 0; i < nmo; ++i) U(i, i) = 1.0;


    default_random_generator.setstate(
            182041 + world.rank() * 10101); // To help with reproducibility for debugging, etc.

    if (world.rank() == 0) {
        //MKL_Set_Num_Threads_Local(16);

        tensorT Q(nmo, natom);
        double breaksym = 0.0; // was 1e-3;
        auto QQ = [&at_to_bf, &at_nbf, &breaksym](const tensorT& C, int i, int j, int a) -> double {
            int lo = at_to_bf[a], nbf = at_nbf[a];
            const double *Ci = &C(i, lo);
            const double *Cj = &C(j, lo);
            double qij = 0.0;
            for (int mu = 0; mu < nbf; ++mu) qij += Ci[mu] * Cj[mu];
            return qij * (1.0 + breaksym * a); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! break symmetry
        };

        auto makeGW = [&Q, &nmo, &natom, &QQ](const tensorT& C, double& W, tensorT& g) -> void {
            W = 0.0;
            for (int i = 0; i < nmo; ++i) {
                for (int a = 0; a < natom; ++a) {
                    Q(i, a) = QQ(C, i, i, a);
                    W += Q(i, a) * Q(i, a);
                }
            }

            for (int i = 0; i < nmo; ++i) {
                for (int j = 0; j < i; ++j) {
                    double Qiiij = 0.0, Qijjj = 0.0;
                    for (int a = 0; a < natom; ++a) {
                        double Qija = QQ(C, i, j, a);
                        Qijjj += Qija * Q(j, a);
                        Qiiij += Qija * Q(i, a);
                    }
                    g(j, i) = Qiiij - Qijjj;
                    g(i, j) = -g(j, i);
                }
            }
        };

        tensorT xprev; // previous search direction
        tensorT gprev; // previous gradient
        bool rprev = true; // if true previous iteration restricted step or did incomplete search (so don't do conjugate)
        const int N = (nmo * (nmo - 1)) / 2; // number of independent variables
        for (int iter = 0; iter < 1200; ++iter) {
            tensorT g(nmo, nmo);
            double W;

            makeGW(C, W, g);

            if (randomize && iter == 0) {
                for (int i = 0; i < nmo; ++i) {
                    for (int j = 0; j < i; ++j) {
                        g(i, j) += 0.1 * (RandomValue<double>() - 0.5);
                        g(j, i) = -g(i, j);
                    }
                }
            }

            double maxg = g.absmax();
            if (doprint) printf("iteration %d W=%.8f maxg=%.2e\n", iter, W, maxg);
            if (maxg < thresh) break;

            // construct search direction using conjugate gradient approach
            tensorT x = copy(g);
            if (!rprev) { // Only apply conjugacy if did LS with real gradient
                double gamma = g.trace(g - gprev) / gprev.trace(gprev);
                if (doprint) print("gamma", gamma);
                x.gaxpy(1.0, xprev, gamma);
            }

            // Perform the line search.
            rprev = false;
            double dxgrad = x.trace(g) * 2.0;  // 2*2 = 4 which should be prefactor on integrals in gradient
            if (dxgrad < 0 || ((iter + 1) % N) == 0) {
                if (doprint) print("resetting since dxgrad -ve or due to dimension", dxgrad, iter, N);
                x = copy(g);
                dxgrad = x.trace(g) * 2.0;
            }
            xprev = x; // Save for next iteration
            gprev = copy(g);

            double mu = 0.05 / std::max(0.1,
                                        maxg); // Take larger inital step (was 0.01), and restrict intial step mu by size of max gradient
            tensorT dU = matrix_exponential(x * mu);
            tensorT newC = inner(dU, C, 0, 0);
            double newW;
            makeGW(newC, newW, g);
            double dxgnew = x.trace(g) * 2.0;

            if (randomize && iter == 0) {
                rprev = true; // since did not use real gradient
            } else { // perform quadratic fit using f(0), df(0)/dx=dxgrad, f(mu) --- actually now use f(0), df(0)/dx, df(mu)/dx for better accuracy
                double f0 = W;
                double f1 = newW;
                //double hess = 2.0*(f1-f0-mu*dxgrad)/(mu*mu);
                double hess = (dxgnew - dxgrad) / mu; // Near convergence this is more accurate
                if (hess >= 0) {
                    if (doprint) print("+ve hessian", hess);
                    hess = -2.0 * dxgrad; // force a bigish step to get out of bad region
                    rprev = true; // since did not do line search
                }
                double mu2 = -dxgrad / hess;
                if (mu2 * maxg > 0.25) {
                    mu2 = 0.25 / maxg; // pi/6 = 0.524, pi/4=0.785
                    rprev = true; // since did not do line search
                }
                double f2p = f0 + dxgrad * mu2 + 0.5 * hess * mu2 * mu2;
                if (doprint) print(f0, f1, f0 - f1, f2p, "dxg", dxgrad, "hess", hess, "mu", mu, "mu2", mu2);
                mu = mu2;
            }

            // if (maxg < 10 * thresh) {
            //     breaksym = 0.0; // was 1e-5;
            //     rprev = true; // since just messed up the gradient
            // }

            dU = matrix_exponential(x * mu);
            U = inner(U, dU, 1, 0);
            C = inner(dU, C, 0, 0);
        }
        bool switched = true;
        while (switched) {
            switched = false;
            for (int i = 0; i < nmo; i++) {
                for (int j = i + 1; j < nmo; j++) {
                    if (set[i] == set[j]) {
                        double sold = U(i, i) * U(i, i) + U(j, j) * U(j, j);
                        double snew = U(i, j) * U(i, j) + U(j, i) * U(j, i);
                        if (snew > sold) {
                            tensorT tmp = copy(U(_, i));
                            U(_, i) = U(_, j);
                            U(_, j) = tmp;
                            switched = true;
                        }
                    }
                }
            }
        }

        // Fix phases.
        for (int i = 0; i < nmo; ++i) {
            if (U(i, i) < 0.0)
                U(_, i).scale(-1.0);
        }
        //MKL_Set_Num_Threads_Local(1);
    }
    //done:
    world.gop.broadcast(U.ptr(), U.size(), 0);

    DistributedMatrix<double> dUT = column_distributed_matrix<double>(world, nmo, nmo);
    dUT.copy_from_replicated(transpose(U));

    // distmatT dUT = distributed_localize_PM(world, mo, ao, set, at_to_bf, at_nbf,
    //                                        thresh, thetamax, randomize, doprint);
    //print(UT);
    return dUT;
}

template<typename T>
std::size_t Localizer::determine_frozen_orbitals(const Tensor<T> fmat) {

    // freezing based on canonical orbital energies
    Tensor<double> eps;
    if (fmat.ndim()==2) {
        auto [eval, evec] = syev(fmat);
        eps=real(eval);
    } else {
        MADNESS_CHECK(fmat.ndim()==1);
        eps=real(fmat);
    }

    std::vector<Function<T,3>> vec(eps.size());
    MolecularOrbitals<T, 3> dummy_mo(vec, eps);
    dummy_mo.recompute_localize_sets();
    auto s = MolecularOrbitals<double, 3>::convert_set_to_slice(dummy_mo.get_localize_sets());

    long nactive = s.back().end - s.back().start + 1;
    long nfrozen = dummy_mo.get_mos().size() - nactive;

    if (not Localizer::check_frozen_consistency(nfrozen,dummy_mo.get_localize_sets())) {
        dummy_mo.pretty_print("mos");
        print("nfrozen", nfrozen);
        MADNESS_EXCEPTION("inconsistent number of frozen orbitals",1);
    };

    return nfrozen;
}

bool Localizer::check_frozen_consistency(const long nfrozen, const std::vector<int>& localize_sets) {

    // fast return
    if (nfrozen==0) return true;

    // check that freeze is consistent with the core/valence block-diagonal structure of the fock matrix
    bool fine = false;
    int ntotal = 0;
    std::vector<Slice> slices = MolecularOrbitals<double, 3>::convert_set_to_slice(localize_sets);
    for (std::size_t i = 0; i < slices.size(); ++i) {
        const Slice& s = slices[i];
        int n_in_set = s.end - s.start + 1;
        ntotal += n_in_set;
        if (nfrozen == ntotal) {
            fine = true;
            break;
        }
    }
    return fine;
}

/// given a unitary transformation matrix undo mere reordering
template<typename T>
void Localizer::undo_reordering(Tensor<T>& U, const Tensor<double>& occ, Tensor<double>& evals) {
    long nmo = U.dim(0);

    // Within blocks with the same occupation number attempt to
    // keep orbitals in the same order (to avoid confusing the
    // non-linear solver).
    // !!!!!!!!!!!!!!!!! NEED TO RESTRICT TO OCCUPIED STATES?
    bool switched = true;
    while (switched) {
        switched = false;
        for (int i = 0; i < nmo; i++) {
            for (int j = i + 1; j < nmo; j++) {
                if (occ(i) == occ(j)) {
//                    double sold = std::abs(U(i, i) * U(i, i) + U(j, j) * U(j, j));
//                    double snew = std::abs(U(i, j) * U(i, j) + U(j, i) * U(j, i));
                    double sold = std::real(U(i,i)*std::conj(U(i,i))) + std::real(U(j,j)*std::conj(U(j,j)));
                    double snew = std::real(U(i,j)*std::conj(U(i,j))) + std::real(U(j,i)*std::conj(U(j,i)));
                    if (snew > sold) {
                        Tensor<T> tmp = copy(U(_, i));
                        U(_, i) = U(_, j);
                        U(_, j) = tmp;
                        std::swap(evals[i], evals[j]);
                        switched = true;
                    }
                }
            }
        }
    }

    // Fix phases.
//    for (long i = 0; i < nmo; ++i)
//        if (std::real(U(i, i)) < 0.0)
//            U(_, i).scale(-1.0);

    for (long i = 0; i < nmo; ++i) U(_, i).scale(conditional_conj(U(i,i))/std::abs(U(i,i)));

}


std::vector<Slice> Localizer::find_degenerate_blocks(const Tensor<double>& eval, const double thresh_degenerate) {
    long nmo = eval.size();
    std::vector<Slice> blocks;
    // Rotations between effectively degenerate states confound
    // the non-linear equation solver ... undo these rotations
    long ilo = 0; // first element of cluster
    while (ilo < nmo - 1) {
        long ihi = ilo;
        while (fabs(eval[ilo] - eval[ihi + 1])
               < thresh_degenerate * 10.0 * std::max(fabs(eval[ilo]), 1.0)) {
            ++ihi;
            if (ihi == nmo - 1)
                break;
        }
        blocks.push_back(Slice(ilo, ihi));
        ilo = ihi + 1;
    }
    return blocks;
}

template<typename T>
Tensor<T> Localizer::undo_rotation(const Tensor<T>& U_in, const std::vector<Slice>& blocks) {
    Tensor<T> U = copy(U_in);
    for (auto block: blocks) {
        long ncluster = block.end - block.start + 1;
        if (ncluster > 1) {
            Tensor<T> q = copy(U(block, block));

            // Polar Decomposition
            Tensor<T> VH(ncluster, ncluster);
            Tensor<T> W(ncluster, ncluster);
            Tensor<double> sigma(ncluster);

            svd(q, W, sigma, VH);
            q = transpose(inner(W, VH)).conj();
            U(_, block) = inner(U(_, block), q);
        }
    }
    return U;
}

template<typename T>
void Localizer::undo_rotations_within_sets(Tensor<T>& U, const std::vector<int>& localized_set) {
    std::vector<Slice> blocks = MolecularOrbitals<T, 3>::convert_set_to_slice(localized_set);
//    for (auto block : blocks) print("block",block);
    U = undo_rotation(U, blocks);
}

/// given a unitary transformation matrix undo rotations between degenerate columns
template<typename T>
void Localizer::undo_degenerate_rotations(Tensor<T>& U, const Tensor<double>& evals, const double thresh_degenerate) {
    MADNESS_CHECK(evals.size() == U.dim(0));
    std::vector<Slice> degenerate_blocks = find_degenerate_blocks(evals, thresh_degenerate);
    U = undo_rotation(U, degenerate_blocks);
}


template<typename T>
Tensor<T> Localizer::matrix_exponential(const Tensor<T>& A) const {
    MADNESS_CHECK(A.dim(0) == A.dim(1));
    typedef Tensor<T> tensorT;

    // Power iteration to estimate the 2-norm of the matrix. Used
    // to use Frobenius or 1-norms but neither were very tight.
    double anorm;
    {
        tensorT x(A.dim(0));
        x.fillrandom();
        x.scale(1.0 / x.normf());
        double prev = 0.0;
        for (int i = 0; i < 100; i++) {
            tensorT xnew = inner(A, inner(A, x, 1, 0), 0, 0);
            anorm = std::sqrt(std::abs((x.trace(xnew))));
            double err = std::abs(prev - anorm) / anorm;
            //print(i,anorm,err,A.normf());
            if (err < 0.01) break; // just need 1-2 digits
            x = xnew.scale(1.0 / xnew.normf());
            prev = anorm;
        }
    }

    // Scale A by a power of 2 until it is "small"
    int n = 0;
    double scale = 1.0;
    while (anorm * scale > 0.089) { // so that 9th order expansion is accurate to 1e-15
        ++n;
        scale *= 0.5;
    }
    tensorT B = scale * A;    // B = A*2^-n

    // Make identity
    tensorT I = tensorT(2, B.dims());
    for (int i = 0; i < I.dim(0); ++i) I(i, i) = 1.0;

    // Compute exp(B) using Taylor series optimized to reduce cost --- Chebyshev is only a minor improvement
    tensorT expB;
    if (anorm > 0.24e-1) {
        tensorT B2 = inner(B, B);
        tensorT B4 = inner(B2, B2);
        tensorT B6 = inner(B4, B2);
        expB = I + inner(B, B6 + 42. * B4 + 840. * B2 + 5040. * I).scale(1. / 5040.) +
               inner(B2, B6 + 56. * B4 + 1680. * B2 + 20160. * I).scale(1. / 40320.);
    } else if (anorm > 0.26e-2) {
        tensorT B2 = inner(B, B);
        tensorT B4 = inner(B2, B2);
        expB = I + inner(B, 42. * B4 + 840. * B2 + 5040. * I).scale(1. / 5040.) +
               inner(B2, 56. * B4 + 1680. * B2 + 20160. * I).scale(1. / 40320.);
    } else if (anorm > 0.18e-4) {
        tensorT B2 = inner(B, B);
        expB = I + inner(B, 840. * B2 + 5040. * I).scale(1. / 5040.) +
               inner(B2, 1680. * B2 + 20160. * I).scale(1. / 40320.);
    } else if (anorm > 4.5e-8) {
        expB = I + B + inner(B, B).scale(0.5);
    } else {
        expB = I + B;
    }

    // // Old algorithm
    // tensorT oldexpB = copy(I);
    // const double tol = 1e-13;
    // int k = 1;
    // tensorT term = B;
    // while (term.normf() > tol) {
    //     oldexpB += term;
    //     term = inner(term, B);
    //     ++k;
    //     term.scale(1.0 / k);
    // }
    // Error check for validation
    // double err = (expB-oldexpB).normf();
    // print("matxerr", anorm, err);

    // Repeatedly square to recover exp(A)
    while (n--) expB = inner(expB, expB);

    return expB;
}


template
MolecularOrbitals<double, 3> Localizer::localize(const MolecularOrbitals<double, 3>& mo_in, bool randomize) const;

template
MolecularOrbitals<double, 3> Localizer::localize(const MolecularOrbitals<double, 3>& mo_in, const Tensor<double>& Fock,
                                                 bool randomize) const;

template
std::size_t Localizer::determine_frozen_orbitals(const Tensor<double> fmat);

template
std::size_t Localizer::determine_frozen_orbitals(const Tensor<double_complex> fmat);

template
MolecularOrbitals<double, 3> Localizer::separate_core_valence(const MolecularOrbitals<double, 3>& mo_in,
                                                            const Tensor<double>& Fock) const;

//template
//MolecularOrbitals<double_complex, 3> Localizer::separate_core_valence(const MolecularOrbitals<double_complex, 3>& mo_in,
//                                                              const Tensor<double_complex>& Fock) const;

template
void
Localizer::undo_degenerate_rotations(Tensor<double>& U, const Tensor<double>& evals, const double thresh_degenerate);

template
void
Localizer::undo_degenerate_rotations(Tensor<double_complex>& U, const Tensor<double>& evals, const double thresh_degenerate);

template
void
Localizer::undo_reordering(Tensor<double>& U, const Tensor<double>& occ, Tensor<double>& eval);

template
void
Localizer::undo_reordering(Tensor<double_complex>& U, const Tensor<double>& occ, Tensor<double>& eval);
}
