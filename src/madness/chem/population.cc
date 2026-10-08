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

/// \file population.cc
/// \brief atomic populations of occupied orbitals: Mulliken, Loewdin and IAO (see population.h)

#include <madness/chem/population.h>
#include <madness/tensor/tensor_lapack.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <utility>

namespace madness {
namespace population {

namespace {

/// the inputs fit together: S square, X with S's rows, one occupation per orbital, one atom per basis function
void check_inputs(const Tensor<double>& S, const Tensor<double>& X, const Tensor<double>& occ,
                  const std::vector<int>& atom_of_bf, const long natom) {
    MADNESS_CHECK_THROW(S.ndim() == 2 and S.dim(0) == S.dim(1), "population: S must be a square matrix");
    MADNESS_CHECK_THROW(X.ndim() == 2 and X.dim(0) == S.dim(0), "population: X must have a row per basis function");
    MADNESS_CHECK_THROW(occ.ndim() == 1 and occ.dim(0) == X.dim(1), "population: one occupation per orbital");
    MADNESS_CHECK_THROW(long(atom_of_bf.size()) == S.dim(0), "population: one atom per basis function");
    for (const int a : atom_of_bf)
        MADNESS_CHECK_THROW(a >= 0 and a < natom, "population: atom index of a basis function out of range");
}

/// electrons per atom from the per-basis-function populations
Tensor<double> sum_on_atoms(const Tensor<double>& per_bf, const std::vector<int>& atom_of_bf, const long natom) {
    Tensor<double> e(natom);
    for (long mu = 0; mu < per_bf.dim(0); ++mu) e(atom_of_bf[mu]) += per_bf(mu);
    return e;
}

} // namespace


const std::vector<std::string>& schemes() {
    static const std::vector<std::string> s = {"mulliken", "lowdin", "iao"};
    return s;
}


Tensor<double> matrix_power(const Tensor<double>& S, const double p, const double tol) {
    Tensor<double> U, e;
    syev(S, U, e);
    const double emax = e.size() ? e.max() : 0.0;
    MADNESS_CHECK_THROW(emax > 0.0, "population: matrix_power needs a positive semidefinite matrix");
    Tensor<double> Ue = copy(U);
    for (long k = 0; k < e.dim(0); ++k) {
        const double f = (e(k) > tol * emax) ? std::pow(e(k), p) : 0.0;
        Ue(_, k).scale(f);
    }
    return inner(Ue, transpose(U));
}


SpinPopulation projection_population(const std::string& scheme, const Tensor<double>& S, const Tensor<double>& X,
                                     const Tensor<double>& occ, const std::vector<int>& atom_of_bf,
                                     const long natom) {
    check_inputs(S, X, occ, atom_of_bf, natom);
    MADNESS_CHECK_THROW(scheme == "mulliken" or scheme == "lowdin",
                        "population: projection schemes are mulliken and lowdin");
    const long nbf = S.dim(0), nocc = X.dim(1);
    SpinPopulation r;
    r.electrons = Tensor<double>(natom);
    r.completeness = Tensor<double>(nocc);
    if (nocc == 0) return r;
    // the fits C = S^-1 X, and W = X^T S^-1 X = C^T S C, the overlap of the fitted orbitals;
    // its diagonal is the completeness c_i of each orbital
    const Tensor<double> C = inner(matrix_power(S, -1.0), X);
    const Tensor<double> W = inner(transpose(X), C);
    for (long i = 0; i < nocc; ++i) r.completeness(i) = W(i, i);
    {
        Tensor<double> U, w;
        syev(W, U, w);
        r.min_w_eigenvalue = w(0L);
    }
    MADNESS_CHECK_THROW(r.min_w_eigenvalue > 1.e-8,
                        "population: the basis misses an occupied orbital (W = X^T S^-1 X is singular)");
    // The fitted orbitals orthonormalized among themselves, Cbar = C W^-1/2: the occupied space
    // projected onto the basis. Renormalizing each fit by its own c_i instead would also give
    // the sum rule, but the c_i depend on the orbital gauge, so localized and canonical
    // orbitals would get different populations.
    const Tensor<double> Cbar = inner(C, matrix_power(W, -0.5));
    Tensor<double> per_bf(nbf);
    if (scheme == "lowdin") {
        const Tensor<double> Y = inner(matrix_power(S, 0.5), Cbar);
        for (long i = 0; i < nocc; ++i)
            for (long mu = 0; mu < nbf; ++mu) per_bf(mu) += occ(i) * Y(mu, i) * Y(mu, i);
    } else {
        const Tensor<double> SCbar = inner(S, Cbar);
        for (long i = 0; i < nocc; ++i)
            for (long mu = 0; mu < nbf; ++mu) per_bf(mu) += occ(i) * Cbar(mu, i) * SCbar(mu, i);
    }
    r.electrons = sum_on_atoms(per_bf, atom_of_bf, natom);
    double ncap = 0.0, nsum = 0.0;
    for (long i = 0; i < nocc; ++i) {
        ncap += occ(i) * r.completeness(i);
        nsum += occ(i);
    }
    r.spilling = (nsum > 0.0) ? 1.0 - ncap / nsum : 0.0;
    return r;
}


SpinPopulation iao_population(const Tensor<double>& S, const Tensor<double>& X, const Tensor<double>& occ,
                              const std::vector<int>& atom_of_bf, const long natom) {
    check_inputs(S, X, occ, atom_of_bf, natom);
    const long nbf = S.dim(0), nocc = X.dim(1);
    SpinPopulation r;
    r.electrons = Tensor<double>(natom);
    if (nocc == 0) return r;
    MADNESS_CHECK_THROW(nocc <= nbf, "population: iao needs at least as many minimal basis functions as orbitals");

    const Tensor<double> SinvX = inner(matrix_power(S, -1.0), X);
    const Tensor<double> W = inner(transpose(X), SinvX);
    {
        Tensor<double> U, w;
        syev(W, U, w);
        r.min_w_eigenvalue = w(0L);
    }
    MADNESS_CHECK_THROW(r.min_w_eigenvalue > 1.e-8,
                        "population: iao: the minimal basis misses an occupied orbital (W is singular)");
    const Tensor<double> D = inner(SinvX, matrix_power(W, -0.5));
    Tensor<double> DDtS = inner(inner(D, transpose(D)), S);
    Tensor<double> M = -1.0 * DDtS;
    Tensor<double> twoDDtS_1 = 2.0 * DDtS;
    for (long mu = 0; mu < nbf; ++mu) {
        M(mu, mu) += 1.0;
        twoDDtS_1(mu, mu) -= 1.0;
    }
    const Tensor<double> K = inner(transpose(X), twoDDtS_1);                 // nocc x nbf
    const Tensor<double> MtX = inner(transpose(M), X);                        // nbf x nocc
    const Tensor<double> MtXK = inner(MtX, K);                                // nbf x nbf
    const Tensor<double> SA = inner(inner(transpose(M), S), M) + MtXK + transpose(MtXK)
                              + inner(transpose(K), K);
    const Tensor<double> AP = MtX + transpose(K);                             // <A|phi>
    const Tensor<double> Q = inner(matrix_power(SA, -0.5), AP);

    Tensor<double> per_bf(nbf);
    for (long i = 0; i < nocc; ++i) {
        double norm = 0.0;
        for (long rho = 0; rho < nbf; ++rho) {
            per_bf(rho) += occ(i) * Q(rho, i) * Q(rho, i);
            norm += Q(rho, i) * Q(rho, i);
        }
        r.max_norm_error = std::max(r.max_norm_error, std::abs(norm - 1.0));
    }
    r.electrons = sum_on_atoms(per_bf, atom_of_bf, natom);
    return r;
}


SpinPopulation spin_population(const std::string& scheme, const Tensor<double>& S, const Tensor<double>& X,
                               const Tensor<double>& occ, const std::vector<int>& atom_of_bf, const long natom) {
    if (scheme == "iao") return iao_population(S, X, occ, atom_of_bf, natom);
    return projection_population(scheme, S, X, occ, atom_of_bf, natom);
}


AtomicPopulations atomic_populations(const std::string& scheme, const std::string& basis, const long nbf,
                                     const std::vector<std::string>& symbols, const std::vector<double>& Z,
                                     const long nalpha, const long nbeta, const SpinPopulation& alpha,
                                     const SpinPopulation& beta) {
    const long natom = long(Z.size());
    MADNESS_CHECK_THROW(long(symbols.size()) == natom and alpha.electrons.dim(0) == natom,
                        "population: one symbol, charge and population per atom");
    MADNESS_CHECK_THROW(nbeta == 0 or beta.electrons.dim(0) == natom, "population: beta populations missing");
    AtomicPopulations r;
    r.scheme = scheme;
    r.basis = basis;
    r.nbf = nbf;
    r.nalpha = nalpha;
    r.nbeta = nbeta;
    r.symbols = symbols;
    r.alpha = alpha;
    r.beta = beta;
    for (long a = 0; a < natom; ++a) {
        const double nb = (nbeta > 0) ? beta.electrons(a) : 0.0;
        r.charge.push_back(Z[a] - alpha.electrons(a) - nb);
        r.spin.push_back(alpha.electrons(a) - nb);
    }
    return r;
}


void AtomicPopulations::print(const double total_charge, const int print_level) const {
    const bool o = open();
    printf("\n population analysis: %s (%s basis %s, %ld functions)\n\n", scheme.c_str(),
           scheme == "iao" ? "minimal" : "projection", basis.c_str(), nbf);
    printf("     atom      charge%s\n", o ? "        spin" : "");
    double qsum = 0.0, ssum = 0.0;
    for (std::size_t a = 0; a < charge.size(); ++a) {
        qsum += charge[a];
        ssum += spin[a];
        if (o) printf("   %3zu %-3s %11.6f %11.6f\n", a, symbols[a].c_str(), charge[a], spin[a]);
        else printf("   %3zu %-3s %11.6f\n", a, symbols[a].c_str(), charge[a]);
    }
    if (o) printf("   sum     %11.6f %11.6f   (expected %.1f and %ld)\n", qsum, ssum, total_charge, nalpha - nbeta);
    else printf("   sum     %11.6f   (expected %.1f)\n", qsum, total_charge);
    const bool b = (nbeta > 0);
    if (scheme == "iao") {
        printf("   min eigenvalue of W %.3e, max |Q column norm - 1| %.1e\n",
               b ? std::min(alpha.min_w_eigenvalue, beta.min_w_eigenvalue) : alpha.min_w_eigenvalue,
               b ? std::max(alpha.max_norm_error, beta.max_norm_error) : alpha.max_norm_error);
    } else if (o) {
        printf("   spilling: alpha %.3e, beta %.3e\n", alpha.spilling, b ? beta.spilling : 0.0);
    } else {
        printf("   spilling: %.3e\n", alpha.spilling);
    }
    if (scheme != "iao" and print_level > 3) {
        // the orbitals the basis describes worst
        for (const auto& [label, p] : {std::make_pair("alpha", &alpha), std::make_pair("beta", &beta)}) {
            if (p->completeness.size() == 0 or (p == &beta and not o)) continue;
            long worst = 0;
            for (long i = 1; i < p->completeness.dim(0); ++i)
                if (p->completeness(i) < p->completeness(worst)) worst = i;
            printf("   %s orbital %ld is the least complete: %.6f\n", label, worst, p->completeness(worst));
        }
    }
}


nlohmann::json AtomicPopulations::to_json() const {
    nlohmann::json j;
    j["basis"] = basis;
    j["charges"] = charge;
    std::vector<double> ea(charge.size()), eb(charge.size());
    for (std::size_t a = 0; a < charge.size(); ++a) {
        ea[a] = alpha.electrons(long(a));
        eb[a] = (nbeta > 0) ? beta.electrons(long(a)) : 0.0;
    }
    j["electrons_alpha"] = ea;
    j["electrons_beta"] = eb;
    if (open()) j["spin"] = spin;
    j["min_w_eigenvalue"] = (nbeta > 0) ? std::min(alpha.min_w_eigenvalue, beta.min_w_eigenvalue)
                                        : alpha.min_w_eigenvalue;
    if (scheme == "iao") {
        j["max_norm_error"] = (nbeta > 0) ? std::max(alpha.max_norm_error, beta.max_norm_error)
                                          : alpha.max_norm_error;
    } else {
        j["spilling_alpha"] = alpha.spilling;
        j["spilling_beta"] = (nbeta > 0) ? beta.spilling : 0.0;
        // each orbital's completeness c_i, in the order of the orbitals given (for canonical
        // orbitals the last is the HOMO, the usual worst case for diffuse anions)
        const auto list = [](const Tensor<double>& c) {
            std::vector<double> v;
            for (long i = 0; i < c.size(); ++i) v.push_back(c(i));
            return v;
        };
        j["completeness_alpha"] = list(alpha.completeness);
        j["completeness_beta"] = (nbeta > 0) ? list(beta.completeness) : std::vector<double>();
    }
    return j;
}

} // namespace population
} // namespace madness
