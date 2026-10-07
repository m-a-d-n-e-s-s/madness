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

/// \file lcao_cholesky.h
/// \brief pivoted Cholesky decomposition of the two-electron integrals over Gaussian basis functions

#ifndef MADNESS_CHEM_LCAO_CHOLESKY_H__INCLUDED
#define MADNESS_CHEM_LCAO_CHOLESKY_H__INCLUDED

#include <madness/chem/lcao_integrals.h>
#include <madness/chem/lcao_scf.h>

#include <memory>
#include <vector>

namespace madness {

class World;

namespace lcao {

/// pivoted, incomplete Cholesky decomposition of the two-electron integral matrix, V ~ L L^T
///
/// V_(mu nu),(lambda sigma) = (mu nu|lambda sigma) runs over the function pairs (FunctionPairs).
/// It is positive semidefinite, and the decomposition stops when the largest residual
/// diagonal is at most tol, so every element of V - L L^T is at most tol in magnitude
/// (Beebe and Linderberg 1977; Koch, Sanchez de Meras and Pedersen 2003):
/// - Pair screening: rows with (mu nu|mu nu) < tol^2/D_max are dropped. Their integrals are
///   bounded by sqrt(D D_max) < tol.
/// - Pivoting by batches: the group pair of the largest residual diagonal D gives the
///   candidates, its rows with D > max(tol, span D_max). Their integral columns are computed
///   for all kept rows (eri_columns), the previous vectors are subtracted (one dgemm per block
///   of rows), and the candidates become vectors in descending D, re-checked after each.
/// - All kernel terms are kept and nothing is Schwarz-skipped, so V stays positive
///   semidefinite. Residual diagonals that rounding makes negative are set to 0.
///
/// Tasks on the thread pool of each rank, with sequential BLAS inside each task. The blocks of
/// rows are consecutive group pairs, cut where they reach 4096 kept rows, independent of the
/// numbers of threads and ranks, and the pivot choice breaks ties by the lowest row.
///
/// With collective, every rank of world takes part (all must construct it, in the same order):
/// - The blocks of rows go to the ranks in contiguous ranges; each rank computes the integral
///   columns, updates and vectors of its rows, so no integral data moves.
/// - Each batch takes two reductions (the pivot) and two broadcasts from the owner of the
///   pivot's group pair: the candidates with their D and rows of L, then their block of V after
///   the update. Every rank then makes the same decisions from the same numbers, so the vectors
///   do not depend on the number of ranks (hash()).
/// - Afterwards the vectors are redistributed in chunks of consecutive vectors, each complete
///   over the kept rows: chunk c is held by rank c % nproc. CholeskyERI works on the chunks.
class CholeskyERIDecomposition {
public:
    /// what the decomposition did, and its wall times in seconds (of this rank)
    struct Stats {
        std::size_t npairs = 0;     ///< function pairs
        std::size_t nkept = 0;      ///< rows kept after pair screening
        std::size_t ncolumns = 0;   ///< integral columns computed (whole group pairs)
        std::size_t nbatches = 0;   ///< pivot batches, one eri_columns call each
        double t_diagonal = 0.0;    ///< pair_diagonal
        double t_integrals = 0.0;   ///< eri_columns, and gathering the candidates' columns
        double t_updates = 0.0;     ///< subtracting the previous vectors, making the new ones
        double t_redistribute = 0.0;    ///< collective: moving the vectors into chunks
    };

    /// @param[in] tol         largest residual diagonal left, which bounds every element of |V - L L^T|
    /// @param[in] span        candidates of a batch need D > span D_max (and > tol)
    /// @param[in] collective  every rank of world takes part and holds a share; otherwise this rank does it all
    CholeskyERIDecomposition(World& world, const SeparatedGaussianIntegrals& ints, double tol, double span = 1.e-2,
                             bool collective = false);

    double tol() const { return tol_; }
    double span() const { return span_; }

    /// true if the vectors are spread over the ranks of world (collective, and more than one rank)
    bool distributed() const { return distributed_; }

    /// the function pairs, all of them: the row numbering of V
    const FunctionPairs& pairs() const { return pairs_; }

    /// the rows kept after pair screening, in ascending order; the vectors are given on these rows
    const std::vector<std::size_t>& rows() const { return rows_; }
    std::size_t nkept() const { return rows_.size(); }

    /// the number of vectors, M
    long nvec() const { return nvec_; }

    /// the vectors come in nchunk() chunks of consecutive vectors, [first, last) of chunk c
    long nchunk() const { return nchunk_; }
    std::pair<long,long> chunk(const long c) const { return {nvec_ * c / nchunk_, nvec_ * (c + 1) / nchunk_}; }

    /// true if this rank holds chunk c
    bool holds(const long c) const { return not distributed_ or c % nproc_ == rank_; }

    /// the vectors of a held chunk one after the other, nkept() values each
    const double* chunk_vectors(const long c) const;

    /// all vectors one after the other, nkept() values each: L_k(rows()[i]) at [k nkept() + i];
    /// only if not distributed
    const double* vectors() const;

    /// an exact hash of all vectors: the XOR of their bit patterns, each rotated by its position;
    /// the same on any number of ranks for the same vectors (collective if distributed)
    unsigned long hash(World& world) const;

    const Stats& stats() const { return stats_; }

private:
    double tol_, span_;
    bool distributed_ = false;
    int nproc_ = 1, rank_ = 0;
    FunctionPairs pairs_;
    std::vector<std::size_t> rows_;
    std::vector<double> L_;                     ///< not distributed: all vectors
    std::vector<std::vector<double>> chunks_;   ///< distributed: the held chunks
    long nvec_ = 0;
    long nchunk_ = 1;
    Stats stats_;
};

/// J and K from the Cholesky vectors of the two-electron integrals, (mu nu|l s) ~ sum_k L_k(mu nu) L_k(l s)
///
/// - J_{mu nu} = sum_k L_k(mu nu) gamma_k, with gamma_k = sum_{l s} L_k(l s) P_{l s} over all
///   index orders.
/// - K = sum_k L_k P L_k, the vectors as symmetric N x N matrices. It is computed as
///   sum_k X_k X_k^T, X_k = L_k Y, from P = Y+ Y+^T - Y- Y-^T, the eigenpairs of P (the
///   negative part is empty for a density matrix, which is positive semidefinite).
/// - One task per held chunk of vectors (several vectors per dgemm for K). The partial sums of
///   all chunks are gathered and added in chunk order on every rank, so J and K do not depend on
///   the scheduling or on the number of ranks. Collective if the decomposition is distributed.
class CholeskyERI : public TwoElectronBuilder {
public:
    /// @param[in] nbf  the number of basis functions, which the density matrices must match
    CholeskyERI(World& world, std::shared_ptr<const CholeskyERIDecomposition> chol, long nbf);

    void jk(const Tensor<double>& Pa, const Tensor<double>& Pb, Tensor<double>& J, Tensor<double>& Ka,
            Tensor<double>& Kb) const override;

    const CholeskyERIDecomposition& decomposition() const { return *chol_; }

private:
    World& world_;
    std::shared_ptr<const CholeskyERIDecomposition> chol_;
    long nbf_;
    std::vector<int> mu_, nu_;      ///< per kept row, its two basis functions
};

} // namespace lcao
} // namespace madness

#endif // MADNESS_CHEM_LCAO_CHOLESKY_H__INCLUDED
