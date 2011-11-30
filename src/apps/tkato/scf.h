/*
   This file is part of MADNESS.

   Copyright (C) 2007,2010 Oak Ridge National Laboratory

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
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


$Id$
*/

/// \file scf.h
/// \brief Molecular HF and DFT code
/// \defgroup moldft The molecular density funcitonal and Hartree-Fock code

#ifndef MADNESS_SCF_H
#define MADNESS_SCF_H

#include <mra/mra.h>
#include <moldft/xcfunctional.h>

struct SCF {
    typedef madness::Function<double,3> functionT;
    typedef std::vector<functionT> vecfuncT;
    typedef std::pair<vecfuncT,vecfuncT> pairvecfuncT;
    typedef std::vector<pairvecfuncT> subspaceT;
    typedef madness::Tensor<double> tensorT;
    typedef madness::SeparatedConvolution<double,3> operatorT;
    typedef std::shared_ptr<operatorT> poperatorT;

    Molecule molecule;
    SCFParameters param;
    XCfunctional xc;
    AtomicBasisSet aobasis;
    functionT vnuc;
    functionT mask;
    vecfuncT amo, bmo;
    std::vector<int> aset, bset;
    vecfuncT ao;
    std::vector<int> at_to_bf, at_nbf;
    tensorT aocc, bocc;
    tensorT aeps, beps;
    poperatorT coulop;
    std::vector< std::shared_ptr<real_derivative_3d> > gradop;
    double vtol;
    double current_energy;

    SCF(World & world, const char *filename);

    void set_protocol(World & world, double thresh);

    void save_mos(World& world);

    void load_mos(World& world);

    void do_plots(World& world);

    void project(World & world);

    void make_nuclear_potential(World & world);

    void project_ao_basis(World & world);

    double PM_q(const tensorT & S, const tensorT & C, int i, int j, int lo, int nbf);

    void localize_PM_task_kernel(tensorT & Q, std::vector<tensorT> & Svec, tensorT & C,
            const bool & doprint, const std::vector<int> & set,
            const double thetamax, tensorT & U, const double thresh);

    tensorT localize_PM(World & world, const vecfuncT & mo, const std::vector<int> & set, const double thresh = 1e-9, const double thetamax = 0.5, const bool randomize = true, const bool doprint = true);

    void analyze_vectors(World & world, const vecfuncT & mo, const tensorT & occ = tensorT(), const tensorT & energy = tensorT(), const std::vector<int> & set = std::vector<int>());

    inline double DIP(const tensorT & dip, int i, int j, int k, int l);

    tensorT localize_boys(World & world, const vecfuncT & mo, const std::vector<int> & set, const double thresh = 1e-9, const double thetamax = 0.5, const bool randomize = true);

    tensorT kinetic_energy_matrix(World & world, const vecfuncT & v);

    vecfuncT core_projection(World & world, const vecfuncT & psi, const bool include_Bc = true);

    double core_projector_derivative(World & world, const vecfuncT & mo, const tensorT & occ, int atom, int axis);

    void initial_guess(World & world);

    void initial_load_bal(World & world);

    functionT make_density(World & world, const tensorT & occ, const vecfuncT & v);

    std::vector<poperatorT> make_bsh_operators(World & world, const tensorT & evals);

    vecfuncT apply_hf_exchange(World & world, const tensorT & occ, const vecfuncT & psi, const vecfuncT & f);

    // Used only for initial guess that is always spin-restricted LDA
    functionT make_lda_potential(World & world, const functionT & arho);


    functionT make_dft_potential(World & world, const vecfuncT& vf, int what);

    double make_dft_energy(World & world, const vecfuncT& vf);

    vecfuncT apply_potential(World & world, const tensorT & occ, const vecfuncT & amo,
            const vecfuncT& vf, const vecfuncT& delrho, const functionT & vlocal, double & exc, int ispin);

    tensorT derivatives(World & world);

    tensorT dipole(World & world);

    void vector_stats(const std::vector<double> & v, double & rms, double & maxabsval);

    vecfuncT compute_residual(World & world, tensorT & occ, tensorT & fock, const vecfuncT & psi, vecfuncT & Vpsi, double & err);

    tensorT make_fock_matrix(World & world, const vecfuncT & psi, const vecfuncT & Vpsi, const tensorT & occ, double & ekinetic);

    /// Compute the two-electron integrals over the provided set of orbitals

    /// Returned is a *replicated* tensor of \f$(ij|kl)\f$ with \f$i>=j\f$
    /// and \f$k>=l\f$.  The symmetry \f$(ij|kl)=(kl|ij)\f$ is enforced.
    Tensor<double> twoint(World& world, const vecfuncT& psi);

    tensorT matrix_exponential(const tensorT& A);

    tensorT diag_fock_matrix(World & world, tensorT& fock, vecfuncT & psi, vecfuncT & Vpsi, tensorT & evals, const tensorT & occ, double thresh);

    void loadbal(World & world, functionT & arho, functionT & brho, functionT & arho_old, functionT & brho_old, subspaceT & subspace);

    void rotate_subspace(World& world, const tensorT& U, subspaceT& subspace, int lo, int nfunc, double trantol);

    void update_subspace(World & world,
            vecfuncT & Vpsia, vecfuncT & Vpsib,
            tensorT & focka, tensorT & fockb,
            subspaceT & subspace, tensorT & Q,
            double & bsh_residual, double & update_residual);

    // For given protocol, solve the DFT/HF/response equations
    void solve(World & world);
};


#endif

