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
*/

/// \file SCFOperators.cc
/// \brief Operators for the molecular HF and DFT code
/// \defgroup chem The molecular density functional and Hartree-Fock code


#include<madness/chem/SCFOperators.h>
#include<madness/chem/SCF.h>
#include<madness/chem/nemo.h>
#include<madness/chem/oep.h>
#include<madness/chem/correlationfactor.h>
#include<madness/chem/xcfunctional.h>
#include<madness/chem/exchangeoperator.h>


using namespace madness;

namespace madness {

template<typename T, std::size_t NDIM>
DistributedMatrix<T> Kinetic<T, NDIM>::kinetic_energy_matrix(World &world,
                                                             const vecfuncT &v) const {
    int n = v.size();
    DistributedMatrix<T> r = column_distributed_matrix<T>(world, n, n);
    reconstruct(world, v);

    // apply the derivative operator on each function for each dimension
    std::vector<vecfuncT> dv(NDIM);
    for (std::size_t i = 0; i < NDIM; ++i) {
        dv[i] = apply(world, *(gradop[i]), v, false);
    }
    world.gop.fence();
    for (std::size_t i = 0; i < NDIM; ++i) {
        compress(world, dv[i], false);
    }
    world.gop.fence();
    for (std::size_t i = 0; i < NDIM; ++i) {
        r += matrix_inner(r.distribution(), dv[i], dv[i], true);
    }
    r *= 0.5;
    return r;
}


template<typename T, std::size_t NDIM>
DistributedMatrix<T> Kinetic<T, NDIM>::kinetic_energy_matrix(World &world,
                                                             const vecfuncT &vbra, const vecfuncT &vket) const {
    int n = vbra.size();
    int m = vket.size();
    DistributedMatrix<T> r = column_distributed_matrix<T>(world, n, m);
    reconstruct(world, vbra);
    reconstruct(world, vket);
    const auto bra_equiv_ket = &vbra == &vket;

    // apply the derivative operator on each function for each dimension
    std::vector<vecfuncT> dvbra(NDIM), dvket(NDIM);
    for (std::size_t i = 0; i < NDIM; ++i) {
        dvbra[i] = apply(world, *(gradop[i]), vbra, false);
        dvket[i] = apply(world, *(gradop[i]), vket, false);
    }
    world.gop.fence();
    for (std::size_t i = 0; i < NDIM; ++i) {
        compress(world, dvbra[i], false);
        compress(world, dvket[i], false);
    }
    world.gop.fence();
    for (std::size_t i = 0; i < NDIM; ++i) {
        r += matrix_inner(r.distribution(), dvbra[i], dvket[i], bra_equiv_ket);
    }
    r *= 0.5;
    return r;
}

// explicit instantiation
template
class Kinetic<double, 1>;

template
class Kinetic<double, 2>;

template
class Kinetic<double, 3>;

template
class Kinetic<double, 4>;

template
class Kinetic<double, 5>;

template
class Kinetic<double, 6>;

template
class Kinetic<double_complex, 1>;

template
class Kinetic<double_complex, 2>;

template
class Kinetic<double_complex, 3>;

template
class Kinetic<double_complex, 4>;

template
class Kinetic<double_complex, 5>;

template
class Kinetic<double_complex, 6>;


template<typename T, std::size_t NDIM>
std::vector<Function<T, NDIM> >
Laplacian<T, NDIM>::operator()(const std::vector<Function<T, NDIM> > &vket) const {

    refine(world, vket);     // for better accuracy
    vecfuncT result = zero_functions_compressed<T, NDIM>(world, vket.size());
    SeparatedConvolution<T, NDIM> smooth = SmoothingOperator<NDIM>(world, eps);


    for (size_t idim = 0; idim < NDIM; ++idim) {
        vecfuncT dvket = apply(world, *gradop[idim].get(), vket);
        refine(world, dvket);
        if (eps > 0.0) dvket = apply(world, smooth, dvket);
        vecfuncT ddvket = apply(world, *gradop[idim].get(), dvket);
        result = add(world, result, ddvket);
    }

    if (eps > 0.0) result = apply(world, smooth, result);

    return result;
}

// explicit instantiation
template
class Laplacian<double, 1>;

template
class Laplacian<double, 2>;

template
class Laplacian<double, 3>;

template
class Laplacian<double, 4>;

template
class Laplacian<double, 5>;

template
class Laplacian<double, 6>;


/// ctor with an SCF calculation providing the MOs and density
template<typename T, std::size_t NDIM>
Coulomb<T, NDIM>::Coulomb(World &world, const Nemo *nemo) : world(world) {
    reset_poisson_operator_ptr(nemo->get_calc()->param.lo(), nemo->get_calc()->param.econv());
    vcoul = compute_potential(nemo);
}

/// ctor with an SCF calculation providing the MOs and density
template<typename T, std::size_t NDIM>
Coulomb<T, NDIM>::Coulomb(World &world, const SCF *calc) : world(world) {
    reset_poisson_operator_ptr(calc->param.lo(), calc->param.econv());
    vcoul = compute_potential(calc);
}

template<typename T, std::size_t NDIM>
void Coulomb<T, NDIM>::reset_poisson_operator_ptr(const double lo, const double econv) {
    poisson.reset(CoulombOperatorPtr(world, lo, econv));
}

template<typename T, std::size_t NDIM>
real_function_3d Coulomb<T, NDIM>::compute_density(const SCF *calc) const {
    real_function_3d density = calc->make_density(world, calc->get_aocc(),
                                                  calc->get_amo());
    if (calc->is_spin_restricted()) {
        density.scale(2.0);
    } else {
        real_function_3d brho = calc->make_density(world, calc->get_bocc(),
                                                   calc->get_bmo());
        density += brho;
    }
    density.truncate();
    return density;
}

template<typename T, std::size_t NDIM>
real_function_3d Coulomb<T, NDIM>::compute_potential(const madness::SCF *calc) const {
    real_function_3d density = compute_density(calc);
    return (*poisson)(density).truncate();
}

/// same as above, but with the additional factor R^2 in the density
template<typename T, std::size_t NDIM>
real_function_3d Coulomb<T, NDIM>::compute_potential(const madness::Nemo *nemo) const {
    real_function_3d density = nemo->make_density(nemo->get_calc()->aocc,
                                                  nemo->get_calc()->amo);
    if (nemo->get_calc()->is_spin_restricted()) {
        density.scale(2.0);
    } else {
        real_function_3d brho = nemo->get_calc()->make_density(world,
                                                               nemo->get_calc()->get_bocc(),
                                                               nemo->get_calc()->get_bmo());
        density += brho;
    }
    density = (density * nemo->R_square).truncate();
    return (*poisson)(density).truncate();
}


template<typename T, std::size_t NDIM>
Nuclear<T, NDIM>::Nuclear(World &world, const SCF *calc) : world(world) {
    ncf = std::shared_ptr<NuclearCorrelationFactor>(
            new PseudoNuclearCorrelationFactor(world,
                                               calc->molecule, calc->potentialmanager, 1.0));
}

template<typename T, std::size_t NDIM>
Nuclear<T, NDIM>::Nuclear(World &world, const NemoBase* nemo) : world(world) {
    ncf = nemo->ncf;
}

template<typename T, std::size_t NDIM>
Nuclear<T, NDIM>::Nuclear(World &world, const Molecule& molecule) : world(world) {
    auto pm_ptr=std::make_shared<PotentialManager>(molecule,"");
    MADNESS_CHECK(molecule.parameters.pure_ae());
    pm_ptr->make_nuclear_potential(world);
    ncf = std::shared_ptr<NuclearCorrelationFactor>(
            new PseudoNuclearCorrelationFactor(world, molecule, pm_ptr, 1.0));
}

template<typename T, std::size_t NDIM>
std::vector<Function<T, NDIM> > Nuclear<T, NDIM>::operator()(const std::vector<Function<T, NDIM> > &vket) const {

    typedef Function<T, NDIM> functionT;
    typedef std::vector<functionT> vecfuncT;

    // shortcut for local nuclear potential (i.e. no correlation factor)
    if (ncf->type() == NuclearCorrelationFactor::None) {
        return truncate(ncf->U2() * vket);
    }

    std::vector<std::shared_ptr<Derivative<T, NDIM> > > gradop =
            gradient_operator<T, NDIM>(world);
    reconstruct(world, vket);
    vecfuncT vresult = zero_functions_compressed<T, NDIM>(world, vket.size());

    // memory-saving algorithm: outer loop over the dimensions
    // apply the derivative operator on each function for each dimension
    for (std::size_t i = 0; i < NDIM; ++i) {
        vecfuncT dv = apply(world, *(gradop[i]), vket, true);
        truncate(world, dv);
        vresult += truncate(ncf->U1(i % 3) * dv);
    }

    return truncate(vresult + ncf->U2() * vket);
}


template<typename T, std::size_t NDIM>
DNuclear<T, NDIM>::DNuclear(World &world, const SCF *calc, const int iatom, const int iaxis)
        : world(world), iatom(iatom), iaxis(iaxis) {
    ncf = std::shared_ptr<NuclearCorrelationFactor>(
            new PseudoNuclearCorrelationFactor(world,
                                               calc->molecule, calc->potentialmanager, 1.0));
}

template<typename T, std::size_t NDIM>
DNuclear<T, NDIM>::DNuclear(World &world, const Nemo *nemo, const int iatom, const int iaxis)
        : world(world), iatom(iatom), iaxis(iaxis) {
    ncf = nemo->ncf;
}

template<typename T, std::size_t NDIM>
std::vector<Function<T, NDIM>> DNuclear<T, NDIM>::operator()(const std::vector<Function<T, NDIM>> &vket) const {

    const double vthresh = FunctionDefaults<NDIM>::get_thresh() * 0.1;

    // compute the U2 potential/ the derivative nuclear potential
    NuclearCorrelationFactor::U2X_functor u2x(ncf.get(), iatom, iaxis);
    real_function_3d u2x_f = real_factory_3d(world).functor(u2x)
            .thresh(vthresh).truncate_on_project();
    std::vector<Function<T, NDIM>> result = mul(world, u2x_f, vket);
    truncate(world, result, vthresh);

    // add U1 and U3 potentials if the nuclear correlation factor exists
    if (ncf->type() != NuclearCorrelationFactor::None) {

        std::vector<std::shared_ptr<Derivative<T, NDIM> > > gradop =
                gradient_operator<T, NDIM>(world);
        reconstruct(world, vket);

        // memory-saving algorithm: outer loop over the dimensions
        // apply the derivative operator on each function for each dimension
        for (std::size_t i = 0; i < NDIM; ++i) {
            std::vector<Function<T, NDIM> > dv = apply(world, *(gradop[i]), vket, true);
            truncate(world, dv);

            // note the two different axis: U1axis (i) and the derivative axis (iaxis)
            // \frac{\partial U1_i}{\partial R_{A,iaxis}}
            // e.g. d/dYA U1x
            NuclearCorrelationFactor::U1X_functor u1x(ncf.get(), iatom, i, iaxis);
            Function<double, 3> U1 = real_factory_3d(world).functor(u1x).truncate_on_project();
            std::vector<Function<T, NDIM> > U1dv = U1 * dv;
            truncate(world, U1dv);
            result = sub(world, result, U1dv);
            truncate(world, result);
        }

        // add the U3X potential
        NuclearCorrelationFactor::U3X_functor u3x(ncf.get(), iatom, iaxis);
        real_function_3d u3x_f = real_factory_3d(world).functor(u3x).truncate_on_project();
        std::vector<Function<T, NDIM> > U3v = mul(world, u3x_f, vket);
        result = sub(world, result, U3v);
        truncate(world, result);
    }
    truncate(world, result);

    return result;
}


/// custom ctor with information about the XC functional
template<typename T, std::size_t NDIM>
XCOperator<T, NDIM>::XCOperator(World &world, std::string xc_data, const bool spin_polarized,
                                const real_function_3d &arho, const real_function_3d &brho, std::string deriv)
        : world(world), dft_deriv(deriv), nbeta(0), ispin(0),
          extra_truncation(FunctionDefaults<3>::get_thresh() * 0.01) {

    nbeta = (brho.norm2() > 0.0);   // does this make sense

    xc = std::shared_ptr<XCfunctional>(new XCfunctional());
    xc->initialize(xc_data, spin_polarized, world);

    xc_args = prep_xc_args(arho, brho);
}

/// custom ctor for the regularized (nemo) path, without a Nemo object
template<typename T, std::size_t NDIM>
XCOperator<T, NDIM>::XCOperator(World &world, std::string xc_data, const bool spin_polarized,
                                const real_function_3d &arho, const real_function_3d &brho,
                                std::shared_ptr<NuclearCorrelationFactor> ncf_,
                                const real_function_3d &arho_reg_, const real_function_3d &brho_reg_,
                                std::string deriv)
        : world(world), dft_deriv(deriv), nbeta(0), ispin(0),
          extra_truncation(FunctionDefaults<3>::get_thresh() * 0.01) {

    nbeta = (brho.norm2() > 0.0);

    xc = std::shared_ptr<XCfunctional>(new XCfunctional());
    xc->initialize(xc_data, spin_polarized, world);

    ncf = ncf_;
    xc_args = prep_xc_args(arho, brho, arho_reg_, brho_reg_);
}

/// custom ctor with the XC functional
template<typename T, std::size_t NDIM>
XCOperator<T, NDIM>::XCOperator(World& world, std::shared_ptr<XCfunctional> xc,
           const bool spin_polarized,
           const int ispin,
           const int nbeta,
           const real_function_3d& arho, const real_function_3d& brho,
           std::string deriv)
    : world(world), dft_deriv(deriv), xc(xc), nbeta(nbeta), ispin(ispin),
      extra_truncation(FunctionDefaults<3>::get_thresh() * 0.01) {
  xc_args = prep_xc_args(arho, brho);
}

template<typename T, std::size_t NDIM>
XCOperator<T, NDIM>::XCOperator(World &world, const SCF *calc, int ispin, std::string deriv)
        : world(world), dft_deriv(deriv), ispin(ispin), extra_truncation(FunctionDefaults<3>::get_thresh() * 0.01) {
    xc = std::shared_ptr<XCfunctional>(new XCfunctional());
    xc->initialize(calc->param.xc(), !calc->param.spin_restricted(), world);
    nbeta = calc->param.nbeta();
    const bool have_beta = xc->is_spin_polarized() && nbeta != 0;

    // compute the alpha and beta densities
    real_function_3d arho, brho;
    arho = calc->make_density(world, calc->aocc, calc->amo);
    if (have_beta) {
        brho = calc->make_density(world, calc->bocc, calc->bmo);
    } else {
        brho = arho;
    }
    xc_args = prep_xc_args(arho, brho);
}

template<typename T, std::size_t NDIM>
XCOperator<T, NDIM>::XCOperator(World &world, const Nemo *nemo, int ispin)
        : world(world), dft_deriv(nemo->get_calc()->param.dft_deriv()), ispin(ispin),
          extra_truncation(FunctionDefaults<3>::get_thresh() * 0.01) {
    xc = std::shared_ptr<XCfunctional>(new XCfunctional());
    xc->initialize(nemo->get_calc()->param.xc(),
                   !nemo->get_calc()->param.spin_restricted(), world);

    ncf = nemo->ncf;

    ncf = nemo->ncf;

    nbeta = nemo->get_calc()->param.nbeta();
    const bool have_beta = xc->is_spin_polarized() && nbeta != 0;

    // compute the alpha and beta densities
    real_function_3d arho, brho;
    real_function_3d arhonemo = nemo->make_density(nemo->get_calc()->aocc, nemo->get_calc()->amo);
    arho = (arhonemo * nemo->R_square).truncate(extra_truncation);
    if (have_beta) {
        real_function_3d brhonemo = nemo->make_density(nemo->get_calc()->bocc, nemo->get_calc()->bmo);
        brho = (brhonemo * nemo->R_square).truncate(extra_truncation);
    } else {
        brho = arho;
    }

    xc_args = prep_xc_args(arho, brho);
}


template<typename T, std::size_t NDIM>
XCOperator<T, NDIM>::XCOperator(World &world, const SCF *calc, const real_function_3d &arho,
                                const real_function_3d &brho, int ispin, std::string deriv)
        : world(world), dft_deriv(deriv), nbeta(calc->param.nbeta()), ispin(ispin),
          extra_truncation(FunctionDefaults<3>::get_thresh() * 0.01) {
    xc = std::shared_ptr<XCfunctional>(new XCfunctional());
    xc->initialize(calc->param.xc(), !calc->param.spin_restricted(), world);
    xc_args = prep_xc_args(arho, brho);
}

template<typename T, std::size_t NDIM>
XCOperator<T, NDIM>::XCOperator(World &world, const Nemo *nemo, const real_function_3d &arho,
                                const real_function_3d &brho, int ispin)
        : world(world), dft_deriv(nemo->get_calc()->param.dft_deriv()),
          nbeta(nemo->get_calc()->param.nbeta()), ispin(ispin),
          extra_truncation(FunctionDefaults<3>::get_thresh() * 0.01) {
    xc = std::shared_ptr<XCfunctional>(new XCfunctional());
    xc->initialize(nemo->get_calc()->param.xc(),
                   not nemo->get_calc()->param.spin_restricted(), world);
    ncf = nemo->ncf;

    xc_args = prep_xc_args(arho, brho);
}

template<typename T, std::size_t NDIM>
std::vector<Function<T, NDIM> > XCOperator<T, NDIM>::operator()(const std::vector<Function<T, NDIM> > &vket) const {
    real_function_3d xc_pot = make_xc_potential();
    double vtol = FunctionDefaults<3>::get_thresh() * 0.1;  // safety
    std::vector<Function<T, NDIM> > result = mul_sparse(world, xc_pot, vket, vtol);
    if (has_tau_term()) result += apply_tau_term(vket);
    return result;
}


template<typename T, std::size_t NDIM>
bool XCOperator<T, NDIM>::has_tau_term() const {
    return xc->needs_tau();
}


template<typename T, std::size_t NDIM>
real_function_3d XCOperator<T, NDIM>::get_tau(const int spin) const {
    return xc_args[spin == 0 ? XCfunctional::enum_taua : XCfunctional::enum_taub];
}


/// compute tau = 1/2 sum_i |grad psi_i|^2 and store it in the intermediates
template<typename T, std::size_t NDIM>
void XCOperator<T, NDIM>::set_tau(const vecfuncT &amo, const Tensor<double> &aocc,
                                  const vecfuncT &bmo, const Tensor<double> &bocc,
                                  const TauU1 u1mode) const {

    MADNESS_CHECK_THROW(is_initialized(), "set_tau called before the intermediates exist");

    // In nemo mode the vectors handed in are the nemos F, with psi = R F, so
    //
    //   grad psi = R (grad F - U1 F),        U1 = -grad(R)/R
    //   |grad psi|^2 = R^2 (|grad F|^2 - 2 F U1.grad F + |U1|^2 F^2)
    //
    // and the whole nuclear cusp sits in U1, which is analytic and precomputed.
    // Only the cusp-free F is differentiated numerically. Forming psi = R F and
    // differentiating that instead would put the cusp straight back under the
    // derivative operator, which is what the regularization exists to avoid.
    // Same decomposition as OEP::compute_total_kinetic_density.
    // Which of the three routes is taken decides what has to be projected: the
    // pointwise route touches no U1 Function at all, and skipping the projection
    // is half its point -- U1vec() and U1_dot_U1 are four functor projections per
    // SCF iteration, each of them deep.
    vecfuncT U1;
    real_function_3d U1dot, R_square;
    const bool u1_as_functions = bool(ncf) and (u1mode == TauU1::mra);
    if (ncf) {
        if (u1_as_functions) {
            U1 = ncf->U1vec();
            NuclearCorrelationFactor::U1_dot_U1_functor u1_dot_u1(ncf.get());
            U1dot = real_factory_3d(world).functor(u1_dot_u1).truncate_on_project();
        }
        R_square = ncf->square();
    }

    const bool have_beta = (xc->is_spin_polarized()) and (nbeta > 0);

    // tau_sigma = 1/2 sum_i occ_{sigma,i} |grad psi_i|^2. The occupations are not
    // decoration: the caller's orbital vectors are sized nmo, not nalpha/nbeta, so
    // they carry virtual orbitals whose occupation is zero. Those add nothing to
    // the density but would inflate an unweighted sum, silently changing the
    // meta-gga energy and potential as soon as virtuals are requested. Fractional
    // occupations need the same weighting, exactly as make_density applies it.
    //
    // Spin bookkeeping: madness stores per-spin occupations, one per spin channel
    // even when spin-restricted (SCF.cc: "madness instead stores 2 identical sets
    // (alpha and beta) with occupation 1"), and make_libxc_args forms the total
    // from the alpha quantities. So the weight is the occupation itself, with no
    // further normalisation, and the usual occ = 1 reproduces the unweighted sum.
    // the smooth ingredients of the product rule, all cusp-free by construction
    struct tau_pieces {
        real_function_3d gradf;      ///< sum_i w_i |grad F_i|^2
        real_function_3d n;          ///< sum_i w_i F_i^2
        real_function_3d G[3];       ///< sum_i w_i F_i dF_i/dx_a = 1/2 grad(n)
    };

    auto compute_tau = [&](const vecfuncT &mo, const Tensor<double> &occ) -> tau_pieces {
        MADNESS_CHECK_THROW(occ.size() >= long(mo.size()),
                            "set_tau: fewer occupation numbers than orbitals");

        // fold the weight into the orbitals as sqrt(w): the derivative is linear,
        // so dot(D(sqrt(w) psi), D(sqrt(w) psi)) is sum_i w_i |grad psi_i|^2 and
        // the vectorised path is preserved. w == 1 is passed through untouched, so
        // integer-occupied cases are bit-identical to an unweighted sum rather
        // than picking up the noise of a redundant scalar multiplication.
        vecfuncT wmo;
        for (size_t i = 0; i < mo.size(); ++i) {
            const double w = occ(long(i));
            if (w == 0.0) continue;             // virtuals carry no density, no tau
            MADNESS_CHECK_THROW(w > 0.0, "set_tau: negative occupation number");
            wmo.push_back(w == 1.0 ? mo[i] : std::sqrt(w) * mo[i]);
        }

        tau_pieces p;
        if (wmo.empty()) {
            p.gradf = real_factory_3d(world).compressed();
            return p;
        }
        p.gradf = real_factory_3d(world).compressed();
        p.n = dot(world, wmo, wmo);
        for (int axis = 0; axis < 3; ++axis) {
            real_derivative_3d D(world, axis);
            if (dft_deriv == "bspline") D.set_bspline1();
            else if (dft_deriv == "ble") D.set_ble1();
            vecfuncT mo_copy = copy(world, wmo);
            refine(world, mo_copy);
            vecfuncT dmo = apply(world, D, mo_copy);
            p.gradf += dot(world, dmo, dmo);
            // G_a = sum_i w_i F_i dF_i/dx_a = 1/2 dn/dx_a. Smooth: F is cusp-free.
            if (ncf) p.G[axis] = dot(world, mo_copy, dmo);
        }
        p.gradf.truncate(extra_truncation);
        if (p.n.is_initialized()) p.n.truncate(extra_truncation);
        for (int axis = 0; axis < 3; ++axis)
            if (p.G[axis].is_initialized()) p.G[axis].truncate(extra_truncation);
        return p;
    };

    // Where the two routes part. The pointwise route stores the *pieces* and lets
    // make_libxc_args contract them against U1 at the quadrature points, so nothing
    // involving U1 is ever projected and tau never exists as a Function at all --
    // which is also why its depth stops taxing every other intermediate through
    // refine_to_common_level. The mra route (and the moldft path, which has no ncf
    // and no U1 to worry about) assembles tau here as before.
    const tau_pieces pa = compute_tau(amo, aocc);
    tau_pieces pb;
    if (have_beta) {
        MADNESS_CHECK_THROW(bmo.size() > 0, "set_tau needs beta orbitals for an "
                                            "open-shell meta-gga calculation");
        pb = compute_tau(bmo, bocc);
    }

    auto assemble = [&](const tau_pieces& p) {
        real_function_3d r = copy(p.gradf);
        if (u1_as_functions) {
            r += U1dot * p.n;
            for (int axis = 0; axis < 3; ++axis) r -= 2.0 * U1[axis] * p.G[axis];
        }
        if (ncf) r = r * R_square;
        // libxc convention: tau = 1/2 sum_i |grad psi_i|^2  (since libxc 2.0.0)
        return (0.5 * r).truncate(extra_truncation);
    };

    if (ncf and u1mode == TauU1::pointwise) {
        xc_args[XCfunctional::enum_nemo_R2] = R_square;
        xc_args[XCfunctional::enum_gradfa] = pa.gradf;
        xc_args[XCfunctional::enum_na] = pa.n;
        xc_args[XCfunctional::enum_Ga_x] = pa.G[0];
        xc_args[XCfunctional::enum_Ga_y] = pa.G[1];
        xc_args[XCfunctional::enum_Ga_z] = pa.G[2];
        if (have_beta) {
            xc_args[XCfunctional::enum_gradfb] = pb.gradf;
            xc_args[XCfunctional::enum_nb] = pb.n;
            xc_args[XCfunctional::enum_Gb_x] = pb.G[0];
            xc_args[XCfunctional::enum_Gb_y] = pb.G[1];
            xc_args[XCfunctional::enum_Gb_z] = pb.G[2];
        }
    } else {
        xc_args[XCfunctional::enum_taua] = assemble(pa);
        if (have_beta) xc_args[XCfunctional::enum_taub] = assemble(pb);
    }
    world.gop.fence();
}


/// apply the non-multiplicative meta-gga term, -1/2 sum_x D_x(vtau D_x psi_i)
template<typename T, std::size_t NDIM>
std::vector<Function<T, NDIM> >
XCOperator<T, NDIM>::apply_tau_term(const std::vector<Function<T, NDIM> > &vket) const {

    MADNESS_CHECK_THROW(has_tau_term(), "apply_tau_term on a non-meta functional");
    MADNESS_CHECK_THROW(vtau.is_initialized(), "apply_tau_term before make_xc_potential");

    const double vtol = FunctionDefaults<3>::get_thresh() * 0.1;

    // 1 + de/dtau is the inverse effective mass of the modified kinetic operator
    // -1/2 nabla.((1 + de/dtau) nabla). Where it goes non-positive the quadratic
    // form loses positive-definiteness and there is no minimum to converge to.
    // Function carries no min/max, so this samples a coarse lattice: it is an
    // indicator, not a bound, hence a warning rather than an assertion.
    if (print_level >= 2) {
        const double L = FunctionDefaults<3>::get_cell_width().max() * 0.5;
        double lo = 1.e10, hi = -1.e10;
        const int n = 7;
        for (int ix = 0; ix < n; ++ix)
            for (int iy = 0; iy < n; ++iy)
                for (int iz = 0; iz < n; ++iz) {
                    const coord_3d r{-L + 2.0 * L * ix / (n - 1),
                                     -L + 2.0 * L * iy / (n - 1),
                                     -L + 2.0 * L * iz / (n - 1)};
                    const double v = vtau(r);
                    lo = std::min(lo, v);
                    hi = std::max(hi, v);
                }
        if (world.rank() == 0) {
            print("meta-gga de/dtau sampled on a", n, "^3 lattice: min", lo, " max", hi);
            if (1.0 + lo <= 0.0)
                print("WARNING: 1 + de/dtau is not positive -- the effective mass "
                      "operator is not elliptic there and the SCF may not converge");
        }
    }

    // In nemo mode the kets are the nemos F and what the caller needs back is the
    // term divided by R, since its equations are for F rather than psi = R F:
    //
    //   R^{-1} [ -1/2 div( v_tau grad(R F) ) ]
    //
    // With grad(R F) = R (grad F - U1 F) and W = v_tau (grad F - U1 F),
    //
    //   div(R W) = R div W + grad(R).W = R (div W - U1.W)
    //
    // so the R factors cancel exactly and the result is -1/2 (div W - U1.W). No
    // division by R anywhere, and the cusp stays in the analytic U1 rather than
    // under a derivative operator. Without an ncf, U1 drops out and W = v_tau
    // grad(psi), which is the plain nested form.
    vecfuncT U1;
    if (ncf) U1 = ncf->U1vec();

    std::vector<Function<T, NDIM> > result =
            zero_functions_compressed<T, NDIM>(world, vket.size());
    for (int axis = 0; axis < 3; ++axis) {
        auto D = make_derivative(axis);
        std::vector<Function<T, NDIM> > vket_copy = copy(world, vket);
        refine(world, vket_copy);
        std::vector<Function<T, NDIM> > W = apply(world, *D, vket_copy);
        if (ncf) W = sub(world, W, mul(world, U1[axis], vket_copy));
        // vtau is only ever multiplied, never differentiated
        W = mul_sparse(world, vtau, W, vtol);
        refine(world, W);
        result = add(world, result, apply(world, *D, W));
        if (ncf) result = sub(world, result, mul(world, U1[axis], W));
    }
    scale(world, result, -0.5);
    truncate(world, result);
    return result;
}


/// gradient operator for the meta-gga term, honouring dft_deriv
template<typename T, std::size_t NDIM>
std::shared_ptr<Derivative<T, NDIM> > XCOperator<T, NDIM>::make_derivative(const int axis) const {
    auto D = std::shared_ptr<Derivative<T, NDIM> >(new Derivative<T, NDIM>(world, axis));
    if (dft_deriv == "bspline") D->set_bspline1();
    else if (dft_deriv == "ble") D->set_ble1();
    return D;
}


/// divergence of a real vector field, honouring dft_deriv

/// vmra.h's div() is div_abgv(), so a caller that has selected
/// `dft_deriv bspline`/`ble` has to name the matching derivative explicitly.
template<typename T, std::size_t NDIM>
real_function_3d XCOperator<T, NDIM>::div_dft_deriv(const vecfuncT& v) const {
    const DerivMethod method = (dft_deriv == "bspline") ? DerivMethod::bspline
                             : (dft_deriv == "ble")     ? DerivMethod::ble
                                                        : DerivMethod::abgv;
    return madness::div_deriv(v, method, true);
}


/// the xc contribution to the Fock matrix

/// See the declaration for the bra convention -- vbra carries R^2, vket does not.
template<typename T, std::size_t NDIM>
Tensor<T> XCOperator<T, NDIM>::operator()(const std::vector<Function<T,NDIM> >& vbra,
                                          const std::vector<Function<T,NDIM> >& vket) const {

    MADNESS_CHECK_THROW(vbra.size() == vket.size(),
                        "XCOperator matrix elements: vbra and vket differ in size");
    MADNESS_CHECK_THROW(vlocal.is_initialized(),
                        "XCOperator matrix elements before make_xc_potential");

    const long nn = long(vket.size());
    Tensor<T> result(nn, nn);
    if (nn == 0) return result;

    const double thresh = FunctionDefaults<3>::get_thresh();
    const double vtol = thresh * 0.1;

    // Verify the bra weighting instead of trusting it. The tau term below
    // distinguishes bra from ket, so a swapped or unweighted bra does not fail --
    // it returns a plausible wrong matrix, which is the worst outcome available.
    // One orbital settles it, so the check is O(1) in the orbital count, and R^2
    // appears here only as an assertion and nowhere in the arithmetic.
    {
        std::vector<Function<T,NDIM> > one(1, vket[0]), ref;
        if (ncf) ref = mul_sparse(world, ncf->square(), one, vtol);
        else     ref = one;
        const double err = (vbra[0] - ref[0]).norm2();
        MADNESS_CHECK_THROW(err < 10.0 * thresh * (1.0 + vbra[0].norm2()),
            ncf ? "XCOperator matrix elements: vbra must be the R^2-weighted "
                  "counterpart of vket -- call it as xcoperator(R2nemo, nemo), the "
                  "way Kinetic is called in Nemo::compute_fock_matrix"
                : "XCOperator matrix elements: without a nuclear correlation factor "
                  "vbra and vket must be the same orbitals");
    }

    // the multiplicative potential is complete here, so all that is missing is the
    // non-multiplicative meta-gga term. apply_tau_term returns it already divided
    // by R, which is exactly what makes an R^2-weighted bra produce the physical
    // <psi_i|.|psi_j>.
    result += matrix_inner(world, mul_sparse(world, vlocal, vbra, vtol), vket);
    if (has_tau_term()) result += matrix_inner(world, vbra, apply_tau_term(vket));
    return result;
}


template<typename T, std::size_t NDIM>
double XCOperator<T, NDIM>::compute_xc_energy() const {

    if (not is_initialized()) {
        MADNESS_EXCEPTION("calling xc energy without intermediates ", 1);
    }
    // same precondition as make_xc_potential(): without it a meta-gga energy is
    // evaluated at the tau floor instead of the orbital tau, which is wrong but
    // finite and therefore easy to miss
    if (has_tau_term() and (not has_tau_args())) {
        MADNESS_EXCEPTION("meta-gga functional without a kinetic energy density: "
                          "call XCOperator::set_tau() with the occupied orbitals "
                          "before compute_xc_energy()", 1);
    }

    refine_to_common_level(world, xc_args);
    real_function_3d vlda = multiop_values<double, xc_functional, 3>
            (xc_functional(*xc, make_u1_functors()), xc_args);
    truncate(world, xc_args);

    return vlda.trace();
}


/// true once set_tau() has supplied tau, by either route
template<typename T, std::size_t NDIM>
bool XCOperator<T, NDIM>::has_tau_args() const {
    if (xc_args[XCfunctional::enum_taua].is_initialized()) return true;      // moldft / mra route
    return xc_args[XCfunctional::enum_gradfa].is_initialized()               // pointwise route
       and xc_args[XCfunctional::enum_nemo_R2].is_initialized();
}


/// the four analytic U1 quantities the xc ops evaluate pointwise
template<typename T, std::size_t NDIM>
nemo_u1_functors XCOperator<T, NDIM>::make_u1_functors() const {
    typedef FunctionFunctorInterface<double, 3> functorT;
    if (not ncf) return nemo_u1_functors();
    if (not xc_args[XCfunctional::enum_gradfa].is_initialized()) return nemo_u1_functors();
    std::vector<std::shared_ptr<functorT> > f;
    for (int axis = 0; axis < 3; ++axis)
        f.push_back(std::shared_ptr<functorT>(
                new NuclearCorrelationFactor::U1_functor(ncf.get(), axis)));
    f.push_back(std::shared_ptr<functorT>(
            new NuclearCorrelationFactor::U1_dot_U1_functor(ncf.get())));
    return nemo_u1_functors(f);
}


/// A thin wrapper whose only job is to stash the result: operator()(vbra,vket)
/// needs whatever the caller got, with no argument to receive it through.
template<typename T, std::size_t NDIM>
real_function_3d XCOperator<T, NDIM>::make_xc_potential() const {
    vlocal = make_xc_potential_impl();
    return vlocal;
}


template<typename T, std::size_t NDIM>
real_function_3d XCOperator<T, NDIM>::make_xc_potential_impl() const {

    if (not is_initialized()) {
        MADNESS_EXCEPTION("calling xc potential without intermediates ", 1);
    }
    if (has_tau_term() and (not has_tau_args())) {
        MADNESS_EXCEPTION("meta-gga functional without a kinetic energy density: "
                          "call XCOperator::set_tau() with the occupied orbitals "
                          "before make_xc_potential()", 1);
    }

    refine_to_common_level(world, xc_args);

    // compute all the contributions to the xc kernel
    xc_potential op(*xc, ispin, make_u1_functors());
    const vecfuncT intermediates = multi_to_multi_op_values(op, xc_args);

    // local part, first term in Yanai2005, Eq. (12)
    real_function_3d dft_pot = intermediates[0];

    // de/dtau -- kept for apply_tau_term, which turns it into the
    // non-multiplicative operator. Comes out of the same pointwise pass.
    if (has_tau_term()) vtau = intermediates[xc->is_spin_polarized() ? 7 : 4];

    if (xc->needs_sigma()) {
        vecfuncT semilocal(3);
        semilocal[0] = intermediates[1];
        semilocal[1] = intermediates[2];
        semilocal[2] = intermediates[3];

        // second term in Yanai2005, Eq. (12)
        real_function_3d gga_pot_same_spin = div_dft_deriv(semilocal);
        dft_pot -= gga_pot_same_spin;

        bool have_beta = xc->is_spin_polarized() && nbeta != 0;

        if (have_beta) {
            semilocal[0] = intermediates[4];
            semilocal[1] = intermediates[5];
            semilocal[2] = intermediates[6];

            // third term in Yanai2005, Eq. (12)
            real_function_3d gga_pot_other_spin = div_dft_deriv(semilocal);
            dft_pot -= gga_pot_other_spin;
        }
    }

    truncate(world, xc_args);
    return dft_pot.truncate();
}


/// apply the xc kernel on a perturbed density

/// cf Eq. (13) of T. Yanai, R. J. Harrison, and N. Handy,
/// “Multiresolution quantum chemistry in multiwavelet bases: time-dependent
/// density functional theory with asymptotically corrected potentials in
/// local density and generalized gradient approximations,”
/// Mol. Phys., vol. 103, no. 2, pp. 413–424, 2005.
///
/// the application of the xc kernel is (RHF only)
/// \f[
///   \frac{\partial^2E_{xc}}{\partial \rho_\alpha^2}\circ\tilde\rho
///      = second_{local} + second_{semilocal} + first_{semilocal}
/// \f]
/// where the second partial derivatives are
/// \f[
///        second_{local} = \frac{\partial^2 f_{xc}}{\partial \rho_\alpha^2}\tilde \rho
///        + 2\frac{\partial^2 f_{xc}}{\partial \rho_\alpha\sigma_{\alpha\alpha}}
///            \left(\vec\nabla \rho_a\cdot \vec \nabla\tilde\rho\right)
/// \f]
///  the second partial derivatives that need to be multiplied with the density gradients
/// \f[
///      second_{semilocal} = -\vec\nabla\cdot\left((\vec\nabla\rho)
///             \left[2\frac{\partial^2 f_{xc}}{\partial\rho_\alpha\partial\sigma_{\alpha\alpha}}\tilde\rho
///             + 4\frac{\partial^2 f_{xc}}{\partial\sigma_{\alpha\alpha}^2}
///                \left(\vec\nabla\rho_\alpha\cdot\vec\nabla\tilde\rho\right)\right]\right)
/// \f]
/// and the first derivatives that need to be multiplied with the density gradients
/// \f[
///      first_{semilocal} =
///        -\vec\nabla\cdot\left(2\frac{\partial f_{xc}}{\partial\sigma_{\alpha\alpha}}\vec\nabla\tilde\rho\right)
/// \f]
template<typename T, std::size_t NDIM>
real_function_3d XCOperator<T, NDIM>::apply_xc_kernel(const real_function_3d &dens_pt,
                                                      const vecfuncT grad_dens_pt) const {

    MADNESS_ASSERT(not xc->is_spin_polarized());    // for now
    MADNESS_ASSERT(ispin == 0);           // for now

    if (not is_initialized()) {
        MADNESS_EXCEPTION("calling apply_xc_kernel without intermediates ", 1);
    }

    vecfuncT ddens_pt = grad_dens_pt;
    prep_xc_args_response(dens_pt, xc_args, ddens_pt);
    refine_to_common_level(world, xc_args);

    // compute all the contributions to the xc kernel
    xc_kernel_apply op(*xc, ispin);
    const vecfuncT intermediates = multi_to_multi_op_values(op, xc_args);

    // lda potential and local parts of the gga potential
    real_function_3d result = intermediates[0];

    // add semilocal gga potentials
    if (xc->is_gga()) {
        // turn intermediates into quantities that can be digested by the div operator
        vecfuncT semilocal(3);
        semilocal[0] = intermediates[1];
        semilocal[1] = intermediates[2];
        semilocal[2] = intermediates[3];

        real_function_3d gga_pot = -1.0 * div_dft_deriv(semilocal);

        result += gga_pot;
    }
    truncate(world, xc_args);
    return result.truncate();
}

/// prepare xc args
template<typename T, std::size_t NDIM>
vecfuncT XCOperator<T, NDIM>::prep_xc_args(const real_function_3d &arho,
                                           const real_function_3d &brho,
                                           const real_function_3d &arho_reg_in,
                                           const real_function_3d &brho_reg_in) const {

    World &world = arho.world();
    vecfuncT xcargs(XCfunctional::number_xc_args);
    const bool have_beta = (xc->is_spin_polarized()) and (nbeta > 0);

    // assign the densities (alpha, beta)
    xcargs[XCfunctional::enum_rhoa] = copy(arho.reconstruct());      // alpha density
    if (have_beta) xcargs[XCfunctional::enum_rhob] = copy(brho.reconstruct());  // beta density
    world.gop.fence();

    // zeta_sigma = grad(ln rho_sigma), so that grad(rho_sigma) = rho_sigma zeta_sigma
    // and sigma_st = rho_s rho_t (zeta_s.zeta_t). Only zeta is stored: the
    // contractions zeta_s.zeta_t are formed pointwise in
    // XCfunctional::make_libxc_args, where they are guaranteed to be the exact Gram
    // matrix of the gradients. Carrying them as their own multiwavelet functions
    // (as this used to) meant the projected product disagreed with the zeta
    // components at the quadrature points, by O(1) near the nuclear cusp -- enough
    // for the sum of squares chi_aa to turn negative and for the total sigma handed
    // to libxc to follow it.
    if (xc->needs_sigma()) {

        auto grad_variant = [&](const real_function_3d& f) {
            if (dft_deriv == "bspline") return grad_bspline_one(f);  // b-spline
            if (dft_deriv == "ble")     return grad_ble_one(f);      // BLE
            return grad(f);                                          // default is abgv
        };

        // zeta = grad log(rho). With a nuclear correlation factor rho = R^2 rho_reg,
        // so zeta = grad log(R^2) + grad log(rho_reg) = -2 U1 + grad log(rho_reg)
        // (U1 = -grad(R)/R). Differentiating log(rho) directly puts the nuclear cusp
        // under the derivative operator, which is what the regularization exists to
        // avoid: the kink forces refinement to the finest level and the O(thresh)
        // noise it leaves there is amplified by 2^n by the derivative. Only the
        // cusp-free rho_reg is differentiated here.
        const bool regularized_zeta = bool(ncf);
        vecfuncT U1;
        if (regularized_zeta) U1 = ncf->U1vec();

        auto make_zeta = [&](const real_function_3d& rho,
                             const real_function_3d& rho_reg) {
            if (regularized_zeta and rho_reg.is_initialized()) {
                real_function_3d logdens = unary_op(rho_reg, logme());
                vecfuncT zeta = grad_variant(logdens);
                for (int axis = 0; axis < 3; ++axis) zeta[axis] -= 2.0 * U1[axis];
                return zeta;
            }
            real_function_3d logdens = unary_op(rho, logme());
            return grad_variant(logdens);
        };

        vecfuncT grada = make_zeta(arho, arho_reg_in);
        xcargs[XCfunctional::enum_zetaa_x] = grada[0];
        xcargs[XCfunctional::enum_zetaa_y] = grada[1];
        xcargs[XCfunctional::enum_zetaa_z] = grada[2];

        if (have_beta) {
            vecfuncT gradb = make_zeta(brho, brho_reg_in);
            xcargs[XCfunctional::enum_zetab_x] = gradb[0];
            xcargs[XCfunctional::enum_zetab_y] = gradb[1];
            xcargs[XCfunctional::enum_zetab_z] = gradb[2];
        }
    }

    world.gop.fence();
    truncate(world, xcargs, extra_truncation);
    return xcargs;
}

/// add intermediates for the response kernels to xc_args
template<typename T, std::size_t NDIM>
void XCOperator<T, NDIM>::prep_xc_args_response(const real_function_3d &dens_pt,
                                                vecfuncT &xc_args, vecfuncT &ddens_pt) const {

    const bool have_beta = (xc->is_spin_polarized()) and (nbeta > 0);
    World &world = dens_pt.world();

    // assign the perturbed density (spin-free)
    xc_args[XCfunctional::enum_rho_pt] = dens_pt;
    world.gop.fence();

    // assign the reduced density gradients with the perturbed density for GGA
    // \sigma_pt   = 2.0 * \nabla \rho_\alpha \cdot \nabla\tilde\rho
    // \sigma_pt_a = \nabla \rho_\alpha \cdot \nabla\tilde\rho
    // \sigma_pt_b = \nabla \rho_\beta \cdot \nabla\tilde\rho
    //
    // using the logarithmic derivatives for rho only we get (alpha and RHF)
    // \sigma_pt = 2.0 * \rho_\alpha (\nabla\zeta_\alpha \cdot \nabla\tilde\rho)
    // \sigma_pt_a = \rho_\alpha (\nabla\zeta_\alpha \cdot \nabla\tilde\rho)
    // we save the functions without multiplying the ground state density rho
    if (xc->is_gga()) {

        if (ddens_pt.size() == 0) ddens_pt = grad(dens_pt);     // spin free
        else print(" using provided ddens_pt in prep_xc_args_response");

        xc_args[XCfunctional::enum_ddens_ptx] = ddens_pt[0];
        xc_args[XCfunctional::enum_ddens_pty] = ddens_pt[1];
        xc_args[XCfunctional::enum_ddens_ptz] = ddens_pt[2];

        std::vector<real_function_3d> zeta(3);
        zeta[0] = xc_args[XCfunctional::enum_zetaa_x];
        zeta[1] = xc_args[XCfunctional::enum_zetaa_y];
        zeta[2] = xc_args[XCfunctional::enum_zetaa_z];
        xc_args[XCfunctional::enum_sigma_pta_div_rho] = dot(world, zeta, ddens_pt);    // sigma_a
        // for RHF add factor 2 on rho; will be done in xcfunctional_libxc::make_libxc_args
        // \sigma_pt = 2 * rho_a * sigma_pta_div_rho
        world.gop.fence();

        if (have_beta) {
            zeta[0] = xc_args[XCfunctional::enum_zetab_x];
            zeta[1] = xc_args[XCfunctional::enum_zetab_y];
            zeta[2] = xc_args[XCfunctional::enum_zetab_z];
            xc_args[XCfunctional::enum_sigma_ptb_div_rho] = dot(world, zeta, ddens_pt);  // sigma_b
        }
        world.gop.fence();
    }
    world.gop.fence();
    truncate(world, xc_args, extra_truncation);
}

/// ctor
template<typename T, std::size_t NDIM>
Exchange<T,NDIM>::Exchange(World& world, const double lo, const double thresh) : impl(new Exchange<T,NDIM>::ExchangeImpl(world,lo,thresh)) {};


/// ctor with a conventional calculation
template<typename T, std::size_t NDIM>
Exchange<T,NDIM>::Exchange(World& world, const SCF *calc, const int ispin) : impl(new Exchange<T,NDIM>::ExchangeImpl(world,calc,ispin)) {};

/// ctor with a nemo calculation
template<typename T, std::size_t NDIM>
Exchange<T,NDIM>::Exchange(World& world, const Nemo *nemo, const int ispin) : impl(new Exchange<T,NDIM>::ExchangeImpl(world,nemo,ispin)) {};

/// apply the exchange operator on a vector of functions

/// note that only one spin is used (either alpha or beta orbitals)
/// @param[in]  vket       the orbitals |i> that the operator is applied on
/// @return     a vector of orbitals  K| i>
template<typename T, std::size_t NDIM>
std::vector<Function<T,NDIM>> Exchange<T,NDIM>::operator()(const std::vector<Function<T,NDIM>>& vket) const {
    impl->set_taskq(this->taskq);
    auto result=impl->operator()(vket);
    this->statistics=impl->get_statistics();
    return result;
};

template<typename T, std::size_t NDIM>
Exchange<T,NDIM>& Exchange<T,NDIM>::set_bra_and_ket(const vecfuncT& bra, const vecfuncT& ket) {
    MADNESS_CHECK(impl);
    impl->set_bra_and_ket(bra, ket);
    return *this;
}

template<typename T, std::size_t NDIM>
bool Exchange<T,NDIM>::is_symmetric() const {
    return impl->is_symmetric();
}

template<typename T, std::size_t NDIM>
Exchange<T,NDIM>& Exchange<T,NDIM>::set_symmetric(const bool flag) {
    impl->symmetric(flag);
    return *this;
}

template<typename T, std::size_t NDIM>
Exchange<T,NDIM>& Exchange<T,NDIM>::set_algorithm(const ExchangeAlgorithm& alg) {
    impl->set_algorithm(alg);
    return *this;
}

template<typename T, std::size_t NDIM>
Exchange<T,NDIM>& Exchange<T,NDIM>::set_macro_task_info(const MacroTaskInfo& info) {
    impl->set_macro_task_info(info);
    return *this;
}

 template<typename T, std::size_t NDIM>
 Exchange<T,NDIM>& Exchange<T,NDIM>::set_printlevel(const long& level) {
    impl->set_printlevel(level);
    return *this;
}

template<typename T, std::size_t NDIM>
Exchange<T,NDIM>& Exchange<T,NDIM>::set_batch_granularity(const long level) {
    impl->set_batch_granularity(level);
    return *this;
}

template<typename T, std::size_t NDIM>
Exchange<T,NDIM>& Exchange<T,NDIM>::set_accumulation_mode(const int mode) {
    impl->set_accumulation_mode(mode);
    return *this;
}

template<typename T, std::size_t NDIM>
Exchange<T,NDIM>& Exchange<T,NDIM>::set_cost_aware_assignment(const bool flag) {
    impl->set_cost_aware_assignment(flag);
    return *this;
}

template<>
Fock<double, 3>::Fock(World &world, const Nemo *nemo) : world(world) {
    auto tmp = nemo->make_fock_operator();
    if (tmp) std::swap(tmp->operators, operators);
    else MADNESS_EXCEPTION("failed to construct fock operator", 1);
}

template<>
Fock<double, 3>::Fock(World &world, const OEP *oep) : world(world) {
    auto tmp = oep->make_fock_operator();
    if (tmp) std::swap(tmp->operators, operators);
    else MADNESS_EXCEPTION("failed to construct fock operator", 1);
}

template<>
Fock<double, 3>::Fock(World &world, const NemoBase *nemobase) : world(world) {
    auto tmp = nemobase->make_fock_operator();
    if (tmp) std::swap(tmp->operators, operators);
    else MADNESS_EXCEPTION("failed to construct fock operator", 1);
}


template class Exchange<double_complex,3>;
template class Exchange<double,3>;

template class Coulomb<double_complex,3>;
template class Coulomb<double,3>;

template class XCOperator<double_complex,3>;
template class XCOperator<double,3>;

template class Nuclear<double_complex,3>;
template class Nuclear<double,3>;

template class DNuclear<double_complex,3>;
template class DNuclear<double,3>;

template class Fock<double_complex,3>;
template class Fock<double,3>;

} // namespace madness


