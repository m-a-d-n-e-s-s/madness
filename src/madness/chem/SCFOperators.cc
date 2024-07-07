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
        : world(world), ispin(ispin), extra_truncation(FunctionDefaults<3>::get_thresh() * 0.01) {
    xc = std::shared_ptr<XCfunctional>(new XCfunctional());
    xc->initialize(nemo->get_calc()->param.xc(),
                   !nemo->get_calc()->param.spin_restricted(), world);

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
        : world(world), nbeta(nemo->get_calc()->param.nbeta()), ispin(ispin), extra_truncation(0.01) {
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
    return mul_sparse(world, xc_pot, vket, vtol);
}

template<typename T, std::size_t NDIM>
double XCOperator<T, NDIM>::compute_xc_energy() const {

    if (not is_initialized()) {
        MADNESS_EXCEPTION("calling xc energy without intermediates ", 1);
    }

    refine_to_common_level(world, xc_args);
    real_function_3d vlda = multiop_values<double, xc_functional, 3>
            (xc_functional(*xc), xc_args);
    truncate(world, xc_args);

    return vlda.trace();
}


template<typename T, std::size_t NDIM>
real_function_3d XCOperator<T, NDIM>::make_xc_potential() const {

    if (not is_initialized()) {
        MADNESS_EXCEPTION("calling xc potential without intermediates ", 1);
    }

    refine_to_common_level(world, xc_args);

    // compute all the contributions to the xc kernel
    xc_potential op(*xc, ispin);
    const vecfuncT intermediates = multi_to_multi_op_values(op, xc_args);

    // local part, first term in Yanai2005, Eq. (12)
    real_function_3d dft_pot = intermediates[0];

    if (xc->is_gga()) {
        vecfuncT semilocal(3);
        semilocal[0] = intermediates[1];
        semilocal[1] = intermediates[2];
        semilocal[2] = intermediates[3];

        // second term in Yanai2005, Eq. (12)
        real_function_3d gga_pot_same_spin = div(semilocal, true);
        dft_pot -= gga_pot_same_spin;

        bool have_beta = xc->is_spin_polarized() && nbeta != 0;

        if (have_beta) {
            semilocal[0] = intermediates[4];
            semilocal[1] = intermediates[5];
            semilocal[2] = intermediates[6];

            // third term in Yanai2005, Eq. (12)
            real_function_3d gga_pot_other_spin = div(semilocal, true);
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

        real_function_3d gga_pot = -1.0 * div(semilocal, true);

        result += gga_pot;
    }
    truncate(world, xc_args);
    return result.truncate();
}

/// prepare xc args
template<typename T, std::size_t NDIM>
vecfuncT XCOperator<T, NDIM>::prep_xc_args(const real_function_3d &arho,
                                           const real_function_3d &brho) const {

    World &world = arho.world();
    vecfuncT xcargs(XCfunctional::number_xc_args);
    const bool have_beta = (xc->is_spin_polarized()) and (nbeta > 0);

    // assign the densities (alpha, beta)
    xcargs[XCfunctional::enum_rhoa] = copy(arho.reconstruct());      // alpha density
    if (have_beta) xcargs[XCfunctional::enum_rhob] = copy(brho.reconstruct());  // beta density
    world.gop.fence();

    // compute the chi quantity such that sigma = rho^2 * chi
    if (xc->is_gga()) {

        real_function_3d logdensa = unary_op(arho, logme());
        vecfuncT grada;
        if (dft_deriv == "bspline") grada = grad_bspline_one(logdensa); // b-spline
        else if (dft_deriv == "ble") grada = grad_ble_one(logdensa);    // BLE
        else grada = grad(logdensa);                                   // Default is abgv
        real_function_3d chi = dot(world, grada, grada);
        xcargs[XCfunctional::enum_chi_aa] = chi;
        xcargs[XCfunctional::enum_zetaa_x] = grada[0];
        xcargs[XCfunctional::enum_zetaa_y] = grada[1];
        xcargs[XCfunctional::enum_zetaa_z] = grada[2];

        if (have_beta) {
            real_function_3d logdensb = unary_op(brho, logme());
            // Bryan's edits for derivatives
            vecfuncT gradb;
            if (dft_deriv == "bspline") gradb = grad_bspline_one(logdensb);  // b-spline
            else if (dft_deriv == "ble") gradb = grad_ble_one(logdensb);     // BLE
            else gradb = grad(logdensb);                                    // Default is abgv
            real_function_3d chib = dot(world, gradb, gradb);
            real_function_3d chiab = dot(world, grada, gradb);
            xcargs[XCfunctional::enum_zetab_x] = gradb[0];
            xcargs[XCfunctional::enum_zetab_y] = gradb[1];
            xcargs[XCfunctional::enum_zetab_z] = gradb[2];
            xcargs[XCfunctional::enum_chi_bb] = chib;
            xcargs[XCfunctional::enum_chi_ab] = chiab;
        }
    }

    world.gop.fence();
    truncate(world, xc_args, extra_truncation);
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
    return impl->operator()(vket);
};

template<typename T, std::size_t NDIM>
Exchange<T,NDIM>& Exchange<T,NDIM>::set_bra_and_ket(const vecfuncT& bra, const vecfuncT& ket) {
    MADNESS_CHECK(impl);
//    if (not impl) {
//        World& world=bra.front().world();
//        impl.reset(new Exchange<T,NDIM>::ExchangeImpl(world,lo1,FunctionDefaults<NDIM>::get_thresh()));
//    }
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
Exchange<T,NDIM>& Exchange<T,NDIM>::set_algorithm(const Algorithm& alg) {
    impl->set_algorithm(alg);
    return *this;
}

template<typename T, std::size_t NDIM>
Exchange<T,NDIM>& Exchange<T,NDIM>::set_printlevel(const long& level) {
    impl->set_printlevel(level);
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


