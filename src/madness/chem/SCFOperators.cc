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
#include<madness/mra/funcplot.h>


using namespace madness;

namespace madness {

namespace {

/// TEMPORARY probe: pointwise min/max of a Function

/// Function carries no min/max, so accumulate over the leaf values. The op is
/// copied into every task, hence the shared, mutex-guarded accumulator.
struct fun_minmax {
    typedef double resultT;
    std::shared_ptr<std::pair<double, double> > acc;
    std::shared_ptr<Mutex> mtx;
    fun_minmax() : acc(new std::pair<double, double>(1.e300, -1.e300)),
                   mtx(new Mutex()) {}
    Tensor<double> operator()(const Key<3>& key, const Tensor<double>& t) const {
        double lo = 1.e300, hi = -1.e300;
        const double* p = t.ptr();
        for (long i = 0; i < t.size(); ++i) {
            lo = std::min(lo, p[i]);
            hi = std::max(hi, p[i]);
        }
        ScopedMutex<Mutex> lock(mtx.get());
        acc->first = std::min(acc->first, lo);
        acc->second = std::max(acc->second, hi);
        return copy(t);
    }
    template<typename Archive> void serialize(Archive& ar) {}
};

/// TEMPORARY probe: one line of census for a Function, gated on MAD_XC_PROBE
void xc_probe(const std::string& tag, const real_function_3d& f) {
    const int level = xc_probe_level();
    if (level < 1) return;
    // f.world() dereferences the impl, so it cannot be called before this check:
    // tau is unassigned for a plain GGA and probing it segfaulted every pbe run
    if (not f.is_initialized()) return;
    World& world = f.world();
    real_function_3d g = copy(f);
    g.reconstruct();
    const double nrm = g.norm2();
    const std::size_t tree = g.tree_size();
    const std::size_t depth = g.max_depth();
    const double gb = double(g.size()) * 8.e-9;
    double lo = 0.0, hi = 0.0;
    if (level >= 2) {
        fun_minmax mm;
        real_function_3d dummy = unary_op(g, mm);
        dummy.clear();
        world.gop.min(mm.acc->first);
        world.gop.max(mm.acc->second);
        lo = mm.acc->first;
        hi = mm.acc->second;
    }
    if (world.rank() == 0) {
        printf("XCPROBE %-24s  norm %12.5e  tree %8lu  depth %3lu  GB %8.4f",
               tag.c_str(), nrm, (unsigned long) tree, (unsigned long) depth, gb);
        if (level >= 2) printf("  min %12.5e  max %12.5e", lo, hi);
        printf("\n");
    }
}


}   // anonymous namespace

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

    nbeta = nemo->get_calc()->param.nbeta();
    const bool have_beta = xc->is_spin_polarized() && nbeta != 0;

    // compute the alpha and beta densities
    real_function_3d arho, brho;
    real_function_3d arhonemo = nemo->make_density(nemo->get_calc()->aocc, nemo->get_calc()->amo);
    arho = (arhonemo * nemo->R_square).truncate(extra_truncation);
    real_function_3d brhonemo;
    if (have_beta) {
        brhonemo = nemo->make_density(nemo->get_calc()->bocc, nemo->get_calc()->bmo);
        brho = (brhonemo * nemo->R_square).truncate(extra_truncation);
    } else {
        brho = arho;
        brhonemo = arhonemo;
    }

    // the regularized densities let prep_xc_args build zeta without putting the
    // nuclear cusp of rho = R^2 rho_reg under a numerical derivative
    xc_args = prep_xc_args(arho, brho, arhonemo, brhonemo);
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
    const int e_tau = (spin == 0) ? XCfunctional::enum_taua : XCfunctional::enum_taub;
    if (xc_args[e_tau].is_initialized()) return xc_args[e_tau];

    // Pointwise route: tau is deliberately not a Function -- it is assembled from
    // the smooth pieces and U1 values inside make_libxc_args, where nothing gets
    // projected. Reconstructing it here reintroduces exactly the projection the
    // route exists to avoid, so this is for diagnostics and tests only and its
    // near-nucleus values are not what the functional saw.
    const int e_grad = (spin == 0) ? XCfunctional::enum_gradfa : XCfunctional::enum_gradfb;
    const int e_n    = (spin == 0) ? XCfunctional::enum_na     : XCfunctional::enum_nb;
    const int e_gx   = (spin == 0) ? XCfunctional::enum_Ga_x   : XCfunctional::enum_Gb_x;
    if (not xc_args[e_grad].is_initialized()) return real_function_3d();
    MADNESS_CHECK_THROW(bool(ncf), "regularized tau pieces without a correlation factor");

    real_function_3d r = copy(xc_args[e_grad]);
    NuclearCorrelationFactor::U1_dot_U1_functor u1_dot_u1(ncf.get());
    const real_function_3d U1dot =
            real_factory_3d(world).functor(u1_dot_u1).truncate_on_project();
    r += U1dot * xc_args[e_n];
    const vecfuncT U1 = ncf->U1vec();
    for (int axis = 0; axis < 3; ++axis) r -= 2.0 * U1[axis] * xc_args[e_gx + axis];
    r = r * ncf->square();
    return (0.5 * r).truncate(extra_truncation);
}


/// compute tau = 1/2 sum_i |grad psi_i|^2 and store it in the intermediates
template<typename T, std::size_t NDIM>
void XCOperator<T, NDIM>::set_tau(const vecfuncT &amo, const Tensor<double> &aocc,
                                  const vecfuncT &bmo, const Tensor<double> &bocc,
                                  const TauU1 u1mode_in) const {

    MADNESS_CHECK_THROW(is_initialized(), "set_tau called before the intermediates exist");

    // resolve the environment default once, so a single run cannot mix routes
    const TauU1 u1mode = (u1mode_in != TauU1::from_env) ? u1mode_in
                         : (xc_tau_no_u1()        ? TauU1::none
                         : (xc_tau_u1_pointwise() ? TauU1::pointwise
                                                  : TauU1::mra));

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
            if (ncf and u1mode != TauU1::none) p.G[axis] = dot(world, mo_copy, dmo);
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


/// inject a kinetic energy density computed elsewhere
template<typename T, std::size_t NDIM>
void XCOperator<T, NDIM>::set_tau(const real_function_3d &taua,
                                  const real_function_3d &taub) const {

    MADNESS_CHECK_THROW(is_initialized(), "set_tau called before the intermediates exist");
    MADNESS_CHECK_THROW(taua.is_initialized(), "set_tau: alpha tau is not initialized");

    const bool have_beta = (xc->is_spin_polarized()) and (nbeta > 0);
    MADNESS_CHECK_THROW(not have_beta or taub.is_initialized(),
                        "set_tau: an open-shell meta-gga needs a beta tau too");

    xc_args[XCfunctional::enum_taua] = copy(taua);
    if (have_beta) xc_args[XCfunctional::enum_taub] = copy(taub);
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


/// divergence of a vector field, honouring dft_deriv
/// see the declaration in SCFOperators.h for why this is shared rather than local
real_function_3d div_dft_deriv(World& world, const std::vector<real_function_3d>& v,
                               const std::string& dft_deriv, bool do_refine) {
    MADNESS_CHECK(v.size() == 3);
    // the operand is a real vector field regardless of any operator's value type,
    // so this cannot go through make_derivative(), which is Derivative<T,NDIM>
    // same preparation as vmra.h's div(v, do_refine=true), which this replaces:
    // refining before differentiating makes the result more accurate
    reconstruct(world, v);
    if (do_refine) refine(world, v);
    // the operators must outlive the un-fenced applies below, so build them all
    // first -- a Derivative destroyed while its deferred tasks are still queued is
    // a use-after-free
    std::vector<std::shared_ptr<Derivative<double, 3> > > D(3);
    for (int axis = 0; axis < 3; ++axis) {
        D[axis].reset(new Derivative<double, 3>(world, axis));
        if (dft_deriv == "bspline") D[axis]->set_bspline1();
        else if (dft_deriv == "ble") D[axis]->set_ble1();
    }
    std::vector<real_function_3d> dv(3);
    for (int axis = 0; axis < 3; ++axis) dv[axis] = apply(*D[axis], v[axis], false);
    world.gop.fence();
    // sum() compresses first -- apply(Derivative,...) returns a reconstructed
    // function, and gaxpy needs both operands compressed
    return sum(world, dv, true);
}


template<typename T, std::size_t NDIM>
real_function_3d XCOperator<T, NDIM>::div_dft_deriv(const vecfuncT& v) const {
    return madness::div_dft_deriv(world, v, dft_deriv);
}


template<typename T, std::size_t NDIM>
bool XCOperator<T, NDIM>::is_weak_form() const {
    return weak_form_ok and xc_weak_gga() and xc->needs_sigma();
}


/// weak-form split of the non-multiplicative xc terms
template<typename T, std::size_t NDIM>
void XCOperator<T, NDIM>::weak_xc_terms(const std::vector<Function<T,NDIM> >& vket,
                                        std::vector<Function<T,NDIM> >& mult,
                                        std::vector<std::vector<Function<T,NDIM> > >& flux) const {

    MADNESS_CHECK_THROW(is_weak_form(), "weak_xc_terms outside weak form");
    MADNESS_CHECK_THROW(semilocal_flux.size() == 3,
                        "weak_xc_terms before make_xc_potential");

    const double vtol = FunctionDefaults<3>::get_thresh() * 0.1;
    const std::size_t n = vket.size();

    vecfuncT U1;
    if (ncf) U1 = ncf->U1vec();

    mult = zero_functions_compressed<T,NDIM>(world, n);

    // accumulate the flux one Cartesian component at a time, as vectors over the
    // orbitals: the vector gaxpy coerces tree states, the scalar Function::gaxpy
    // throws on a mismatch instead
    std::vector<std::vector<Function<T,NDIM> > > Y(3);
    for (int axis = 0; axis < 3; ++axis)
        Y[axis] = zero_functions_compressed<T,NDIM>(world, n);

    // ---- semilocal part: Y_i += X F_i,  mult_i += X.grad(F_i)
    //
    // grad(F_i) is the gradient of the *nemo*, which is cusp-free, and X is only
    // ever multiplied. Both products inherit X's jump at the nucleus, which is
    // harmless: the first is convolved with G, the second enters V psi and is
    // convolved with G as well. Nothing differentiates X.
    for (int axis = 0; axis < 3; ++axis) {
        auto D = make_derivative(axis);
        std::vector<Function<T,NDIM> > dket = apply(world, *D, vket, false);
        world.gop.fence();
        mult += mul_sparse(world, semilocal_flux[axis], dket, vtol);
        Y[axis] += mul_sparse(world, semilocal_flux[axis], vket, vtol);
    }

    // ---- meta-gga part: W_i = v_tau (grad F_i - U1 F_i),
    //      Y_i += 1/2 W_i,   mult_i += 1/2 U1.W_i
    //
    // Same R cancellation as apply_tau_term: with psi = R F,
    //   R^{-1}[-1/2 div(v_tau grad(R F))] = -1/2 (div W - U1.W),
    // so the divergence piece is -1/2 div W (through G) and the rest is
    // +1/2 U1.W (multiplicative). v_tau is never differentiated either way --
    // and unlike apply_tau_term, neither is the product v_tau grad(psi).
    if (has_tau_term()) {
        MADNESS_CHECK_THROW(vtau.is_initialized(),
                            "weak_xc_terms before make_xc_potential");
        for (int axis = 0; axis < 3; ++axis) {
            auto D = make_derivative(axis);
            std::vector<Function<T,NDIM> > W = apply(world, *D, vket, false);
            world.gop.fence();
            if (ncf) W = sub(world, W, mul(world, U1[axis], vket));
            W = mul_sparse(world, vtau, W, vtol);
            Y[axis] += 0.5 * W;
            if (ncf) mult += 0.5 * mul(world, U1[axis], W);
        }
    }

    truncate(world, mult);
    for (int axis = 0; axis < 3; ++axis) truncate(world, Y[axis]);

    // transpose into the per-orbital vector fields the caller wants
    flux.assign(n, std::vector<Function<T,NDIM> >(3));
    for (std::size_t i = 0; i < n; ++i)
        for (int axis = 0; axis < 3; ++axis) flux[i][axis] = Y[axis][i];
}


/// xc contribution to the Fock matrix in weak form
template<typename T, std::size_t NDIM>
Tensor<T> XCOperator<T, NDIM>::weak_xc_matrix(const std::vector<Function<T,NDIM> >& vket,
                                              const real_function_3d& v_local) const {

    MADNESS_CHECK_THROW(is_weak_form(), "weak_xc_matrix outside weak form");
    MADNESS_CHECK_THROW(semilocal_flux.size() == 3,
                        "weak_xc_matrix before make_xc_potential");

    const double vtol = FunctionDefaults<3>::get_thresh() * 0.1;
    const std::size_t n = vket.size();

    // the bra carries R^2 when there is a nuclear correlation factor, so that the
    // integrals below are over the physical orbitals psi = R F
    real_function_3d R2;
    if (ncf) R2 = ncf->square();
    auto weight = [&](const real_function_3d& f) {
        return ncf ? (R2 * f).truncate() : f;
    };

    vecfuncT U1;
    if (ncf) U1 = ncf->U1vec();

    const long nn = long(n);
    Tensor<T> result(nn, nn);

    // local term: <psi_i| de/drho |psi_j>
    {
        std::vector<Function<T,NDIM> > bra =
                mul_sparse(world, weight(v_local), vket, vtol);
        result += matrix_inner(world, bra, vket);
    }

    // semilocal term: int X.grad(psi_i psi_j)
    //   = sum_a int R^2 X_a [ F_j dF_i/da + F_i dF_j/da ] - 2 int R^2 (U1.X) F_i F_j
    // The first bracket is A + A^T with A_ij = int (R^2 X_a dF_i/da) F_j, so the
    // result is symmetric by construction -- which the divergence form is not.
    for (int axis = 0; axis < 3; ++axis) {
        auto D = make_derivative(axis);
        std::vector<Function<T,NDIM> > dket = apply(world, *D, vket, false);
        world.gop.fence();
        std::vector<Function<T,NDIM> > bra =
                mul_sparse(world, weight(semilocal_flux[axis]), dket, vtol);
        Tensor<T> A = matrix_inner(world, bra, vket);
        result += A;
        result += transpose(A);
    }
    if (ncf) {
        real_function_3d u1x = dot(world, U1, semilocal_flux).truncate();
        std::vector<Function<T,NDIM> > bra =
                mul_sparse(world, weight(u1x), vket, vtol);
        result -= 2.0 * matrix_inner(world, bra, vket);
    }

    // meta-gga term: 1/2 int v_tau grad(psi_i).grad(psi_j)
    //   = 1/2 sum_a int R^2 v_tau (dF_i/da - U1_a F_i)(dF_j/da - U1_a F_j)
    if (has_tau_term()) {
        MADNESS_CHECK_THROW(vtau.is_initialized(),
                            "weak_xc_matrix before make_xc_potential");
        for (int axis = 0; axis < 3; ++axis) {
            auto D = make_derivative(axis);
            std::vector<Function<T,NDIM> > g = apply(world, *D, vket, false);
            world.gop.fence();
            if (ncf) g = sub(world, g, mul(world, U1[axis], vket));
            truncate(world, g);
            std::vector<Function<T,NDIM> > bra =
                    mul_sparse(world, weight(vtau), g, vtol);
            result += 0.5 * matrix_inner(world, bra, g);
        }
    }

    return result;
}


/// true once set_tau() has supplied tau, by either route
template<typename T, std::size_t NDIM>
bool XCOperator<T, NDIM>::has_tau_args() const {
    if (xc_args[XCfunctional::enum_taua].is_initialized()) return true;      // moldft / mra route
    return xc_args[XCfunctional::enum_gradfa].is_initialized()               // pointwise route
       and xc_args[XCfunctional::enum_nemo_R2].is_initialized();
}


/// the four analytic U1 quantities the xc ops evaluate pointwise

/// Empty unless the pointwise route is in use. They are handed to the op rather
/// than projected into xc_args precisely because a projected product with U1 is
/// what rings; see nemo_u1_functors.
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


template<typename T, std::size_t NDIM>
real_function_3d XCOperator<T, NDIM>::make_xc_potential() const {

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

    // 2 de/dsigma, diagnostic, from the same pointwise pass. Index from the shared
    // helper so producer and consumer cannot disagree about the layout.
    {
        const int islot = xc_dfds_slot(xc->is_spin_polarized(), xc->needs_tau(),
                                       xc->needs_sigma());
        if (islot >= 0) {
            MADNESS_CHECK_THROW(islot < int(intermediates.size()),
                                "MAD_XC_EXPORT_DFDS slot past the end of vxc's result");
            dfdsigma = intermediates[islot];
        }
    }

    xc_probe("rho_a", xc_args[XCfunctional::enum_rhoa]);
    xc_probe("tau_a", xc_args[XCfunctional::enum_taua]);
    xc_probe("local(de/drho)", dft_pot);
    if (has_tau_term()) xc_probe("vtau", vtau);

    if (xc->needs_sigma()) {

        const bool factored = xc_factored_gga_potential();
        const bool have_beta = xc->is_spin_polarized() && nbeta != 0;

        // the density that the factored form divides out of the flux: the total
        // density when spin-restricted (that is what vxc used), the density of the
        // spin this operator acts on for the same-spin term
        real_function_3d rho_same, rho_other;
        if (factored) {
            if (xc->is_spin_polarized()) {
                rho_same = xc_args[ispin == 0 ? XCfunctional::enum_rhoa
                                              : XCfunctional::enum_rhob];
                rho_other = xc_args[ispin == 0 ? XCfunctional::enum_rhob
                                               : XCfunctional::enum_rhoa];
            } else {
                rho_same = 2.0 * xc_args[XCfunctional::enum_rhoa];
            }
        }

        vecfuncT semilocal(3);
        semilocal[0] = intermediates[1];
        semilocal[1] = intermediates[2];
        semilocal[2] = intermediates[3];
        xc_probe("flux_same[x]", semilocal[0]);

        // MAD_XC_PLOT: dump the semilocal ingredients on ONE grid, as columns of
        // a single plot_plane file, so they are sampled identically and can be
        // compared point by point. Column order is documented here and nowhere
        // else, so keep the two in step:
        //   3 rho   4-6 grad(rho)_{x,y,z}   7-9 zeta_{x,y,z}
        //   10-12 X_{x,y,z}   13 df/drho
        // (columns 1-2 are the plane coordinates). rho and grad(rho) are the
        // TOTAL density: closed shell, so rho = 2 rho_alpha and grad rho = rho zeta.
        // All of these are post-refine_to_common_level, i.e. exactly what the
        // functional and the divergence actually see.
        if (not xc_plot_tag().empty()) {
            const real_function_3d& ra = xc_args[XCfunctional::enum_rhoa];
            std::vector<real_function_3d> cols;
            auto push = [&](const real_function_3d& f) {
                real_function_3d g = copy(f);
                g.reconstruct();
                cols.push_back(g);
            };
            const real_function_3d rho_tot = (2.0 * ra).truncate();
            push(rho_tot);
            for (int ax = 0; ax < 3; ++ax)   // grad(rho) = rho * zeta
                push((rho_tot * xc_args[XCfunctional::enum_zetaa_x + ax]).truncate());
            for (int ax = 0; ax < 3; ++ax)
                push(xc_args[XCfunctional::enum_zetaa_x + ax]);
            for (int ax = 0; ax < 3; ++ax)
                push(semilocal[ax]);
            push(intermediates[0]);          // df/drho, before any divergence
            plot_plane<3>(world, cols, "diag_" + xc_plot_tag());
        }

        if (is_weak_form()) {
            // keep the flux, take no divergence at all. The same-spin and
            // cross-spin contributions enter the potential as one divergence, so
            // they are summed here and the caller never needs them apart.
            MADNESS_CHECK_THROW(not xc_factored_gga_potential(),
                                "MAD_XC_WEAK_GGA and MAD_XC_FACTORED_GGA are "
                                "mutually exclusive: the weak form needs the flux "
                                "X, the factored form stores X/rho");
            semilocal_flux = copy(world, semilocal);
            if (have_beta) {
                semilocal_flux[0] += intermediates[4];
                semilocal_flux[1] += intermediates[5];
                semilocal_flux[2] += intermediates[6];
            }
            truncate(world, semilocal_flux);
            xc_probe("weak flux[x]", semilocal_flux[0]);
            truncate(world, xc_args);
            xc_probe("xc_pot (local only)", dft_pot);
            return dft_pot.truncate();
        }

        // ---- algorithm A'': div(X) = grad(u).grad(rho) + u laplacian(rho),
        //      with laplacian(rho) = div(rho zeta), i.e. one divergence of the
        //      semi-analytic gradient. See XC-POTENTIAL-ALGORITHMS.md 1G.3.1.
        //      The point of this branch is to measure what the expansion costs:
        //      div_dft_deriv() satisfies int div(X) = 0 structurally, a sum of two
        //      separately truncated terms does not.
        if (xc_onelevel_div()) {
            MADNESS_CHECK_THROW(not xc->is_spin_polarized(),
                                "MAD_XC_ONELEVEL_DIV is spin-restricted only");
            MADNESS_CHECK_THROW(not xc_factored_gga_potential(),
                                "MAD_XC_ONELEVEL_DIV needs the unfactored flux");

            const std::size_t nr = intermediates.size();
            const real_function_3d u = intermediates[nr - 4];

            // grad(rho) as the flux itself uses it, i.e. the *munged* rho*zeta.
            // Rebuilding it here from the MRA rho and zeta gave a field that does
            // not vanish at the box wall (zeta is garbage where rho is at noise
            // level), so div() picked up a surface term and int laplacian(rho)
            // came out at -17 instead of 0.
            vecfuncT drho(3);
            const bool use_ddens = bool(ncf) and arho_reg.is_initialized()
                                   and not xc_a2_munged_grad();
            if (use_ddens) {
                // make_ddensity's route: grad(rho) = R^2 (grad n - 2 U1 n), with n
                // the total regularized density. No zeta anywhere, so nothing
                // carries tail garbage to the box wall.
                const real_function_3d n_tot = 2.0 * arho_reg;
                const real_function_3d R2 = ncf->square();
                const vecfuncT U1 = ncf->U1vec();
                real_function_3d n_ref = copy(n_tot).refine();
                for (int axis = 0; axis < 3; ++axis) {
                    real_derivative_3d D = free_space_derivative<double, 3>(world, axis);
                    drho[axis] = (R2 * (D(n_ref) - 2.0 * U1[axis] * n_tot)).truncate();
                }
            } else {
                for (int axis = 0; axis < 3; ++axis)
                    drho[axis] = intermediates[nr - 3 + axis];
            }
            truncate(world, drho);

            const real_function_3d lap = div_dft_deriv(drho);          // conserving
            vecfuncT gu = grad(copy(u).refine(), true);                // one derivative of u
            truncate(world, gu);

            const real_function_3d term1 = dot(world, gu, drho).truncate();
            const real_function_3d term2 = (u * lap).truncate();
            const real_function_3d divX = (term1 + term2).truncate();

            // the reference: the same quantity as a genuine divergence
            const real_function_3d divX_ref = div_dft_deriv(semilocal);

            if (xc_probe_level() >= 1) {
                // Every integral here goes through a reconstructed copy: trace()
                // on a function fresh out of truncate()/sum()/div()/an operator
                // apply or refine_to_common_level() returns a silently wrong
                // value. The N_elec line below is the validation -- if it is not
                // the electron count, no other number in this block is usable.
                auto itrace = [&](const real_function_3d& f) {
                    if (not f.is_initialized()) return 0.0;
                    real_function_3d g = copy(f);
                    g.reconstruct();
                    return g.trace();
                };

                const real_function_3d& rho_a = xc_args[XCfunctional::enum_rhoa];
                const double i_rho_a = itrace(rho_a);
                const double i_reg   = itrace(arho_reg);
                real_function_3d rho_rebuilt;
                if (use_ddens) rho_rebuilt = (ncf->square() * (2.0 * arho_reg)).truncate();

                // the difference relative to the POTENTIAL, not just to the
                // divergence term, plus where it sits pointwise
                real_function_3d dp = copy(dft_pot); dp.compress();
                real_function_3d dref = copy(divX_ref); dref.compress();
                real_function_3d da2 = copy(divX); da2.compress();
                const real_function_3d vxc_ref = (dp - dref).truncate();
                const real_function_3d diff = (da2 - dref).truncate();
                const double nv = vxc_ref.norm2();
                xc_probe("A2: divX - reference", diff);
                xc_probe("A2: v_xc (reference)", vxc_ref);

                const double t1 = itrace(term1), t2 = itrace(term2);
                const double tl = itrace(lap);
                const double te = itrace(divX), tr = itrace(divX_ref);
                const double dn = (divX - divX_ref).norm2(), rn = divX_ref.norm2();

                if (world.rank() == 0) {
                    print("A2: grad(rho) route =", use_ddens ? "make_ddensity" : "rho*zeta");
                    printf("A2CHK int 2*rho_alpha        %14.8f   <== must equal N_elec\n",
                           2.0 * i_rho_a);
                    if (use_ddens)
                        printf("A2CHK int R^2*2*arho_reg     %14.8f   <== must equal N_elec\n",
                               itrace(rho_rebuilt));
                    printf("A2CHK int arho_reg           %14.8f   (= N_elec/2 / <R^2>)\n", i_reg);
                    printf("A2TRACE int laplacian(rho)   %14.6e   (exact 0)\n", tl);
                    printf("A2TRACE int grad(u).grad(rho)%14.6e\n", t1);
                    printf("A2TRACE int u laplacian(rho) %14.6e\n", t2);
                    printf("A2TRACE int div(X) expanded  %14.6e   <-- the cost of A''\n", te);
                    printf("A2TRACE int div(X) reference %14.6e   (structural zero)\n", tr);
                    printf("A2TRACE ||expanded-reference|| %12.4e   rel(divX) %10.4e"
                           "   rel(v_xc) %10.4e\n", dn, dn / rn, dn / nv);
                    printf("A2TRACE grad(rho) route = %s\n",
                           use_ddens ? "make_ddensity (unmunged)" : "munged rho*zeta (== X)");
                }
            }

            dft_pot -= divX;
            truncate(world, xc_args);
            xc_probe("xc_pot", dft_pot);
            if (not xc_plot_tag().empty()) {
                real_function_3d vplot = copy(dft_pot);
                vplot.reconstruct();
                plot_plane<3>(world, std::vector<real_function_3d>{vplot},
                              "vxc_" + xc_plot_tag());
            }
            return dft_pot.truncate();
        }

        // drop the noise the common refinement level added, before the derivative
        // can amplify it by 2^n -- see xc_flux_truncation()
        const double ftrunc = xc_flux_truncation();
        if (ftrunc > 0.0) {
            truncate(world, semilocal, ftrunc * FunctionDefaults<3>::get_thresh());
            xc_probe("flux_same[x] trunc", semilocal[0]);
        }

        // second term in Yanai2005, Eq. (12)
        real_function_3d gga_pot_same_spin = div_dft_deriv(semilocal);
        xc_probe("div(flux_same)", gga_pot_same_spin);
        if (factored) {
            // div(rho B) = rho div B + rho (zeta.B); the second piece was folded
            // into the local term by vxc, so only the multiplication is left. rho
            // damps the noise of the divergence exactly where the term is small,
            // which is what a bare div() of the full flux cannot do.
            gga_pot_same_spin = (rho_same * gga_pot_same_spin).truncate();
            xc_probe("rho*div(flux_same)", gga_pot_same_spin);
        }
        dft_pot -= gga_pot_same_spin;

        if (have_beta) {
            semilocal[0] = intermediates[4];
            semilocal[1] = intermediates[5];
            semilocal[2] = intermediates[6];

            // third term in Yanai2005, Eq. (12)
            if (ftrunc > 0.0)
                truncate(world, semilocal, ftrunc * FunctionDefaults<3>::get_thresh());
            real_function_3d gga_pot_other_spin = div_dft_deriv(semilocal);
            if (factored) {
                gga_pot_other_spin = (rho_other * gga_pot_other_spin).truncate();
            }
            dft_pot -= gga_pot_other_spin;
        }
    }

    truncate(world, xc_args);
    xc_probe("xc_pot", dft_pot);
    if (not xc_plot_tag().empty()) {
        // plot_plane evaluates pointwise, which requires a reconstructed function
        real_function_3d vplot = copy(dft_pot);
        vplot.reconstruct();
        plot_plane<3>(world, std::vector<real_function_3d>{vplot},
                      "vxc_" + xc_plot_tag());
    }
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
                                           const real_function_3d &brho,
                                           const real_function_3d &arho_reg_in,
                                           const real_function_3d &brho_reg) const {

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
            else if (dft_deriv == "ble") return grad_ble_one(f);     // BLE
            else return grad(f);                                    // Default is abgv
        };

        // zeta = grad log(rho). With a nuclear correlation factor the density is
        // rho = R^2 rho_reg, so
        //
        //   zeta = grad log(R^2) + grad log(rho_reg) = -2 U1 + grad log(rho_reg)
        //
        // (U1 = -grad(R)/R). Differentiating log(rho) directly puts the nuclear
        // cusp under the derivative operator, which is what the regularization
        // exists to avoid: the kink forces refinement to the finest level and the
        // O(thresh) noise it leaves there is amplified by 2^n by the derivative.
        // Only the cusp-free rho_reg is differentiated here, exactly as in set_tau.
        arho_reg = arho_reg_in;          // for algorithm A'' (may be unassigned)
        const bool regularized_zeta = bool(ncf) and xc_nemo_zeta();
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
            vecfuncT gradb = make_zeta(brho, brho_reg);
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


