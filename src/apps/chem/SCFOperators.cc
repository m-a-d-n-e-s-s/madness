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


#include <apps/chem/SCFOperators.h>
#include <apps/chem/SCF.h>
#include <apps/chem/nemo.h>
#include <apps/chem/correlationfactor.h>
#include <apps/chem/xcfunctional.h>

using namespace madness;

namespace madness {


/// munge the lhs quantity under the condition that the rhs quantity > tol
struct binary_munging{

    binary_munging() : TOL(1.e-7) {}
    binary_munging(const double t) : TOL(t) {}
    double TOL;

    void operator()(const Key<3>& key, Tensor<double>& U, const Tensor<double>& t,
            const Tensor<double>& rho) const {
        ITERATOR(U,
            double d = t(IND);
            double p = rho(IND);
            U(IND) = p > TOL ? d : 0.0;
        );
    }
    template <typename Archive>
    void serialize(Archive& ar) {}

};

template<typename T, std::size_t NDIM>
DistributedMatrix<T> Kinetic<T,NDIM>::kinetic_energy_matrix(World & world,
        const vecfuncT & v) const {
    int n = v.size();
    DistributedMatrix<T> r = column_distributed_matrix<T>(world, n, n);
    reconstruct(world, v);

    // apply the derivative operator on each function for each dimension
    std::vector<vecfuncT> dv(NDIM);
    for (std::size_t i=0; i<NDIM; ++i) {
        dv[i]=apply(world, *(gradop[i]), v, false);
    }
    world.gop.fence();
    for (std::size_t i=0; i<NDIM; ++i) {
        compress(world,dv[i],false);
    }
    world.gop.fence();
    for (std::size_t i=0; i<NDIM; ++i) {
        r += matrix_inner(r.distribution(), dv[i], dv[i], true);
    }
    r *= 0.5;
    return r;
}


template<typename T, std::size_t NDIM>
DistributedMatrix<T> Kinetic<T,NDIM>::kinetic_energy_matrix(World & world,
        const vecfuncT & vbra, const vecfuncT & vket) const {
    int n = vbra.size();
    int m = vket.size();
    DistributedMatrix<T> r = column_distributed_matrix<T>(world, n, m);
    reconstruct(world, vbra);
    reconstruct(world, vket);

    // apply the derivative operator on each function for each dimension
    std::vector<vecfuncT> dvbra(NDIM),dvket(NDIM);
    for (std::size_t i=0; i<NDIM; ++i) {
        dvbra[i]=apply(world, *(gradop[i]), vbra, false);
        dvket[i]=apply(world, *(gradop[i]), vket, false);
    }
    world.gop.fence();
    for (std::size_t i=0; i<NDIM; ++i) {
        compress(world,dvbra[i],false);
        compress(world,dvket[i],false);
    }
    world.gop.fence();
    for (std::size_t i=0; i<NDIM; ++i) {
        r += matrix_inner(r.distribution(), dvbra[i], dvket[i], false);
    }
    r *= 0.5;
    return r;
}

// explicit instantiation
template class Kinetic<double,1>;
template class Kinetic<double,2>;
template class Kinetic<double,3>;
template class Kinetic<double,4>;
template class Kinetic<double,5>;
template class Kinetic<double,6>;


template<typename T, std::size_t NDIM>
std::vector<Function<T,NDIM> >
Laplacian<T,NDIM>::operator()(const std::vector<Function<T,NDIM> >& vket) const {

    refine(world,vket);     // for better accuracy
    vecfuncT result=zero_functions_compressed<T,NDIM>(world,vket.size());
    SeparatedConvolution<T,NDIM> smooth=SmoothingOperator<NDIM>(world,eps);


    for (int idim=0; idim<NDIM; ++idim) {
        vecfuncT dvket=apply(world,*gradop[idim].get(),vket);
        refine(world,dvket);
        if (eps>0.0) dvket=apply(world,smooth,dvket);
        vecfuncT ddvket=apply(world,*gradop[idim].get(),dvket);
        result=add(world,result,ddvket);
    }

    if (eps>0.0) result=apply(world,smooth,result);

    return result;
}

// explicit instantiation
template class Laplacian<double,1>;
template class Laplacian<double,2>;
template class Laplacian<double,3>;
template class Laplacian<double,4>;
template class Laplacian<double,5>;
template class Laplacian<double,6>;


Coulomb::Coulomb(World& world, const Nemo* nemo) : world(world),
        R_square(nemo->R_square) {
    vcoul=compute_potential(nemo);
}


real_function_3d Coulomb::compute_density(const SCF* calc) const {
    real_function_3d density = calc->make_density(world, calc->get_aocc(),
            calc->get_amo());
    if (calc->is_spin_restricted()) {
        density.scale(2.0);
    } else {
        real_function_3d brho = calc->make_density(world, calc->get_bocc(),
                calc->get_bmo());
        density+=brho;
    }
    density.truncate();
    return density;
}

real_function_3d Coulomb::compute_potential(const madness::SCF* calc) const {
    real_function_3d density=compute_density(calc);
    return calc->make_coulomb_potential(density);
}

/// same as above, but with the additional factor R^2 in the density
real_function_3d Coulomb::compute_potential(const madness::Nemo* nemo) const {
    real_function_3d density=nemo->make_density(nemo->get_calc()->aocc,
            nemo->get_calc()->amo);
    if (nemo->get_calc()->is_spin_restricted()) {
        density.scale(2.0);
    } else {
        real_function_3d brho = nemo->get_calc()->make_density(world,
                nemo->get_calc()->get_bocc(),nemo->get_calc()->get_bmo());
        density+=brho;
    }
    density=density*R_square;
    density.truncate();
    return nemo->get_calc()->make_coulomb_potential(density);
}

real_function_3d Coulomb::compute_potential(const real_function_3d& density,
        double lo, double econv) const {
    real_convolution_3d poisson = CoulombOperator(world, lo, econv);
    return poisson(density).truncate();
}


Nuclear::Nuclear(World& world, const SCF* calc) : world(world) {
    ncf=std::shared_ptr<NuclearCorrelationFactor>(
            new PseudoNuclearCorrelationFactor(world,
            calc->molecule,calc->potentialmanager,1.0));
}

Nuclear::Nuclear(World& world, const Nemo* nemo) : world(world) {
    ncf=nemo->nuclear_correlation;
}

vecfuncT Nuclear::operator()(const vecfuncT& vket) const {

    const static std::size_t NDIM=3;

    // shortcut for local nuclear potential (i.e. no correlation factor)
    if (ncf->type()==NuclearCorrelationFactor::None) {
        vecfuncT result=mul(world,ncf->U2(),vket);
        truncate(world,result);
        return result;
    }

    std::vector< std::shared_ptr<Derivative<double,NDIM> > > gradop =
            gradient_operator<double,NDIM>(world);
    reconstruct(world, vket);
    vecfuncT vresult=zero_functions_compressed<double,NDIM>(world,vket.size());

    // memory-saving algorithm: outer loop over the dimensions
    // apply the derivative operator on each function for each dimension
    for (std::size_t i=0; i<NDIM; ++i) {
        std::vector<Function<double,NDIM> > dv=apply(world, *(gradop[i]), vket, true);
        truncate(world,dv);
        real_function_3d U1=ncf->U1(i%3);
        std::vector<Function<double,NDIM> > U1dv=mul(world,U1,dv);
        truncate(world,U1dv);
        vresult=add(world,vresult,U1dv);
    }

    real_function_3d U2=ncf->U2();
    std::vector<Function<double,NDIM> > U2v=mul(world,U2,vket);
    vresult=add(world,vresult,U2v);
    truncate(world,vresult);
    return vresult;
}


DNuclear::DNuclear(World& world, const SCF* calc, const int iatom, const int iaxis)
    : world(world), iatom(iatom), iaxis(iaxis) {
    ncf=std::shared_ptr<NuclearCorrelationFactor>(
            new PseudoNuclearCorrelationFactor(world,
            calc->molecule,calc->potentialmanager,1.0));
}

DNuclear::DNuclear(World& world, const Nemo* nemo, const int iatom, const int iaxis)
           : world(world), iatom(iatom), iaxis(iaxis) {
    ncf=nemo->nuclear_correlation;
}

vecfuncT DNuclear::operator()(const vecfuncT& vket) const {

    const static std::size_t NDIM=3;
    const double vthresh=FunctionDefaults<NDIM>::get_thresh()*0.1;

    // compute the U2 potential/ the derivative nuclear potential
    NuclearCorrelationFactor::U2X_functor u2x(ncf.get(),iatom,iaxis);
    real_function_3d u2x_f=real_factory_3d(world).functor(u2x)
            .thresh(vthresh).truncate_on_project();
    vecfuncT result=mul(world,u2x_f,vket);
    truncate(world,result,vthresh);

    // add U1 and U3 potentials if the nuclear correlation factor exists
    if (ncf->type() != NuclearCorrelationFactor::None) {

        std::vector< std::shared_ptr<Derivative<double,NDIM> > > gradop =
                gradient_operator<double,NDIM>(world);
        reconstruct(world, vket);

        // memory-saving algorithm: outer loop over the dimensions
        // apply the derivative operator on each function for each dimension
        for (std::size_t i=0; i<NDIM; ++i) {
            std::vector<Function<double,NDIM> > dv=apply(world, *(gradop[i]), vket, true);
            truncate(world,dv);

            // note the two different axis: U1axis (i) and the derivative axis (iaxis)
            // \frac{\partial U1_i}{\partial R_{A,iaxis}}
            // e.g. d/dYA U1x
            NuclearCorrelationFactor::U1X_functor u1x(ncf.get(),iatom,i,iaxis);
            real_function_3d U1=real_factory_3d(world).functor(u1x).truncate_on_project();
            std::vector<Function<double,NDIM> > U1dv=mul(world,U1,dv);
            truncate(world,U1dv);
            result=sub(world,result,U1dv);
            truncate(world,result);
        }

        // add the U3X potential
        NuclearCorrelationFactor::U3X_functor u3x(ncf.get(),iatom,iaxis);
        real_function_3d u3x_f=real_factory_3d(world).functor(u3x).truncate_on_project();
        std::vector<Function<double,NDIM> > U3v=mul(world,u3x_f,vket);
        result=sub(world,result,U3v);
        truncate(world,result);
    }
    truncate(world,result);

    return result;
}


Exchange::Exchange(World& world, const SCF* calc, const int ispin)
        : world(world), small_memory_(true), same_(false) {
    if (ispin==0) { // alpha spin
        mo_ket=calc->amo;
        occ=calc->aocc;
    } else if (ispin==1) {  // beta spin
        mo_ket=calc->bmo;
        occ=calc->bocc;
    }
    mo_bra=mo_ket;
    poisson = std::shared_ptr<real_convolution_3d>(
            CoulombOperatorPtr(world, calc->param.lo, calc->param.econv));
}

Exchange::Exchange(World& world, const Nemo* nemo, const int ispin)
    : world(world), small_memory_(true), same_(false) {

    if (ispin==0) { // alpha spin
        mo_ket=nemo->get_calc()->amo;
        occ=nemo->get_calc()->aocc;
    } else if (ispin==1) {  // beta spin
        mo_ket=nemo->get_calc()->bmo;
        occ=nemo->get_calc()->bocc;
    }

    mo_bra=mul(world,nemo->nuclear_correlation->square(),mo_ket);
    truncate(world,mo_bra);
    poisson = std::shared_ptr<real_convolution_3d>(
            CoulombOperatorPtr(world, nemo->get_calc()->param.lo,
                    nemo->get_calc()->param.econv));

}

void Exchange::set_parameters(const vecfuncT& bra, const vecfuncT& ket,
        const Tensor<double>& occ1, const double lo, const double econv) {
    mo_bra=copy(world,bra);
    mo_ket=copy(world,ket);
    occ=copy(occ1);
    poisson = std::shared_ptr<real_convolution_3d>(
            CoulombOperatorPtr(world, lo, econv));
}

vecfuncT Exchange::operator()(const vecfuncT& vket) const {
    const bool same = this->same();
    int nocc = mo_bra.size();
    int nf = vket.size();
    double tol = FunctionDefaults < 3 > ::get_thresh(); /// Important this is consistent with Coulomb
    vecfuncT Kf = zero_functions_compressed<double, 3>(world, nf);
    reconstruct(world, mo_bra);
    norm_tree(world, mo_bra);
    reconstruct(world, mo_ket);
    norm_tree(world, mo_ket);
    if (!same) {
        reconstruct(world, vket);
        norm_tree(world, vket);
    }

    if (small_memory_) {     // Smaller memory algorithm ... possible 2x saving using i-j sym
        for(int i=0; i<nocc; ++i){
            if(occ[i] > 0.0){
                vecfuncT psif = mul_sparse(world, mo_bra[i], vket, tol); /// was vtol
                truncate(world, psif);
                psif = apply(world, *poisson.get(), psif);
                truncate(world, psif);
                psif = mul_sparse(world, mo_ket[i], psif, tol); /// was vtol
                gaxpy(world, 1.0, Kf, occ[i], psif);
            }
        }
    } else {    // Larger memory algorithm ... use i-j sym if psi==f
        vecfuncT psif;
        for (int i = 0; i < nocc; ++i) {
            int jtop = nf;
            if (same)
                jtop = i + 1;
            for (int j = 0; j < jtop; ++j) {
                psif.push_back(mul_sparse(mo_bra[i], vket[j], tol, false));
            }
        }

        world.gop.fence();
        truncate(world, psif);
        psif = apply(world, *poisson.get(), psif);
        truncate(world, psif, tol);
        reconstruct(world, psif);
        norm_tree(world, psif);
        vecfuncT psipsif = zero_functions<double, 3>(world, nf * nocc);
        int ij = 0;
        for (int i = 0; i < nocc; ++i) {
            int jtop = nf;
            if (same)
                jtop = i + 1;
            for (int j = 0; j < jtop; ++j, ++ij) {
                psipsif[i * nf + j] = mul_sparse(psif[ij], mo_ket[i], false);
                if (same && i != j) {
                    psipsif[j * nf + i] = mul_sparse(psif[ij], mo_ket[j], false);
                }
            }
        }
        world.gop.fence();
        psif.clear();
        world.gop.fence();
        compress(world, psipsif);
        for (int i = 0; i < nocc; ++i) {
            for (int j = 0; j < nf; ++j) {
                Kf[j].gaxpy(1.0, psipsif[i * nf + j], occ[i], false);
            }
        }
        world.gop.fence();
        psipsif.clear();
        world.gop.fence();
    }
    truncate(world, Kf, tol);
    return Kf;

}

/// custom ctor with information about the XC functional
XCOperator::XCOperator(World& world, std::string xc_data, const bool spin_polarized,
        const real_function_3d& arho, const real_function_3d& brho)
    : world(world), nbeta(0), ispin(0), extra_truncation(FunctionDefaults<3>::get_thresh()*0.01) {
    xc=std::shared_ptr<XCfunctional> (new XCfunctional());
    xc->initialize(xc_data, spin_polarized, world);
    const bool have_beta=xc->is_spin_polarized() && nbeta != 0;

    std::vector<real_function_3d> darho(3),dbrho(3);
    if (xc->is_gga()) {
        for (int iaxis=0; iaxis<3; ++iaxis) {
            Derivative<double,3> D = free_space_derivative<double,3>(world, iaxis);
            darho[iaxis]=D(arho);
            if (have_beta) dbrho[iaxis]=D(brho);
        }
    }

    xc_args=prep_xc_args(arho,brho,darho,dbrho);

}


XCOperator::XCOperator(World& world, const SCF* calc, int ispin) : world(world),
        ispin(ispin), extra_truncation(FunctionDefaults<3>::get_thresh()*0.01) {
    xc=std::shared_ptr<XCfunctional> (new XCfunctional());
    xc->initialize(calc->param.xc_data, !calc->param.spin_restricted, world);
    nbeta=calc->param.nbeta;
    const bool have_beta=xc->is_spin_polarized() && nbeta != 0;

    // compute the alpha and beta densities
    real_function_3d arho,brho;
    arho=calc->make_density(world,calc->aocc,calc->amo);
    if (have_beta) {
        brho=calc->make_density(world,calc->bocc,calc->bmo);
    } else {
        brho=arho;
    }
    std::vector<real_function_3d> darho,dbrho;
    if (xc->is_gga()) {
        darho=calc->nabla(arho);
        if (have_beta) dbrho=calc->nabla(brho);
    }

    xc_args=prep_xc_args(arho,brho,darho,dbrho);
}

XCOperator::XCOperator(World& world, const Nemo* nemo, int ispin) : world(world),
        ispin(ispin), extra_truncation(FunctionDefaults<3>::get_thresh()*0.01) {
    xc=std::shared_ptr<XCfunctional> (new XCfunctional());
    xc->initialize(nemo->get_calc()->param.xc_data,
            !nemo->get_calc()->param.spin_restricted, world);

    nbeta=nemo->get_calc()->param.nbeta;
    const bool have_beta=xc->is_spin_polarized() && nbeta != 0;

    // compute the alpha and beta densities
    real_function_3d arho,brho;
    real_function_3d arhonemo=nemo->make_density(nemo->get_calc()->aocc,nemo->get_calc()->amo);
    arho=(arhonemo*nemo->R_square).truncate(extra_truncation);
    if (have_beta) {
        real_function_3d brhonemo=nemo->make_density(nemo->get_calc()->bocc,nemo->get_calc()->bmo);
        brho=(brhonemo*nemo->R_square).truncate(extra_truncation);
    } else {
        brho=arho;
    }

    std::vector<real_function_3d> darho(3),dbrho(3);
    if (xc->is_gga()) {
        darho[0]=nemo->make_ddensity(arhonemo,0);
        darho[1]=nemo->make_ddensity(arhonemo,1);
        darho[2]=nemo->make_ddensity(arhonemo,2);
        if (have_beta) {
            real_function_3d brhonemo=nemo->make_density(nemo->get_calc()->bocc,nemo->get_calc()->bmo);
            dbrho[0]=nemo->make_ddensity(brhonemo,0);
            dbrho[1]=nemo->make_ddensity(brhonemo,1);
            dbrho[2]=nemo->make_ddensity(brhonemo,2);
        }
    }
    xc_args=prep_xc_args(arho,brho,darho,dbrho);
}


XCOperator::XCOperator(World& world, const SCF* calc, const real_function_3d& arho,
        const real_function_3d& brho, int ispin)
        : world(world), nbeta(calc->param.nbeta), ispin(ispin),
          extra_truncation(FunctionDefaults<3>::get_thresh()*0.01) {
    xc=std::shared_ptr<XCfunctional> (new XCfunctional());
    xc->initialize(calc->param.xc_data, !calc->param.spin_restricted, world);
    const bool have_beta=xc->is_spin_polarized() && nbeta != 0;


    std::vector<real_function_3d> darho,dbrho;
    if (xc->is_gga()) {
        darho=calc->nabla(arho);
        if (have_beta) dbrho=calc->nabla(brho);
    }

    xc_args=prep_xc_args(arho,brho,darho,dbrho);
}

XCOperator::XCOperator(World& world, const Nemo* nemo, const real_function_3d& arho,
        const real_function_3d& brho, int ispin) : world(world),
        nbeta(nemo->get_calc()->param.nbeta), ispin(ispin), extra_truncation(0.01) {
    xc=std::shared_ptr<XCfunctional> (new XCfunctional());
    xc->initialize(nemo->get_calc()->param.xc_data,
            not nemo->get_calc()->param.spin_restricted, world);
    const bool have_beta=xc->is_spin_polarized() && nbeta != 0;

    std::vector<real_function_3d> darho,dbrho;
    if (xc->is_gga()) {
        darho=nemo->get_calc()->nabla(arho);
        if (have_beta) dbrho=nemo->get_calc()->nabla(brho);
    }

    xc_args=prep_xc_args(arho,brho,darho,dbrho);
}

vecfuncT XCOperator::operator()(const vecfuncT& vket) const {
    real_function_3d xc_pot=make_xc_potential();
    double vtol = FunctionDefaults<3>::get_thresh() * 0.1;  // safety
    return mul_sparse(world, xc_pot, vket, vtol);
}

double XCOperator::compute_xc_energy() const {

    if (not is_initialized()) {
        MADNESS_EXCEPTION("calling xc energy without intermediates ",1);
    }

    refine_to_common_level(world,xc_args);
    real_function_3d vlda=multiop_values<double, xc_functional, 3>
            (xc_functional(*xc), xc_args);
    truncate(world,xc_args);

    return vlda.trace();
}


real_function_3d XCOperator::make_xc_potential() const {

    if (not is_initialized()) {
        MADNESS_EXCEPTION("calling xc potential without intermediates ",1);
    }

    refine_to_common_level(world,xc_args);

    // LDA/GGA local part
    real_function_3d dft_pot=multiop_values<double, xc_potential, 3>
                (xc_potential(*xc, ispin, XCfunctional::potential_rho), xc_args);
//    save(dft_pot,"lda_pot");

    // GGA part
    //
    // What = potential_rho : Vrho
    // What = potential_same_spin : Vsigma_ss
    // What = potential_mixed_spin : Vsigma_ab
    //
    // close shell
    //       v_xc = vrho - Div( 2Vsig_aa * Grad(rho_a))
    // open shell
    //       v_xc = vrho - Div( 2*Vsig_aa*Grad(rho)+ Vsig_ab*Grad(rho_b)
    //                   + Vsig_ba*Grad(rho_a) + 2*Vsig_bb*Grad(rho_b))
    //

    if (xc->is_gga() ) {
        const real_function_3d& rho=xc_args[XCfunctional::enum_rhoa];
        double tol=xc->get_ggatol();

        bool have_beta=xc->is_spin_polarized() && nbeta != 0;

        real_function_3d gga_pot=real_factory_3d(world).compressed();

        if (not xc->is_spin_polarized()) {      // RHF case
            MADNESS_ASSERT(ispin==0);
            // get Vsigma_aa
            functionT vsigaa = multiop_values<double, xc_potential, 3>
                (xc_potential(*xc, ispin, XCfunctional::potential_same_spin), xc_args); //.truncate();
//            save(vsigaa,"vsigaa");

            for (int axis=0; axis<3; axis++) {
                functionT gradn_alpha = xc_args[XCfunctional::enum_drhoa_x+axis];
//                functionT ddel = 4.0*vsigaa*gradn_alpha;    // fac 2 for formula, fac 2 to rho=2rho_\alpha
                functionT ddel1 = 4.0* vsigaa*gradn_alpha;    // fac 2 for formula, fac 2 to rho=2rho_\alpha
                real_function_3d ddel=binary_op(ddel1,rho,binary_munging(tol));
                Derivative<double,3> D = free_space_derivative<double,3>(world, axis);

//                save(ddel,"ddel"+stringify(axis));
//                Derivative<double,3> D = free_space_derivative<double,3>(world, axis);
                functionT vxc2=D(ddel);
//                (vxc2,"vxc2"+stringify(axis));
                gga_pot-=vxc2;//.truncate();
            }
//            save(gga_pot,"gga_pot");

        } else if (have_beta) {                                // UHF case

            // get Vsigma_aa or Vsigma_bb
            functionT vsigaa = multiop_values<double, xc_potential, 3>
                    (xc_potential(*xc, ispin, XCfunctional::potential_same_spin), xc_args); //.truncate();
            // get Vsigma_ab
            functionT vsigab= multiop_values<double, xc_potential, 3>
                    (xc_potential(*xc, ispin, XCfunctional::potential_mixed_spin), xc_args); //.truncate();

            for (int axis=0; axis<3; axis++) {
                real_function_3d drho_same_spin, drho_other_spin;
                if (ispin==0) {
                    drho_same_spin=xc_args[XCfunctional::enum_drhoa_x+axis];
                    drho_other_spin=xc_args[XCfunctional::enum_drhob_x+axis];
                } else {
                    drho_same_spin=xc_args[XCfunctional::enum_drhob_x+axis];
                    drho_other_spin=xc_args[XCfunctional::enum_drhoa_x+axis];
                }

                functionT ddel1 = 2.0* vsigaa*drho_same_spin + vsigab*drho_other_spin;
                real_function_3d ddel=binary_op(ddel1,rho,binary_munging(tol));
                Derivative<double,3> D = free_space_derivative<double,3>(world, axis);
                functionT vxc2=D(ddel);
                gga_pot-=vxc2;//.truncate();
            }
        }
//        save(gga_pot,"gga_pot");

        real_function_3d gga_pot_munged=binary_op(gga_pot,rho,binary_munging(tol));
//        save(gga_pot_munged,"gga_pot_munged");

        dft_pot+=gga_pot_munged;
    } //is gga

    truncate(world,xc_args);

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
//  the second partial derivatives that need to be multiplied with the density gradients
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
real_function_3d XCOperator::apply_xc_kernel(const real_function_3d& dens_pt) const {

    MADNESS_ASSERT(not xc->is_spin_polarized());    // for now
    MADNESS_ASSERT(ispin==0);           // for now

    if (not is_initialized()) {
        MADNESS_EXCEPTION("calling apply_xc_kernel without intermediates ",1);
    }

    vecfuncT ddens_pt;
    prep_xc_args_response(dens_pt, xc_args, ddens_pt);
    refine_to_common_level(world,xc_args);

    // compute the various terms from the xc kernel

    // compute the local terms: second_{local}
    real_function_3d result=multiop_values<double, xc_kernel_apply, 3>
            (xc_kernel_apply(*xc, ispin, XCfunctional::kernel_second_local), xc_args);
//    save(result,"local_apply");

    if (xc->is_gga()) {
        // compute the semilocal terms, second partial derivatives
        real_function_3d semilocal2a=multiop_values<double, xc_kernel_apply, 3>
                (xc_kernel_apply(*xc, ispin, XCfunctional::kernel_second_semilocal), xc_args);
        save(semilocal2a,"semilocal2a");

        const real_function_3d& rho=xc_args[XCfunctional::enum_rhoa];
        double tol=xc->get_ggatol();
        real_function_3d semilocal2=binary_op(semilocal2a,rho,binary_munging(tol));
//        real_function_3d semilocal2=semilocal2a;
        save(semilocal2,"semilocal2_ddel");

        for (int idim=0; idim<3; ++idim) {
            real_function_3d ddens = xc_args[XCfunctional::enum_drhoa_x+idim];
            real_function_3d ddel = 2.0*semilocal2*ddens;
            save(ddel,"ddel"+stringify(idim));
            Derivative<double,3> D = free_space_derivative<double,3>(world, idim);
            real_function_3d vxc2=D(ddel);
            save(vxc2,"vxc2"+stringify(idim));
            result-=vxc2;
        }

        // compute the semilocal terms, first partial derivative
        real_function_3d semilocal1a=multiop_values<double, xc_kernel_apply, 3>
                (xc_kernel_apply(*xc, ispin, XCfunctional::kernel_first_semilocal), xc_args);
        real_function_3d semilocal1=binary_op(semilocal1a,rho,binary_munging(tol));
//        real_function_3d semilocal1=semilocal1a;
//        save(semilocal1a,"semilocal1a");

        for (int idim=0; idim<3; ++idim) {
            real_function_3d ddel = semilocal1*ddens_pt[idim];
            Derivative<double,3> D = free_space_derivative<double,3>(world, idim);
            real_function_3d vxc2=D(ddel);
            result-=vxc2;
        }
        real_function_3d result1=binary_op(result,rho,binary_munging(tol));
        result=result1;
    }

    truncate(world,xc_args);
    return result.truncate();
}

/// prepare xc args
vecfuncT XCOperator::prep_xc_args(const real_function_3d& arho,
        const real_function_3d& brho, const std::vector<real_function_3d>& darho,
        const std::vector<real_function_3d>& dbrho) const {

    World& world=arho.world();
    vecfuncT xcargs(XCfunctional::number_xc_args);
    const bool have_beta=(xc->is_spin_polarized()) and (nbeta>0);

    // assign the densities (alpha, beta)
    arho.reconstruct();
    xcargs[XCfunctional::enum_rhoa]=copy(arho);      // alpha density
//    save(xcargs[XCfunctional::enum_rhoa],"rho");

    if (have_beta) {
        brho.reconstruct();
        xcargs[XCfunctional::enum_rhob]=copy(brho);  // beta density
    }
    world.gop.fence();

    if (xc->is_gga()) {
        // compute the gradients of the densities
        std::vector< std::shared_ptr<Derivative<double,3> > > gradop =
                gradient_operator<double,3>(world);

        // assign the gradients of the densities
//        xcargs[XCfunctional::enum_drhoa_x]=(*gradop[0])(arho, false);
//        xcargs[XCfunctional::enum_drhoa_y]=(*gradop[1])(arho, false);
//        xcargs[XCfunctional::enum_drhoa_z]=(*gradop[2])(arho, false);
        xcargs[XCfunctional::enum_drhoa_x]=darho[0];
        xcargs[XCfunctional::enum_drhoa_y]=darho[1];
        xcargs[XCfunctional::enum_drhoa_z]=darho[2];

        if (have_beta) {
//            xcargs[XCfunctional::enum_drhob_x]=(*gradop[0])(brho, false);
//            xcargs[XCfunctional::enum_drhob_y]=(*gradop[1])(brho, false);
//            xcargs[XCfunctional::enum_drhob_z]=(*gradop[2])(brho, false);
            xcargs[XCfunctional::enum_drhob_x]=dbrho[0];
            xcargs[XCfunctional::enum_drhob_y]=dbrho[1];
            xcargs[XCfunctional::enum_drhob_z]=dbrho[2];
        }
        world.gop.fence();

        // autorefine before squaring
        real_function_3d drhoa_x=copy(xcargs[XCfunctional::enum_drhoa_x]).refine();
        real_function_3d drhoa_y=copy(xcargs[XCfunctional::enum_drhoa_y]).refine();
        real_function_3d drhoa_z=copy(xcargs[XCfunctional::enum_drhoa_z]).refine();

        // assign the reduced densities sigma
        xcargs[XCfunctional::enum_saa]=        // sigma_aa
                (drhoa_x * drhoa_x + drhoa_y * drhoa_y + drhoa_z * drhoa_z);
        if (have_beta) {

            // autorefine before squaring
            real_function_3d drhob_x=copy(xcargs[XCfunctional::enum_drhob_x]).refine();
            real_function_3d drhob_y=copy(xcargs[XCfunctional::enum_drhob_y]).refine();
            real_function_3d drhob_z=copy(xcargs[XCfunctional::enum_drhob_z]).refine();

            xcargs[XCfunctional::enum_sab]=    // sigma_ab
                    (drhoa_x * drhob_x + drhoa_y * drhob_y + drhoa_z * drhob_z);
            xcargs[XCfunctional::enum_sbb]=    // sigma_bb
                    (drhob_x * drhob_x + drhob_y * drhob_y + drhob_z * drhob_z);

            // this is needed for proper munging of sigma_ab
            xcargs[XCfunctional::enum_sigtot]= xcargs[XCfunctional::enum_saa]
                           +2.0*xcargs[XCfunctional::enum_sab]+xcargs[XCfunctional::enum_sbb];
        }
        world.gop.fence();
//        save(xcargs[XCfunctional::enum_saa],"sigma");

    }

    world.gop.fence();

    truncate(world,xc_args,extra_truncation);

    return xcargs;
}

/// add intermediates for the response kernels to xc_args
void XCOperator::prep_xc_args_response(const real_function_3d& dens_pt,
        vecfuncT& xc_args, vecfuncT& ddens_pt) const {

    const bool have_beta=(xc->is_spin_polarized()) and (nbeta>0);
    World& world=dens_pt.world();

    const real_function_3d& arho=xc_args[XCfunctional::enum_rhoa];
    const real_function_3d& brho=xc_args[XCfunctional::enum_rhob];
    ddens_pt=vecfuncT(3);

    // assign the perturbed density (spin-free ??)
    xc_args[XCfunctional::enum_rho_pt]=dens_pt;
    world.gop.fence();

    // assign the reduced density gradients with the perturbed density for GGA
    // \sigma_pt_a = \nabla \rho_\alpha \cdot \nabla\tilde\rho
    // \sigma_pt_b = \nabla \rho_\beta \cdot \nabla\tilde\rho
    if (xc->is_gga()) {

        std::vector< std::shared_ptr<Derivative<double,3> > > gradop =
                gradient_operator<double,3>(world);
        arho.reconstruct(false);
        brho.reconstruct(false);
        dens_pt.reconstruct(false);
        world.gop.fence();

        std::vector<real_function_3d> ddensa(3),ddensb(3);
        for (std::size_t i=0; i<3; ++i) {
            ddensa[i]=(*gradop[i])(arho, false);
            if (have_beta) ddensb[i]=(*gradop[i])(brho, false);
            ddens_pt[i]=(*gradop[i])(dens_pt, false);
        }
        world.gop.fence();

        // autorefine before squaring
        for (std::size_t i=0; i<3; ++i) {
            ddensa[i].refine(false);
            if (have_beta) ddensb[i].refine(false);
            ddens_pt[i].refine(false);
        }
        world.gop.fence();

        xc_args[XCfunctional::enum_sigma_pta]=dot(world,ddensa,ddens_pt);    // sigma_a
        if (have_beta) xc_args[XCfunctional::enum_sigma_ptb]= dot(world,ddensb,ddens_pt);
//        save(xc_args[XCfunctional::enum_sigma_pta],"sigma_pt");
        world.gop.fence();
    }
    world.gop.fence();
    truncate(world,xc_args,extra_truncation);
}


void XCOperator::prep_xc_args_old(const real_function_3d& arho,
        const real_function_3d& brho, vecfuncT& delrho, vecfuncT& vf) const {

    delrho.clear();
    vf.clear();

    arho.reconstruct();
    if (nbeta != 0 && xc->is_spin_polarized()) brho.reconstruct();

    vf.push_back(arho);
    if (xc->is_spin_polarized()) vf.push_back(brho);

    if (xc->is_gga()) {

        std::vector< std::shared_ptr<Derivative<double,3> > > gradop =
                gradient_operator<double,3>(world);

        for (int axis = 0; axis < 3; ++axis)
            delrho.push_back((*gradop[axis])(arho, false)); // delrho
        if (xc->is_spin_polarized() && nbeta != 0)
            for (int axis = 0; axis < 3; ++axis)
                delrho.push_back((*gradop[axis])(brho, false));

        world.gop.fence(); // NECESSARY

        vf.push_back(delrho[0] * delrho[0] + delrho[1] * delrho[1]
                  + delrho[2] * delrho[2]);     // sigma_aa

        if (xc->is_spin_polarized()) {
            if (nbeta != 0) {
                vf.push_back(delrho[0] * delrho[3] + delrho[1] * delrho[4]
                           + delrho[2] * delrho[5]); // sigma_ab
                vf.push_back(delrho[3] * delrho[3] + delrho[4] * delrho[4]
                           + delrho[5] * delrho[5]); // sigma_bb
            } else {
                vf.push_back(real_function_3d());
                vf.push_back(real_function_3d());
            }
        }

        world.gop.fence(); // NECESSARY
    }
    if (vf.size()) {
        reconstruct(world, vf);
        refine_to_common_level(world,vf); // Ugly but temporary (I hope!)
    }

//    // this is a nasty hack, just adding something so that make_libxc_args
//    // receives 5 arguments has to be here or refine_to_common_level(vf) above
//    // hangs, but we really need a better solution for when nbeta=0
//    if (xc->is_spin_polarized() && nbeta == 0 && xc->is_gga()){
//        vf.push_back(brho);
//        vf.push_back(brho);
//    }
}

bool XCOperator::is_initialized() const {
    return (xc_args.size()>0);
}


Fock::Fock(World& world, const SCF* calc,
        std::shared_ptr<NuclearCorrelationFactor> ncf )
    : world(world),
      J(world,calc),
      K(world,calc,0),
      T(world),
      V(world,ncf) {
}



} // namespace madness


