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
        r += matrix_inner(r.distribution(), dvbra[i], dvket[i], true);
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
    real_function_3d u2x_f=real_factory_3d(world).functor2(u2x)
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
            real_function_3d U1=real_factory_3d(world).functor2(u1x).truncate_on_project();
            std::vector<Function<double,NDIM> > U1dv=mul(world,U1,dv);
            truncate(world,U1dv);
            result=sub(world,result,U1dv);
            truncate(world,result);
        }

        // add the U3X potential
        NuclearCorrelationFactor::U3X_functor u3x(ncf.get(),iatom,iaxis);
        real_function_3d u3x_f=real_factory_3d(world).functor2(u3x).truncate_on_project();
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


XCOperator::XCOperator(World& world, const SCF* calc, int ispin) : world(world),
        ispin(ispin) {
    xc=std::shared_ptr<XCfunctional> (new XCfunctional());
    xc->initialize(calc->param.xc_data, !calc->param.spin_restricted, world);
    nbeta=calc->param.nbeta;

    // compute the alpha and beta densities
    real_function_3d arho,brho;
    arho=calc->make_density(world,calc->aocc,calc->amo);
    if (xc->is_spin_polarized() && nbeta != 0) {
        brho=calc->make_density(world,calc->bocc,calc->bmo);
    } else {
        brho=arho;
    }
    prep_xc_args(arho,brho,delrho,vf);
}

XCOperator::XCOperator(World& world, const Nemo* nemo, int ispin) : world(world),
        ispin(ispin) {
    xc=std::shared_ptr<XCfunctional> (new XCfunctional());
    xc->initialize(nemo->get_calc()->param.xc_data,
            !nemo->get_calc()->param.spin_restricted, world);
    nbeta=nemo->get_calc()->param.nbeta;

    // compute the alpha and beta densities
    real_function_3d arho,brho;
    arho=nemo->make_density(nemo->get_calc()->aocc,nemo->get_calc()->amo);
    arho=(arho*nemo->R_square).truncate();
    if (xc->is_spin_polarized() && nbeta != 0) {
        brho=nemo->make_density(nemo->get_calc()->bocc,nemo->get_calc()->bmo);
        brho=(brho*nemo->R_square).truncate();
    } else {
        brho=arho;
    }
    prep_xc_args(arho,brho,delrho,vf);
}


XCOperator::XCOperator(World& world, const SCF* calc, const real_function_3d& arho,
        const real_function_3d& brho, int ispin)
        : world(world), nbeta(calc->param.nbeta), ispin(ispin) {
    xc=std::shared_ptr<XCfunctional> (new XCfunctional());
    xc->initialize(calc->param.xc_data, !calc->param.spin_restricted, world);
    prep_xc_args(arho,brho,delrho,vf);
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

    // actually independent of ispin!! (legacy calling list?)
    real_function_3d vlda=multiop_values<double, xc_functional, 3>
            (xc_functional(*xc, ispin), vf);
    return vlda.trace();
}


real_function_3d XCOperator::make_xc_potential() const {

    if (not is_initialized()) {
        MADNESS_EXCEPTION("calling xc potential without intermediates ",1);
    }
    // LDA part
    real_function_3d dft_pot=multiop_values<double, xc_potential, 3>
                (xc_potential(*xc, ispin, 0), vf);

    // GGA part
    //
    // What = 0 : Vrho
    // What = 1 : Vsigma_ss
    // What = 2 : Vsigma_ab
    //
    // close shell
    //       v_xc = vrho - Div( 2Vsig_aa * Grad(rho_a))
    // open shell
    //       v_xc = vrho - Div( 2*Vsig_aa*Grad(rho)+ Vsig_ab*Grad(rho_b)
    //                   + Vsig_ba*Grad(rho_a) + 2*Vsig_bb*Grad(rho_b))
    //

    if (xc->is_gga() ) {
        // get Vsigma_aa (if it is the case and Vsigma_bb)
        functionT vsigaa = multiop_values<double, xc_potential, 3>
        (xc_potential(*xc, ispin, 1), vf); //.truncate();
        functionT vsigab;
        if (xc->is_spin_polarized() && nbeta != 0)// V_ab
            vsigab = multiop_values<double, xc_potential, 3>
                    (xc_potential(*xc, ispin, 2), vf); //.truncate();

        for (int axis=0; axis<3; axis++) {
            functionT gradn = delrho[axis + 3*ispin];
            functionT ddel = vsigaa*gradn;
            if (xc->is_spin_polarized() && nbeta != 0) {
                functionT vsab = vsigab*delrho[axis + 3*(1-ispin)];
                ddel = ddel + vsab;
            }
            ddel.scale(xc->is_spin_polarized() ? 2.0 : 4.0);
            Derivative<double,3> D = free_space_derivative<double,3>(world, axis);
            functionT vxc2=D(ddel);
            dft_pot-=vxc2;//.truncate();
        }
    } //is gga
    world.gop.fence();
    return dft_pot;
}

real_function_3d XCOperator::make_xc_kernel() const {
    MADNESS_EXCEPTION("no make_xc_kernel yet",1);
    return multiop_values<double, xc_kernel, 3>(xc_kernel(*xc, ispin, 0), vf);
}

void XCOperator::prep_xc_args(const real_function_3d& arho,
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

        if (xc->is_spin_polarized() && nbeta != 0) {
            vf.push_back(delrho[0] * delrho[3] + delrho[1] * delrho[4]
                       + delrho[2] * delrho[5]); // sigma_ab
            vf.push_back(delrho[3] * delrho[3] + delrho[4] * delrho[4]
                       + delrho[5] * delrho[5]); // sigma_bb
        }

        world.gop.fence(); // NECESSARY
    }
    if (vf.size()) {
        reconstruct(world, vf);
        vf[0].refine_to_common_level(vf); // Ugly but temporary (I hope!)
    }

    // this is a nasty hack, just adding something so that make_libxc_args
    // receives 5 arguments has to be here or refine_to_common_level(vf) above
    // hangs, but we really need a better solution for when nbeta=0
    if (xc->is_spin_polarized() && nbeta == 0 && xc->is_gga()){
        vf.push_back(brho);
        vf.push_back(brho);
    }
}

bool XCOperator::is_initialized() const {
    bool cond=(vf.size()>0);
    if (not xc->is_lda()) cond=(cond and (delrho.size()>0));
    return cond;
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


