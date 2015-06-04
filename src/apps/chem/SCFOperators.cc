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
#include <apps/chem/correlationfactor.h>


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

Coulomb::Coulomb(World& world, const madness::SCF* calc) : world(world) {

    real_function_3d density = calc->make_density(world, calc->get_aocc(),
            calc->get_amo());
    if (calc->is_spin_restricted()) {
        density.scale(2.0);
    } else {
        real_function_3d brho = calc->make_density(world, calc->get_bocc(),
                calc->get_bmo());
        density+=brho;
    }
    vcoul=calc->make_coulomb_potential(density);
}

real_function_3d Coulomb::compute_potential(const real_function_3d& density,
        double lo, double econv) const {
    real_convolution_3d poisson = CoulombOperator(world, lo, econv);
    return poisson(density);
}


real_function_3d Nuclear::operator()(const real_function_3d& ket) const {
    return ncf->apply_U(ket);
}

Exchange::Exchange(World& world, const vecfuncT& mo,
        const Tensor<double> occ, const SCF* calc,
        const real_function_3d& R2 ) : world(world), small_memory(true),
                mo(mo), occ(occ), R2(R2) {
    poisson = std::shared_ptr<real_convolution_3d>(
            CoulombOperatorPtr(world, calc->param.lo, calc->param.econv));
}

vecfuncT Exchange::operator()(const vecfuncT& vket) const {
    const bool same = (&mo == &vket);
    int nocc = mo.size();
    int nf = vket.size();
    double tol = FunctionDefaults < 3 > ::get_thresh(); /// Important this is consistent with Coulomb
    vecfuncT Kf = zero_functions_compressed<double, 3>(world, nf);
    reconstruct(world, mo);
    norm_tree(world, mo);
    if (!same) {
        reconstruct(world, vket);
        norm_tree(world, vket);
    }

    if (small_memory) {     // Smaller memory algorithm ... possible 2x saving using i-j sym
        for(int i=0; i<nocc; ++i){
            if(occ[i] > 0.0){
                vecfuncT psif = mul_sparse(world, mo[i], vket, tol); /// was vtol
                truncate(world, psif);
                psif = apply(world, *poisson.get(), psif);
                truncate(world, psif);
                psif = mul_sparse(world, mo[i], psif, tol); /// was vtol
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
                psif.push_back(mul_sparse(mo[i], vket[j], tol, false));
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
                psipsif[i * nf + j] = mul_sparse(psif[ij], mo[i], false);
                if (same && i != j) {
                    psipsif[j * nf + i] = mul_sparse(psif[ij], mo[j], false);
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


Fock::Fock(World& world, const SCF* calc,
        std::shared_ptr<NuclearCorrelationFactor> ncf )
    : world(world),
      J(world,calc),
      K(world,calc->amo,calc->aocc,calc,ncf->square()),
      T(world),
      V(world,ncf) {
}



} // namespace madness


