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

/// \file SCFOperators.h
/// \brief Operators for the molecular HF and DFT code
/// \defgroup chem The molecular density functional and Hartree-Fock code


#ifndef MADNESS_CHEM_SCFOPERATORS_H_
#define MADNESS_CHEM_SCFOPERATORS_H_

#include <madness.h>
using namespace madness;

namespace madness {

// forward declaration
class SCF;
class NuclearCorrelationFactor;

typedef std::vector<real_function_3d> vecfuncT;

template<typename T, std::size_t NDIM>
class Kinetic {
    typedef DistributedMatrix<T> distmatT;
    typedef Function<T,NDIM> functionT;
    typedef std::vector<functionT> vecfuncT;
    typedef Tensor<T> tensorT;

public:
    Kinetic(World& world) : world(world) {
        gradop = gradient_operator<T,NDIM>(world);
    }

    functionT operator()(const functionT& ket) const {
        MADNESS_EXCEPTION("do not apply the kinetic energy operator on a function!",1);
        return ket;
    }

    T operator()(const functionT& bra, const functionT ket) const {
        T ke = 0.0;
        for (std::size_t axis = 0; axis < NDIM; axis++) {
            Derivative<T,NDIM> D = free_space_derivative<T,NDIM>(world,axis);
            const functionT dket = D(ket);
            const functionT dbra = D(bra);
            ke += 0.5 * (inner(dket, dbra));
        }
        return ke;
    }

    tensorT operator()(const vecfuncT& vbra, const vecfuncT& vket) const {
        distmatT dkinetic;
        if (&vbra==&vket) {
            dkinetic = kinetic_energy_matrix(world,vbra);
        } else {
            dkinetic = kinetic_energy_matrix(world,vbra,vket);
        }
        tensorT kinetic(vbra.size(),vket.size());
        dkinetic.copy_to_replicated(kinetic);
        return kinetic;
    }

private:
    World& world;
    std::vector< std::shared_ptr<Derivative<T,NDIM>> > gradop;

    distmatT kinetic_energy_matrix(World & world, const vecfuncT & v) const;
    distmatT kinetic_energy_matrix(World & world, const vecfuncT & vbra,
            const vecfuncT & vket) const;

};


class Coulomb {
public:

    Coulomb(World& world) : world(world) {};

    Coulomb(World& world, const SCF* calc);

    real_function_3d operator()(const real_function_3d& ket) const {
        return (vcoul*ket).truncate();
    }

    vecfuncT operator()(const vecfuncT& vket) const {
        vecfuncT tmp=mul(world,vcoul,vket);
        truncate(world,tmp);
        return tmp;
    }

    double operator()(const real_function_3d& bra, const real_function_3d ket) const {
        return inner(bra,vcoul*ket);
    }

    Tensor<double> operator()(const vecfuncT& vbra, const vecfuncT& vket) const {
        vecfuncT vJket;
        for (std::size_t i=0; i<vket.size(); ++i) {
            vJket.push_back(this->operator()(vket[i]));
        }
        return matrix_inner(world,vbra,vJket);
    }

    const real_function_3d& potential() const {return vcoul;}
    real_function_3d& potential() {return vcoul;}
    real_function_3d compute_potential(const real_function_3d& density,
            double lo=1.e-4, double econv=FunctionDefaults<3>::get_thresh()) const;
    real_function_3d compute_potential(const SCF* calc) const;

private:
    real_function_3d vcoul; ///< the coulomb potential
    World& world;
};


class Nuclear {
public:
    Nuclear(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf)
        : world(world), ncf(ncf) {}

    real_function_3d operator()(const real_function_3d& ket) const;

    vecfuncT operator()(const vecfuncT& vket) const {
        vecfuncT result(vket.size());
        for (std::size_t i=0; i<vket.size(); ++i) result[i]=this->operator()(vket[i]);
        truncate(world,result);
        return result;
    }

    double operator()(const real_function_3d& bra, const real_function_3d& ket) const {
        return inner(bra,this->operator()(ket));
    }

    Tensor<double> operator()(const vecfuncT& vbra, const vecfuncT& vket) const {
        vecfuncT vVket;
        for (std::size_t i=0; i<vket.size(); ++i) {
            vVket.push_back(this->operator()(vket[i]));
        }
        return matrix_inner(world,vbra,vVket);
    }

private:
    World& world;
    std::shared_ptr<NuclearCorrelationFactor> ncf;
};


class Exchange {
public:
    Exchange(World& world, const vecfuncT& amo, const Tensor<double> occ,
            const SCF* calc, const real_function_3d& R2 );

    real_function_3d operator()(const real_function_3d& ket) const {
        real_function_3d result = real_factory_3d(world).compressed(true);
        real_function_3d R2ket=R2*ket;
        for (std::size_t k = 0; k < mo.size(); ++k) {
            real_function_3d ik = mo[k] * R2ket;
            result += mo[k] * (*poisson)(ik);
        }
        return result;
    }

    /// apply the exchange operator on a vector of functions

    /// note that only one spin is used (either alpha or beta orbitals)
    /// @param[in]  vket       the orbitals |i> that the operator is applied on
    /// @return     a vector of orbitals  K| i>
    vecfuncT operator()(const vecfuncT& vket) const;

//    vecfuncT operator()(const vecfuncT& vket) const {
//        vecfuncT result(vket.size());
//        for (std::size_t i=0; i<vket.size(); ++i) result[i]=this->operator()(vket[i]);
//        truncate(world,result);
//        return result;
//    }

    double operator()(const real_function_3d& bra, const real_function_3d ket) const {
        return inner(bra,this->operator()(ket));
    }

    Tensor<double> operator()(const vecfuncT& vbra, const vecfuncT& vket) const {
        vecfuncT vKket;
        for (std::size_t i=0; i<vket.size(); ++i) {
            vKket.push_back(this->operator()(vket[i]));
        }
        return matrix_inner(world,vbra,vKket);
    }

private:
    World& world;
    bool small_memory;
    const vecfuncT mo;
    const Tensor<double> occ;
    const real_function_3d R2;
    std::shared_ptr<real_convolution_3d> poisson;
};




class Fock {
public:
    Fock(World& world, const SCF* calc, std::shared_ptr<NuclearCorrelationFactor> ncf);

    real_function_3d operator()(const real_function_3d& ket) const {
        real_function_3d result;
        return result;
    }
    double operator()(const real_function_3d& bra, const real_function_3d ket) const {
        double J_00 = J(bra,ket);
        double K_00 = K(bra,ket);
        double T_00 = T(bra,ket);
        double V_00 = V(bra,ket);
        return T_00 + J_00 - K_00 + V_00;
    }

    Tensor<double> operator()(const vecfuncT& vbra, const vecfuncT& vket) const {
        double wtime=-wall_time(); double ctime=-cpu_time();
        Tensor<double> kmat=K(vbra,vket);
        Tensor<double> jmat=J(vbra,vket);
        Tensor<double> tmat=T(vbra,vket);
        Tensor<double> vmat=V(vbra,vket);
        Tensor<double> fock=tmat+jmat-kmat+vmat;
        wtime+=wall_time(); ctime+=cpu_time();
        if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", "fock matrix", wtime, ctime);
        return fock;
    }


private:
    World& world;
    Coulomb J;
    Exchange K;
    Kinetic<double,3> T;
    Nuclear V;
};

}
#endif /* MADNESS_CHEM_SCFOPERATORS_H_ */
