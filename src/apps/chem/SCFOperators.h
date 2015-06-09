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
class Nemo;
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

    vecfuncT operator()(const vecfuncT& vket) const {
        MADNESS_EXCEPTION("do not apply the kinetic energy operator on a function!",1);
        return vket;
    }

    T operator()(const functionT& bra, const functionT ket) const {
        vecfuncT vbra(1,bra), vket(1,ket);
        Tensor<T> tmat=this->operator()(vbra,vket);
        return tmat(0l,0l);
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
    std::vector< std::shared_ptr<Derivative<T,NDIM> > > gradop;

    distmatT kinetic_energy_matrix(World & world, const vecfuncT & v) const;
    distmatT kinetic_energy_matrix(World & world, const vecfuncT & vbra,
            const vecfuncT & vket) const;

};


class Coulomb {
public:

    /// default empty ctor
    Coulomb(World& world) : world(world) {};

    /// ctor with an SCF calculation providing the MOs and density
    Coulomb(World& world, const SCF* calc) : world(world) {
        vcoul=compute_potential(calc);
    }

    /// ctor with an SCF calculation providing the MOs and density
    Coulomb(World& world, const Nemo* nemo);

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

    /// getter for the Coulomb potential
    const real_function_3d& potential() const {return vcoul;}

    /// setter for the Coulomb potential
    real_function_3d& potential() {return vcoul;}

    real_function_3d compute_density(const SCF* calc) const;

    /// given a density compute the Coulomb potential

    /// this function uses a newly constructed Poisson operator. Note that
    /// the accuracy parameters must be consistent with the exchange operator.
    real_function_3d compute_potential(const real_function_3d& density,
            double lo=1.e-4, double econv=FunctionDefaults<3>::get_thresh()) const;

    /// given a set of MOs in an SCF calculation, compute the Coulomb potential

    /// this function uses the Poisson operator of the SCF calculation
    real_function_3d compute_potential(const SCF* calc) const;

    /// given a set of MOs in an SCF calculation, compute the Coulomb potential

    /// this function uses the Poisson operator of the SCF calculation
    real_function_3d compute_potential(const Nemo* nemo) const;

private:
    World& world;
    real_function_3d vcoul; ///< the coulomb potential
    real_function_3d R_square;    ///< square of the nuclear correlation factor, if any
};


class Nuclear {
public:

    Nuclear(World& world, const SCF* calc);

    Nuclear(World& world, const Nemo* nemo);

    Nuclear(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf)
        : world(world), ncf(ncf) {}


    real_function_3d operator()(const real_function_3d& ket) const {
        vecfuncT vket(1,ket);
        return this->operator()(vket)[0];
    }

    vecfuncT operator()(const vecfuncT& vket) const;

    double operator()(const real_function_3d& bra, const real_function_3d& ket) const {
        return inner(bra,this->operator()(ket));
    }

    Tensor<double> operator()(const vecfuncT& vbra, const vecfuncT& vket) const {
        vecfuncT vVket=this->operator()(vket);
        return matrix_inner(world,vbra,vVket);
    }

private:
    World& world;
    std::shared_ptr<NuclearCorrelationFactor> ncf;

};


class Exchange {
public:

    /// default ctor
    Exchange(World& world) : world(world), small_memory_(true), same_(false) {};

    /// ctor with a conventional calculation
    Exchange(World& world, const SCF* calc, const int ispin);

    /// ctor with a nemo calculation
    Exchange(World& world, const Nemo* nemo, const int ispin);

    void set_parameters(const vecfuncT& bra, const vecfuncT& ket,
            const Tensor<double>& occ, const double lo=1.e-4,
            const double econv=FunctionDefaults<3>::get_thresh());

    real_function_3d operator()(const real_function_3d& ket) const {
        vecfuncT vket(1,ket);
        vecfuncT vKket=this->operator()(vket);
        return vKket[0];
    }

    /// apply the exchange operator on a vector of functions

    /// note that only one spin is used (either alpha or beta orbitals)
    /// @param[in]  vket       the orbitals |i> that the operator is applied on
    /// @return     a vector of orbitals  K| i>
    vecfuncT operator()(const vecfuncT& vket) const;

    /// compute the matrix element <bra | K | ket>

    /// @param[in]  bra    real_funtion_3d, the bra state
    /// @param[in]  ket    real_funtion_3d, the ket state
    double operator()(const real_function_3d& bra, const real_function_3d ket) const {
        return inner(bra,this->operator()(ket));
    }

    /// compute the matrix < vbra | K | vket >

    /// @param[in]  vbra    vector of real_funtion_3d, the set of bra states
    /// @param[in]  vket    vector of real_funtion_3d, the set of ket states
    /// @return K_ij
    Tensor<double> operator()(const vecfuncT& vbra, const vecfuncT& vket) const {
        vecfuncT vKket=this->operator()(vket);
        return matrix_inner(world,vbra,vKket);
    }

    bool& small_memory() {return small_memory_;}
    bool small_memory() const {return small_memory_;}
    Exchange& small_memory(const bool flag) {
        small_memory_=flag;
        return *this;
    }

    bool& same() {return same_;}
    bool same() const {return same_;}
    Exchange& same(const bool flag) {
        same_=flag;
        return *this;
    }

private:

    World& world;
    bool small_memory_;
    bool same_;
    vecfuncT mo_bra, mo_ket;    ///< MOs for bra and ket
    Tensor<double> occ;
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
