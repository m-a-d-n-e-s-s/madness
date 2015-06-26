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

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES

/*!
 \file examples/nemo.h
 \brief solve the HF equations using numerical exponential MOs

 The source is
 <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
 /trunk/src/apps/examples/nemo.h>here</a>.

 */

#ifndef NEMO_H_
#define NEMO_H_

#include <madness/mra/mra.h>
#include <madness/mra/funcplot.h>
#include <madness/mra/operator.h>
#include <madness/mra/lbdeux.h>
#include <chem/SCF.h>
#include <chem/correlationfactor.h>
#include <chem/molecular_optimizer.h>
#include <examples/nonlinsol.h>
#include <madness/mra/vmra.h>

namespace madness {

class PNO;

// this class needs to be moved to vmra.h !!

// This class is used to store information for the non-linear solver
template<typename T, std::size_t NDIM>
class vecfunc {
public:
	World& world;
	std::vector<Function<T, NDIM> > x;

	vecfunc(World& world, const std::vector<Function<T, NDIM> >& x1) :
			world(world), x(x1) {
	}

	vecfunc(const std::vector<Function<T, NDIM> >& x1) :
			world(x1[0].world()), x(x1) {
	}

	vecfunc(const vecfunc& other) :
			world(other.world), x(other.x) {
	}

	vecfunc& operator=(const vecfunc& other) {
		x = other.x;
		return *this;
	}

	vecfunc operator-(const vecfunc& b) const {
		return vecfunc(world, sub(world, x, b.x));
	}

	vecfunc operator+=(const vecfunc& b) { // Operator+= necessary
		x = add(world, x, b.x);
		return *this;
	}

	vecfunc operator*(double a) const { // Scale by a constant necessary

		PROFILE_BLOCK(Vscale);
		std::vector<Function<T,NDIM> > result(x.size());
		for (unsigned int i=0; i<x.size(); ++i) {
			result[i]=mul(a,x[i],false);
		}
		world.gop.fence();

//		scale(world, x, a);
		return result;
	}
};

/// the non-linear solver requires an inner product
template<typename T, std::size_t NDIM>
T inner(const vecfunc<T, NDIM>& a, const vecfunc<T, NDIM>& b) {
	Tensor<T> i = inner(a.world, a.x, b.x);
	return i.sum();
}

// The default constructor for functions does not initialize
// them to any value, but the solver needs functions initialized
// to zero for which we also need the world object.
template<typename T, std::size_t NDIM>
struct allocator {
	World& world;
	const int n;

	/// @param[in]	world	the world
	/// @param[in]	nn		the number of functions in a given vector
	allocator(World& world, const int nn) :
			world(world), n(nn) {
	}

	/// allocate a vector of n empty functions
	vecfunc<T, NDIM> operator()() {
		return vecfunc<T, NDIM>(world, zero_functions<T, NDIM>(world, n));
	}
};

/// The Nemo class
class Nemo: public MolecularOptimizationTargetInterface {
	typedef std::shared_ptr<real_convolution_3d> poperatorT;
	friend class PNO;
public:

	/// ctor

	/// @param[in]	world1	the world
	/// @param[in]	calc	the SCF
	Nemo(World& world1, std::shared_ptr<SCF> calc) :
			world(world1), calc(calc), ttt(0.0), sss(0.0), coords_sum(-1.0) {}

	double value() {return value(calc->molecule.get_all_coords());}

	double value(const Tensor<double>& x);

	/// compute the nuclear gradients
	Tensor<double> gradient(const Tensor<double>& x);

	bool provides_gradient() const {return true;}

	/// returns the molecular hessian matrix at structure x
	Tensor<double> hessian(const Tensor<double>& x);

	/// solve the CPHF equations for the nuclear displacements

	/// @param[in]  iatom   the atom A to be moved
	/// @param[in]  iaxis   the coordinate X of iatom to be moved
	/// @return     \frac{\partial}{\partial X_A} \varphi
	vecfuncT cphf(const int iatom, const int iaxis,
	        const vecfuncT& guess=vecfuncT()) const;

	/// solve teh CPHF equation for all displacements
	std::vector<vecfuncT> compute_all_cphf() const;

    /// solve the CPHF equations for the nuclear displacements

	/// works if no nuclear correlation factor is employed
    /// @param[in]  iatom   the atom A to be moved
    /// @param[in]  iaxis   the coordinate X of iatom to be moved
    /// @return     \frac{\partial}{\partial X_A} \varphi
    vecfuncT cphf_no_ncf(const int iatom, const int iaxis,
            const vecfuncT& guess=vecfuncT()) const;

    /// compute the perturbed density for CPHF
    real_function_3d compute_perturbed_density(const vecfuncT& mo,
            const vecfuncT& xi) const;

	/// returns the vibrational frequencies

	/// @param[in]  hessian the hessian matrix
    /// @param[in]  project whether to project out translation and rotation
    /// @param[in]  print_hessian   whether to print the hessian matrix
	/// @return the frequencies in atomic units
	Tensor<double> compute_frequencies(const Tensor<double>& hessian,
	        const bool project, const bool print_hessian) const;

	/// compute the mass-weight the hessian matrix

	/// use as: hessian.emul(massweights);
	/// @param[in]  molecule    for getting access to the atomic masses
	/// @return the mass-weighting matrix for the hessian
	Tensor<double> massweights(const Molecule& molecule) const;

	std::shared_ptr<SCF> get_calc() const {return calc;}

	/// compute the Fock matrix from scratch
	tensorT compute_fock_matrix(const vecfuncT& nemo, const tensorT& occ) const;

	/// return a reference to the molecule
	Molecule& molecule() {return calc->molecule;}

    /// return a reference to the molecule
    Molecule& molecule() const {
        return calc->molecule;
    }

    /// make the density (alpha or beta)
    real_function_3d make_density(World& world, const Tensor<double>& occ,
            const vecfuncT& nemo) const;

private:

	/// the world
	World& world;

	std::shared_ptr<SCF> calc;

	mutable double ttt, sss;
	void START_TIMER(World& world) const {
	    world.gop.fence(); ttt=wall_time(); sss=cpu_time();
	}

	void END_TIMER(World& world, const std::string msg) const {
	    END_TIMER(world,msg.c_str());
	}

	void END_TIMER(World& world, const char* msg) const {
	    ttt=wall_time()-ttt; sss=cpu_time()-sss;
	    if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
	}

public:

	/// the nuclear correlation factor
	std::shared_ptr<NuclearCorrelationFactor> nuclear_correlation;

	/// the nuclear correlation factor
	real_function_3d R;

	/// the inverse nuclear correlation factor
	real_function_3d R_inverse;

    /// the square of the nuclear correlation factor
    real_function_3d R_square;

private:

	/// sum of square of coords at last solved geometry
	mutable double coords_sum;

	/// a poisson solver
	std::shared_ptr<real_convolution_3d> poisson;

	void print_nuclear_corrfac() const;

	/// solve the HF equations
	double solve();

	/// given nemos, compute the HF energy
	double compute_energy(const vecfuncT& psi, const vecfuncT& Jpsi,
			const vecfuncT& Kpsi) const;

    /// given nemos, compute the HF energy using the regularized expressions for T and V
    double compute_energy_regularized(const vecfuncT& nemo, const vecfuncT& Jnemo,
            const vecfuncT& Knemo, const vecfuncT& Unemo) const;

	/// compute the reconstructed orbitals, and all potentials applied on nemo

	/// to use these potentials in the fock matrix computation they must
	/// be multiplied by the nuclear correlation factor
	/// @param[in]	nemo	the nemo orbitals
	/// @param[out]	psi		the reconstructed, full orbitals
	/// @param[out]	Jnemo	Coulomb operator applied on the nemos
	/// @param[out]	Knemo	exchange operator applied on the nemos
	/// @param[out]	Vnemo	nuclear potential applied on the nemos
	/// @param[out]	Unemo	regularized nuclear potential applied on the nemos
	void compute_nemo_potentials(const vecfuncT& nemo, vecfuncT& psi,
			vecfuncT& Jnemo, vecfuncT& Knemo, vecfuncT& Vnemo,
			vecfuncT& Unemo) const;

	/// normalize the nemos
	void normalize(vecfuncT& nemo) const;

    /// orthonormalize the vectors

    /// @param[in,out]	nemo	the vectors to be orthonormalized
    void orthonormalize(vecfuncT& nemo) const;

	/// return the Coulomb potential
	real_function_3d get_coulomb_potential(const vecfuncT& psi) const;

	bool is_dft() const {return calc->xc.is_dft();}

	/// localize the nemo orbitals
	vecfuncT localize(const vecfuncT& nemo) const;

	vecfuncT apply_exchange(const vecfuncT& nemo, const vecfuncT& psi) const;

	template<typename solverT>
	void rotate_subspace(World& world, const tensorT& U, solverT& solver,
			int lo, int nfunc, double trantol) const;

    /// save a function
    template<typename T, size_t NDIM>
    void save_function(const Function<T,NDIM>& f, const std::string name) const;

    /// load a function
    template<typename T, size_t NDIM>
    void load_function(Function<T,NDIM>& f, const std::string name) const;

    /// save a function
    template<typename T, size_t NDIM>
    void save_function(const std::vector<Function<T,NDIM> >& f, const std::string name) const;

    /// load a function
    template<typename T, size_t NDIM>
    void load_function(std::vector<Function<T,NDIM> >& f, const std::string name) const;

};

}

#endif /* NEMO_H_ */

