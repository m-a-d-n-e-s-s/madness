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
#include <chem/SCFProtocol.h>
#include <chem/correlationfactor.h>
#include <chem/molecular_optimizer.h>
#include <examples/nonlinsol.h>
#include <madness/mra/vmra.h>
#include <chem/pcm.h>

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

	void construct_nuclear_correlation_factor() {
		// construct the nuclear potential
		// Make the nuclear potential, initial orbitals, etc.
		calc->make_nuclear_potential(world);
		calc->potentialmanager->vnuclear().print_size("vnuc");
		calc->project_ao_basis(world);
//		save_function(calc->potentialmanager->vnuclear(),"vnuc");
	    // construct the nuclear correlation factor:
	    nuclear_correlation=create_nuclear_correlation_factor(world,*calc);
	    R = nuclear_correlation->function();
	    R.set_thresh(FunctionDefaults<3>::get_thresh());
	    R_square = nuclear_correlation->square();
	    R_square.set_thresh(FunctionDefaults<3>::get_thresh());
	}

	double value() {return value(calc->molecule.get_all_coords());}

	double value(const Tensor<double>& x);

	/// compute the nuclear gradients
	Tensor<double> gradient(const Tensor<double>& x);

	bool provides_gradient() const {return true;}

	/// returns the molecular hessian matrix at structure x
	Tensor<double> hessian(const Tensor<double>& x);

	/// purify and symmetrize the hessian

	/// The hessian should be symmetric, but it is not, because
	/// \f[
	///  \langle i^{Y_B}|H^{X_A}|i\rangle \neq \langle i|H^{X_A}|i^{Y_B}\rangle
	/// \f]
	/// does holds analytically, but not numerically. If the two numbers
	/// differ, pick the more trustworthy, which is the one with a heavy
	/// atom causing the perturbed density and the light atom being the
	/// nuclear singularity.
	/// @param[in]  hessian the raw hessian
	/// @return     a symmetrized hessian
	Tensor<double> purify_hessian(const Tensor<double>& hessian) const;


	/// solve the CPHF equations for the nuclear displacements

	/// this function computes that part of the orbital response that is
	/// orthogonal to the occupied space. If no NCF's are used this
	/// corresponds to the normal response. If NCF's are used the part
	/// parallel to the occupied space must be added!
	/// \f[
	///     F^X = F^\perp + F^\parallel
	/// \f]
	/// cf parallel_CPHF()
	/// @param[in]  iatom   the atom A to be moved
	/// @param[in]  iaxis   the coordinate X of iatom to be moved
	/// @return     \ket{i^X} or \ket{F^\perp}
	vecfuncT solve_cphf(const int iatom, const int iaxis, const Tensor<double> fock,
	        const vecfuncT& guess, const vecfuncT& rhsconst,
	        const Tensor<double> incomplete_hessian, const vecfuncT& parallel,
	        const SCFProtocol& p, const std::string& xc_data) const;

	/// solve the CPHF equation for all displacements

	/// this function computes the nemo response F^X
    /// \f[
    ///     F^X = F^\perp + F^\parallel
    /// \f]
	/// To reconstruct the unregularized orbital response (not recommended):
	/// \f[
	///   i^X   = R^X F + R F^X
	/// \f]
	/// The orbital response i^X constructed in this way is automatically
	/// orthogonal to the occupied space because of the parallel term F^\parallel
	/// @return a vector of the nemo response F^X for all displacements
	std::vector<vecfuncT> compute_all_cphf();

    /// this function computes that part of the orbital response that is
    /// parallel to the occupied space.
    /// \f[
    ///     F^X = F^\perp + F^\parallel
    /// \f]
    /// If no NCF's are used F^\parallel vanishes.
    /// If NCF's are used this term does not vanish because the derivatives of
    /// the NCF does not vanish, and it is given by
    /// \f[
    ///  F_i^\parallel = -\frac{1}{2}\sum_k|F_k ><F_k | (R^2)^X | F_i>
    /// \f]
    vecfuncT compute_cphf_parallel_term(const int iatom, const int iaxis) const;

    /// compute the IR intensities in the double harmonic approximation

    /// use the projected normal modes; units are km/mol
    /// @param[in]  normalmodes the normal modes
    /// @param[in]  dens_pt the perturbed densities for each nuclear displacement
    Tensor<double> compute_IR_intensities(const Tensor<double>& normalmodes,
            const vecfuncT& dens_pt) const;

	std::shared_ptr<SCF> get_calc() const {return calc;}

	PCM get_pcm()const{return pcm;}

	/// compute the Fock matrix from scratch
	tensorT compute_fock_matrix(const vecfuncT& nemo, const tensorT& occ) const;

	/// return a reference to the molecule
	Molecule& molecule() {return calc->molecule;}

    /// return a reference to the molecule
    Molecule& molecule() const {
        return calc->molecule;
    }

    /// make the density (alpha or beta)
    real_function_3d make_density(const Tensor<double>& occ,
            const vecfuncT& nemo) const;

    /// make the density using different bra and ket vectors

    /// e.g. for computing the perturbed density \sum_i \phi_i \phi_i^X
    /// or when using nemos: \sum_i R2nemo_i nemo_i
    real_function_3d make_density(const tensorT & occ,
            const vecfuncT& bra, const vecfuncT& ket, const bool refine=false) const;

    /// make the derivative of the density

    /// \f$ \nabla\rho = 2R^X R \rho_R + R^2\nabla \rho_R \f$
    /// @param[in]  rhonemo    the regularized density
    /// @param[in]  axis       the component of the nabla operator
    /// @return     the gradient of the *reconstructed* density
    real_function_3d make_ddensity(const real_function_3d& rhonemo,
            const int axis) const;

    /// compute the reduced densities sigma (gamma) for GGA functionals
    real_function_3d make_sigma(const real_function_3d& rho1,
            const real_function_3d& rho2) const;


    /// the Laplacian of the density

    /// The Laplacian should currently only be used for subsequent convolution
    /// with a Green's function (which is reasonably stable), but not on its own!
    ///
    /// The Laplacian of the cuspy density is numerically fairly unstable:
    ///  - a singular term may be rewritten using the nuclear potential (see below)
    ///  - the Laplacian of the regularized density is still very noisy
    ///
    /// It may be computed as
    /// \f[
    ///   \Delta \rho = \Delta (R^2 \rho_R)
    ///          = \Delta (R^2) \rho_R + 2\nabla R \nabla \rho_R + R^2 \Delta \rho_R
    ///          = 2 R^2 U1^2 \rho_R -4 R^2 ( U-V ) \rho_R + R^2 \Delta\rho_R
    /// \f]
    /// where we can use the identity
    /// \f[
    ///   U=V + R^{-1}[T,R]
    ///   -2 R (U-V) = \Delta R + 2\nabla R\dot \nabla
    /// \f]
    /// first term comes from the definition of the U potential as the commutator
    /// over the kinetic energy (aka the Laplacian)
    /// @param[in]  rhonemo    the regularized density \rho_R
    /// @return     the laplacian of the reconstructed density \Delta (R^2\rho_R)
    real_function_3d make_laplacian_density(const real_function_3d& rhonemo) const;

    /// compute the kinetic energy potential using Eq. (16) of
    /// R. A. King and N. C. Handy, “Kinetic energy functionals from the Kohn–Sham potential,”
    /// Phys. Chem. Chem. Phys., vol. 2, no. 22, pp. 5049–5056, 2000.
    real_function_3d kinetic_energy_potential(const vecfuncT& nemo) const;


    /// smooth a function by projecting it onto k-1 and then average with k

    /// kept it here for further testing
    static void smoothen(real_function_3d& f) {
        int k=f.get_impl()->get_k();
        real_function_3d fproj=project(f,k-1);
        real_function_3d freproj=project(fproj,k);
        f=0.5*(f+freproj);
    }

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

	struct timer {
        World& world;
	    double ttt,sss;
	    timer(World& world) : world(world) {
	        world.gop.fence();
	        ttt=wall_time();
	        sss=cpu_time();
	    }

	    void tag(const std::string msg) {
            world.gop.fence();
	        double tt1=wall_time()-ttt;
	        double ss1=cpu_time()-sss;
	        if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg.c_str(), ss1, tt1);
	        ttt=wall_time();
	        sss=cpu_time();
	    }

	    void end(const std::string msg) {
            world.gop.fence();
            double tt1=wall_time()-ttt;
            double ss1=cpu_time()-sss;
            if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg.c_str(), ss1, tt1);
        }

	};

public:

	/// the nuclear correlation factor
	std::shared_ptr<NuclearCorrelationFactor> nuclear_correlation;

	/// the nuclear correlation factor
	real_function_3d R;

    /// the square of the nuclear correlation factor
    real_function_3d R_square;

private:

	/// sum of square of coords at last solved geometry
	mutable double coords_sum;

	/// a poisson solver
	std::shared_ptr<real_convolution_3d> poisson;

	/// polarizable continuum model
	PCM pcm;

	void print_nuclear_corrfac() const;

	/// adapt the thresholds consistently to a common value
    void set_protocol(const double thresh) {

        calc->set_protocol<3>(world,thresh);

        // (re) construct nuclear potential and correlation factors
        construct_nuclear_correlation_factor();

        // (re) construct the Poisson solver
        poisson = std::shared_ptr<real_convolution_3d>(
                CoulombOperatorPtr(world, calc->param.lo, FunctionDefaults<3>::get_thresh()));

        // set thresholds for the MOs
        set_thresh(world,calc->amo,thresh);
        set_thresh(world,calc->bmo,thresh);

    }

	/// solve the HF equations
	double solve(const SCFProtocol& proto);

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
	/// @param[out]	pcmnemo	PCM (solvent) potential applied on the nemos
	/// @param[out]	Unemo	regularized nuclear potential applied on the nemos
	void compute_nemo_potentials(const vecfuncT& nemo, vecfuncT& psi,
			vecfuncT& Jnemo, vecfuncT& Knemo, vecfuncT& pcmnemo,
			vecfuncT& Unemo) const;

	/// normalize the nemos
	void normalize(vecfuncT& nemo) const;

    /// orthonormalize the vectors

    /// @param[in,out]	nemo	the vectors to be orthonormalized
    void orthonormalize(vecfuncT& nemo) const;

	/// return the Coulomb potential
	real_function_3d get_coulomb_potential(const vecfuncT& psi) const;

	/// compute the incomplete hessian

	/// incomplete hessian is the nuclear-nuclear contribution, and the
	/// contribution from the second derivative of the nuclear potential,
	/// and also the derivative of the nuclear correlation factor.
	/// i.e. all contributions that *do not* contain the regularized perturbed
	/// density, but it will contain parts of the perturbed density
	Tensor<double> make_incomplete_hessian() const;

    /// compute the complementary incomplete hessian

	/// @param[in]  xi the response functions including the parallel part
    Tensor<double> make_incomplete_hessian_response_part(
            const std::vector<vecfuncT>& xi) const;

	/// compute the constant term for the CPHF equations

	/// mainly all terms with the nuclear correlation factor's derivatives
	vecfuncT make_cphf_constant_term(const int iatom, const int iaxis,
	        const vecfuncT& R2nemo, const real_function_3d& rhonemo) const;

	public:

	bool is_dft() const {return calc->xc.is_dft();}

	bool do_pcm() const {return calc->param.pcm_data != "none";}
	
	private:

	/// localize the nemo orbitals
    vecfuncT localize(const vecfuncT& nemo, const double dconv, const bool randomize) const;

	/// return the threshold for vanishing elements in orbital rotations
    double trantol() const {
        return calc->vtol / std::min(30.0, double(get_calc()->amo.size()));
    }

	template<typename solverT>
	void rotate_subspace(World& world, const tensorT& U, solverT& solver,
			int lo, int nfunc) const;

    tensorT Q2(const tensorT& s) const {
        tensorT Q = -0.5*s;
        for (int i=0; i<s.dim(0); ++i) Q(i,i) += 1.5;
        return Q;
    }

    void make_plots(const real_function_3d &f,const std::string &name="function")const{
        double width = FunctionDefaults<3>::get_cell_min_width()/2.0 - 1.e-3;
        plot_plane(world,f,name);
        coord_3d start(0.0); start[0]=-width;
        coord_3d end(0.0); end[0]=width;
        plot_line(("line_"+name).c_str(),1000,start,end,f);
    }

    /// save a function
    template<typename T, size_t NDIM>
    void save_function(const std::vector<Function<T,NDIM> >& f, const std::string name) const;

    /// load a function
    template<typename T, size_t NDIM>
    void load_function(std::vector<Function<T,NDIM> >& f, const std::string name) const;

};

/// rotate the KAIN subspace (cf. SCF.cc)
template<typename solverT>
void Nemo::rotate_subspace(World& world, const tensorT& U, solverT& solver,
        int lo, int nfunc) const {
    std::vector < vecfunc<double, 3> > &ulist = solver.get_ulist();
    std::vector < vecfunc<double, 3> > &rlist = solver.get_rlist();
    for (unsigned int iter = 0; iter < ulist.size(); ++iter) {
        vecfuncT& v = ulist[iter].x;
        vecfuncT& r = rlist[iter].x;
        vecfuncT vnew = transform(world, vecfuncT(&v[lo], &v[lo + nfunc]), U,
                trantol(), false);
        vecfuncT rnew = transform(world, vecfuncT(&r[lo], &r[lo + nfunc]), U,
                trantol(), true);

        world.gop.fence();
        for (int i=0; i<nfunc; i++) {
            v[i] = vnew[i];
            r[i] = rnew[i];
        }
    }
    world.gop.fence();
}

/// save a function
template<typename T, size_t NDIM>
void Nemo::save_function(const std::vector<Function<T,NDIM> >& f, const std::string name) const {
    if (world.rank()==0) print("saving vector of functions",name);
    archive::ParallelOutputArchive ar(world, name.c_str(), 1);
    ar & f.size();
    for (const Function<T,NDIM>& ff:f)  ar & ff;
}

/// save a function
template<typename T, size_t NDIM>
void Nemo::load_function(std::vector<Function<T,NDIM> >& f, const std::string name) const {
    if (world.rank()==0) print("loading vector of functions",name);
    archive::ParallelInputArchive ar(world, name.c_str(), 1);
    std::size_t fsize=0;
    ar & fsize;
    f.resize(fsize);
    for (std::size_t i=0; i<fsize; ++i) ar & f[i];
}

}

#endif /* NEMO_H_ */

