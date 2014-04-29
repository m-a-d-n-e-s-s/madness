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

#define WORLD_INSTANTIATE_STATIC_TEMPLATES


/*!
  \file examples/tdhf.cc
  \brief compute the time-dependent HF equations (currently CIS approximation)

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
  /trunk/src/apps/examples/tdhf.cc>here</a>.

*/

#include <examples/projector.h>
//#include <examples/mp2.h>

#include<examples/nonlinsol.h>
#include<moldft/moldft.h>
#include <mra/operator.h>
#include <mra/mra.h>
#include <mra/vmra.h>
#include <mra/lbdeux.h>

#include<iomanip>
#include<iostream>

using namespace madness;

const static double bsh_eps=1.e-6;

namespace madness {

typedef std::vector<Function<double,3> > vecfuncT;

// This class is used to store information for the non-linear solver
struct F {
    World& world;
    vecfuncT x;

    F(World& world, const vecfuncT& x1) : world(world), x(x1) {}

    F(const F& other) : world(other.world), x(other.x) {}

    F& operator=(const F& other) {
    	x=other.x;
    	return *this;
    }

    F operator-(const F& b) {
        return F(world,sub(world,x,b.x));
    }

    F operator+=(const F& b) { // Operator+= necessary
    	x=add(world,x,b.x);
    	return *this;
    }

    F operator*(double a) { // Scale by a constant necessary
    	scale(world,x,a);
        return *this;
    }
};

/// the non-linear solver requires an inner product
double inner(const F& a, const F& b) {
    Tensor<double> i=inner(a.world,a.x,b.x);
    return i.sum();
}

// The default constructor for functions does not initialize
// them to any value, but the solver needs functions initialized
// to zero for which we also need the world object.
struct allocator {
    World& world;
    const int n;

    /// @param[in]	world	the world
    /// @param[in]	nn		the number of functions in a given vector
    allocator(World& world, const int nn) : world(world), n(nn) {}

    /// allocate a vector of n empty functions
    F operator()() {
        return F(world,zero_functions<double,3>(world,n));
    }
};


/// POD holding excitation energy and response vector for a single excitation
struct root {
	root() : omega(0.0) {}
	root(vecfuncT& x, double omega) : x(x), omega(omega) {}
	vecfuncT x;
	double omega;
	std::vector<double> amplitudes_;
};

/// for convenience
double inner(const root& a, const root& b) {
	if (a.x.size()==0) return 0.0;
	return inner(a.x[0].world(),a.x,b.x).sum();
}

double rfunction(const coord_3d& r) {
    return sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
}

/// helper struct for computing the moments
struct xyz {
	int direction;
	xyz(int direction) : direction(direction) {}
	double operator()(const coord_3d& r) const {
		return r[direction];
	}
};

/// The CIS class holds all machinery to compute excited state properties
class CIS {
	typedef SeparatedConvolution<double,3> operatorT;
	typedef std::shared_ptr<operatorT> poperatorT;
	typedef std::shared_ptr< FunctionFunctorInterface<double,3> > functorT;

public:

	/// ctor

	/// @param[in]	world	the world
	/// @param[in]	hf		the HartreeFock reference state
	/// @param[in]	input	the input file name
	CIS(World& world, const Calculation& calc, const std::string input)
		: world(world), calc_(calc), guess_("all_orbitals"), nroot_(5)
		, nfreeze_(0), econv_(calc.param.econv)
		, dconv_(calc.param.dconv), print_grid_(false),
		fixed_point_(false) {


		omega_=std::vector<double>(9,100.0);

        std::ifstream f(input.c_str());
        position_stream(f, "CIS");
        std::string s, tag;
        while (std::getline(f,s)) {
            std::istringstream ss(s);
            ss >> tag;
            if (tag == "end") break;
            else if (tag == "guess") ss >> guess_;
            else if (tag == "nroot") ss >> nroot_;
            else if (tag == "freeze") ss >> nfreeze_;
            else if (tag == "econv") ss >> econv_;
            else if (tag == "dconv") ss >> dconv_;
            else if (tag == "fixed_point") fixed_point_=true;
            else if (tag == "print_grid") print_grid_=true;
            else if (tag == "omega0") ss >> omega_[0];
            else if (tag == "omega1") ss >> omega_[1];
            else if (tag == "omega2") ss >> omega_[2];
            else if (tag == "omega3") ss >> omega_[3];
            else if (tag == "omega4") ss >> omega_[4];
            else if (tag == "omega5") ss >> omega_[5];
            else if (tag == "omega6") ss >> omega_[6];
            else if (tag == "omega7") ss >> omega_[7];
            else if (tag == "omega8") ss >> omega_[8];
            else continue;
        }

        if (world.rank() == 0) {
            madness::print("\n ======= CIS info =======\n");
    		if (nfreeze_==0) madness::print("   # frozen orbitals ","none");
    		if (nfreeze_>0) madness::print("   # frozen orbitals ",0, " to ",nfreeze_-1);
			madness::print("     active orbitals ", nfreeze_," to ",get_calc().param.nalpha-1);

            madness::print("          guess from ", guess_);
            madness::print("        threshold 3D ", FunctionDefaults<3>::get_thresh());
            madness::print("  energy convergence ", econv_);
            madness::print("max residual (dconv) ", dconv_);
            madness::print("     number of roots ", nroot_);
            madness::print(" omega ", omega_[0],omega_[1],omega_[2]);

        }

        lo=get_calc().param.lo;
	}

	/// return the HF reference
	const Calculation& get_calc() const {return calc_;}

	/// solve the CIS equations for n roots
	void solve() {

		// for convenience
	    const vecfuncT& amo=get_calc().amo;
	    const std::size_t nmo=amo.size();
		const int noct=nmo-nfreeze_;	// # of active orbitals

	    // this is common for all roots in all iterations
	    exchange_intermediate_=make_exchange_intermediate(active_mo(),active_mo());

    	// guess the amplitudes for all roots
	    for (std::size_t iroot=0; iroot<nroot_; ++iroot) {
	    	root root1=guess_amplitudes(iroot);
	        roots().push_back(root1);
	    }

	    // one KAIN solver for each root
	    typedef  XNonlinearSolver<F,double,allocator> solverT;
	    solverT onesolver(allocator(world,noct));
	    onesolver.set_maxsub(3);
	    std::vector<solverT> solver(nroot_,onesolver);

	    bool converged=iterate_all_CIS_roots(world,solver,roots());
	    if (converged) {
	    	analyze(roots());
	    	if (world.rank()==0) print(" CIS iterations converged ");
	    } else {
	    	if (world.rank()==0) print(" CIS iterations not converged ");
	    }
	}

	/// return the roots of the response equation
	std::vector<root>& roots() {return roots_;}

	/// are we supposed to print the grid for an external guess
	bool print_grid() const {return print_grid_;}

private:

	/// the world
	World& world;

	/// the HartreeFock reference state
	const Calculation& calc_;

	/// the excited states aka the roots of the response equation
    std::vector<root> roots_;

    /// intermediate for the two-electron interaction term, Eq. (8)

	/// the intermediate is the same for all roots:
	/// \[
	///   int[p,i] = \int 1/r12 \phi_i(1) * \phi_p(1)
	/// \]
    /// with p \in noct, i \in nocc
    std::vector<vecfuncT> exchange_intermediate_;

    /// the coulomb potential
    mutable real_function_3d coulomb_;

    /// where we get our guess from (all_virtual, koala)
    std::string guess_;

    /// the phases of the guess and ours might differ
    Tensor<double> guess_phases_;

    /// number of roots we are supposed to solve
    int nroot_;

    /// number of frozen orbitals
    int nfreeze_;

    /// guess for the excitation energies
    std::vector<double> omega_;

    /// energy convergence threshold
    double econv_;

    /// density convergence threshold (=residual)
    double dconv_;

    /// flag if the grid for the density should be printed

    /// external programs (e.g. Koala) need this grid to export the guess roots
    bool print_grid_;

    /// perform a fixed-point iteration of the given excitation energies
    bool fixed_point_;

    double lo;

    /// guess amplitudes

    /// note that the orbitals from an external guess might be rotated wrt the
    /// MRA orbitals, which would spoil the guess. Therefore we rotate the
    /// guess with the overlap matrix so that the guess conforms with the MRA
    /// phases. The overlap to the external orbitals (if applicable) is computed
    /// only once
    /// @param[in]	iroot	guess for root iroot
    root guess_amplitudes(const int iroot) {

    	// for convenience
    	const std::size_t nmo=get_calc().amo.size();	// all orbitals
    	const int noct=nmo-nfreeze_;						// active orbitals

    	// default empty root
    	root root;
		root.amplitudes_=std::vector<double>(nmo,1.0);

        // guess an excitation energy: 0.9* HOMO or use input from file
        root.omega=-0.9*get_calc().aeps(nmo-1);
    	if (omega_[iroot]!=100.0) root.omega=omega_[iroot];

    	// check if there's a root on disk, if so return those amplitudes
    	root.x.resize(noct);
    	if (load_root(world,iroot,root)) return root;
	    root.x.clear();

	    // get the guess from koala
	    if (guess_=="koala") {

			vecfuncT koala_amo;

			// first we need to determine the rotation of the external orbitals
	    	// to the MRA orbitals, so that we can subsequently rotate the guess
	    	if (not guess_phases_.has_data()) {

				// read koala's orbitals from disk
				for (std::size_t i=nfreeze_; i<nmo; ++i) {
					real_function_3d x_i=real_factory_3d(world).empty();
					const std::string valuefile="grid.koala.orbital"+stringify(i);
					x_i.get_impl()->read_grid2<3>(valuefile,functorT());
					koala_amo.push_back(x_i);
				}
				// this is the transformation matrix for the rotation
				guess_phases_=matrix_inner(world,koala_amo,get_calc().amo);
				guess_phases_=guess_phases_(_,Slice(nfreeze_,nmo-1));

	    	}

	    	// read the actual external guess from file
	    	for (std::size_t i=nfreeze_; i<nmo; ++i) {

    	    	// this is the file where the guess is on disk
        	    const std::string valuefile="grid.koala.orbital"+stringify(i)
        	    		+".excitation"+stringify(iroot);
        		real_function_3d x_i=real_factory_3d(world).empty();
        		x_i.get_impl()->read_grid2<3>(valuefile,functorT());
    	    	root.x.push_back(x_i);
    	    }

	    	// compute the inverse of the overlap matrix
	    	Tensor<double> S=(guess_phases_+transpose(guess_phases_)).scale(0.5);
	    	Tensor<double> U, eval;
	    	syev(S,U,eval);
	    	Tensor<double> Sinv=copy(U);
	    	for (int i=0; i<U.dim(0); ++i) {
	    		for (int j=0; j<U.dim(1); ++j) {
	    			Sinv(i,j)/=eval(j);
	    		}
	    	}
	    	Sinv=inner(Sinv,U,-1,-1);

	    	// now rotate the active orbitals from the guess to conform with
    	    // the MRA orbitals
//	    	Sinv=Sinv(Slice(nfreeze_,nmo-1),Slice(nfreeze_,nmo-1));
    	    root.x=transform(world,root.x,Sinv);
    	    save_root(world,iroot,root);


    	} else if (guess_=="all_orbitals") {
    		// Take a linear combination of all orbitals as guess, because
    		// there are not enough virtuals in the minimal basis set for all
    		// possible symmetries
    		if (world.rank()==0) {
    			print("taking as guess all orbitals");
    		}
    		real_function_3d all_orbitals=real_factory_3d(world);
    		for (std::size_t ivir=0; ivir<get_calc().ao.size(); ++ivir) {
    			all_orbitals+=get_calc().ao[ivir];
			}

    		// construct the projector on the occupied space
    		const vecfuncT& amo=get_calc().amo;
    		Projector<double,3> rho0(amo);

    	    real_function_3d r=real_factory_3d(world).f(rfunction);

    		// multiply the guess with r and project out the occupied space
    		real_function_3d r_all_orbitals=all_orbitals*r;
    		r_all_orbitals-=rho0(r_all_orbitals);

    		for (std::size_t iocc=nfreeze_; iocc<nmo; ++iocc) {
    			root.x.push_back(copy(r_all_orbitals));
    		}


    	} else {
          	MADNESS_EXCEPTION("unknown source to guess CIS amplitudes",1);
    	}
    	MADNESS_ASSERT(root.x.size()==noct);

    	// normalize this root
    	normalize(world,root);

    	return root;
    }


    /// solve the CIS equations for all roots

    /// @param[in]	world	the world
    /// @param[in]	solver	the KAIN solver (unused right now)
    /// @param[inout]	roots	on entry: guess for the roots
    ///                         on successful exit: converged root
    /// @return	convergence reached or not
	template<typename solverT>
    bool iterate_all_CIS_roots(World& world, std::vector<solverT>& solver,
    		std::vector<root>& roots) const {

		// check orthogonality of the roots
		orthonormalize(world,roots);
		Tensor<double> ovlp=overlap(roots,roots);
		for (int i=0; i<ovlp.dim(0); ++i) ovlp(i,i)-=1.0;
		if (ovlp.normf()/ovlp.size()>econv_) {
			print(ovlp);
			MADNESS_EXCEPTION("non-orthogonal roots",1);
		}

		// construct the projector on the occupied space
		const vecfuncT& amo=get_calc().amo;
	    vecfuncT act=this->active_mo();
		const std::size_t noct=act.size();
		Projector<double,3> rho0(amo);

		// start iterations
		for (int iter=0; iter<50; ++iter) {

			Tensor<double> error(roots.size());
			// print progress on the computation
			if (world.rank()==0) {
				print("starting iteration ", iter, " at time ", wall_time());
				print(" root   excitation energy   energy correction         error");
			}

			// apply the BSH operator on each root
			START_TIMER(world);
			for (std::size_t iroot=0; iroot<roots.size(); ++iroot) {
				std::cout << std::setw(4) << iroot << " ";
				error(iroot)=iterate_one_CIS_root(world,solver[iroot],roots[iroot]);
			}
			orthonormalize(world,roots);
			END_TIMER(world,"BSH step");


			// compute the Fock matrix elements
			Tensor<double> Fock_pt=make_perturbed_fock_matrix(roots,act,rho0);

			// add the term with the orbital energies:
			//  F_pq -= \sum_i\epsilon_i <x^p_i | x^q_i>
			for (std::size_t i=0; i<roots.size(); ++i) {
				for (std::size_t j=0; j<roots.size(); ++j) {
					Tensor<double> eij=inner(world,roots[i].x,roots[j].x);
					for (int ii=0; ii<noct; ++ii) {
						Fock_pt(i,j)-=get_calc().aeps[ii+nfreeze_]*eij[ii];
					}
				}
			}

			// diagonalize the roots in their subspace
//			print("perturbed Fock matrix ");
//			print(Fock_pt);
			Tensor<double> U, evals;
			syev(Fock_pt,U,evals);
			print("eval(F)",evals);
//			print("U");
//			print(U);

			// check energy and density convergence (max residual norm)
			bool dconverged=(error.max()<dconv_);
			bool econverged=true;
			for (std::size_t i=0; i<roots.size(); ++i) {
				if (std::fabs(roots[i].omega - evals(i))>econv_) econverged=false;
			}
			if (econverged and dconverged) return true;

			// diagonalize the amplitudes in the space of the perturbed Fock
			// matrix
	        std::vector<vecfuncT> vc(nroot_);
			for (std::size_t iroot=0; iroot<nroot_; ++iroot) {
				vc[iroot]=zero_functions<double,3>(world,noct);
		        compress(world, vc[iroot]);
		        compress(world,roots[iroot].x);
			}

	        for (int i=0; i<roots.size(); ++i) {
	            for (int j=0; j<roots.size(); ++j) {
	               	gaxpy(world,1.0,vc[i],U(j,i),roots[j].x);
	            }
	        }

			for (std::size_t iroot=0; iroot<roots.size(); ++iroot) {
		        roots[iroot].x=vc[iroot];
		        normalize(world,roots[iroot]);
		        roots[iroot].omega=evals[iroot];
			}

			for (std::size_t i=0; i<roots.size(); ++i) save_root(world,i,roots[i]);
		}
		return false;
    }

    /// iterate the TDHF or CIS equations

	/// follow Eq (4) of
	/// T. Yanai, R. J. Harrison, and N. Handy,
	/// ÒMultiresolution quantum chemistry in multiwavelet bases: time-dependent
	/// density functional theory with asymptotically corrected potentials in
	/// local density and generalized gradient approximations,Ó
	/// Mol. Phys., vol. 103, no. 2, pp. 413Ð424, 2005.
	///
	/// The convergence criterion is that the excitation amplitudes don't change
    /// @param[in]		world	the world
    /// @param[in]		solver	the KAIN solver (not used right now..)
    /// @param[inout]	root	the current root that we solve
    /// @return			the residual error
	template<typename solverT>
	double iterate_one_CIS_root(World& world, solverT& solver, root& thisroot) const {

		// for convenience
		const vecfuncT& amo=get_calc().amo;
		const int nmo=amo.size();		// # of orbitals in the HF calculation
		const int noct=nmo-nfreeze_;	// # of active orbitals

		vecfuncT& x=thisroot.x;
		double& omega=thisroot.omega;

	    // these are the active orbitals in a vector (shallow-copied)
	    vecfuncT active_mo=this->active_mo();

		// the projector on the unperturbed density
		Projector<double,3> rho0(amo);

		// apply the unperturbed potential
		const vecfuncT V0=apply_fock_potential(x);

		// apply the gamma potential on the active MOs; Eq. (4-6)
		const vecfuncT Gamma=apply_gamma(x,active_mo,rho0);

		// add the local potential and the electron interaction term
		vecfuncT Vphi=add(world,V0,Gamma);
		scale(world,Vphi,-2.0);
		truncate(world,Vphi);

		// the bound-state helmholtz function for omega < orbital energy
	    std::vector<poperatorT> bsh(noct);
	    for(int p = 0; p<noct; ++p){
	        double eps = get_calc().aeps[p+nfreeze_] + omega;	// orbital energy
	        if(eps > 0){
	            if(world.rank() == 0)
	            	print("bsh: warning: positive eigenvalue", p+nfreeze_, eps);
	            eps = -0.03;
	        }
	        bsh[p] = poperatorT(BSHOperatorPtr3D(world, sqrt(-2.0 * eps), lo, bsh_eps));
	    }
	    const vecfuncT GVphi=apply(world,bsh,Vphi);

	    // compute the residual aka the difference between to old and the new solution vector
	    const vecfuncT residual=sub(world,x,GVphi);

		// update the excitation energy omega
		Tensor<double> t1=inner(world,Vphi,residual);
		double t2=0.0;
		for (int i=0; i<noct; ++i) {
			double n=GVphi[i].norm2();
			t2+=n*n;
		}
		// remove factor 2 from Vphi coming from the BSH application
		double delta=0.5*t1.sum()/t2;

		const double error=sqrt(inner(F(world,residual),F(world,residual)));

		// some update on the progress for the user
		if (world.rank()==0) {
			std::cout << std::scientific << std::setprecision(10);
			std::cout << std::setw(20) << omega	<< std::setw(20)<< delta
					<< std::setw(19) << error << std::endl;
        }

		// update the energy
		if (not fixed_point_) omega+=delta;

	    // update the x vector: orthogonalize against the occupied space
	    if (std::fabs(delta)>1.e-3) {
	    	x=GVphi;
	    } else {
	    	F ff=solver.update(F(world,x),F(world,residual));
	    	x=ff.x;
	    }

		x=GVphi;
		for (int p=0; p<noct; ++p) x[p] -= rho0(GVphi[p]);

		return error;
	}

	/// iterate the TDHF or CIS equations -- deprecated !!

	/// follow Eq (4) of
	/// T. Yanai, R. J. Harrison, and N. Handy,
	/// ÒMultiresolution quantum chemistry in multiwavelet bases: time-dependent
	/// density functional theory with asymptotically corrected potentials in
	/// local density and generalized gradient approximations,Ó
	/// Mol. Phys., vol. 103, no. 2, pp. 413Ð424, 2005.
	///
	/// The convergence criterion is that the excitation amplitudes don't change
    /// @param[in]		world	the world
    /// @param[in]		solver	the KAIN solver
    /// @param[inout]	root	the current root that we solve
	/// @param[in]		excited_states all lower-lying excited states for orthog.
	/// @param[in]		iteration	the current iteration
    /// @return			if the iteration on this root has converged
	template<typename solverT>
	bool iterate_CIS(World& world, solverT& solver, root& thisroot,
			const std::vector<root >& excited_states, const int iteration) const {

		// for convenience
		const vecfuncT& amo=get_calc().amo;
		const int nmo=amo.size();		// # of orbitals in the HF calculation
		const int noct=nmo-nfreeze_;	// # of active orbitals

		vecfuncT& x=thisroot.x;
		double& omega=thisroot.omega;

	    // these are the active orbitals in a vector (shallow-copied)
	    vecfuncT active_mo=this->active_mo();

		// the projector on the unperturbed density
		Projector<double,3> rho0(amo);

		// apply the unperturbed potential
		const vecfuncT V0=apply_fock_potential(x);

		// apply the gamma potential on the active MOs; Eq. (4-6)
		const vecfuncT Gamma=apply_gamma(x,active_mo,rho0);

		// add the local potential and the electron interaction term
		vecfuncT Vphi=add(world,V0,Gamma);
		scale(world,Vphi,-2.0);
		truncate(world,Vphi);

		// the bound-state helmholtz function for omega < orbital energy
	    std::vector<poperatorT> bsh(noct);
	    for(int p = 0; p<noct; ++p){
	        double eps = get_calc().aeps[p+nfreeze_] + omega;
	        if(eps > 0){
	            if(world.rank() == 0)
	            	print("bsh: warning: positive eigenvalue", p+nfreeze_, eps);
	            eps = -0.03;
	        }
	        bsh[p] = poperatorT(BSHOperatorPtr3D(world, sqrt(-2.0 * eps), lo, bsh_eps));
	    }
	    vecfuncT GVphi=apply(world,bsh,Vphi);

	    // compute the residual aka the difference between to old and the new solution vector
	    vecfuncT residual=sub(world,x,GVphi);

		// update the excitation energy omega
		Tensor<double> t1=inner(world,Vphi,residual);
		double t2=0.0;
		for (int i=0; i<noct; ++i) {
			double n=GVphi[i].norm2();
			t2+=n*n;
		}
		// remove factor 2 from Vphi coming from the BSH application
		double delta=0.5*t1.sum()/t2;

		// update the energy; damp the first few iterations
		if (iteration<3) ;
		else if (iteration<5) omega+=0.3*delta;
		else if (iteration<10) omega+=0.7*delta;
		else omega+=delta;

		// some update on the progress for the user
		const double error=sqrt(inner(F(world,residual),F(world,residual)));
        if (world.rank()==0) {
			std::cout << std::setw(6) << iteration << "    ";
			std::cout << std::scientific << std::setprecision(10);
			std::cout << std::setw(20) << omega	<< std::setw(20)<< delta
					<< std::setw(19) << error << " " ;
			std::cout.width(8); std::cout.precision(1);
			std::cout << std::fixed << wall_time()<< std::endl;
        }

		// generate a new trial solution vector
	    if (iteration<50) {
	    	x=GVphi;
	    } else {
	    	F ff=solver.update(F(world,x),F(world,residual));
	    	x=ff.x;
	    }

	    // update the x vector: orthogonalize against the occupied space
		for (int p=0; p<noct; ++p) x[p] -= rho0(x[p]);

		// orthogonalize the roots wrt each other
		for (std::size_t r=0; r<excited_states.size(); ++r) {

			const vecfuncT& lower=excited_states[r].x;
			// skip self
			if (&lower==&x) continue;

			Tensor<double> ovlp=inner(world,lower,x);

        	compress(world,lower,false);
        	compress(world,x,true);

            for (unsigned int p=0; p<x.size(); ++p) {
                x[p].gaxpy(1.0, lower[p], -ovlp(p), false);
            }
            world.gop.fence();
		}

		// normalize the transition density for all occupied orbitals
		root dummy(x,omega);
		normalize(world,dummy);
		truncate(world,x);
		const std::vector<double> norms=norm2s(world,x);

		// check energy, residual, and individual amplitude convergence
		bool converged=true;
		if (error>dconv_) converged=false;
		if (std::fabs(delta)>econv_) converged=false;

		// this is not converged if the amplitudes change by a lot (> 1 percent)
		for (std::size_t i=0; i<norms.size(); ++i) {
			const double change=std::fabs(thisroot.amplitudes_[i]/norms[i]-1.0);
			if (change>0.01) converged=false;
		}
		thisroot.amplitudes_=norms;

		// if convergence is reached print out weights and leave
		if (converged) {
			if (world.rank()==0) {
				std::cout.width(10); std::cout.precision(6);
				print("\nconverged:  ",omega);
		    	std::cout << std::fixed;
				for (int p=nfreeze_; p<nmo; ++p) {
					std::cout << "norm(x_"<<p<<") **2  "
							<< thisroot.amplitudes_[p-nfreeze_] << std::endl;
				}
			}
			return true;
		}
		return false;
	}

	/// apply the gamma potential of Eq. (6) on the x vector

	/// @param[in]	x		the response amplitudes
	/// @param[in]	act		the active orbitals p
	/// @param[in]	rho0	the projector on all (frozen & active) MOs
	/// @return		(1-\rho^0) \Gamma \phi_p
	vecfuncT apply_gamma(const vecfuncT& x, const vecfuncT& act,
			const Projector<double,3>& rho0) const {

		// for convenience
		const std::size_t noct=act.size();

		// now construct the two-electron contribution Gamma, cf Eqs. (7,8)
		real_function_3d rhoprime=real_factory_3d(world);

		// a poisson solver
	    std::shared_ptr<real_convolution_3d> poisson
	    	=std::shared_ptr<real_convolution_3d>
	    	(CoulombOperatorPtr(world,lo,bsh_eps));

		// the Coulomb part Eq. (7) and Eq. (3)
		for (int i=0; i<noct; ++i) rhoprime+=x[i]*act[i];
		real_function_3d pp=2.0*(*poisson)(rhoprime);

		// the Coulomb part Eq. (7) and Eq. (3)
		vecfuncT Gamma=mul(world,pp,act);

		// the exchange part Eq. (8)
		for (std::size_t p=0; p<noct; ++p) {

			// this is x_i * \int 1/r12 \phi_i \phi_p
			vecfuncT x_Ppi=mul(world,x,exchange_intermediate_[p]);
			for (std::size_t i=0; i<noct; ++i) Gamma[p]-=x_Ppi[i];

			// project out the zeroth-order density Eq. (4)
			Gamma[p]-=rho0(Gamma[p]);
		}
		return Gamma;
	}

	/// make the zeroth-order potential J-K+V_nuc

	/// note the use of all (frozen & active) orbitals in the computation
	/// @param[in]	x	the response vector
	/// @return		(J-K+V_nuc) x
	vecfuncT apply_fock_potential(const vecfuncT& x) const {

	    // the local potential V^0 of Eq. (4)
		real_function_3d coulomb;
	    real_function_3d vlocal = get_calc().potentialmanager->vnuclear() +
	    		get_coulomb_potential();

	    // make the potential for V0*xp
		vecfuncT Vx=mul(world,vlocal,x);

		// and the exchange potential is K xp
		vecfuncT Kx=get_calc().apply_hf_exchange(world,get_calc().aocc,
				get_calc().amo,x);

		// sum up: V0 xp = V_loc xp - K xp
		vecfuncT V0=sub(world,Vx,Kx);

		return V0;

	}

    /// return the Coulomb potential
    real_function_3d get_coulomb_potential() const {
        MADNESS_ASSERT(get_calc().param.spin_restricted);
        if (coulomb_.is_initialized()) return copy(coulomb_);
        functionT rho = get_calc().make_density(world, get_calc().aocc,
        		get_calc().amo).scale(2.0);
        coulomb_=get_calc().make_coulomb_potential(rho);
        return copy(coulomb_);
    }

	/// make the 2-electron interaction intermediate

	/// the intermediate is the same for all roots:
	/// \f[
	///   Int[i,p] = \int \frac{1}{r_{12}} \phi_i(1) * \phi_p(1)
	/// \f]
	/// both i and p are active MOs
	/// @param[in]	active_mo	active orbitals in the CIS computation
	/// @param[in]	amo			all MOs of the HF calculation
	/// @return		a vector of vectors of functions: [noct][nocc]
	std::vector<vecfuncT> make_exchange_intermediate(const vecfuncT& active_mo,
			const vecfuncT& amo) const {
		// a poisson solver
	    std::shared_ptr<real_convolution_3d> poisson
	    	=std::shared_ptr<real_convolution_3d>
	    	(CoulombOperatorPtr(world,lo,bsh_eps));

	    std::vector<vecfuncT> intermediate(amo.size());

	    for (std::size_t p=0; p<active_mo.size(); ++p) {
			intermediate[p]=apply(world,(*poisson),mul(world,active_mo[p],amo));
	    }
	    return intermediate;
	}

	/// compute the perturbed fock matrix

	/// the matrix is given by
	/// \f[
	///   F_{pq} = < x^p_i | F | x^q_i > + < x^p_i | \Gamma^q | \phi_i >
	/// \f]
	/// where an amplitude is given by its components x_i, such that and
	/// summation over the occupied orbitals i is implied
	/// \f[
	///   < \vec x^p | \vec x^q > = \sum_i x^p_i x^q_i
	/// \f]
	/// and similar for the Fock matrix elements.
	/// Note this is NOT the response matrix A
	Tensor<double> make_perturbed_fock_matrix(const std::vector<root>& roots,
			const vecfuncT& act, const Projector<double,3>& rho0) const {

		const std::size_t nroot=roots.size();
		Tensor<double> Fock_pt(nroot,nroot);

		START_TIMER(world);
		// apply the unperturbed Fock operator and the Gamma potential on all
		// components of all roots
		for (std::size_t iroot=0; iroot<nroot; ++iroot) {

			const vecfuncT& xp=roots[iroot].x;
			// apply the unperturbed potential and the Gamma potential
			// on the response amplitude x^p
			const vecfuncT V0=apply_fock_potential(xp);
			const vecfuncT Gamma=apply_gamma(xp,act,rho0);

			const vecfuncT V=add(world,V0,Gamma);

			// compute the matrix element V_pq
			// <p | V | q> = \sum_i <x^p_i | V | x^q_i>
			for (std::size_t jroot=0; jroot<nroot; ++jroot) {
				const vecfuncT& xq=roots[jroot].x;
				Tensor<double> xpi_Vxqi=inner(world,xq,V);
				Fock_pt(iroot,jroot)=xpi_Vxqi.sum();
			}
		}
		END_TIMER(world,"Fock_pt, potential");
		START_TIMER(world);

		// add kinetic part
    	for (std::size_t iroot=0; iroot<nroot; ++iroot) {
			const vecfuncT& xp=roots[iroot].x;
	        reconstruct(world, xp);
    	}

    	std::vector< std::shared_ptr<real_derivative_3d> > gradop;
        gradop = gradient_operator<double,3>(world);

		for(int axis = 0;axis < 3;++axis) {

			std::vector<vecfuncT> dxp;
        	for (std::size_t iroot=0; iroot<nroot; ++iroot) {

        		const vecfuncT& xp=roots[iroot].x;
	            const vecfuncT d = apply(world, *(gradop[axis]), xp);
	            dxp.push_back(d);
	        }
        	for (std::size_t iroot=0; iroot<nroot; ++iroot) {
            	for (std::size_t jroot=0; jroot<nroot; ++jroot) {
            		Tensor<double> xpi_Txqi=inner(world,dxp[iroot],dxp[jroot]);
            		Fock_pt(iroot,jroot)+=0.5*xpi_Txqi.sum();
            	}
        	}
		}
		END_TIMER(world,"Fock_pt, kinetic");
		return Fock_pt;
	}

	/// load a converged root from disk

	/// @param[in]	world 	the world
	/// @param[in]	iroot	the i-th root
	/// @param[inout]	x	the x-vector for the i-th root
	/// @param[out]	omega	the excitation energy
	/// @return	successfully loaded a root or not
	bool load_root(World& world, const int i, root& root)  const {
		std::string filename="root_"+stringify(i);
		bool exists=archive::ParallelInputArchive::exists(world,filename.c_str());
		if (not exists) return false;

		archive::ParallelInputArchive ar(world, filename.c_str(), 1);
		ar & root.omega;
		for (std::size_t i=0; i<root.x.size(); ++i) ar & root.x[i];
		return true;
	}

	/// save a converged root to disk

	/// @param[in]	world 	the world
	/// @param[in]	iroot	the i-th root
	/// @param[inout]	x	the x-vector for the i-th root
	/// @param[in]	omega	the excitation energy
	void save_root(World& world, const int i, const root& root) const {
		std::string filename="root_"+stringify(i);
		archive::ParallelOutputArchive ar(world, filename.c_str(), 1);
		ar & root.omega;
		for (std::size_t i=0; i<root.x.size(); ++i) ar & root.x[i];
	}

	/// normalize the excitation amplitudes

	/// normalize the set of excitation amplitudes, such that the sum of square
	/// of all amplitudes equals 1.
	/// @param[in]		world the world
	/// @param[inout]	x	the excitation vector
	void normalize(World& world, root& x) const {
		const double n2=inner(x,x);
		scale(world,x.x,1.0/sqrt(n2));
	}

	/// orthonormalize all roots using Gram-Schmidt
	void orthonormalize(World& world, std::vector<root>& roots) const {

		// first normalize
		for (std::size_t r=0; r<roots.size(); ++r) {
			normalize(world,roots[r]);
		}

		// orthogonalize the roots wrt each other
		for (std::size_t r=0; r<roots.size(); ++r) {
			vecfuncT& x=roots[r].x;
			for (std::size_t rr=0; rr<r; ++rr) {
				const vecfuncT& lower=roots[rr].x;

				const double ovlp=inner(world,lower,x).sum();
				compress(world,lower,false);
				compress(world,x,true);

				for (unsigned int p=0; p<x.size(); ++p) {
					x[p].gaxpy(1.0, lower[p], -ovlp, false);
				}
				world.gop.fence();
			}
			normalize(world,roots[r]);
		}
	}

	/// compute the overlap between 2 sets of roots
	Tensor<double> overlap(const std::vector<root>& r1,
			const std::vector<root>& r2) const {

		Tensor<double> ovlp(r1.size(),r2.size());
		for (std::size_t i=0; i<r1.size(); ++i) {
			const vecfuncT& x1=r1[i].x;
			for (std::size_t j=0; j<r2.size(); ++j) {
				const vecfuncT& x2=r2[j].x;
				ovlp(i,j)=inner(world,x1,x2).sum();
			}
		}
		return ovlp;
	}


	/// compute the oscillator strength in the length representation

	/// the oscillator strength is given by
	/// \f[
	/// f = 2/3 * \omega |<x | \vec \mu | i >| ^2 * 2
	/// \f]
	/// where \f$ x \f$ is the excited state, and \f$ i \f$ is the ground state
	/// @param[in]	root	a converged root
	double oscillator_strength_length(const root& root) const {
		Tensor<double> mu_if(3);
		for (int idim=0; idim<3; idim++) {
		    real_function_3d ri = real_factory_3d(world).functor2(xyz(idim));
		    vecfuncT amo_times_x=mul(world,ri,active_mo());
			Tensor<double> a=inner(world,amo_times_x,root.x);
			mu_if(idim)=a.sum();
		}
		const double f= 2.0/3.0 * root.omega * mu_if.sumsq() * 2.0;
		return f;
	}

	/// compute the oscillator strength in the velocity representation

	/// the oscillator strength is given by
	/// \f[
	/// f = 2/(3 * \omega) |<x | \vec p | i >| ^2 * 2
	/// \f]
	/// where \f$ x \f$ is the excited state, and \f$ i \f$ is the ground state
	/// @param[in]	root	a converged root
	double oscillator_strength_velocity(const root& root) const {
		Tensor<double> p_if(3);
		// compute the derivatives of the MOs in all 3 directions
		for (int idim=0; idim<3; idim++) {
            real_derivative_3d D = free_space_derivative<double,3>(world, idim);
		    vecfuncT Damo=apply(world,D,active_mo());
			Tensor<double> a=inner(world,Damo,root.x);
			p_if(idim)=a.sum();
		}
		const double f= 2.0/(3.0 * root.omega) * p_if.sumsq() * 2.0;
		return f;
	}

	/// analyze the root: oscillator strength and contributions from occ
	void analyze(const std::vector<root>& roots) const {

		const int noct=active_mo().size();

		std::vector<root>::const_iterator it;
		int iroot=0;
		for (it=roots.begin(); it!=roots.end(); ++it, ++iroot) {
			std::vector<double> norms=norm2s(world,it->x);

			// compute the oscillator strengths
			double osl=this->oscillator_strength_length(*it);
			double osv=this->oscillator_strength_velocity(*it);

			std::cout << std::scientific << std::setprecision(10) << std::setw(20);
			if (world.rank()==0) {
				std::cout.width(10); std::cout.precision(8);
	    		print("excitation energy for root ",iroot,": ",it->omega);
				print("oscillator strength (length)    ", osl);
				print("oscillator strength (velocity)  ", osv);

				// print out the most important amplitudes
				print("\n  dominant contributions ");

				for (std::size_t p=0; p<noct; ++p) {
					const double amplitude=norms[p]*norms[p];
					if (amplitude > 0.1) {
						std::cout << "  norm(x_"<<p<<") **2  ";
						std::cout.width(10); std::cout.precision(6);
						std::cout << amplitude << std::endl;
					}
				}
			}
		}
	}

	/// return the active MOs only, note the shallow copy
	const vecfuncT active_mo() const {
		const vecfuncT& amo=get_calc().amo;
	    vecfuncT actmo;
	    for (std::size_t i=nfreeze_; i<amo.size(); ++i) actmo.push_back(amo[i]);
	    return actmo;
	}

};

}




int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    if (world.rank() == 0) {
    	print("\n  TDHF -- time-dependent Hartree-Fock in the CIS approximation  \n");
    	printf("starting at time %.1f\n", wall_time());
       	print("\nmain() compiled at ",__TIME__," on ",__DATE__);

    }
    startup(world,argc,argv);
    std::cout.precision(6);
    typedef std::vector<functionT> vecfuncT;

    // take the HF orbitals to start
    const std::string input="input";
	Calculation calc(world,input.c_str());
    calc.molecule.print();
    print("\n");
    calc.param.print(world);

    // solve the ground state energy; calc is a reference
    MolecularEnergy me(world,calc);
    double hf_energy=me.value(calc.molecule.get_all_coords());

    if (world.rank()==0) print("MRA hf energy", hf_energy);
    if (world.rank()==0) {
    	printf("\n\n starting TDHF section at time %.1f\n",wall_time());
    	print("nuclear repulsion: ", calc.molecule.nuclear_repulsion_energy());
    	print("hf energy:         ", hf_energy);
    	print("orbital energies:  ");
    	for (std::size_t i=0; i<calc.amo.size(); ++i)
    		print("     ",calc.aeps[i]);
    }



    // construct the CIS solver, it requires a converged HF reference
    CIS cis(world,calc,input);

    // print grid information to file to get a better guess from external
    if (cis.print_grid()) {
    	if (world.rank()==0) print("printing grid for koala\n");
    	real_function_3d density=cis.get_calc().make_density(world,calc.aocc,
    		calc.amo);
    	density.get_impl()->print_grid("grid");
    } else {

    	// solve the response equation
    	cis.solve();
    }

    if (world.rank() == 0) printf("finished at time %.1f\n", wall_time());
    finalize();
    return 0;
}
