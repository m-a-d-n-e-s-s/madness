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
  \brief compute the time-dependent HF equations (currently in the CCS approximation)

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local
  /trunk/src/apps/examples/tdhf.cc>here</a>.

*/

#include <examples/mp2.h>


#include <mra/operator.h>
#include <mra/mra.h>
#include <mra/operator.h>
#include <mra/lbdeux.h>

using namespace madness;

namespace madness {



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


/// POD holding the excitation energy and the response vector for a single excitation
struct root {
	root(vecfuncT& x, double omega) : x(x), omega(omega) {}
	vecfuncT x;
	double omega;
};


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

public:

	/// ctor

	/// @param[in]	hf	the HartreeFock reference state
	CIS(World& world, const HartreeFock& hf) : world(world), hf_(hf) {
	}

	/// return the HF reference
	const HartreeFock& hf() const {return hf_;}

	/// solve the CIS equations for n roots

	/// @param[in] nroot	the number of roots to solve
	void solve(const std::size_t& nroot) {

		// for convenience
	    const std::size_t nmo=hf().get_calc().amo.size();

	    // this is common for all roots in all iterations
	    exchange_intermediate_=make_exchange_intermediate(hf().get_calc().amo);

	    // loop over all roots
	    for (std::size_t iroot=0; iroot<nroot; ++iroot) {

	    	if (world.rank()==0) print("\nworking on root ",iroot,"\n");
	        // take virtual orbitals as a start guess for the transition density
	        vecfuncT x(nmo);
	        for (std::size_t i=0; i<nmo; ++i) x[i]=real_factory_3d(world);

	        // guess an excitation energy: 0.9* HOMO
	        double omega=-0.9*hf().get_calc().aeps(nmo-1);

	        // check if there's a root on disk. Otherwise take a linear
	        // combination of all orbitals as guess, because there are
	        // not enough virtuals in the minimal basis set for all
	        // possible symmetries
	    	if (not load_root(world,iroot,x,omega)) {
	    		real_function_3d all_virtuals=real_factory_3d(world);
	    		for (std::size_t ivir=0; ivir<hf().get_calc().ao.size(); ++ivir) {
	        		all_virtuals+=hf().get_calc().ao[ivir];
	    		}
	    		const double norm=all_virtuals.norm2();
	    		all_virtuals.scale(1.0/(norm*norm));
	    		for (std::size_t iocc=0; iocc<nmo; ++iocc) {
	        		x[iocc]=all_virtuals;
	    		}
	    	}

	    	// KAIN solver
	        XNonlinearSolver<F,double,allocator> solver =
	        		XNonlinearSolver<F,double,allocator>(allocator(world,nmo));
	        solver.set_maxsub(3);

	        // print progress on the computation
	        if (world.rank()==0) {
	        	print(" iteration    excitation energy   energy correction         error");
	        }
	    	// the excitation amplitudes for each occupied orbital
	    	std::vector<double> amplitudes(nmo,1.0);
	        for (int iter=0; iter<50; ++iter) {
	        	bool converged=iterate_CCS(world,solver,omega,x,
	        			roots(),iter,hf().get_calc().param.econv,amplitudes);
	        	// save for restart
	        	save_root(world,iroot,x,omega);
	        	if (converged) break;
	        }

	        roots().push_back(root(x,omega));
	        double osl=this->oscillator_strength_length(roots().back());
	        double osv=this->oscillator_strength_velocity(roots().back());

			std::cout << std::scientific << std::setprecision(10) << std::setw(20);
	        if (world.rank()==0) {
	        	print("excitation energy              ", omega);
	        	print("oscillator strength (length)   ", osl);
	        	print("oscillator strength (velocity) ", osv);
	        }

	    }
	}

	/// return the roots of the response equation
	std::vector<root>& roots() {return roots_;}

private:

	/// the world
	World& world;

	/// the HartreeFock reference state
	const HartreeFock& hf_;

	/// the excited states aka the roots of the response equation
    std::vector<root> roots_;

    /// intermediate for the two-electron interaction term, Eq. (8)
    std::vector<vecfuncT> exchange_intermediate_;

	/// iterate the TDHF or CCS equations

	/// follow Eq (4) of
	/// T. Yanai, R. J. Harrison, and N. Handy,
	/// ÒMultiresolution quantum chemistry in multiwavelet bases: time-dependent density
	/// functional theory with asymptotically corrected potentials in local density and
	/// generalized gradient approximations,Ó Mol. Phys., vol. 103, no. 2, pp. 413Ð424, 2005.
	///
	/// The convergence criterion is that the excitation amplitudes don't change
	/// @param[in]		amo	the unperturbed orbitals
	/// @param[in]		omega	the lowest excitation energy
	/// @param[inout]	x	the response function
	/// @param[in]		hf 	the HF reference object
	/// @param[in]		excited_states all lower-lying excited states for orthogonalization
	/// @param[in]	iteration	the current iteration
	template<typename solverT>
	bool iterate_CCS(World& world, solverT& solver, double& omega,
			vecfuncT& x, const std::vector<root >& excited_states,
			const int iteration, const double& thresh, std::vector<double>& amplitudes) const {

		// for convenience
		const vecfuncT& amo=hf().get_calc().amo;
		const int nmo=amo.size();
		// the projector on the unperturbed density
		Projector<double,3> rho0(amo);

		// a poisson solver
	    std::shared_ptr<real_convolution_3d> poisson
	    	=std::shared_ptr<real_convolution_3d>
	    	(CoulombOperatorPtr(world,lo,bsh_eps));
	    // the local potential V^0 of Eq. (4)
	    real_function_3d vlocal = hf().get_nuclear_potential() + hf().get_coulomb_potential();

	    // make the potential for V0*xp
		vecfuncT Vx=mul(world,vlocal,x);

		// and the exchange potential is K xp
		vecfuncT Kx=hf().get_calc().apply_hf_exchange(world,hf().get_calc().aocc,amo,x);

		// sum up: V0 xp = V_loc xp - K xp
		vecfuncT V0=sub(world,Vx,Kx);

		// now construct the two-electron contribution Gamma, cf Eqs. (7,8)
		real_function_3d rhoprime=real_factory_3d(world);

		// the Coulomb part Eq. (7) and Eq. (3)
		for (int i=0; i<nmo; ++i) rhoprime+=x[i]*amo[i];
		real_function_3d pp=2.0*(*poisson)(rhoprime);

		// the Coulomb part Eq. (7) and Eq. (3)
		vecfuncT Gamma=mul(world,pp,amo);

		// the exchange part Eq. (8)
		for (int p=0; p<nmo; ++p) {

			// this is x_i * \int 1/r12 \phi_i \phi_p
			vecfuncT x_Ppi=mul(world,x,exchange_intermediate_[p]);
			for (int i=0; i<nmo; ++i) Gamma[p]-=x_Ppi[i];

			// project out the zeroth-order density Eq. (4)
			Gamma[p]-=rho0(Gamma[p]);
		}

		// add the local potential and the electron interaction term
		vecfuncT Vphi=add(world,V0,Gamma);
		scale(world,Vphi,-2.0);
		truncate(world,Vphi);

		// the bound-state helmholtz function for omega < orbital energy
	    std::vector<poperatorT> bsh(nmo);
	    for(int p = 0; p<nmo; ++p){
	        double eps = hf().orbital_energy(p) + omega;
	        if(eps > 0){
	            if(world.rank() == 0)
	            	print("bsh: warning: positive eigenvalue", p, eps);
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
		for (std::size_t p=0; p<amo.size(); ++p) {
			double n=GVphi[p].norm2();
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
		double error=sqrt(inner(F(world,residual),F(world,residual)));
        if (world.rank()==0) {
			std::cout << std::setw(6) << iteration << "     ";
			std::cout << std::scientific << std::setprecision(10);
			std::cout << std::setw(20) << omega	<< std::setw(20)<< delta
					<< std::setw(20) << error << " " << wall_time()<< std::endl;
        }

		// generate a new trial solution vector
	    if (iteration<17) {
	    	x=GVphi;
	    } else {
	    	F ff=solver.update(F(world,x),F(world,residual));
	    	x=ff.x;
	    }

	    // update the x vector: orthogonalize against the occupied space
		for (int p=0; p<nmo; ++p) {
			x[p] -= rho0(x[p]);
			for (std::size_t r=0; r<excited_states.size(); ++r) {
				const real_function_3d& x_rp=excited_states[r].x[p];
				x[p] -= x_rp*inner(x_rp,x[p]);
			}
		}

		// normalize the transition density for all occupied orbitals
		const std::vector<double> norms=normalize(world,x);
		truncate(world,x);

		// convergence flag
		bool converged=true;

		// this is not converged if the amplitudes change by a lot (> 1 percent)
		if (error>sqrt(thresh)) converged=false;

		// this is not converted if the energy changes a lot (> thresh)
		if (std::fabs(delta)>thresh) converged=false;

		// this is not converged if the amplitudes change by a lot (> 1 percent)
		for (std::size_t i=0; i<norms.size(); ++i) {
			if (std::fabs(amplitudes[i]/norms[i] - 1.0)>0.01) {
				converged=false;
			}
		}
		amplitudes=norms;

		// if convergence is reached print out weights and leave
		if (converged) {
			if (world.rank()==0) {
				std::cout.width(10); std::cout.precision(6);
				print("\nconverged:  ",omega);
		    	std::cout << std::fixed;
				for (int p=0; p<nmo; ++p) {
					std::cout << "norm(x_"<<p<<") **2  " << amplitudes[p] << std::endl;
				}
			}
			return true;
		}
		return false;
	}

	/// make the 2-electron interaction intermediate

	/// the intermediate is teh same for all roots:
	/// \[
	///   int[i,p] = \int 1/r12 \phi_i(1) * \phi_p(1)
	/// \]
	/// @param[in]	amo	active orbitals in the CIS computation
	/// @return		a vector of vectors of functions
	std::vector<vecfuncT> make_exchange_intermediate(const vecfuncT& amo) const {
		// a poisson solver
	    std::shared_ptr<real_convolution_3d> poisson
	    	=std::shared_ptr<real_convolution_3d>
	    	(CoulombOperatorPtr(world,lo,bsh_eps));

	    std::vector<vecfuncT> intermediate(amo.size());
	    for (std::size_t i=0; i<amo.size(); ++i) {
			intermediate[i]=apply(world,(*poisson),mul(world,amo[i],amo));
	    }
	    return intermediate;
	}

	/// load a converged root from disk

	/// @param[in]	world 	the world
	/// @param[in]	iroot	the i-th root
	/// @param[inout]	x	the x-vector for the i-th root
	/// @param[out]	omega	the excitation energy
	/// @return	successfully loaded a root or not
	bool load_root(World& world, const int i, vecfuncT& x, double& omega) {
		std::string filename="root_"+stringify(i);
		bool exists=archive::ParallelInputArchive::exists(world,filename.c_str());
		if (not exists) return false;

		archive::ParallelInputArchive ar(world, filename.c_str(), 1);
		ar & omega;
		for (std::size_t i=0; i<x.size(); ++i) ar & x[i];
		return true;
	}

	/// save a converged root to disk

	/// @param[in]	world 	the world
	/// @param[in]	iroot	the i-th root
	/// @param[inout]	x	the x-vector for the i-th root
	/// @param[in]	omega	the excitation energy
	void save_root(World& world, const int i, const vecfuncT& x, const double omega) {
		std::string filename="root_"+stringify(i);
		archive::ParallelOutputArchive ar(world, filename.c_str(), 1);
		ar & omega;
		for (std::size_t i=0; i<x.size(); ++i) ar & x[i];
	}

	/// normalize the excitation amplitudes

	/// normalize the set of excitation amplitudes, such that the sum of square
	/// of all amplitudes equals 1.
	/// @param[in]		world the world
	/// @param[inout]	x	the excitation vector
	/// @return			the squared amplitudes for each x
	std::vector<double> normalize(World& world, vecfuncT& x) const {
		std::vector<double> n=norm2s(world,x);
		double totalnorm2=0.0;
		for (std::size_t i=0; i<x.size(); ++i) totalnorm2+=n[i]*n[i];
		totalnorm2=sqrt(totalnorm2);
		for (std::size_t i=0; i<x.size(); ++i) {
			x[i].scale(1.0/totalnorm2);
			n[i]=n[i]*n[i]/totalnorm2;
		}
		return n;
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
		    real_function_3d ri  = real_factory_3d(world).functor2(xyz(idim));
		    vecfuncT amo_times_x=mul(world,ri,hf().get_calc().amo);
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
		    vecfuncT Damo=apply(world,D,hf().get_calc().amo);
			Tensor<double> a=inner(world,Damo,root.x);
			p_if(idim)=a.sum();
		}
		const double f= 2.0/(3.0 * root.omega) * p_if.sumsq() * 2.0;
		return f;
	}

};

}




int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    if (world.rank() == 0) {
    	print("\n  TDHF -- time-dependent Hartree-Fock in the CCS approximation  \n");
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

    HartreeFock hf(world,calc);
    const double hf_energy=hf.value();
    if (world.rank()==0) {
    	printf("\n\n starting TDHF section at time %.1f\n",wall_time());
    	print("nuclear repulsion: ", hf.get_calc().molecule.nuclear_repulsion_energy());
    	print("hf energy:         ", hf_energy);
    	print("orbital energies:  ");
    	for (std::size_t i=0; i<hf.get_calc().amo.size(); ++i)
    		print("     ",hf.get_calc().aeps[i]);
    }


    // construct the CIS solver, it requires a converged HF reference
    CIS cis(world,hf);

    // solve the response equation for n roots;
    cis.solve(3);

    if (world.rank()==0) {
    	print("CIS solver ended");
    	std::cout << std::fixed;
    	for (std::size_t i=0; i<cis.roots().size(); ++i) {

			std::cout.width(10); std::cout.precision(6);
    		print("excitation energy for root ",i,": ",cis.roots()[i].omega);

    		// print out the most important amplitudes
    		print("  dominant contributions ");
    		for (std::size_t p=0; p<hf.nocc(); ++p) {
    		functionT xp=cis.roots()[i].x[p];
    			double amplitude=xp.norm2();
    			amplitude*=amplitude;
    			if (amplitude > 0.1) {
    				if (world.rank()==0) {
    					std::cout << "  norm(x_"<<p<<") **2  ";
    					std::cout.width(10); std::cout.precision(6);
    					std::cout << amplitude << std::endl;
    				}
    			}
    		}
    	}
    }



    if (world.rank() == 0) printf("finished at time %.1f\n", wall_time());
    finalize();
    return 0;
}
