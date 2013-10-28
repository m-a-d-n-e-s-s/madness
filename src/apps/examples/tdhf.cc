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


/// load a converged root from disk

/// @param[in]	world 	the world
/// @param[in]	iroot	the i-th root
/// @param[inout]	x	the x-vector for the i-th root
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

/// load a converged root from disk

/// @param[in]	world 	the world
/// @param[in]	iroot	the i-th root
/// @param[inout]	x	the x-vector for the i-th root
/// @return	successfully loaded a root or not
bool save_root(World& world, const int i, const vecfuncT& x, const double omega) {
	std::string filename="root_"+stringify(i);
	archive::ParallelOutputArchive ar(world, filename.c_str(), 1);
	ar & omega;
	for (std::size_t i=0; i<x.size(); ++i) ar & x[i];
	return true;
}

/// normalize the excitation amplitudes

/// @param[in]		world the world
/// @param[inout]	x	the excitation vector
/// @return			the squared amplitudes for each x
std::vector<double> normalize(World& world, vecfuncT& x) {
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
bool iterate_CCS(World& world, const vecfuncT& amo, double& omega,
		vecfuncT& x, const HartreeFock& hf, const std::vector<vecfuncT >& excited_states,
		const int iteration, const double& thresh, std::vector<double>& amplitudes) {

	// for convenience
	const int nmo=amo.size();
	// the projector on the unperturbed density
	Projector<double,3> rho0(amo);

	// a poisson solver
    std::shared_ptr<real_convolution_3d> poisson
    	=std::shared_ptr<real_convolution_3d>
    	(CoulombOperatorPtr(world,lo,bsh_eps));
    // the local potential V^0 of Eq. (4)
    real_function_3d vlocal = hf.get_nuclear_potential() + hf.get_coulomb_potential();

    // make the potential for V0*xp
	vecfuncT Vx=mul(world,vlocal,x);

	// and the exchange potential is K xp
	vecfuncT Kx=hf.get_calc().apply_hf_exchange(world,hf.get_calc().aocc,amo,x);

	// sum up: V0 xp = V_loc xp - K xp
	vecfuncT V0=sub(world,Vx,Kx);

	// now construct the two-electron contribution Gamma, cf Eqs. (7,8)
	vecfuncT Gamma(nmo);
	real_function_3d rhoprime=real_factory_3d(world);

	// the Coulomb part Eq. (7) and Eq. (3)
	for (int i=0; i<nmo; ++i) rhoprime+=x[i]*amo[i];
	real_function_3d pp=2.0*(*poisson)(rhoprime);

	for (int p=0; p<nmo; ++p) {

		// the Coulomb part Eq. (7) and Eq. (3)
		Gamma[p]=pp*amo[p];

		// the exchange part Eq. (8)
		for (int i=0; i<nmo; ++i) Gamma[p]-=x[i]*(*poisson)(amo[i]*amo[p]);

		// project out the zeroth-order density Eq. (4)
		Gamma[p]-=rho0(Gamma[p]);
	}

	// add the local potential and the electron interaction term
	vecfuncT Vphi=add(world,V0,Gamma);
	scale(world,Vphi,-2.0);

	// the bound-state helmholtz function for omega < orbital energy
    std::vector<poperatorT> bsh(nmo);
    for(int p = 0; p<nmo; ++p){
        double eps = hf.orbital_energy(p) + omega;
        if(eps > 0){
            if(world.rank() == 0)
            	print("bsh: warning: positive eigenvalue", p, eps);
            eps = -0.03;
        }
        bsh[p] = poperatorT(BSHOperatorPtr3D(world, sqrt(-2.0 * eps), lo, bsh_eps));
    }
    vecfuncT GVphi=apply(world,bsh,Vphi);


	// update the excitation energy omega
	double t1=0.0;
	double t2=0.0;
	for (std::size_t p=0; p<amo.size(); ++p) {
		// remove factor 2 from Vphi coming from the BSH application
		t1+=0.5*(inner(Vphi[p],GVphi[p]) - inner(Vphi[p],x[p]));
		double n=GVphi[p].norm2();
		t2+=n*n;
	}
	double delta=-t1/t2;
	// do some damping for the first few iterations
	if (iteration<5) ;
	else if (iteration<8) omega+=0.1*delta;
	else if (iteration<10) omega+=0.3*delta;
	else if (iteration<15) omega+=0.5*delta;
	else omega+=delta;
	if (world.rank()==0) {
		std::cout << " iteration " << std::setw(3) << iteration << "  ";
		std::cout.width(16);
		std::cout.precision(10);
		std::cout << std::scientific;
		std::cout << omega << " " << delta << std::endl;
	}

	// update the x vector: orthogonalize against the occupied space
	for (int p=0; p<nmo; ++p) {
		GVphi[p] -= rho0(GVphi[p]);
		for (std::size_t r=0; r<excited_states.size(); ++r) {
			const real_function_3d& x_rp=excited_states[r][p];
			GVphi[p] -= x_rp*inner(x_rp,GVphi[p]);
		}
		x[p]=GVphi[p];
	}

	// normalize the transition density for all occupied orbitals
	const std::vector<double> norms=normalize(world,x);

	// convergence flag
	bool converged=true;

	// this is not converted if the energy changes a lot (> thresh)
	if (std::fabs(delta)>thresh) converged=false;

	// this is not converged if the amplitudes change by a lot (> 1 percent)
	for (std::size_t i=0; i<norms.size(); ++i) {
		if (std::fabs(amplitudes[i]/norms[i] - 1.0)>0.01) {
			print("not converged: ",i,amplitudes[i],norms[i]);
			converged=false;
		}
	}
	amplitudes=norms;

	// print out the most important amplitudes
	for (int p=0; p<nmo; ++p) {
		if (norms[p] > 0.3) {
			if (world.rank()==0) {
				std::cout << "norm(x_"<<p<<") **2  " << std::fixed;
				std::cout.width(10); std::cout.precision(6);
				std::cout << norms[p] << std::endl;
			}
		}
	}

	// if convergence is reached print out weights and leave
	if (converged) {
		if (world.rank()==0) {
			print("\nconverged:  ",omega);
			for (int p=0; p<nmo; ++p) {
				std::cout << "norm(x_"<<p<<") **2  " << amplitudes[p] << std::endl;
			}
		}
		return true;
	}
	return false;
}



int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    if (world.rank() == 0) {
    	print("\n  TDHF -- time-dependent Hartree-Fock in the CCS approximation  \n");
    	printf("starting at time %.1f\n", wall_time());
    }
    startup(world,argc,argv);
    std::cout.precision(6);
    typedef std::vector<functionT> vecfuncT;

    // take the HF orbitals to start
    const std::string input="input";
	Calculation calc(world,input.c_str());
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

    // for convenience
    const std::size_t nmo=hf.get_calc().amo.size();
    const std::size_t nroot=3;
    std::vector<vecfuncT> excited_states;

    // loop over all roots
    for (std::size_t iroot=0; iroot<nroot; ++iroot) {

    	if (world.rank()==0) print("\nworking on root ",iroot,"\n");
        // take virtual orbitals as a start guess for the transition density
        vecfuncT x(nmo);
        for (std::size_t i=0; i<nmo; ++i) x[i]=real_factory_3d(world);

        // guess an excitation energy
        double omega=0.4;

        // guess a HOMO-LUMO excitation, unless there's a root on disk
    	if (not load_root(world,iroot,x,omega)) {
    		real_function_3d all_virtuals=real_factory_3d(world);
    		for (std::size_t ivir=nmo; ivir<hf.get_calc().ao.size(); ++ivir) {
        		all_virtuals+=hf.get_calc().ao[ivir];
    		}
    		const double norm=all_virtuals.norm2();
    		all_virtuals.scale(1/(norm*norm));
    		for (std::size_t iocc=0; iocc<nmo; ++iocc) {
        		x[iocc]=all_virtuals;
    		}
    	}

    	// the excitation amplitudes for each occupied orbital
    	std::vector<double> amplitudes(nmo,1.0);
        for (int iter=0; iter<30; ++iter) {
        	bool converged=iterate_CCS(world,hf.get_calc().amo,omega,x,hf,
        			excited_states,iter,hf.get_calc().param.econv,amplitudes);
        	if (converged) break;
        }
        save_root(world,iroot,x,omega);

        excited_states.push_back(x);

        if (world.rank()==0) print("excitation energy ", omega,"\n");

    }



    if (world.rank() == 0) printf("finished at time %.1f\n", wall_time());
    finalize();
    return 0;
}
