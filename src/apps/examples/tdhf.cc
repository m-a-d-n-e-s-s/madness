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
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local/trunk/src/apps/examples/tdhf.cc>here</a>.

*/

#include <examples/mp2.h>


#include <mra/operator.h>
#include <mra/mra.h>
#include <mra/operator.h>
#include <mra/lbdeux.h>

using namespace madness;

/// iterate the TDHF or CCS equations

/// follow Eq (4) of
/// T. Yanai, R. J. Harrison, and N. Handy,
/// ÒMultiresolution quantum chemistry in multiwavelet bases: time-dependent density
/// functional theory with asymptotically corrected potentials in local density and
/// generalized gradient approximations,Ó Mol. Phys., vol. 103, no. 2, pp. 413Ð424, 2005.
/// @param[in]		amo	the unperturbed orbitals
/// @param[in]		omega	the lowest excitation energy
/// @param[inout]	x	the response function
/// @param[in]		hf 	the HF reference object
/// @param[in]	iteration	the current iteration
double iterate_CCS(World& world, const vecfuncT& amo, double& omega, vecfuncT& x,
		HartreeFock& hf, const int iteration) {

	// for convenience
	const int nmo=amo.size();
	// the projector on the unperturbed density
	Projector<double,3> rho0(amo);
	// a poisson solver
    std::shared_ptr<real_convolution_3d> poisson=std::shared_ptr<real_convolution_3d>
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
            if(world.rank() == 0) print("bsh: warning: positive eigenvalue", p, eps);
            eps = -0.03;
        }
        bsh[p] = poperatorT(BSHOperatorPtr3D(world, sqrt(-2.0 * eps),  lo, bsh_eps));
    }
    vecfuncT GVphi=apply(world,bsh,Vphi);


	// update the excitation energy omega
	double t1=0.0;
	double t2=0.0;
	for (int p=0; p<amo.size(); ++p) {
		// remove factor 2 from Vphi coming from the BSH application
		t1+=0.5*(inner(Vphi[p],GVphi[p]) - inner(Vphi[p],x[p]));
		double n=GVphi[p].norm2();
		t2+=n*n;
	}
	double delta=-t1/t2;
	// do some damping for the first few iterations
	if (iteration<5) ;
	else if (iteration<8) omega+=0.05*delta;
	else if (iteration<10) omega+=0.1*delta;
	else omega+=0.3*delta;
	if (world.rank()==0) print("t1, t2, delta, omega",t1,t2,delta,omega);

	// update the x vector: orthogonalize against the occupied space
	for (int p=0; p<nmo; ++p) GVphi[p] -= rho0(GVphi[p]);

	// normalize the transition density for all
	double norm2=0.0;
	for (int p=0; p<nmo; ++p) {
		x[p]=GVphi[p];
		double norm=x[p].norm2();
		norm2+=norm*norm;
		std::cout << "norm(x_"<<p<<")  " << norm << std::endl;
	}
	for (int p=0; p<nmo; ++p) x[p].scale(1.0/sqrt(norm2));


	return omega;
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
    	print("\n\n starting TDHF section \n");
    	print("nuclear repulsion: ", hf.get_calc().molecule.nuclear_repulsion_energy());
    	print("hf energy:         ", hf_energy);
    	print("orbital energies:  ");
    	for (int i=0; i<hf.get_calc().amo.size(); ++i) print("     ",hf.get_calc().aeps[i]);
    }

    // get a virtual orbital as a start guess for the transition density
    const int nmo=hf.get_calc().amo.size();
    vecfuncT xp(nmo);
    for (int i=0; i<nmo; ++i) xp[i]=real_factory_3d(world);
//    for (int i=0; i<hf.get_calc().ao.size(); ++i) xp[0]+=hf.get_calc().ao[i];

    // guess a HOMO-LUMO excitation
    xp[nmo-1]=hf.get_calc().ao[nmo];

    // guess an excitation energy for the h2 molecule
    double omega=0.2;
    for (int iter=0; iter<30; ++iter) {
    	iterate_CCS(world,hf.get_calc().amo,omega,xp,hf,iter);
    }

    if (world.rank()==0) print("excitation energy ", omega);

    if (world.rank() == 0) printf("finished at time %.1f\n", wall_time());
    finalize();
    return 0;
}
