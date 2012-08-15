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
  \file examples/oep.cc
  \brief optimized effective potentials for DFT

  The source is <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local/trunk/src/apps/examples/oep.cc>here</a>.

*/

#include <examples/mp2.h>

#include <mra/operator.h>
#include <mra/mra.h>
#include <mra/mraimpl.h>
#include <mra/operator.h>
#include <mra/lbdeux.h>

using namespace madness;

static const double   rcut = 0.01; // Smoothing distance in 1e potential
static const double d12cut = 0.01; // Smoothing distance in wave function

typedef Tensor<double> tensorT;

static double rr(const coord_3d& r) {
	return sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
}

// Smoothed 1/r potential (c is the smoothing distance)
static double u(double r, double c) {
    r = r/c;
    double r2 = r*r, pot;
    if (r > 6.5){
        pot = 1.0/r;
    } else if (r > 1e-2) {
        pot = erf(r)/r + exp(-r2)*0.56418958354775630;
    } else{
        pot = 1.6925687506432689-r2*(0.94031597257959381-r2*(0.39493270848342941-0.12089776790309064*r2));
    }
    
    return pot/c;
}


struct LBCost {
    double leaf_value;
    double parent_value;
    LBCost(double leaf_value=1.0, double parent_value=1.0) 
        : leaf_value(leaf_value)
        , parent_value(parent_value) 
    {}

    double operator()(const Key<6>& key, const FunctionNode<double,6>& node) const {
        if (key.level() <= 1) {
            return 100.0*(leaf_value+parent_value);
        }
        else if (node.is_leaf()) {
            return leaf_value;
        }
        else {
            return parent_value;
        }
    }
};


static Tensor<double> read_grid(int k, std::string filename) {

    FILE* file = fopen(filename.c_str(),"r");

    Tensor<double> a(k);
    for (int i=0; i<k; ++i) {
    	double c;
    	if (fscanf(file,"%lf",&c) != 1) {
    		print("read_grid: failed reading data from file", filename);
    		MADNESS_EXCEPTION("",0);
    	}
    	a(i) = c;
    }
    return a;
}

struct coulomb {
	coulomb() {}
	double operator()(const coord_3d& xyz) const {
		return -4.0*u(rr(xyz),0.000001);
	}
};


/// wrapper class for the optimized effective potential

/// The potential is given on a radial grid. Interpolate if necessary
struct recpot {

	/// the number of grid points
	long k;

	/// the effective potential tabulated
	tensorT potential;

	/// the nuclear charge
	double Z;

	/// ctor with the nuclear charge to subtract for better interpolation

	/// @param[in]	file_grid 	name of the file with grid points
	/// @param[in]	file_pot1 	name of the file with the potential on the grid points
	/// @param[in]	file_pot2 	name of the other file with the potential on the grid points
	/// @param[in]	Z			nuclear charge (unused at the moment)
	recpot(const std::string file_grid, const std::string file_pot1, const std::string file_pot2, const long Z) : Z(double(Z)) {

		// number of grid points
		k=400;

		// read the grid points and the potential values
		potential=tensorT(2,k);
		// this is an assignment of a slice
		potential(0,_)=read_grid(k,file_grid);
		potential(1,_)=read_grid(k,file_pot1);
		potential(1,_)+=read_grid(k,file_pot2);

		// subtract the nuclear potential; must be consistent with operator()
		for (int i=0; i<k; ++i) {
			double r=potential(0,i);
			potential(1,i)+=Z*u(r,1.e-6);
		}
	}

	/// return the value of the potential; interpolate if necessary

	/// @param[in]	xyz		cartesian coordinates of a point
	/// @return		the value of the potential at point xyz
	double operator()(const coord_3d& xyz) const {
		double r=rr(xyz);
		return (interpolate(r)- Z* u(r,1.e-6));
	}

	/// interpolate the radial potential from the grid
	double interpolate(const double& r) const {

		int i=0;
		// upon loop exit i will be the index with the grid point right behind r
		if (r<1.0) i=195;
		if (r<0.3) i=145;
		if (r<0.1) i=100;
		if (r<0.02) i=45;
		for (i=0; i<k; ++i) {
			if (potential(0,i)<r) continue;
			break;
		}

		// fit a linear curve: f(x)=a*x + b
		double delta_y=potential(1,i) - potential(1,i-1);
		double delta_x=potential(0,i) - potential(0,i-1);
		double a=delta_y/delta_x;

		double dx=r-potential(0,i-1);
		double val=potential(1,i-1) + dx*a;
		return val;
	}

};

void plot(const real_function_3d& f, const std::string filename, const long k) {
    FILE* file = fopen(filename.c_str(),"w");

    for (int i=0; i<k; ++i) {
    	double z=0.001+double(i)*0.01;
    	coord_3d r=vec(0.0,0.0,z);
    	double c=f(r);
    	fprintf(file,"%lf %lf\n",z,c);
    }
    fclose(file);
}

// print the radial density: r^2 rho
void plot_radial_density(const real_function_3d& rho, const std::string filename, const tensorT& grid) {
	FILE* file = fopen(filename.c_str(),"w");
	for (int i=0; i<grid.dim(0); ++i) {
		double r=grid(i);
		coord_3d xyz=vec(0.0,0.0,r);
		double c=r*r*rho(xyz);
		fprintf(file,"%lf %lf\n",r,c);
	}
	fclose(file);
}

void compute_energy(World& world, const real_function_3d& psi, const real_function_3d& pot, double& ke, double& pe) {

	double kinetic_energy = 0.0;
    for (int axis=0; axis<3; axis++) {
    	real_derivative_3d D = free_space_derivative<double,3>(world, axis);
    	real_function_3d dpsi = D(psi);
    	kinetic_energy += 0.5*inner(dpsi,dpsi);
    }
    ke=kinetic_energy;

    pe=inner(psi,pot*psi);
    if(world.rank() == 0) {
        printf("compute the energy at time   %.1fs\n", wall_time());
        printf("kinetic energy      %12.8f\n", ke);
        printf("potential energy    %12.8f\n", pe);
        printf("total energy        %12.8f\n", pe+ke);
    }
}

/// apply the Green's function on V*psi, update psi and the energy
void iterate(World& world, const real_function_3d& V, real_function_3d& psi, double& eps) {

	real_function_3d Vpsi = (V*psi);
    Vpsi.scale(-2.0).truncate();

    real_convolution_3d op = BSHOperator3D(world, sqrt(-2*eps), 1.e-6, 1e-6);
    real_function_3d tmp=op(Vpsi).truncate();

    double norm = tmp.norm2();
    real_function_3d r = tmp-psi;
    double rnorm = r.norm2();
    double eps_new = eps - 0.5*inner(Vpsi,r)/(norm*norm);
    if (world.rank() == 0) {
        print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
        print("eps_new",eps_new);
    }
    psi = tmp.scale(1.0/norm);
    if (eps_new<0.0) eps = eps_new;
}

/// orthogonalize orbital i against all other orbitals
void orthogonalize(std::vector<real_function_3d>& orbitals, const int ii) {
	MADNESS_ASSERT(ii<orbitals.size());

	real_function_3d& phi=orbitals[ii];

	// loop over all other orbitals
	for (int i=0; i<ii; ++i) {
		const real_function_3d orbital=orbitals[i];
		double ovlp=inner(orbital,phi);
		double norm=orbital.norm2();
		phi-=(ovlp/norm/norm)*orbital;

	}

	double n=phi.norm2();
	phi.scale(1.0/n);
}

/// solve the residual equations

/// @param[in]		potential	the effective potential
/// @param[inout]	eps			guesses for the orbital energies
/// @param[inout]	orbitals	the first n roots of the equation
void solve(World& world, const real_function_3d& potential, tensorT& eps,
		std::vector<real_function_3d>& orbitals) {

	const long nroots=eps.size();
	print("solving for",nroots,"roots of the effective potential");

	// loop over all roots
	for (int i=0; i<nroots; ++i) {

		real_function_3d& phi=orbitals[i];
		double eiger=eps(i);

		double ke=0.0,pe=0.0;
		for (int j=0; j<10; ++j) {
			compute_energy(world,phi,potential,ke,pe);
			iterate(world,potential,phi,eiger);
			orthogonalize(orbitals,i);
		}
	}
}


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    if (world.rank() == 0) {
    	print("\n  OEP -- optimized effective potentials for DFT  \n");
    	printf("starting at time %.1f\n", wall_time());
    }
    startup(world,argc,argv);
    std::cout.precision(6);

    // take as a guess the HF orbitals
    const std::string input="input";
	Calculation calc(world,input.c_str());
	HartreeFock hf(world,calc);
    hf.value();


    recpot pot_functor("grid.txt","Be_recpot_num_dzp.txt","Be_startpot.txt",4);
    real_function_3d potential=real_factory_3d(world).functor2(pot_functor).truncate_on_project();;

    tensorT grid(pot_functor.k);
    grid=pot_functor.potential(0,_);

    {
		real_function_3d rho=real_factory_3d(world);
		for (const real_function_3d& orbital : hf.get_calc().amo) {
			rho+=orbital*orbital;
		}
		rho.scale(2.0);
		plot_radial_density(rho,"r2rho_hf",grid);
    }

    potential.print_size("potential");

    // solve the residual equations
    std::vector<real_function_3d> orbitals=hf.get_calc().amo;
    tensorT eps=hf.get_calc().aeps;

    solve(world,potential,eps,orbitals);

    {
		real_function_3d rho=real_factory_3d(world);
		for (const real_function_3d& orbital : orbitals) {
			rho+=orbital*orbital;
		}
		rho.scale(2.0);
		plot_radial_density(rho,"r2rho",grid);
    }

    if (world.rank() == 0) printf("finished at time %.1f\n", wall_time());
    finalize();
    return 0;
}
