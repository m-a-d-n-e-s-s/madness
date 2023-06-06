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

/*!
  \file examples/gygi_slution.cc
  \brief compute the dielectric cavity and the electrostatic potential of hydrogen atom in water
  \defgroup examplegygi compute the dielectric cavity and the electrostatic potential of hydrogen atom in water
  \ingroup examples

  The source is <a href=https://github.com/m-a-d-n-e-s-s/madness/blob/master/src/examples/gygi_slution.cc>here</a>.

  \par Points of interest
  - compute the dielectric functional (of density)
  - compute the electrostatic potential by convolving the free space Green's function with the effective charge and the induced surface charge
  - compute the derivatives of the dielectric functional and the electrosstatic potential with respect to the density
  \par Background
  - This program is an implementation of the solvation model in THE JOURNAL OF CHEMICAL PHYSICS 124, 074103 (2006)
  - The DFT equation is not solved
  - The test system isa hydrogen atom (1s orbital)
*/
 //We will test this for a hydrogen atom
#include <madness/mra/mra.h>
#include <madness/constants.h>
#include <ctime>
#include <madness/tensor/solvers.h>
#include <madness/mra/funcplot.h>
#include <madness/mra/nonlinsol.h>
using namespace madness;
typedef real_function_3d realfunc;
//define parameter
const double epsilon = 78.304;
const double beta = 1.300;
const double rho_0 = 0.0004;
const double Gamma = 0.07197;
const int k= 8;
const double L = 10;
const double thresh = 1e-6;
const int maxiter = 20;
//Macro for timing
double XXstart;
#define TIME(MSG,X) XXstart=wall_time();          \
                    X; \
		    if (world.rank() == 0)print("timer: ", MSG,"used",wall_time() - XXstart) \

//nucleu charge density
double nuc_func(const coord_3d& r) {
  const double expnt = 100.0;
  const double coeff = -1.0*pow(1.0/constants::pi*expnt,0.5*3);
  return coeff*exp(-expnt*(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));
}
//define the electronic density function of a hydrogen atom
double rho_func(const coord_3d& r){
  double R = sqrt(r[0]*r[0] +r[1]*r[1] + r[2]*r[2]);
  double coef = -1.0/constants::pi;
  return coef*exp(-2.0*R);
}
//define a class for the Gygi potential and the surface cavity
class GygiPot {
private:
  const realfunc& rho;    //density for the poisson's equation
  const realfunc& rhot;    //density for the dielectric functional
  const double& rho_0;  // the density at the solut-solvent interface
  const double& epsilon; //assyptotic value of the dielectric functional
  const int& maxiter;  // maximum iterations in the solver
  World& world;
  const double& Gamma;   //the surface tention of the solvent
  static const double cutrho = 1e-8;  //cutoff value of the density
  //utility function

  //unary operator to determine the reciprocal of a madness function
  template<typename T,int DIM>
  struct Reciprocal {
    void operator()(const Key<DIM>& key, Tensor<T>& t) const {
      UNARY_OPTIMIZED_ITERATOR(T, t,*_p0=1.0/(*_p0));
    }
    template <typename Archive>void serialize(Archive& ar) {}
  };
  //unary operator to raise electronic density to the power 2beta
  template<typename T,int DIM>
  struct Pow {
    void operator()(const Key<DIM>& key,Tensor<T>& t) const {
      UNARY_OPTIMIZED_ITERATOR(T,t,
			       if(*_p0 <cutrho)
			       	 *_p0 = 0.0;
			       else
				 *_p0 = pow(*_p0,2*beta);
			       );
    }
    template <typename Archive>void serialize(Archive& ar) {}
  };
  //unary operator to raise electronic density to the power 2beta - 1
  template<typename T,int DIM>
  struct Pow_beta_one {
    void operator()(const Key<DIM>& key,Tensor<T>& t) const {
      UNARY_OPTIMIZED_ITERATOR(T,t,
			       if (*_p0 < cutrho)
			       	 *_p0 = 0.0;
			       else
				 *_p0 = pow(*_p0,2*beta-1);
			       );
    }
    template <typename Archive>void serialize(Archive& ar) {}
  };
  //density ratio raised to the power 2beta
  realfunc ratio_rho()const {
    double rerho_0 =-1.0/rho_0;
    realfunc rr = rerho_0*rhot;
    rr.unaryop(Pow<double,3>());
    return rr;
  }
  //density ratio raised to the power 2*beta -1
  realfunc rho_beta_one() const {
    double rerho_0 =-1.0/rho_0 ;
    realfunc rr = rerho_0*rhot;
    rr.unaryop(Pow_beta_one<double,3>());
    return rr;
  }
  //reciprocal of 1 + ratio_rho()
  realfunc re_one_plus_ratio_rho() const{
    realfunc value = 1.0 + ratio_rho();
    value.unaryop(Reciprocal<double,3>());
    return value;
  }
  //gradient of density
  realfunc grad_rho(const realfunc& dens) const {
    real_derivative_3d Dx = free_space_derivative<double,3>(dens.world(), 0);
    real_derivative_3d Dy = free_space_derivative<double,3>(dens.world(), 1);
    real_derivative_3d Dz = free_space_derivative<double,3>(dens.world(), 2);
    return (Dx(dens) + Dy(dens) + Dz(dens));
  }
  //reciprocal of gradient of density
  realfunc re_grad_rho(const realfunc& dens) const{
    realfunc value = grad_rho(dens);
    value.unaryop(Reciprocal<double,3>());
    return value;
  }
public:
  //define methods of the class
  //define the dielectric functional
  realfunc epsilon_rho()const {
    return 1.0 + ((epsilon - 1.0)/2.0)*(1.0 +((1.0 -ratio_rho())*re_one_plus_ratio_rho()));
    // return ratio_rho();
  }
  //reciprocal of the dielectric functional
  realfunc re_epsilon_rho() const{
    realfunc value = epsilon_rho();
    value.unaryop(Reciprocal<double,3>());
    return value;
  }
  //derivative of the dielectric functional
  realfunc depsilon_drho() const {
    //this is the derivative w.r.t to rho
    double fac = (1.0 - epsilon)/rho_0;
    realfunc nume =2.0*beta*rho_beta_one();
    realfunc denom = re_one_plus_ratio_rho()*re_one_plus_ratio_rho();
    return fac*nume*denom;
    // return rho_beta_one();
  }
  //derivative of the dielectric functional
  double cavitation_energy() const {
    double fac =1.0/rho_0;
    double quantum_surface = 0.5*fac*beta*grad_rho(rho).norm2();
    double convfact = 6.423049507e-4; // 1N/m = 6.423049507e-4a.u
    return convfact*Gamma*quantum_surface;
  }
  //cavitation potential to be included in the KS equation
  realfunc dfree_drho() const {
    double fac = beta/(4.0*rho_0);
    real_derivative_3d Dx = free_space_derivative<double,3>(rhot.world(), 0);
    real_derivative_3d Dy = free_space_derivative<double,3>(rhot.world(), 1);
    real_derivative_3d Dz = free_space_derivative<double,3>(rhot.world(), 2);
    realfunc gradsq = Dx(rho)*Dx(rho) + Dy(rho)*Dy(rho) + Dz(rho)*Dz(rho);
    realfunc nume = fac*grad_rho(gradsq);
    realfunc denom =(1.0/grad_rho(rho).norm2())*re_grad_rho(rho);
    double convfact = 6.423049507e-4; // 1N/m = 6.423049507e-4a.u
    return nume*denom*convfact*Gamma;
  }
  /*defining the quantum surface with a Gaussian
  double expnt_surface() const{
    real_derivative_3d Dx = free_space_derivative<double,3>(rhot.world(), 0);
    real_derivative_3d Dy = free_space_derivative<double,3>(rhot.world(), 1);
    real_derivative_3d Dz = free_space_derivative<double,3>(rhot.world(), 2);
    realfunc grad_rho = Dx(rho) + Dy(rho) + Dz(rho);
    double fac = 3.0/(rho_0*sqrt(constants::pi));
    return fac*grad_rho.norm2();
  }
  */
//Define the guess potential for the Poisson's equation
  realfunc GuessPotential(World& world) const {
    double tol = madness::FunctionDefaults<3>::get_thresh();
    real_convolution_3d op = madness::CoulombOperator(world, tol*10.0, tol*0.1);
    //return op(rho);
    return op(re_epsilon_rho()*rho);   //U_0
}
//Molecular potential i.e potential due to the molecular charge distribution
  realfunc MolecularPotential()const {
    return (re_epsilon_rho()*rho);
  }
  //compute the gradient of epsilon[rho]
  realfunc depsilon_dr() const {
    return depsilon_drho()*grad_rho(rhot);
  }
  //compute the surface charge
  realfunc make_surfcharge(const realfunc& u) const {
    realfunc gradU = grad_rho(u);
    //  realfunc eprho = epsilon_rho();
    const double rfourpi = 1.0/(4.0*constants::pi);
    return (rfourpi*re_epsilon_rho()*depsilon_dr()*gradU).truncate();
  }

  //Define the electrostatic potential
  realfunc ESP(World& world)const {
    const bool USE_SOLVER = true;
    double tol = std::max(1e-3,FunctionDefaults<3>::get_thresh());
    real_convolution_3d op = madness::CoulombOperator(world, tol*10.0, tol*0.1);
    realfunc U = GuessPotential(world);
    double unorm = U.norm2();
    print("U.norm2: ", unorm);
    if (USE_SOLVER) {
      madness::NonlinearSolver solver;
      // This section employs a non-linear equation solver from solvers.h
      //  http://onlinelibrary.wiley.com/doi/10.1002/jcc.10108/abstract
      if (world.rank() == 0){
        print("\n\n");//for formating output
	madness::print("    Computing the Electrostatic Potential   ");
	madness::print("           ______________________           \n ");

	madness::print("iteration          residue norm2                soln(10.0)            used \n");
      }
      realfunc uvec, rvec;
      for (int iter=0; iter<maxiter; iter++) {
	double start = wall_time();
        uvec = U;
        realfunc W = MolecularPotential();
        realfunc Scharge = make_surfcharge(U);
        rvec = U + op(W + Scharge);
        realfunc U_new = solver.update(uvec,rvec);
        double err = rvec.norm2();
        if (world.rank()==0)
	  madness::print("  ", iter,"             " , err,"           ",U(coord_3d(10.0)),"   ", wall_time()-start );
        if (err >0.3*unorm){ U = 0.5*U + 0.5*U_new;
        }
        else
	  U = U_new;
        if(err < 10.0*tol) break;
      }
    }
    return U ;
  }
  //Define the dielectric potential
  realfunc V_epsilon(const realfunc& u) const {
    realfunc grad_ESP = grad_rho(u);
    double coef = -1.0/(8.0*constants::pi);
    return coef*grad_ESP*grad_ESP*depsilon_drho();
  }
  //Defining the gradient of the ESP w.r.t rho
  realfunc dESP_drho(const realfunc u) const {
    return u +V_epsilon(u);
  }
  //Defining the Constructor
  GygiPot(const realfunc& rho,
	  const realfunc& rhot,
	  const double& rho_0,
	  const double& epsilon,
	  const int& maxiter,
	  World& world,
	  const double Gamma)
    :rho(rho),
    rhot(rhot),
    rho_0(rho_0),
    epsilon(epsilon),
    maxiter(maxiter),
    world(world),
    Gamma(Gamma){}
  };
int main(int argc, char **argv){
  initialize(argc,argv);
  World world(SafeMPI::COMM_WORLD);
  startup(world,argc,argv);
  //Function default
  FunctionDefaults<3>::set_k(k);
  FunctionDefaults<3>::set_thresh(thresh);
  FunctionDefaults<3>::set_cubic_cell(-L, L);
  FunctionDefaults<3>::set_initial_level(4);
  FunctionDefaults<3>::set_truncate_on_project(false);
  FunctionDefaults<3>::set_bc(BC_FREE);

  //Print parameters
  print("k            ",k);
  print("thresh       ",thresh);
  print("rho_0        ", rho_0);
  print("beta         ",beta);
  print("epsilon      ",epsilon);
  print("maxiter      ",maxiter);
  print("L            ",L);
  TIME("make eletronic density:  ",const realfunc rho = real_factory_3d(world).f(rho_func));
  TIME("make nuclear density:  ",const realfunc rhon = real_factory_3d(world).f(nuc_func));
  print("total charge: ", rho.trace());
  print("total charge: ", rhon.trace());
  //total charge
  const realfunc rhot = rho + rhon;
  //Make a GygiPot object
  GygiPot gygi(rho,rhot,rho_0,epsilon,maxiter,world,Gamma);
  TIME("make ESP:                         ", realfunc U = gygi.ESP(world));
  //derivative of dielectric function w.r.t rho
  TIME("make depsilon_drho:               ", realfunc depdrho = gygi.depsilon_drho());
  //total charge
  TIME("make epsilon[rho]:                ",realfunc eprho = gygi.epsilon_rho());
  TIME("make reciprocal of dielectric:    ",realfunc reprho = gygi.re_epsilon_rho());
  TIME("make charge:                      ",realfunc charge = gygi.make_surfcharge(U));
  TIME("make V_epsilon:                   ",realfunc Ve = gygi.V_epsilon(U));
  TIME("make dESP_drho:                   ",realfunc desp_drho = gygi.dESP_drho(U));
  TIME("make cavitation energy:           ",double cavE = gygi.cavitation_energy());
  //TIME("make exponent surface:            ",double expsurf = gygi.expnt_surface());
  TIME("make deriv of free energy:        ",realfunc dfree = gygi.dfree_drho());
  //make the dielectric functional
  print("plotting ");
  coord_3d lo(0.0), hi(0.0); // Range for line plotting
  lo[0] =-10.0;
  hi[0] = 10.0;
  plot_line("ESP.dat", 401, lo, hi, U);
  plot_line("dielec_cavity.dat", 401, lo, hi, eprho);
  plot_line("surfcharge.dat", 401, lo, hi, charge);
  plot_line("Ueps.dat", 401, lo,hi, Ve );
  plot_line("rhot.dat", 401, lo, hi,rhot );
  plot_line("desp_drho.dat", 401, lo, hi, desp_drho );
  plot_line("dfree_drho.dat", 401, lo, hi, dfree );
  plotdx(eprho,"dielec_cavity.dx");
  print("\ncompute quantum surface \n");
  print(" dfree_drho  ",dfree.norm2());
  print("cav energy:   ",cavE);
  finalize();
  return 0;
}
