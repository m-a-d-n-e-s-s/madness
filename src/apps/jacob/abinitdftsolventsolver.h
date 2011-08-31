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
                                                                                                                    
  The source is <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local/trunk/src/apps/jacob/abinitdftsolventsolver>here</a>.  

  \par Points of interest
  - compute the dielectric functional (of density)
  - compute the electrostatic potential by convolving the free space Green's function with the effective charge and the induced surface charge
  - compute the derivatives of the dielectric functional and the electrosstatic potential with respect to the density
  \par Background 
  - This program is an implementation of the solvation model in THE JOURNAL OF CHEMICAL PHYSICS 124, 074103 (2006) 
*/
#ifndef MADNESS_ABINITDFTSOLVENTSOLVER_H
#define MADNESS_ABINITDFTSOLVENTSOLVER_H
#include <mra/mra.h>
#include <constants.h>
#include <ctime>
#include <linalg/solvers.h>
#include <mra/funcplot.h>
#include <examples/nonlinsol.h>

typedef real_function_3d realfunc;

//define a class for the Gygi potential and the surface cavity
class DFTSolventSolver{
private:
  const realfunc& rho;    //density for the dielectric functional
  const realfunc& rhot;    //density for the dielectric functional
  const double& rho_0;  // the density at the solut-solvent interface 
  const double& epsilon; //assyptotic value of the dielectric functional
  const int& maxiter;  // maximum iterations in the solver
  World& world;
  const double& Gamma;   //the surface tention of the solvent
  const double& beta;
  const double& minlen; // minimul in the coulomb operator
  const double& thresh;
  real_convolution_3d op; //< Coulomb operator
  static const double cutrho = 1e-12;  //cutoff value of the density
  //utility function
  
  //unary operator to determine the reciprocal of a madness function
  template<typename T,int DIM>
  struct Reciprocal {
    void operator()(const Key<DIM>& key, Tensor<T>& t) const {
      UNARY_OPTIMIZED_ITERATOR(T, t,*_p0=1.0/(*_p0));
    }
    template <typename Archive>void serialize(Archive& ar) {}
  };
  //Compute the dielectric cavity on the go
  template<typename T,int DIM>
  struct Epsilon_rho {
    double rho0;
    double beta;
    double epsilon;
    double cutrho;
    Epsilon_rho(double rho0, double beta, double epsilon, double cutrho)
      : rho0(rho0), beta(beta), epsilon(epsilon), cutrho(cutrho)
    {}

    void operator()(const Key<DIM>& key,Tensor<T>& t) const {
      UNARY_OPTIMIZED_ITERATOR(T,t,
			       if(std::fabs(*_p0) < cutrho)
				 *_p0 = epsilon;
			       else {
				 T ratio = std::pow(std::fabs(*_p0)/rho0,2.0*beta);
				 *_p0 = (epsilon + ratio)/(1.0 + ratio);
			       }
			       );
    }
    template <typename Archive>void serialize(Archive& ar) {}
  };
  //Compute the normalization constant of  the density on the go
  template<typename T,int DIM>
  struct normconst{
    double rho0;
    double beta;
    double epsilon;
    double cutrho;
    normconst(double rho0, double beta, double epsilon, double cutrho)
      : rho0(rho0), beta(beta), epsilon(epsilon), cutrho(cutrho)
    {}

    void operator()(const Key<DIM>& key,Tensor<T>& t) const {
      UNARY_OPTIMIZED_ITERATOR(T,t,
			       if(std::fabs(*_p0) < cutrho)
				 *_p0 = 0.0;
			       else {
				 T ratio = std::pow(std::fabs(*_p0)/rho0,2.0*beta);
				 *_p0 = (epsilon - 1.0)/(epsilon + ratio);
			       }
			       );
    }
    template <typename Archive>void serialize(Archive& ar) {}
  };
  //Compute the derivative of the dielectric cavity on the go
  template<typename T,int DIM>
  struct dEpsilon_drho{
    double rho0;
    double beta;
    double epsilon;
    double cutrho;
    dEpsilon_drho(double rho0, double beta, double epsilon, double cutrho)
      : rho0(rho0), beta(beta), epsilon(epsilon), cutrho(cutrho)
    {}

    void operator()(const Key<DIM>& key,Tensor<T>& t) const {
      UNARY_OPTIMIZED_ITERATOR(T,t,
			       if(std::fabs(*_p0) < cutrho)
				*_p0 = 0.0;
			       else {
				 double fac = (1.0 - epsilon)/rho0;
				 T ratone = std::pow(std::fabs(*_p0)/rho0,2.0*beta -1.0);
				 T ratio = std::pow(std::fabs(*_p0)/rho0,2.0*beta);
				 *_p0 = fac*2.0*beta*ratone/((1.0 + ratio)*(1.0 + ratio));
			       }
			       );
    }
    template <typename Archive>void serialize(Archive& ar) {}
  };
  //Compute the ratio of the derivative of epsilon by epsilon on the go
  template<typename T,int DIM>
  struct ratioepsilon{
    double rho0;
    double beta;
    double epsilon;
    double cutrho;
    ratioepsilon(double rho0, double beta, double epsilon, double cutrho)
      : rho0(rho0), beta(beta), epsilon(epsilon), cutrho(cutrho)
    {}

    void operator()(const Key<DIM>& key,Tensor<T>& t) const {
      UNARY_OPTIMIZED_ITERATOR(T,t,
			       if(std::fabs(*_p0) < cutrho)
				 *_p0 = 0.0;
			       else {
				 double fac = (1.0 - epsilon)/rho0;
				 T ratone = std::pow(std::fabs(*_p0)/rho0,2.0*beta -1.0);
				 T ratio = std::pow(std::fabs(*_p0)/rho0,2.0*beta);
			       *_p0 = (fac*2.0*beta*ratone)/((1.0 + ratio)*(epsilon + ratio));
			       }
			       );
    }
    template <typename Archive>void serialize(Archive& ar) {}
  };
  //Compute the reciprocal of dielectric cavity on the go
  template<typename T,int DIM>
  struct repsilon_rho {
    double rho0;
    double beta;
    double epsilon;
    double cutrho;
    repsilon_rho(double rho0, double beta, double epsilon, double cutrho)
      : rho0(rho0), beta(beta), epsilon(epsilon), cutrho(cutrho)
    {}

    void operator()(const Key<DIM>& key,Tensor<T>& t) const {
      UNARY_OPTIMIZED_ITERATOR(T,t,
			       if(std::fabs(*_p0) < cutrho)
				*_p0 = 1.0/epsilon;
			       else {
				 T ratio = std::pow(std::fabs(*_p0)/rho0,2*beta);
				 *_p0 = (1.0 + ratio)/(epsilon + ratio);
			       }
			       );
    }
    template <typename Archive>void serialize(Archive& ar) {}
  };
  //gradient of density
  realfunc grad_of(const realfunc& dens) const {
    real_derivative_3d Dx = free_space_derivative<double,3>(dens.world(), 0);
    real_derivative_3d Dy = free_space_derivative<double,3>(dens.world(), 1);
    real_derivative_3d Dz = free_space_derivative<double,3>(dens.world(), 2);
    return (Dx(dens) + Dy(dens) + Dz(dens));
  }

  //reciprocal of gradient of density
  realfunc re_grad_rho(const realfunc& dens) const{
    realfunc value = grad_of(dens);
    value.unaryop(Reciprocal<double,3>());
    return value;
  }
public:
  //make epsilon
  realfunc make_epsilon() const {
    realfunc value = copy(rhot);
    value.unaryop(Epsilon_rho<double,3>(rho_0, beta, epsilon,cutrho));
    return value;
  }
  //make normalization constant
  realfunc make_normconst() const {
    realfunc value = copy(rhot);
    value.unaryop(normconst<double,3>(rho_0, beta, epsilon,cutrho));
    return value; 
  }
  //make reciprocal of epsilon
  realfunc make_repsilon() const {
    realfunc value = copy(rhot);
    value.unaryop(repsilon_rho<double,3>(rho_0, beta, epsilon,cutrho));
    return value;
  }
  
  //make derivative of epsilon
  realfunc make_depsilon_drho() const {
    realfunc value = copy(rhot);
    value.unaryop(dEpsilon_drho<double,3>(rho_0, beta,epsilon,cutrho));
    return value;
  }
  //make ratio of the derivative of epsilon w.r.t rho by epsilon
  realfunc make_ratioepsilon() const {
    realfunc value = copy(rhot);
    value.unaryop(ratioepsilon<double,3>(rho_0, beta,epsilon,cutrho));
    return value; 
  }
  //cavitation energy
  double cavitation_energy() const {
    double fac =1.0/rho_0;
    double quantum_surface = 0.25*fac*beta*(grad_of(rhot).norm2());
    double convfact = 6.423049507e-4; // 1N/m = 6.423049507e−4a.u 
    return convfact*Gamma*quantum_surface;
  }
  /*
  //cavitation potential to be included in the KS equation
  realfunc dfree_drho() const {
    double fac = beta()/(4.0*rho_0);
    real_derivative_3d Dx = free_space_derivative<double,3>(rhot.world(), 0);
    real_derivative_3d Dy = free_space_derivative<double,3>(rhot.world(), 1);
    real_derivative_3d Dz = free_space_derivative<double,3>(rhot.world(), 2);
    realfunc gradsq = Dx(rhot)*Dx(rhot) + Dy(rhot)*Dy(rhot) + Dz(rhot)*Dz(rhot);
    realfunc nume = fac*grad_of(gradsq);
    realfunc denom =(1.0/grad_of(rhot).norm2())*re_grad_rho(rhot);
    double convfact = 6.423049507e-4; // 1N/m = 6.423049507e−4a.u 
    return nume*denom*convfact*Gamma;
  }
  // defining the quantum surface with a Gaussian
  double expnt_surface() const{
    real_derivative_3d Dx = free_space_derivative<double,3>(rhot.world(), 0);
    real_derivative_3d Dy = free_space_derivative<double,3>(rhot.world(), 1);
    real_derivative_3d Dz = free_space_derivative<double,3>(rhot.world(), 2);
    realfunc grad_rho = Dx(rhot) + Dy(rhot) + Dz(rhot);
    double fac = 3.0/(rho_0*sqrt(constants::pi));
    return fac*grad_rho.norm2();
  }
  */
//Define the guess potential for the Poisson's equation
  realfunc GuessPotential() const {
    return op(make_repsilon()*rhot);   //U_0  
}
//Molecular potential i.e potential due to the molecular charge distribution                                                                             
  realfunc MolecularPotential()const {
    return op(rhot);
  }
  //compute the surface charge                                                                                                                            
  realfunc make_surfcharge(const realfunc& u) const {
    real_derivative_3d Dx = free_space_derivative<double,3>(rhot.world(), 0);
    real_derivative_3d Dy = free_space_derivative<double,3>(rhot.world(), 1);
    real_derivative_3d Dz = free_space_derivative<double,3>(rhot.world(), 2);
    real_derivative_3d Dxx = free_space_derivative<double,3>(u.world(), 0);
    real_derivative_3d Dyy = free_space_derivative<double,3>(u.world(), 1);
    real_derivative_3d Dzz = free_space_derivative<double,3>(u.world(), 2);
    realfunc pgrad = (Dx(rhot)*Dxx(u) + Dy(rhot)*Dyy(u) + Dz(rhot)*Dzz(u));
    //    realfunc gradU = grad_of(u);
    const double rfourpi = 1.0/(4.0*constants::pi);
    return (rfourpi*make_ratioepsilon()*pgrad).truncate();
  }
  
  //Define the electrostatic potential
  realfunc ESP()const {
    const bool USE_SOLVER = true;
    double tol = std::max(1e-3,FunctionDefaults<3>::get_thresh());
    realfunc charge = make_repsilon()*rhot;
    realfunc tcharge = make_normconst()*rhot; //total molecular charge in solvent
    realfunc U = op(charge);
    realfunc Uvac = op(rhot);
    realfunc U0 = op(tcharge);
    realfunc Ur = U - Uvac; 
    double unorm = Ur.norm2();
    //    print("U.norm2: ", unorm);
    //start plots
    coord_3d lo(0.0), hi(0.0);
    lo[0] = -20.0;
    hi[0] = 20.0;
    plot_line("epsilon.dat", 401, hi, lo,make_epsilon());
    plot_line("ratioepsilon.dat", 401, hi, lo,make_ratioepsilon());
    // plot_line("repsilon.dat", 401, hi, lo,make_repsilon());
    // plot_line("depsilondrho.dat", 401, hi, lo,make_depsilon_drho());
    // plot_line("depsilondr.dat", 401, hi, lo,depsilon_dr());
    // plot_line("guesspot.dat", 401, hi, lo,U0);
    // plot_line("surfacecharge.dat", 401, hi, lo,make_surfcharge(U));
    //end plots
    if (USE_SOLVER) {
      madness::NonlinearSolver solver;
      // This section employs a non-linear equation solver from solvers.h                                                                                  
      //  http://onlinelibrary.wiley.com/doi/10.1002/jcc.10108/abstract                                                                               
      if (world.rank() == 0){
        print("\n\n");//for formating output                                                                                                              
	madness::print("    Computing the Electrostatic Potential   ");
	madness::print("           ______________________           \n ");
	
	madness::print("iteration          residue norm2            soln(10.0)  ");
      }
      realfunc uvec, rvec;
      for (int iter=0; iter<maxiter; iter++) {
        uvec = Ur;
        realfunc Scharge = make_surfcharge(Ur + Uvac);
        rvec = (Ur - U0 + op(Scharge)).truncate();
	plot_line("surfacepot.dat", 401, hi, lo,op(Scharge));
        realfunc U_new = solver.update(uvec,rvec);
        double err = rvec.norm2();
        if (world.rank()==0)
	  madness::print("  ", iter,"             " , err,"           ",Ur(coord_3d(10.0)));
        if (err >0.3*unorm) Ur = 0.5*Ur + 0.5*U_new;
	else
	  Ur = U_new;
        if(err < 10.0*tol) break;
      }
    }
    return Ur;
  }
  //Define the dielectric potential
  realfunc V_epsilon(const realfunc& u) const {
    realfunc grad_ESP = grad_of(u);
    double coef = -1.0/(8.0*constants::pi);
    return coef*grad_ESP*grad_ESP*make_depsilon_drho();
  }
  //Defining the gradient of the ESP w.r.t rho
  realfunc dESP_drho(const realfunc u) const {
    realfunc dep = u +V_epsilon(u);
    return dep;
  }
  //Defining the Constructor
 DFTSolventSolver(const realfunc& rho,
		  const realfunc& rhot,
		  const double& rho_0,
		  const double& epsilon,
		  const int& maxiter,
		  World& world,
		  const double& Gamma,
		  const double& beta,
		  const double& minlen)
   :rho(rho),
    rhot(rhot),
    rho_0(rho_0),
    epsilon(epsilon),
    maxiter(maxiter),
    world(world),
    Gamma(Gamma),
    beta(beta),
    minlen(minlen),
    thresh(FunctionDefaults<3>::get_thresh()),
    op(CoulombOperator(world, minlen, thresh)){}
};
#endif
