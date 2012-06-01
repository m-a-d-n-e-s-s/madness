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
  \file trunk/src/apps/jacob/abinitdftsolventsolver.h       
  \brief ab initio computation of the solvent-solute interaction potential
  \defgroup examplegygi compute the dielectric cavity and the electrostatic potential of solute in solvent
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
  const realfunc& rho;    //electronic density 
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
 //Binary operator to multiply grad U and depsilondrho
  struct Bop {
      void operator()(const Key<3>& key,
                      real_tensor U,
                      const real_tensor& gradu,
                      const real_tensor& dedrho) const {
            ITERATOR(U,
                     double d = gradu(IND);
                     double p = dedrho(IND);
                     if (std::abs(p)<cutrho || std::abs(d)<cutrho)
                         U(IND) = 0.0;
                     else
                         U(IND) = d*d*p;
                     );
        }

        template <typename Archive>
        void serialize(Archive& ar) {}
    };

  //Compute the dielectric cavity on the go
  template<typename T,int DIM>
  struct Epsilon_rho {
    double rho0;
    double beta;
    double epsilon;
    double cutrho;

    Epsilon_rho() {};

    Epsilon_rho(double rho0, double beta, double epsilon, double cutrho)
      : rho0(rho0), beta(beta), epsilon(epsilon), cutrho(cutrho)
    {}

    void operator()(const Key<DIM>& key,Tensor<T>& t) const {
      UNARY_OPTIMIZED_ITERATOR(T,t,
                               T rho = std::fabs(*_p0);
                               if(rho < cutrho)
				 *_p0 = epsilon;
			       else {
				 T ratio = std::pow(rho/rho0,2.0*beta);
				 T result  = (epsilon + ratio)/(1.0 + ratio);
				 //if (result > epsilon) {
				 //  print("!!!", *_p0, rho0, beta, epsilon, ratio, result);
				 //  throw "done";
                                 // }
				 *_p0 = result;
			       }
			       );
    }

    template <typename Archive>void serialize(Archive& ar) {
        ar & rho0 & beta & epsilon & cutrho;
    }
  };
  //Compute the normalization constant of  the density on the go
  template<typename T,int DIM>
  struct normconst{
    double rho0;
    double beta;
    double epsilon;
    double cutrho;
      normconst(){}
      normconst(double rho0, double beta, double epsilon, double cutrho)
          : rho0(rho0), beta(beta), epsilon(epsilon), cutrho(cutrho)
      {}

      void operator()(const Key<DIM>& key,Tensor<T>& t) const {
          UNARY_OPTIMIZED_ITERATOR(T,t,
                                   T rho = std::fabs(*_p0);
                                   if(rho < cutrho)
                                       *_p0 = 0.0;
                                   else {
                                       T ratio = std::pow(rho/rho0,2.0*beta);
                                       *_p0 = (epsilon - 1.0)/(epsilon + ratio);
                                   }
                                   );
      }
      template <typename Archive>void serialize(Archive& ar) {
          ar & rho0 & beta & epsilon & cutrho;
      }
  };
  //Compute the derivative of the dielectric cavity on the go
  template<typename T,int DIM>
  struct dEpsilon_drho{
    double rho0;
    double beta;
    double epsilon;
    double cutrho;
    dEpsilon_drho(){}
    dEpsilon_drho(double rho0, double beta, double epsilon, double cutrho)
      : rho0(rho0), beta(beta), epsilon(epsilon), cutrho(cutrho)
    {}

    void operator()(const Key<DIM>& key,Tensor<T>& t) const {
      UNARY_OPTIMIZED_ITERATOR(T,t,
                               T rho = std::fabs(*_p0);
			       if(rho < cutrho)
				*_p0 = 0.0;
			       else {
				 double fac = (1.0 - epsilon)/rho0;
				 T ratone = std::pow(rho/rho0,2.0*beta -1.0);
				 T ratio = std::pow(rho/rho0,2.0*beta);
				 *_p0 = (fac*2.0*beta*ratone)/((1.0 + ratio)*(1.0 + ratio));
			       }
			       );
    }
    template <typename Archive>void serialize(Archive& ar) {
        ar & rho0 & beta & epsilon & cutrho;
    }
  };
    //Compute the surface of the cavity on the go
  template<typename T,int DIM>
  struct dielectric_surface{
    double rho0;
    double beta;
    double epsilon;
    double cutrho;
      dielectric_surface(){}
      dielectric_surface(double rho0, double beta, double epsilon, double cutrho)
          : rho0(rho0), beta(beta), epsilon(epsilon), cutrho(cutrho)
      {}

    void operator()(const Key<DIM>& key,Tensor<T>& t) const {
      UNARY_OPTIMIZED_ITERATOR(T,t,
                               T rho = std::abs(*_p0)-rho0;
			       if(rho < cutrho)
				*_p0 = 0.0;
			       else {
                                   double fac = 2.0*beta;
                                   //T ratone = std::pow(rho/rho0,2.0*beta-1);
                                   T ratio = std::pow(rho,2.0*beta);
                                   *_p0 = ((fac*ratio)/((1.0 + ratio)*rho*(1.0 + ratio)));
                               }
                               );
    }
      template<typename Archive>void serialize(Archive& ar) {
          ar & rho0 & beta & epsilon & cutrho;
    }
  };
  //Compute the ratio of the derivative of epsilon by epsilon on the go
    template<typename T,int DIM>
    struct ratioepsilon{
        double rho0;
        double beta;
        double epsilon;
        double cutrho;
        ratioepsilon(){}
        ratioepsilon(double rho0, double beta, double epsilon, double cutrho)
            : rho0(rho0), beta(beta), epsilon(epsilon), cutrho(cutrho)
        {}
        
        void operator()(const Key<DIM>& key,Tensor<T>& t) const {
            UNARY_OPTIMIZED_ITERATOR(T,t,
                                     T rho = std::fabs(*_p0);
                                     if(rho < cutrho)
                                         *_p0 = 0.0;
                                     else {
                                         //   double fac = (1.0 - epsilon)/rho0;
                                         //   T ratone = std::pow(std::fabs(*_p0)/rho0,2.0*beta -1.0);
                                         //   
                                         // *_p0 = (fac*2.0*beta*ratone)/((1.0 + ratio)*(epsilon + ratio));
                                         
                                         T r = std::pow(rho/rho0,2.0*beta);
                                         T rs = std::pow(rho/rho0,4.0*beta);
                                         *_p0 = 2.0*beta*(1.0-epsilon)*(r + rs)/(rho*(1+r)*(1+r)*(epsilon+r));
                                     }
                                     );
        }
        template <typename Archive>void serialize(Archive& ar) {
            ar & rho0 & beta & epsilon & cutrho;
        }
    };
    //Compute the reciprocal of dielectric cavity on the go
  template<typename T,int DIM>
  struct repsilon_rho {
    double rho0;
    double beta;
    double epsilon;
    double cutrho;
      repsilon_rho(){}
    repsilon_rho(double rho0, double beta, double epsilon, double cutrho)
      : rho0(rho0), beta(beta), epsilon(epsilon), cutrho(cutrho)
    {}

    void operator()(const Key<DIM>& key,Tensor<T>& t) const {
      UNARY_OPTIMIZED_ITERATOR(T,t,
                               T rho = std::fabs(*_p0);
			       if(rho < cutrho)
				*_p0 = 1.0/epsilon;
			       else {
				 T ratio = std::pow(rho/rho0,2*beta);
				 *_p0 = (1.0 + ratio)/(epsilon + ratio);
			       }
			       );
    }
    template <typename Archive>void serialize(Archive& ar) {
        ar & rho0 & beta & epsilon & cutrho;
    }
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
    realfunc value = copy(rho);
    value.unaryop(Epsilon_rho<double,3>(rho_0, beta, epsilon,cutrho));
    return value;
  }
  //make normalization constant
  realfunc make_normconst() const {
    realfunc value = copy(rho);
    value.unaryop(normconst<double,3>(rho_0, beta, epsilon,cutrho));
    return value;
  }
  //make dielectric surface
  realfunc make_surface() const {
      realfunc value = copy(rho);
      value.unaryop(dielectric_surface<double,3>(rho_0, beta,epsilon,cutrho));
    return value;
  }
  //make reciprocal of epsilon
  realfunc make_repsilon() const {
    realfunc value = copy(rho);
    value.unaryop(repsilon_rho<double,3>(rho_0, beta, epsilon,cutrho));
    return value;
  }

  //make derivative of epsilon
  realfunc make_depsilon_drho() const {
    realfunc value = copy(rho);
    value.unaryop(dEpsilon_drho<double,3>(rho_0, beta,epsilon,cutrho));
    return value;
  }
  //make ratio of the derivative of epsilon w.r.t rho by epsilon
  realfunc make_ratioepsilon() const {
    realfunc value = copy(rho);
    value.unaryop(ratioepsilon<double,3>(rho_0, beta,epsilon,cutrho));
    return value;
  }

  //cavitation energy
  double cavitation_energy() const {
      double quantum_surface =(make_surface()).norm2()*(grad_of(rho).norm2());
      double convfact = 6.423049507e-4; // 1N/m = 6.423049507e−4a.u 
      return convfact*Gamma*quantum_surface;
  }
 
    //cavitation potential to be included in the KS equation
  realfunc dfree_drho() const {
      double rnorm = 1.0/grad_of(rho).norm2();
      double fac = 2.0*beta*Gamma;
      real_derivative_3d Dx = free_space_derivative<double,3>(rho.world(), 0);
      real_derivative_3d Dy = free_space_derivative<double,3>(rho.world(), 1);
      real_derivative_3d Dz = free_space_derivative<double,3>(rho.world(), 2);
      realfunc gradxrho = Dx(rho);
      realfunc gradyrho = Dy(rho);
      realfunc gradzrho = Dz(rho);
      realfunc gxGgx = gradxrho*(Dx(gradxrho) + Dx(gradyrho) + Dx(gradzrho));
      realfunc gyGgy = gradyrho*(Dy(gradxrho) + Dy(gradyrho) + Dy(gradzrho));
      realfunc gzGgz = gradzrho*(Dz(gradxrho) + Dz(gradyrho) + Dz(gradzrho));
      realfunc gradnormrho = (gxGgx + gyGgy + gzGgz).scale(rnorm);
      realfunc scoef = (make_surface().norm2())*re_grad_rho(rho);
      realfunc ddelG = fac*scoef*gradnormrho;
      double convfact = 6.423049507e-4; // 1N/m = 6.423049507e−4a.u 
      return ddelG.scale(convfact);
  }
    //compute the gradient of epsilon[rho] 
    realfunc depsilon_dr() const {
        real_derivative_3d Dx = free_space_derivative<double,3>(rho.world(), 0);
        real_derivative_3d Dy = free_space_derivative<double,3>(rho.world(), 1);
        real_derivative_3d Dz = free_space_derivative<double,3>(rho.world(), 2);
        realfunc grad = (Dx(rho) + Dy(rho) + Dz(rho));
        realfunc depdrho = copy(rho).unaryop(dEpsilon_drho<double,3>(rho_0, beta,epsilon,cutrho));
        //    return make_depsilon_drho()*grad_of(rhot);
        return grad*depdrho;
  }
  //compute the surface charge                                                                                                                            
  realfunc make_surfcharge(const realfunc& u) const {
    real_derivative_3d Dx = free_space_derivative<double,3>(rho.world(), 0);
    real_derivative_3d Dy = free_space_derivative<double,3>(rho.world(), 1);
    real_derivative_3d Dz = free_space_derivative<double,3>(rho.world(), 2);
    realfunc pgrad = (Dx(rho)*Dx(u) + Dy(rho)*Dy(u) + Dz(rho)*Dz(u));
    //    realfunc gradU = grad_of(u);
    const double rfourpi = 1.0/(4.0*constants::pi);
    return (make_ratioepsilon()*pgrad).scale(rfourpi);
  }
  
  //Define the electrostatic potential
  realfunc ESP()const {
    const bool USE_SOLVER = true;
    double tol = std::max(1e-7,FunctionDefaults<3>::get_thresh());
    realfunc charge = make_repsilon()*rhot;
    realfunc tcharge = make_normconst()*rhot; //total molecular charge in solvent
    realfunc U0 = op(charge);  //U
    //    double einf = -1.0/epsilon;
    realfunc U = op(rhot); //Uvac
    realfunc Ug = op(tcharge);//U0
    realfunc Ur = U;// - Uvac; 
    double unorm = U.norm2();
    //print("U.norm2: ", unorm);
    /*coord_3d lo(0.0), hi(0.0);
    lo[0] = -20.0;
    hi[0] = 20.0;
    plot_line("iso_rho.dat", 10001, hi, lo,rho);
    plot_line("iso_dfree_drho.dat", 10001, hi, lo,dfree_drho());
    plot_line("potential.dat", 10001, hi, lo,U0);
    plot_line("iso_epsilon.dat", 10001, hi, lo,make_epsilon());
    plotdx(make_epsilon(),"iso_epsilon_pot.dx");
    plot_line("iso_ratioepsilon.dat", 10001, hi, lo,make_ratioepsilon());
    //plot_line("repsilon.dat", 10001, hi, lo,make_repsilon());
    plot_line("iso_depsilondrho.dat", 10001, hi, lo,make_depsilon_drho());
    plot_line("iso_depsilondr.dat", 10001, hi, lo,depsilon_dr());
    //    plot_line("iso_pot.dat", 10001, hi, lo,U);
    plot_line("iso_surfacecharge.dat", 10001, hi, lo,make_surfcharge(U));
    // throw "done";*/
    if (USE_SOLVER) {
        madness::NonlinearSolver solver;//(5);
      // This section employs a non-linear equation solver from solvers.h                                                                                  
      //  http://onlinelibrary.wiley.com/doi/10.1002/jcc.10108/abstract                                                                               
      if (world.rank() == 0){
        print("\n\n");//for formating output                                                                                                              
	madness::print("    Computing the Electrostatic Potential   ");
	madness::print("           ______________________           \n ");
	
	madness::print("iteration          residue norm2            soln(10.0)  ");
      }

      for (int iter=0; iter<maxiter; iter++) {
          realfunc uvec = U;
          realfunc Scharge = make_surfcharge(U);// - Uvac);
          realfunc rvec = (U - U0 - op(Scharge)).truncate();
          //plotdx(op(Scharge),"iso_surface_pot.dx");
          realfunc U_new = solver.update(uvec,rvec);
          double err = rvec.norm2();
          if (world.rank()==0)
              // madness::print("  ", iter,"             " , err,"           ",U(coord_3d(10.0)));
              std::printf("%8d %22.10f %22.10f \n", iter,err,Ur(coord_3d(10.0)));
          if (err >0.3*unorm) U = 0.5*U + 0.5*U_new;
          //if (err >0.3*unorm) Ur = 0.5*Ur + 0.5*U_new;
          else
              U = U_new;
          //Ur = U_new;
          if(err < 10.0*tol) break;
      }
    }
    //plot_line("iso_total_pot.dat", 10001, hi, lo,U);
    //plotdx(U,"iso_total_pot.dx");
    // throw "done";
    realfunc rxtnpot = U - op(rhot);
    return rxtnpot;
  }
    
  //Defining the polarization of the dielectric continuum in the presence of 
 // an external electric field. Used in the response of solvated molecule  
    realfunc Laplace_ESP(const realfunc& uguess)const {
        const bool USE_SOLVER = true;
        double tol = std::max(1e-7,FunctionDefaults<3>::get_thresh());
        realfunc U = uguess;
        double unorm = U.norm2();
        print("U.norm2: ", unorm);
        //start plots
        coord_3d lo(0.0), hi(0.0);
        lo[0] = -20.0;
        hi[0] = 20.0;
        if (USE_SOLVER) {
            madness::NonlinearSolver solver;//(5);
            if (world.rank() == 0){
                print("\n\n");//for formating output
                madness::print("    Computing the Continuum-Field Interaction Potential   ");
                madness::print("           ______________________           \n ");
                
                madness::print("iteration          residue norm2            soln(10.0)  ");
            }
            
            for (int iter=0; iter<maxiter; iter++) {
                realfunc uvec = U;
                realfunc Scharge = make_surfcharge(U);// - Uvac);
                realfunc rvec = (U - op(Scharge)).truncate();
                plot_line("continuum_surface_pot.dat", 10001, hi, lo,op(Scharge));
                plotdx(op(Scharge),"continuum_surface_pot.dx");
                realfunc U_new = solver.update(uvec,rvec);
                double err = rvec.norm2();
                if (world.rank()==0)
                    std::printf("%8d %22.10f %22.10f \n", iter,err,U(coord_3d(10.0)));
                //madness::print("  ", iter,"             " , err,"           ",U(coord_3d(10.0)));
                if (err >0.3*unorm) U = 0.5*U + 0.5*U_new;
                else
                    U = U_new;
                if(err < 10.0*tol) break;
            }
        }
        // plot_line("continuum_field_pot.dat", 10001, hi, lo,U);
        //plotdx(U,"continuum_field_pot.dx");
        return U;
    }
//Defining the derivative of the ESP w.r.t rho
// this function is very noisy and is not called
  realfunc dESP_drho(const realfunc& u)const {
    realfunc gu = grad_of(u);
    realfunc dedrho= make_depsilon_drho();
    realfunc ues = binary_op(gu,dedrho, Bop());
    double coef = -1.0/(8.0*constants::pi);
    realfunc dep = u + ues.scale(coef);
    realfunc gusq = gu*gu;;
    coord_3d lo(0.0), hi(0.0);
    lo[0] = -20.0;
    hi[0] = 20.0;
    plot_line("iso_square_gradU.dat", 10001, hi, lo,gusq);
    return dep;
  }
    //computes components of the the electric field due to the surface charge 
    realfunc make_electric_field(const realfunc& u) const {
        real_derivative_3d Dx = free_space_derivative<double,3>(u.world(), 0);
        real_derivative_3d Dy = free_space_derivative<double,3>(u.world(), 1);
        real_derivative_3d Dz = free_space_derivative<double,3>(u.world(), 2);
        double fac = -1.0/(4.0*constants::pi);
        realfunc Sigma =(Dx(u) + Dy(u) + Dz(u)).scale(fac); //excess charge on colloid surface
        realfunc uxc  = op(Sigma); //coulomb potential due to excess charge                                                                         
        realfunc dx = Dx(uxc) ;
        realfunc dy = Dy(uxc) ;
        realfunc dz = Dz(uxc) ;
        return (dx + dy + dz).scale(-1.0);
    }
    //calculate the average reaction field(\int C(r)F_r(r)d \tau/\int C(r)d\tau  
    //the mask is that of the molecule because the average field is that felt by the molecule 
    double ave_rxn_field(const real_function_3d& F_r,const real_function_3d& mask)const {
        real_function_3d  pdt = mask*F_r;
        double numerator = pdt.trace();
        double denominator = mask.trace();
        return (numerator/denominator);
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
