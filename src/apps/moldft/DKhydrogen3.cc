//Solves for the ground state energy of the hydrogenic atom
//Using first-order Douglas Kroll
//
//TODO: 
//1. Update operators. Current ones are only accurate to ~1e-10 (fixed, I think)
//2. Update rele, calculating the expectation value of the Hamiltonian: Using the straight potential is wrong.

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/constants.h>
#include <stdlib.h> //for atoi
#include <vector>
#include "DKops.h"

using namespace madness;

static const long k = 8;             //wavelet order in madness
static const double thresh = 1e-6;   //desired precision when refining
static const double speedoflight = 137.0359895; //speed of light
static const double PI = constants::pi;  //mmm, pie
static const double small = 1e-8; // shortest length scale we want to resolve

static double ttt;
static void START_TIMER(World& world){
     world.gop.fence();
     ttt=wall_time();
}

static void END_TIMER(World& world, char* message){
     world.gop.fence();
     ttt = wall_time() - ttt;
     if(world.rank()==0) print(message, ttt);
}

//Functor to make the nuclear potential, slightly perturbed to remove singularity
class PotentialFunctor : public FunctionFunctorInterface<double,3> {
     private:
          int m_Z;
     public:
          PotentialFunctor(int Z){
               m_Z = Z;
          }
          double operator() (const coord_3d&r) const {
               return -m_Z/(sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]+small*small)); 
          }
};

//norelativistic ground state functor
class NRGFunctor : public FunctionFunctorInterface<double,3> {
     private:
          int m_Z;

     public:
          NRGFunctor(int Z){
               m_Z = Z;
          }

          double operator() (const coord_3d& r) const {
               return std::exp(-m_Z*std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]+small*small)); 
          }

};

double exact_energy(double Z) {
     double c = speedoflight;
     double d = sqrt(1 - Z*Z/(c*c));
     return c*c / sqrt(1.0 + Z*Z/(c*c*d*d)) - c*c; // may lose a few digits
}


double rele(World& world, const real_function_3d& psi, const real_function_3d Vpsi, real_convolution_3d& Pbarop, real_function_3d& Vnuc){
     if(world.rank() == 0) print("Now in rele");
     real_function_3d result = real_factory_3d(world);
     double answer(0.0);
     for(int axis = 0; axis < 3; axis++){
          real_derivative_3d D = free_space_derivative<double,3>(world,axis);
          real_function_3d dpsi = D(psi);
          dpsi = D(dpsi);
          result -= 0.5*dpsi;
     }
     
     answer = inner(psi,result)/inner(psi,psi);
     if(world.rank()==0) print("   Nonrelativistic KE: ", answer);
     
     result = apply(Pbarop, result);
     result.scale(2.0*speedoflight);
     answer = inner(psi,result)/inner(psi,psi);

     answer = answer/pow(2.0*PI,1.5); //I can't find this factor in my notes but I know it's here.....

     if(world.rank()==0) print("   Relativistic KE: ", answer);

     double PE = inner(psi, Vnuc*psi)/inner(psi,psi);

     if(world.rank()==0) print("Non-adjusted potential energy: ", PE);

     //double PE = inner(psi,Vpsi)/inner(psi,psi); //not really the potential energy

     //if(world.rank()==0) print("   Adjusted potential energy: ", PE);
     
     return answer + PE;
}

real_function_3d makerightside(World& world, real_function_3d& V, real_function_3d& psi, real_convolution_3d& PbarAop, real_convolution_3d& Aop){

     if(world.rank()==0) print("now in makerightside");
     real_function_3d Vpsi = V*psi;
     Vpsi.truncate();
     real_function_3d AVpsi = apply(Aop,Vpsi);
     real_function_3d VApsi = V*apply(Aop,psi);
     real_function_3d AVApsi = apply(Aop,VApsi);
     
     real_function_3d result = 4.0*pow(PI,3.0)*Vpsi + 2.0*pow(PI,1.5)*(AVpsi + VApsi) + AVApsi;
     //result = pow(2.0*PI,-1.5)*result;
     //real_function_3d result = 0.5*Vpsi + (1.0/sqrt(2.0))*(AVpsi + VApsi) + AVApsi;


     //double ttemp = inner(psi,result)/inner(psi,psi);
     //if(world.rank()==0) print("AVApsi expectation value:", ttemp);
     if(world.rank()==0) print(" made functions");

     result = pow(2*PI,-3.0)*result; //I can't find this factor in my notes either, but it has to be in here....
     result.truncate();

     real_function_3d temporary(world);
     for(int axis=0; axis<3; axis++){
          real_derivative_3d D = free_space_derivative<double,3>(world,axis);
          real_function_3d t = D(psi);
          t = apply(PbarAop,t);
          t = V*t;
          t = D(t);
          t = apply(PbarAop,t);
          temporary += t; 
     }
     temporary = pow(2.0*PI,-1.5)*temporary;

     //ttemp = inner(temporary,psi)/inner(psi,psi);
     //if(world.rank()==0) print("APVPApsi expectation value: ", ttemp);
     if(world.rank()==0) print("made result");

     result = result + temporary;

     //ttemp = result.norm2();
     //if(world.rank()==0) print("Made a rightside! Norm is: ", ttemp);

     return result;
}

double iterate(World& world, real_function_3d Vpsi, real_function_3d& psi, real_convolution_3d& Ebarop){

     Vpsi.truncate();  //Vpsi should already contain the result of makerightside
     Vpsi.scale(-1.0); //move Vpsi to the right hand side of the equation

     real_function_3d tmp = apply(Ebarop, Vpsi);
     double norm = tmp.norm2();
     tmp.scale(1.0/norm);
     real_function_3d residual = psi-tmp;
     double resid = residual.norm2();
     if(world.rank()==0) print("Residual: ", resid);
     psi = tmp;

     return resid;


}

int main(int argc, char** argv){
     //set up parallel environment (this is boiler-plate)
     initialize(argc, argv);
     World world(SafeMPI::COMM_WORLD);
     startup(world,argc,argv);
     std::cout.precision(10);
     int Z = atoi(argv[1]);
     
     //some of our madness parameters depend on Z:
     const double L = 46.0/Z;
     const double initialguessenergy = -Z*Z/2.0;


     //Set necessary madness parameters
     FunctionDefaults<3>::set_k(k);
     FunctionDefaults<3>::set_thresh(thresh);
     FunctionDefaults<3>::set_truncate_mode(1);
     FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);
     //FunctionDefaults<3>::set_initial_level(5);


     //Starting guess is the nonrelativistic exact solution
     NRGFunctor nonrelguess(Z);
     real_function_3d psi = real_factory_3d(world).functor(nonrelguess);
     double normconst = psi.norm2();
     psi.scale(1.0/normconst);

     //make potential
     PotentialFunctor myV(Z);
     real_function_3d Vnuc = real_factory_3d(world).functor(myV).truncate_mode(0);


     //make operators
     real_convolution_3d Pbarop = Pbar(world);
     real_convolution_3d Aop = A(world);
     real_convolution_3d PbarAop = PbarA(world);
     real_function_3d rightside = makerightside(world, Vnuc, psi, PbarAop, Aop);


     double energy = initialguessenergy;
     const double exact = exact_energy(Z);
     double residual = 1.0/thresh;
     int iter = 1;
     if(world.rank()==0) {
          print("\n            Box Width: ", L);
          print("        Wavelet Order: ", k);
          print(" Refinement Threshold: ", thresh);
          print("                    Z: ", Z);
          print("   Exact Dirac Energy: ", exact);
          print("      Starting Energy: ", energy, "\n");
     }

     
     //iterate solution
     //while(residual > sqrt(thresh)){
     while(residual > sqrt(thresh)/10.0 or iter < 4){
          psi.truncate();
          if(world.rank()==0) print("Iteration: ", iter++);
          real_convolution_3d Ebarop = Ebar(world,energy);
          residual = iterate(world, rightside, psi, Ebarop);
          rightside = makerightside(world, Vnuc, psi, PbarAop, Aop);
          energy = rele(world, psi, rightside, Pbarop, Vnuc);
          if(world.rank()==0) print("Energy: ", energy);
          //energy = exact;
     }

     if(world.rank()==0) {
          print("\nFinal Energy: ", energy);
          print("Calculated Relativistic Correction: ", initialguessenergy - energy);
          print("Scaled Energy (-eps/Z^2): ", -energy/(Z*Z), "\n");
     }

     world.gop.fence();
     finalize();
     return 0;
}

//kthxbye
