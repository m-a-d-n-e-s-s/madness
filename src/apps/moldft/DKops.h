#ifndef MADNESS_APPS_MOLDFT_DKOPS_H_INCLUDED
#define MADNESS_APPS_MOLDFT_DKOPS_H_INCLUDED


#include <madness/mra/mra.h>
#include <madness/constants.h>
#include <madness/mra/operator.h>
#include <iostream>
#include <fstream>
#include <exception>

//#include "../../../../mpfrc++-3/mpreal.h"

using namespace madness;

static const double opthresh = 1e-16;

double q(double t){
     return exp(-t);
}

double w(double t, double eps){
     double c = 137.0359895;
     double PI = constants::pi;
     return 1.0/(c*sqrt(PI))*exp(-1.0/2.0*t-c*c*exp(-t))-(1.0+eps/(c*c))*exp((2.0*eps+(eps*eps)/(c*c))*exp(-t)-t)*(erfc((c+eps/c)*(exp(-t/2.0)))-2.0);
}

double Ebark(double t, double eps){
     double R = 1e-8;
     return pow(1.0/(2.0*q(t)), 1.5)*w(t, eps)*exp(-1.0/(4.0*q(t))*R*R);
}

real_convolution_3d Ebar(World& world, double eps){
     Vector<double,3> args = vec(0.0,0.0,0.0);


     std::vector<double> cvec;
     std::vector<double> tvec;

     double fac = pow(2.0*constants::pi,-1.5);
     double dt = 1.0/8.0;
     double tval = -10.0;
     while(tval < 100.0){

          if(Ebark(tval, eps) > opthresh){
               //if(world.rank()==0) print("   Ebark(",tval,",",eps," = ",Ebark(tval,eps));
               cvec.push_back(dt/pow(2.0*q(tval),1.5)*w(tval,eps)*fac);
               tvec.push_back(1.0/(4.0*q(tval)));
          }
          //else{
               //if(world.rank()==0) print("Ebark(",tval,",",eps," = ",Ebark(tval,eps));

          //}
          tval = tval + dt;
     }
     
     Tensor<double> ctensor(cvec.size());
     Tensor<double> ttensor(tvec.size());
     for(int i = 0; i < cvec.size(); i++){
          ctensor[i] = cvec[i];
          ttensor[i] = tvec[i];
     }

     int n = cvec.size();
     //if(world.rank()==0) print("Made an Ebar! n = ",n);
     //if(world.rank()==0) print("Made an Ebar! Here's what's inside:\n\nc:\n",ctensor,"\nt:\n",ttensor);

     //test ebar?
     /*
     if(world.rank()==0) print("Testing Ebar. Value should be ~22.680321967196");
     double testvalue = 0.0;
     for(int i = 0; i < cvec.size(); i++){
          testvalue += ctensor[i]*exp(-1e-2*ttensor[i]);
     }
     if(world.rank()==0) print("And the result is: ", testvalue);
     */

     //return real_convolution_3d(world, args, ctensor, ttensor);
     return real_convolution_3d(world, ctensor, ttensor);
}

//real_convolution_3d Ebar_fixed(World& world){
//     //Tensor<double> c(305l), t(305l);
//     std::vector<double> c(365), t(365);
//     Vector<double,3> args = vec(0.0,0.0,0.0);
//     #include "RelCoeffs/Ebar_coeffs.dat"
//     int n = 305;
//     
//     for(int i=0; i < n; i++){
//          if(c[i]*exp(-t[i]*1e-16) < opthresh){
//               c.erase(c.begin()+i);
//               t.erase(t.begin()+i);
//               i--;
//               n--;
//          }
//     }
//     
//     //print("n = ", n);
//
//     Tensor<double> ctens(n), ttens(n);
//     for(int i = 0; i < n; i++){
//          ctens[i] = c[i];
//          ttens[i] = t[i];
//     }
//     //if(world.rank()==0) print("Made a Ebar (fixed)! n = ", n);
//     //if(world.rank()==0) print("Made a Pbar! Here's what's inside:\n\nc:\n",ctens,"\nt:\n",ttens);
//     //return real_convolution_3d(world, args, ctens, ttens);
//     return real_convolution_3d(world, ctens, ttens);
//}

real_convolution_3d Pbar(World& world){
     //Tensor<double> c(305l), t(305l);
     std::vector<double> c;
     std::vector<double> t;
     Vector<double,3> args = vec(0.0,0.0,0.0);

     std::ifstream inf("/gpfs/home/jscanderson/DKproject/Pbar_t.csv");
     if(!inf){
          if(world.rank() == 0) std::cerr << "Unable to open Pbar_t.csv" << std::endl;
          exit(1);
     }
     std::string strInput;
     getline(inf, strInput);
     while(inf){
          t.push_back(std::stod(strInput));
          getline(inf, strInput);
     }
     inf.close();
     //if(world.rank()==0) print("Done Reading Pbar_t");

     inf.open("/gpfs/home/jscanderson/DKproject/Pbar_c.csv");
     if(!inf){
          if(world.rank() == 0) std::cerr << "Unable to open Pbar_c.csv" << std::endl;
          exit(1);
     }
     getline(inf, strInput);
     while(inf){
          //if(world.rank()==0) print("a", strInput, "b");
          c.push_back(std::stod(strInput));
          getline(inf, strInput);
     }
     inf.close();
     //if(world.rank()==0) print("Done Reading Pbar_c");
     

     int n = c.size();
     /*
     for(int i=0; i < n; i++){
          if(c[i]*exp(-t[i]*1e-30) < opthresh){
               c.erase(c.begin()+i);
               t.erase(t.begin()+i);
               i--;
               n--;
          }
     }
     */
     //print("n = ", n);

     Tensor<double> ctens(n), ttens(n);
     for(int i = 0; i < n; i++){
          ctens[i] = c[i];
          ttens[i] = t[i];
     }
     //if(world.rank()==0) print("Made a Pbar! n = ", n);
     //if(world.rank()==0) print("Made a Pbar! Here's what's inside:\n\nc:\n",ctens,"\nt:\n",ttens);
     //return real_convolution_3d(world, args, ctens, ttens);
     return real_convolution_3d(world, ctens, ttens);
}

real_convolution_3d Tbar(World& world){
     //Tensor<double> c(305l), t(305l);
     std::vector<double> c;
     std::vector<double> t;
     Vector<double,3> args = vec(0.0,0.0,0.0);

     std::ifstream inf("/gpfs/home/jscanderson/DKproject/Tbar_t.csv");
     if(!inf){
          if(world.rank() == 0) std::cerr << "Unable to open Tbar_t.csv" << std::endl;
          exit(1);
     }
     std::string strInput;
     getline(inf, strInput);
     while(inf){
          t.push_back(std::stod(strInput));
          getline(inf, strInput);
     }
     inf.close();
     //if(world.rank()==0) print("Done Reading Tbar_t");

     inf.open("/gpfs/home/jscanderson/DKproject/Tbar_c.csv");
     if(!inf){
          if(world.rank() == 0) std::cerr << "Unable to open Tbar_c.csv" << std::endl;
          exit(1);
     }
     getline(inf, strInput);
     while(inf){
          //if(world.rank()==0) print("a", strInput, "b");
          c.push_back(std::stod(strInput));
          getline(inf, strInput);
     }
     inf.close();
     //if(world.rank()==0) print("Done Reading Tbar_c");
     

     int n = c.size();
     /*
     for(int i=0; i < n; i++){
          if(c[i]*exp(-t[i]*1e-30) < opthresh){
               c.erase(c.begin()+i);
               t.erase(t.begin()+i);
               i--;
               n--;
          }
     }
     */
     //print("n = ", n);

     Tensor<double> ctens(n), ttens(n);
     for(int i = 0; i < n; i++){
          ctens[i] = c[i];
          ttens[i] = t[i];
     }
     //if(world.rank()==0) print("Made a Tbar! n = ", n);
     //if(world.rank()==0) print("Made a Pbar! Here's what's inside:\n\nc:\n",ctens,"\nt:\n",ttens);
     //return real_convolution_3d(world, args, ctens, ttens);
     return real_convolution_3d(world, ctens, ttens);
}


real_convolution_3d A(World& world){
     //Tensor<double> c(301l), t(301l);
     std::vector<double> c, t;
     Vector<double,3> args = vec(0.0,0.0,0.0);
     
     std::ifstream inf("/gpfs/home/jscanderson/DKproject/A_t.csv");
     if(!inf){
          if(world.rank() == 0) std::cerr << "Unable to open A_t.csv" << std::endl;
          exit(1);
     }
     std::string strInput;
     getline(inf, strInput);
     while(inf){
          t.push_back(std::stod(strInput));
          getline(inf, strInput);
     }
     inf.close();
     //if(world.rank()==0) print("Done Reading A_t");

     inf.open("/gpfs/home/jscanderson/DKproject/A_c.csv");
     if(!inf){
          if(world.rank() == 0) std::cerr << "Unable to open A_c.csv" << std::endl;
          exit(1);
     }
     getline(inf, strInput);
     while(inf){
          c.push_back(std::stod(strInput));
          getline(inf, strInput);
     }
     inf.close();
     //if(world.rank()==0) print("Done Reading A_c");

     int n = c.size();
     /*
     for(int i=0; i < n; i++){
          if(c[i]*exp(-t[i]*1e-30) < opthresh){
               c.erase(c.begin()+i);
               t.erase(t.begin()+i);
               i--;
               n--;
          }
     }
     */
     Tensor<double> ctens(n), ttens(n);
     for(int i = 0; i < n; i++){
          ctens[i] = c[i];
          ttens[i] = t[i];
     }
     //if(world.rank()==0) print("Made an A!, n = ", n);
     //if(world.rank()==0) print("Made an A! Here's what's inside:\n\nc:\n",ctens,"\nt:\n",ttens);
     //return real_convolution_3d(world, args, ctens, ttens);
     return real_convolution_3d(world, ctens, ttens);
}

//Same as A(World& world) above but adds 1/sqrt(2) to the largest coefficient
real_convolution_3d A2(World& world){
     //Tensor<double> c(301l), t(301l);
     std::vector<double> c, t;
     Vector<double,3> args = vec(0.0,0.0,0.0);
     
     std::ifstream inf("/gpfs/home/jscanderson/DKproject/A_t.csv");
     if(!inf){
          if(world.rank() == 0) std::cerr << "Unable to open A_t.csv" << std::endl;
          exit(1);
     }
     std::string strInput;
     getline(inf, strInput);
     while(inf){
          t.push_back(std::stod(strInput));
          getline(inf, strInput);
     }
     inf.close();
     //if(world.rank()==0) print("Done Reading A_t");

     inf.open("/gpfs/home/jscanderson/DKproject/A_c.csv");
     if(!inf){
          if(world.rank() == 0) std::cerr << "Unable to open A_c.csv" << std::endl;
          exit(1);
     }
     getline(inf, strInput);
     while(inf){
          c.push_back(std::stod(strInput));
          getline(inf, strInput);
     }
     inf.close();
     //if(world.rank()==0) print("Done Reading A_c");

     int n = c.size();
     /*
     for(int i=0; i < n; i++){
          if(c[i]*exp(-t[i]*1e-30) < opthresh){
               c.erase(c.begin()+i);
               t.erase(t.begin()+i);
               i--;
               n--;
          }
     }
     */
     double fac = pow(2.0*constants::pi,1.5);
     Tensor<double> ctens(n), ttens(n);
     for(int i = 0; i < n; i++){
          ctens[i] = c[i]*fac;
          ttens[i] = t[i];
     }
     //c[n-1] += 1/sqrt(2)*fac*std::pow(t[n-1]*constants::pi,1.5); //largest exponent is at the end, so add 1/sqrt(2) to the coefficient there
     //if(world.rank()==0) print("Made an A!, n = ", n);
     //if(world.rank()==0) print("Made an A! Here's what's inside:\n\nc:\n",ctens,"\nt:\n",ttens);
     //return real_convolution_3d(world, args, ctens, ttens);
     return real_convolution_3d(world, ctens, ttens);
}

real_convolution_3d PbarA(World& world){
     //Tensor<double> c(441l), t(441l);
     std::vector<double> c, t;
     Vector<double,3> args = vec(0.0,0.0,0.0);
     
     std::ifstream inf("/gpfs/home/jscanderson/DKproject/PbarA_t.csv");
     if(!inf){
          if(world.rank() == 0) std::cerr << "Unable to open PbarA_t.csv" << std::endl;
          exit(1);
     }
     std::string strInput;
     getline(inf, strInput);
     while(inf){
          t.push_back(std::stod(strInput));
          getline(inf, strInput);
     }
     inf.close();
     //if(world.rank()==0) print("Done Reading PbarA_t");

     inf.open("/gpfs/home/jscanderson/DKproject/PbarA_c.csv");
     if(!inf){
          if(world.rank() == 0) std::cerr << "Unable to open PbarA_c.csv" << std::endl;
          exit(1);
     }
     getline(inf, strInput);
     while(inf){
          c.push_back(std::stod(strInput));
          getline(inf, strInput);
     }
     inf.close();
     //if(world.rank()==0) print("Done Reading PbarA_c");

     int n = c.size();
     /*
     for(int i=0; i < n; i++){
          if(c[i]*exp(-t[i]*1e-30) < opthresh){
               c.erase(c.begin()+i);
               t.erase(t.begin()+i);
               i--;
               n--;
          }
     }
     */
     double fac = pow(2.0*constants::pi,1.5);
     Tensor<double> ctens(n), ttens(n);
     for(int i = 0; i < n; i++){
          ctens[i] = c[i]*fac;
          ttens[i] = t[i];
     }
     //if(world.rank()==0) print("Made a PbarA!, n = ", n);
     //if(world.rank()==0) print("Made a PbarA! Here's what's inside:\n\nc:\n",ctens,"\nt:\n",ttens);
     //return real_convolution_3d(world, args, ctens, ttens);
     return real_convolution_3d(world, ctens, ttens);
}

std::vector<real_convolution_3d_ptr> gradPbarA(World& world){
     const double pi = constants::pi;
     std::vector<double> c, t;
     
     std::ifstream inf("/gpfs/home/jscanderson/DKproject/PbarA_t.csv");
     if(!inf){
          if(world.rank() == 0) std::cerr << "Unable to open PbarA_t.csv" << std::endl;
          exit(1);
     }
     std::string strInput;
     getline(inf, strInput);
     while(inf){
          t.push_back(std::stod(strInput));
          getline(inf, strInput);
     }
     inf.close();
     //if(world.rank()==0) print("Done Reading PbarA_t");

     inf.open("/gpfs/home/jscanderson/DKproject/PbarA_c.csv");
     if(!inf){
          if(world.rank() == 0) std::cerr << "Unable to open PbarA_c.csv" << std::endl;
          exit(1);
     }
     getline(inf, strInput);
     while(inf){
          c.push_back(std::stod(strInput));
          getline(inf, strInput);
     }
     inf.close();
     //if(world.rank()==0) print("Done Reading PbarA_c");

     int n = c.size();
     /*
     for(int i=0; i < n; i++){
          if(c[i]*exp(-t[i]*1e-30) < opthresh){
               c.erase(c.begin()+i);
               t.erase(t.begin()+i);
               i--;
               n--;
          }
     }
     */
     double fac = pow(2.0*pi,1.5);
     for(int i = 0; i < n; i++){
          c[i] = c[i]*fac;
     }

     const Tensor<double> width = FunctionDefaults<3>::get_cell_width();
     const int k = FunctionDefaults<3>::get_k();
     double hi = width.normf(); // Diagonal width of cell

     std::vector<real_convolution_3d_ptr> gradG(3);

     for (int dir=0; dir<3; dir++) {
          std::vector< ConvolutionND<double,3> > ops(n);
          for (int mu=0; mu<n; mu++) {
          // We cache the normalized operator so the factor is the value we must multiply
          // by to recover the coeff we want.
          double cc = std::pow(t[mu]/pi,1.5); // Normalization coeff
          ops[mu].setfac(c[mu]/cc/width[dir]);

               for (int d=0; d<3; d++) {        
                    if (d != dir) ops[mu].setop(d,GaussianConvolution1DCache<double>::get(k, t[mu]*width[d]*width[d], 0, false));
               }
          ops[mu].setop(dir,GaussianConvolution1DCache<double>::get(k, t[mu]*width[dir]*width[dir], 1, false));
          }
          gradG[dir] = real_convolution_3d_ptr(new SeparatedConvolution<double,3>(world, ops));
     }

     return gradG;

}

#endif
