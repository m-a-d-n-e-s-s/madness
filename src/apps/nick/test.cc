/*********************
 * This is a 1D test of the ScatteringWF for spline checking.
 * I'm combining wavef.h wavef.cc and projPsi.cc into one file 
 * to keep test this outside of the madness framework
 * Derivations of the variable mesh are in Notebook D pg 58-60
 *********************/
#include <complex>
#include <iostream>
using std::scientific;
using std::cout;
using std::endl;
#include <fstream>
using std::ofstream;
#include <float.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <nick/hyp.h>
typedef std::complex<double> complexd;
const complexd I(0,1);
const complexd one(1,0);
const double PI = M_PI;
double Z = 1.0;
double k = 0.15;
//The furthest out we will evaluate this function is sqrt(3)*boxSize
//And the kr + kDOTr will bring along another factor of 2
//2*sqrt(3) = 3.46
double dx = 0.1;
double range  = 100.0;
const int n = floor(2*range/dx) + 1;
double* aR = new double[n];
double* bR = new double[n];
double* cR = new double[n];
double* dR = new double[n];
double* aI = new double[n];
double* bI = new double[n];
double* cI = new double[n];
double* dI = new double[n];
double* h  = new double[n-1];
double* x  = new double[n];
double toX(int i) {
    return dx*i*i/n/2;
}
int fromX( double xx ) {
    return floor( sqrt(2*n*xx/dx) );
}
complexd splined1F1(double xx) {
    // Horner's scheme is suppose to be an efficeint method of evaluating
    // polynomials. See Wikipedia for details.
    int j = fromX(xx);     // variable mesh
    //int j = floor( xx/dx ); // constant mesh
    double y = (xx - x[j]);
    return complexd(aR[j],aI[j]) + y*(complexd(bR[j],bI[j]) + 
                                      y*(complexd(cR[j],cI[j]) +
                                         y*complexd(dR[j],dI[j])
                                        )
                                     );
}
complexd gamma(double re, double im)
{
    gsl_sf_result lnr;
    gsl_sf_result arg;
    int status = gsl_sf_lngamma_complex_e(re, im, &lnr, &arg);
    if(status != 0) throw "Error: gsl_sf_lngamma: " + status;
    complexd ANS(exp(lnr.val)*cos(arg.val), exp(lnr.val)*sin(arg.val) );
    return ANS;
}
complexd expmPI_k   = exp(-PI/k);
complexd expPI_2k   = exp(PI/(2*k));
complexd gamma1pI_k = gamma(1.0,1/k);
complexd gammamI_k  = gamma(0.0,-1/k);
/****************************************************************
 * The asymptotic form of the hypergeometric function given by
 * Abramowitz and Stegun 13.5.1
 * **************************************************************/
complexd aForm3(complexd ZZ) {
    //cout << scientific;
    //cout.precision(15);
    complexd cA2 = pow(ZZ,I/k);
    complexd cB1 = exp(ZZ);
    complexd cB2 = pow(ZZ,-one-I/k);
    complexd cA = expmPI_k*cA2/gamma1pI_k;
    complexd cB = cB1*cB2/gammamI_k;
    //print("cA = ", cA,"cB = ", cB);
    complexd termA(0,0);
    complexd termB(0,0);
    int maxTerms = 20;
    complexd zrn = 1;
    complexd mzrn = 1;
    complexd zr = 1.0/ZZ;
    double nFact = 1.0;            //0! = 1
    complexd pochAA(1.0,0.0);      //Pochhammer is the counting up factorial (A)_0 = 1
    complexd poch1mAA(1.0,0.0);   //(BB-AA)_n
    for(int n=0; n<=maxTerms; n++) {
        complexd contribA = pochAA*pochAA*mzrn/nFact;
        termA += contribA;
        complexd contribB = poch1mAA*poch1mAA*zrn/nFact;
        termB += contribB;
        mzrn     *= -zr;
        zrn      *=  zr;
        nFact    *= n+1;  //(n+1) is the number to be used in the next iteration
        pochAA   *= complexd(n,-1.0/k); //(x)_n = x(x+1)(x+2)..(x+n-1)
        poch1mAA *= complexd(1+n,1.0/k);        
    }
    return cA*termA + cB*termB;
}

/*********************************************************************************
 * Here is where I splice together my two representations of the hypergeometric
 * function. See personal journal D page 65 for a derivation of the cut off
 * and the "Finding the series - aForm boundary" section of Splines.nb for the 
 * emperical results
 *********************************************************************************/
complexd f11(complexd ZZ) {
        if(fabs(imag(ZZ)) <= 40 ) return conhyp(-I/k,one,ZZ);
    else return aForm3(ZZ);
    //return  conhyp(-I/k,one,ZZ);
}
double diffR(double x) {
    return real(splined1F1(x) - f11(-I*x));
}
double diffI(double x) {
    return imag(splined1F1(x) - f11(-I*x));
}
int main(int argc, char**argv) {
    /***********************************************************************
     * By evaluating the scattering function on a 1D grid we can do a time 
     * saving table lookup when evaluating in 3D.
     * I'm following the natural cubic spline by following algorithm 3.4
     * in Numerical Analysis by Burden and Faires 7th edition to interpolate
     * between the points.
     **********************************************************************/ 
    double l[n];
    double mu[n-1];
    double zR[n];
    double zI[n];
    double alphaR[n-1];
    double alphaI[n-1];
    complexd value;
    for(int i=0; i<=n; i++) {
        //x[i] = i*dx; //constant mesh
        x[i] = toX(i); //variable mesh
        //cout << "x[" << i << "] = " << toX(i) << endl;
    }
    for(int i=0; i<=n; i++) {
        value = f11(-I*x[i]);
        aR[i] = real(value);
        aI[i] = imag(value);
        //cout << "ar[" << i << "] = " << aR[i] << endl;
        //cout << "ai[" << i << "] = " << aI[i] << endl;
    }
    mu[0] = 0.5;
    h[0]  = x[1] - x[0];
    l[0]  = 2*h[0];
    alphaR[0] = 3*(aR[1] - aR[0])/h[0] + 3.0/k;
    alphaI[0] = 3*(aI[1] - aI[0])/h[0];
    zR[0]  = alphaR[0]/l[0];
    zI[0]  = alphaI[0]/l[0];
    for(int i=1; i<=n-1; i++) {
        h[i] = x[i+1]-x[i];
        alphaR[i] = 3.0/h[i]*(aR[i+1] - aR[i]) - 3.0/h[i-1]*(aR[i] - aR[i-1]);
        alphaI[i] = 3.0/h[i]*(aI[i+1] - aI[i]) - 3.0/h[i-1]*(aI[i] - aI[i-1]);
        l[i] = 2*(x[i+1] - x[i-1] - h[i-1]*mu[i-1]);
        mu[i] = h[i]/l[i];
        zR[i] = (alphaR[i] - h[i-1]*zR[i-1])/l[i];
        zI[i] = (alphaI[i] - h[i-1]*zI[i-1])/l[i];
    }
    l[n] = 1.0;
    zR[n] = 0.0;
    zI[n] = 0.0;
    cR[n] = 0.0;
    cI[n] = 0.0;
    for(int j=n-1; j>=0; j--) {
        cR[j] = zR[j] - mu[j]*cR[j+1];
        cI[j] = zI[j] - mu[j]*cI[j+1];
        bR[j] = (aR[j+1] - aR[j])/h[j] - h[j]*(cR[j+1] + 2*cR[j])/3;
        bI[j] = (aI[j+1] - aI[j])/h[j] - h[j]*(cI[j+1] + 2*cI[j])/3;
        dR[j] = (cR[j+1] - cR[j])/(3*h[j]);
        dI[j] = (cI[j+1] - cI[j])/(3*h[j]);
    }
    ofstream F11, splined, diff, fdiff, CONHYP;
    F11.open("f11.dat");
    F11.precision(15);
    splined.open("splined.dat");
    splined.precision(15);
    diff.open("diff.dat");
    diff.precision(15);
    fdiff.open("fdiff.dat");
    fdiff.precision(15);
    CONHYP.open("conhyp.dat");
    CONHYP.precision(15);
    // DANGER DANGER
    // I'm not sure why x[0] is being set to 1, but it is.
    // Rather than find out why, I'm going to just set it back to 0.0;
    //cout << "x[0] = " << x[0] << endl; 
    //cout << "toX(0) = " << toX(0) << endl;
    x[0] = 0.0;
    for(int i=0; i<n; i++) {
        F11     << x[i] << "\t" << aR[i]  << "\t" << aI[i]  << endl;
    }
    //The spline file is sampled at approximately 7 times as many points
    //to make error calculations easy

//     for(double r=0; r < range-dx; r += dx/sqrt(50) ) { // Constant Mesh
//         splined << r << "\t" << real(splined1F1(r)) << "\t" << imag(splined1F1(r)) << endl;
//         diff    << r << "\t" << diffR(r)            << "\t" << diffI(r)            << endl;
//     }
    // Variable Mesh
    double r = 0.0;
    for(int i=0; r<range-dx; i++) {
        r = toX(i)/sqrt(50); //Now r scales the same way as the grid points do
        splined << r << "\t" << real(splined1F1(r)) << "\t" << imag(splined1F1(r)) << endl;
        diff    << r << "\t" << diffR(r)            << "\t" << diffI(r)            << endl;
        fdiff   << r << "\t" << real(conhyp(-I/k,one,-I*r) -  aForm3(-I*r) ) 
                << "\t"      << imag(conhyp(-I/k,one,-I*r) -  aForm3(-I*r) )      << endl;
    }
    delete aR;
    delete bR;
    delete cR;
    delete dR;
    delete aI;
    delete bI;
    delete cI;
    delete dI;
    delete x;
    return 0;
}
