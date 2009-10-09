//\file wavef.cc
//\brief The hydrogenic bound and continuum states
/************************************************************************
 * Here is a madness representation of the hydrogenic wave functions.
 * The bound states come from the Gnu Scientific Library. The unbound
 * states are generated with the confluent hypergeometric function which
 * uses gmp and mpfr for extended precision
 * 
 * Using: Gnu Scientific Library          http://www.gnu.org/software/gsl/
 *        GNU Multiple Precision library  http://gmplib.org/
 *        mpfr                            http://www.mpfr.org/
 * By:    Nick Vence
 ************************************************************************/

#include "wavef.h"
#include "hyp.h"
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <float.h>
#include <time.h>

using namespace madness;

time_t before, after;
//MPI printing macros
#define PRINTLINE(str) if(world.rank()==0) cout << str << endl;
//#define START_TIMER  tt=cpu_time()
//#define START_TIMER world.gop.fence(); tt=cpu_time()
#define END_TIMER_C(msg,cplx) tt=cpu_time()-tt; cout << "Timer: " << msg << "(" << real(cplx) << " + " << imag(cplx) << "I) took " << tt << " seconds" << endl
#define END_TIMER(msg) tt=cpu_time()-tt; printf("timer: %24.24s    took%8.2f seconds\n", msg, tt)
//#define END_TIMER(msg) tt=cpu_time()-tt;  if (world.rank()==0) printf("timer: %24.24s    took%8.2f seconds\n", msg, tt)
#define PRINT_COMPLEX(msg,re,im) if(world.rank()==0) printf("%34.34s %9.6f + %9.6fI\n", msg, re, im)


/*****************************************
 *Exp[ I*(k.r) ]
 *****************************************/
Expikr::Expikr( const vector3D& kVec) : kVec(kVec)
{
    double sum = 0.0;
    for(int i=0; i<NDIM; i++) { sum += kVec[i]*kVec[i]; }
    k = sqrt(sum);
    costhK = kVec[2]/k;
}
complexd Expikr::operator()(const vector3D& rVec) const
{
    double kDOTr = 0.0;
    for(int i=0; i<NDIM; i++) {
        kDOTr += kVec[i]*rVec[i];
    }
    return exp(I*kDOTr);
}

/*****************************************
 *Exp[ -I*(kr + k.r) ]
 *****************************************/
Expikr2::Expikr2( const vector3D& kVec) : kVec(kVec)
{
    double sum = 0.0;
    for(int i=0; i<NDIM; i++) { sum += kVec[i]*kVec[i]; }
    k = sqrt(sum);
    costhK = kVec[2]/k;
}
complexd Expikr2::operator()(const vector3D& rVec) const
{
    double kDOTr = 0.0;
    double r2 = 0.0;
    for(int i=0; i<NDIM; i++) { r2 += rVec[i]*rVec[i]; }
    double r = sqrt(r2);
    for(int i=0; i<NDIM; i++) {
        kDOTr += kVec[i]*rVec[i];
    }
    return exp(-I*(k*r + kDOTr));
}

/******************************************
 * BoundWF
 ******************************************/
BoundWF::BoundWF(double Z, int nn, int ll, int mm ) : Z(Z)
{
    gsl_set_error_handler_off();
    if(nn < 1) {
	cerr << "Thou shalt not have negative n!" << endl;
	exit(1);
    }
    if(ll<0 || ll>=nn) {
    cerr << "n = " << nn << "\tl = " << ll << endl;
	cerr << "l has broken the quantum commandments!" << endl;
	exit(1);
    }
    if(abs(mm) > ll) {
    cerr << "n = " << nn << "\tl = " << ll << "\tm = " << mm << endl;
	cerr << "m out of bounds error!" << endl;
	exit(1);
    }
    n=nn;
    l=ll;
    m=mm;
}
complexd BoundWF::operator()(const vector3D& rVec) const {
    double sum = 0.0;
    for(int i=0; i<NDIM; i++) { sum += rVec[i]*rVec[i]; }
    double r = sqrt(sum);
    double cosTH;
    if(r==0.0) 
        cosTH = 1.0;
    else
        cosTH = rVec[2]/r;
    gsl_sf_result Rnl;
    int status = gsl_sf_hydrogenicR_e(n, l, Z, r, &Rnl);
    //    gsl_set_error_handler(NULL);     //Turns on the default error handler
    if(status == GSL_EUNDRFLW) { return complexd(0,0); }
    else if(status != 0)            MADNESS_EXCEPTION("gsl_ERROR: ",status);
    if(m==0) { return complexd(Rnl.val * gsl_sf_legendre_sphPlm(l, m, cosTH), 0.0); }
    else {
	gsl_sf_result rPhi;
	gsl_sf_result phi ; 
	gsl_sf_rect_to_polar(rVec[0], rVec[1], &rPhi, &phi);
	return complexd(   Rnl.val 
                         * gsl_sf_legendre_sphPlm(l, abs(m), cosTH)
                         * exp(complexd(0,m*phi.val))
                       );
    }
}

struct MemberFuncPtr {
    ScatteringWF* thisObj;
    MemberFuncPtr(ScatteringWF* obj) : thisObj(obj) {}
    complexd operator()(double x) { return thisObj->f11(x); }
};

/***********************************************************************
 * The Scattering Wave Function
 * See Landau and Lifshitz Quantum Mechanics Volume 3
 * Third Edition Formula (136.9)
 **********************************************************************/
ScatteringWF::ScatteringWF(World& world, double Z, const vector3D& kVec)
    :Z(Z)
    ,kVec(kVec)
    ,k(sqrt(kVec[0]*kVec[0] + kVec[1]*kVec[1] + kVec[2]*kVec[2]))
    ,domain(k*sqrt(3)*pow(FunctionDefaults<NDIM>::get_cell_volume(),1.0/3.0))    
{
    expmPI_k   = exp(-PI/k);
    expPI_2k   = exp(PI/(2*k));
    gamma1pI_k = gamma(1.0,1/k);
    gammamI_k  = gamma(0.0,-1/k);
    expPI_2kXgamma1pI_k = expPI_2k * gamma1pI_k ;
    one = complexd(1.0, 0.0);
    TOL = 1e-12;
    /***********************************************************************
     * How far must we tabulate our 1F1 to cover the domain?
     * V^(1/3) gives us the length of the box
     * Our spherical function needs only one octant of the box    V^(1/3)/2
     * sqrt(3) allows us to reach the corner of the cube  sqrt(3)*V^(1/3)/2
     * kr + kDOTr brings along another factor of 2k     k*sqrt(3)*V^(1/3)
     **********************************************************************/
    dx = 4e-3;   //Mesh spacing
    ra = 5.0; //boundary cutoff
    n = floor(0.5*domain/dx*(1 + sqrt(1 + 4*ra/domain))) + 1;
    time( &before );
    MemberFuncPtr p1F1(this);
    fit1F1 = CubicInterpolationTable<complexd>(world, 0.0, domain, n, p1F1);
    time( &after );
    cout << "Computing the CubicInterpolationTable took " << after - before 
         << " seconds." << endl;
    //cout << "Max Error = " << fit1F1.err(p1F1) << "\t" << endl;
    //cout << "DONE with err\n";
//     boundary = ra + 0.5*dx;
//     xi = 2;   //density speedup (should be 2 or greater)
//     //cout << "\ndomain = " << domain << "\t\t ra = " << ra << "\t\tdx = " << dx << endl;
//     /*********************************************************
//      * Calculating n for a variable mesh
//      * (0)     x'[1] = dx (the coarsest mesh spacing)
//      * (1)     xI[i] = alpha*i^2          region    |   region  
//      * (2)    xII[i] = dx*i + beta           I      |     II
//      * (3)   xI'[n1] = xII'[n1]          quadratic  |   uniform
//      *    2*alpha*n1 = dx                   mesh    |    mesh
//      *         alpha = dx/2/n1          ____________|____________
//      * (4)    xI[n1] = xII[n1]          012  ...   n1    ...    n
//      *(dx/2/n1)n1*n1 = dx*n1 + beta    |<-   n1   ->|<-   n2   ->|
//      *       dx*n1/2 = dx*n1 + beta
//      *          beta = -dx*n1/2
//      *            n2 =  distanceII*densityII
//      *            n2 = (domain-ra)*1/dx
//      *            n1 =   distanceI*densityI
//      *            n1 =          ra*xi/dx
//      *             n = n1 + n2
//      ********************************************************/
//     n1 = floor(ra*xi/dx) + 1;
//     n  = floor(n1 + (domain - ra)/dx) + 1;
//     alpha = dx*dx/(2*ra*xi);
//     beta  = -dx*n1/2;
//     dx_2n1 = dx/2/n1;
//     dxn1_2 = dx*n1/2;
//     //cout << "n = " << n << "\t\t n1 = " << n1 << "\t\txi = " << xi << endl;
//     h  = new double[n];
//     x  = new double[n+1];
//     // Initializing the Variable mesh
//     for(int i=0; i<=n; i++) {
//         x[i] = i*dx; //constant mesh
//         //x[i] = toX(i); //variable mesh
//         //cout << "x[" << i << "] = " << toX(i) << endl;
//     }
//     // Spacing
//     for(int i=0; i<n; i++) {
//         h[i] = x[i+1] - x[i];
//     }
//     taylor3Fit();
    /*****************************************************
     * Below is the quality control for my splines
     ****************************************************/
//     ofstream F11, splined, diff, fdiff,Xfile;
//     F11.open("f11.dat");
//     F11.precision(15);
//     splined.open("splined.dat");
//     splined.precision(15);
//     diff.open("diff.dat");
//     diff.precision(15);
//     fdiff.open("fdiff.dat");
//     fdiff.precision(15);
//     for(int i=0; i<n; i++) {
//         //        F11     << x[i] << "\t" << aR[i]  << "\t" << aI[i]  << endl;
//         F11     << x[i] << "\t" << fR[i]  << "\t" << fI[i]  << endl;
//     }
//     //The spline file is sampled at approximately 7 times as many points
//     //to make error calculations easy
//     //     for(double r=0; r < range-dx; r += dx/sqrt(50) ) { // Constant Mesh
//     //         splined << r << "\t" << real(approx1F1(r)) << "\t" << imag(approx1F1(r)) << endl;
//     //         diff    << r << "\t" << diffR(r)            << "\t" << diffI(r)            << endl;
//     //     }
//     // Variable Mesh
//     int ii=0;
//     for(double r=0; r<domain; r +=dx/sqrt(200) ) {
//         //r scales the same way as the grid points do
//         ii++;
//         splined << r << "\t" << real(fit1F1(r)) << "\t" << imag(fit1F1(r)) << endl;
//         diff    << r << "\t" << diffR(r)            << "\t" << diffI(r)            << endl;
//         fdiff   << r << "\t" << real(conhyp(-I/k,one,-I*r) -  aForm3(-I*r) ) 
//                 << "\t"      << imag(conhyp(-I/k,one,-I*r) -  aForm3(-I*r) )      << endl;
//     }
}

ScatteringWF::ScatteringWF(double Z, const vector3D& kVec)
    :Z(Z)
    ,kVec(kVec)
    ,k(sqrt(kVec[0]*kVec[0] + kVec[1]*kVec[1] + kVec[2]*kVec[2]))
    ,domain(k*sqrt(3)*pow(FunctionDefaults<NDIM>::get_cell_volume(),1.0/3.0))    
{
    expmPI_k   = exp(-PI/k);
    expPI_2k   = exp(PI/(2*k));
    gamma1pI_k = gamma(1.0,1/k);
    gammamI_k  = gamma(0.0,-1/k);
    expPI_2kXgamma1pI_k = expPI_2k * gamma1pI_k ;
    one = complexd(1.0, 0.0);
    TOL = 1e-12;
    /***********************************************************************
     * How far must we tabulate our 1F1 to cover the domain?
     * V^(1/3) gives us the length of the box
     * Our spherical function needs only one octant of the box    V^(1/3)/2
     * sqrt(3) allows us to reach the corner of the cube  sqrt(3)*V^(1/3)/2
     * kr + kDOTr brings along another factor of 2k     k*sqrt(3)*V^(1/3)
     **********************************************************************/
    dx = 4e-3;   //Mesh spacing
    ra = 5.0; //boundary cutoff
    n = floor(0.5*domain/dx*(1 + sqrt(1 + 4*ra/domain))) + 1;
    time( &before );
    MemberFuncPtr p1F1(this);
    fit1F1 = CubicInterpolationTable<complexd>(0.0, domain, n, p1F1);
    time( &after );
    cout << "Computing the CubicInterpolationTable took " << after - before 
         << " seconds." << endl;
    //cout << "Max Error = " << fit1F1.err(p1F1) << "\t" << endl;
    //cout << "DONE with err\n";
//     boundary = ra + 0.5*dx;
//     xi = 2;   //density speedup (should be 2 or greater)
//     //cout << "\ndomain = " << domain << "\t\t ra = " << ra << "\t\tdx = " << dx << endl;
//     /*********************************************************
//      * Calculating n for a variable mesh
//      * (0)     x'[1] = dx (the coarsest mesh spacing)
//      * (1)     xI[i] = alpha*i^2          region    |   region  
//      * (2)    xII[i] = dx*i + beta           I      |     II
//      * (3)   xI'[n1] = xII'[n1]          quadratic  |   uniform
//      *    2*alpha*n1 = dx                   mesh    |    mesh
//      *         alpha = dx/2/n1          ____________|____________
//      * (4)    xI[n1] = xII[n1]          012  ...   n1    ...    n
//      *(dx/2/n1)n1*n1 = dx*n1 + beta    |<-   n1   ->|<-   n2   ->|
//      *       dx*n1/2 = dx*n1 + beta
//      *          beta = -dx*n1/2
//      *            n2 =  distanceII*densityII
//      *            n2 = (domain-ra)*1/dx
//      *            n1 =   distanceI*densityI
//      *            n1 =          ra*xi/dx
//      *             n = n1 + n2
//      ********************************************************/
//     n1 = floor(ra*xi/dx) + 1;
//     n  = floor(n1 + (domain - ra)/dx) + 1;
//     alpha = dx*dx/(2*ra*xi);
//     beta  = -dx*n1/2;
//     dx_2n1 = dx/2/n1;
//     dxn1_2 = dx*n1/2;
//     //cout << "n = " << n << "\t\t n1 = " << n1 << "\t\txi = " << xi << endl;
//     h  = new double[n];
//     x  = new double[n+1];
//     // Initializing the Variable mesh
//     for(int i=0; i<=n; i++) {
//         x[i] = i*dx; //constant mesh
//         //x[i] = toX(i); //variable mesh
//         //cout << "x[" << i << "] = " << toX(i) << endl;
//     }
//     // Spacing
//     for(int i=0; i<n; i++) {
//         h[i] = x[i+1] - x[i];
//     }
//     taylor3Fit();
    /*****************************************************
     * Below is the quality control for my splines
     ****************************************************/
//     ofstream F11, splined, diff, fdiff,Xfile;
//     F11.open("f11.dat");
//     F11.precision(15);
//     splined.open("splined.dat");
//     splined.precision(15);
//     diff.open("diff.dat");
//     diff.precision(15);
//     fdiff.open("fdiff.dat");
//     fdiff.precision(15);
//     for(int i=0; i<n; i++) {
//         //        F11     << x[i] << "\t" << aR[i]  << "\t" << aI[i]  << endl;
//         F11     << x[i] << "\t" << fR[i]  << "\t" << fI[i]  << endl;
//     }
//     //The spline file is sampled at approximately 7 times as many points
//     //to make error calculations easy
//     //     for(double r=0; r < range-dx; r += dx/sqrt(50) ) { // Constant Mesh
//     //         splined << r << "\t" << real(approx1F1(r)) << "\t" << imag(approx1F1(r)) << endl;
//     //         diff    << r << "\t" << diffR(r)            << "\t" << diffI(r)            << endl;
//     //     }
//     // Variable Mesh
//     int ii=0;
//     for(double r=0; r<domain; r +=dx/sqrt(200) ) {
//         //r scales the same way as the grid points do
//         ii++;
//         splined << r << "\t" << real(fit1F1(r)) << "\t" << imag(fit1F1(r)) << endl;
//         diff    << r << "\t" << diffR(r)            << "\t" << diffI(r)            << endl;
//         fdiff   << r << "\t" << real(conhyp(-I/k,one,-I*r) -  aForm3(-I*r) ) 
//                 << "\t"      << imag(conhyp(-I/k,one,-I*r) -  aForm3(-I*r) )      << endl;
//     }
}

ScatteringWF::~ScatteringWF() {
//     delete aR;
//     delete bR;
//     delete cR;
//     delete dR;
//     delete aI;
//     delete bI;
//     delete cI;
//     delete dI;

//     delete fR;
//     delete fI;
//     delete fpR;
//     delete fpI;
//     delete fppR;
//     delete fppI;
//     delete h;
//     delete x;
}

complexd ScatteringWF::approx1F1(double xx) const {
    int j = fromX(xx);
    double y = (xx - x[j]);
    //    cout << "y = " << y << "\t fpR[" << j << "] = " << fpR[j] << endl;
    // NEEDED FOR: splineFit()
//     return complexd(aR[j],aI[j]) + y*(complexd(bR[j],bI[j]) + 
//                                       y*(complexd(cR[j],cI[j]) +
//                                          y*complexd(dR[j],dI[j]) ) );
//      return complexd( fR[j] + y*(fpR[j] + y*(fppR[j]/2 + y*fpppR[j]/6)),
//                      fI[j] + y*(fpI[j] + y*(fppI[j]/2 + y*fpppI[j]/6)) );
    return complexd( fR[j] + y*(fpR[j] + y*fppR[j]/2 ),
                     fI[j] + y*(fpI[j] + y*fppI[j]/2) );
}

double   ScatteringWF::diffR(double x)    const {
    return real(fit1F1(x) - f11(x));
}
double   ScatteringWF::diffI(double x)    const {
    return imag(fit1F1(x) - f11(x));
}
double   ScatteringWF::toX(double s)      const {
    return dx*s*s/(ra/dx + s);
}
double   ScatteringWF::toS(double x)      const {
    return (x + sqrt(x*x + 4*ra*x))/2/dx;
}
double   ScatteringWF::toX(int i)         const {
    if( i <= n1 ) return alpha*i*i;
    return dx*i - dxn1_2;
}
int      ScatteringWF::fromX( double xx ) const {
    return floor( xx/dx );
    if( xx < boundary ) {
        return floor( sqrt(2*xi*ra*xx/dx/dx) + 0.5);
    }
    int index =  floor( (xx - beta)/dx + 0.5);
    if(index > n) {
        cout << "index = " << index << endl;
        cout << "n = "         << n << endl;
        cout << "x = "        << xx << endl;
        throw "ScatteringWF: index out of bounds\n increase domain";
    }
    return index;
}
complexd ScatteringWF::f11(double xx) const {
    complexd ZZ(0.0,-xx);
    if(xx <= 40 ) return conhyp(-I/k,one,ZZ);
    else return aForm3(ZZ);
}

complexd ScatteringWF::operator()(const vector3D& rVec) const {
    double kDOTr = kVec[0]*rVec[0] + kVec[1]*rVec[1] + kVec[2]*rVec[2];
    double r     = sqrt(rVec[0]*rVec[0] + rVec[1]*rVec[1] + rVec[2]*rVec[2]);
    //print("expPI_2kXgamma1pI_k = ", expPI_2kXgamma1pI_k, "exp(I*kDOTr) =",
    //       exp(I*kDOTr)," f11(-I*(k*r + kDOTr)) =",fit1F1(k*r + kDOTr));
    return 0.0634936359342 //  = (2PI)^-(3/2)
         * expPI_2kXgamma1pI_k
         * exp(I*kDOTr)
         * fit1F1(k*r + kDOTr);
}

/***********************************************************************
 * Method 3:
 * Lagrange interpolating polynomials are unique and they should
 * accomodate the variable mesh.
 * ********************************************************************/
void ScatteringWF::polyFit() {
    //MUST UNCOMMENT ARRAY DECLARATIONS IN wavef.h
    //MUST CHANGE RETURN IN approx1F1()
    //MUST UNCOMMENT delete statements in ~ScatteringWF()
    //choose const int order
    //make an array of gsl_interp_accel
    //make an array of gsl_interp
}

/***********************************************************************
 * Method 2:
 * Using a Third order Taylor Polynomial which doesn't exhibit propper
 * convergence behavior (error is 3rd order)
 ***********************************************************************/
void ScatteringWF::taylor3Fit() {
    //MUST UNCOMMENT ARRAY DECLARATIONS IN wavef.h
    //MUST CHANGE RETURN IN approx1F1()
    //MUST UNCOMMENT delete statements in ~ScatteringWF()
    fR   = new double[n+1];
    fpR  = new double[n+1];
    fppR = new double[n+1];
    fpppR = new double[n+1];
    fI   = new double[n+1];
    fpI  = new double[n+1];
    fppI = new double[n+1];
    fpppI = new double[n+1];
    complexd value;
    for(int i=0; i<=n; i++) {
        value = f11(x[i]);
        fR[i] = real(value);
        fI[i] = imag(value);
        //cout << "fR[" << i << "] = " << fR[i] << "\t fI[" << i << "] = " << fI[i] << endl;
    }
    // 1st derivative
    fpR[0] = -1.0/k; //CHECK
    fpI[0] = 0.0;    
    for(int i=1; i<n; i++) {
        fpR[i] = 0.5*((fR[i] - fR[i-1])/h[i-1] + (fR[i+1] - fR[i])/h[i]);
        fpI[i] = 0.5*((fI[i] - fI[i-1])/h[i-1] + (fI[i+1] - fI[i])/h[i]);
//         cout << "f'L[" << i << "] = " << (fR[i] - fR[i-1])/h[i-1] 
//              << "\t f'R[" << i << "] = " << (fR[i+1] - fR[i])/h[i]
//              << "\t fpR[" << i << "] = " << fpR[i] <<  endl;
    }
    fpR[n] = (fR[n] - fR[n-1])/h[n-1];
    fpI[n] = (fI[n] - fI[n-1])/h[n-1];
    // 2nd derivative
    fppR[0] = (fpR[1] - fpR[0])/h[0];
    fppI[0] = (fpI[1] - fpI[0])/h[0];
    for(int i=1; i<n-1; i++) {
        fppR[i] = 0.5*( (fpR[i+1] - fpR[i])/h[i] + (fpR[i] - fpR[i-1])/h[i-1]);
        fppI[i] = 0.5*( (fpI[i+1] - fpI[i])/h[i] + (fpI[i] - fpI[i-1])/h[i-1]);
        //cout << "fppR[" << i << "] = " << fppR[i] << "\t fppI[" << i << "] = " << fppI[i] << endl;
    }
    fppR[n] = (fpR[n] - fpR[n-1])/h[n-1];
    fppI[n] = (fpI[n] - fpI[n-1])/h[n-1];
    //3nd derivatives
    fpppR[0] = (fppR[1] - fppR[0])/h[0];
    fpppI[0] = (fppI[1] - fppI[0])/h[0];
    for(int i=1; i<n-1; i++) {
        fpppR[i] = 0.5*( (fppR[i+1] - fppR[i])/h[i] + (fppR[i] - fppR[i-1])/h[i-1]);
        fpppI[i] = 0.5*( (fppI[i+1] - fppI[i])/h[i] + (fppI[i] - fppI[i-1])/h[i-1]);
        //cout << "fpppR[" << i << "] = " << fpppR[i] << "\t fpppI[" << i << "] = " << fpppI[i] << endl;
    }
    fpppR[n] = (fppR[n] - fppR[n-1])/h[n-1];
    fpppI[n] = (fppI[n] - fppI[n-1])/h[n-1];
}
/**********************************************************************
 * Method 1:
 * By evaluating the scattering function on a 1D grid we can do a time 
 * saving table lookup when evaluating in 3D.
 * I'm following the natural cubic spline (algorithm 3.4)
 * and the clamped cubic spline (algorithm 3.5)
 * in Numerical Analysis by Burden and Faires 7th edition to interpolate
 * between the points.
 ***********************************************************************/
void ScatteringWF::splineFit() {
    //MUST UNCOMMENT ARRAY DECLARATIONS IN wavef.h
    //MUST CHANGE RETURN IN approx1F1()
    //MUST UNCOMMENT delete statements in ~ScatteringWF()
//     aR = new double[n+1];
//     bR = new double[n+1];
//     cR = new double[n+1];
//     dR = new double[n+1];
//     aI = new double[n+1];
//     bI = new double[n+1];
//     cI = new double[n+1];
//     dI = new double[n+1];
//     double* l = new double[n+1];
//     double* mu = new double[n+1];
//     double* zR = new double[n+1];
//     double* zI = new double[n+1];
//     double* alphaR = new double[n+1];
//     double* alphaI = new double[n+1];

//     for(int i=0; i<=n; i++) {
//         value = f11(x[i]);
//         aR[i] = real(value);
//         aI[i] = imag(value);
//         //cout << "ar[" << i << "] = " << aR[i] << "\t ai[" << i << "] = " << aI[i] << endl;
//     }
//     cout << "a[i] is initialized" << endl;
//     mu[0] = 0.5;
//     h[0]  = x[1] - x[0];
//     l[0]  = 2*h[0];
//     alphaR[0] = 3*(aR[1] - aR[0])/h[0] + 3.0/k;
//     alphaI[0] = 3*(aI[1] - aI[0])/h[0];
//     zR[0]  = alphaR[0]/l[0];
//     zI[0]  = alphaI[0]/l[0];
//     for(int i=1; i<=n-1; i++) {
//         h[i] = x[i+1]-x[i];
//         alphaR[i] = 3.0/h[i]*(aR[i+1] - aR[i]) - 3.0/h[i-1]*(aR[i] - aR[i-1]);
//         alphaI[i] = 3.0/h[i]*(aI[i+1] - aI[i]) - 3.0/h[i-1]*(aI[i] - aI[i-1]);
//         l[i] = 2*(x[i+1] - x[i-1] - h[i-1]*mu[i-1]);
//         mu[i] = h[i]/l[i];
//         zR[i] = (alphaR[i] - h[i-1]*zR[i-1])/l[i];
//         zI[i] = (alphaI[i] - h[i-1]*zI[i-1])/l[i];
//     }
//     l[n] = 1.0;
//     zR[n] = 0.0;
//     zI[n] = 0.0;
//     cR[n] = 0.0;
//     cI[n] = 0.0;
//     for(int j=n-1; j>=0; j--) {
//         cR[j] = zR[j] - mu[j]*cR[j+1];
//         cI[j] = zI[j] - mu[j]*cI[j+1];
//         bR[j] = (aR[j+1] - aR[j])/h[j] - h[j]*(cR[j+1] + 2*cR[j])/3;
//         bI[j] = (aI[j+1] - aI[j])/h[j] - h[j]*(cI[j+1] + 2*cI[j])/3;
//         dR[j] = (cR[j+1] - cR[j])/(3*h[j]);
//         dI[j] = (cI[j+1] - cI[j])/(3*h[j]);
//     }    
//     delete l;
//     delete mu;
//     delete zR;
//     delete zI;
//     delete alphaR;
//     delete alphaI;
}

/****************************************************************
 * The asymptotic form of the hypergeometric function given by
 * Abramowitz and Stegun 13.5.1
 * **************************************************************/
complexd ScatteringWF::aForm(complexd AA, complexd BB, complexd ZZ) const {
     complexd coeffA = gamma(BB)* exp(-1.0*I*PI*AA) * pow(ZZ,-AA)/gamma(BB-AA);
     complexd coeffB = gamma(BB)* exp(ZZ) * pow(ZZ,AA-BB)/gamma(AA);
     complexd termA(0,0);
     complexd termB(0,0);
     int maxTerms = 20;
     for(int n=0; n<=maxTerms; n++)
         {
             termA += pochhammer(AA,n)*pochhammer(1.0+AA-BB,n)*pow(-1.0*ZZ,-1.0*n)
                 /(double)fact(n);
         }
     for(int n=0; n<=maxTerms; n++)
         {
             termB += pochhammer(BB-AA,n)*pochhammer(1.0-AA,n)*pow(ZZ,-1.0*n)
                 /(double)fact(n);
         }
     return coeffA*termA + coeffB*termB;
}
complexd ScatteringWF::aForm1(complexd AA, complexd BB, complexd ZZ) const {
     complexd coeffA = gamma(BB)* exp(-1.0*I*PI*AA) * pow(ZZ,-AA)/gamma(BB-AA);
     complexd coeffB = gamma(BB)* exp(ZZ) * pow(ZZ,AA-BB)/gamma(AA);
     complexd termA(0,0);
     complexd termB(0,0);
     int maxTerms = 20;
     complexd zrn = 1;
     complexd mzrn = 1;
     complexd zr = 1.0/ZZ;
     
     for(int n=0; n<=maxTerms; n++) {
         termA += pochhammer(AA,n)*pochhammer(1.0+AA-BB,n)*mzrn/(double)fact(n);
         termB += pochhammer(BB-AA,n)*pochhammer(1.0-AA,n)*zrn /(double)fact(n);
         mzrn  *= -zr;
         zrn   *=  zr;
     }
     return coeffA*termA + coeffB*termB;
}
complexd ScatteringWF::aForm2(complexd AA, complexd BB, complexd ZZ) const {
    cout << scientific;
    cout.precision(15);
    //    START_TIMER;
    complexd coeffA1 = gamma(BB);
    //END_TIMER_C("gamma",BB);
    //START_TIMER;
    complexd coeffA2 = exp(-I*PI*AA);
    //END_TIMER_C("exp",-I*PI*AA);
    //START_TIMER;
    complexd coeffA3 = pow(ZZ,-AA);
    //END_TIMER_C("pow",ZZ);
    //START_TIMER;
    complexd coeffA_4= gamma(BB-AA);
    //END_TIMER_C("gamma",BB-AA);
    //START_TIMER;
    complexd coeffB1 = gamma(BB);
    //END_TIMER_C("gamma",BB);
    //START_TIMER;
    complexd coeffB2 = exp(ZZ);
    //END_TIMER_C("exp",ZZ);
    //START_TIMER;
    complexd coeffB3 = pow(ZZ,AA-BB);
    //END_TIMER_C("pow",ZZ);
    //START_TIMER;
    complexd coeffB_4= gamma(AA);
    //END_TIMER_C("gamma",AA);
    complexd coeffA = coeffA1*coeffA2*coeffA3/coeffA_4;
    complexd coeffB = coeffB1*coeffB2*coeffB3/coeffB_4;
    //Print("coeffA = ", coeffA,"coeffB = ", coeffB);
    complexd termA(0,0);
    complexd termB(0,0);
    int maxTerms = 20;
    complexd zrn = 1;
    complexd mzrn = 1;
    complexd zr = 1.0/ZZ;
    double nFact = 1.0;            //0! = 1
    complexd pochAA(1.0,0.0);      //Pochhammer is the counting up factorial (A)_0 = 1
    complexd poch1pAAmBB(1.0,0.0); //(1+AA-BB)_n
    complexd pochBBmAA(1.0,0.0);   //(BB-AA)_n
    complexd poch1mAA(1.0,0.0);   //(BB-AA)_n
    for(int n=0; n<=maxTerms; n++) {
        complexd contribA = pochAA*poch1pAAmBB*mzrn/nFact;
        termA += contribA;
        complexd contribB = pochBBmAA*poch1mAA*zrn /nFact;
        termB += contribB;
        //print("contribA = ",contribA,"\tcontribB = ",contribB, "termA =", termA, "termB =", termB);
        mzrn  *= -zr;
        zrn   *= zr;
        nFact *= n+1;  //(n+1) is the number to be used in the next iteration
        pochAA*= complexd(n,0)+AA; //(x)_n = x(x+1)(x+2)..(x+n-1)
        poch1pAAmBB*= complexd(1+n,0)+AA-BB;
        pochBBmAA  *= complexd(n,0)  +BB-AA;
        poch1mAA   *= complexd(1+n,0)   -AA;        
    }
    return coeffA*termA + coeffB*termB;
}
complexd ScatteringWF::aForm3(complexd ZZ) const {
    complexd cA2 = pow(ZZ,I/k);
    complexd cB1 = exp(ZZ);
    complexd cB2 = pow(ZZ,-one-I/k);
    complexd cA = expmPI_k*cA2/gamma1pI_k;
    complexd cB = cB1*cB2/gammamI_k;
    complexd termA(0,0);
    complexd termB(0,0);
    int maxTerms = 10;
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
        //if(abs(contribA)<TOL && abs(contribB)<TOL) break;
        //print("contribA = ",contribA,"\tcontribB = ",contribB, "termA =", termA, "termB =", termB);
        mzrn     *= -zr;
        zrn      *=  zr;
        nFact    *= n+1;  //(n+1) is the number to be used in the next iteration
        pochAA   *= complexd(n,-1.0/k); //(x)_n = x(x+1)(x+2)..(x+n-1)
        poch1mAA *= complexd(1+n,1.0/k);        
    }
    return cA*termA + cB*termB;
}
int fact(int n) {
    int result = 1;
    if( n<0 )
        {
            cout << "Factorial out of bounds error" << endl;
            cerr << "Factorial out of bounds error" << endl;
            return 0;
        }
    for(int i=n; i>1; i--) { result *= i; }
    return result;
}
void testFact(World& world) {
    if(world.rank() == 0) 
	cout << "Testing Factorial:============================================" << endl;
    if(world.rank() == 0) cout << fact(0) << endl;
    if(world.rank() == 0) cout << fact(1) << endl;
    if(world.rank() == 0) cout << fact(2) << endl;
    if(world.rank() == 0) cout << fact(3) << endl;
    if(world.rank() == 0) cout << fact(-1) << endl;
}
complexd gamma(double re, double im) {
    gsl_sf_result lnr;
    gsl_sf_result arg;
    int status = gsl_sf_lngamma_complex_e(re, im, &lnr, &arg);
    if(status != 0) throw "Error: gsl_sf_lngamma: " + status;
    complexd ANS(exp(lnr.val)*cos(arg.val), exp(lnr.val)*sin(arg.val) );
    return ANS;
}
complexd gamma(complexd AA) {
    gsl_sf_result lnr;
    gsl_sf_result arg;
    int status = gsl_sf_lngamma_complex_e(real(AA), imag(AA), &lnr, &arg);
    if(status != 0) throw "Error: gsl_sf_lngamma: " + status;
    complexd ANS(exp(lnr.val)*cos(arg.val), exp(lnr.val)*sin(arg.val) );
    return ANS;
}
void testGamma(World& world) {
    if(world.rank() == 0) cout << "Testing Gamma:================================================" << endl;
    if(world.rank() == 0) cout << "gamma(3.0,0.0) = " << gamma(3.0,0.0) << endl;
    if(world.rank() == 0) cout << "gamma(0.0,3.0) = " << gamma(0.0,3.0) << endl;
    if(world.rank() == 0) cout << "gamma(3.0,1.0) = " << gamma(3.0,1.0) << endl;
    if(world.rank() == 0) cout << "gamma(1.0,3.0) = " << gamma(1.0,3.0) << endl;
}
	
complexd pochhammer(complexd AA, int n) {
    complexd VAL(1.0,0.0);
    for(int i=0; i<n; i++)
        {
            complexd II(i,0.0);
            VAL *= AA + II;
        }
    return VAL;
}
void testPochhammer(World& world) {
    complexd ZERO(0.0,0.0);
    complexd ONE(1.0,0.0);
    complexd _ONE(-1.0,0.0);
    complexd Iplus3(3,1);
    if(world.rank() == 0) cout << "Testing Pochhammer:===========================================" << endl;
    if(world.rank() == 0) cout << "Pochhammer[" << ZERO << ", " << 0 <<"] = " << pochhammer(ZERO,0) << endl;
    if(world.rank() == 0) cout << "Pochhammer[" << ZERO << ", " << 1 <<"] = " << pochhammer(ZERO,1) << endl;
    if(world.rank() == 0) cout << "Pochhammer[" << ONE <<  ", " << 0 <<"] = " << pochhammer(ONE,0) << endl;
    if(world.rank() == 0) cout << "Pochhammer[" << ONE <<  ", " << 1 <<"] = " << pochhammer(ONE,1) << endl;
    if(world.rank() == 0) cout << "Pochhammer[" << I <<    ", " << 2 <<"] = " << pochhammer(I,2) << endl;
    if(world.rank() == 0) cout << "Pochhammer[" << I <<    ", " << 2 <<"] = " << pochhammer(I,2) << endl;
    if(world.rank() == 0) cout << "Pochhammer[" << Iplus3 <<", "<< 2 <<"] = " << pochhammer(Iplus3,2) << endl;
    if(world.rank() == 0) cout << "Pochhammer[" << _ONE << ", " << 3 <<"] = " << pochhammer(ZERO,3) << endl;
}
//The Coulomb potential
complexd V(const vector3D& r) {
    return -1.0/sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2] + 1e-8);
}

// void testWF(World& world) {
//     double dARR[3] = {0, 0, 0.5};
//     const vector3D kVec(dARR);
//     ScatteringWF phi(1.0,kVec);
//     if(world.rank() == 0)
// 	{
//             cout.precision(3);
//             testPochhammer(world);
//             testFact(world);
//             cout.precision(12);
//             if(world.rank() == 0) cout << 
//             "Testing the wave functions:===================================" << endl;
// 	}

//     //  Bound States
//     START_TIMER;
//     Function<complexd,NDIM> psi_100 = FunctionFactory<complexd,NDIM>(world).functor( 
// 	functorT(new BoundWF(1.0, 1, 0, 0)));  
//     END_TIMER("Projecting |100>            ");
//     // A second way to declare a functor
//     START_TIMER;
//     functorT psiA(new BoundWF(1.0, 2, 1, 1));  
//     Function<complexd,NDIM> psi_211 = FunctionFactory<complexd,NDIM>(world).functor( psiA );
//     END_TIMER("Projecting |211>    ");  
    
//     //  Checking orthonormality
//     START_TIMER;
//     complexd printMe = psi_100.inner(psi_100);
//     PRINT_COMPLEX("<100|100> =                        ",real(printMe),imag(printMe));
//     END_TIMER("<100|100>                            ");  
//     printMe = psi_100.inner(psi_211);
//     PRINT_COMPLEX("<100|211> =                        ",real(printMe),imag(printMe));
//     END_TIMER("<100|211>                            "); 

//     //  z-axis k vector
//     const double zHat[NDIM] = {0.0, 0.0, 1.0};
//     vector3D k1Vec(zHat);

//     //	Calculate energy of |100>
//     START_TIMER;
//     Function<complexd,NDIM> dPsi_100 = diff(psi_100,0);
//     END_TIMER("Differentiating Psi_100             ");
//     START_TIMER;
//     Function<complexd,NDIM> v = FunctionFactory<complexd,NDIM>(world).f(V);
//     END_TIMER("Projecting 1/r                      ");
//     START_TIMER;
//     complexd KE = 3*0.5*(dPsi_100.inner(dPsi_100));
//     END_TIMER("<dPsi_100|dPsi_100>                 ");
//     PRINT_COMPLEX("KE =                              ",real(KE), imag(KE));
//     START_TIMER;
//     complexd PE = psi_100.inner(v*psi_100);
//     END_TIMER("<Psi_100|V|Psi_100>                 ");
//     PRINT_COMPLEX("PE =                              ",real(PE), imag(PE));
//     PRINT_COMPLEX("The total energy is:              ",real(PE+KE), imag(PE+KE));

//     //	Calculate energy of |211>
//     START_TIMER;
//     Function<complexd,NDIM> dPsi_211 = diff(psi_211,0);
//     END_TIMER("Projecting dPsi_211                ");
//     START_TIMER;
//     KE = 3*0.5*(dPsi_211.inner(dPsi_211));
//     END_TIMER("<dPsi_211|dPsi_211>                 ");
//     PRINT_COMPLEX("KE =                              ",real(KE), imag(KE));
//     START_TIMER;
//     PE = psi_211.inner(v*psi_211);
//     END_TIMER("<Psi_211|V|Psi_211>                 ");
//     PRINT_COMPLEX("PE =                              ",real(PE), imag(PE));
//     PRINT_COMPLEX("The total energy is:              ",real(PE+KE), imag(PE+KE));

//     //  Scattering States
//     START_TIMER;
//     Function<complexd,NDIM> psi_k1 = FunctionFactory<complexd,NDIM>(world).functor( 
//                                      functorT( new ScatteringWF(1.0, k1Vec)));
//     END_TIMER("Projecting |psi_k1>                 ");

//     //  Checking orthogonality
//     START_TIMER;
//     printMe = psi_100.inner(psi_k1);
//     END_TIMER("Projecting dPsi_211                 ");
//     PRINT_COMPLEX("<Psi_k1|100>                      ",real(printMe), imag(printMe));

//     //  Linearly polarized laser operator 
//     vector3D FVec(zHat);
//     START_TIMER;
//     Function<complexd,NDIM> laserOp = FunctionFactory<complexd,NDIM>(world).functor( 
//                                                         functorT( new Expikr(FVec) ));
//     END_TIMER("Projecting ExpIkr                    ");
//     START_TIMER;
//     Function<complexd,NDIM> ket = laserOp*psi_100;
//     END_TIMER("ExpIkr * Psi_100       ");
//     START_TIMER;
//     printMe =  psi_k1.inner(ket);
//     END_TIMER("<psi_k1|Exp[ik.r]|100> =             ");
//     PRINT_COMPLEX("<psi_k1|Exp[ik.r]|100> =           ",real(printMe), imag(printMe));

//     //  Off-axis k vector
//     double k   = 1.0; 
//     double thK = PI/4;
//     const double temp2[NDIM] = {0.0, k*sin(thK), k*cos(thK)};
//     vector3D k1_45Vec(temp2);

//     //  Off-Axis Positive Energy State
//     START_TIMER;
//     Function<complexd,NDIM> psi_k1_45 = FunctionFactory<complexd,NDIM>(world).functor( 
//                                         functorT( new ScatteringWF(1.0, k1_45Vec) ));
//     END_TIMER("Projecting |psi_k1_45>               ");
//     START_TIMER;
//     printMe =  psi_k1_45.inner(laserOp*psi_100);
//     END_TIMER("<psi_k1_45|Exp[ik.r]|100> =             ");
//     PRINT_COMPLEX("<psi_k1_45|Exp[ik.r]|100> =           ",real(printMe),imag(printMe)); 
// }

