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

//MPI printing macros
time_t before, after;
#define PRINTLINE(str) if(world.rank()==0) cout << str << endl;
#define END_TIMER_C(msg,cplx) tt=cpu_time()-tt; cout << "Timer: " << msg << "(" << real(cplx) << " + " << imag(cplx) << "I) took " << tt << " seconds" << endl
#define END_TIMER(msg) tt=cpu_time()-tt; printf("timer: %24.24s    took%8.2f seconds\n", msg, tt)

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

/**********************************************************************
 * The Scattering Wave Function
 * See Landau and Lifshitz Quantum Mechanics Volume 3
 * Third Edition Formula (136.9)
 ********************************************************************** 
 * How far must we tabulate our 1F1 to cover the domain?
 * V^(1/3) gives us the length of the box
 * Our spherical function needs only one octant of the box    V^(1/3)/2
 * sqrt(3) allows us to reach the corner of the cube  sqrt(3)*V^(1/3)/2
 * kr + kDOTr brings along another factor of 2k     k*sqrt(3)*V^(1/3)
 **********************************************************************/
///World is needed for timing the length of the CubicInterpolationTable
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
    dx = 4e-3;   //Mesh spacing <- OPTIMIZE
    ra = 5.0; //boundary cutoff <- Variable Mesh
    n = floor(0.5*domain/dx*(1 + sqrt(1 + 4*ra/domain))) + 1;
//     time( &before );
    MemberFuncPtr p1F1(this);
    fit1F1 = CubicInterpolationTable<complexd>(world, 0.0, domain, n, p1F1);
//     time( &after );
//     PRINTLINE("Computing the CubicInterpolationTable took " << after - before 
//          << " seconds.");
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
    dx = 4e-3;   //Mesh spacing
    ra = 5.0; //boundary cutoff
    n = floor(0.5*domain/dx*(1 + sqrt(1 + 4*ra/domain))) + 1;
    //time( &before );
    MemberFuncPtr p1F1(this);
    fit1F1 = CubicInterpolationTable<complexd>(0.0, domain, n, p1F1);
    //time( &after );
    //PRINT("Computing the CubicInterpolationTable took " << after - before);
    //PRINTLINE( " seconds.");
}

complexd ScatteringWF::approx1F1(double xx) const {
    int j = fromX(xx);
    double y = (xx - x[j]);
    return complexd( fR[j] + y*(fpR[j] + y*fppR[j]/2 ),
                     fI[j] + y*(fpI[j] + y*fppI[j]/2) );
}

double   ScatteringWF::diffR(double x)    const {
    complexd ZZ(0.0,-x);
    return real(aForm3(ZZ) - conhyp(-I/k,one,ZZ));
}
double   ScatteringWF::diffI(double x)    const {
    complexd ZZ(0.0,-x);
    return imag(aForm3(ZZ) - conhyp(-I/k,one,ZZ));
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
    //The cutoff was done by finding the minimum difference between
    //conhyp(k,r) - aForm(k,r) for different k values
    //20 + 7exp(-6k) is the emperical fit
    //if(xx <= 20 + 7*std::exp(-6*k)) return conhyp(-I/k,one,ZZ);
    if(xx <= 4.5/k/k + 1) return conhyp(-I/k,one,ZZ);
    else return aForm3(ZZ);
}

/****************************************************
 * The Scattering Wave Function
 * See Landau and Lifshitz Quantum Mechanics Volume 3
 * Third Edition Formula (136.9)
 ****************************************************/
complexd ScatteringWF::operator()(const vector3D& rVec) const {
    double kDOTr = kVec[0]*rVec[0] + kVec[1]*rVec[1] + kVec[2]*rVec[2];
    double r     = sqrt(rVec[0]*rVec[0] + rVec[1]*rVec[1] + rVec[2]*rVec[2]);
    return 0.0634936359342 //  = (2PI)^-(3/2)
         * expPI_2kXgamma1pI_k
         * exp(I*kDOTr)
         * fit1F1(k*r + kDOTr);
}

/****************************************************************
 * The asymptotic form of the hypergeometric function given by
 * Abramowitz and Stegun 13.5.1
 * **************************************************************/
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
        //print("contribA = ",contribA,"\tcontribB = ",contribB, "termA =", termA, "termB =", termB);
        mzrn     *= -zr;
        zrn      *=  zr;
        nFact    *= n+1;  //(n+1) is the number to be used in the next iteration
        pochAA   *= complexd(n,-1.0/k); //(x)_n = x(x+1)(x+2)..(x+n-1)
        poch1mAA *= complexd(1+n,1.0/k);        
    }
    return cA*termA + cB*termB;
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
