/**********************
 * By:Nick Vence
 * Testing the Hypergeometric 1F1 
 * ********************/

#include <nick/wavef.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <ctime>
#include <string>

clock_t breakTime(string message, clock_t started);


/******************************************
 * ScatteringWF
 ******************************************/
ScatteringWF::ScatteringWF(double M, double Z, const coordT& kVec) : 
    WaveFunction(M,Z), kVec(kVec) 
{
    double sum = 0.0;
    for(int i=0; i<NDIM; i++) { sum += kVec[i]*kVec[i]; }
    k = sqrt(sum);
    //z = r*cos(th)
    //cos(th) = z/r
    costhK = kVec[2]/k;
}

complexd ScatteringWF::operator()(const coordT& rVec) const
{
    double sum = 0.0;
    for(int i=0; i<NDIM; i++) { sum += rVec[i]*rVec[i]; }
    double r = sqrt(sum);
    double kDOTr = 0.0;
    for(int i=0; i<NDIM; i++) 
    { 
	kDOTr += rVec[i]*kVec[i]; 
    }
    return exp(PI/(2*k))
	 * gamma(1.0+I/k)
 	 * exp(I*kDOTr)
	 * f11(-1.0*I/k, 1.0, -1.0*I*(k*r + kDOTr) );
}


/******************************************
 * BoundWF
 ******************************************/
BoundWF::BoundWF(double M, double Z, int nn, int ll, int mm ) : WaveFunction(M,Z)
{
    if(nn < 1) {
	cerr << "Thou shalt not have negative n!" << endl;
	exit(1);
    }
    if(ll<0 || ll>=nn) {
	cerr << "l has broken the quantum commandments!" << endl;
	exit(1);
    }
    if(abs(mm) > ll) {
	cerr << "m out of bounds error!" << endl;
	exit(1);
    }
    n=nn;
    l=ll;
    m=mm;
}

complexd BoundWF::operator()(const coordT& rVec) const
{
    double sum = 0.0;
    for(int i=0; i<NDIM; i++) { sum += rVec[i]*rVec[i]; }
    double r = sqrt(sum);
    //z = r*cos(th)
    //cos(th) = z/r
    double costh = rVec[2]/r;
    if(m==0) {
	return complexd(gsl_sf_hydrogenicR(n, l, 1.0*_Z, r) 
		      * gsl_sf_legendre_sphPlm(l, m, costh), 0.0);
    }
    else {
	gsl_sf_result rPhi;
	gsl_sf_result phi ; 
	gsl_sf_rect_to_polar(rVec[0], rVec[1], &rPhi, &phi);
	return gsl_sf_hydrogenicR(n, l, _Z, r)
	    * gsl_sf_legendre_sphPlm(l, m, costh)
	    * exp(complexd(0.0, m*phi.val));
    }
}

/*************************
 *WaveFuction Constructor
 *************************/
WaveFunction::WaveFunction(double M, double Z) 
{
    _M = abs(M);
    _Z = Z;
}


/*****************************************
 *Exp[ I*(k.r) ]
 *****************************************/
Expikr::Expikr( const coordT& kVec)
{
    double sum = 0.0;
    for(int i=0; i<NDIM; i++) { sum += kVec[i]*kVec[i]; }
    k = sqrt(sum);
    //z = r*cos(th)
    //cos(th) = z/r
    costhK = kVec[2]/k;
}

complexd Expikr::operator()(const coordT& rVec) const
{
    double kDOTr = 0.0;
    for(int i=0; i<NDIM; i++) {
	kDOTr += kVec[i]*rVec[i];
    }
    return exp(I*kDOTr);
}

    
/*******************************************************
 * Here is where I splice together my two functions
 * See personal journal C page 31 for a derivation of the cut off
 *******************************************************/
complexd f11(complexd AA, complexd BB, complexd ZZ)
{
    double k = 1.0/abs(imag(AA));
    if(abs(imag(ZZ)) <= 11.0 + 1.0/k + 0.5/(k*k) ) return conHyp(AA,BB,ZZ);
    else return aForm(AA,BB,ZZ);
}

/*********************************************************
 *The function from Barnett's code
 *Currently I'm not using it, but it may be faster in some regions of 1F1
 *I'll have to check it out later
 *********************************************************/
complexd hypergf(complexd AA, complexd BB, complexd XX)
{
    double EPS = 5E-15;
    int LIMIT  = 20000;
    int KIND   = 1;      	//Using normal precision in complex artithmetic 
    double ERR = 1E-18;
    int NITS   = 0;
    double FPMAX = 1E200;
    double ACC8  = 2E-18;
    double ACC16 = 2E-200;
    return   hypergf_(&AA, &BB, &XX, &EPS, &LIMIT, &KIND,
		      &ERR, &NITS, &FPMAX, &ACC8, &ACC16);
}

/*********************************************************
 *Tom707 Confluent Hypergeometric function
 ********************************************************/
complexd conHyp(complexd AA, complexd BB, complexd XX)
{
    int IP = 0;
    int LNCHF = 0;
    return conhyp_(&AA, &BB, &XX, &LNCHF, &IP);
}


/****************************************************************
 * The asymptotic form of the hypergeometric function given by
 * Abramowitz and Stegun 13.5.1
 * **************************************************************/
complexd aForm(complexd AA, complexd BB, complexd ZZ)
{
    complexd coeffA = gamma(BB)* exp(-1.0*I*PI*AA) * pow(ZZ,-AA)/gamma(BB-AA);
    complexd coeffB = gamma(BB)* exp(ZZ) * pow(ZZ,AA-BB)/gamma(AA);
    complexd termA(0,0);
    complexd termB(0,0);
    int maxTerms = 10;
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

void test1F1(complexd (*func1F1)(complexd,complexd,complexd), char *fName)
{
    cout << "Testing 1F1:==================================================" << endl;
    double r[] = { 1.0, 10.0, 100.0, 1000.0, 10000.0 };
    double k[] = { 1.0, 0.1, 0.01, 0.001};
    double theta = 0.0;
    double thetaK = PI-0.1;
    int rMAX   = sizeof(r)/sizeof(r[0]);
    int kMAX   = sizeof(k)/sizeof(k[0]);
    complexd BB(1.0,0.0);
    for(int i=0; i<kMAX; i++)
    {
	complexd AA(0.0, -1.0/k[i]);
	cout << "***********" << endl;
	for(int j=0; j<rMAX; j++)
	{
	    complexd ZZ(0.0, -1*k[i]*r[j]*(1 + cos(theta)*cos(thetaK) ));
	    printf( "%s[ %5.0fI, 1, %9.3f]= %+.10e, %+.10eI\n",fName, 
		    imag(AA), imag(ZZ), real((*func1F1)(AA, BB, ZZ)),
		    imag((*func1F1)(AA, BB, ZZ))
		);
	}
    }
}

int fact(int n)
{
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

void testFact()
{
    cout << "Testing Factorial:============================================" << endl;
    cout << fact(0) << endl;
    cout << fact(1) << endl;
    cout << fact(2) << endl;
    cout << fact(3) << endl;
    cout << fact(-1) << endl;
}

complexd gamma(complexd AA)
{
    gsl_sf_result lnr;
    gsl_sf_result arg;
    gsl_sf_lngamma_complex_e(real(AA), imag(AA), &lnr, &arg);
    complexd ANS(exp(lnr.val)*cos(arg.val), exp(lnr.val)*sin(arg.val) );
    return ANS;
}

void testGamma(double re, double im)
{
    cout << "Testing Gamma:================================================" << endl;
    gsl_sf_result lnr;
    gsl_sf_result arg;
    gsl_sf_lngamma_complex_e(re, im, &lnr, &arg);
    cout << "lnr = " <<  lnr.val << endl;
    cout << "arg = " <<  arg.val << endl;
    cout << "re = "  <<  exp(lnr.val)*cos(arg.val) << endl;
    cout << "im = "  <<  exp(lnr.val)*sin(arg.val) << endl;
}
	
complexd pochhammer(complexd AA, int n)
{
    complexd VAL(1.0,0.0);
    for(int i=0; i<n; i++)
    {
	complexd II(i,0.0);
	VAL *= AA + II;
    }
    return VAL;
}

void testPochhammer()
{
    complexd ZERO(0.0,0.0);
    complexd ONE(1.0,0.0);
    complexd _ONE(-1.0,0.0);
    complexd Iplus3(3,1);
    cout << "Testing Pochhammer:===========================================" << endl;
    cout << "Pochhammer[" << ZERO << ", " << 0 <<"] = " << pochhammer(ZERO,0) << endl;
    cout << "Pochhammer[" << ZERO << 1 <<"] = " << pochhammer(ZERO,1) << endl;
    cout << "Pochhammer[" << ONE <<  0 <<"] = " << pochhammer(ONE,0) << endl;
    cout << "Pochhammer[" << ONE <<  1 <<"] = " << pochhammer(ONE,1) << endl;
    cout << "Pochhammer[" << I <<    2 <<"] = " << pochhammer(I,2) << endl;
    cout << "Pochhammer[" << I <<    2 <<"] = " << pochhammer(I,2) << endl;
    cout << "Pochhammer[" << Iplus3 <<  2 <<"] = " << pochhammer(Iplus3,2) << endl;
    cout << "Pochhammer[" << _ONE << 3 <<"] = " << pochhammer(ZERO,3) << endl;
}

clock_t breakTime(string message, clock_t start) 
{
    cout.precision(1);
    clock_t finish = clock();
    cout << message << " took " << ((double)finish - start)/CLOCKS_PER_SEC << " seconds " << endl;
    cout.precision(12);
    return clock();
}

void f11Tester(World& world) 
{
    cout.precision(3);
    testPochhammer();
    testGamma(1.0,0.0);
    testFact();
    cout.precision(12);
    test1F1(f11   ,"   f11");
    clock_t time = clock();
    cout << "Testing the wave functions:==================================" << endl;

//  Bound States
    Function<complexd,NDIM> psi_100 = FunctionFactory<complexd,NDIM>(world).functor( 
	functorT(new BoundWF(1.0, 1.0, 1, 0, 0)));  
    time = breakTime("Projecting |100>            ",time);
    // A second way to declare a functor
    functorT psiA(new BoundWF(1.0, 1.0, 2, 1, 1));  
    Function<complexd,NDIM> psi_211 = FunctionFactory<complexd,NDIM>(world).functor( psiA );
    time = breakTime("Projecting |211>            ",time);  

//  Checking orthonormality
    cout << "<100|100> = " << psi_100.inner(psi_100) << endl;
    time = breakTime("<100|100>                   ",time);  
    cout << "<100|211> = " << psi_100.inner(psi_211) << endl;
    time = breakTime("<100|211>                   ",time);  

//  z-axis k vector
    const double zHat[NDIM] = {0.0, 0.0, 1.0};
    coordT k1Vec(zHat);

//  Scattering States
    Function<complexd,NDIM> psi_k1 = FunctionFactory<complexd,NDIM>(world).functor( 
	functorT( new ScatteringWF(1.0, 1.0, k1Vec)));
    time = breakTime("Projecting |k=1.0z^Hat>     ",time);

//  Checking orthogonality
    cout << "<100|psi_k1> = " << psi_100.inner(psi_k1) << endl;
    time = breakTime("<100|k=1.0z^Hat>            ",time);  

//  Linearly polarized laser operator 
    coordT FVec(zHat);
    Function<complexd,NDIM> laserOp = FunctionFactory<complexd,NDIM>(world).functor( 
	functorT( new Expikr(FVec) ));
    time = breakTime("Projecting laserOp          ",time);  
    cout << "<psi_k1|Exp[ik.r]|100> = " << psi_k1.inner(laserOp*psi_100) << endl;
    time = breakTime("<psi_k1|Exp[ik.r]|100>      ",time);  

//  Off-axis k vector
    double k   = 1.0; 
    double thK = PI/4;
    const double temp2[NDIM] = {0.0, k*sin(thK), k*cos(thK)};
    coordT k1_45Vec(temp2);

//  Off-Axis Positive Energy State
    Function<complexd,NDIM> psi_k1_45 = FunctionFactory<complexd,NDIM>(world).functor( 
	functorT( new ScatteringWF(1.0, 1.0, k1_45Vec) ));
    time = breakTime("Projecting |psi_k1_45>      ",time);  
    cout << "<psi_k1|Exp[ik.r]|100> = " << psi_k1_45.inner(laserOp*psi_100) << endl;
    time = breakTime("<psi_k1_45|Exp[iF.r]|100>      ",time);  

}

			     
