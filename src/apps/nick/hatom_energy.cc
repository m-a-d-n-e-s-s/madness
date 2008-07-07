/// \file hatom_energy.cc
/// \brief Compute the energy of the hydrogen atom ground state

#include <mra/mra.h>
#include <complex>
#include <gsl/gsl_sf_legendre.h>

typedef std::complex<double> complexd;

using namespace madness;
const int NL = 2;
extern "C" void coulcc_(complexd* XX, complexd* ETA1, complexd* ZLMIN,
			const int* NL, complexd* FC, complexd* GC,
			complexd* FCP, complexd* GCP, complexd* SIG,
			int* MODE1, int* KFN, int* IFAIL);

// exp[iF.r] = exp[iF*z] 
complexd expifr(const Vector<double,3>& r) {
  double g = 1.0;			//The momentum transfer from the field
  complexd i(0.0,1.0);
	//F.r = z when the field is aligned along the z axis
  return exp(i*g*r[2]);
}

// exp[-r]
complexd psi_ground(const Vector<double,3>& rr) {
  double r = sqrt(rr[0]*rr[0]+rr[1]*rr[1]+rr[2]*rr[2]);
  return 2.0*exp(-r);
}

// Positive energy radial eigenstate
complexd psi_test(double r) 
{
  if( r < 1E-15) return 0.0;
	double   Z = 1.0;
	double   k = 1.0;
	double   eta = Z/k;
  complexd XX(k*r,0.0);
  complexd ETA1(eta,0.0); // Corresponds to Z=1 and k=0.1 -> E=1/2*k^2
  complexd ZLMIN(0.0,0.0);
  complexd FC[NL];
  complexd GC[NL];
  complexd FCP[NL];
  complexd GCP[NL];
  complexd SIG[NL];
  int MODE1 = 4;			//return only F
  int KFN   = 0;		  //return complex Coulomb function
  int IFAIL = 1;			//Print Error messages 0 -> no error messages
  print("XX = ",XX);
  print("eta = ",eta);
  coulcc_(&XX, &ETA1, &ZLMIN, &NL, FC, GC, FCP, GCP, SIG, &MODE1, &KFN, 
	  &IFAIL);
  return FC[0];
}

// Positive energy radial eigenstate
complexd psi_k(const Vector<double,3>& r) 
{
  double rr = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
  if( rr < 1E-15) return 0.0;
	double   Z = 1.0;
	double   k = 1.0;
	double   eta = Z/k;
  complexd XX(k*rr,0.0);
  complexd ETA1(eta,0.0); // Corresponds to Z=1 and k=0.1 -> E=1/2*k^2
  complexd ZLMIN(0.0,0.0);
  complexd FC[NL];
  complexd GC[NL];
  complexd FCP[NL];
  complexd GCP[NL];
  complexd SIG[NL];
  int MODE1 = 4;			//return only F
  int KFN   = 0;		  //return complex Coulomb function
  int IFAIL = 1;			//Print Error messages 0 -> no error messages
  //  print("XX = ",XX);
  //print("In psi_r before coulcc_:  XX = ",XX);
  coulcc_(&XX, &ETA1, &ZLMIN, &NL, FC, GC, FCP, GCP, SIG, &MODE1, &KFN, &IFAIL);
  return FC[0]/rr;
}

int main(int argc, char**argv) {
  // Initialize the parallel programming environment
  MPI::Init(argc, argv);
  World world(MPI::COMM_WORLD);
  
  // Load info for MADNESS numerical routines
  startup(world,argc,argv);
  
  // Setup defaults for numerical functions
	double thresh = 1e-2;
  FunctionDefaults<3>::k = 11;                // Wavelet order
  FunctionDefaults<3>::thresh = thresh;         // Accuracy
  FunctionDefaults<3>::refine = true;         // Enable adaptive refinement
  FunctionDefaults<3>::initial_level = 2;     // Initial projection level
  for (int i=0; i<3; i++) {
    FunctionDefaults<3>::cell(i,0) = -20.0;   // User-defined volume
    FunctionDefaults<3>::cell(i,1) =  20.0;
  }
	/*
	for(int i=0; i<10; i++) {
		 double r = 1.0*i;
		 print("integrand = ", psi_test(r)*psi_ground_test(r)*expifr_test(r));
	}
	*/
  print("I am starting u");
  Function<complexd,3> u = FunctionFactory<complexd,3>(world).f(psi_ground);
  print("I am finished with u");
  Function<complexd,3> v = FunctionFactory<complexd,3>(world).f(expifr);
  print("I am finished with v");
  Function<complexd,3> w = FunctionFactory<complexd,3>(world).f(psi_k);
  print("I am finished with w");
	print("Polynomials up to order k = ", 11);
	print("Accuracy = ", thresh);
  complexd AFint = (w.inner(v*u));
  print("the Atomic Form integral is",AFint);
  MPI::Finalize();				//FLAG
  return 0;
}
