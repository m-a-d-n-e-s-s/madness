/// \file mra/test.cc

#include <iostream>
using std::cout;
using std::endl;
#include <fstream>
using std::ifstream;
#include <serialize/archive.h>
#include <mra/mra.h>
#include <misc/misc.h>
#include <misc/communicator.h>
#include <mra/twoscale.h>
#include <mra/legendre.h>
#include <tensor/tensor.h>

using namespace madness;

double *Cut;

const double PI           = 3.1415926535897932384;
const double THREE_SQRTPI = 5.31736155271654808184;
const double thresh       = 1e-5;
const double antoau       = 1.889725989;

struct input_data {
  int ncent;		// the number of atoms;
  double *coords;	// coordinates of atoms;
  double *charge;	// charge on each atoms;
  long k; 		// wavelet order; 
  double L;		// box size;
  double eacc;		// energy precision;
  long initial_level;	// initial level; 
};

static input_data* Geom;

double min(double a, double b)
{
    if (a < b)
	return a;
    else
	return b;
}

void set_cut(input_data& data)
{
  Geom = &data;
  Cut = new double[Geom->ncent]; 
  for(long i=0; i<Geom->ncent;i++){
//    cout << " xyz = " << i << " " << Geom->coords[3*i] << " " << Geom->coords[3*i+1] << " " << Geom->coords[3*i+2] << endl;
    Cut[i] = 
    pow(min(1e-3, Geom->eacc*0.3)/2.0/0.00435/pow(Geom->charge[i],5.0),(1.0/3.0));
  }
};

void unset_cut(){
  delete [] Cut; 
};

static double u(double r)
{
  /* Regularized 1/r */

  double r2 = r*r, pot;
  if (r > 6.5){
    pot = 1.0/r;
  }else if (r > 1e-8){
    pot = erf(r)/r + (exp(-r2) + 16*exp(-4.0*r2))/(THREE_SQRTPI);
  }else{
    pot = (2.0 + 17.0/3.0)/sqrt(PI);
  }

  return pot;
};

double V(double x, double y, double z)
{
  /* Regularized nuclear potential */
  double sum = 0.0;
  int i;

  /* Convert from [0,1] into [-L/2,L/2] and atomic units */

//  cout << " Geom->L = " << Geom->L << endl;
  x = (x-0.5)*(Geom->L);
  y = (y-0.5)*(Geom->L);
  z = (z-0.5)*(Geom->L);

  for (i=0; i<Geom->ncent; i++) {
    double xx = x-(Geom->coords[3*i]);
    double yy = y-(Geom->coords[3*i+1]);
    double zz = z-(Geom->coords[3*i+2]);
    double r = sqrt(xx*xx + yy*yy + zz*zz);

    sum += -1.0*(Geom->charge[i])*u(r/Cut[i])/Cut[i];
  }
  return sum;
};
 
double NuclearRep()
{
  double sum = 0.0, rsq;
  for (int i=0; i<Geom->ncent; i++) {
    for (int j=0; j<i; j++) {
      rsq = 0.0;
      for (int k=0; k<3; k++) {
        rsq += pow((Geom->coords[3*i+k] - Geom->coords[3*j+k]),2.0);
      }
      sum += (Geom->charge[i]*Geom->charge[j])/sqrt(rsq);
    }
  }
  return sum;
};

double fred(double x, double y, double z) {
    double fac = pow(2.0*65.0/PI,0.75);
    x-=0.5; y-=0.5; z-=0.5;
    return fac*exp(-65.0*(x*x+y*y+z*z));
}
//double fred(double x, double y, double z) {
//    return x*x+y*y*z*z;
//}

//double_complex cfred(double x, double y, double z) {
//    return x*x+y*y*z*z;
//}

int main(int argc, char* argv[]) {
    char ss[256];
    double x, y, z;
    long charge;
    // The following should be used to setup all calculations
    // 1) Initialize parallel environment
    // 2) Setup communication information
    // 3) Redirect standard output+err for parallel processes
    // 4) Load coeffs and quadrature information from file
    // 5) Setup default OctTreeLayout
    // 6) Sanity check
    MPI::Init(argc, argv);
    Communicator comm;
    redirectio(comm);
    comm.print();
    input_data testdata;
    if(comm.rank()==0){
      ifstream fin(argv[1]);
//      ifstream fin;
//      fin.open(argv[1]);
      if (!fin) {
         cout << argv[1] << " The input file cannot be opened!" << endl;
      }
      if (fin.getline(ss,256)) {
        testdata.ncent = atoi(ss); 
        cout << testdata.ncent << endl;
        testdata.coords = new double[3*testdata.ncent];
        testdata.charge = new double[testdata.ncent];
      }
      for(int i=0; i<testdata.ncent; i++) {
        if (fin.getline(ss,256)) {
          sscanf(ss,"%d %lf %lf %lf",&charge,&x,&y,&z);
          testdata.charge[i] = static_cast<double>(charge);
          testdata.coords[3*i] = x*antoau;
          testdata.coords[3*i+1] = y*antoau;
          testdata.coords[3*i+2] = z*antoau;
          //cout << testdata.charge[i] << " " << testdata.coords[3*i] << " " << testdata.coords[3*i+1] << " " << testdata.coords[3*i+2] << endl;
        }
          cout << testdata.charge[i] << " " << testdata.coords[3*i] << " " << testdata.coords[3*i+1] << " " << testdata.coords[3*i+2] << endl;
      }
      if (fin.getline(ss,256)) {
        testdata.k = atol(ss); 
        cout << testdata.k << endl;
      }
      if (fin.getline(ss,256)) {
        testdata.L = atof(ss); 
        cout << testdata.L << endl;
      }
      if (fin.getline(ss,256)) {
        testdata.eacc = atof(ss); 
        cout << testdata.eacc << endl;
      }
      if (fin.getline(ss,256)) {
        testdata.initial_level = atol(ss); 
        cout << testdata.initial_level << endl;
      }
      fin.close();
    }
    cout << "Done getting the data" << endl;
    comm.Bcast(&testdata.ncent, 1, 0);
    if (comm.rank() != 0)
    {
	testdata.coords = new double[3*testdata.ncent];
	testdata.charge = new double[testdata.ncent];
    }
    for(int i=0; i<testdata.ncent; i++) {
      cout << "in bcast loop, i = " << i << endl;
      comm.Bcast(testdata.charge[i], 0);
      comm.Bcast(testdata.coords[3*i], 0);
      comm.Bcast(testdata.coords[3*i+1], 0);
      comm.Bcast(testdata.coords[3*i+2], 0);
      cout << "end of bcast loop, i = " << i << endl;
    }
    cout << "Done broadcasting the charges and coordinates" << endl;
    comm.Bcast(testdata.ncent, 0);
    comm.Bcast(testdata.k, 0);
    comm.Bcast(testdata.L, 0);
    comm.Bcast(testdata.eacc, 0);
    comm.Bcast(testdata.initial_level, 0);
/*
*/
    cout << "Before Loading the coeffs and quadrature" << endl;
    load_coeffs(comm);
    load_quadrature(comm);
    cout << "Loaded the coeffs and quadrature" << endl;
    FunctionDefaults::tree = new FunctionOctTree(OctTree<FunctionNode>::create_default(comm,2));
    cout << "Created new tree" << endl;
    if (!gauss_legendre_test()) comm.Abort();
    cout << "Performed gauss-Legendre test" << endl;
    if (!test_two_scale_coefficients()) comm.Abort();
    cout << "Performed two-scale coefficients test" << endl;
    
    (void)set_cut(testdata);
    cout << "did set_cut of testdata" << endl;
    FunctionDefaults::k=testdata.k;
    FunctionDefaults::initial_level=testdata.initial_level;
    Function<double> f = FunctionFactory<double>(V).refine(1).compress(0).initial_level(testdata.initial_level).thresh(1e-9);
    cout << "Refined compressed whatever" << endl;
/*
*/
/*
*/
    // The follwing should be used for succesful termination
    (void)unset_cut();
    cout << "did unset_cut thingy" << endl;
    comm.close(); 
    MPI::Finalize();
    return 0;
}

