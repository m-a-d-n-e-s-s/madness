#include <iostream>
using std::cout;
using std::endl;
#include <cmath>
using namespace std;
#include <fstream>
using namespace std;

#include <algorithm>

/// \file mra/test.cc


#include <mra/mra.h>
#include <misc/misc.h>
#include <misc/communicator.h>
#include <mra/twoscale.h>
#include <mra/legendre.h>
#include <tensor/tensor.h>


#include <misc/madexcept.h>
#include <misc/print.h>
#include <misc/madpapi.h>

#define endl "\n";


using namespace madness;


const double PI           = 3.1415926535897932384;
const double THREE_SQRTPI = 5.31736155271654808184;
const double antoau       = 1.889725989;

const int MAXPROC = 12000;
const int MAXCENT = 4096;
double Cut[MAXCENT];
struct input_data {
  double coords[MAXCENT*3];	// coordinates of atoms;
  double charge[MAXCENT];	// charge on each atoms;
  int ncent;		// the number of atoms;
  long k; 		// wavelet order; 
  double Lx;		// box size;
  double Ly;		// box size;
  double Lz;		// box size;
  double eacc;		// energy precision;
  long initial_level;	// initial level; 
};

static input_data* Geom;

void set_cut(input_data& data)
{
  Geom = &data;
  for(long i=0; i<Geom->ncent;i++){
//    cout << " xyz = " << i << " " << Geom->coords[3*i] << " " << Geom->coords[3*i+1] << " " << Geom->coords[3*i+2] << endl;
    Cut[i] = 
    pow(min(1e-3, 1e-5*0.3)/2.0/0.00435/pow(Geom->charge[i],5.0),(1.0/3.0));
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

inline void tosim(double& x, double& y, double& z) {x= x/Geom->Lx + 0.5;y= y/Geom->Ly + 0.5;z= z/Geom->Lz + 0.5;};
inline void tomol(double& x, double& y, double& z) {x= (x-0.5)*Geom->Lx;y= (y-0.5)*Geom->Ly;z= (z-0.5)*Geom->Lz;};

double V(double x, double y, double z)
{
  /* Regularized nuclear potential */
  double sum = 0.0;
  int i;

  /* Convert from [0,1] into [-L/2,L/2] and atomic units */

  tomol(x,y,z);

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
        rsq += pow((Geom->coords[3*i+k] - Geom->coords[3*j+k]),2);
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
    Communicator& comm = startup(argc,argv);
/*
    MPI::Init(argc, argv);
    Communicator comm;
    redirectio(comm);
    comm.print();
*/
    input_data testdata;
    if(comm.rank()==0){
      ifstream fin;
      fin.open(argv[1]);
      if (!fin) {
         cout << argv[1] << " The input file cannot be opened!" << endl;
      }
      if (fin.getline(ss,256)) {
        testdata.ncent = atoi(ss); 
        cout << testdata.ncent << endl;
      }
      double xx=0.0, yy=0.0, zz=0.0;
      for(int i=0; i<testdata.ncent; i++) {
        if (fin.getline(ss,256)) {
          sscanf(ss,"%ld %lf %lf %lf",&charge,&x,&y,&z);
          testdata.charge[i] = charge;
          testdata.coords[3*i] = x*antoau;
          testdata.coords[3*i+1] = y*antoau;
          testdata.coords[3*i+2] = z*antoau;
          xx += x*antoau;
          yy += y*antoau;
          zz += z*antoau;
        }
        else {
            MADNESS_EXCEPTION("BAD INPUT",i);
        }
      }
      // Recenter the coordinates and also compute bounding box
      xx/=testdata.ncent; yy/=testdata.ncent; zz/=testdata.ncent;
      double bx=0.0, by=0.0, bz=0.0;
      for(int i=0; i<testdata.ncent; i++) {
          testdata.coords[3*i  ] -= xx;
          testdata.coords[3*i+1] -= yy;
          testdata.coords[3*i+2] -= zz;
          bx = std::max(bx,std::abs(testdata.coords[3*i]));
          by = std::max(by,std::abs(testdata.coords[3*i+1]));
          bz = std::max(bz,std::abs(testdata.coords[3*i+2]));

          //cout << testdata.charge[i] << " " << testdata.coords[3*i] << " " << testdata.coords[3*i+1] << " " << testdata.coords[3*i+2] << endl;
      }
      if (fin.getline(ss,256)) {
        testdata.k = atol(ss); 
        cout << testdata.k << endl;
      }
      if (fin.getline(ss,256)) {
        double L = atof(ss); 
        cout << L << endl;
      }
      bx*=2.0; by*=2.0; bz*=2.0;
      print("bounding box",bx,by,bz);
      testdata.Lx = bx;
      testdata.Ly = by;
      testdata.Lz = bz;
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
    comm.Bcast(testdata.ncent, 0);
    for(int i=0; i<testdata.ncent; i++) {
      comm.Bcast(testdata.charge[i], 0);
      comm.Bcast(testdata.coords[3*i], 0);
      comm.Bcast(testdata.coords[3*i+1], 0);
      comm.Bcast(testdata.coords[3*i+2], 0);
    }
    comm.Bcast(testdata.k, 0);
    comm.Bcast(testdata.Lx, 0);
    comm.Bcast(testdata.Ly, 0);
    comm.Bcast(testdata.Lz, 0);
    comm.Bcast(testdata.eacc, 0);  
    comm.Bcast(testdata.initial_level, 0);
    
    (void)set_cut(testdata);

    FunctionDefaults::k=testdata.k;
    FunctionDefaults::initial_level=testdata.initial_level;
    FunctionDefaults::k = testdata.k;
    FunctionDefaults::thresh = testdata.eacc;

    print("INITIAL LEVEL",testdata.initial_level,testdata.k,testdata.eacc,V(0.5,0.6,0.7));
    

    comm.Barrier();

    papi::start();
    Function<double> f = FunctionFactory<double>(V).refine().nocompress();
    papi::stop("Initial projection+refine");

    int natom[MAXPROC];
    for(int i=0; i<comm.nproc(); i++) natom[i] = 0;

    double minx, miny, minz;
    double maxx, maxy, maxz;
    minx = miny = minz = 999999.0;
    maxx = maxy = maxz = -99999.0;
    for(int i=0; i<testdata.ncent; i++) {
        double x = testdata.coords[3*i];
        double y = testdata.coords[3*i+1];
        double z = testdata.coords[3*i+2];
        minx = std::min(minx,x);
        miny = std::min(miny,y);
        minz = std::min(minz,z);
        maxx = std::max(maxx,x);
        maxy = std::max(maxy,y);
        maxz = std::max(maxz,z);
    }
    print("MINMAX mol",minx, miny, minz, maxx, maxy, maxz);
    minx = miny = minz = 999999.0;
    maxx = maxy = maxz = -99999.0;
    for(int i=0; i<testdata.ncent; i++) {
        double x = testdata.coords[3*i];
        double y = testdata.coords[3*i+1];
        double z = testdata.coords[3*i+2];
        tosim(x,y,z);
        minx = std::min(minx,x);
        miny = std::min(miny,y);
        minz = std::min(minz,z);
        maxx = std::max(maxx,x);
        maxy = std::max(maxy,y);
        maxz = std::max(maxz,z);
    }
    print("MINMAX sim",minx, miny, minz, maxx, maxy, maxz);

    print("process that owns each atom");
    for(int i=0; i<testdata.ncent; i++) {
        double value;
        double x = testdata.coords[3*i];
        double y = testdata.coords[3*i+1];
        double z = testdata.coords[3*i+2];
        tosim(x,y,z);
        bool flag = f.eval_local(x,y,z,&value);
        ProcessID owner = flag ? comm.rank() : 0;
        owner = comm.global_sum(owner);
        natom[owner]++;
        value = comm.global_sum(value);
        print(i,"->",owner,x,y,z,value);
    }
    print("no. of atoms per process");
    int nnn = 0;
    for (int i=0; i<comm.nproc(); i++) {
        print(i,"->",natom[i]);
        nnn += natom[i];
    }
    print("CONSERVATION OF ATOMS?",nnn);
    cout.flush();

    if( std::abs( f(0.33,0.44,0.55) - V(0.33,0.44,0.55)) > 3e-3 ) {
      cout << " test1 failed" << endl;
    }
    cout << " test1 passed" << endl;

    papi::start();
    f.compress();
    papi::stop("Compress");
    print("normsq after compression   ",f.norm2sq());

    papi::start();
    f.reconstruct();
    papi::stop("Reconstruct");
    print("normsq after reconstruction",f.norm2sq());

    papi::start();
    f*f;
    papi::stop("Multiplication");

    papi::start();
    f.square();
    papi::stop("Squaring");

    papi::start();
    f.diff(0);
    papi::stop("Differentiation");

    f.compress();
    f.pnorms();

    print("BARRIER at end");
    comm.Barrier();
    // The follwing should be used for succesful termination
    comm.close(); 
    //(void)unset_cut();
    MPI::Finalize();
    return 0;
}
