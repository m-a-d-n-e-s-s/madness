#include<vector>
using namespace std;

#include <cstdlib>

#include <iostream>
using std::cout;
using std::endl;

#include <cmath>

#include <misc/communicator.h>
using madness::Communicator;

#include <mra/twoscale.h>
#include <mra/legendre.h>

/// \file serialize/test.cc
/// \brief Tests serialization by some of the archives

#define ARCHIVE_REGISTER_TYPE_INSTANTIATE_HERE

#include <serialize/textfsar.h>
using madness::archive::TextFstreamInputArchive;
using madness::archive::TextFstreamOutputArchive;

#include <serialize/binfsar.h>
using madness::archive::BinaryFstreamInputArchive;
using madness::archive::BinaryFstreamOutputArchive;

#include <serialize/vecar.h>
using madness::archive::VectorInputArchive;
using madness::archive::VectorOutputArchive;

#include <mra/mra.h>
using madness::Function;
//using madness::Function::save_local;
using madness::FunctionDefaults;
using madness::FunctionFactory;
using madness::FunctionOctTree;
using madness::FunctionNode;

#include <octtree/octtree.h>
using madness::OctTree;
using madness::OctTreeT;

#include <misc/misc.h>
using madness::redirectio;

const double PI = 3.1415926535897932384;

double fred(double x, double y, double z) {
    double fac = std::pow(2.0*65.0/PI,0.75);
    x-=0.5;
    y-=0.5;
    z-=0.5;
    return fac*std::exp(-65.0*(x*x+y*y+z*z));
}


double expnt, xxx, yyy, zzz;
double ranfred(double x, double y, double z) {
    double fac = std::pow(2.0*expnt/PI,0.75);
    x=x-xxx;
    y=y-yyy;
    z=z-zzz;
    return fac*std::exp(-expnt*(x*x+y*y+z*z));
}

double drand() {
  return random()/double(RAND_MAX);
}

void ranfn() {
  expnt = drand()*10000;
  xxx = drand();
  yyy = drand();
  zzz = drand();
  print("RANFN",expnt,xxx,yyy,zzz);
}


int main(int argc, char* argv[]) {
  
  Communicator& comm = startup(argc,argv);

  try {

    ranfn();
    Function<double> ftest = FunctionFactory<double>(ranfred).k(7).thresh(1e-5);
    Function<double> ftest2 = FunctionFactory<double>().k(7).thresh(1e-5);
    //Function<double> ftest2 = FunctionFactory<double>();
    const char* f = "tserialize.dat";
    //TextFstreamOutputArchive oar(f);
    ftest.compress();
    ftest.truncate();
    //ftest.save_local(oar);
    long partLevel = 2;
    ftest.save(f, partLevel);
    //TextFstreamInputArchive iar(f);
    //ftest.load_local(iar);
    ftest2.load(f);
    cout << " class subtraction test " << (ftest - ftest2).norm2sq() << endl;

  } catch (char const* msg) {
    std::cerr << "Exception (string): " << msg << std::endl;
    comm.Abort();
  } catch (std::exception& e) {
    std::cerr << "Exception (std): " << e.what() << std::endl;
    comm.Abort();
  } catch (MadnessException& e) {
    std::cerr << e << std::endl;
    comm.Abort();
  } catch (TensorException& e) {
    std::cerr << e << std::endl;
    comm.Abort();
  } catch (MPI::Exception& e) {
    std::cerr << "Exception (mpi): code=" << e.Get_error_code()
	      << ", class=" << e.Get_error_class()
	      << ", string=" << e.Get_error_string() << std::endl;
    comm.Abort();
  } catch (...) {
    std::cerr << "Exception (general)" << std::endl;
    comm.Abort();
  }

  // The following should be used for succesful termination
 done:
  
  print("TENSOR COUNT OUTSIDE OF SCOPE",BaseTensor::get_instance_count());
  
  comm.close();
  MPI::Finalize();
  
  return 1;
}
