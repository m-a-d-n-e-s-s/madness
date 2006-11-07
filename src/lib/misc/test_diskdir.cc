#include <iostream>
using std::cout;
using std::endl;

#include <cstring>
using std::strcmp;

#include <vector>
using std::vector;

/// \file mra/test.cc

#include <mra/mra.h>
#include <misc/misc.h>
#include <misc/communicator.h>
#include <mra/twoscale.h>
#include <mra/legendre.h>
#include <tensor/tensor.h>
#include <misc/madexcept.h>
#include <misc/diskdir.h>
#include <serialize/vecar.h>

#include <unistd.h>

using namespace madness;

const double PI = 3.1415926535897932384;

double fred(double x, double y, double z) {
    double fac = pow(2.0*65.0/PI,0.75);
    x-=0.5;
    y-=0.5;
    z-=0.5;
    return fac*exp(-65.0*(x*x+y*y+z*z));
}

int main(int argc, char* argv[]) {
  
  Communicator& comm = startup(argc,argv);
  FunctionDefaults::k=15;
  FunctionDefaults::initial_level=0;

  print("INITIAL TENSOR COUNT",BaseTensor::get_instance_count());
  
  // To ensure reliable cleanup catch all C++ exceptions here
  try {

    Function<double> f = FunctionFactory<double>(fred).thresh(1e-3);
    Function<double> f2 = FunctionFactory<double>();
    DiskDir<double> test("./tdiskdir.dat", comm);
    
    //test["mo"];
    //test("mo", 1) << f;
    test("mo",1);
    test << f;
    test >> f2;
    cout << " class subtraction test " << (f - f2).norm2sq() << endl;
    vector< Function<double> > mo_test;
    mo_test.resize(5);
    for(int i = 0; i < 5; i++) {
      mo_test[i] = FunctionFactory<double>(fred); 
    };
    test << mo_test;
    vector< Function<double> > mo_test2;
    test >> mo_test2;
    unsigned int nmo_test2 = mo_test2.size();
    int nsize_mo_test2 = nmo_test2;
    for(int i = 0; i < nsize_mo_test2; i++) {
      cout << " class subtraction test " << i << " " << (mo_test[i] - mo_test2[i]).norm2sq() << endl;
    }

    DiskDir<double> test2("./tothertest.dat", comm);
    //double double_test = 10.0;
    long long_test = 10;
    test2 << long_test;
    long_test = 1;
    test2 >> long_test;
    cout << "long_test=" << long_test << endl;
    double double_test = 10.1;
    test2 << double_test;
    double_test = 1.5;
    test2 >> double_test;
    cout << "double_test=" << double_test << endl;
    bool  bool_test = true;
    test2 << bool_test;
    bool_test = false;
    test2 >> bool_test;
    cout << "bool_test=" << bool_test << endl;
    //    char char_test[256];
    /*
      char* char_test;
      test2 << "abc";
      test2 >> char_test;
      cout << "char_test=" << char_test << endl;
    */
    
    comm.Barrier();
    vector<int> vectorint_test;
    vectorint_test.resize(5);
    for(int i=0; i<5; i++) {
      vectorint_test[i] = i;
    }
    test2 << vectorint_test;
    vector<int> vectorint_test2;
    test2 >> vectorint_test2;
    int nsize = vectorint_test2.size();
    for(int i = 0; i < nsize; i++) {
      cout << vectorint_test[i] << endl;
    }
    
    print("TENSOR COUNT AT END OF SCOPE",BaseTensor::get_instance_count());
    
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
  return 0;
}

