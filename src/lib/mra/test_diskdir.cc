#include <vector>
using namespace std;
#include <iostream>
using std::cout;
using std::endl;

#include <mra/mra.h>
#include <mra/diskdir.h>
#include <misc/misc.h>
#include <misc/communicator.h>
#include <mra/twoscale.h>
#include <mra/legendre.h>
#include <tensor/tensor.h>

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
// This is the test program for DiskDir class.
    // Before treatment
    MPI::Init(argc, argv);
    Communicator comm;
    madness::comm_default = &comm;
    redirectio(comm);
    comm.print();
    load_coeffs(comm);
    load_quadrature(comm);
    FunctionDefaults::tree = new FunctionOctTree(OctTree<FunctionNode>::create_default(comm,2));
    if (!gauss_legendre_test()) comm.Abort();
    if (!test_two_scale_coefficients()) comm.Abort();

    // test field ->
    FunctionDefaults::k=7;
    FunctionDefaults::initial_level=1;
    Function<double> f = FunctionFactory<double>(fred).k(3).refine(1).compress(1).initial_level(2).thresh(1e-1);
    f.compress();
    Function<double> g;
    DiskDir<double> testDiskDir("./tdiskdir.dat", "keep", "write", comm);
    cout << " before << operator " << endl;
    testDiskDir << f;
    cout << " after << operator " << endl;
    //testDiskDir.classClose();
    DiskDir<double> testDiskDir2("./tdiskdir.dat", "keep", "read", comm);
    testDiskDir2 >> g;
    //vector< Function<double> > mos[2];
    //testDiskDir << mos;
    //testDiskDir >> mos;
    // <-test field 

    // Finalize
    comm.close();
    MPI::Finalize();
    return 0;
}
