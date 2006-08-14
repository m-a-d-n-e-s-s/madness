/// \file mra/test.cc

#include <fstream>
using std::ifstream;
#include <iostream>
using std::cout;
using std::endl;
#include <mra/mraX.h>
#include <misc/misc.h>
#include <misc/communicator.h>
#include <mra/twoscale.h>
#include <mra/legendre.h>
#include <tensor/tensor.h>
//#include <octtree/sendrecv.h>

using namespace madness;

double *Cut;

const double PI           = 3.1415926535897932384;
const double THREE_SQRTPI = 5.31736155271654808184;
const double thresh       = 1e-5;
const double antoau       = 1.889725989;

struct input_data {
    int ncent;            // the number of atoms;
    double *coords;       // coordinates of atoms;
    double *charge;       // charge on each atoms;
    long k;               // wavelet order;
    double L;             // box size;
    double eacc;          // energy precision;
    long initial_level;   // initial level;
};

static input_data* Geom;

#define min(a,b) ((a) < (b) ? (a) : (b))

void set_cut(input_data& data)
{
    Geom = &data;
    Cut = new double[Geom->ncent];
    for(long i = 0; i < Geom->ncent; i++) {
        Cut[i] = pow(min(1e-3, Geom->eacc*0.3)/2.0/0.00435/pow(Geom->charge[i],5.0),(1.0/3.0));
    }
};

void unset_cut() {
  delete [] Cut;
};

static double u(double r)
{
  /* Regularized 1/r */

    double r2 = r*r, pot;
    if (r > 6.5) {
        pot = 1.0/r;
    } else if (r > 1e-8) {
        pot = erf(r)/r + (exp(-r2) + 16*exp(-4.0*r2))/(THREE_SQRTPI);
    } else {
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

    x = (x-0.5)*(Geom->L);
    y = (y-0.5)*(Geom->L);
    z = (z-0.5)*(Geom->L);

    for (i = 0; i < Geom->ncent; i++) {
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
    for (int i = 0; i < Geom->ncent; i++) {
        for (int j = 0; j < i; j++) {
            rsq = 0.0;
            for (int k = 0; k < 3; k++) {
                rsq += pow((Geom->coords[3*i+k] - Geom->coords[3*j+k]),2.0);
            }
            sum += (Geom->charge[i]*Geom->charge[j])/sqrt(rsq);
        }
    }
    return sum;
};

double fred(double x, double y, double z) {
    const double PI = 3.14;
    double fac = pow(2.0*65.0/PI,0.75);
    x-=0.5; y-=0.5; z-=0.5;
    return fac*exp(-65.0*(x*x+y*y+z*z));
}

int main(int argc, char* argv[]) {
    // The following should be used to setup all calculations
    // 1) Initialize parallel environment
    // 2) Setup communication information
    // 3) Redirect standard output+err for parallel processes
    // 4) Load coeffs and quadrature information from file
    // 5) Setup default OctTreeLayout
    // 6) Sanity check

    MPI::Init(argc, argv);

    Communicator comm;
    madness::comm_default = &comm;
    redirectio(comm);
    comm.print();
//    comm.set_debug(true);
    cout << "beginning of main" << endl;

    char ss[256];
    double x, y, z;
    long charge;

    input_data testdata;
    if (comm.rank() == 0)
    {
	double s0, e0;
        s0 = MPI::Wtime();
      	ifstream fin(argv[1]);
      	if (!fin) {
            cout << argv[1] << " The input file cannot be opened!" << endl;
      	}
        if (fin.getline(ss,256)) {
            testdata.ncent = atoi(ss);
            testdata.coords = new double[3*testdata.ncent];
            testdata.charge = new double[testdata.ncent];
         }
      	for (int i=0; i<testdata.ncent; i++) {
            if (fin.getline(ss,256)) {
          	sscanf(ss,"%d %lf %lf %lf",&charge,&x,&y,&z);
          	testdata.charge[i] = static_cast<double>(charge);
          	testdata.coords[3*i] = x*antoau;
          	testdata.coords[3*i+1] = y*antoau;
          	testdata.coords[3*i+2] = z*antoau;
            }
      	}
      	if (fin.getline(ss,256)) {
            testdata.k = atol(ss);
      	}
      	if (fin.getline(ss,256)) {
            testdata.L = atof(ss);
      	}
      	if (fin.getline(ss,256)) {
            testdata.eacc = atof(ss);
      	}
      	if (fin.getline(ss,256)) {
            testdata.initial_level = atol(ss);
      	}
      	fin.close();
      	e0 = MPI::Wtime();
        cout << "Time to load data: " << e0 - s0 << endl;
    }

    cout << "Done getting the data" << endl;

    // Now, broadcast the data to all processors

    double s1, e1;

    s1 = MPI::Wtime();

    comm.Bcast(&testdata.ncent, 1, 0);

    if (comm.rank() != 0)
    {
        testdata.coords = new double[3*testdata.ncent];
        testdata.charge = new double[testdata.ncent];
    }
    for(int i=0; i<testdata.ncent; i++) {
      	comm.Bcast(testdata.charge[i], 0);
      	comm.Bcast(testdata.coords[3*i], 0);
      	comm.Bcast(testdata.coords[3*i+1], 0);
      	comm.Bcast(testdata.coords[3*i+2], 0);
    }
    comm.Bcast(testdata.ncent, 0);
    comm.Bcast(testdata.k, 0);
    comm.Bcast(testdata.L, 0);
    comm.Bcast(testdata.eacc, 0);
    comm.Bcast(testdata.initial_level, 0);
    e1 = MPI::Wtime();
    cout << "Time to bcast data: " << e1 - s1 << endl;

    set_cut(testdata);


    load_coeffs(comm);
    load_quadrature(comm);

    // Now, create the original tree

    OctTree<FunctionNode> *t = new OctTree<FunctionNode>();

    std::vector<SharedPtr<OctTree<FunctionNode> > > treeList;
    std::cout << "created tree list" << std::endl;
    treeList.push_back(SharedPtr<OctTree<FunctionNode> >(t->create_default(comm,1)));
    std::cout << "pushed back tree list" << std::endl;
    FunctionDefaults::tree = SharedPtr<FunctionOctTree> (new FunctionOctTree(treeList));
    std::cout << "made tree" << std::endl;

    // Test coefficients and stuff

    if (!gauss_legendre_test()) comm.Abort();
    std::cout << "tested gauss-legendre" << std::endl;
    if (!test_two_scale_coefficients()) comm.Abort();
    std::cout << "tested two scale coeffs" << std::endl;

    FunctionDefaults::k=7;
    FunctionDefaults::initial_level=2;
    std::cout << "set function defaults k and initial_level" << std::endl;

    MPI::COMM_WORLD.Barrier();

//    xterm_debug(comm,0,0);

    try {
	std::cout << "about to create function f" << std::endl;

	MPI::COMM_WORLD.Barrier();
	double t1 = MPI::Wtime();
//	Function<double> f = FunctionFactory<double>(fred).thresh(1e-5).compress(0);
//	Function<double> f = FunctionFactory<double>(fred).norefine().thresh(1e-2).compress(0);
//	Function<double> f = FunctionFactory<double>(fred).thresh(1e-2).compress(0);

//	Function<double> f = FunctionFactory<double>(V).refine(1).compress(0).initial_level(testdata.initial_level).thresh(testdata.eacc);
//	Function<double> f = FunctionFactory<double>(V).compress(0).thresh(testdata.eacc);
	Function<double> f = FunctionFactory<double>(V).compress(0).thresh(1e-4);
	MPI::COMM_WORLD.Barrier();
	double t2 = MPI::Wtime();
	std::cout << "created function" << std::endl;

	print("Tree in scaling function form");
	f.pnorms();

	std::cout << "about to compress" << std::endl;
	MPI::COMM_WORLD.Barrier();
	double t3 = MPI::Wtime();
	f.compress();
	MPI::COMM_WORLD.Barrier();
	double t4 = MPI::Wtime();

	print("Tree in wavelet form");
	f.pnorms();

	std::cout << "about to reconstruct" << std::endl;
	MPI::COMM_WORLD.Barrier();
	double t5 = MPI::Wtime();
	f.reconstruct();
	MPI::COMM_WORLD.Barrier();
	double t6 = MPI::Wtime();

	//Now, do load balancing
	std::cout << "before load balancing: " << std::endl;
    	FunctionDefaults::tree->depthFirstTraverse();

	std::cout << "about to load balance" << std::endl;
	MPI::COMM_WORLD.Barrier();
	double t7 = MPI::Wtime();
	balanceFunctionOctTree(FunctionDefaults::tree);
	MPI::COMM_WORLD.Barrier();
	double t8 = MPI::Wtime();

//	f.data->trees->depthFirstTraverse();
	double cost = f.computeMaxCost();
	std::cout << "cost = " << cost << std::endl;
//	Function<double> g = FunctionFactory<double>(V).compress(0).thresh(testdata.eacc);
	Function<double> g = FunctionFactory<double>(V).compress(0).thresh(1e-4);
	MPI::COMM_WORLD.Barrier();
	double t85 = MPI::Wtime();
//	g.data->trees->depthFirstTraverse();

	std::cout << "about to compress" << std::endl;
	MPI::COMM_WORLD.Barrier();
	double t9 = MPI::Wtime();
	g.compress();
	MPI::COMM_WORLD.Barrier();
	double t10 = MPI::Wtime();

	print("Tree in wavelet form");
	g.pnorms();

	std::cout << "about to reconstruct" << std::endl;
	MPI::COMM_WORLD.Barrier();
	double t11 = MPI::Wtime();
	g.reconstruct();
	MPI::COMM_WORLD.Barrier();
	double t12 = MPI::Wtime();

	// Now, set up load balancing weighting the level exponentially
	FunctionDefaults::tree->initDepthCost(2.0);
//	FunctionDefaults::tree->initDepthCost(8.0);
	std::cout << "reinitialized cost function" << std::endl;
	MPI::COMM_WORLD.Barrier();
	double t13 = MPI::Wtime();
	Function<double> h = FunctionFactory<double>(V).compress(0).thresh(1e-4);
	MPI::COMM_WORLD.Barrier();
	double t14 = MPI::Wtime();

	MPI::COMM_WORLD.Barrier();
	double t15 = MPI::Wtime();
	balanceFunctionOctTree(FunctionDefaults::tree);
	MPI::COMM_WORLD.Barrier();
	double t16 = MPI::Wtime();
	std::cout << "finished load balance" << std::endl;

/*
	FunctionDefaults::tree->initDepthCost(2.0);
	std::cout << "doing second load balance just for kicks" << std::endl;
	balanceFunctionOctTree(FunctionDefaults::tree);
*/

	std::cout << "about to compress" << std::endl;
	MPI::COMM_WORLD.Barrier();
	double t17 = MPI::Wtime();
	h.compress();
	MPI::COMM_WORLD.Barrier();
	double t18 = MPI::Wtime();

	std::cout << "about to reconstruct" << std::endl;
	MPI::COMM_WORLD.Barrier();
	double t19 = MPI::Wtime();
	h.reconstruct();
	MPI::COMM_WORLD.Barrier();
	double t20 = MPI::Wtime();

	// Now, set up load balancing weighting the level exponentially
	FunctionDefaults::tree->initDepthCost(8.0);
	std::cout << "reinitialized cost function" << std::endl;
	MPI::COMM_WORLD.Barrier();
	double t21 = MPI::Wtime();
	Function<double> s = FunctionFactory<double>(V).compress(0).thresh(1e-4);
	MPI::COMM_WORLD.Barrier();
	double t22 = MPI::Wtime();

	MPI::COMM_WORLD.Barrier();
	double t23 = MPI::Wtime();
	balanceFunctionOctTree(FunctionDefaults::tree);
	MPI::COMM_WORLD.Barrier();
	double t24 = MPI::Wtime();
	std::cout << "finished load balance" << std::endl;

	std::cout << "about to compress" << std::endl;
	MPI::COMM_WORLD.Barrier();
	double t25 = MPI::Wtime();
	s.compress();
	MPI::COMM_WORLD.Barrier();
	double t26 = MPI::Wtime();

	std::cout << "about to reconstruct" << std::endl;
	MPI::COMM_WORLD.Barrier();
	double t27 = MPI::Wtime();
	s.reconstruct();
	MPI::COMM_WORLD.Barrier();
	double t28 = MPI::Wtime();

	std::cout << "GAME OVER!!!" << std::endl;

	std::cout << std::endl << "Timings: " << std::endl;
	std::cout << "           Initialize    Load Balance    Compress     Reconstruct" << std::endl;
	std::cout << "-----------------------------------------------------------------" << std::endl;
	std::cout << "Default       " << t2-t1 << "    ------         " << t4-t3 << "       " << t6-t5 
		<< std::endl;
	std::cout << "Balanced (1)  " << t85-t8 << "    " << t8-t7 << "        " << t10-t9 << "       " 
		<< t12-t11 << std::endl;
	std::cout << "Balanced (2)  " << t14-t13 << "    " << t16-t15 << "        " << t18-t17 << "       " 
		<< t20-t19 << std::endl;
	std::cout << "Balanced (8)  " << t22-t21 << "    " << t24-t23 << "        " << t26-t25 << "       " 
		<< t28-t27 << std::endl;

    }
    catch (char const* msg) {
        std::cerr << "Exception (string): " << msg << std::endl;
        comm.Abort();
    }
    catch (std::exception& e) {
        std::cerr << "Exception (std): " << e.what() << std::endl;
        comm.Abort();
    }
    catch (TensorException& e) {
        std::cerr << e << std::endl;
        comm.Abort();
    }
    catch (MPI::Exception& e) {
        std::cerr << "Exception (mpi): code=" << e.Get_error_code() 
                  << ", class=" << e.Get_error_class() 
                  << ", string=" << e.Get_error_string() << std::endl;
        comm.Abort();
    }
    catch (...) {
        std::cerr << "Exception (general)" << std::endl;
        comm.Abort();
    }
    
    comm.close(); 
    MPI::Finalize();
    return 0;
}

