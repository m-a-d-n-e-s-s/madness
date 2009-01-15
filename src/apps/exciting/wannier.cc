#include <mra/mra.h>
#include <iostream>
#include <cmath>

extern "C" void readinput_();
extern "C" void wann_init1_();
extern "C" void wann_unk_(int* n,double* vpl,double* vrc,double* val);

using namespace madness;

typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;
typedef FunctionFactory<double,3> factoryT;
typedef Function<double,3> functionT;
typedef Vector<double,3> coordT;

double func(int n, const Vector<int,3>& nk, double xx, double yy, double zz)
{
  // initialize value to zero
  double val[2];
  val[0] = 0.0; val[1] = 0.0;
  // number of k-points
  int nkpt = nk[0] * nk[1] * nk[2];
  //compute \sum_{k}u_{nk}(r)
  double vkl[3];
  // translate point
  double vr[3];
  vr[0] = xx; vr[1] = yy; vr[2] = zz;
  for (int j1 = 0; j1 < nk[0]; j1++)
  {
    for (int j2 = 0; j2 < nk[1]; j2++)
    {
      for (int j3 = 0; j3 < nk[2]; j3++)
      {
        // temp variable
        double val0[2];
        //k-point in lattice coordinates
        vkl[0] = 1.0 * j1 / nk[0];
        vkl[1] = 1.0 * j2 / nk[1];
        vkl[2] = 1.0 * j3 / nk[2];
        //get value of u_{nk}(r)
        wann_unk_(&n, vkl, vr, val0);
        val[0] += val0[0];
        val[1] += val0[1];
      }
    }
  }
  val[0]=val[0]/nkpt;
  val[1]=val[1]/nkpt;

  return sqrt(val[0]*val[0] + val[1]*val[1]);
}

template<typename T, int NDIM>
class Wannier: public FunctionFunctorInterface<T, NDIM>
{
public:
  typedef Vector<double, NDIM> coordT;
  int _n;
  coordT _center;
  Vector<int,3> _nk;

  Wannier(int n, const coordT& center, const Vector<int,3>& nk) :
    _n(n), _center(center), _nk(nk) {}

  T operator()(const coordT& x) const
  {
    return func(_n, _nk, x[0], x[1], x[2]);
  }
};

void test_wannier(World& world)
{
    //k-mesh division
    Vector<int,3> nk(4);
    //cener point of the box
    Vector<double,3> center(0.0);
    //index of Wannier function
    int n=3;

    readinput_();
    wann_init1_();

    // Function defaults
    int funck = 5;
    double thresh = 1e-3;
    double bsize = 4.0;
    FunctionDefaults<3>::set_k(funck);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_cubic_cell(-bsize, bsize);

    functorT functor(new Wannier<double,3>(n, center, nk));
    functionT w = factoryT(world).functor(functor);

    {
      w.reconstruct();
      if (world.rank() == 0)  printf("\n");
      double L = bsize / 2;
      double bstep = L / 100.0;
      for (int i = 0; i < 101; i++)
      {
        coordT p(-L / 2 + i * bstep);
        double fval = func(n, nk, p[0], p[1], p[2]);
        double fdiff = w(p) - fval;
        if (world.rank() == 0)
          printf("%10.2f%15.8f%15.8f%15.8f\n", p[0], w(p), fval, fdiff);
      }
    }

    // Plot to OpenDX
    vector<long> npt(3,101);
    plotdx(w, "wannier.dx", FunctionDefaults<3>::get_cell(), npt);

}

#define TO_STRING(s) TO_STRING2(s)
#define TO_STRING2(s) #s

//*****************************************************************************
int main(int argc, char** argv)
{
  initialize(argc, argv);
  World world(MPI::COMM_WORLD);
  if (world.rank() == 0)
  {
    print("");
    print("--------------------------------------------");
    print("   MADNESS", " multiresolution testsuite");
    print("--------------------------------------------");
    print("");
    print("   number of processors ...", world.size());
    print("    processor frequency ...", cpu_frequency());
    print("            host system ...", TO_STRING(HOST_SYSTEM));
    print("             byte order ...", TO_STRING(MADNESS_BYTE_ORDER));
    print("          configured by ...", MADNESS_CONFIGURATION_USER);
    print("          configured on ...", MADNESS_CONFIGURATION_HOST);
    print("          configured at ...", MADNESS_CONFIGURATION_DATE);
    print("                    CXX ...", MADNESS_CONFIGURATION_CXX);
    print("               CXXFLAGS ...", MADNESS_CONFIGURATION_CXXFLAGS);
#ifdef WORLD_WATCHDOG
    print("               watchdog ...", WATCHDOG_BARK_INTERVAL,
        WATCHDOG_TIMEOUT);
#endif
#ifdef OPTERON_TUNE
    print("             tuning for ...", "opteron");
#elif defined(CORE_DUO_TUNE)
    print("             tuning for ...", "core duo");
#else
    print("             tuning for ...", "core2");
#endif
#ifdef BOUNDS_CHECKING
    print(" tensor bounds checking ...", "enabled");
#endif
#ifdef TENSOR_INSTANCE_COUNT
    print("  tensor instance count ...", "enabled");
#endif
    print(" ");
  }

  try
  {
    printf("WSTHORNTON: Starting up the world ... \n");

    startup(world,argc,argv);
    if (world.rank() == 0) print("Initial tensor instance count", BaseTensor::get_instance_count());
    test_wannier(world);
  }
  catch (const MPI::Exception& e)
  {
    //        print(e);
    error("caught an MPI exception");
  }
  catch (const madness::MadnessException& e)
  {
    print(e);
    error("caught a MADNESS exception");
  }
  catch (const madness::TensorException& e)
  {
    print(e);
    error("caught a Tensor exception");
  }
  catch (const char* s)
  {
    print(s);
    error("caught a string exception");
  }
  catch (const std::string& s)
  {
    print(s);
    error("caught a string (class) exception");
  }
  catch (const std::exception& e)
  {
    print(e.what());
    error("caught an STL exception");
  }
  catch (...)
  {
    error("caught unhandled exception");
  }

  if (world.rank() == 0)
    print("entering final fence");
  world.gop.fence();
  if (world.rank() == 0)
    print("done with final fence");
  if (world.rank() == 0)
    print("Final tensor instance count", BaseTensor::get_instance_count());

  finalize();
  return 0;
}

