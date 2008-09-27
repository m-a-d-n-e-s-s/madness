#include <mra/mra.h>
#include <iostream>

#include "dft.h"
//#include "hartreefock.h"

using std::cout;
using std::endl;

using namespace madness;

const double PI = 3.1415926535897932384;

typedef Vector<double,3> coordT;

//  //***************************************************************************
//  template <typename T, int NDIM>
//  class NuclearChargeDensityOp : public EigSolverOp<T,NDIM>
//  {
//  public:
//    typedef Function<T,NDIM> funcT;
//    //*************************************************************************
//    // Constructor
//    NuclearChargeDensityOp(World& world, double coeff, double thresh);
//    //*************************************************************************
//
//    //*************************************************************************
//    // Is there an orbitally-dependent term?
//    virtual bool is_od() {return false;}
//    //*************************************************************************
//
//    //*************************************************************************
//    // Is there a density-dependent term?
//    virtual bool is_rd() {return true;}
//    //*************************************************************************
//
//    //*************************************************************************
//    virtual funcT op_r(const funcT& rho, const funcT& rhon, const funcT& psi)
//    {
//
//    }
//    //*************************************************************************
//
//  private:
//    //*************************************************************************
//    funcT _rhon;
//    //*************************************************************************
//  };
//  //***************************************************************************

//*****************************************************************************
static double rho_func_hydro(const coordT& r)
{
  const double x=r[0], y=r[1], z=r[2];
  double e1 = 100.0;
  double coeff = pow(e1/PI, 1.5);
  return -1.0 * coeff * exp(-e1 * (x*x + y*y + z*z));
}
//*****************************************************************************

//*****************************************************************************
double psi_func_hydro(const Vector<double,3>& r)
{
  double x = r[0];
  double y = r[1];
  double z = r[2];
  return exp(-x*x + y*y + z*z);
}
//*****************************************************************************

////*****************************************************************************
//void test_hydro(int argc, char** argv)
//{
//  MPI::Init(argc, argv);
//  World world(MPI::COMM_WORLD);
//  startup(world,argc,argv);
//
//  // Box size
//  double L = 10.0;
//
//  // Function defaults
//  int funck = 8;
//  double thresh = 1e-6;
//  FunctionDefaults<3>::set_k(funck);
//  FunctionDefaults<3>::set_thresh(thresh);
//  FunctionDefaults<3>::set_refine(true);
//  FunctionDefaults<3>::set_initial_level(2);
//  FunctionDefaults<3>::set_truncate_mode(1);
//  FunctionDefaults<3>::set_cubic_cell(-L/2, L/2);
//
//  // Nuclear potential (Be)
//  //const coordT origin(0.0);
//  if (world.rank() == 0) madness::print("Creating Function object for nuclear charge density ...");
//  Function<double,3> rhon = FunctionFactory<double,3>(world).f(rho_func_hydro).thresh(thresh);
//  rhon.truncate();
//
//  // Guess for the wavefunctions
//  if (world.rank() == 0) madness::print("Creating wavefunction's ...");
//  Function<double,3> psi = FunctionFactory<double,3>(world).f(psi_func_hydro);
//  psi.truncate();
//  psi.scale(1.0/psi.norm2());
//  std::vector<Function<double,3> > phis;
//  phis.push_back(psi);
//  // Create list of eigenvalues
//  std::vector<double> eigs;
//  eigs.push_back(-0.9);
//  // Create eigensolver
//  if (world.rank() == 0) madness::print("Creating Eigensolver object...");
//  // Constructor for non-periodic system
//  EigSolver<double,3> solver(world, phis, eigs, ops, thresh);
//  if (world.rank() == 0) madness::print("Diagonalizing Hamiltonian ...");
//  solver.solve(15);
//
//  double eval = solver.get_eig(0);
//  Function<double,3> func = solver.get_phi(0);
//  if (world.rank() == 0) printf("reconstructing func ...\n\n");
//
//  MPI::Finalize();
//}
////*****************************************************************************

#define TO_STRING(s) TO_STRING2(s)
#define TO_STRING2(s) #s

//*****************************************************************************
int main(int argc, char** argv)
{
  MPI::Init(argc, argv);
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
//    test_hf_he(world);
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
  MPI::Finalize();

  return 0;
}
//*****************************************************************************
