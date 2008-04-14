#include <mra/mra.h>
#include <iostream>

#include "hartreefock.h"

using std::cout;
using std::endl;

using namespace madness;

const double PI = 3.1415926535897932384;

/// Returns radius for smoothing nuclear potential with energy precision eprec
static double smoothing_parameter(double Z, double eprec) {
    // The min is since asymptotic form not so good at low acc.
    // The 2 is from two electrons in 1s closed shell.
    if (Z == 0.0) return 1.0;
    double Z5 = Z*Z*Z*Z*Z;
    double c = pow(std::min(1e-3,eprec)/2.0/0.00435/Z5,1.0/3.0);
    return c;
}


/// Regularized 1/r potential.

/// Invoke as \c u(r/c)/c where \c c is the radius of the
/// smoothed volume.
static double smoothed_potential(double r) {
    const double THREE_SQRTPI = 5.31736155271654808184;
    double r2 = r*r, pot;
    if (r > 6.5){
        pot = 1.0/r;
    } else if (r > 1e-8){
        pot = erf(r)/r + (exp(-r2) + 16.0*exp(-4.0*r2))/(THREE_SQRTPI);
    } else{
        pot = (2.0 + 17.0/3.0)/sqrt(PI);
    }
    
    return pot;
}

//*****************************************************************************
static double psi_func_be1(const coordT& rr)
{
  const double x=rr[0], y=rr[1], z=rr[2];
  double r = sqrt(x*x+y*y+z*z);
  return exp(-4.0*r+1e-4);
}
//*****************************************************************************

//*****************************************************************************
static double psi_func_be2(const coordT& rr)
{
  const double x=rr[0], y=rr[1], z=rr[2];
  double r = sqrt(x*x+y*y+z*z);
  return (1.0 - 2.0*r*exp(-2.0*r));
}
//*****************************************************************************

//*****************************************************************************
static double psi_func_he(const coordT& r)
{
  const double x=r[0], y=r[1], z=r[2];
  return 6.0*exp(-2.0*sqrt(x*x+y*y+z*z)+1e-4);
}
//*****************************************************************************

//*****************************************************************************
static double V_func_he(const coordT& r)
{
  const double x=r[0], y=r[1], z=r[2];
  return -2.0/(sqrt(x*x+y*y+z*z)+1e-6);
}
//*****************************************************************************

//*****************************************************************************
static double V_func_be(const coordT& r)
{
  const double x=r[0], y=r[1], z=r[2];
  double rr = sqrt(x*x + y*y + z*z);
  double c = smoothing_parameter(4.0, 1e-7);
  return -4.0 * smoothed_potential(rr/c) / c;
}
//*****************************************************************************

//*****************************************************************************
double V_func_h2(const Vector<double,3>& r)
{
  double x = r[0];
  double y = r[1];
  double z = r[2];
  return -1.0/(sqrt(x*x + y*y + (z-0.7)*(z-0.7)) + 1e-8) + 
    -1.0/(sqrt(x*x + y*y + (z+0.7)*(z+0.7)) + 1e-8);
}
//*****************************************************************************

//*****************************************************************************
double psi_func_h2(const Vector<double,3>& r)
{
  double x = r[0];
  double y = r[1];
  double z = r[2];
//  return exp(-0.5*(x*x + y*y + (z-0.7)*(z-0.7))) + 
//    exp(-0.5*(x*x + y*y + (z+0.7)*(z+0.7)));
  return exp(-sqrt(x*x + y*y + (z-0.7)*(z-0.7))) + 
    exp(-sqrt(x*x + y*y + (z+0.7)*(z+0.7)));
}
//*****************************************************************************

//*****************************************************************************
template<typename T, int NDIM> class Gaussian :
  public FunctionFunctorInterface<T,NDIM>
{
public:
  typedef Vector<double,NDIM> coordT;
  const coordT center;
  const double exponent;
  const T coefficient;

  Gaussian(const coordT& center, double exponent, T coefficient) :
    center(center), exponent(exponent), coefficient(coefficient)
  {
  }
  
  T operator()(const coordT& x) const
  {
    double sum = 0.0;
    for (int i=0; i<NDIM; i++)
    {
      double xx = center[i]-x[i];
      sum += xx*xx;
    };
    return coefficient*exp(-exponent*sum);
  }
  
};
//*****************************************************************************

//*****************************************************************************
class H2Potential :
  public FunctionFunctorInterface<double,3>
{
public:
  typedef Vector<double,3> coordT;
  
  H2Potential() {}

  double operator()(const coordT& x) const
  {
    double xx = x[0];
    double yy = x[1];
    double zz = x[2];
    return -1.0/(sqrt(xx*xx + yy*yy + (zz-0.7)*(zz-7.0)) + 1e-08) +
      -1.0/(sqrt(xx*xx + yy*yy + (zz+0.7)*(zz+7.0)) + 1e-08);
  }
};
//*****************************************************************************

//*****************************************************************************
template<typename T, int NDIM> class HarmonicOsc3D :
  public FunctionFunctorInterface<T,NDIM>
{
public:
  typedef Vector<double,NDIM> coordT;
  const coordT _center;
  const double _coefficient;
  const double _offset;

  //***************************************************************************
  HarmonicOsc3D(const coordT& center, double coefficient, double offset) :
    _center(center), _coefficient(coefficient), _offset(offset)
  {
  }
  //***************************************************************************

  //***************************************************************************
  T operator()(const coordT& x) const
  {
    double sum = 0.0;
    for (int i=0; i<NDIM; i++)
    {
      double xx = _center[i]-x[i];
      sum += xx*xx;
    };
    return (_coefficient * sum) + _offset;
  }
  //***************************************************************************

};
//*****************************************************************************

//*****************************************************************************
void test_hf_ho(World& world)
{
  cout << "Running test application HartreeFock ..." << endl;
  
  typedef Vector<double,3> coordT;
  typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;

  // Dimensions of the bounding box
  double bsize = 10.0;
  for (int i=0; i<3; i++)
  {
    FunctionDefaults<3>::cell(i,0) = -bsize;
    FunctionDefaults<3>::cell(i,1) = bsize;
  }
  // Function defaults
  FunctionDefaults<3>::k = 5;
  FunctionDefaults<3>::thresh = 1e-3;
  FunctionDefaults<3>::refine = true;
  FunctionDefaults<3>::initial_level = 2;
  
  // Nuclear potential (harmonic oscillator)
  const coordT origin(0.0);
  const double coeff = 0.5;
  const double offset = -50.0;
  functorT Vnuc_functor(new HarmonicOsc3D<double,3>(origin, coeff, offset));
  Function<double,3> Vnuc = FunctionFactory<double,3>(world).functor(Vnuc_functor);
  
  // Guess for the wavefunction
  functorT wavefunc_functor(new Gaussian<double,3>(origin, -0.5, 100.0));
  Function<double,3> psi = FunctionFactory<double,3>(world).functor(Vnuc_functor);
  psi.scale(1.0/psi.norm2());
  printf("Norm of psi = %.5f\n\n", psi.norm2());
  // Create HartreeFock object
  cout << "Creating HartreeFock object..." << endl;
  HartreeFock hf(world, Vnuc, psi, -42.5, false, false, 1e-5);
  cout << "Running HartreeFock object..." << endl;
  hf.hartree_fock(10);
  printf("Ground state is: %.5f\n", hf.get_eig(0));

}
//*****************************************************************************

//*****************************************************************************
void test_hf_h2(World& world)
{
  cout << "Running test application HartreeFock ..." << endl;
  
  typedef Vector<double,3> coordT;
  typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;

  // Dimensions of the bounding box
  double bsize = 30.0;
  for (int i=0; i<3; i++)
  {
    FunctionDefaults<3>::cell(i,0) = -bsize;
    FunctionDefaults<3>::cell(i,1) = bsize;
  }
  // Function defaults
  FunctionDefaults<3>::k = 7;
  FunctionDefaults<3>::thresh = 1e-5;
  FunctionDefaults<3>::refine = true;
  FunctionDefaults<3>::initial_level = 2;
  
  // Nuclear potential (harmonic oscillator)
  const coordT origin(0.0);
  cout << "Creating Function object for nuclear potential ..." << endl;
  Function<double,3> Vnuc = FunctionFactory<double,3>(world).f(V_func_h2);
 
  // Guess for the wavefunction
  cout << "Creating wavefunction psi ..." << endl;
  Function<double,3> psi = FunctionFactory<double,3>(world).f(psi_func_h2);
  psi.scale(1.0/psi.norm2());
  printf("Norm of psi = %.5f\n\n", psi.norm2());
  // Create HartreeFock object
  cout << "Creating HartreeFock object..." << endl;
  HartreeFock hf(world, Vnuc, psi, -0.6, true, true, 1e-5);
  cout << "Running HartreeFock object..." << endl;
  hf.hartree_fock(10);
//  double ke = 2.0 * hf.calculate_tot_ke_sp();
//  double pe = 2.0 * hf.calculate_tot_pe_sp();
//  double ce = hf.calculate_tot_coulomb_energy();
//  double ee = hf.calculate_tot_exchange_energy();
//  double ne = 1.0/1.4;
//  printf("Kinetic energy:\t\t\t %.8f\n", ke);
//  printf("Potential energy:\t\t %.8f\n", pe);
//  printf("Two-electron energy:\t\t %.8f\n", 2.0*ce - ee);
//  printf("Total energy:\t\t\t %.8f\n", ke + pe + 2.0*ce - ee + ne);
}
//*****************************************************************************

//*****************************************************************************
void test_hf_he(World& world)
{
  cout << "Running test application HartreeFock ..." << endl;
  
  typedef Vector<double,3> coordT;
  typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;

  // Dimensions of the bounding box
  double bsize = 22.4;
  for (int i=0; i<3; i++)
  {
    FunctionDefaults<3>::cell(i,0) = -bsize;
    FunctionDefaults<3>::cell(i,1) = bsize;
  }
  // Function defaults
  FunctionDefaults<3>::k = 10;
  FunctionDefaults<3>::thresh = 1e-8;
  FunctionDefaults<3>::refine = true;
  FunctionDefaults<3>::initial_level = 2;
  
  // Nuclear potential (harmonic oscillator)
  const coordT origin(0.0);
  cout << "Creating Function object for nuclear potential ..." << endl;
  Function<double,3> Vnuc = FunctionFactory<double,3>(world).f(V_func_he);
 
  // Guess for the wavefunction
  cout << "Creating wavefunction psi ..." << endl;
  Function<double,3> psi = FunctionFactory<double,3>(world).f(psi_func_he);
  psi.scale(1.0/psi.norm2());
  printf("Norm of psi = %.5f\n\n", psi.norm2());
  // Create HartreeFock object
  cout << "Creating HartreeFock object..." << endl;
  HartreeFock hf(world, Vnuc, psi, -0.6, true, true, 1e-5);
  cout << "Running HartreeFock object..." << endl;
  hf.hartree_fock(10);
//  double ke = 2.0 * hf.calculate_tot_ke_sp();
//  double pe = 2.0 * hf.calculate_tot_pe_sp();
//  double ce = hf.calculate_tot_coulomb_energy();
//  double ee = hf.calculate_tot_exchange_energy();
//  printf("Kinetic energy:\t\t\t %.8f\n", ke);
//  printf("Potential energy:\t\t %.8f\n", pe);
//  printf("Two-electron energy:\t\t %.8f\n", 2.0*ce - ee);
//  printf("Total energy:\t\t\t %.8f\n", ke + pe + 2.0*ce - ee);
}
//*****************************************************************************

//*****************************************************************************
void test_hf_be(World& world)
{
  //if (world.rank() == 0) cout << "Running test application HartreeFock ..." << endl;
  
  typedef Vector<double,3> coordT;
  typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;

  // Dimensions of the bounding box
  double bsize = 40.0;
  for (int i=0; i<3; i++)
  {
    FunctionDefaults<3>::cell(i,0) = -bsize;
    FunctionDefaults<3>::cell(i,1) = bsize;
  }
  // Function defaults
  int funck = 8;
  double thresh = 1e-6;
  FunctionDefaults<3>::k = funck;
  FunctionDefaults<3>::thresh = thresh;
  FunctionDefaults<3>::refine = true;
  FunctionDefaults<3>::initial_level = 2;
  FunctionDefaults<3>::truncate_mode = 1;
  
  // Nuclear potential (Be)
  const coordT origin(0.0);
  if (world.rank() == 0) cout << "Creating Function object for nuclear potential ..." << endl;
  cout << "Creating Function object for nuclear potential ..." << endl;
  Function<double,3> Vnuc = FunctionFactory<double,3>(world).f(V_func_be).thresh(thresh);
 
  // Guess for the wavefunctions
  if (world.rank() == 0) cout << "Creating wavefunction's ..." << endl;
  cout << "Creating wavefunction's ..." << endl;
  Function<double,3> psi1 = FunctionFactory<double,3>(world).f(psi_func_be1);
  psi1.scale(1.0/psi1.norm2());
  //if (world.rank() == 0) printf("Norm of psi1 = %.5f\n\n", psi1.norm2());
  printf("Norm of psi1 = %.5f\n\n", psi1.norm2());
  Function<double,3> psi2 = FunctionFactory<double,3>(world).f(psi_func_be2);
  psi2.scale(1.0/psi2.norm2());
  //if (world.rank() == 0) printf("Norm of psi2 = %.5f\n\n", psi2.norm2());
  printf("Norm of psi2 = %.5f\n\n", psi2.norm2());
  // Create list of wavefunctions
  std::vector<funcT> phis;
  phis.push_back(psi1);
  phis.push_back(psi2);
  // Creat list of eigenvalues
  std::vector<double> eigs;
  eigs.push_back(-5.0);
  eigs.push_back(-0.5);
  // Create HartreeFock object
  if (world.rank() == 0) cout << "Creating HartreeFock object..." << endl;
  cout << "Creating HartreeFock object..." << endl;
  HartreeFock hf(world, Vnuc, phis, eigs, true, true, thresh);
  if (world.rank() == 0) cout << "Running HartreeFock object..." << endl;
  cout << "Running HartreeFock object..." << endl;
  hf.hartree_fock(20);
//  double ke = 2.0 * hf.calculate_tot_ke_sp();
//  double pe = 2.0 * hf.calculate_tot_pe_sp();
//  double ce = hf.calculate_tot_coulomb_energy();
//  double ee = hf.calculate_tot_exchange_energy();
//  printf("Kinetic energy:\t\t\t %.8f\n", ke);
//  printf("Potential energy:\t\t %.8f\n", pe);
//  printf("Two-electron energy:\t\t %.8f\n", 2.0*ce - ee);
//  printf("Total energy:\t\t\t %.8f\n", ke + pe + 2.0*ce - ee);
}
//*****************************************************************************

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
    print("   MADNESS", MADNESS_VERSION, "multiresolution testsuite");
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
    test_hf_be(world);
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
