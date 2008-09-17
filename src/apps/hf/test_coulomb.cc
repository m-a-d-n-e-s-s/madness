#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <constants.h>
#include "eigsolver.h"
#include "dft.h"
//#include "outputwriter.h"

using namespace madness;
using std::endl;

#define WST_PI madness::constants::pi

//#define DEBUG_STREAM *(OutputWriter::instance()->debug_stream())
//#define LOG_STREAM *(OutputWriter::instance()->log_stream())

#define DEBUG_STREAM std::cout
#define LOG_STREAM std::cout

typedef Vector<double,3> coordT3d;
typedef Vector<double,1> coordT1d;

const double L = 30.0;
const double alpha = 13.5;

/// Returns radius for smoothing nuclear potential with energy precision eprec
//*****************************************************************************
static double smoothing_parameter(double Z, double eprec) {
    // The min is since asymptotic form not so good at low acc.
    // The 2 is from two electrons in 1s closed shell.
    if (Z == 0.0) return 1.0;
    double Z5 = Z*Z*Z*Z*Z;
    double c = pow(std::min(1e-3,eprec)/2.0/0.00435/Z5,1.0/3.0);
    return c;
}
//*****************************************************************************


/// Regularized 1/r potential.
//*****************************************************************************
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
        pot = (2.0 + 17.0/3.0)/sqrt(WST_PI);
    }

    return pot;
}
//*****************************************************************************

//*****************************************************************************
static double V_cosine1D(const coordT1d& r)
{
  const double x=r[0];
  double twopi = 2 * WST_PI;
  return  -alpha*(cos(twopi*x/L)+1);
}
//*****************************************************************************

//*****************************************************************************
static double V_cosine3D(const coordT3d& r)
{
  const double x=r[0]; const double y = r[1]; const double z = r[2];
  double twopi = 2 * WST_PI;
  return  -alpha*(cos(twopi*x/L)*cos(twopi*y/L)*cos(twopi*z/L) + 1);
}
//*****************************************************************************

//*****************************************************************************
static double V_cosine3D_sep(const coordT3d& r)
{
  const double x=r[0]; const double y = r[1]; const double z = r[2];
  double twopi = 2 * WST_PI;
  return  -alpha*(cos(twopi*x/L) + cos(twopi*y/L) + cos(twopi*z/L) + 3);
}
//*****************************************************************************

//*****************************************************************************
template <int NDIM>
static double gaussian_func(const Vector<double, NDIM> r)
{
  double sum = 0.0;
  for (int i = 0; i < NDIM; i++)
    sum += r[i]*r[i];
  return exp(-2*L*sum);
}
//*****************************************************************************

//*****************************************************************************
static double phi_initial_guess1D(const coordT1d& r)
{
  return gaussian_func<1>(r);
}
//*****************************************************************************

//*****************************************************************************
static double phi_initial_guess3D(const coordT3d& r)
{
  return gaussian_func<3>(r);
}
//*****************************************************************************

//*****************************************************************************
static double psi_func_he(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
  return 6.0*exp(-2.0*sqrt(x*x+y*y+z*z)+1e-4);
}
//*****************************************************************************

//*****************************************************************************
static double V_func_he(const coordT3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
  double rr = sqrt(x*x + y*y + z*z);
  double c = smoothing_parameter(2.0, 1e-7);
  return -2.0 * smoothed_potential(rr/c) / c;
}
//*****************************************************************************


//*****************************************************************************
void test_cosine_1D(int argc, char** argv)
{
  MPI::Init(argc, argv);
  World world(MPI::COMM_WORLD);
  startup(world,argc,argv);

  // Function defaults
  int funck = 6;
  double thresh = 1e-4;
  FunctionDefaults<1>::set_k(funck);
  FunctionDefaults<1>::set_thresh(thresh);
  FunctionDefaults<1>::set_refine(true);
  FunctionDefaults<1>::set_initial_level(2);
  FunctionDefaults<1>::set_truncate_mode(1);
  FunctionDefaults<1>::set_cubic_cell(-L/2, L/2);

  // Nuclear potential (Be)
  //const coordT origin(0.0);
  if (world.rank() == 0) madness::print("Creating Function object for nuclear potential ...");
  Function<double,1> Vnuc = FunctionFactory<double,1>(world).f(V_cosine1D).thresh(thresh);

  // Guess for the wavefunctions
  if (world.rank() == 0) madness::print("Creating wavefunction's ...");
  Function<double,1> psi = FunctionFactory<double,1>(world).f(phi_initial_guess1D);
  psi.scale(1.0/psi.norm2());
  std::vector<Function<double,1> > phis;
  phis.push_back(psi);
  // Create list of eigenvalues
  std::vector<double> eigs;
  eigs.push_back(-5.0);
  // Create list of operators
  std::vector<EigSolverOp<double,1>*> ops;
  ops.push_back(new DFTNuclearPotentialOp<double,1>(world, Vnuc, 1.0, thresh));
  // Create eigensolver
  if (world.rank() == 0) madness::print("Creating Eigensolver object...");
  // Constructor for non-periodic system
  EigSolver<double,1> solver(world, phis, eigs, ops, thresh);
  if (world.rank() == 0) madness::print("Diagonalizing Hamiltonian ...");
  solver.solve(11);

  MPI::Finalize();
}
//*****************************************************************************

//*****************************************************************************
void test_cosine_3D(int argc, char** argv)
{
  MPI::Init(argc, argv);
  World world(MPI::COMM_WORLD);
  startup(world,argc,argv);

  // Function defaults
  int funck = 8;
  double thresh = 1e-6;
  FunctionDefaults<3>::set_k(funck);
  FunctionDefaults<3>::set_thresh(thresh);
  FunctionDefaults<3>::set_refine(true);
  FunctionDefaults<3>::set_initial_level(2);
  FunctionDefaults<3>::set_truncate_mode(1);
  FunctionDefaults<3>::set_cubic_cell(-L/2, L/2);

  // Nuclear potential (Be)
  //const coordT origin(0.0);
  if (world.rank() == 0) madness::print("Creating Function object for nuclear potential ...");
  Function<double,3> Vnuc = FunctionFactory<double,3>(world).f(V_cosine3D_sep).thresh(thresh);

  // Guess for the wavefunctions
  if (world.rank() == 0) madness::print("Creating wavefunction's ...");
  Function<double,3> psi = FunctionFactory<double,3>(world).f(phi_initial_guess3D);
  psi.scale(1.0/psi.norm2());
  std::vector<Function<double,3> > phis;
  phis.push_back(psi);
  // Create list of eigenvalues
  std::vector<double> eigs;
  eigs.push_back(-50.0);
  // Create list of operators
  std::vector<EigSolverOp<double,3>*> ops;
  ops.push_back(new DFTNuclearPotentialOp<double,3>(world, Vnuc, 1.0, thresh));
  // Create eigensolver
  if (world.rank() == 0) madness::print("Creating Eigensolver object...");
  // Constructor for non-periodic system
  EigSolver<double,3> solver(world, phis, eigs, ops, thresh);
  if (world.rank() == 0) madness::print("Diagonalizing Hamiltonian ...");
  solver.solve(45);

  double eval = solver.get_eig(0);
  Function<double,3> func = solver.get_phi(0);
  if (world.rank() == 0) printf("reconstructing func ...\n\n");
  func.reconstruct();
  coordT3d pt1(0.3);
  coordT3d pt2(1.3);
  if (world.rank() == 0) printf("evaluating points ...\n\n");
  double funcpt1 = func(pt1);
  double funcpt2 = func(pt2);
  if (world.rank() == 0) printf("eval = %.8f\n\n", eval);
  if (world.rank() == 0) printf("func(pt1) = %.8f\tfunc(pt2) = %.8f\n\n", funcpt1, funcpt2);

  MPI::Finalize();
}
//*****************************************************************************

//*****************************************************************************
void test_he(int argc, char** argv)
{
  MPI::Init(argc, argv);
  World world(MPI::COMM_WORLD);
  startup(world,argc,argv);

  // Setup output files
//  OutputWriter::instance()->init_debug("test_he_debug.txt");
//  OutputWriter::instance()->init_log("test_he_log.txt");

  // Function defaults
  int funck = 6;
  double thresh = 1e-4;
  FunctionDefaults<3>::set_k(funck);
  FunctionDefaults<3>::set_thresh(thresh);
  FunctionDefaults<3>::set_refine(true);
  FunctionDefaults<3>::set_initial_level(2);
  FunctionDefaults<3>::set_truncate_mode(1);
  FunctionDefaults<3>::set_cubic_cell(-L/2, L/2);

  // Nuclear potential (Be)
  //const coordT origin(0.0);
  if (world.rank() == 0) LOG_STREAM << "Creating Function object for nuclear potential ...\n" << endl;
  Function<double,3> Vnuc = FunctionFactory<double,3>(world).f(V_func_he).thresh(thresh);

  // Guess for the wavefunctions
  if (world.rank() == 0) LOG_STREAM << "Creating wavefunction's ...\n" << endl;
  Function<double,3> psi = FunctionFactory<double,3>(world).f(psi_func_he);
  psi.scale(1.0/psi.norm2());
  std::vector<Function<double,3> > phis;
  phis.push_back(psi);
  // Create list of eigenvalues
  std::vector<double> eigs;
  eigs.push_back(-5.0);
  // Create dft object
  if (world.rank() == 0) LOG_STREAM << "Creating DFT object..." << endl;
  // Constructor for non-periodic system
  DFT<double,3> dft(world, Vnuc, phis, eigs, thresh, true);
  if (world.rank() == 0) LOG_STREAM << "Running DFT object..." << endl;
  dft.solve(31);

  MPI::Finalize();
}
//*****************************************************************************

int main(int argc, char** argv)
{
  //test_he(argc, argv);
  test_cosine_3D(argc, argv);
  return 0;
}
