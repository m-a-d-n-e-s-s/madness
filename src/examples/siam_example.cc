#include <madness/mra/mra.h>
#include <madness/constants.h>
#include "nonlinsol.h"
using namespace madness;
using namespace std;

const double eps_int = 1.0; ///< Interior dielectric
const double eps_ext =10.0; ///< Exterior dielectric
const double R = 2.0;       ///< Radius of sphere
const double xi = 100.0;    ///< Exponent for delta function approx.
const double sigma = 0.1;   ///< Surface "width"

/// Functor representing the log-derivative of epsilon along axis
class DSphere : public FunctionFunctorInterface<double, 3> {
  const double fac;
  const int axis;

protected:
  /// Derivative of the mask w.r.t. s=sdf/sigma
  double dmask(double s) const {
    if (fabs(s) > 6.0) return 0.0;
    return -exp(-s*s) / sqrt(constants::pi);
  }
  
public:
  DSphere(double Vint, double Vext, int axis) 
    : fac(log(Vint/Vext)), axis(axis)  {}
  
  double operator()(const coord_3d &r) const {
    double d = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    double sdf = (d-R) / sigma;
    return fac*r[axis]*dmask(sdf)/(sigma*d);
  }
};

double rho_function(const coord_3d &r) {
  const double coeff = pow(xi / constants::pi, 1.5);
  return coeff*exp(-xi * (r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));
}

int main(int argc, char **argv) {
  const double rfourpi = 1.0/(4.0*constants::pi);
  initialize(argc, argv);
  World world(SafeMPI::COMM_WORLD);
  startup(world, argc, argv);

  // Function defaults
  FunctionDefaults<3>::set_k(6);
  FunctionDefaults<3>::set_thresh(1e-4);
  FunctionDefaults<3>::set_cubic_cell(-5, 5);
  FunctionDefaults<3>::set_initial_level(3);
  FunctionDefaults<3>::set_truncate_on_project(true);

  // Make integral and derivative operators
  real_convolution_3d op = CoulombOperator(world, 1e-4, FunctionDefaults<3>::get_thresh());
  real_derivative_3d Dx = free_space_derivative<double,3>(world, 0);
  real_derivative_3d Dy = free_space_derivative<double,3>(world, 1);
  real_derivative_3d Dz = free_space_derivative<double,3>(world, 2);

  // Make functors for dielectric related quantities
  real_functor_3d epsx_functor(new DSphere(eps_int, eps_ext, 0));
  real_functor_3d epsy_functor(new DSphere(eps_int, eps_ext, 1));
  real_functor_3d epsz_functor(new DSphere(eps_int, eps_ext, 2));

  // Make the actual numerical functions
  real_function_3d rho = real_factory_3d(world).f(rho_function);
  real_function_3d eps_x = real_factory_3d(world).functor(epsx_functor);
  real_function_3d eps_y = real_factory_3d(world).functor(epsy_functor);
  real_function_3d eps_z = real_factory_3d(world).functor(epsz_functor);

  real_function_3d u = op(rho).truncate(); // Initial guess
  NonlinearSolver solver;
  for (int iter=0; iter<20; iter++) {
    real_function_3d surf_rho = rfourpi*(eps_x*Dx(u) + eps_y*Dy(u) + eps_z*Dz(u));
    real_function_3d r = (u - op(rho + surf_rho)).truncate(); // residual
    real_function_3d unew = solver.update(u, r);
    double change = (unew - u).norm2();
    u = unew;

    print("iter", iter, "change", change, "surf charge", surf_rho.trace());
    
    if (change < 10.*FunctionDefaults<3>::get_thresh()) break;
  }
  
  coord_3d lo{0.0,0.0,0.0}, hi{5.0,0.0,0.0}; // Range for line plotting
  plot_line("testpot.dat", 501, lo, hi, u);
  
  finalize();
  return 0;
}
