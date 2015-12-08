#include <madness/mra/mra.h>
#include <madness/constants.h>
#include "nonlinsol.h"
using namespace madness;
using namespace std;

const double eps_int = 1.0; // Interior dielectric
const double eps_ext =10.0; // Exterior dielectric
const double R = 2.0; // Radius of cavity
const double expnt = 100.0; // Exponent of Gaussian approx to delta function
const double sigma = 0.1; // Surface width

class DSphere : public FunctionFunctorInterface<double,3>  {
  double R, sigma, fac;
  int axis;
protected:
  double radius(const coord_3d& r) const {
    return sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
  }
  
  double dmask(double s) const { // Derivative of the mask w.r.t. s=sdf/sigma
    const double rsqrtpi = 1.0/sqrt(madness::constants::pi);
    if (fabs(s) > 5.5) return 0.0;
    return -exp(-s*s)*rsqrtpi;
  }
  
  coord_3d gradient(const coord_3d& r) const { // Gradient of mask
    double d = radius(r);
    double sdf = (d-R)/sigma;
    return r*(dmask(sdf)/(sigma*d));
  }
  
public:
  DSphere(double radius, double sigma, double Vint, double Vext, int axis) 
    : R(radius), sigma(sigma), fac(log(Vint/Vext)), axis(axis)
  {}
  
  double operator()(const madness::coord_3d& r) const {
    return fac*gradient(r)[axis];
  }
};


double charge_function(const coord_3d& r) {
    const double coeff = pow(1.0/constants::pi*expnt,1.5);
    return coeff*exp(-expnt*(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));
}

double exact_function(const coord_3d& x) {
    double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    if (r > R) return 1.0/(eps_ext*r);
    else return erf(sqrt(expnt)*r)/(eps_int*r) + (1.0/eps_ext - 1.0/eps_int)/R;
}

int main(int argc, char **argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);

    // Function defaults
    FunctionDefaults<3>::set_k(10);
    FunctionDefaults<3>::set_thresh(1e-8);
    FunctionDefaults<3>::set_cubic_cell(-5, 5);
    FunctionDefaults<3>::set_initial_level(3);
    FunctionDefaults<3>::set_truncate_on_project(true);

    // Integral and derivative operators
    real_convolution_3d op = CoulombOperator(world, 1e-4, FunctionDefaults<3>::get_thresh());
    real_derivative_3d Dx = free_space_derivative<double,3>(world, 0);
    real_derivative_3d Dy = free_space_derivative<double,3>(world, 1);
    real_derivative_3d Dz = free_space_derivative<double,3>(world, 2);

    // Functors for dielectric related quantities
    real_functor_3d gradx_functor(new DSphere(R, sigma, eps_int, eps_ext, 0));
    real_functor_3d grady_functor(new DSphere(R, sigma, eps_int, eps_ext, 1));
    real_functor_3d gradz_functor(new DSphere(R, sigma, eps_int, eps_ext, 2));

    // Make the actual functions
    real_function_3d exact = real_factory_3d(world).f(exact_function);
    real_function_3d charge = real_factory_3d(world).f(charge_function);
    real_function_3d eps_x = real_factory_3d(world).functor(gradx_functor);
    real_function_3d eps_y = real_factory_3d(world).functor(grady_functor);
    real_function_3d eps_z = real_factory_3d(world).functor(gradz_functor);

    const double rfourpi = 1.0/(4.0*constants::pi);

    real_function_3d u = op(charge).truncate();
    NonlinearSolver solver;
    for (int iter=0; iter<20; iter++) {
      real_function_3d surf_charge = rfourpi*(eps_x*Dx(u) + eps_y*Dy(u) + eps_z*Dz(u));
      real_function_3d r = (u - op(charge + surf_charge)).truncate();
      real_function_3d unew = solver.update(u, r);
      double change = (unew-u).norm2();
      u = unew;

      print("iter", iter, "change", change, "surf charge", surf_charge.trace());
      
      if (change < 10.0*FunctionDefaults<3>::get_thresh()) break;
    }
    
    coord_3d lo{-5.0,0.0,0.0}, hi{0.5,0.0,0.0}; // Range for line plotting
    plot_line("testpot.dat", 301, lo, hi, u, exact);
    
    finalize();
    return 0;
}
