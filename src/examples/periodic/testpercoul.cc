#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

const double L = 10;

// Normalized Gaussian exponent a produces potential erf(sqrt(a)*r)/r
// This function is a difference of a tight and diffuse Gaussian, has no net moments, so decay is determined
// by the diffuse Gaussian
static double rho_gaussian_func3d(const madness::coord_3d& r)
{
    const double x=r[0], y=r[1], z=r[2];
    /* diffuse Gaussian ~ electronic density */
    const double expnt1 = 1.0;
    const double coeff1 = pow(expnt1/M_PI, 1.5);
    double piece1 = coeff1 * exp(-expnt1 * (x*x + y*y + z*z));
    /* tight Gaussian ~ nuclear density */
    const double expnt2 = 10000.0;
    const double coeff2 = pow(expnt2/M_PI, 1.5);
    double piece2 = coeff2 * exp(-expnt2 * (x*x + y*y + z*z));

    return piece1-piece2;
}

// This function test both the periodic and non-periodic versions of the Coulomb
// operator. In order to make this test valid set L to a high value so that
// charge distribution should not be able to see its neighbor.
int main(int argc, char**argv) {
    using namespace madness;
    {
      auto &world = initialize(argc, argv);
      startup(world, argc, argv, true);
      
      // Function defaults
      int k = 10;
      double thresh = 1e-8;
      double eps = 1e-8;
      FunctionDefaults<3>::set_k(k);
      FunctionDefaults<3>::set_cubic_cell(-L / 2, L / 2);
      FunctionDefaults<3>::set_thresh(thresh);
      FunctionDefaults<3>::set_refine(true);
      FunctionDefaults<3>::set_initial_level(2);
      FunctionDefaults<3>::set_truncate_mode(1);
      
      Function<double, 3> V_nonperiodic(world), V_periodic;
      
      {
	FunctionDefaults<3>::set_bc(BC_PERIODIC);
	printf("building gaussian diff charge distribution ...\n\n");
	auto special_pt = std::vector<coord_3d>{coord_3d(std::array<double,3>{0.0, 0.0, 0.0})};
	Function<double, 3> rho = FunctionFactory<double, 3>(world).special_points(special_pt).initial_level(4).f(rho_gaussian_func3d);
	rho.truncate();
	
	BoundaryConditions<3> bc_periodic(BC_PERIODIC);
	SeparatedConvolution<double, 3> pop = CoulombOperator(world, 1e-4, eps, bc_periodic);
	printf("applying periodic operator ...\n\n");
	V_periodic = apply(pop, rho);
	V_periodic.truncate();
	
	V_periodic.reconstruct();
	V_nonperiodic.reconstruct();
      }
      
      double bstep = L / 100.0;
      printf("   z\t\tV_per[z]\tV_nonper[z]\terror\n");
      for (int i = 0; i < 101; i++) {
	double x = -L / 2 + i * bstep;
        coord_3d p(std::array<double, 3>{x,x,x});
        double error = fabs(V_periodic(p) - V_nonperiodic(p));
        printf("%.2f\t\t%.8f\t%.8f\t%.8f\n", p[0], V_periodic(p), V_nonperiodic(p), error);
      }
    }
    madness::finalize();
}
