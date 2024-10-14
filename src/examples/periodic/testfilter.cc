#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

const double L = 10;

static double gaussian_3d(const madness::coord_3d& r, const double expnt) {
    const double x=r[0], y=r[1], z=r[2];
    const double coeff = pow(expnt/M_PI, 1.5);
    return coeff * exp(-expnt * (x*x + y*y + z*z));
}

static double rho_electronic_3d(const madness::coord_3d& r) { return gaussian_3d(r, 1); }
static double rho_nuclear_3d(const madness::coord_3d& r) { return -gaussian_3d(r, 10000); }

static double rho_gaussian_func_3d(const madness::coord_3d& r) {
    return rho_electronic_3d(r) + rho_nuclear_3d(r);
}

template<typename T, std::size_t NDIM>
void filter_moments_inplace(madness::Function<T,NDIM>& f, const int k, const bool fence=true) {
    auto filter_op = [&](madness::Tensor<T> &coeff) {
        if constexpr (NDIM == 3) {
            if (k >= 0) {
                coeff(0, 0, 0) = 0.;
            }
            if (k >= 1) {
                coeff(1, 0, 0) = 0.;
                coeff(0, 1, 0) = 0.;
                coeff(0, 0, 1) = 0.;
            }
            if (k >= 2) {
                coeff(2, 0, 0) = 0.;
                coeff(0, 2, 0) = 0.;
                coeff(0, 0, 2) = 0.;
                coeff(1, 1, 0) = 0.;
                coeff(1, 0, 1) = 0.;
                coeff(0, 1, 1) = 0.;
            }
            if (k >= 3)
                abort();// TODO implement support for higher moments
        } else
            static_assert("filter_moments_inplace at the moment is implemented for NDIM=3");
    };

    // on [-L,L] normalized scaling function k is Sqrt[(2 k + 1)/(2 L)] LegendreP[k, x/L]
    // the coefficient of x^k in LegendreP[k, x] is 2^(-k) Binomial[2 k, k], hence
    // coefficient of x^k in scaling function k is Sqrt[(2 k + 1)/(2 L)] Binomial[2 k, k]/(2 L)^k
    // thus if all moments of up to k-1 vanish, k-th moment (=expectation value of x^k) of
    // scaling function k is its coefficient times Sqrt[(2 L)/ (2 k + 1)] (2 L)^k / Binomial[2 k, k]
    f.unaryop_coeff([&](const madness::Key<NDIM> &key, madness::Tensor<T> &coeff) {
        if (f.is_reconstructed()) {
            filter_op(coeff);
        } else if (f.is_compressed()) {
            // for compressed form only need 1 node, but public interface does not allow mutable access to the data
            if (key == f.get_impl()->key0()) {
                filter_op(coeff);
            }
        } else {
            MADNESS_EXCEPTION("filter_moments_inplace(f): f must be either compressed or reconstructed", 1);
        }
    },
                    fence);
}

int main(int argc, char**argv) {
    using namespace madness;
    {
      auto &world = initialize(argc, argv);
      startup(world, argc, argv, true);
      
      // Function defaults
      int k = 11;
      double thresh = 1e-9;
      double eps = 1e-9;
      FunctionDefaults<3>::set_k(k);
      FunctionDefaults<3>::set_cubic_cell(-L / 2, L / 2);
      FunctionDefaults<3>::set_thresh(thresh);
      FunctionDefaults<3>::set_refine(true);
      FunctionDefaults<3>::set_initial_level(2);
      FunctionDefaults<3>::set_truncate_mode(1);
      
      Function<double, 3> rho, rhoE, rhoN, V_periodic, V_periodicE, V_periodicN;
      
      {
        constexpr double filter_level = 0;
	printf("building gaussian diff charge distribution ...\n\n");
	std::vector<coord_3d> special_pt{coord_3d(std::array<double,3>{0.0, 0.0, 0.0})};
	rho = FunctionFactory<double, 3>(world).special_points(special_pt).initial_level(4).f(rho_gaussian_func_3d);
        filter_moments_inplace(rho, filter_level);
	rho.truncate();
	rhoE = FunctionFactory<double, 3>(world).special_points(special_pt).initial_level(4).f(rho_electronic_3d);
        filter_moments_inplace(rhoE, filter_level);
	rhoE.truncate();
	rhoN = FunctionFactory<double, 3>(world).special_points(special_pt).initial_level(4).f(rho_nuclear_3d);
        filter_moments_inplace(rhoN, filter_level);
	rhoN.truncate();
	
	BoundaryConditions<3> bc_periodic(BC_PERIODIC);
	SeparatedConvolution<double, 3> pop = CoulombOperator(world, 1e-4, eps, bc_periodic);
	printf("applying periodic operator ...\n\n");
	V_periodic = apply(pop, rho);
	V_periodic.truncate();
	V_periodicE = apply(pop, rhoE);
	V_periodicE.truncate();
	V_periodicN = apply(pop, rhoN);
	V_periodicN.truncate();
	
	V_periodic.reconstruct();
	V_periodicE.reconstruct();
	V_periodicN.reconstruct();
      }
      
      double bstep = L / 1000.0;
      printf("   z\t\tV_(E+N)[z]\tV_E[z]\t\tV_N[z]\t\tV_E[z]+V_n[z]\trho\t\terror\n");
      for (int i = 0; i < 1001; i++) {
	double x = -L / 2 + i * bstep;
        coord_3d p(std::array<double, 3>{0,0,x});
        printf("%.3f\t\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", p[2], V_periodic(p), V_periodicE(p), V_periodicN(p), V_periodicE(p) + V_periodicN(p), rho(p), V_periodic(p) - V_periodicN(p) - V_periodicE(p));
      }
    }
    madness::finalize();
}
