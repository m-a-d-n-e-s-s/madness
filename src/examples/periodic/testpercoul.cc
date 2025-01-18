#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

const int L = 20;
const int Lx = L;
const int Ly = L;
const int Lz = L;

// N.B. Normalized Gaussian exponent a produces potential erf(sqrt(a)*r)/r

/// unit-normalized constant
static double unit_func3d(const madness::coord_3d& r)
{
  return 1/(L*L*L);
}

/// diffuse gaussian (a=1) at 0,0,0
static double gdiffuse_000_func3d(const madness::coord_3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
  /* diffuse Gaussian ~ electronic density */
  const double expnt1 = 1.0;
  const double coeff1 = pow(expnt1/M_PI, 1.5);
  return coeff1 * exp(-expnt1 * (x*x + y*y + z*z));
}

/// tight gaussian (a=1e4) at 0,0,0
static double gtight_000_func3d(const madness::coord_3d& r)
{
  const double x=r[0], y=r[1], z=r[2];
  /* diffuse Gaussian ~ electronic density */
  const double expnt1 = 1.e4;
  const double coeff1 = pow(expnt1/M_PI, 1.5);
  return coeff1 * exp(-expnt1 * (x*x + y*y + z*z));
}

/// diffuse gaussian (a=1) at 0,0,1
static double gdiffuse_001_func3d(const madness::coord_3d& r)
{
  const double x=r[0], y=r[1], z=r[2]-1.;
  /* diffuse Gaussian ~ electronic density */
  const double expnt1 = 1.0;
  const double coeff1 = pow(expnt1/M_PI, 1.5);
  return coeff1 * exp(-expnt1 * (x*x + y*y + z*z));
}

/// diffuse gaussian (a=1) at 0,0,-1
static double gdiffuse_00m1_func3d(const madness::coord_3d& r)
{
  const double x=r[0], y=r[1], z=r[2]+1.;
  /* diffuse Gaussian ~ electronic density */
  const double expnt1 = 1.0;
  const double coeff1 = pow(expnt1/M_PI, 1.5);
  return coeff1 * exp(-expnt1 * (x*x + y*y + z*z));
}

/// tight gaussian (a=1e4) at 0,0,-1
static double gtight_00m1_func3d(const madness::coord_3d& r)
{
  const double x=r[0], y=r[1], z=r[2]+1.;
  /* diffuse Gaussian ~ electronic density */
  const double expnt1 = 1.e4;
  const double coeff1 = pow(expnt1/M_PI, 1.5);
  return coeff1 * exp(-expnt1 * (x*x + y*y + z*z));
}

/// `K`-th moment function along `Axis`, w.r.t to the center of the cell
template <std::size_t Axis, std::size_t K = 1>
static double kth_moment_func3d(const madness::coord_3d& r) {
  MADNESS_ASSERT(Axis < 3);
  auto cell = madness::FunctionDefaults<3>::get_cell();
  auto o_axis = (cell(Axis, 1) - cell(Axis, 0)) / 2.;
  return pow(r[Axis] - o_axis, K);
}

/// spheropole (spherical second moment) function w.r.t the center of the cell
static double r2_func3d(const madness::coord_3d& r) {
  auto cell = madness::FunctionDefaults<3>::get_cell();
  const auto x = r[0] - (cell(0, 1) - cell(0, 0)) / 2.;
  const auto y = r[1] - (cell(1, 1) - cell(1, 0)) / 2.;
  const auto z = r[2] - (cell(2, 1) - cell(2, 0)) / 2.;
  return x * x + y * y + z * z;
}

/// \brief Filter the moments of a function
/// \param[in,out] f on entry: a function to filter; on exit: function with Cartesian moments up to (and including) order \p k zeroed out
/// \param[in] k the maximum order of Cartesian moments to zero out
/// \param[in] fence if true, synchronize the function's world after filtering (default: true); since this is an in-place operation, only set `fence=false` when you know there are no pending readers of \p f
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

// This function test both the periodic and non-periodic versions of the Coulomb
// operator. In order to make this test valid set L to a high value so that
// charge distribution should not be able to see its neighbor.
int main(int argc, char**argv) {
  using namespace madness;
  auto &world = initialize(argc, argv);
  startup(world, argc, argv, true);

  int nerrors = 0;

  {

    // Function defaults
    int k = 7;
    double eps = std::pow(10., -k+2);
    FunctionDefaults<3>::set_k(k);
    Tensor<double> cell(3, 2);
    cell(0, 0) = -Lx / 2;
    cell(0, 1) = Ly / 2;
    cell(1, 0) = -Ly / 2;
    cell(1, 1) = Ly / 2;
    cell(2, 0) = -Lz / 2;
    cell(2, 1) = Lz / 2;
    FunctionDefaults<3>::set_cell(cell);
    FunctionDefaults<3>::set_thresh(eps);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_initial_level(8);
    FunctionDefaults<3>::set_truncate_mode(0);
    BoundaryConditions<3> bc_open(BC_FREE);
    BoundaryConditions<3> bc_periodic(BC_PERIODIC);
    FunctionDefaults<3>::set_bc(bc_periodic);

    // Create test charge density and the exact solution to Poisson's equation
    // with said charge density
    auto r000 =
        std::vector<coord_3d>{coord_3d(std::array<double, 3>{0.0, 0.0, 0.0})};
    auto r001 =
        std::vector<coord_3d>{coord_3d(std::array<double, 3>{0.0, 0.0, 1.0})};
    auto r00m1 =
        std::vector<coord_3d>{coord_3d(std::array<double, 3>{0.0, 0.0, -1.0})};
    Function<double, 3> unit =
        FunctionFactory<double, 3>(world).initial_level(4).f(unit_func3d);
    printf("projected unit\n");
    Function<double, 3> gdiffuse_000 = FunctionFactory<double, 3>(world)
                                           .initial_level(4)
                                           .special_points(r000)
                                           .f(gdiffuse_000_func3d);
    printf("projected gdiffuse_000\n");
    Function<double, 3> gtight_000 = FunctionFactory<double, 3>(world)
                                         .initial_level(4)
                                         .special_points(r000)
                                         .f(gtight_000_func3d);
    printf("projected gtight_000\n");
    Function<double, 3> gdiffuse_001 = FunctionFactory<double, 3>(world)
                                           .initial_level(4)
                                           .special_points(r001)
                                           .f(gdiffuse_001_func3d);
    printf("projected gdiffuse_001\n");
    Function<double, 3> gdiffuse_00m1 = FunctionFactory<double, 3>(world)
                                            .initial_level(4)
                                            .special_points(r00m1)
                                            .f(gdiffuse_00m1_func3d);
    printf("projected gdiffuse_00m1\n");
    Function<double, 3> gtight_00m1 = FunctionFactory<double, 3>(world)
                                          .initial_level(4)
                                          .special_points(r00m1)
                                          .f(gtight_00m1_func3d);
    printf("projected gtight_00m1\n");
    Function<double, 3> r2 =
        FunctionFactory<double, 3>(world).initial_level(4).f(r2_func3d);
    printf("projected r2\n");
    Function<double, 3> x =
        FunctionFactory<double, 3>(world).initial_level(4).f(
            kth_moment_func3d<0, 1>);
    Function<double, 3> y =
        FunctionFactory<double, 3>(world).initial_level(4).f(
            kth_moment_func3d<1, 1>);
    Function<double, 3> z =
        FunctionFactory<double, 3>(world).initial_level(4).f(
            kth_moment_func3d<2, 1>);
    printf("projected dipoles\n");
    unit.truncate();
    gdiffuse_000.truncate();
    gtight_000.truncate();
    gdiffuse_001.truncate();
    gdiffuse_00m1.truncate();
    gtight_00m1.truncate();
    r2.truncate();
    x.truncate();
    y.truncate();
    z.truncate();
    printf("truncated\n");

    printf("\n");

    print("gdiffuse = diffuse Gaussian at (0,0,-1)");
    auto gdiffuse = gdiffuse_00m1;
    print("gtight = tight Gaussian at (0,0,-1)");
    auto gtight = gtight_00m1;
    auto rho = gdiffuse - gtight;
    // filter_moments_inplace(rho,2);
    printf("rho = gdiffuse - gtight\n\n");

    // Average value of the test function
    printf("Average value of gdiffuse is %.8f\n", gdiffuse.trace());
    printf("Average value of gtight is %.8f\n", gtight.trace());
    printf("Average value of rho is %.8f\n", rho.trace());
    printf("Spheropole of rho is %.8f\n", rho.inner(r2));
    printf("Average value of unit density is %.8f\n", unit.trace());
    printf("Spheropole of unit density is %.8f\n", unit.inner(r2));

    printf("\n");

    // Create operators and apply
    SeparatedConvolution<double, 3> op =
        CoulombOperator(world, 1e-10, eps, bc_open.is_periodic());
    SeparatedConvolution<double, 3> pop =
        CoulombOperator(world, 1e-10, eps, bc_periodic.is_periodic());
    // N.B. non-periodic Coulomb with range restriction to [-L/2,L/2]
    SeparatedConvolution<double, 3> op_rr(
        world,
        madness::OperatorInfo(0.0, 1e-10, eps, madness::OT_G12,
                              /* truncate? */ false,
                              /* range restriction? */
                              bc_periodic.make_range_vector(1 /*, 1./L*/)
                              ),
        madness::no_lattice_sum<3>(), madness::FunctionDefaults<3>::get_k());
    // N.B. Coulomb with range restriction to [-L/2,L/2]
    SeparatedConvolution<double, 3> pop_rr(
        world,
        madness::OperatorInfo(0.0, 1e-10, eps, madness::OT_G12,
                              /* truncate? */ false,
                              /* range restriction? */
                              bc_periodic.make_range_vector(1 /*, 1./L*/)
                              ),
        madness::lattice_sum<3>(), madness::FunctionDefaults<3>::get_k());

    // print out norms vs displacement length
    //    std::cout << "RP operator norms\n";
    //    for (int n = 0; n <= 10; ++n) {
    //      for (int l = 0; l <= ((1 << n) + 5); ++l) {
    //        std::cout << "n=" << n << " l={" << l << ",0,0} ||op_{n,{l,0,0}}||="
    //                  << op_rr.norm(n, Key<3>(n, Vector<Translation, 3>({l, 0, 0})),
    //                                Key<3>(n, Vector<Translation, 3>({0,0,0})))
    //                  << "\n";
    //      }
    //    }
    //    std::cout << std::endl;

    const std::vector<std::reference_wrapper<const SeparatedConvolution<double, 3>>> test_operators = {op, pop, op_rr, pop_rr};
    const std::vector<std::string> test_operator_names = {"NP", "P", "RNP", "RP"};
    std::map<std::string, std::reference_wrapper<const SeparatedConvolution<double, 3>>> str2op;
    for(int i=0; i!=test_operators.size(); ++i) {
      str2op.emplace(test_operator_names[i], test_operators[i]);
    }
    const std::vector<Function<double, 3>> test_functions = {rho, gdiffuse, gtight};
    const std::vector<std::string> test_function_names = {"rho", "gdiffuse", "gtight"};  // must be this order since we assume 0 = 1 - 2

    std::vector<std::vector<Function<double, 3>>> V_op_f;
    int oi = 0;
    for(const auto& o: test_operators) {
      int fi = 0;

      V_op_f.push_back({});
      V_op_f.back().reserve(test_functions.size());

      for(const auto& f: test_functions) {

        printf("applying operator %s to %s ... ", test_operator_names[oi].c_str(), test_function_names[fi].c_str());
        const auto tstart = cpu_time();
        V_op_f[oi].emplace_back(apply(o.get(), f).truncate());
        const auto tstop = cpu_time();
        printf("%5.2e sec\n", tstop - tstart);

        printf("  <V_%s[%s]> = %.8f\n",
               test_operator_names[oi].c_str(), test_function_names[fi].c_str(), V_op_f[oi].back().trace());

        ++fi;
      }

      printf("\n");
      ++oi;
    }
    printf("\n");

    std::map<std::string, std::vector<Function<double, 3>>> str2V;
    for(int i=0; i!=test_operators.size(); ++i) {
      str2V.emplace(test_operator_names[i], V_op_f[i]);
    }

    for(int oi=0; oi != test_operators.size(); ++oi) {
      std::string ostr = test_operator_names[oi];
      std::string vstr = "V_" + ostr;
      const auto error =(str2V[ostr][0] - str2V[ostr][1] + str2V[ostr][2]).norm2();
      const auto tol = k <=10 ? 1e2 * eps : 5e-7;
      const bool success = error <= 1e2 * eps;
      if (!success) ++nerrors;
      std::cout << "||" << vstr << "(rho)-" << vstr << "(gdiffuse)+" << vstr
                << "(gtight)||="
                << error << (success ? " (PASS)" : " (FAIL)")
                << std::endl;
    }
    printf("\n");

    const int step = 1;
    const std::vector<std::string> axis_names = {"X", "Z"};
    for(const auto& axis_name: axis_names) {

      auto make_coord = [&](int r) {
        return axis_name=="X" ? coord_3d{double(r), 0, 0} : (axis_name=="Y" ? coord_3d{0, double(r), 0} : coord_3d{0, 0, double(r)});
      };

      for(int fi=0; fi != test_functions.size(); ++fi) {
        const auto& f = test_functions[fi];
        const auto L = axis_name=="X" ? Lx : (axis_name=="Y" ? Ly : Lz);

        printf("Scan along %s axis\n", axis_name.c_str());
        printf("%10c%18s%18s%18s%18s%18s%18s%18s%18s\n", std::tolower(axis_name[0]),
               test_function_names[fi].c_str(), "V_NP", "V_P", "V_RNP", "V_RP",
               "V_P-V_NP", "V_P-V_RP",
               "V_RP-V_NP");
        for (int r = -L / 2; r <= L / 2;
             r += step) {
          auto p = make_coord(r);
          printf("%10.2f%18.8f%18.8f%18.8f%18.8f%18.8f%18.8f%18.8f%18.8f\n", (double)r,
                 f(p), str2V["NP"][fi](p), str2V["P"][fi](p), str2V["RNP"][fi](p), str2V["RP"][fi](p),
                 str2V["P"][fi](p) - str2V["NP"][fi](p),
                 str2V["P"][fi](p) - str2V["RP"][fi](p),
                 str2V["RP"][fi](p) - str2V["NP"][fi](p)
                 );
        }

        // check periodicity of potential
        for (std::string opstr : {"P", "RP"}) {
          const auto ptol = k <= 10 ? 1e2 * eps : 1e-5;
          std::string result_str = "(PASS)";
          if (std::abs(str2V[opstr][fi](make_coord(L / 2)) -
                       str2V[opstr][fi](make_coord(-L / 2))) > ptol) {
            ++nerrors;
            result_str = "(FAIL)";
          }
          printf("check if V_%s(%s) is periodic along %s: %s\n", opstr.c_str(),
                 test_function_names[fi].c_str(),
                 axis_name.c_str(), result_str.c_str());
        }
      }
    }

  }

  madness::finalize();

  printf("nerrors = %d\n", nerrors);

  return nerrors;
}
