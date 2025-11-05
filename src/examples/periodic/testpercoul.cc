#include <madness/mra/mra.h>
#include <madness/mra/mw.h>
#include <madness/mra/operator.h>
#include <madness/chem/potentialmanager.h>

const int L = 18;
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
    int k = 10;
    double eps = std::pow(10., -k+2);
    FunctionDefaults<3>::set_k(k);
    Tensor<double> cell(3, 2);
    cell(0, 0) = -Lx / 2;
    cell(0, 1) = Ly / 2;
    cell(1, 0) = -Ly / 2;
    cell(1, 1) = Ly / 2;
    cell(2, 0) = -Lz / 2;
    cell(2, 1) = Lz / 2;
    BoundaryConditions<3> bc_open(BC_FREE);
    BoundaryConditions<3> bc_periodic(BC_PERIODIC);
    BoundaryConditions<3> bc_mixed({BC_FREE, BC_FREE, BC_FREE, BC_FREE, BC_PERIODIC, BC_PERIODIC});
    FunctionDefaults<3>::set_bc(bc_mixed);
    const auto bc = FunctionDefaults<3>::get_bc();
    Displacements<3>::reset_periodic_axes(bc.is_periodic());

    FunctionDefaults<3>::set_cell(cell);
    FunctionDefaults<3>::set_thresh(eps);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_initial_level(2);
    FunctionDefaults<3>::set_truncate_mode(0);

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
        CoulombOperator(world, 1e-10, eps, bc.is_periodic());

    auto range = bc.make_range<3>(1, .5/L);

    // N.B. non-periodic Coulomb with range restriction to [-L/2,L/2]
    SeparatedConvolution<double, 3> op_rr(
        world,
        madness::OperatorInfo(0.0, 1e-10, eps, madness::OT_G12,
                              /* truncate? */ false,
                              /* range restriction? */
                              range
                              ),
        madness::no_lattice_sum<3>());
    // N.B. Coulomb with range restriction to [-L/2,L/2]
    SeparatedConvolution<double, 3> pop_rr(
        world,
        madness::OperatorInfo(0.0, 1e-10, eps, madness::OT_G12,
                              /* truncate? */ false,
                              /* range restriction? */
                              range
                              ),
        bc.is_periodic());
    // N.B. Coulomb with range restriction to [-L/2,L/2]
    SeparatedConvolution<double, 3> pop2_rr(
        world,
        madness::OperatorInfo(0.0, 1e-10, eps, madness::OT_G12,
                              /* truncate? */ false,
                              /* range restriction? */
                              range
                              ),
        madness::no_lattice_sum<3>());
    pop2_rr.set_domain_periodicity(bc.is_periodic());

    // check operator norms
    {
      int nerrors_local = 0;

      constexpr auto ND = 3;
      for (int n = 1; n <= 5; ++n) {
        const auto twonm1 = 1 << (n - 1);
        const auto twonm3 = 1 << std::max((n - 1), 0);
        Key<ND> key(n, {twonm3, twonm1, 0});

        std::size_t disp_count = 0;
        std::array<std::optional<std::int64_t>, ND> box_radius;
        std::array<std::optional<std::int64_t>, ND> surface_thickness;
        auto &range = op_rr.get_range();
        for (int d = 0; d != ND; ++d) {
          if (range[d]) {
            box_radius[d] = range[d].N();
            surface_thickness[d] = range[d].finite_soft() ? 2 : 0;
          }
        }

        // surface displacements
        {
          BoxSurfaceDisplacementRange<ND> range_boundary_face_displacements(
              key, box_radius, surface_thickness, madness::no_lattice_sum<3>(),
              [](const auto level, const auto &dest,
                 const auto &displacement) -> bool { return true; });
          // check that all displacements are unique:
          {
            std::vector disps(range_boundary_face_displacements.begin(),
                              range_boundary_face_displacements.end());
            std::sort(disps.begin(), disps.end());
            auto it = std::unique(disps.begin(), disps.end());

            if (it != disps.end()) {
              std::cout << "Duplicates found!!" << std::endl;
              abort();
            }
          }

          auto process_displacement = [&](const Key<3>& disp) {

            auto rp2 = pop2_rr.get_ops();
            auto rp = pop_rr.get_ops();
            auto rnp = op_rr.get_ops();
            const auto twon = 1 << n;
            MADNESS_ASSERT(rp2.size() == rp.size());
            for(size_t mu=0; mu != rp.size(); ++mu) {
              for(int d=0; d!=3; ++d) {
                const auto& rp_R = rp[mu].getop(d)->nonstandard(n, disp[d])->R;
                const auto& rp2_R = rp2[mu].getop(d)->nonstandard(n, disp[d])->R;
                if (!rp_R.has_data()) {
                  MADNESS_ASSERT(!rp2_R.has_data());
                  continue;
                }

                if (bc.is_periodic()[d]) {
                  auto rp_mu_gau =
                      std::dynamic_pointer_cast<GaussianConvolution1D<double>>(
                          rp[mu].getop(d));

                  const auto &rnp_d0_R =
                      rnp[mu].getop(d)->nonstandard(n, disp[d])->R;
                  const auto &rnp_d1_R =
                      rnp[mu].getop(d)->nonstandard(n, disp[d] + twon)->R;
                  const auto &rnp_dm1_R =
                      rnp[mu].getop(d)->nonstandard(n, disp[d] - twon)->R;
                  MADNESS_ASSERT(rnp_d0_R.has_data() || rnp_d1_R.has_data() ||
                                 rnp_dm1_R.has_data());
                  Tensor<double> rp_R_recomputed;
                  if (rnp_d0_R.has_data())
                    rp_R_recomputed = copy(rnp_d0_R);
                  if (rnp_d1_R.has_data())
                    rp_R_recomputed = rp_R_recomputed.has_data()
                                          ? rp_R_recomputed + rnp_d1_R
                                          : copy(rnp_d1_R);
                  if (rnp_dm1_R.has_data())
                    rp_R_recomputed = rp_R_recomputed.has_data()
                                          ? rp_R_recomputed + rnp_dm1_R
                                          : copy(rnp_dm1_R);
                  const auto error_rp = (rp_R - rp_R_recomputed).normf();
                  const auto error_rp2 = (rp2_R - rnp_d0_R).normf();

                  if (error_rp > 1e-14 || error_rp2 > 1e-14) {
                    if (true) {
                      std::cout << "||RP||{n=" << n << ",l=0} -> {n=" << n
                                << ",l=" << disp.translation()
                                << "} = " << std::scientific
                                << pop_rr.norm(
                                       n, Key<3>(n, disp.translation()),
                                       Key<3>(n, Vector<Translation, 3>{0, 0, 0}))
                                << "\n";
                      std::cout << "||RP2||{n=" << n << ",l=0} -> {n=" << n
                                << ",l=" << disp.translation()
                                << "} = " << std::scientific
                                << pop2_rr.norm(
                                       n, Key<3>(n, disp.translation()),
                                       Key<3>(n, Vector<Translation, 3>{0, 0, 0}))
                                << "\n";
                    }
                    ++nerrors_local;
                  }
                }
              }
            }
          };

          for (int l = 0; l != (1 << (n - 1)) + 5; ++l) {
            process_displacement(Key<3>(n, Vector<Translation, 3>({l, 0, 0})));
            ++disp_count;
          }

          if (n <= 3) {
            for (auto &&disp : range_boundary_face_displacements) {
              process_displacement(disp);
              ++disp_count;
            }
          }

          process_displacement(Key<3>(n, Vector<Translation, 3>({-4,0,7})));
          process_displacement(Key<3>(n, Vector<Translation, 3>({-4,1,7})));

        } // box displacements

      }
      std::cout << "RNP operator norms check ... " << (nerrors_local > 0 ? "(FAIL)" : "(PASS)") << std::endl;
      nerrors += nerrors_local;
    }

    const std::vector<std::reference_wrapper<const SeparatedConvolution<double, 3>>> test_operators = {op, pop, op_rr, pop_rr, pop2_rr};
    const std::vector<std::string> test_operator_names = {"NP", "P", "RNP", "RP", "RP2"};
    std::map<std::string, std::reference_wrapper<const SeparatedConvolution<double, 3>>> str2op;
    for(size_t i=0; i!=test_operators.size(); ++i) {
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
    for(size_t i=0; i!=test_operators.size(); ++i) {
      str2V.emplace(test_operator_names[i], V_op_f[i]);
    }

    for(size_t oi=0; oi != test_operators.size(); ++oi) {
      std::string ostr = test_operator_names[oi];
      std::string vstr = "V_" + ostr;
      const auto error =(str2V[ostr][0] - str2V[ostr][1] + str2V[ostr][2]).norm2();
      const auto tol = k <=10 ? 2e2 * eps : 5e-7;
      const bool success = error <= tol;
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

      for(size_t fi=0; fi != test_functions.size(); ++fi) {
        const auto& f = test_functions[fi];
        const auto L = axis_name=="X" ? Lx : (axis_name=="Y" ? Ly : Lz);

        printf("Scan along %s axis\n", axis_name.c_str());
        printf("%10c%18s%18s%18s%18s%18s%18s%18s%18s%18s%18s\n", std::tolower(axis_name[0]),
               test_function_names[fi].c_str(), "V_NP", "V_P", "V_RNP", "V_RP", "V_RP2",
               "V_P-V_NP", "V_P-V_RP",
               "V_RP-V_NP", "V_RP-V_RP2");
        for (int r = -L / 2; r <= L / 2;
             r += step) {
          auto p = make_coord(r);
          printf("%10.2f%18.8f%18.8f%18.8f%18.8f%18.8f%18.8f%18.8f%18.8f%18.8f%18.8f\n", (double)r,
                 f(p), str2V["NP"][fi](p), str2V["P"][fi](p), str2V["RNP"][fi](p), str2V["RP"][fi](p), str2V["RP2"][fi](p),
                 str2V["P"][fi](p) - str2V["NP"][fi](p),
                 str2V["P"][fi](p) - str2V["RP"][fi](p),
                 str2V["RP"][fi](p) - str2V["NP"][fi](p),
                 str2V["RP"][fi](p) - str2V["RP2"][fi](p)
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
    }  // gaussian density test

    // test superfine structure
    if (false) {  // tiny MRA box
      const auto n=15;
      const auto twon = (1<<n);
      const auto one_over_twon = 1./twon;
      const Key<3> key(15, {(1<<14) - 1, (1<<14)-1, 18932});
      const auto ff =
          ScalingFunctionFunctor<3>({Lx/2, Ly/2, Lz/2}, key, {0,0,0});
      std::cout << "ff(ff.special_points()[0]) = " << ff(ff.special_points()[0]) << std::endl;
      madness::real_function_3d f =
          madness::real_factory_3d(world)
              .functor(ff)
              .truncate_mode(0)
          //.special_points(std::vector{coord_3d{-0.5*one_over_twon,-0.5*one_over_twon,-L + (18932*2*L+0.5) * one_over_twon}})
          ;
      f.truncate();
      f.print_tree(std::cout);

      auto V_np = op(f);
      auto V_rnp = op_rr(f);

      std::string axis_name = "Z";
      printf("Scan along %s axis\n", axis_name.c_str());
      printf("%10c%18s%18s%18s\n", std::tolower(axis_name[0]),
             "ρn", "V_NP", "V_RNP");
      const auto step = 0.18;
      for (double r = -L / 2; r <= L / 2;
           r += step) {
        coord_3d p = {-0.5*one_over_twon, -0.5*one_over_twon, (double)r};
        printf("%10.2f%18.8f%18.8f%18.8f\n", (double)r,
               f(p), V_np(p), V_rnp(p)
        );
      }

    }  // superfine box test

    if (true) { // sum of point charges
      int nerrors_local = 0;
      auto check = [&](double error, double tol) {
        bool success = error < tol;
        if (!success)
          ++nerrors_local;
      };
      auto scheck = [&](const auto& str, double error, double tol) {
        const bool success = error < tol;
        if (!success)
          ++nerrors_local;
        std::cout << str << " = " << error
                  << (success ? "(PASS)" : "(FAIL)") << std::endl;
      };

      constexpr auto x0 = -1.0;
      const auto xs = {-9., -7.9, -4.7, -4.3, x0};

      std::vector<madness::Atom> mad_atoms;
      // uncomment all charges to stress test
      for (const auto& atom : std::vector<std::array<double, 3>>{
//               {x0, 0., 1.4},
//               {x0, 0., -0.4},
//               {x0, 0., 3.2},
//               {x0, 0., -2.2},
//               {x0, 0., 5.0},
//               {x0, 0., -4.0},
//               {x0, 0., 6.8},
//               {x0, 0., -5.8},
//               {x0, 0., 8.6},
               {x0, 0., -7.6}
           }) {
        mad_atoms.emplace_back(atom[0], atom[1], atom[2], +1., 1);
      }
      madness::Molecule mol(std::move(mad_atoms),
                            eps/10.);
      auto ρn_func = std::make_shared<madness::NuclearDensityFunctor>(mol);
      using fi3d = madness::FunctionFunctorInterface<double, 3>;
      madness::real_function_3d ρnuc =
          madness::real_factory_3d(world)
              .functor(std::static_pointer_cast<fi3d>(ρn_func))
              .truncate_mode(0)
              .truncate_on_project();
      ρnuc.truncate();

      auto V_np = op(ρnuc);
      auto V_np_exact_func = std::make_shared<madness::MolecularPotentialFunctor>(mol);
      madness::real_function_3d  V_np_exact =
          madness::real_factory_3d(world)
              .functor(std::static_pointer_cast<fi3d>(V_np_exact_func))
              .truncate_mode(0)
              .truncate_on_project();
      V_np_exact *= -1;
      scheck("||NP - NP(exact)||", (V_np - V_np_exact).norm2(), 1e2 * eps);

      {
        for (auto x: xs) {
          printf("Scan along {%f,0,z} ray\n", x);
          printf("%10s%18s%18s%18s%18s%18s\n", "z", "V_NP",
                 "V_NP(exact)", "δV_NP", "V_NP(exact,mw)",
                 "V_NP(exact,δmw)");
          const auto step = 0.18;
          for (double z = -L / 2; z <= L / 2; z += step) {
            coord_3d p = {x, 0, (double)z};
            const auto v = V_np(p);
            const auto v_ex = -(*V_np_exact_func)(p);
            const auto v_ex_mw = V_np_exact(p);

            printf("%10.2f%18.8f%18.8f%18.8f%18.8f%18.8f\n", (double)z, v,
                   v_ex, v-v_ex, v_ex, v_ex_mw - v_ex);

            // near the point charges the errors due to projection can be large, but operator application reduces them
            check(v-v_ex_mw, (x == x0 ? 5e2 : 5) * eps);
            check(v-v_ex, 5 * eps);
          }
        }
      }

      auto V_rp = pop_rr(ρnuc);
      auto V_rp2 = pop2_rr(ρnuc);

      scheck("||RP - RP2||", (V_rp - V_rp2).norm2(), 1e2 * eps);

      auto V_rp_exact_func = std::make_shared<madness::WignerSeitzPotentialFunctor>(mol, FunctionDefaults<3>::get_cell(), FunctionDefaults<3>::get_bc(), range);
      madness::real_function_3d  V_rp_exact =
          madness::real_factory_3d(world)
              .functor(std::static_pointer_cast<fi3d>(V_rp_exact_func))
              .truncate_mode(0)
              .truncate_on_project();
      V_rp_exact *= -1;
      scheck("||RP - RP(exact)||", (V_rp - V_rp_exact).norm2(), 1e2 * eps);
      scheck("||RP2 - RP(exact)||", (V_rp2 - V_rp_exact).norm2(), 1e2 * eps);

      {
        for (auto x: xs) {
          printf("Scan along {%f,0,z} ray\n", x);
          printf("%10s%18s%18s%18s%18s%18s\n", "z", "V_RP",
                 "V_RP(exact)", "δV_RP", "V_RP(exact,mw)",
                 "V_RP(exact,δmw)");
          const auto step = 0.18;
          for (double z = -L / 2; z <= L / 2; z += step) {
            coord_3d p = {x, 0, (double)z};
            const auto v = V_rp(p);
            const auto v_ex_mw = V_rp_exact(p);
            const auto v_ex = -(*V_rp_exact_func)(p);

            printf("%10.2f%18.8f%18.8f%18.8f%18.8f%18.8f\n", (double)z, v,
                   v_ex, v - v_ex, v_ex_mw, v_ex_mw - v_ex);

            // near the point charges the errors due to projection can be large, but operator application reduces them
            check(v-v_ex_mw, (x == x0 ? 5e2 : 5) * eps);
            check(v-v_ex, 5 * eps);
          }
        }
      }

      auto V_rnp = op_rr(ρnuc);
      auto V_rnp_exact_func = std::make_shared<madness::WignerSeitzPotentialFunctor>(mol, FunctionDefaults<3>::get_cell(), FunctionDefaults<3>::get_bc(), range, std::array{0,0,0});
      madness::real_function_3d  V_rnp_exact =
          madness::real_factory_3d(world)
              .functor(std::static_pointer_cast<fi3d>(V_rnp_exact_func))
              .truncate_mode(1)
              .truncate_on_project();
      V_rnp_exact *= -1;
      scheck("||RNP - RNP(exact)||", (V_rnp - V_rnp_exact).norm2(), 1e2 * eps);

      {
        for (auto x: xs) {
          printf("Scan along {%f,0,z} ray\n", x);
          printf("%10s%18s%18s%18s%18s%18s\n", "z", "V_RNP",
                 "V_RNP(exact)", "δV_RNP", "V_RNP(exact, mw)",
                 "V_RNP(exact,δmw)");
          const auto step = 0.18;
          for (double z = -L / 2; z <= L / 2; z += step) {
            coord_3d p = {x, 0, (double)z};
            const auto v = V_rnp(p);
            const auto v_ex_mw = V_rnp_exact(p);
            const auto v_ex = -(*V_rnp_exact_func)(p);

            printf("%10.2f%18.8f%18.8f%18.8f%18.8f%18.8f\n", (double)z, v,
                   v_ex, v - v_ex, v_ex_mw, v_ex_mw - v_ex);

            // near the point charges the errors due to projection can be large, but operator application reduces them
            check(v-v_ex_mw, (x == x0 ? 5e2 : 5) * eps);
            check(v-v_ex, 5 * eps);
          }
        }
      }

      plot_plane(world, std::vector{V_np - V_np_exact,
                                    V_rnp - V_rnp_exact, V_rp2 - V_rp_exact, V_rp - V_rp2}, "input");

      auto V_p = pop(ρnuc);
      {
        std::string axis_name = "Z";
        printf("Scan along %s axis\n", axis_name.c_str());
        printf("%10c%18s%18s%18s%18s%18s%18s\n", std::tolower(axis_name[0]),
               "ρn", "V_NP", "V_RNP", "V_P", "V_RP", "V_RP2");
        const auto step = 0.18;
        for (double r = -L / 2; r <= L / 2; r += step) {
          coord_3d p = {x0, 0, (double)r};
          printf("%10.2f%18.8f%18.8f%18.8f%18.8f%18.8f%18.8f\n", (double)r,
                 ρnuc(p), V_np(p), V_rnp(p), V_p(p), V_rp(p), V_rp2(p));
        }
      }

      // plot contributions from individual Gaussians in the separated representation
      // do this for every 10th term in the separated representation
      if (false) {
        const auto rank = op.get_rank();
        MADNESS_ASSERT(rank == op_rr.get_rank());
        for (int mu = 0; mu < rank; mu += 10) {
          std::cout << "=== op[mu].alpha = "
                    << std::static_pointer_cast<GaussianConvolution1D<double>>(
                           op.get_ops().at(mu).getop(0))
                           ->expnt
                    << std::endl;
          SeparatedConvolution<double, 3> op_mu(
              world, std::vector{op.get_ops().at(mu)});
          SeparatedConvolution<double, 3> op_rr_mu(
              world, std::vector{op_rr.get_ops().at(mu)});

          auto V_np = op_mu(ρnuc);
          auto V_rnp = op_rr_mu(ρnuc);

          std::string axis_name = "Z";
          printf("Scan along %s axis\n", axis_name.c_str());
          printf("%10c%18s%18s%18s\n", std::tolower(axis_name[0]), "ρn", "V_NP",
                 "V_RNP");
          const auto step = 0.18;
          for (double r = -L / 2; r <= L / 2; r += step) {
            coord_3d p = {0, 0, (double)r};
            printf("%10.2f%18.8f%18.8f%18.8f\n", (double)r, ρnuc(p), V_np(p),
                   V_rnp(p));
          }
        }
      }

      nerrors += nerrors_local;
    }  // point charge test

  }

  madness::finalize();

  printf("nerrors = %d\n", nerrors);

  return nerrors;
}
