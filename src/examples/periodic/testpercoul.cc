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

  {

    // Function defaults
    int k = 10;
    double eps = 1e-8;
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
                                            .initial_level(6).notruncate_on_project()
                                            .special_points(r00m1)
                                            .f(gdiffuse_00m1_func3d);
    printf("projected gdiffuse_00m1\n");
    Function<double, 3> gtight_00m1 = FunctionFactory<double, 3>(world)
                                          .initial_level(6).notruncate_on_project()
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
//    gdiffuse_00m1.truncate();
//    gtight_00m1.truncate();
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
        madness::OperatorInfo(0.0, 1e-10, eps,
                              madness::OT_G12, /* truncate? */ false, /* range restriction? */ bc_periodic.make_range_vector(1/*, 1./L*/)),
        madness::no_lattice_sum<3>(),
        madness::FunctionDefaults<3>::get_k());
    // N.B. Coulomb with range restriction to [-L/2,L/2]
    SeparatedConvolution<double, 3> pop_rr(
                world,
                madness::OperatorInfo(0.0, 1e-10, eps,
                                      madness::OT_G12, /* truncate? */ false, /* range restriction? */ bc_periodic.make_range_vector(1/*, 1./L*/)),
                madness::lattice_sum<3>(),
                madness::FunctionDefaults<3>::get_k());

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

    printf("applying non-periodic operator to rho ...\n");
    Function<double, 3> V_nonperiodic = apply(op, rho);
    V_nonperiodic.truncate();
    printf("applying non-periodic operator to gdiffuse ...\n");
    Function<double, 3> V_nonperiodic_d = apply(op, gdiffuse);
    V_nonperiodic_d.truncate();
    printf("applying non-periodic operator to gtight ...\n");
    Function<double, 3> V_nonperiodic_t = apply(op, gtight);
    V_nonperiodic_t.truncate();

    printf("applying periodic operator to rho ...\n");
    Function<double, 3> V_periodic = apply(pop, rho);
    V_periodic.truncate();
    printf("applying periodic operator to gdiffuse ...\n");
    Function<double, 3> V_periodic_d = apply(pop, gdiffuse);
    V_periodic_d.truncate();
    printf("applying periodic operator to gtight ...\n");
    Function<double, 3> V_periodic_t = apply(pop, gtight);
    V_periodic_t.truncate();

    printf("applying range-restricted periodic operator to rho ...\n");
    Function<double, 3> V_rrperiodic = apply(pop_rr, rho);
    V_rrperiodic.truncate();
    printf("applying range-restricted periodic operator to gdiffuse ...\n");
    Function<double, 3> V_rrperiodic_d = apply(pop_rr, gdiffuse);
    V_rrperiodic_d.truncate();
    printf("applying range-restricted periodic operator to gtight ...\n");
    Function<double, 3> V_rrperiodic_t = apply(pop_rr, gtight);
    V_rrperiodic_t.truncate();

    printf("\n");

    printf("Average value of V_nonperiodic(rho) is %.8f\n",
           V_nonperiodic.trace());
    printf("Average value of V_periodic(rho) is %.8f\n", V_periodic.trace());
    printf("Average value of V_periodic(gdiffuse) is %.8f\n",
           V_periodic_d.trace());
    printf("Average value of V_periodic(gtight) is %.8f\n",
           V_periodic_t.trace());
    printf("Average value of V_nonperiodic(gdiffuse) is %.8f\n",
           V_nonperiodic_d.trace());
    printf("Average value of V_nonperiodic(gtight) is %.8f\n",
           V_nonperiodic_t.trace());
    printf("1st x moment of V_nonperiodic(rho) is %.8f\n",
           V_nonperiodic.inner(x));
    printf("1st y moment of V_nonperiodic(rho) is %.8f\n",
           V_nonperiodic.inner(y));
    printf("1st z moment of V_nonperiodic(rho) is %.8f\n",
           V_nonperiodic.inner(z));
    printf("1st x moment of V_periodic(rho) is %.8f\n", V_periodic.inner(x));
    printf("1st y moment of V_periodic(rho) is %.8f\n", V_periodic.inner(y));
    printf("1st z moment of V_periodic(rho) is %.8f\n", V_periodic.inner(z));

    auto V_periodic_regauged =
        V_periodic + unit * (V_nonperiodic.trace() - V_periodic.trace());

    int step = 1;
    printf("Scan along X axis\n");
    printf("%10s%18s%18s%18s%18s%18s%18s%18s\n", "x", "rho",
           "V_np[gdiffuse]", "V_np[gtight]",
           "V_p[gdiffuse]", "V_p[gtight]",
           "V_rp[gdiffuse]", "V_rp[gtight]");
    for (int rx = std::max(-100, -Lx / 2); rx <= std::min(100, Lx / 2);
         rx += step) {
      coord_3d p(std::array<double, 3>{double(rx), 0, 0});
      printf("%10.2f%18.8f%18.8f%18.8f%18.8f%18.8f%18.8f%18.8f\n", p[0],
             rho(p), V_periodic_d(p), V_periodic_t(p),
             V_nonperiodic_d(p), V_nonperiodic_t(p),
             V_rrperiodic_d(p), V_rrperiodic_t(p));
    }
    printf("Scan along X axis\n");
    printf("%10s%18s%18s%18s%18s%18s%18s%18s\n", "x", "rho",
           "V_np[rho]",
           "V_p[rho]", "V_rp[rho]",
           "V_p[rho]-V_np[rho]",
           "V_p[rho]-V_rp[rho]",
           "V_rp[rho]-V_np[rho]");
    for (int rx = std::max(-100, -Lx / 2); rx <= std::min(100, Lx / 2);
         rx += step) {
      coord_3d p(std::array<double, 3>{double(rx), 0, 0});
      printf("%10.2f%18.8f%18.8f%18.8f%18.8f%18.8f%18.8f%18.8f\n", p[0],
             rho(p), V_nonperiodic(p), V_periodic(p), V_rrperiodic(p),
             V_periodic(p) - V_nonperiodic(p), V_periodic(p) - V_rrperiodic(p), V_rrperiodic(p) - V_nonperiodic(p));
    }
    printf("Scan along X axis\n");
    printf("%10s%36s%36s%36s\n", "x", "V_np([rho]-[gdiffuse]+[gtight])",
           "V_p([rho]-[gdiffuse]+[gtight])", "V_rp([rho]-[gdiffuse]+[gtight])");
    for (int rx = std::max(-100, -Lx / 2); rx <= std::min(100, Lx / 2);
         rx += step) {
      coord_3d p(std::array<double, 3>{double(rx), 0, 0});
      printf("%10.2f%36.8f%36.8f%36.8f\n", p[0],
             V_nonperiodic(p)-V_nonperiodic_d(p)+V_nonperiodic_t(p),
             V_periodic(p)-V_periodic_d(p)+V_periodic_t(p),
             V_rrperiodic(p)-V_rrperiodic_d(p)+V_rrperiodic_t(p));
    }
    printf("Scan along Z axis\n");
    printf("%10s%18s%18s%18s%18s%18s%18s%18s\n", "x", "rho",
           "V_np[gdiffuse]", "V_np[gtight]",
           "V_p[gdiffuse]", "V_p[gtight]",
           "V_rp[gdiffuse]", "V_rp[gtight]");
    for (int rz = std::max(-100, -Lz / 2); rz <= std::min(100, Lz / 2);
         rz += step) {
      coord_3d p(std::array<double, 3>{0, 0, double(rz)});
      printf("%10.2f%18.8f%18.8f%18.8f%18.8f%18.8f%18.8f%18.8f\n", p[2],
             rho(p), V_periodic_d(p), V_periodic_t(p),
             V_nonperiodic_d(p), V_nonperiodic_t(p),
             V_rrperiodic_d(p), V_rrperiodic_t(p));
    }
    printf("Scan along Z axis\n");
    printf("%10s%18s%18s%18s%18s%18s%18s%18s\n", "z", "rho",
           "V_np[rho]",
           "V_p[rho]", "V_rp[rho]",
           "V_p[rho]-V_np[rho]",
           "V_p[rho]-V_rp[rho]",
           "V_rp[rho]-V_np[rho]");
    for (int rz = std::max(-100, -Lz / 2); rz <= std::min(100, Lz / 2);
         rz += step) {
      coord_3d p(std::array<double, 3>{0, 0, double(rz)});
      printf("%10.2f%18.8f%18.8f%18.8f%18.8f%18.8f%18.8f%18.8f\n", p[2],
             rho(p), V_nonperiodic(p), V_periodic(p), V_rrperiodic(p),
             V_periodic(p) - V_nonperiodic(p), V_periodic(p) - V_rrperiodic(p), V_rrperiodic(p) - V_nonperiodic(p));
    }
    printf("Scan along Z axis\n");
    printf("%10s%36s%36s%36s\n", "x", "V_np([rho]-[gdiffuse]+[gtight])",
           "V_p([rho]-[gdiffuse]+[gtight])", "V_rp([rho]-[gdiffuse]+[gtight])");
    for (int rz = std::max(-100, -Lz / 2); rz <= std::min(100, Lz / 2);
         rz += step) {
      coord_3d p(std::array<double, 3>{0, 0, double(rz)});
      printf("%10.2f%36.8f%36.8f%36.8f\n", p[2],
             V_nonperiodic(p)-V_nonperiodic_d(p)+V_nonperiodic_t(p),
             V_periodic(p)-V_periodic_d(p)+V_periodic_t(p),
             V_rrperiodic(p)-V_rrperiodic_d(p)+V_rrperiodic_t(p));
    }

    std::cout << "||V_np(rho)-V_np(gdiffuse)+V_np(gtight)||=" << (V_nonperiodic-V_nonperiodic_d+V_nonperiodic_t).trace() << std::endl;
    std::cout << "||V_p(rho)-V_p(gdiffuse)+V_p(gtight)||=" << (V_periodic-V_periodic_d+V_periodic_t).trace() << std::endl;
    std::cout << "||V_rp(rho)-V_rp(gdiffuse)+V_rp(gtight)||=" << (V_rrperiodic-V_rrperiodic_d+V_rrperiodic_t).trace() << std::endl;

    if (false) {
      const auto rank = pop_rr.get_rank();
      for (int mu = 0; mu != rank; ++mu) {
        SeparatedConvolution<double, 3> pop_mu(world, std::vector<ConvolutionND<double, 3>>{1, pop.ops_mu(mu)}, madness::lattice_sum<3>(),
                                                  madness::FunctionDefaults<3>::get_k());
        SeparatedConvolution<double, 3> pop_rr_mu(world, std::vector<ConvolutionND<double, 3>>{1, pop_rr.ops_mu(mu)}, madness::lattice_sum<3>(),
                                              madness::FunctionDefaults<3>::get_k());

        std::cout << "expr="
                  << std::dynamic_pointer_cast<
                         madness::GaussianConvolution1D<double>>(pop_rr.ops_mu(mu).getop(0))
                         ->expnt << std::endl
                  << "||V_p(rho)-V_p(gdiffuse)+V_p(gtight)||="
                  << (apply(pop_mu, rho) - apply(pop_mu, gdiffuse_00m1) +
                      apply(pop_mu, gtight_00m1))
                         .trace()
                  << std::endl
                  << "||V_rp(rho)-V_rp(gdiffuse)+V_rp(gtight)||="
                  << (apply(pop_rr_mu, rho) - apply(pop_rr_mu, gdiffuse_00m1) +
                      apply(pop_rr_mu, gtight_00m1))
                         .trace()
                  << std::endl
            ;
      }
    }
  }

  madness::finalize();
}
