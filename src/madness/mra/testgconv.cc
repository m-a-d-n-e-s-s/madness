/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


  $Id$
*/

/// \file mra/testgconv.cc
/// \brief Test convolution with Gaussian * polyn

#include <madness/mra/mra.h>
#include <madness/mra/mw.h>
#include <madness/mra/operator.h>
#include <madness/constants.h>

using namespace madness;

// WARNING even legacy tests prior to range-restricted op additions would not pass with k=8! Need k=10.
static const int k = 10;
static const double thresh = 1e-6;
static constexpr double L = 50;
static constexpr double Lx = L;
static constexpr double Ly = L;
static constexpr double Lz = L;
const auto npts = 10001;

constexpr std::size_t NDIM = 1;
static_assert(NDIM >= 1 && NDIM <= 3);
using coord_t = Vector<double,NDIM>;
using dim_t = Vector<double,NDIM>;
constexpr coord_t one(1.0);  // {1,1,1...}
constexpr auto make_rx = [](double x) -> coord_t {
  static_assert(NDIM >= 1 && NDIM <= 3);
  if constexpr (NDIM==1)
    return {x};
  else if constexpr (NDIM==2)
    return {x, 0.0};
  else
    return {x, 0.0, 0.0};
};
constexpr coord_t onex = make_rx(1.);  // {1,1,1...}
using real_factory_t = FunctionFactory<double, NDIM>;
using real_function_t = Function<double, NDIM>;
using real_convolution_t = SeparatedConvolution<double, NDIM>;

/// \brief Return the size of the simulation cell in user coordinates
constexpr dim_t cell_extent() {
  static_assert(NDIM >= 1 && NDIM <= 3);
  if constexpr (NDIM == 1)
    return {2*Lx};
  else if constexpr (NDIM == 2) {
    return {2*Lx, 2*Ly};
  }
  else { // NDIM == 3
    return {2*Lx, 2*Ly, 2*Lz};
  }
}

double cell_volume() {
  const auto size = cell_extent();
  return std::accumulate(
      size.begin(), size.end(), 1.,
      [](const auto &acc, const auto &v) { return acc * v; });
}

// normalized identity function = 1/volume
double unit(const coord_t& r) {
  const double fac = 1/cell_volume();
  return fac;
}

// normalized identity function in the right half of the box = 1/L along first dimension
constexpr double rhunit(const coord_t& r) {
  return r[0] < 0 ? 0 : 1/(2 * cell_volume());
}

// 0th scaling function = in box {13,4177}
constexpr double s_13_4177_0(const coord_t& r) {
  const auto n = 13;
  const auto ootwonm1 = 1. / (1 << (n - 1));
  const auto δ = ootwonm1 * L;
  const auto l = 4177;
  const auto xmin = -L + l * δ;
  const auto xmax = xmin + δ;
  return r[0] <= xmin || r[0] > xmax ? 0 : pow(2.0,n/2.)/sqrt(2.0*L);
}

// exp(-|r|^2) / (sqrt(pi))^d = normalized unit gaussian at the origin
double g(const coord_t& r) {
    static const double fac = std::pow(constants::inv_sqrt_pi, coord_t::static_size);
    return exp(-inner(r,r)) * fac;
}

// exp(-|r-1|^2) / sqrt(pi) = normalized unit gaussian at 1
double g1(const coord_t& r) {
  static const double fac = std::pow(constants::inv_sqrt_pi, coord_t::static_size);
  return exp(-inner(r-onex,r-onex)) * fac;
}

// exp(-(|r|/(2L))^2) / (2L sqrt(pi)) = diffuse gaussian at origin, normalized over infinity
double gd(const coord_t& r) {
  static const double fac = std::pow(constants::inv_sqrt_pi, coord_t::static_size) / cell_volume();
  const auto rt = r/cell_extent();
  return exp(-inner(rt,rt)) * fac;
}

// exp(-(|r-1|/L)^2) / (L sqrt(pi)) = diffuse gaussian at 1, normalized over infinity
double gd1(const coord_t& r) {
  static const double fac = std::pow(constants::inv_sqrt_pi, coord_t::static_size) / cell_volume();
  const auto rt = (r-onex)/cell_extent();
  return exp(-inner(rt,rt)) * fac;
}

// exp(-alpha * |r|^2) * sqrt(alpha/pi)^d = tight gaussian at origin, normalized over infinity
double gt(const coord_t& r) {
  const auto alpha = 1e4;
  static const double fac = std::pow(sqrt(alpha) * constants::inv_sqrt_pi, coord_t::static_size);
  return exp(-alpha * inner(r,r)) * fac;
}

// exp(-alpha * |r-1|^2) * sqrt(alpha/pi)^d = tight gaussian at 1, normalized over infinity
double gt1(const coord_t& r) {
  const auto alpha = 1e4;
  static const double fac = std::pow(sqrt(alpha) * constants::inv_sqrt_pi, coord_t::static_size);
  const auto rt = r-onex;
  return exp(-alpha * inner(rt,rt)) * fac;
}

// g' = exp(-r^2) / sqrt(pi) = derivative of a normalized gaussian
double gprime(const coord_1d& r) {
    static const double fac = 1.0/sqrt(constants::pi);
    return -2.* r[0] * exp(-r[0]*r[0]) * fac;
}

// conv g()
double convg(const coord_t& r) {
    static const double fac = pow(1.0/sqrt(2.0*constants::pi), NDIM);
    return exp(-0.5*inner(r,r)) * fac;
}

// conv g1()
double convg1(const coord_1d& r) {
  static const double fac = 1.0/sqrt(2.0*constants::pi);
  return exp(-(r[0]-1)*(r[0]-1)*0.5) * fac;
}

// sqrt(8)*x*exp(-x^2)
double h(const coord_1d& r) {
    static const double fac = sqrt(8.0);
    return exp(-r[0]*r[0]) * fac * r[0];
}

// g() conv h() == h() conv g()
double gconvh(const coord_1d& r) {
    return exp(-0.5*r[0]*r[0]) * r[0];
}

// D(conv) (g())
double conv_prime_g(const coord_1d& r) {
    static const double fac = 1.0/sqrt(2.0*constants::pi);
    return -exp(-0.5*r[0]*r[0]) * r[0] *fac;
}


// 1d Gaussian on [-L,L] centered at R, with exponent a, normalized over real axis: sqrt(a/pi) exp(-a (x-R)^2)
struct G : public FunctionFunctorInterface<double, 1> {
  double a, R;
  G(double a, double R) : a(a), R(R) {
  }
  double operator()(const coord_1d& r) const {
    const auto x = r[0];
    const auto result = sqrt(a / constants::pi) * exp(-a * pow(R - x, 2));
    return result;
  }
};

// convolution of 1-d Gaussian kernel sqrt(a/pi) exp(-a (x-y)^2) with 1d Gaussian on [-L,L] centered at R sqrt(b/pi) exp(-b (x-R)^2)
// in Wolfram: Integrate[ Exp[-a (x - y)^2] Sqrt[a/Pi] Exp[-b (y - R)^2]  Sqrt[b/Pi], {y, -L, L}]
struct GConvGNP : public FunctionFunctorInterface<double, 1> {
    double a, b, R, L;
    GConvGNP(double a, double b, double R, double L) : a(a), b(b), R(R), L(L) {
    }
    double operator()(const coord_1d& r) const {
      const auto x = r[0];
      const auto result =
          (sqrt(a) * sqrt(b) * exp(-a * b * pow(R - x, 2) / (a + b)) *
           (erf(((a + b) * L - b * R - a * x) / sqrt(a + b)) -
            erf((-(b * (L + R)) - a * (L + x)) / sqrt(a + b)))) /
          (2. * sqrt(a + b) * sqrt(constants::pi));
      return result;
    }
};

// convolution of 1-d Gaussian kernel sqrt(a/pi) exp(-a (x-y)^2) restricted to |x-y|<=L
// with 1d Gaussian on [-L,L] centered at R sqrt(b/pi) exp(-b (x-R)^2)
// in Wolfram: Integrate[ Exp[-a (x - y)^2] Sqrt[a/Pi] Exp[-b (y - R)^2]  Sqrt[b/Pi], {y, Max[x - L, -L], Max[Min[x + L, L], Max[x - L, -L]]}]
struct GConvGRNP : public FunctionFunctorInterface<double, 1> {
  double a, b, R, L;
  GConvGRNP(double a, double b, double R, double L) : a(a), b(b), R(R), L(L) {
  }
  double operator()(const coord_1d& r) const {
    const auto x = r[0];
    using std::max;
    using std::min;
    const auto result =
        (sqrt(a) * sqrt(b) * exp(-a * b * pow(R - x, 2) / (a + b)) *
         (erf((-(b*R) - a*x + (a + b)* max(max(-L,-L + x), min(L,L + x)))/sqrt(a + b))
          -erf((-(b*R) - a*x + (a + b)*max(-L,-L + x))/sqrt(a + b)))) /
        (2. * sqrt(a + b) * sqrt(constants::pi));
    return result;
  }
};

// periodic version of GConvGRNP
struct GConvGRP : public FunctionFunctorInterface<double, 1> {
  GConvGRNP rnp;
  coord_1d period;
  GConvGRP(double a, double b, double R, double L)
      : rnp(a, b, R, L), period{2 * L} {}
  double operator()(const coord_1d &r) const {
    return rnp(r) + rnp(r - period) + rnp(r + period);
  }
};

int test_gconv(World& world) {
    constexpr coord_t origin(0.0), lo = make_rx(-L), hi = make_rx(L);
    double width = 2.0*L;
    int success=0;

    if (world.rank() == 0) print("Test gconv operation");

    real_function_t f = real_factory_t(world).f(g);
    double error=f.trace()-1.0;
    print("error in integral(g) ", error);
    if (abs(error)>FunctionDefaults<NDIM>::get_thresh()) success++;
    print("success 0 ", success);

    real_function_t f1 = real_factory_t(world).f(g1);
    error=f1.trace()-1.0;
    print("error in integral(g1) ", error);
    if (abs(error)>FunctionDefaults<NDIM>::get_thresh()) success++;
    print("success 0 ", success);

    // convolve with a normalized Gaussian kernel
    std::vector< std::shared_ptr< Convolution1D<double> > > ops(1);
    ops[0].reset(new GaussianConvolution1D<double>(k, width/sqrt(constants::pi),
            width*width, 0, false));
    real_convolution_t op(world, ops);

    real_function_t opf = op(f);
    error=opf.trace()-1.0;
    print("error in integral(op(g)) ", error);
    if (abs(error)>FunctionDefaults<NDIM>::get_thresh()) success++;
    print("success 1 ", success);

    real_function_t opf1 = op(f1);
    error=opf1.trace()-1.0;
    print("error in integral(op(g1)) ", error);
    if (abs(error)>FunctionDefaults<NDIM>::get_thresh()) success++;
    print("success 1b ", success);

    real_function_t exact = real_factory_t(world).f(convg);
    print("norm2(g conv g - exact)", (opf - exact).norm2());
    error = (opf - exact).norm2();
    if (abs(error) > FunctionDefaults<NDIM>::get_thresh())
      success++;
    print("success 2 ", success);

    if constexpr (NDIM == 1) {

      real_function_1d exact1 = real_factory_1d(world).f(convg1);
      print("norm2(g1 conv g1 - exact)", (opf1 - exact1).norm2());
      error = (opf1 - exact1).norm2();
      if (abs(error) > FunctionDefaults<NDIM>::get_thresh())
        success++;
      print("success 2b ", success);

      real_function_1d q = real_factory_1d(world).f(h);
      error = q.trace();
      print("error in integral(h) ", error);
      if (abs(error) > FunctionDefaults<NDIM>::get_thresh())
        success++;
      print("success 3 ", success);

      error = q.norm2() - sqrt(sqrt(2.0 * constants::pi));
      print("error in norm2(h)", error);
      if (abs(error) > FunctionDefaults<NDIM>::get_thresh())
        success++;
      print("success 4 ", success);

      real_function_1d opq = op(q);
      exact = real_factory_1d(world).f(gconvh);
      error = (opq - exact).norm2();
      print("norm2(g conv h - exact)", error);
      if (abs(error) > FunctionDefaults<NDIM>::get_thresh())
        success++;
      print("success 5 ", success);

      // test the convolution with a derivative Gaussian:
      // result(y) = \int g'(x-y) f(x) = \int g(x-y) f'(x)
      // where we use
      // f(x)      = exp(-x^2) / sqrt(pi)
      // f'(x)     = -2 x exp(-x^2) / sqrt(pi)
      // result(y) = -y exp(-y/2) / sqrt(2 pi)

      // the derivative Gaussian convolution kernel:
      // note the scaling of the coeffs because the derivative operator brings
      // down the scaling factor of the exponent
      ops[0].reset(new GaussianConvolution1D<double>(
          k, 1.0 / sqrt(constants::pi), width * width, 1, false));

      real_convolution_1d oph(world, ops);

      // this is the result hardwired
      const real_function_1d convpg = real_factory_1d(world).f(conv_prime_g);

      // apply the derivative Gaussian on f
      opq = oph(f);

      // apply the Gaussian on the derivative of f
      const real_function_1d fp = real_factory_1d(world).f(gprime);
      real_function_1d opfp = op(fp);

      // the error
      const double error1 = (opq - convpg).norm2();
      const double error2 = (opfp - convpg).norm2();
      print("norm2(conv' g - exact)", error1);
      print("norm2(conv g' - exact)", error2);
      if (abs(error1) > FunctionDefaults<NDIM>::get_thresh())
        success++;
      print("success 6a ", success);
      if (abs(error2) > FunctionDefaults<NDIM>::get_thresh())
        success++;
      print("success 6b ", success);

      plot_line("opf.dat", npts, lo, hi, q, opq, convpg);
    }

    ///////////////////////////////////////////////////
    // test convolution with range restricted Gaussians
    ///////////////////////////////////////////////////

    // validate BoxSurfaceDisplacementRange by making sure
    // - it does not produce duplicate displacements
    {
      constexpr auto ND = 2;
      const int n = 4;
      Key<ND> key(n, Vector<Translation, ND>(1<<(n-1)));

      std::size_t disp_count = 0;
      const auto process_surface_displacements =
          [&]() {
            std::array<std::optional<std::int64_t>, ND> box_radius;
            std::array<std::optional<std::int64_t>, ND> surface_thickness;
            for (int d = 0; d != ND; ++d) {
              box_radius[d] = 1;
              surface_thickness[d] = 1;
            }
            BoxSurfaceDisplacementRange<ND> range_boundary_face_displacements(
                key, box_radius, surface_thickness, array_of_bools<ND>{false},
                [&](const auto level, const auto &dest, const auto &displacement) -> bool {
                  // skip displacements not in domain
                  const auto twon = (1 << level);
                  for (auto d = 0; d != ND; ++d) {
                    if (dest[d].has_value()) {
                      if (dest[d] < 0 ||
                          dest[d] >= twon)
                        return false;
                    }
                  }
                  return true;
                });
            std::vector<Key<ND>> disps;
            for (auto &&disp : range_boundary_face_displacements) {
//              std::cout << "disp = " << disp << " dest = " << key.neighbor(disp)
//                        << std::endl;
              disps.push_back(disp);
              ++disp_count;
            }
            std::sort(disps.begin(), disps.end());
            auto it = std::unique(disps.begin(), disps.end());

            if (it != disps.end()) {
              std::cout << "Duplicates found!!" << std::endl;
              success++;
            }
          };
      process_surface_displacements();
    }

    // now test the convolutions
    if constexpr (NDIM == 1) {

      int nerrors = 0;

      // gaussian exponents in *user* coordinates
      const std::vector<double> gaussian_exponents = {1 / (width * width), 1., 1e4, 1e8};
      std::vector<std::pair<real_function_t, real_function_t>> gaussians_01;
      for (const auto gaussian_exponent : gaussian_exponents) {
        auto factory = real_factory_t(world);
        auto make_gaussian = [&](const auto O) {
          real_function_1d result =
              factory.special_points(std::vector<coord_t>{coord_t(O)}).initial_level(15)
                  .functor(G(gaussian_exponent, O));
          result.truncate();
          return result;
        };
        gaussians_01.emplace_back(make_gaussian(0), make_gaussian(1));
      }

      const std::vector<double> kernel_exponents = {1/(width*width), 1., 25., 1e4};

      // some tests do not pass with default tolerances
      using α_β_O_opstr_t = std::tuple<double, double, std::size_t,std::string>;
      const std::map<α_β_O_opstr_t,double> custom_thresholds =
          {
              {std::tuple{kernel_exponents[1], gaussian_exponents[2], 0, "RNP"}, 2e-4},
              {std::tuple{kernel_exponents[3], gaussian_exponents[0], 0, "NP"}, 2e-4},
              {std::tuple{kernel_exponents[3], gaussian_exponents[0], 1, "NP"}, 2e-4},
              {std::tuple{kernel_exponents[3], gaussian_exponents[0], 0, "RNP"}, 2e-4},
              {std::tuple{kernel_exponents[3], gaussian_exponents[0], 1, "RNP"}, 2e-4}
          };

      // kernel exponents in *simulation* coordinates
      for (const auto kernel_exponent : kernel_exponents) {
        for (const auto sigma : {0}) {
          KernelRange range(1, sigma / (2 * L));

          // nonperiodic range-unlimited kernel
          ops[0].reset(new GaussianConvolution1D<double>(
              k, sqrt(kernel_exponent / constants::pi), kernel_exponent,
              /* deriv */ 0,
              /* lattice summed? */ false,
              /* bloch_k */ 0.)); // range unrestricted
          real_convolution_t opnp(world, ops);

          // nonperiodic range-restricted kernel
          ops[0].reset(new GaussianConvolution1D<double>(
              k, sqrt(kernel_exponent / constants::pi), kernel_exponent,
              /* deriv */ 0,
              /* lattice summed? */ false, /* bloch_k */ 0.,
              range)); // range restricted
          real_convolution_t oprnp(world, ops);
          oprnp.set_domain_periodicity(array_of_bools<NDIM>{false});
          // periodic range-restricted kernel implemented by expicit lattice sum
          real_convolution_t oprp2(world, ops);
          oprp2.set_domain_periodicity(array_of_bools<NDIM>{true});

          // periodic range-restricted kernel
          ops[0].reset(new GaussianConvolution1D<double>(
              k, sqrt(kernel_exponent / constants::pi), kernel_exponent,
              /* deriv */ 0,
              /* lattice summed? */ true, /* bloch_k */ 0.,
              range)); // range restricted
          real_convolution_t oprp(world, ops);
          FunctionDefaults<NDIM>::set_bc(BC_PERIODIC);

          // print out norms vs displacement length
          //      std::cout << "RNP norms\n";
          //      for (int n = 0; n <= 10; ++n) {
          //        for (int l = 0; l <= ((1 << n) + 5); ++l) {
          //          std::cout << "n=" << n << " l=" << l << " ||op_{n,l}||="
          //                    << oprnp.norm(n, Key<NDIM>(n, Vector<Translation, NDIM>(l)),
          //                                  Key<NDIM>(n, Vector<Translation, NDIM>(0)))
          //                    << "\n";
          //        }
          //      }
          //      std::cout << std::endl;

          for (size_t ig = 0; ig != gaussians_01.size(); ++ig) {
            const auto &[g, g1] = gaussians_01[ig];
            const auto gaussian_exponent = gaussian_exponents[ig];

            auto check = [&](const auto &op, const auto &f,
                             const std::size_t R) {
              const auto rr = op.range_restricted();
              const auto lattice_summed = op.lattice_summed().any();
              const auto periodic = lattice_summed || op.domain_is_periodic().any();
              if (!rr)
                MADNESS_ASSERT(!periodic);
              const std::string opstr = !rr ? "NP" : (periodic ? (lattice_summed ? "RP" : "RP2") : "RNP");
              const auto V = op(f);
              real_function_t V_exact;
              auto factory = real_factory_t(world).special_points(
                  std::vector<coord_t>{coord_t(R)});
              if (!rr) {
                V_exact = factory.functor(GConvGNP(kernel_exponent / (width * width), gaussian_exponent, R, L));
              } else {
                if (!periodic) {
                  V_exact = factory.functor(GConvGRNP(kernel_exponent / (width * width), gaussian_exponent, R, L));
                } else {
                  V_exact = factory.functor(GConvGRP(kernel_exponent / (width * width), gaussian_exponent, R, L));
                }
              }
              V_exact.truncate();

              const auto error2 = (V - V_exact).norm2();
              std::cout << opstr << "[α=" << (kernel_exponent/(width*width));
              if (op.range_restricted())
                std::cout << ", σ=" << sigma;
              std::cout << "](|x-y|) * g[" << gaussian_exponent << "](x";
              if (R > 0)
                std::cout << "-" << R;
              std::cout << "): ||mra-exact|| = " << std::scientific << error2;
              const auto default_error2_tol = thresh * 1e2;
              // look for overriden tolerances
              auto it = custom_thresholds.find(std::tuple{kernel_exponent, gaussian_exponent, R, opstr});
              const auto error2_tol = (it != custom_thresholds.end()) ? it->second : default_error2_tol;
              if (error2 > error2_tol) {
                ++nerrors;
                std::cout << " (FAIL)";
              }
              std::cout << std::endl;

              // optional plotting
              if constexpr (false) {
                const std::string ndstr =
                    (std::string("-") + std::to_string(NDIM) + "d");
                const std::string Rstr = "-R=" + std::to_string(R);
                const std::string ostr = "-" + opstr;
                const std::string sigmastr =
                    sigma > 0 ? (std::string("-σ=") + std::to_string(sigma)) : "";
                const std::string suffix = ndstr + ostr + Rstr + sigmastr;

                auto to_string = [](const double &x) {
                  if (double(int64_t(x)) == x) {
                    return std::to_string((int64_t)x);
                  }
                  else {
                    char buf[100];
                    if (x >= 1) {
                      snprintf(buf, 100, "%f", x);
                    }
                    else {
                      const auto width = (int)(-floor(log10(x)));
                      std::string fmt = "%." + std::to_string(width) + "f";
                      snprintf(buf, 100, fmt.c_str(), x);
                    }
                    return std::string(buf);
                  }
                };
                plot_line((std::string("gconv-") + to_string(kernel_exponent) +
                           "-g-" + to_string(gaussian_exponent*width*width) + suffix + ".dat")
                              .c_str(),
                          npts, lo, hi, V, V_exact);
              }
            };

            auto check_rp_vs_rp2 = [&](const auto &oprp, const auto &oprp2, const auto &f,
                                       const std::size_t R) {
              const auto V = oprp(f);
              const auto V2 = oprp2(f);

              const auto error2 = (V - V2).norm2();
              std::cout << "op[α=" << (kernel_exponent / (width * width));
              if (op.range_restricted())
                std::cout << ", σ=" << sigma;
              std::cout << "](|x-y|) * g[" << gaussian_exponent << "](x";
              if (R > 0)
                std::cout << "-" << R;
              std::cout << "): ||RP-RP2|| = " << std::scientific << error2;
              const auto default_error2_tol = thresh * 1e2;
              const auto error2_tol = default_error2_tol;
              if (error2 > error2_tol) {
                ++nerrors;
                std::cout << " (FAIL)";
              }
              std::cout << std::endl;
            };

            check(opnp, g, 0);
            check(opnp, g1, 1);
            check(oprnp, g, 0);
            check(oprnp, g1, 1);
            check(oprp, g, 0);
            check(oprp, g1, 1);
            check(oprp2, g, 0);
            check(oprp2, g1, 1);

            check_rp_vs_rp2(oprp, oprp2, g, 0);
            check_rp_vs_rp2(oprp, oprp2, g1, 1);

          }

          // linearity test
          {
            std::vector<std::reference_wrapper<const real_convolution_1d>> ops{
                opnp, oprnp, oprp, oprp2};
            const std::vector<std::string> op_labels = {"NP", "RNP", "RP", "RP2"};
            for(size_t iop=0; iop!=ops.size(); ++iop) {
              const auto &op = ops[iop].get();
              const auto &opstr = op_labels[iop];
              // evaluate op(f1+ f2 + ...) - op(f1) - op(f2) - ...
              real_function_1d error_0, fsum_0;
              real_function_1d error_1, fsum_1;
              for (size_t i = 0; i != gaussians_01.size(); ++i) {
                const auto &[g, g1] = gaussians_01[i];
                if (i == 0) {
                  fsum_0 = copy(g);
                  fsum_1 = copy(g1);
                  error_0 = -1 * op(g);
                  error_1 = -1 * op(g1);
                } else {
                  fsum_0 += g;
                  fsum_1 += g1;
                  error_0 -= op(g);
                  error_1 -= op(g1);
                }
              }
              error_0 += op(fsum_0);
              error_1 += op(fsum_1);

              const auto error_tol = thresh * 3e2;
              const auto error_0_norm2 = error_0.norm2();
              std::cout << opstr
                        << " linearity test 0: ||op(f1 + f2 + ...) - op(f1) - op(f2) - ...|| = "
                        << error_0_norm2;
              if (error_0_norm2 > error_tol) {
                ++nerrors;
                std::cout << " (FAIL)";
              }
              std::cout << std::endl;
              const auto error_1_norm2 = error_1.norm2();
              std::cout << opstr
                        << " linearity test 1: ||op(f1 + f2 + ...) - op(f1) - op(f2) - ...|| = "
                        << error_1_norm2;
              if (error_1_norm2 > error_tol) {
                ++nerrors;
                std::cout << " (FAIL)";
              }
              std::cout << std::endl;
            }
          }

        }
      }

      success += nerrors;
    }

    world.gop.fence();
    return success;
}


int main(int argc, char**argv) {
    auto&& world = initialize(argc,argv);

    int success=0;
    try {
        startup(world,argc,argv);

        FunctionDefaults<NDIM>::set_cubic_cell(-L,L);
        FunctionDefaults<NDIM>::set_k(k);
        FunctionDefaults<NDIM>::set_thresh(thresh);
        FunctionDefaults<NDIM>::set_refine(true);
        FunctionDefaults<NDIM>::set_initial_level(5);
        FunctionDefaults<NDIM>::set_truncate_mode(1);
        if (world.rank()==0) {
        	print(" threshold  ", thresh);
        	print(" polynomial ", k,"\n");
        }

        FunctionDefaults<2>::set_cubic_cell(-L,L);
        FunctionDefaults<2>::set_k(k);
        FunctionDefaults<2>::set_thresh(thresh);
        FunctionDefaults<2>::set_refine(true);
        FunctionDefaults<2>::set_initial_level(5);
        FunctionDefaults<2>::set_truncate_mode(1);

        success+=test_gconv(world);

    }
    catch (const SafeMPI::Exception& e) {
        print(e);
        error("caught an MPI exception");
    }
    catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    }
    catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    }
    catch (char* s) {
        print(s);
        error("caught a c-string exception");
    }
//    catch (const char* s) {
//        print(s);
//        error("caught a c-string exception");
//    }
    catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    }
    catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    }
    catch (...) {
        error("caught unhandled exception");
    }

    world.gop.fence();
    finalize();

    return success;
}

