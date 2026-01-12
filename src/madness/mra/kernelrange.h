//
// Created by Eduard Valeyev on 1/3/25.
//

#ifndef MADNESS_MRA_KERNELRANGE_H__INCLUDED
#define MADNESS_MRA_KERNELRANGE_H__INCLUDED

#include <madness/world/madness_exception.h>
#include <madness/world/worldhash.h>

#include <cmath>
#include <optional>

namespace madness {

/// Class to control the maximum range of lattice summation.
/// Only a single integer is needed.
class LatticeRange {
  int range_ = 0; // Initialize to OBC
public:
  LatticeRange() = default;
  explicit LatticeRange(int range) {
    set_range(range);
  }
  explicit LatticeRange(bool lattice_sum) {
    if (lattice_sum) set_range_inf();
  }

  void set_range(int range) {
    if (range == std::numeric_limits<int>::max()) range_ = range;
    else {
      MADNESS_ASSERT(range_ >= 0 && range % 2 == 1);
      range_ = (range - 1) / 2;
    }
  }

  // Convenience function to make lattice range inactive.
  void set_range_inf() {
    set_range(std::numeric_limits<int>::max());
  }

  [[nodiscard]] int get_range() const {return range_;}

  [[nodiscard]] bool infinite() const {return range_ == std::numeric_limits<int>::max();}

  /// @return true if range is limited
  explicit operator bool() const { return static_cast<bool>(range_); }
};

/// To limit the range of kernel K(x-y) it is multiplied (in user coordinates) by a restrictor function
/// \f$ r(N/2 - |x-y|) \f$, where \f$ r(x) \f$ is identity for unrestricted kernel or one of the choices
/// encoded by KernelRange::Type
class KernelRange {
public:

  /// types of range restrictor functions
  enum Type {
    /// Heavyside (step) function: \f$ r(x) \equiv \theta(x)\f$
    Hard,
    /// erf approximation to Heavyside function: \f$ r(x) \equiv \theta_\sigma(x) \equiv (1 + \erf(x/\sigma))/2 \f$;
    /// in the limit $\sigma\to 0$  this becomes the Heaviside function.
    SoftErf
  };

  /// restrictor function
  struct Restrictor {
    //Restrictor() = default;
    //Restrictor(Type type) : type_(type) { MADNESS_ASSERT(type == Hard); }
      Restrictor(Type type = Hard, double sigma=1.0) : type_(type) { // Avoid uninitialized warning on sigma if type is Hard 
      MADNESS_ASSERT(sigma >= 0);
      // if (sigma == 0)
      //   MADNESS_ASSERT(type == Hard);
      // else {
      //   MADNESS_ASSERT(type != Hard);
        this->sigma_w_inverse_ = {sigma, 1./sigma};
      // }
    }

    Type type() const { return type_; }
    double sigma() const {
      MADNESS_ASSERT(sigma_w_inverse_);
      return sigma_w_inverse_->first;
    }
    double sigma_inverse() const {
      MADNESS_ASSERT(sigma_w_inverse_);
      return sigma_w_inverse_->second;
    }

    /// @return value of restrictor function r(x) at x
    double value(double x) const {
      auto f = [](double x, double ooσ) {
        return (1 + std::erf(x * ooσ)) * 0.5;
      };

      switch (type_) {
      case Hard:
        // use regularized step function to ensure that value at x=0 is 1/2
        return f(x, 1e8);
      case SoftErf:
        return f(x, this->sigma_inverse());
      default:
        MADNESS_EXCEPTION("KernelRange: unknown restrictor type",1);
      }
    }

    /// @return true if \p r1 and \p r2 are equal
    friend bool operator==(const Restrictor& r1, const Restrictor& r2) {
      return r1.type_ == r2.type_ && (r1.type_ == KernelRange::Type::Hard ||
                                      (r1.type_ != KernelRange::Type::Hard &&
                                       r1.sigma() == r2.sigma()));
    }

    hashT hash() const {
      hashT result = hash_value((int)type_);
      if (sigma())
        hash_combine(result, sigma());
      return result;
    }

  private:
    Type type_ = Hard;
    std::optional<std::pair<double,double>> sigma_w_inverse_;
  };

  /// constructs a null (i.e., infinite) kernel range
  /// @post `this->infinite()==true`
  KernelRange() = default;
  /// constructs a finite soft (`sigma > 0`) or hard (`sigma==0`) kernel range
  /// @param sigma regularization parameter (lengthscale in simulation [0,1] coordinate units) controls the softness of the range restrictor
  /// @pre `sigma>=0`
  /// @post `this->soft()==true`
  KernelRange(unsigned int N, double sigma) {
    MADNESS_ASSERT(sigma >= 0);
    if (sigma == 0)
      data.emplace(N, Restrictor{Hard});
    else
      data.emplace(N, Restrictor{SoftErf, sigma});
  }
  /// constructs a finite (hard) kernel range
  /// @post `this->hard()==true`
  KernelRange(unsigned int N) { data.emplace(N, Restrictor{}); }

  unsigned int N() const { return data.value().N; }
  Type type() const { return restrictor().type(); }
  double sigma() const { return restrictor().sigma(); }

  /// @return true if range is limited
  bool finite() const { return data.has_value(); }
  /// @return true if range is limited using a hard window function (Heaviside)
  bool finite_hard() const { return finite() && type() == Hard; }
  /// @return true if range is limited using a soft window function (sigmoid)
  bool finite_soft() const { return finite() && type() != Hard; }
  /// @return true if range is unlimited
  bool infinite() const { return !finite(); }

  /// @return true if range is limited
  explicit operator bool() const { return finite(); }

  /// @return true if \p r1 and \p r2 are equal
  friend bool operator==(const KernelRange& r1, const KernelRange& r2) {
    if (r1.finite() != r2.finite())
      return false;
    if (r1.finite())
      return r1.N() == r2.N() && r1.restrictor() == r2.restrictor();
    return true;
  }

  /// @return value of restrictor function at N/2 - abs(r)
  double value(double r) const {
    if (infinite())
      return 1.;
    else {
      auto x = N() * 0.5 - std::abs(r);
      return restrictor().value(x);
    }
  }

  hashT hash() const {
    hashT result = 0;
    if (!infinite()) {
      result = hash_value(N());
      hash_combine(result, restrictor());
    }
    return result;
  }

  static constexpr double extent_default_epsilon = std::numeric_limits<double>::epsilon();

  /// @return max value of `|x-y|` (rounded up, in units of 1/2)
  /// for which `r(N/2 - |x-y|)` is greater than @p epsilon
  int iextent_x2(double epsilon = extent_default_epsilon) const {
    MADNESS_ASSERT(epsilon > 0);
    if (infinite())
      return std::numeric_limits<int>::max();
    else {
      return data.value().iextent_x2(epsilon);
    }
  }

private:

  struct Data {
    unsigned int N;
    Restrictor restrictor;
    int iextent_x2_default;  // memoized extent(int_extent_default_epsilon)

    Data(unsigned int n, Restrictor r) : N(n), restrictor(r), iextent_x2_default(compute_iextent_x2(N,restrictor, extent_default_epsilon)) {}

    /// @return max value of `|x-y|` (rounded up, in units of 1/2)
    /// for which `r(N/2 - |x-y|)` is greater than @p epsilon
    int iextent_x2(double epsilon) const {
      if (epsilon == extent_default_epsilon)
        return iextent_x2_default;
      else
        return compute_iextent_x2(N, restrictor, epsilon);
    }

    /// @return max value of `|x-y|` (rounded up, in units of 1/2)
    /// for which `r(N/2 - |x-y|)` is greater than @p epsilon
    static int compute_iextent_x2(int N, const Restrictor& restrictor, double epsilon) {
      MADNESS_ASSERT(epsilon > 0);
      switch (restrictor.type()) {
      case Hard:
        return N;
      case SoftErf: {
        constexpr bool use_newton_solver = false;
        if (use_newton_solver) {
          auto f01 = [&, ooσ = 1. / restrictor.sigma(),
                      oosqrtpiσ =
                          1. / (sqrt(M_PI) * restrictor.sigma())](double x) {
            const auto arg = (0.5 * N - x) * ooσ;
            const auto fx = (1 + std::erf(arg)) * 0.5;
            const auto dfx = -std::exp(-arg * arg) * oosqrtpiσ;
            return std::make_pair(fx, dfx);
          };
          auto iter = 0;
          const auto max_iter = 50;
          auto x = 0.5;
          auto dx = 1.;
          double fx = 1.;
          // step-restricted newton iteration
          while ((std::abs(dx) > 0.01 || std::abs(fx) > epsilon) &&
                 iter < max_iter) {
            double dfx;
            std::tie(fx, dfx) = f01(x);
            dx = -(fx - epsilon) / dfx;
            auto restricted_step = [](const double step,
                                      const double maxabs_step) {
              auto sign = [](const auto x) { return x >= 0 ? 1. : -1.; };
              return std::abs(step) > maxabs_step ? sign(step) * maxabs_step
                                                  : step;
            };
            x += restricted_step(dx, 0.2);
            ++iter;
          }
          return static_cast<int>(ceil(x * 2.));
        } else {
          // keep increasing x until f(x) falls below epsilon
          int x = N;
          auto value = [&](const double r) {
            auto x = N * 0.5 - std::abs(r);
            return restrictor.value(x);
          };
          double fx = value(x * 0.5);
          while (fx >= epsilon) {
            ++x;
            fx = value(x * 0.5);
          }
          return x;
        }
      }
      default:
        MADNESS_EXCEPTION("KernelRange:: unknown restrictor type",1);
      }
    }

  };  // Data

  std::optional<Data> data;
  const Restrictor &restrictor() const { return data.value().restrictor; }
};

}  // madness

#endif // MADNESS_MRA_KERNELRANGE_H__INCLUDED
