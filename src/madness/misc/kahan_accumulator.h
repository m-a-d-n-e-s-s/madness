//
// Created by Eduard Valeyev on 10/02/21.
//

#ifndef MADNESS_MISC_KAHAN_ACCUMULATOR_H_
#define MADNESS_MISC_KAHAN_ACCUMULATOR_H_

#include <iosfwd>
#include <type_traits>

namespace madness {

template <typename T, typename Enabler = void>
struct KahanAccumulator;

/// implements Kahan summation for real numbers
template <typename Real>
struct KahanAccumulator<Real,
    std::enable_if_t<std::is_floating_point_v<Real>>> {
  KahanAccumulator() = default;
  KahanAccumulator(const KahanAccumulator&) = default;
  KahanAccumulator& operator=(const KahanAccumulator&) = default;

  template <typename Real_,
            typename = std::enable_if_t<std::is_floating_point_v<Real_>>>
  KahanAccumulator(Real_ v) : value_(v) {}

  KahanAccumulator(Real v, Real c) : value_(v), correction_(c) {}

  template <typename Real_>
  KahanAccumulator(KahanAccumulator<Real_> v)
      : value_(v.value_), correction_(v.correction_) {}

  explicit operator Real() const { return value_ + correction_; }

  template <typename Real_,
            typename = std::enable_if_t<std::is_floating_point_v<Real_>>>
  KahanAccumulator& operator+=(Real_ v) {
    volatile auto y = v - correction_;
    volatile auto t = value_ + y;
    correction_ = (t - value_) - y;
    value_ = t;
    return *this;
  }

  template <typename Real_,
            typename = std::enable_if_t<std::is_floating_point_v<Real_>>>
  KahanAccumulator& operator-=(Real_ v) {
    volatile auto minus_y = v + correction_;
    volatile auto t = value_ - minus_y;
    correction_ = (t - value_) + minus_y;
    value_ = t;
    return *this;
  }

  template <typename Real_>
  KahanAccumulator& operator+=(const KahanAccumulator<Real_>& v) {
    *this += v.correction_;
    *this += v.value_;
    return *this;
  }

  template <typename Real_>
  KahanAccumulator& operator-=(const KahanAccumulator<Real_>& v) {
    *this -= v.correction_;
    *this -= v.value_;
    return *this;
  }

  KahanAccumulator operator-() const {
    return KahanAccumulator(-value_, -correction_);
  }

  auto value() const { return value_; }
  auto correction() const { return correction_; }

  template <typename Archive>
  void serialize(Archive& ar) {
    ar& value_& correction_;
  }

 private:
  Real value_ = Real{0};
  Real correction_ = Real{0};
};

template <typename Real1, typename Real2>
auto operator+(KahanAccumulator<Real1> v1, Real2 v2) {
  KahanAccumulator<decltype(std::declval<Real1>() + std::declval<Real2>())>
      result(v1);
  result += v2;
  return result;
}

template <typename Real1, typename Real2>
auto operator+(Real2 v2, KahanAccumulator<Real1> v1) {
  KahanAccumulator<decltype(std::declval<Real1>() + std::declval<Real2>())>
      result(v1);
  result += v2;
  return result;
}

template <typename Real1, typename Real2>
auto operator+(KahanAccumulator<Real1> v1, KahanAccumulator<Real2> v2) {
  KahanAccumulator<decltype(std::declval<Real1>() + std::declval<Real2>())>
      result(v1);
  result += v2;
  return result;
}

template <typename Real1, typename Real2>
auto operator-(KahanAccumulator<Real1> v1, Real2 v2) {
  KahanAccumulator<decltype(std::declval<Real1>() - std::declval<Real2>())>
      result(v1);
  result -= v2;
  return result;
}

template <typename Real1, typename Real2>
auto operator-(Real2 v2, KahanAccumulator<Real1> v1) {
  KahanAccumulator<decltype(std::declval<Real2>() - std::declval<Real1>())>
      result(v2);
  result -= v1;
  return result;
}

template <typename Real1, typename Real2>
auto operator-(KahanAccumulator<Real1> v1, KahanAccumulator<Real2> v2) {
  KahanAccumulator<decltype(std::declval<Real1>() - std::declval<Real2>())>
      result(v1);
  result -= v2;
  return result;
}

template <typename Char, typename Real>
std::basic_ostream<Char>& operator<<(std::basic_ostream<Char>& os,
                                     const KahanAccumulator<Real>& v) {
  os << "{" << v.value() << "," << v.correction() << "}";
  return os;
}

/// implements Kahan summation for complex numbers
template <typename Complex>
struct KahanAccumulator<Complex,
                        std::enable_if_t<!std::is_floating_point_v<Complex>>> {
  using Real = typename Complex::value_type;
  using RealAccumulator = KahanAccumulator<Real>;

  KahanAccumulator() = default;
  KahanAccumulator(const KahanAccumulator&) = default;
  KahanAccumulator& operator=(const KahanAccumulator&) = default;

  template <typename Complex_,
            typename = std::enable_if_t<!std::is_floating_point_v<Complex_>>>
  KahanAccumulator(const Complex_& v) : real_(v.real()), imag_(v.imag()) {}

  template <typename Complex_>
  KahanAccumulator(const KahanAccumulator<Complex_>& v)
      : real_(v.real_), imag_(v.imag_) {}

  explicit operator Complex() const { return Complex(static_cast<Real>(real_), static_cast<Real>(imag_)); }

  template <typename Complex_,
            typename = std::enable_if_t<!std::is_floating_point_v<Complex_>>>
  KahanAccumulator& operator+=(const Complex_& v) {
    real_ += v.real();
    imag_ += v.imag();
    return *this;
  }

  template <typename Complex_,
            typename = std::enable_if_t<!std::is_floating_point_v<Complex_>>>
  KahanAccumulator& operator-=(const Complex_& v) {
    real_ -= v.real();
    imag_ -= v.imag();
    return *this;
  }

  template <typename Complex_>
  KahanAccumulator& operator+=(const KahanAccumulator<Complex_>& v) {
    real_ += v.real();
    imag_ += v.imag();
    return *this;
  }

  template <typename Complex_>
  KahanAccumulator& operator-=(const KahanAccumulator<Complex_>& v) {
    real_ -= v.real();
    imag_ -= v.imag();
    return *this;
  }

  template <typename Archive>
  void serialize(Archive& ar) {
    ar& real_& imag_;
  }

private:
  RealAccumulator real_ = {};
  RealAccumulator imag_ = {};
};

}  // namespace madness

#endif  // MADNESS_MISC_KAHAN_ACCUMULATOR_H_
