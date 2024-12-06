//
// Created by Eduard Valeyev on 12/4/24.
//

#ifndef MADNESS_ARRAY_OF_BOOLS_H
#define MADNESS_ARRAY_OF_BOOLS_H

#include <cstddef>
#include <array>
#include <numeric>

namespace madness {

/// syntactic sugar for std::array<bool, N>

/// reason to exist: so can apply logical operations to the entire pack of bools, e.g. `a && b`,
/// and perform queries like `a.any()`, `a.none()`, etc.
template <std::size_t N>
class array_of_bools : public std::array<bool, N> {
public:
  using base_type = std::array<bool, N>;

  const auto& as_array() const { return static_cast<const base_type&>(*this); }
  auto& as_array() { return static_cast<base_type&>(*this); }

  /// default ctor public only for `N==0`
  template <std::size_t NN = N, typename = std::enable_if_t<NN==0>>
  array_of_bools() {}

  /// constructs an array filled with \p v
  explicit array_of_bools(bool v) { as_array().fill(v); }

  /// constructs an array filled with `{v, vs...}`
  template <typename ... Bools, typename = std::enable_if_t<sizeof...(Bools)+1==N>>
  explicit array_of_bools(bool v, Bools... vs) : base_type{{v, static_cast<bool>(vs)...}} {}

  /// @return true if any element is true
  bool any() const { return std::accumulate(this->begin(), this->end(), false, std::logical_or{}); }

  /// @return true if no elements are true
  bool none() const { return !any(); }

  /// @return true if all elements are true
  bool all() const { return std::accumulate(this->begin(), this->end(), false, std::logical_and{}); }

  /// @return first \p C elements
  template <std::size_t C, typename = std::enable_if_t<C <= N>>
  array_of_bools<C> front() const {
    array_of_bools<C> result;
    std::copy(this->begin(), this->begin()+C, result.begin());
    return result;
  }

  /// assigns \p a to the first \p C elements of this
  /// @param a the array to assign to the front of this
  /// @return reference to this object
  template <std::size_t C, typename = std::enable_if_t<C <= N>>
  array_of_bools& assign_front(const array_of_bools<C>& a) {
    std::copy(a.begin(), a.end(), this->begin());
    return *this;
  }

  /// @return array with first \p C elements obtained by logical AND between
  /// \p a and the first \p C elements of this, the rest filled with the
  /// remainder of this
  template <std::size_t C, typename = std::enable_if_t<C <= N>>
  array_of_bools and_front(const array_of_bools<C>& a) const {
    array_of_bools result;
    const auto it = std::transform(a.begin(), a.end(), this->begin(), result.begin(),
                             std::logical_and{});
    std::copy(this->begin() + C, this->end(), it);
    return result;
  }

  /// @return array with first \p C elements obtained by logical OR between
  /// \p a and the first \p C elements of this, the rest filled with the
  /// remainder of this
  template <std::size_t C, typename = std::enable_if_t<C <= N>>
  array_of_bools or_front(const array_of_bools<C>& a) const {
    array_of_bools result;
    const auto it = std::transform(a.begin(), a.end(), this->begin(), result.begin(),
                             std::logical_or{});
    std::copy(this->begin() + C, this->end(), it);
    return result;
  }

  /// @return last \p C elements
  template <std::size_t C, typename = std::enable_if_t<C <= N>>
  array_of_bools<C> back() const {
    array_of_bools<C> result;
    std::copy(this->begin() + (N - C), this->end(), result.begin());
    return result;
  }

  /// assigns \p a to the last \p C elements of this
  /// @param a the array to assign to the back of this
  /// @return reference to this object
  template <std::size_t C, typename = std::enable_if_t<C <= N>>
  array_of_bools& assign_back(const array_of_bools<C>& a) {
    std::copy(a.begin(), a.end(), this->begin() + (N - C));
    return *this;
  }

  /// @return array with last \p C elements obtained by logical AND between
  /// \p a and the last \p C elements of this, the rest filled with the
  /// remainder of this
  template <std::size_t C, typename = std::enable_if_t<C <= N>>
  array_of_bools and_back(const array_of_bools<C>& a) const {
    array_of_bools result;
    const auto it =
        std::copy(this->begin(), this->begin() + (N - C), result.begin());
    std::transform(a.begin(), a.end(), this->begin() + (N - C), it,
                   std::logical_and{});
    return result;
  }

  /// @return array with last \p C elements obtained by logical OR between
  /// \p a and the last \p C elements of this, the rest filled with the
  /// remainder of this
  template <std::size_t C, typename = std::enable_if_t<C <= N>>
  array_of_bools or_back(const array_of_bools<C>& a) const {
    array_of_bools result;
    const auto it =
        std::copy(this->begin(), this->begin() + (N - C), result.begin());
    std::transform(a.begin(), a.end(), this->begin() + (N - C), it,
                   std::logical_or{});
    return result;
  }

  friend array_of_bools<N> operator&&(const array_of_bools<N>& a, const array_of_bools<N>& b) {
    array_of_bools<N> result;
    std::transform(a.begin(), a.end(), b.begin(), result.begin(), std::logical_and{});
    return result;
  }
  friend array_of_bools<N> operator||(const array_of_bools<N>& a, const array_of_bools<N>& b) {
    array_of_bools<N> result;
    std::transform(a.begin(), a.end(), b.begin(), result.begin(), std::logical_or{});
    return result;
  }
  friend array_of_bools<N> operator^(const array_of_bools<N>& a, const array_of_bools<N>& b) {
    array_of_bools<N> result;
    std::transform(a.begin(), a.end(), b.begin(), result.begin(), [](auto b1, auto b2) {
      return b1 ^ b2;
    });
    return result;
  }
  friend array_of_bools<N> operator!(const array_of_bools<N>& a) {
    array_of_bools<N> result;
    std::transform(a.begin(), a.end(), result.begin(), std::logical_not{});
    return result;
  }

  friend array_of_bools<N> operator&&(const array_of_bools<N>& a, bool b) {
    array_of_bools<N> result;
    std::transform(a.begin(), a.end(), result.begin(), [b](auto a_v) {
      return a_v && b;
    });
    return result;
  }
  friend array_of_bools<N> operator&&(bool b, const array_of_bools<N>& a) {
    return a && b;
  }
  friend array_of_bools<N> operator||(const array_of_bools<N>& a, bool b) {
    array_of_bools<N> result;
    std::transform(a.begin(), a.end(), result.begin(), [b](auto a_v) {
      return a_v || b;
    });
    return result;
  }
  friend array_of_bools<N> operator||(bool b, const array_of_bools<N>& a) {
    return a || b;
  }
  friend array_of_bools<N> operator^(const array_of_bools<N>& a, bool b) {
    array_of_bools<N> result;
    std::transform(a.begin(), a.end(), result.begin(), [b](auto a_v) {
      return a_v ^ b;
    });
    return result;
  }
  friend array_of_bools<N> operator^(bool b, const array_of_bools<N>& a) {
    return a ^ b;
  }

private:
  // "default" ctor for N!=0 only for internal use, need tagged dispatch
  struct nonempty_default_ctor_tag{};
  template <std::size_t NN = N, typename = std::enable_if_t<NN!=0>>
  explicit array_of_bools(nonempty_default_ctor_tag = {}) {}

  template <std::size_t NN>
  friend class array_of_bools;
};

static_assert(std::is_same_v<array_of_bools<1>::value_type, bool>);
static_assert(std::is_same_v<decltype(array_of_bools<1>{true}.data()), bool*>);
static_assert(std::is_same_v<decltype(array_of_bools<2>{true, false}[0]), bool&>);
static_assert(std::is_same_v<decltype(array_of_bools<0>{} && array_of_bools<0>{}), array_of_bools<0>>);
static_assert(std::is_same_v<decltype(array_of_bools<0>{} || array_of_bools<0>{}), array_of_bools<0>>);
static_assert(std::is_same_v<decltype(array_of_bools<0>{} ^ array_of_bools<0>{}), array_of_bools<0>>);
static_assert(std::is_same_v<decltype(!array_of_bools<0>{}), array_of_bools<0>>);
static_assert(std::is_same_v<decltype(array_of_bools<0>{} && true), array_of_bools<0>>);
static_assert(std::is_same_v<decltype(true && array_of_bools<0>{}), array_of_bools<0>>);
static_assert(std::is_same_v<decltype(array_of_bools<0>{} || true), array_of_bools<0>>);
static_assert(std::is_same_v<decltype(true || array_of_bools<0>{}), array_of_bools<0>>);
static_assert(std::is_same_v<decltype(array_of_bools<0>{} ^ true), array_of_bools<0>>);
static_assert(std::is_same_v<decltype(true ^ array_of_bools<0>{}), array_of_bools<0>>);

}

#endif // MADNESS_ARRAY_OF_BOOLS_H
