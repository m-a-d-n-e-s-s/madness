/*
 *   Small class to hold response functions and to interact with KAIN solver.
 */

#ifndef SRC_APPS_ADRIAN_RESPONSEFUNCTION2_H_
#define SRC_APPS_ADRIAN_RESPONSEFUNCTION2_H_

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <vector>

// real_function_3d == Function<double,3>
namespace madness {
class responseVector {
public:
  unsigned int g_states;

  vector_real_function_3d vec;

public:
  // default constructor
  responseVector() : g_states(0) {}
  // zero constructor
  responseVector(World& world, unsigned int num_orbitals)
      : g_states(num_orbitals),
        vec(zero_functions<double, 3>(world, g_states)) {}
  // copy constructor
  responseVector(const responseVector& y)
      : g_states(y.g_states), vec(madness::copy(y[0].world(), y.vec)) {
    vec[0].world().gop.fence();
  }
  // this time we are going to copy a vector_real_function_3d
  explicit responseVector(const vector_real_function_3d vec)
      : g_states(vec.size()), vec(madness::copy(vec[0].world(), vec)) {
    vec[0].world().gop.fence();
  };
  // 1D accessor for vec  returns the function at position i
  real_function_3d& operator[](int64_t i) { return vec.at(i); }
  const real_function_3d& operator[](int64_t i) const { return vec.at(i); }
  bool same_size(const responseVector& y) const {
    // teseting asdfasdfad  iasdfasdf
    //
    //
    //
    // asdf   asd  a a
    //
    return (g_states == y.g_states);
  }

  // Kain
  responseVector operator+(const responseVector& b) {
    // assert that b is the same size of our vector vec
    MADNESS_ASSERT(same_size(b));
    responseVector result(vec[0].world(), g_states);
    result.vec = add(vec[0].world(), vec, b.vec);
    result[0].world().gop.fence();
    return result;
  }
  // assert that b is the same size of our vector x
  responseVector operator-(const responseVector& b) {
    MADNESS_ASSERT(same_size(b));
    responseVector result(vec[0].world(), g_states);
    result.vec = sub(vec[0].world(), vec, b.vec);
    result[0].world().gop.fence();
    return result;
  }
  //+=
  responseVector operator+=(const responseVector& b) {
    MADNESS_ASSERT(same_size(b));
    vec = add(vec[0].world(), vec, b.vec);
    return *this;
  }
  responseVector operator+=(responseVector& b) {
    MADNESS_ASSERT(same_size(b));
    vec = add(vec[0].world(), vec, b.vec);
    return *this;
  }
  // mul * a*V
  responseVector operator*(double a) {
    responseVector result(*this);
    madness::scale(vec[0].world(), result.vec, a, false);
    return result;
    ;
  }
  responseVector operator*=(double a) {
    madness::scale(vec[0].world(), this->vec, a, false);
    return *this;
  }
  // mul f(x)*V
  responseVector operator*(real_function_3d f) {
    responseVector result(*this);
    /// Multiplies a function against a vector of functions --- q[i] = a * v[i]
    result.vec = mul(f.world(), f, result.vec, false);
    return result;
  }
  responseVector copy() const {
    responseVector result(vec[0].world(), g_states);
    result.vec = madness::copy(vec[0].world(), vec, false);
    return result;
  }
  // return size
  unsigned int size() const { return g_states; }
  void zero() { vec = zero_functions<double, 3>(vec[0].world(), g_states); }
  void compress_vec() { compress(vec[0].world(), vec, true); };
  void reconstruct_vec() { reconstruct(vec[0].world(), vec, true); };
  void truncate_vec() { truncate(vec[0].world(), vec); };
  double norm2() { return sqrt(inner(vec, vec)); }
};

inline double inner(responseVector& a, responseVector& b) {
  MADNESS_ASSERT(a.size() > 0);
  MADNESS_ASSERT(a.size() == b.size());
  double value = 0.0;
  // inner of two sets of vector functions
  value = inner(a, b);
  return value;
}
// real_function_3d dot(rvec a, rvec vf) { return dot(a[0].world(), a, vf); }

typedef std::vector<vector_real_function_3d> response_matrix;

class ResponseVectors {
  // Member variables
public:
  int r_states; // Num. of resp. states
  int g_states; // Num. of ground states
  response_matrix x;

  // Member functions
public:
  // Default constructor
  ResponseVectors() : r_states(0), g_states(0) {}

  // Initializes functions to zero
  // m = number of response states
  // n = number of ground state orbitals

  // Zero Constructor
  ResponseVectors(World& world, int num_states, int num_orbitals)
      : r_states(num_states), g_states(num_orbitals) {
    for (unsigned int i = 0; i < num_states; i++) {
      x.push_back(zero_functions<double, 3>(world, g_states, true));
    }
    x[0][0].world().gop.fence();
  }

  // Copy constructor
  ResponseVectors(const ResponseVectors& rf_copy)
      : r_states(rf_copy.r_states), g_states(rf_copy.g_states) {
    for (unsigned int i = 0; i < r_states; i++) {
      x.push_back(madness::copy(rf_copy[0][0].world(), rf_copy[i]));
    }
    x[0][0].world().gop.fence();
  }

  // Determines if two ResponseFunctions are the same size
  bool same_size(const ResponseVectors& rf_copy) const {
    return (r_states == rf_copy.r_states && g_states == rf_copy.g_states);
  }

  // 1D accessor for x
  // std::vector<Function<double, 3>> &operator[](int64_t i) { return x[i]; }
  vector_real_function_3d& operator[](int64_t i) { return x.at(i); }

  const vector_real_function_3d& operator[](int64_t i) const { return x.at(i); }

  // KAIN must have this
  // element wise addition.  we add each vector separatly
  ResponseVectors operator+(const ResponseVectors& b) {
    MADNESS_ASSERT(same_size(b));

    ResponseVectors result(x[0][0].world(), r_states, g_states);

    for (unsigned int i = 0; i < r_states; i++) {
      result[i] = x[i] + b[i];
    }

    result[0][0].world().gop.fence();
    return result;
  }

  friend ResponseVectors operator-(const ResponseVectors& a,
                                   const ResponseVectors& b) {
    MADNESS_ASSERT(a.same_size(b));

    ResponseVectors result(a[0][0].world(), a.r_states, a.g_states);

    for (unsigned int i = 0; i < a.r_states; i++) {
      result[i] = a[i] - b[i];
    }
  }

  // KAIN must have this
  // Scaling by a constant
  ResponseVectors operator*(double a) {
    ResponseVectors result(*this);

    for (unsigned int i = 0; i < r_states; i++) {
      x[i] = x[i] * a;
    }

    x[0][0].world().gop.fence();
    return result;
  }

  // Scaling all internal functions by an external function
  // g[i][j] = x[i][j] * f
  ResponseVectors operator*(const Function<double, 3>& f) {
    ResponseVectors result(*this);

    for (unsigned int i = 0; i < r_states; i++) {
      // Using vmra.h funciton
      result[i] = x[i] * f;
    }

    f.world().gop.fence();
    return result;
  }

  // KAIN must have this
  ResponseVectors operator+=(const ResponseVectors b) {
    MADNESS_ASSERT(same_size(b));

    for (unsigned int i = 0; i < r_states; i++) {
      x[i] += b[i];
    }

    return *this;
  }

  // Returns a deep copy
  ResponseVectors copy() const {
    ResponseVectors result(x[0][0].world(), r_states, g_states);

    for (unsigned int i = 0; i < r_states; i++) {
      result.x[i] = madness::copy(x[0][0].world(), x[i]);
    }
    x[0][0].world().gop.fence();

    return result;
  }

  // Mimicking std::vector with these 4
  void push_back(const vector_real_function_3d& f) {
    x.push_back(f);
    r_states++;

    // Be smart with g_states
    if (g_states > 0) {
      MADNESS_ASSERT(g_states = f.size());
    } else { // g_states == 0 (empty vector)
      g_states = f.size();
    }
  }
  void pop_back() {
    MADNESS_ASSERT(r_states >= 1);
    x.pop_back();
    r_states--;

    // Be smart with g_states
    if (r_states == 0) { // removed last item
      g_states = 0;
    }
  }
  void clear() {
    x.clear();
    r_states = 0;
    g_states = 0;
  }
  unsigned int size() const { return r_states; }

  // Mimicing standard madness calls with these 3
  void zero() {
    for (unsigned int k = 0; k < r_states; k++) {
      x[k] = zero_functions<double, 3>(x[0][0].world(), g_states);
    }
  }

  void compress_rf() {
    for (unsigned int k = 0; k < r_states; k++) {
      compress(x[0][0].world(), x[k], true);
    }
  }

  void reconstruct_rf() {
    for (unsigned int k = 0; k < r_states; k++) {
      reconstruct(x[0][0].world(), x[k], true);
    }
  }

  void truncate_rf() {
    for (unsigned int k = 0; k < r_states; k++) {
      truncate(x[0][0].world(), x[k], FunctionDefaults<3>::get_thresh(), true);
    }
  }
  void truncate_rf(double tol) {
    for (unsigned int k = 0; k < r_states; k++) {
      truncate(x[0][0].world(), x[k], tol, true);
    }
  }

  // Returns norms of each state
  Tensor<double> norm2() {
    Tensor<double> answer(r_states);
    for (unsigned int i = 0; i < r_states; i++) {
      answer(i) = madness::norm2(x[0][0].world(), x[i]);
    }
    return answer;
  }

  // Scales each state (read: entire row) by corresponding vector element
  //     new[i] = old[i] * mat[i]
  void scale(Tensor<double>& mat) {
    for (unsigned int i = 0; i < r_states; i++)
      x[i] = x[i] * mat[i];
    x[0][0].world().gop.fence();
  }
};
// Final piece for KAIN
inline double inner(ResponseVectors& a, ResponseVectors& b) {
  MADNESS_ASSERT(a.size() > 0);
  MADNESS_ASSERT(a.size() == b.size());
  MADNESS_ASSERT(a[0].size() > 0);
  MADNESS_ASSERT(a[0].size() == b[0].size());

  double value = 0.0;

  for (unsigned int i = 0; i < a.size(); i++) {
    // vmra.h function
    value += inner(a[i], b[i]);
  }
  // how do i two sets inner two response vectors?
  // The response vectors hold N response vectors for perturbation
  // To take the inner between 2 response vectors we need to sum
  // the components and then inner.  we get left with N inner values
  // Do we sum each of them? no I don't believe so.

  return value;
}

} // End namespace madness
#endif // SRC_APPS_ADRIAN_RESPONSEFUNCTION2_H_

// Deuces
