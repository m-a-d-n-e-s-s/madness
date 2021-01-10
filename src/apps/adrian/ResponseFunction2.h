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


typedef std::vector<vector_real_function_3d> response_matrix;

class response_space {
  // Member variables
public:
  size_t r_states; // Num. of resp. states
  size_t g_states; // Num. of ground states
  response_matrix x;

  // Member functions
public:
  // Default constructor
  response_space() : r_states(0), g_states(0), x(NULL) {}

  // Initializes functions to zero
  // m = number of response states
  // n = number of ground state orbitals

  // Zero Constructor
  response_space(World& world, size_t num_states, size_t num_orbitals)
      : r_states(num_states), g_states(num_orbitals), x(NULL) {
    for (size_t i = 0; i < r_states; i++) {
      x.push_back(zero_functions<double, 3>(world, g_states, true));
    }
    x[0][0].world().gop.fence();
  }
  // Conversion from  Constructor
  explicit response_space(const response_matrix& x)
      : r_states(x.size()), g_states(x[0].size()), x(NULL) {

    this->x = x;
    x[0][0].world().gop.fence();
  }

  // Copy constructor
  response_space(const response_space& other) {
    r_states = other.r_states;
    g_states = other.g_states;
    x.clear();
    //*this = ResponseVectors(other.x[0][0].world(), r_states, g_states);
    for (unsigned int i = 0; i < r_states; i++) {
      x.push_back(other.x[i]);
    }
  }

  response_space& operator=(const response_space& other) {
    //
    if (this != &other) {
    }
    r_states = other.r_states;
    g_states = other.g_states;
    x.clear();
    for (unsigned int i = 0; i < r_states; i++) {
      // this->x.push_back(madness::copy(other[0][0].world(), other[i]));
      x.push_back(other.x[i]);
    }
    return *this;
  }

  // Determines if two ResponseFunctions are the same size
  bool same_size(const response_space& rf_copy) const {
    return (r_states == rf_copy.r_states && g_states == rf_copy.g_states);
  }

  // 1D accessor for x
  // std::vector<Function<double, 3>> &operator[](int64_t i) { return x[i]; }
  vector_real_function_3d& operator[](int64_t i) { return x.at(i); }

  const vector_real_function_3d& operator[](int64_t i) const { return x.at(i); }

  // KAIN must have this
  // element wise addition.  we add each vector separatly
  response_space operator+(const response_space& b) {
    MADNESS_ASSERT(same_size(b));

    response_space result(x[0][0].world(), r_states, g_states);

    for (unsigned int i = 0; i < r_states; i++) {
      result[i] = x[i] + b[i];
    }

    result[0][0].world().gop.fence();
    return result;
  }

  friend response_space operator-(const response_space& a,
                                   const response_space& b) {
    MADNESS_ASSERT(a.same_size(b));

    response_space result(a[0][0].world(), a.r_states, a.g_states);

    for (unsigned int i = 0; i < a.r_states; i++) {
      result[i] = a[i] - b[i];
    }
    return result;
  }

  // KAIN must have this
  // Scaling by a constant
  response_space operator*(double a) {
    response_space result(*this);

    for (unsigned int i = 0; i < r_states; i++) {
      x[i] = x[i] * a;
    }

    x[0][0].world().gop.fence();
    return result;
  }

  // Scaling all internal functions by an external function
  // g[i][j] = x[i][j] * f
  friend response_space operator*(const response_space& X,
                                   const Function<double, 3>& f) {
    response_space result(X);

    for (unsigned int i = 0; i < X.r_states; i++) {
      // Using vmra.h funciton
      result[i] = X.x[i] * f;
    }

    f.world().gop.fence();
    return result;
  }
  // Scaling all internal functions by an external function
  // g[i][j] = x[i][j] * f
  friend response_space operator*(const Function<double, 3>& f,
                                   const response_space& X) {
    response_space result(X);

    for (unsigned int i = 0; i < X.r_states; i++) {
      // Using vmra.h funciton
      result[i] = X.x[i] * f;
    }

    f.world().gop.fence();
    return result;
  }
  response_space operator*(const Function<double, 3>& f) {
    response_space result(*this);

    for (unsigned int i = 0; i < r_states; i++) {
      // Using vmra.h funciton
      result[i] = mul(f.world(), f, x[i]);
    }

    f.world().gop.fence();
    return result;
  }

  // KAIN must have this
  response_space operator+=(const response_space b) {
    MADNESS_ASSERT(same_size(b));

    for (unsigned int i = 0; i < r_states; i++) {
      x[i] += b[i];
    }

    return *this;
  }

  // Returns a deep copy
  response_space copy() const {
    response_space result(x[0][0].world(), r_states, g_states);

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
inline double inner(response_space& a, response_space& b) {
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
