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
  size_t num_vectors;  // Num. of resp. states
  size_t num_orbitals; // Num. of ground states
  response_matrix x;

  // Member functions
public:
  // Default constructor
  response_space() : num_vectors(0), num_orbitals(0), x(NULL) {}

  // Copy constructor
  response_space(const response_space& other) {
    num_vectors = other.num_vectors;
    num_orbitals = other.num_orbitals;
    x.clear();
    //*this = ResponseVectors(other.x[0][0].world(), r_states, g_states);
    for (unsigned int i = 0; i < num_vectors; i++) {
      x.push_back(other.x[i]);
    }
  }
  // Initializes functions to zero
  // m = number of response states
  // n = number of ground state orbitals

  // Zero Constructor
  response_space(World& world, size_t num_states, size_t num_orbitals)
      : num_vectors(num_states), num_orbitals(num_orbitals), x(NULL) {
    for (size_t i = 0; i < num_vectors; i++) {
      x.push_back(zero_functions<double, 3>(world, num_orbitals, true));
    }
    x[0][0].world().gop.fence();
  }
  // Conversion from  Constructor
  explicit response_space(const response_matrix& x)
      : num_vectors(x.size()), num_orbitals(x[0].size()), x(NULL) {

    this->x = x;
    x[0][0].world().gop.fence();
  }
  friend size_t space_size(const response_space& x) { return x.num_vectors; }
  friend size_t vector_size(const response_space& x) { return x.num_orbitals; }
  // Determines if two ResponseFunctions are the same size
  friend bool same_size(const response_space& x, const response_space& y) {
    return ((space_size(x) == space_size(x)) &&
            (vector_size(y) == vector_size(x)));
  }
  // assignment
  response_space& operator=(const response_space& x) {
    //
    if (this != &x) {
      if (same_size(*this, x)) {
      } else {
        this->~response_space();      // deconstruct
        new (this) response_space(x); // copy constructor
      }
    }
    return *this;
  }

  // 1D accessor for x
  // std::vector<Function<double, 3>> &operator[](int64_t i) { return x[i]; }
  vector_real_function_3d& operator[](int64_t i) { return x.at(i); }

  const vector_real_function_3d& operator[](int64_t i) const { return x.at(i); }

  // KAIN must have this
  // element wise addition.  we add each vector separatly
  response_space operator+(const response_space& b) {
    MADNESS_ASSERT(same_size(*this, b));

    response_space result(x[0][0].world(), num_vectors, num_orbitals);

    for (unsigned int i = 0; i < num_vectors; i++) {
      result[i] = x[i] + b[i];
    }

    result[0][0].world().gop.fence();
    return result;
  }

  friend response_space operator-(const response_space& a,
                                  const response_space& b) {
    MADNESS_ASSERT(same_size(a, b));

    response_space result(a[0][0].world(), a.num_vectors, a.num_orbitals);

    for (unsigned int i = 0; i < a.num_vectors; i++) {
      result[i] = a[i] - b[i];
    }
    return result;
  }

  // KAIN must have this
  // Scaling by a constant
  response_space operator*(double a) {
    response_space result(*this);

    for (unsigned int i = 0; i < num_vectors; i++) {
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

    for (unsigned int i = 0; i < X.num_vectors; i++) {
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

    for (unsigned int i = 0; i < X.num_vectors; i++) {
      // Using vmra.h funciton
      result[i] = X.x[i] * f;
    }

    f.world().gop.fence();
    return result;
  }
  response_space operator*(const Function<double, 3>& f) {
    response_space result(*this);

    for (unsigned int i = 0; i < num_vectors; i++) {
      // Using vmra.h funciton
      result[i] = mul(f.world(), f, x[i]);
    }

    f.world().gop.fence();
    return result;
  }

  // KAIN must have this
  response_space operator+=(const response_space b) {
    MADNESS_ASSERT(same_size(*this, b));

    for (unsigned int i = 0; i < num_vectors; i++) {
      x[i] += b[i];
    }

    return *this;
  }

  // Returns a deep copy
  response_space copy() const {
    response_space result(x[0][0].world(), num_vectors, num_orbitals);

    for (unsigned int i = 0; i < num_vectors; i++) {
      result.x[i] = madness::copy(x[0][0].world(), x[i]);
    }
    x[0][0].world().gop.fence();

    return result;
  }

  // Mimicking std::vector with these 4
  void push_back(const vector_real_function_3d& f) {
    x.push_back(f);
    num_vectors++;

    // Be smart with g_states
    if (num_orbitals > 0) {
      MADNESS_ASSERT(num_orbitals = f.size());
    } else { // g_states == 0 (empty vector)
      num_orbitals = f.size();
    }
  }
  void pop_back() {
    MADNESS_ASSERT(num_vectors >= 1);
    x.pop_back();
    num_vectors--;

    // Be smart with g_states
    if (num_vectors == 0) { // removed last item
      num_orbitals = 0;
    }
  }
  void clear() {
    x.clear();
    num_vectors = 0;
    num_orbitals = 0;
  }
  unsigned int size() const { return num_vectors; }

  // Mimicing standard madness calls with these 3
  void zero() {
    for (unsigned int k = 0; k < num_vectors; k++) {
      x[k] = zero_functions<double, 3>(x[0][0].world(), num_orbitals);
    }
  }

  void compress_rf() {
    for (unsigned int k = 0; k < num_vectors; k++) {
      compress(x[0][0].world(), x[k], true);
    }
  }

  void reconstruct_rf() {
    for (unsigned int k = 0; k < num_vectors; k++) {
      reconstruct(x[0][0].world(), x[k], true);
    }
  }

  void truncate_rf() {
    for (unsigned int k = 0; k < num_vectors; k++) {
      truncate(x[0][0].world(), x[k], FunctionDefaults<3>::get_thresh(), true);
    }
  }
  void truncate_rf(double tol) {
    for (unsigned int k = 0; k < num_vectors; k++) {
      truncate(x[0][0].world(), x[k], tol, true);
    }
  }

  // Returns norms of each state
  Tensor<double> norm2() {
    Tensor<double> answer(num_vectors);
    for (unsigned int i = 0; i < num_vectors; i++) {
      answer(i) = madness::norm2(x[0][0].world(), x[i]);
    }
    return answer;
  }

  // Scales each state (read: entire row) by corresponding vector element
  //     new[i] = old[i] * mat[i]
  void scale(Tensor<double>& mat) {
    for (unsigned int i = 0; i < num_vectors; i++)
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
