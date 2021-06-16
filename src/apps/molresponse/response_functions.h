/* Copyright 2021 Adrian Hurtado
 *   Small class to hold response functions and to interact with KAIN solver.
 */

#ifndef SRC_APPS_MOLRESPONSE_RESPONSE_FUNCTIONS_H_
#define SRC_APPS_MOLRESPONSE_RESPONSE_FUNCTIONS_H_

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <vector>

// real_function_3d == Function<double,3>
namespace madness {
/**
 * @brief response vector class holds a single x or y
 *
 */
struct response_vector {
private:
  size_t num_orbitals;
  vector_real_function_3d x; // n response_functions

public:
  // default constructor
  response_vector() : num_orbitals(size_t(0)), x() {}
  // copy constructor
  response_vector(const response_vector& x);
  // zero function constructor
  explicit response_vector(World& world, size_t num_orbs) : num_orbitals(num_orbs) {
    x = zero_functions<double, 3>(world, num_orbitals);
  }
  ~response_vector() {}
  response_vector& operator=(const response_vector& x);
  friend size_t size(const response_vector& x) { return x.num_orbitals; }
  real_function_3d& operator[](size_t n) {
    assert(n < num_orbitals);
    return x.at(n);
  }
  const real_function_3d& operator[](size_t n) const {
    assert(n < num_orbitals);
    return x[n];
  }
};

/**
 * @brief response matrix holds response vectors for response state
 *
 */
struct response_space {
  // Member variables
  /**
   * @brief vector of vector of real 3d functions
   *
   */
  typedef std::vector<vector_real_function_3d> response_matrix;

public:
  size_t num_states;   // Num. of resp. states
  size_t num_orbitals; // Num. of ground states
  response_matrix x;

  // Member functions
public:
  /**
   * @brief default Construct a new response space object
   * num_states(0)
   * num_orbitals(0)
   * x() default constructor of std::vector
   */
  response_space() : num_states(0), num_orbitals(0), x() {}

  // Copy constructor
  /**
   * @brief copy construct a new response space object
   * we are using copying defined by std:vector
   * we copy madness functions therefore we are copying pointers to function
   * implementations
   * @param y
   */
  response_space(const response_space& y)
      : num_states(y.size()), num_orbitals(y.size_orbitals()), x(y.x) {}
  // assignment
  response_space& operator=(const response_space& y) {
    //
    if (this != &y) { // is it the same object?
      this->num_states = y.size();
      this->num_orbitals = y.size_orbitals();
      this->x = y.x;
    }
    return *this; //
  }
  // Initializes functions to zero
  // m = number of response states
  // n = number of ground state orbitals

  // Zero Constructor constructs m vectors
  /**
   * @brief Construct a new response space with zero functions
   *
   * @param world
   * @param num_states
   * @param num_orbitals
   */
  response_space(World& world, size_t num_states, size_t num_orbitals)
      : num_states(num_states), num_orbitals(num_orbitals), x() {
    for (size_t i = 0; i < num_states; i++) {
      this->x.emplace_back(zero_functions<double, 3>(world, num_orbitals, true));
    }
  }
  // Conversion from respones_matrix
  /**
   * @brief Construct a new response space object from vector of functions
   *
   * @param x
   */
  explicit response_space(const response_matrix& x)
      : num_states(x.size()), num_orbitals(x[0].size()), x(x) {}
  // Determines if two ResponseFunctions are the same size
  friend bool same_size(const response_space& x, const response_space& y) {
    return ((x.size() == y.size()) && (x.size_orbitals() == y.size_orbitals()));
  }

  // 1D accessor for x
  // std::vector<Function<double, 3>> &operator[](int64_t i) { return x[i]; }
  /**
   * @brief access vector of functions with std::vector.at()
   *
   * @param i
   * @return vector_real_function_3d&
   */
  vector_real_function_3d& operator[](size_t i) { return x.at(i); }
  /**
   * @brief access vector of functions const
   *
   * @param i
   * @return const vector_real_function_3d&
   */
  const vector_real_function_3d& operator[](size_t i) const { return x.at(i); }

  // KAIN must have this
  // element wise addition.  we add each vector separatly
  // addition c = this.x+b
  // we need a new function
  /**
   * @brief elementwise addition of response_space
   *
   * @param rhs_y
   * @return response_space
   */
  response_space operator+(const response_space& rhs_y) const {
    MADNESS_ASSERT(size() > 0);
    MADNESS_ASSERT(same_size(*this, rhs_y)); // assert that same size

    World& world = this->x[0][0].world();

    response_space result(world, num_states, num_orbitals); // create zero_functions

    for (size_t i = 0; i < num_states; i++) {
      result[i] = add(world, x[i], rhs_y[i]);
    }
    return result;
  }
  /*
  friend response_space operator+(const response_space& a,
                                  const response_space& b) {
    return a.operator+(b);
  }
  */

  /**
   * @brief elementwise subtraction of response space
   *
   * @param rhs_y
   * @return response_space
   */
  response_space operator-(const response_space& rhs_y) const {
    MADNESS_ASSERT(size() > 0);
    MADNESS_ASSERT(same_size(*this, rhs_y)); // assert that same size

    World& world = this->x[0][0].world();

    response_space result(world, num_states, num_orbitals); // create zero_functions

    for (size_t i = 0; i < num_states; i++) {
      result[i] = sub(world, x[i], rhs_y[i]);
    }
    world.gop.fence();
    return result;
  }
  /*
    friend response_space operator-(const response_space& a,
                                    const response_space& b) {
      return a.operator-(b);
    }
    */

  // KAIN must have this
  // Scaling by a constant
  // return response_space scaled
  /*
  response_space operator*(double a) {
    World& world = x[0][0].world();
    // response_space result(*this);// this a problem
    response_space result = this->copy();  // deep copy

    for (size_t i = 0; i < num_states; i++) {
      madness::scale(world, result[i], a);
    }

    return result;
  }
  */

  /**
   * @brief multiplication by scalar
   *
   * @param y
   * @param a
   * @return response_space
   */
  friend response_space operator*(response_space y, double a) {
    World& world = y.x.at(0).at(0).world();
    response_space result = y.copy(); // deep copy

    for (unsigned int i = 0; i < y.num_states; i++) {
      madness::scale(world, result[i], a);
    }

    return result;
  }
  friend response_space operator*(double a, response_space y) {
    World& world = y.x.at(0).at(0).world();
    response_space result = y.copy(); // deep copy

    for (unsigned int i = 0; i < y.num_states; i++) {
      madness::scale(world, result[i], a);
    }

    return result;
  }
  response_space& operator*=(double a) {
    World& world = this->x[0][0].world();

    for (size_t i = 0; i < num_states; i++) {
      madness::scale(world, this->x[i], a);
    }

    return *this;
  }

  // Scaling all internal functions by an external function
  // g[i][j] = x[i][j] * f
  friend response_space operator*(const response_space& a, const Function<double, 3>& f) {
    World& world = a.x.at(0).at(0).world();
    response_space result(world, a.num_states, a.num_orbitals); // create zero_functions

    for (unsigned int i = 0; i < a.num_states; i++) {
      // Using vmra.h funciton
      result[i] = mul(f.world(), f, a[i]);
    }

    return result;
  }
  // Scaling all internal functions by an external function
  // g[i][j] = x[i][j] * f
  friend response_space operator*(const Function<double, 3>& f, const response_space& a) {
    // commutative property
    return a * f;
  }
  response_space operator*(const Function<double, 3>& f) {
    World& world = x[0][0].world();
    response_space result(world, num_states, num_orbitals); // create zero_functions

    for (size_t i = 0; i < num_states; i++) {
      // Using vmra.h funciton
      result[i] = mul(f.world(), f, x[i]);
    }

    return result;
  }

  friend response_space operator*(const response_space& a, const Tensor<double>& b) {
    MADNESS_ASSERT(a.size() > 0);
    MADNESS_ASSERT(a[0].size() > 0);
    World& world = a[0][0].world();
    response_space result(world, a.num_states, a.num_orbitals);

    for (unsigned int i = 0; i < a.size(); i++) {
      // Using vmra.h definitions
      result[i] = transform(world, a[i], b, false);
    }

    return result;
  }
  // KAIN must have this
  response_space& operator+=(const response_space b) {
    MADNESS_ASSERT(same_size(*this, b));
    World& world = x[0][0].world();
    for (size_t i = 0; i < num_states; i++) {
      this->x[i] = add(world, this->x[i], b[i]);
    }

    return *this;
  }

  // Returns a deep copy
  response_space copy() const {
    response_space result(x[0][0].world(), num_states, num_orbitals);

    for (size_t i = 0; i < num_states; i++) {
      result.x[i] = madness::copy(x[0][0].world(), x[i]);
    }

    return result;
  }

  // Mimicking std::vector with these 4
  void push_back(const vector_real_function_3d& f) {
    x.push_back(f);
    num_states++;

    // Be smart with g_states
    if (num_orbitals > 0) {
      MADNESS_ASSERT(num_orbitals = f.size());
    } else { // g_states == 0 (empty vector)
      num_orbitals = f.size();
    }
  }
  void pop_back() {
    MADNESS_ASSERT(num_states >= 1);
    x.pop_back();
    num_states--;

    // Be smart with g_states
    if (num_states == 0) { // removed last item
      num_orbitals = 0;
    }
  }
  void clear() {
    x.clear();
    num_states = 0;
    num_orbitals = 0;
  }
  size_t size() const { return num_states; }
  size_t size_orbitals() const { return num_orbitals; }

  // Mimicing standard madness calls with these 3
  void zero() {
    for (size_t k = 0; k < num_states; k++) {
      x[k] = zero_functions<double, 3>(x[0][0].world(), num_orbitals);
    }
  }

  void compress_rf() {
    for (size_t k = 0; k < num_states; k++) {
      compress(x[0][0].world(), x[k], true);
    }
  }

  void reconstruct_rf() {
    for (size_t k = 0; k < num_states; k++) {
      reconstruct(x[0][0].world(), x[k], true);
    }
  }

  void truncate_rf() {
    for (size_t k = 0; k < num_states; k++) {
      truncate(x[0][0].world(), x[k], FunctionDefaults<3>::get_thresh(), true);
    }
  }
  void truncate_rf(double tol) {
    for (size_t k = 0; k < num_states; k++) {
      truncate(x[0][0].world(), x[k], tol, true);
    }
  }

  // Returns norms of each state
  Tensor<double> norm2() {
    Tensor<double> answer(num_states);
    for (size_t i = 0; i < num_states; i++) {
      answer(i) = madness::norm2(x[0][0].world(), x[i]);
    }
    return answer;
  }

  // Scales each state (read: entire row) by corresponding vector element
  //     new[i] = old[i] * mat[i]
  void scale(Tensor<double>& mat) {
    for (size_t i = 0; i < num_states; i++)
      madness::scale(x[0][0].world(), x[i], mat[i], false);
    // x[i] = x[i] * mat[i];
  }
  friend bool operator==(const response_space& x, const response_space& y) {
    if (!same_size(x, y))
      return false;
    for (size_t b = 0; b < x.size(); ++b) {
      for (size_t k = 0; b < x.size_orbitals(); ++k) {
        if ((x[b][k] - y[b][k]).norm2() > FunctionDefaults<3>::get_thresh()) // this may be strict
          return false;
      }
    }
    return true;
  }
  friend Tensor<double> response_space_inner(response_space& a, response_space& b) {
    MADNESS_ASSERT(a.size() > 0);
    MADNESS_ASSERT(a.size() == b.size());
    MADNESS_ASSERT(a[0].size() > 0);
    MADNESS_ASSERT(a[0].size() == b[0].size());

    World& world = a[0][0].world();

    size_t dim_1 = a.size();
    size_t dim_2 = b[0].size();
    // Need to take transpose of each input ResponseFunction
    response_space aT(world, dim_2, dim_1);
    response_space bT(world, dim_2, dim_1);
    for (size_t i = 0; i < dim_1; i++) {
      for (size_t j = 0; j < dim_2; j++) {
        aT[j][i] = a[i][j];
        bT[j][i] = b[i][j];
      }
    }
    // Container for result
    Tensor<double> result(dim_1, dim_1);
    for (size_t p = 0; p < dim_2; p++) {
      result += matrix_inner(world, aT[p], bT[p]);
    }
    return result;
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
// Final piece for KAIN

} // End namespace madness
#endif // SRC_APPS_MOLRESPONSE_RESPONSE_FUNCTIONS_H_

// Deuces
