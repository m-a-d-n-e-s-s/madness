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
/**
 * @brief response vector class holds a single x or y
 * 
 */
struct response_vector {
 private:
  size_t num_orbitals;
  vector_real_function_3d x;  // n response_functions

 public:
  // default constructor
  response_vector() : num_orbitals(size_t(0)), x() {}
  // copy constructor
  response_vector(const response_vector& x);
  // zero function constructor
  explicit response_vector(World& world, size_t num_orbs)
      : num_orbitals(num_orbs) {
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
  size_t num_states;    // Num. of resp. states
  size_t num_orbitals;  // Num. of ground states
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
   * we copy madness functions therefore we are copying pointers to function implementations
   * @param y 
   */
  response_space(const response_space& y)
      : num_states(y.size()), num_orbitals(y.size_orbitals()), x(y.x) {}
  // assignment
  response_space& operator=(const response_space& y) {
    //
    if (this != &y) {  // is it the same object?
      this->num_states = y.size();
      this->num_orbitals = y.size_orbitals();
      this->x = y.x;
    }
    return *this;  //
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
      this->x.emplace_back(
          zero_functions<double, 3>(world, num_orbitals, true));
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
    MADNESS_ASSERT(same_size(*this, rhs_y));  // assert that same size

    World& world = this->x[0][0].world();

    response_space result(
        world, num_states, num_orbitals);  // create zero_functions

    for (unsigned int i = 0; i < num_states; i++) {
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
    MADNESS_ASSERT(same_size(*this, rhs_y));  // assert that same size

    World& world = this->x[0][0].world();

    response_space result(
        world, num_states, num_orbitals);  // create zero_functions

    for (unsigned int i = 0; i < num_states; i++) {
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

    for (unsigned int i = 0; i < num_states; i++) {
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
    response_space result = y.copy();  // deep copy

    for (unsigned int i = 0; i < y.num_states; i++) {
      madness::scale(world, result[i], a);
    }

    return result;
  }
  friend response_space operator*(double a, response_space y) {
    World& world = y.x.at(0).at(0).world();
    response_space result = y.copy();  // deep copy

    for (unsigned int i = 0; i < y.num_states; i++) {
      madness::scale(world, result[i], a);
    }

    return result;
  }
  response_space& operator*=(double a) {
    World& world = this->x[0][0].world();

    for (unsigned int i = 0; i < num_states; i++) {
      madness::scale(world, this->x[i], a);
    }

    return *this;
  }

  // Scaling all internal functions by an external function
  // g[i][j] = x[i][j] * f
  friend response_space operator*(const response_space& a,
                                  const Function<double, 3>& f) {
    World& world = a.x.at(0).at(0).world();
    response_space result(
        world, a.num_states, a.num_orbitals);  // create zero_functions

    for (unsigned int i = 0; i < a.num_states; i++) {
      // Using vmra.h funciton
      result[i] = mul(f.world(), f, a[i]);
    }

    return result;
  }
  // Scaling all internal functions by an external function
  // g[i][j] = x[i][j] * f
  friend response_space operator*(const Function<double, 3>& f,
                                  const response_space& a) {
    // commutative property
    return a * f;
  }
  response_space operator*(const Function<double, 3>& f) {
    World& world = x[0][0].world();
    response_space result(
        world, num_states, num_orbitals);  // create zero_functions

    for (size_t i = 0; i < num_states; i++) {
      // Using vmra.h funciton
      result[i] = mul(f.world(), f, x[i]);
    }

    return result;
  }

  friend response_space operator*(const response_space& a,
                                  const Tensor<double>& b) {
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
    /*
    for (size_t i = 0; i < num_states; i++) {
      for (size_t j = 0; j < num_orbitals; j++) {
        this->x[i][j] += b[i][j];
      }
    }
    */
for (unsigned int i = 0; i < num_states; i++) {
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
    } else {  // g_states == 0 (empty vector)
      num_orbitals = f.size();
    }
  }
  void pop_back() {
    MADNESS_ASSERT(num_states >= 1);
    x.pop_back();
    num_states--;

    // Be smart with g_states
    if (num_states == 0) {  // removed last item
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
    for (unsigned int k = 0; k < num_states; k++) {
      x[k] = zero_functions<double, 3>(x[0][0].world(), num_orbitals);
    }
  }

  void compress_rf() {
    for (unsigned int k = 0; k < num_states; k++) {
      compress(x[0][0].world(), x[k], true);
    }
  }

  void reconstruct_rf() {
    for (unsigned int k = 0; k < num_states; k++) {
      reconstruct(x[0][0].world(), x[k], true);
    }
  }

  void truncate_rf() {
    for (unsigned int k = 0; k < num_states; k++) {
      truncate(x[0][0].world(), x[k], FunctionDefaults<3>::get_thresh(), true);
    }
  }
  void truncate_rf(double tol) {
    for (unsigned int k = 0; k < num_states; k++) {
      truncate(x[0][0].world(), x[k], tol, true);
    }
  }

  // Returns norms of each state
  Tensor<double> norm2() {
    Tensor<double> answer(num_states);
    for (unsigned int i = 0; i < num_states; i++) {
      answer(i) = madness::norm2(x[0][0].world(), x[i]);
    }
    return answer;
  }

  // Scales each state (read: entire row) by corresponding vector element
  //     new[i] = old[i] * mat[i]
  void scale(Tensor<double>& mat) {
    for (unsigned int i = 0; i < num_states; i++)
      madness::scale(x[0][0].world(), x[i], mat[i], false);
    // x[i] = x[i] * mat[i];
  }
  friend bool operator==(const response_space& x, const response_space& y) {
    if (!same_size(x, y)) return false;
    for (size_t b = 0; b < x.size(); ++b) {
      for (size_t k = 0; b < x.size_orbitals(); ++k) {
        if ((x[b][k] - y[b][k]).norm2() >
            FunctionDefaults<3>::get_thresh())  // this may be strict
          return false;
      }
    }
    return true;
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
struct X_space {
 private:
  size_t num_states;    // Num. of resp. states
  size_t num_orbitals;  // Num. of ground states
 public:
  response_space X, Y;

 public:
  // default constructor
  X_space() : num_states(0), num_orbitals(0), X(), Y() {}
  // Copy constructor
  X_space(const X_space& A)
      : num_states(size_states(A)),
        num_orbitals(size_orbitals(A)),
        X(A.X),
        Y(A.Y) {}
  // assignment
  X_space& operator=(const X_space& B) {
    if (this != &B) {  // is it the same object?
      this->num_states = size_states(B);
      this->num_orbitals = size_orbitals(B);
      this->X = B.X;
      this->Y = B.Y;
    }
    return *this;  // shallow copy
  }
  // access...reshaper
  /*
  X_space operator[](size_t n) {
    MADNESS_ASSERT(n < num_states);
    World& world = X[0][0].world();
    X_space newX(world, size_t(1), num_orbitals);
    response_space single_X(world, 1, num_states);
    response_space single_Y(world, 1, num_states);
    single_X.push_back(X[n]);
    single_Y.push_back(Y[n]);
    newX.X = response_space(single_X);
    newX.Y = response_space(single_Y);
    return newX;
  }
  */
  // Zero Constructor
  X_space(World& world, size_t num_states, size_t num_orbitals)
      : num_states(num_states),
        num_orbitals(num_orbitals),
        X(world, num_states, num_orbitals),
        Y(world, num_states, num_orbitals) {}
  // explicit constructor from 2 resonse_space
  explicit X_space(response_space& X, response_space& Y) {
    MADNESS_ASSERT(X.size() == Y.size());
    MADNESS_ASSERT(X[0].size() == Y[0].size());
    this->num_states = X.size();
    this->num_orbitals = X[0].size();
    this->X = X.copy();
    this->Y = Y.copy();
  }

  X_space operator+(const X_space B) {
    MADNESS_ASSERT(same_size(*this, B));
    World& world = this->X[0][0].world();
    X_space result(world, num_states, num_orbitals);
    result.X = X + B.X;
    result.Y = Y + B.Y;
    return result;
  }
  X_space& operator+=(const X_space B) {
    MADNESS_ASSERT(same_size(*this, B));
    this->X += B.X;
    this->Y += B.Y;
    return *this;
  }

  friend X_space operator+(const X_space& A, const X_space& B) {
    MADNESS_ASSERT(same_size(A, B));

    World& world = A.X[0][0].world();
    X_space result(
        world, A.num_states, A.num_orbitals);  // create zero_functions

    result.X = A.X + B.X;
    result.Y = A.Y + B.Y;
    return result;
  }

  X_space operator-(const X_space B) {
    MADNESS_ASSERT(same_size(*this, B));
    World& world = this->X[0][0].world();
    X_space result(world, num_states, num_orbitals);
    result.X = X - B.X;
    result.Y = Y - B.Y;
    return result;
  }

  friend X_space operator-(const X_space& A, const X_space& B) {
    MADNESS_ASSERT(same_size(A, B));

    World& world = A.X[0][0].world();
    X_space result(
        world, A.num_states, A.num_orbitals);  // create zero_functions

    result.X = A.X - B.X;
    result.Y = A.Y - B.Y;
    return result;
  }

  friend X_space operator*(const X_space& A, const double& b) {
    World& world = A.X[0][0].world();
    X_space result(
        world, A.num_states, A.num_orbitals);  // create zero_functions

    result.X = A.X * b;
    result.Y = A.Y * b;
    return result;
  }
  friend X_space operator*(const double& b, const X_space& A) {
    World& world = A.X[0][0].world();
    X_space result(
        world, A.num_states, A.num_orbitals);  // create zero_functions

    result.X = A.X * b;
    result.Y = A.Y * b;
    return result;
  }
  X_space operator*(const double& b) {
    this->X *= b;
    this->Y *= b;
    return *this;
  }

  friend X_space operator*(const X_space& A, const Function<double, 3>& f) {
    World& world = A.X[0][0].world();
    X_space result(
        world, A.num_states, A.num_orbitals);  // create zero_functions

    result.X = A.X * f;
    result.Y = A.Y * f;
    return result;
  }
  friend X_space operator*(const Function<double, 3>& f, const X_space& A) {
    World& world = A.X[0][0].world();
    X_space result(
        world, A.num_states, A.num_orbitals);  // create zero_functions

    result.X = A.X * f;
    result.Y = A.Y * f;
    return result;
  }

  friend X_space operator*(const X_space& A, const Tensor<double>& b) {
    MADNESS_ASSERT(size_states(A) > 0);
    MADNESS_ASSERT(size_orbitals(A) > 0);

    World& world = A.X[0][0].world();
    X_space result(world, A.num_states, A.num_orbitals);
    result.X = A.X * b;
    result.Y = A.Y * b;

    return result;
  }
  inline friend Tensor<double> inner(X_space& A, X_space& B) {
    MADNESS_ASSERT(size_states(A) > 0);
    MADNESS_ASSERT(size_orbitals(A) > 0);
    MADNESS_ASSERT(same_size(A, B));
    Tensor<double> G(A.num_states, A.num_states);

    World& world = A.X[0][0].world();
    response_space Collapse(world, A.num_states, A.num_states);

    for (size_t i(0); i < A.num_states; i++) {
      for (size_t j(0); j < A.num_states; j++) {
        Collapse[i][j] =
            dot(world, A.X[i], B.X[j]) + dot(world, A.Y[i], B.Y[j]);
        G(i, j) = Collapse[i][j].trace();
      }
    }

    return G;
  }
  friend size_t size_states(const X_space& x) { return x.num_states; }
  friend size_t size_orbitals(const X_space& x) { return x.num_orbitals; }
  friend bool same_size(const X_space& A, const X_space& B) {
    return ((size_states(A) == size_states(B) &&
             size_orbitals(A) == size_orbitals(B)));
  }
};
// The default constructor for functions does not initialize them to nahy value,
// but the solver needs the functions initialized to zero for which we also need
// the world object.

struct X_vector : public X_space {
  X_vector(World& world, size_t num_orbitals)
      : X_space(world, size_t(1), num_orbitals) {}

  X_vector(X_space A, size_t b)
      : X_space(A.X[0][0].world(), size_t(1), size_orbitals(A)) {
    X[0] = A.X[b];
    Y[0] = A.Y[b];
    // this->X[0].assign(A.X[b].begin(), A.X[b].end());  // = single_X;
    // this->Y[0].assign(A.Y[b].begin(), A.Y[b].end());  // = single_X;
  }
  friend X_vector operator-(const X_vector& A, const X_vector& B) {
    MADNESS_ASSERT(same_size(A, B));

    World& world = A.X[0][0].world();
    X_vector result(world, size_orbitals(A));  // create zero_functions
    result.X = A.X - B.X;
    result.Y = A.Y - B.Y;
    return result;
  }
  friend X_vector operator*(const X_vector& A, const double& c) {
    World& world = A.X[0][0].world();
    X_vector result(world, size_orbitals(A));  // create zero_functions
    result.X = A.X * c;
    result.Y = A.Y * c;
    return result;
  }
  X_vector& operator+=(const X_vector& B) {
    MADNESS_ASSERT(same_size(*this, B));

    this->X += B.X;
    this->Y += B.Y;

    return *this;
  }
  inline friend double inner(X_vector& A, X_vector& B) {
    MADNESS_ASSERT(size_states(A) == 1);
    MADNESS_ASSERT(size_orbitals(A) > 0);
    MADNESS_ASSERT(same_size(A, B));

    World& world = A.X[0][0].world();

    real_function_3d density =
        dot(world, A.X[0], B.X[0]) + dot(world, A.Y[0], B.Y[0]);

    return density.trace();
  }
};
struct X_space_allocator {
  World& world;
  const size_t num_states;
  const size_t num_orbitals;
  X_space_allocator(World& world, size_t num_orbitals)
      : world(world), num_states(size_t(1)), num_orbitals(num_orbitals) {}
  // overloading the default constructor () operator
  X_vector operator()() { return X_vector(world, num_orbitals); }
  // Copy constructor

  X_space_allocator operator=(const X_space_allocator& other) {
    return X_space_allocator(world, other.num_orbitals);
  }
};

}  // End namespace madness
#endif  // SRC_APPS_ADRIAN_RESPONSEFUNCTION2_H_

// Deuces
