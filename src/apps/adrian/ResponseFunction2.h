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

typedef std::vector<vector_real_function_3d> response_matrix;

struct response_space {
  // Member variables
 public:
  size_t num_states;    // Num. of resp. states
  size_t num_orbitals;  // Num. of ground states
  response_matrix x;

  // Member functions
 public:
  // Default constructor
  response_space() : num_states(0), num_orbitals(0), x() {}

  // Copy constructor
  response_space(const response_space& y)
      : num_states(y.size()), num_orbitals(y.vec_size()), x() {
    x = y.x;
  }
  // Initializes functions to zero
  // m = number of response states
  // n = number of ground state orbitals

  // Zero Constructor
  response_space(World& world, size_t num_states, size_t num_orbitals)
      : num_states(num_states), num_orbitals(num_orbitals), x() {
    for (size_t i = 0; i < num_states; i++) {
      this->x.emplace_back(
          zero_functions<double, 3>(world, num_orbitals, true));
    }
    x[0][0].world().gop.fence();
  }
  // Conversion from  Constructor
  explicit response_space(const response_matrix& x)
      : num_states(x.size()), num_orbitals(x[0].size()), x() {
    this->x = x;
    x[0][0].world().gop.fence();
  }
  // Determines if two ResponseFunctions are the same size
  friend bool same_size(const response_space& x, const response_space& y) {
    return ((x.size() == y.size()) && (x.vec_size() == y.vec_size()));
  }
  // assignment
  response_space& operator=(const response_space& y) {
    //
    if (this != &y) {             // is it the same object?
      if (same_size(*this, y)) {  // is it same size?
        this->x = y.x;
        /*
        for (size_t b = 0; b < x.size(); b++) {
          (*this)[b] = x[b];
        }
        */
      } else {                         // if not the same size
        this->~response_space();       // deconstruct response_space
        new (this) response_space(y);  //  call copy constructor
      }
    }
    return *this;  // shallow copy
  }

  // 1D accessor for x
  // std::vector<Function<double, 3>> &operator[](int64_t i) { return x[i]; }
  vector_real_function_3d& operator[](int64_t i) { return x.at(i); }

  const vector_real_function_3d& operator[](int64_t i) const { return x.at(i); }

  // KAIN must have this
  // element wise addition.  we add each vector separatly
  // addition c = this.x+b
  // we need a new function
  response_space operator+(const response_space& rhs_y) {
    MADNESS_ASSERT(same_size(*this, rhs_y));  // assert that same size

    World& world = this->x[0][0].world();

    response_space result(
        world, num_states, num_orbitals);  // create zero_functions

    for (unsigned int i = 0; i < num_states; i++) {
      result[i] = add(world, x[i], rhs_y[i]);
    }
    return result;
  }
  friend response_space operator+(const response_space& a,
                                  const response_space& b) {
    MADNESS_ASSERT(same_size(a, b));

    World& world = a.x[0][0].world();
    response_space result(
        world, a.num_states, a.num_orbitals);  // create zero_functions

    for (unsigned int i = 0; i < a.num_states; i++) {
      result[i] = add(world, a[i], b[i]);
    }
    return result;
  }
  response_space operator-(const response_space& rhs_y) {
    MADNESS_ASSERT(same_size(*this, rhs_y));  // assert that same size

    World& world = this->x[0][0].world();

    response_space result(
        world, num_states, num_orbitals);  // create zero_functions

    for (unsigned int i = 0; i < num_states; i++) {
      result[i] = sub(world, x[i], rhs_y[i]);
    }
    return result;
  }

  friend response_space operator-(const response_space& a,
                                  const response_space& b) {
    MADNESS_ASSERT(same_size(a, b));

    World& world = a.x[0][0].world();
    response_space result(
        world, a.num_states, a.num_orbitals);  // create zero_functions

    for (unsigned int i = 0; i < a.num_states; i++) {
      result[i] = sub(world, a[i], b[i]);
    }
    return result;
  }

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
  friend response_space operator*(response_space y, double a) {
    World& world = y[0][0].world();
    response_space result = y.copy();  // deep copy

    for (unsigned int i = 0; i < y.num_states; i++) {
      madness::scale(world, result[i], a);
    }

    return result;
  }
  friend response_space operator*(double a, response_space y) {
    World& world = y[0][0].world();
    response_space result = y.copy();  // deep copy

    for (unsigned int i = 0; i < y.num_states; i++) {
      madness::scale(world, result[i], a);
    }

    return result;
  }
  response_space operator*=(double a) {
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
    World& world = a.x[0][0].world();
    response_space result(
        world, a.num_states, a.num_orbitals);  // create zero_functions

    for (unsigned int i = 0; i < a.num_states; i++) {
      // Using vmra.h funciton
      result[i] = mul(f.world(), f, a[i]);
    }

    f.world().gop.fence();
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

    f.world().gop.fence();
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
  response_space operator+=(const response_space b) {
    MADNESS_ASSERT(same_size(*this, b));
    World& world = x[0][0].world();
    for (unsigned int i = 0; i < num_states; i++) {
      add(world, this->x[i], b[i]);
    }
    return *this;
  }

  // Returns a deep copy
  response_space copy() const {
    response_space result(x[0][0].world(), num_states, num_orbitals);

    for (unsigned int i = 0; i < num_states; i++) {
      result.x[i] = madness::copy(x[0][0].world(), x[i]);
    }
    x[0][0].world().gop.fence();

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
  size_t vec_size() const { return num_orbitals; }

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
    x[0][0].world().gop.fence();
  }
  friend bool operator==(const response_space& x, const response_space& y) {
    if (!same_size(x, y)) return false;
    for (size_t b = 0; b < x.size(); ++b) {
      for (size_t k = 0; b < x.vec_size(); ++k) {
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
  // Zero Constructor
  X_space(World& world, size_t num_states, size_t num_orbitals)
      : num_states(num_states),
        num_orbitals(num_orbitals),
        X(world, num_states, num_orbitals),
        Y(world, num_states, num_orbitals) {}
  // assignment
  X_space& operator=(const X_space& B) {
    if (this != &B) {             // is it the same object?
      if (same_size(*this, B)) {  // is it same size?
        this->X = B.X;
        this->Y = B.Y;
      } else {                  // if not the same size
        this->~X_space();       // deconstruct response_space
        new (this) X_space(B);  //  call copy constructor
      }
    }
    return *this;  // shallow copy
  }

  X_space operator+(const X_space B) {
    MADNESS_ASSERT(same_size(*this, B));
    World& world = this->X[0][0].world();
    X_space result(world, num_states, num_orbitals);
    result.X = X + B.X;
    result.Y = Y + B.Y;
    return result;
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

}  // End namespace madness
#endif  // SRC_APPS_ADRIAN_RESPONSEFUNCTION2_H_

// Deuces
