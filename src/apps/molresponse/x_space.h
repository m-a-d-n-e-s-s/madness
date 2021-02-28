
#ifndef SRC_APPS_molresponse_x_space_H_
#define SRC_APPS_molresponse_x_space_H_

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <vector>


namespace madness {
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
}

#endif  // SRC_APPS_molresponse_x_space_H_