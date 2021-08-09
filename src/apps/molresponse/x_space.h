// Copyright 2021 Adrian Hurtado
#ifndef SRC_APPS_MOLRESPONSE_X_SPACE_H_
#define SRC_APPS_MOLRESPONSE_X_SPACE_H_

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <vector>

#include "molresponse/response_functions.h"

namespace madness {
struct X_space {
 private:
  size_t n_states;    // Num. of resp. states
  size_t n_orbtials;  // Num. of ground states

 public:
  response_space X, Y;

 public:
  size_t num_states() { return n_states; }
  size_t num_orbitals() { return n_orbtials; }
  // default constructor
  X_space() : n_states(0), n_orbtials(0), X(), Y() {}
  // Copy constructor
  X_space(const X_space& A)
      : n_states(size_states(A)),
        n_orbtials(size_orbitals(A)),
        X(A.X),
        Y(A.Y) {}
  X_space copy() const {
    X_space copyX(X[0][0].world(), n_states, n_orbtials);
    copyX.X = X.copy();
    copyX.Y = Y.copy();
    return copyX;
  }
  /// Create a new copy of the function with different distribution and optional
  /// fence

  /// Works in either basis.  Different distributions imply
  /// asynchronous communication and the optional fence is
  /// collective.
  X_space copy(const std::shared_ptr<WorldDCPmapInterface<Key<3> > >& pmap,
               bool fence = false) const {
    X_space copyX(X[0][0].world(), n_states, n_orbtials);
    copyX.X = X.copy(pmap, fence);
    copyX.Y = Y.copy(pmap, fence);
    return copyX;
  }
  // assignment
  X_space& operator=(const X_space& B) {
    if (this != &B) {  // is it the same object?
      this->n_states = size_states(B);
      this->n_orbtials = size_orbitals(B);
      this->X = B.X;
      this->Y = B.Y;
    }
    return *this;  // shallow copy
  }
  // Zero Constructor
  X_space(World& world, size_t n_states, size_t n_orbtials)
      : n_states(n_states),
        n_orbtials(n_orbtials),
        X(world, n_states, n_orbtials),
        Y(world, n_states, n_orbtials) {}
  // explicit constructor from 2 resonse_space
  explicit X_space(response_space& X, response_space& Y) {
    MADNESS_ASSERT(X.size() == Y.size());
    MADNESS_ASSERT(X[0].size() == Y[0].size());
    this->n_states = X.size();
    this->n_orbtials = X[0].size();
    this->X = X.copy();
    this->Y = Y.copy();
  }
  void clear() {
    X.clear();
    Y.clear();
  }
  X_space operator+(const X_space B) {
    MADNESS_ASSERT(same_size(*this, B));
    World& world = this->X[0][0].world();
    X_space result(world, n_states, n_orbtials);
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
    X_space result(world, A.n_states, A.n_orbtials);  // create zero_functions

    result.X = A.X + B.X;
    result.Y = A.Y + B.Y;
    return result;
  }

  X_space operator-(const X_space B) {
    MADNESS_ASSERT(same_size(*this, B));
    World& world = this->X[0][0].world();
    X_space result(world, n_states, n_orbtials);
    result.X = X - B.X;
    result.Y = Y - B.Y;
    return result;
  }

  friend X_space operator-(const X_space& A, const X_space& B) {
    MADNESS_ASSERT(same_size(A, B));

    World& world = A.X[0][0].world();
    X_space result(world, A.n_states, A.n_orbtials);  // create zero_functions

    result.X = A.X - B.X;
    result.Y = A.Y - B.Y;
    return result;
  }

  friend X_space operator*(const X_space& A, const double& b) {
    World& world = A.X[0][0].world();
    X_space result(world, A.n_states, A.n_orbtials);  // create zero_functions

    result.X = A.X * b;
    result.Y = A.Y * b;
    return result;
  }
  friend X_space operator*(const double& b, const X_space& A) {
    World& world = A.X[0][0].world();
    X_space result(world, A.n_states, A.n_orbtials);  // create zero_functions

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
    X_space result(world, A.n_states, A.n_orbtials);  // create zero_functions

    result.X = A.X * f;
    result.Y = A.Y * f;
    return result;
  }
  friend X_space operator*(const Function<double, 3>& f, const X_space& A) {
    World& world = A.X[0][0].world();
    X_space result(world, A.n_states, A.n_orbtials);  // create zero_functions

    result.X = A.X * f;
    result.Y = A.Y * f;
    return result;
  }

  friend X_space operator*(const X_space& A, const Tensor<double>& b) {
    MADNESS_ASSERT(size_states(A) > 0);
    MADNESS_ASSERT(size_orbitals(A) > 0);

    World& world = A.X[0][0].world();
    X_space result(world, A.n_states, A.n_orbtials);
    result.X = A.X * b;
    result.Y = A.Y * b;

    return result;
  }
  inline friend Tensor<double> inner(X_space& A, X_space& B) {
    MADNESS_ASSERT(size_states(A) > 0);
    MADNESS_ASSERT(size_orbitals(A) > 0);
    MADNESS_ASSERT(same_size(A, B));
    Tensor<double> G(A.n_states, A.n_states);
    Tensor<double> G1(A.n_states, A.n_states);
    Tensor<double> G2(A.n_states, A.n_states);
    G1 = response_space_inner(A.X, B.X);
    G2 = response_space_inner(A.Y, B.Y);
    // TODO find a way to print seperate pieces based on a flag
    print("inner <AX|BX>");
    print(G1);
    print("inner <AY|BY>");
    print(G2);

    G = G1 + G2;
    return G;
  }

  void truncate() {
    X.truncate_rf();
    Y.truncate_rf();
  }

  void print_norm2() {
    for (size_t i = 0; i < num_states(); i++) {
      std::cout << "state " << i;
      std::cout << " X: ";

      for (size_t j = 0; j < num_orbitals(); j++) {
        std::cout << " " << X[i][j].norm2() << " ";
      }
      std::cout << " Y: ";
      for (size_t j = 0; j < num_orbitals(); j++) {
        std::cout << " " << Y[i][j].norm2() << " ";
      }
      std::cout << std::endl;
    }
  }

  friend size_t size_states(const X_space& x) { return x.n_states; }
  friend size_t size_orbitals(const X_space& x) { return x.n_orbtials; }
  friend bool same_size(const X_space& A, const X_space& B) {
    return ((size_states(A) == size_states(B) &&
             size_orbitals(A) == size_orbitals(B)));
  }
};
// The default constructor for functions does not initialize them to nahy value,
// but the solver needs the functions initialized to zero for which we also need
// the world object.

struct X_vector : public X_space {
  X_vector(World& world, size_t n_orbtials)
      : X_space(world, size_t(1), n_orbtials) {}

  X_vector(X_space A, size_t b)
      : X_space(A.X[0][0].world(), size_t(1), size_orbitals(A)) {
    X[0] = A.X[b];
    Y[0] = A.Y[b];
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
    return inner(world, A.X[0], B.X[0]).sum() +
           inner(world, A.Y[0], B.X[0]).sum();
  }
};
struct X_space_allocator {
  World& world;
  const size_t n_states;
  const size_t n_orbtials;
  X_space_allocator(World& world, size_t n_orbtials)
      : world(world), n_states(size_t(1)), n_orbtials(n_orbtials) {}
  // overloading the default constructor () operator
  X_vector operator()() { return X_vector(world, n_orbtials); }
  // Copy constructor

  X_space_allocator operator=(const X_space_allocator& other) {
    return X_space_allocator(world, other.n_orbtials);
  }
};
}  // namespace madness

#endif  // SRC_APPS_MOLRESPONSE_X_SPACE_H_
