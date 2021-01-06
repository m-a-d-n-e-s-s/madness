

#ifndef SRC_APPS_ADRIAN_RESPONSE_SPACE_H_
#define SRC_APPS_ADRIAN_RESPONSE_SPACE_H_

#include "adrian/response_vector.h"

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <vector>

namespace madness {

template <typename T> struct response_space {

  typedef Function<T, 3> function_T;
  typedef T value_type;
  typedef response_vector<T> response_vector_T;

  size_t num_orbitals;
  size_t num_vectors;
  std::vector<response_vector_T> V;

  // conversion from 2 vector functions to response_vector
  explicit response_space(const std::vector<response_vector_T>& W) {
    World world = W.x[0].world();
    num_vectors = W.size();
    num_orbitals = W[0].size();
    for (size_t b; b < num_vectors; b++) {
      V.push_back(response_vector(W[b]));
    }
  }

  // copy constructor
  explicit response_space(const response_space& other) {
    *this = response_space(other.V);
  }

  // default constructor
  response_space() {}
  // destructor
  ~response_space() {}
  // zero_constructor
  response_space(World& world, size_t num_vecs, size_t num_orbs)
      : num_vectors(num_vecs), num_orbitals(num_orbs) {
    for (size_t b; b < num_vectors; b++) {
      V.push_back(response_vector_T(world, num_orbitals));
    }
  }
  // assignment operator
  response_space& operator=(const response_space& W) {
    num_orbitals = W.num_orbitals;
    num_vectors = W.num_vectors;
    V = W.V;
    return *this;
  }

  response_vector_T& operator[](size_t b) { return V.at(b); }
  const response_vector_T& operator[](size_t b) const { return V.at(b); }

  // same size predicate (requirement)
  bool same_size(const response_space& other) const {
    return (num_vectors == other.num_vectors &&
            other.num_orbitals == other.num_orbitals);
  }
  friend bool same_size(const response_space V, const response_space& W) {
    return (V.num_vectors == W.num_vectors && W.num_orbitals == W.num_orbitals);
  }
  response_space operator+(const response_space W) {
    MADNESS_ASSERT(same_size(W));
    World world = W.x[0].world();
    response_space result(world, num_vectors, num_orbitals);
    for (size_t b = 0; b < num_vectors; b++) {
      result[b] = V[b] + W[b];
    }
    return result;
  }
  friend response_space operator+(const response_space V,
                                  const response_space W) {
    MADNESS_ASSERT(V.same_size(W));
    World world = V.x[0].world();
    response_space result(world, V.num_vectors, V.num_orbitals);
    for (size_t b = 0; b < V.num_vectors; b++) {
      result[b] = V[b] + W[b];
    }
    return result;
  }
  response_space operator-(const response_space W) {
    MADNESS_ASSERT(same_size(W));
    World world = W.x[0].world();
    response_space result(world, num_vectors, num_orbitals);
    for (size_t b = 0; b < num_vectors; b++) {
      result[b] = V[b] + W[b];
    }
    return result;
  }
  friend response_space operator-(const response_space V,
                                  const response_space W) {
    MADNESS_ASSERT(V.same_size(W));
    World world = V.x[0].world();
    response_space result(world, V.num_vectors, V.num_orbitals);
    for (size_t b = 0; b < V.num_vectors; b++) {
      result[b] = V[b] + W[b];
    }
    return result;
  }
  // KAIN must have this
  // Scaling by a constant
  response_space operator*(double a) {
    response_space result(*this);

    for (size_t b = 0; b < num_vectors; b++) {
      V[b] = V[b] * a;
    }

    V[0][0].world().gop.fence();
    return result;
  }

  // Scaling all internal functions by an external function
  // g[i][j] = x[i][j] * f
  response_space operator*(const function_T& f) {
    response_space result(*this);

    for (unsigned int b = 0; b < num_vectors; b++) {
      // Using vmra.h funciton
      result[b] = V[b] * f;
    }

    f.world().gop.fence();
    return result;
  }

  // KAIN must have this
  response_space operator+=(const response_space W) {
    MADNESS_ASSERT(same_size(W));
    for (unsigned int b = 0; b < num_vectors; b++) {
      V[b] += W.V[b];
    }
    return *this;
  }
  void zero() {
    for (size_t b = 0; b < num_vectors; b++) {
      V[b].zero();
    }
  }

  void compress() {
    for (size_t b = 0; b < num_vectors; b++) {
      V[b].compress();
    }
  }
  void reconstruct() {
    for (size_t b = 0; b < num_vectors; b++) {
      V[b].reconstruct();
    }
  }
  void truncate() {
    for (size_t b = 0; b < num_vectors; b++) {
      V[b].truncate();
    }
  }

  // Returns norms of each state
  Tensor<double> norm2() {
    Tensor<double> norms(num_vectors);
    for (size_t b = 0; b < num_orbitals; b++) {
      norms(b) = V[b].norm2();
    }
    return norms;
  }

  // Scales each state (read: entire row) by corresponding vector element
  //     new[i] = old[i] * mat[i]
  void scale(Tensor<double>& mat) {
    for (size_t b = 0; b < num_vectors; b++)
      V[b] = V[b] * mat[b];
    V[0].x.world().gop.fence();
  }

  Tensor<double> inner(response_space& V, response_space& W) {
    // assert the two spaces same size
    MADNESS_ASSERT(same_size(V, W));
    Tensor<double> values(num_vectors);
    for (size_t b = 0; b < num_vectors; b++) {
      values(b) = inner(V[b], W[b]);
    }
  }
};

// response allocator
template <typename T> struct response_space_allocator {
  typedef Function<T, 3> function_T;
  typedef response_vector<T> response_vector_T;
  World& world;
  const size_t num_vectors;
  const size_t num_orbitals;
  response_space_allocator(World& world, const size_t num_vectors,
                           const size_t num_orbitals)
      : world(world), num_vectors(num_vectors), num_orbitals(num_orbitals) {}
  // overloading ()
  response_vector_T operator()() {
    return response_vector_T(world, num_vectors, num_orbitals);
  }
  // Copy constructor
  response_space_allocator operator=(const response_space_allocator& other) {
    return response_space_allocator(world, other.num_orbitals);
  }
};

} // namespace madness

#endif // SRC_APPS_ADRIAN_RESPONSE_SPACE_H_
