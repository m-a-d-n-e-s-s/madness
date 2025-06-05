/*
 * Copyright 2021 Adrian Hurtado
 * Some basic operators for ResponseFunction objects
 */
#include "basic_operators.h"

namespace madness {

// Returns a shallow copy of the transpose of a vector of vector of functions
response_space transpose(response_space &f) {
  MADNESS_ASSERT(f.size() > 0);
  MADNESS_ASSERT(f[0].size() > 0);

  // Get sizes
  size_t m = f.size();
  size_t n = f[0].size();

  // Return container
  response_space g(f[0][0].world(), m, n);

  // Now do shallow copies
  for (size_t i = 0; i < m; i++)
    for (size_t j = 0; j < n; j++)
      g[j][i] = f[i][j];

  // Done
  return g;
}

// Truncate a vector of vector of functions
void truncate(World &world, response_space &v, double tol, bool fence) {
  MADNESS_ASSERT(v.size() > 0);
  MADNESS_ASSERT(v[0].size() > 0);

  for (unsigned int i = 0; i < v.size(); i++) {
    truncate(world, v[i], tol, fence);
  }
}

// Apply a vector of vector of operators to a vector of vector of functions
// g[i][j] = op[i][j](f[i][j])
response_space
apply(World &world,
      std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> &op,
      response_space &f) {
  MADNESS_ASSERT(f.size() > 0);
  MADNESS_ASSERT(f.size() == op.size());
  MADNESS_ASSERT(f[0].size() == op[0].size());

  response_space result(f[0][0].world(), f.size(), f[0].size());

  for (const auto &i : f.active) {
    // Using vmra.h function, line 889
    result[i] = apply(world, op[i], f[i]);
  }

  return result;
}
// Apply a vector of operators to a set of response states
//
response_space apply(World &world,
                     std::vector<std::shared_ptr<real_convolution_3d>> &op,
                     response_space &f) {
  MADNESS_ASSERT(f.size() > 0);
  MADNESS_ASSERT(f[0].size() == op.size());
  // response_space result = f.copy();
  response_space result =
      response_space(f[f.active.front()][0].world(), f.size(), f[0].size());
  result.active = f.active;
  for (const auto &i : f.active) {
    result.x[i] = apply(world, op, f[i]);
  }
  world.gop.fence();

  return result;
}

// Apply the derivative operator to a vector of vector of functions
response_space apply(World &world, real_derivative_3d &op, response_space &f) {
  MADNESS_ASSERT(f.size() > 0);

  response_space result =
      response_space(f[f.active.front()][0].world(), f.size(), f[0].size());
  result.active = f.active;
  for (const auto &i : f.active) {
    result.x[i] = apply(world, op, f[i]);
  }
  world.gop.fence();

  return result;
}
} // namespace madness

// Deuces
