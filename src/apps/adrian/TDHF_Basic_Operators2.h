/*
 * Some basic operators for ResponseFunction objects
 */
#include <ResponseFunction2.h>
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <memory>
#include <vector>

#ifndef MADNESS_APPS_TDHF_OPS_H_INCLUDED
#define MADNESS_APPS_TDHF_OPS_H_INCLUDED

using namespace madness;

// Returns a shallow copy of the transpose of a vector of vector of functions
ResponseFunction transpose(ResponseFunction &f) {
  MADNESS_ASSERT(f.size() > 0);
  MADNESS_ASSERT(f[0].size() > 0);

  // Get sizes
  int m = f.size();
  int n = f[0].size();

  // Return container
  ResponseFunction g(f[0][0].world(), m, n);

  // Now do shallow copies
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      g[j][i] = f[i][j];

  // Done
  return g;
}

// Multiplication of a vector of vectors by a matrix,
//  *  g[i][k] = \sum_{j} a[i][j] * b(j,k)
// !  NOTE: NO BOUNDS CHECKING ON THE TENSOR b!!!!
ResponseFunction scale_2d(World &world, ResponseFunction &a,
                          Tensor<double> &b) {
  MADNESS_ASSERT(a.size() > 0);
  MADNESS_ASSERT(a[0].size() > 0);

  ResponseFunction result;

  for (unsigned int i = 0; i < a.size(); i++) {
    // Using vmra.h definitions
    std::vector<real_function_3d> tmp = transform(world, a[i], b, false);
    result.push_back(tmp);
  }

  // Does this need to be here?
  world.gop.fence();

  return result;
}

// Multiplication of a vector of vectors by a scalar g[i][j] = a[i][j] * b
ResponseFunction scale(ResponseFunction a, double b) {
  MADNESS_ASSERT(a.size() > 0);
  MADNESS_ASSERT(a[0].size() > 0);

  ResponseFunction result;

  for (unsigned int i = 0; i < a.size(); i++) {
    // Using vmra.h definitions
    // std::vector<real_function_3d> temp = a[i] * b;
    // result.push_back(temp);
    result.push_back(a[i] * b);
  }

  return result;
}

// Truncate a vector of vector of functions
void truncate(World &world, ResponseFunction &v, double tol = 0.0,
              bool fence = true) {
  MADNESS_ASSERT(v.size() > 0);
  MADNESS_ASSERT(v[0].size() > 0);

  for (unsigned int i = 0; i < v.size(); i++) {
    truncate(world, v[i], tol, fence);
  }
}

// Apply a vector of vector of operators to a vector of vector of functions
// g[i][j] = op[i][j](f[i][j])
ResponseFunction
apply(World &world,
      std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> &op,
      ResponseFunction &f) {
  MADNESS_ASSERT(f.size() > 0);
  MADNESS_ASSERT(f.size() == op.size());
  MADNESS_ASSERT(f[0].size() == op[0].size());

  ResponseFunction result(f[0][0].world(), f.size(), f[0].size());

  for (unsigned int i = 0; i < f.size(); i++) {
    // Using vmra.h function, line 889
    result[i] = apply(world, op[i], f[i]);
  }

  return result;
}

// Apply the derivative operator to a vector of vector of functions
ResponseFunction apply(World &world, real_derivative_3d &op,
                       ResponseFunction &f) {
  MADNESS_ASSERT(f.size() > 0);

  ResponseFunction result;

  for (unsigned int i = 0; i < f.size(); i++) {
    std::vector<real_function_3d> temp = apply(world, op, f[i]);
    result.push_back(temp);
  }

  return result;
}

#endif

// Deuces
