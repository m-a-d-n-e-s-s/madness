/*
 * Some basic operators for ResponseFunction objects
 */
#include "TDHF_Basic_Operators2.h"

namespace madness {

// Returns a shallow copy of the transpose of a vector of vector of functions
response_space transpose(response_space &f) {
  MADNESS_ASSERT(f.size() > 0);
  MADNESS_ASSERT(f[0].size() > 0);

  // Get sizes
  int m = f.size();
  int n = f[0].size();

  // Return container
  response_space g(f[0][0].world(), m, n);

  // Now do shallow copies
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++) g[j][i] = f[i][j];

  // Done
  return g;
}

// Multiplication of a vector of vectors by a matrix,
//  *  g[i][k] = \sum_{j} a[i][j] * b(j,k)
// !  NOTE: NO BOUNDS CHECKING ON THE TENSOR b!!!!
response_space scale_2d(World &world, const response_space &a,
                          const Tensor<double> &b) {
  MADNESS_ASSERT(a.size() > 0);
  MADNESS_ASSERT(a[0].size() > 0);

  response_space result;

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
response_space scale(response_space a, double b) {
  MADNESS_ASSERT(a.size() > 0);
  MADNESS_ASSERT(a[0].size() > 0);

//  print("double b",b);
  response_space result;

  for (unsigned int i = 0; i < a.size(); i++) {
    // Using vmra.h definitions
    // std::vector<real_function_3d> temp = a[i] * b;
    // result.push_back(temp);
    result.push_back(a[i] * b);
    
  }
 // print("norms of scale result");
 // print(result.norm2());

  return result;
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
response_space apply(
    World &world,
    std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> &op,
    response_space &f) {
  MADNESS_ASSERT(f.size() > 0);
  MADNESS_ASSERT(f.size() == op.size());
  MADNESS_ASSERT(f[0].size() == op[0].size());

  response_space result(f[0][0].world(), f.size(), f[0].size());

  for (unsigned int i = 0; i < f.size(); i++) {
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

  response_space result(f[0][0].world(), f.size(), f[0].size());

  for (unsigned int i = 0; i < f.size(); i++) {
    // Using vmra.h function, line 889
    /// Applies a vector of operators to a vector of functions --- q[i] =
    /// apply(op[i],f[i])
    result[i] = apply(world, op, f[i]);
  }

  return result;
}

// Apply the derivative operator to a vector of vector of functions
response_space apply(World &world, real_derivative_3d &op,
                       response_space &f) {
  MADNESS_ASSERT(f.size() > 0);

  response_space result;

  for (unsigned int i = 0; i < f.size(); i++) {
    std::vector<real_function_3d> temp = apply(world, op, f[i]);
    result.push_back(temp);
  }

  return result;
}
}  // namespace madness

// Deuces
