/*
 * Some basic operators for ResponseFunction objects
 */
#ifndef SRC_APPS_ADRIAN_TDHF_BASIC_OPERATORS2_H_
#define SRC_APPS_ADRIAN_TDHF_BASIC_OPERATORS2_H_

#include <ResponseFunction2.h>
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

#include <algorithm>
#include <memory>
#include <vector>

namespace madness {
// Returns a shallow copy of the transpose of a vector of vector of functions
response_space transpose(response_space& f);

// Multiplication of a vector of vectors by a matrix,
//  *  g[i][k] = \sum_{j} a[i][j] * b(j,k)
// !  NOTE: NO BOUNDS CHECKING ON THE TENSOR b!!!!
response_space scale_2d(World& world, const response_space& a,
                         const Tensor<double>& b);

// Multiplication of a vector of vectors by a scalar g[i][j] = a[i][j] * b
response_space scale(response_space a, double b);

// Truncate a vector of vector of functions
void truncate(World& world, response_space& v,
              double tol = FunctionDefaults<3>::get_thresh(),
              bool fence = true);

// Apply a vector of vector of operators to a vector of vector of functions
// g[i][j] = op[i][j](f[i][j])
response_space
apply(World& world,
      std::vector<std::vector<std::shared_ptr<real_convolution_3d>>>& op,
      response_space& f);
// frequecy case
response_space apply(World& world,
                      std::vector<std::shared_ptr<real_convolution_3d>>& op,
                      response_space& f);

// Apply the derivative operator to a vector of vector of functions
response_space apply(World& world, real_derivative_3d& op, response_space& f);
} // namespace madness

#endif // SRC_APPS_ADRIAN_TDHF_BASIC_OPERATORS2_H_
