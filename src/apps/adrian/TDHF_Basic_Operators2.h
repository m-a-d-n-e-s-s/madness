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
ResponseFunction transpose(ResponseFunction& f);

// Multiplication of a vector of vectors by a matrix,
//  *  g[i][k] = \sum_{j} a[i][j] * b(j,k)
// !  NOTE: NO BOUNDS CHECKING ON THE TENSOR b!!!!
ResponseFunction scale_2d(World& world, const ResponseFunction& a,
                          const Tensor<double>& b);

// Multiplication of a vector of vectors by a scalar g[i][j] = a[i][j] * b
ResponseFunction scale(ResponseFunction a, double b);

// Truncate a vector of vector of functions
void truncate(World& world, ResponseFunction& v,
              double tol = FunctionDefaults<3>::get_thresh(),
              bool fence = true);

// Apply a vector of vector of operators to a vector of vector of functions
// g[i][j] = op[i][j](f[i][j])
ResponseFunction apply(
    World& world,
    std::vector<std::vector<std::shared_ptr<real_convolution_3d>>>& op,
    ResponseFunction& f);
// frequecy case
ResponseFunction apply(World& world,
                       std::vector<std::shared_ptr<real_convolution_3d>>& op,
                       ResponseFunction& f);

// Apply the derivative operator to a vector of vector of functions
ResponseFunction apply(World& world, real_derivative_3d& op,
                       ResponseFunction& f);
}  // namespace madness

#endif  // SRC_APPS_ADRIAN_TDHF_BASIC_OPERATORS2_H_
