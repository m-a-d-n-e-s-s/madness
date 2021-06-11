#ifndef SRC_APPS_molresponse_PROPERTY_FUNCTIONS_H_
#define SRC_APPS_molresponse_PROPERTY_FUNCTIONS_H_

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "molresponse/property_operators.h"

typedef Tensor<double> TensorT;
typedef Function<double, 3> FunctionT;
typedef std::shared_ptr<FunctionFunctorInterface<double, 3>> FunctorT;
typedef FunctionFactory<double, 3> FactoryT;
typedef Vector<double, 3> CoordinateT;
typedef std::vector<real_function_3d> VectorFunction3DT;

// returns tensor alpha_ij=inner(rho_b,op_c)
Tensor<double> ComputeSecondOrderPropertyTensor(World &world,
                                                const VectorFunction3DT &rho_b,
                                                const Property &op_c);

void PrintSecondOrderAnalysis(World &world, const Tensor<double> alpha_tensor,
                              const Tensor<double> omega,
                              const ResponseParameters r_params);

#endif  // SRC_APPS_molresponse_PROPERTY_FUNCTIONS_H_
