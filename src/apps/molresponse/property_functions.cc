

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
                                                const Property &op_c) {
  return matrix_inner(world, rho_b, op_c.operator_vector, true);
}

void PrintSecondOrderAnalysis(World &world, const Tensor<double> alpha_tensor,
                              const Tensor<double> omega,
                              const ResponseParameters Rparams) {
  Tensor<double> V, epolar;
  syev(alpha_tensor, V, epolar);
  double Dpolar_average = 0.0;
  double Dpolar_iso = 0.0;
  for (unsigned int i = 0; i < 3; ++i)
    Dpolar_average = Dpolar_average + epolar[i];
  Dpolar_average = Dpolar_average / 3.0;
  Dpolar_iso =
      sqrt(.5) * sqrt(std::pow(alpha_tensor(0, 0) - alpha_tensor(1, 1), 2) +
                      std::pow(alpha_tensor(1, 1) - alpha_tensor(2, 2), 2) +
                      std::pow(alpha_tensor(2, 2) - alpha_tensor(0, 0), 2));

size_t   num_states = Rparams.states;

  if (world.rank() == 0) {
    print("\nTotal Dynamic Polarizability Tensor");
    printf("\nFrequency  = %.6f a.u.\n\n", omega(0, 0));
    // printf("\nWavelength = %.6f a.u.\n\n", Rparams.omega * ???);
    print(alpha_tensor);
    printf("\tEigenvalues = ");
    printf("\t %.6f \t %.6f \t %.6f \n", epolar[0], epolar[1], epolar[2]);
    printf("\tIsotropic   = \t %.6f \n", Dpolar_average);
    printf("\tAnisotropic = \t %.6f \n", Dpolar_iso);
    printf("\n");

    for (long i = 0; i < num_states; i++) {
      print(epolar[i]);
    }
  }
}
