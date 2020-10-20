#ifndef SRC_APPS_ADRIAN_DENSITY_FREQUENCY_RESPONSE_FUNCTIONS_H_
#define SRC_APPS_ADRIAN_DENSITY_FREQUENCY_RESPONSE_FUNCTIONS_H_

#include <ResponseFunction2.h>
#include <TDDFT.h>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "../../madness/mra/funcplot.h"

namespace madness {
// base class for a density
// operator used to create it
// homogeneous sol----x and y functions
// particular sol --- depends on lower order functions used to create it
// it also needs an xc functional
// The Rparams and Gparmas used to create the density
//
class FirstOrderDensity {
 private:
  Tensor<double> omega;
  std::string property;

  int num_response_states;
  int num_ground_states;

  XCfunctional xcf;

  ResponseParameters Rparams;
  GroundParameters Gparams;

  ResponseFunction x;
  ResponseFunction y;

  // first order frequency response densities
  std::vector<real_function_3d> rho_omega;

 public:
  // Collective constructor
  FirstOrderDensity(World &world, const char *filename)
      : FirstOrderDensity(
            world,
            (world.rank() == 0 ? std::make_shared<std::ifstream>(filename)
                               : nullptr)) {}

  FirstOrderDensity(World &world, std::shared_ptr<std::istream> density_input) {
    TDHF calc(world, density_input);
    if (calc.Rparams.property) {
      calc.ComputeFrequencyResponse(world);
    } else {
      calc.solve(world);
    }
    // right now everything uses copy
    property = calc.Rparams.response_type;
    omega = calc.Rparams.omega;

    x = calc.GetResponseFunctions("x");
    y = calc.GetResponseFunctions("y");

    num_response_states = x.size();
    num_ground_states = x[0].size();
    // stuff i want to save
    Rparams = calc.GetResponseParameters();
    Gparams = calc.GetGroundParameters();
    // get the response densities for our states
    rho_omega = calc.transition_density(world, Gparams.orbitals, x, y);
  }
  void PrintDensityInformation() {
    // print
    //
    print("Response Density Information");
    print(property, " response at", omega(0, 0), "frequency using ", Rparams.xc,
          " exchange functional");
    print("Number of Response States : ", num_response_states);
    print("Number of Ground States : ", num_ground_states);
  }

  void PlotResponseDensity(World &world) {
    // Doing line plots along each axis
    // Doing line plots along each axis
    if (world.rank() == 0) print("\n\nStarting plots");
    coord_3d lo, hi;
    char plotname[500];
    double Lp = std::min(Gparams.L, 24.0);
    if (world.rank() == 0) print("x:");
    // x axis
    lo[0] = 0.0;
    lo[1] = 0.0;
    lo[2] = 0.0;
    hi[0] = Lp;
    hi[1] = 0.0;
    hi[2] = 0.0;

    for (int i = 0; i < num_response_states; i++) {
      std::snprintf(plotname, sizeof(plotname),
                    "plot_transition_density_%d_%d_x.plt",
                    FunctionDefaults<3>::get_k(), i);
      plot_line(plotname, 5001, lo, hi, rho_omega[i]);
    }
  }
};

typedef Tensor<double> tensorT;
typedef Function<double, 3> functionT;
typedef std::shared_ptr<FunctionFunctorInterface<double, 3>> functorT;
typedef FunctionFactory<double, 3> factoryT;

typedef Vector<double, 3> coordT;

class MolecularDerivativeFunctor : public FunctionFunctorInterface<double, 3> {
  typedef Vector<double, 3> coordT;

 private:
  const Molecule &molecule;
  const int atom;
  const int axis;

 public:
  MolecularDerivativeFunctor(const Molecule &molecule, int atom, int axis)
      : molecule(molecule), atom(atom), axis(axis) {}

  double operator()(const coordT &x) const {
    return molecule.nuclear_attraction_potential_derivative(atom, axis, x[0],
                                                            x[1], x[2]);
  }

  std::vector<coordT> special_points() const {
    return std::vector<coordT>(1, molecule.get_atom(atom).get_coords());
  }
};  // namespace madness
class Property {
 public:
  std::string property;
  std::vector<real_function_3d> operator_vector;

  Property(World &world, std::string property_type) : operator_vector(3) {
    property = property_type;
    MADNESS_ASSERT(property.compare("dipole") == 0);
    for (int i = 0; i < 3; i++) {
      std::vector<int> f(3, 0);
      f[i] = 1;
      operator_vector.push_back(real_factory_3d(world).functor(
          real_functor_3d(new BS_MomentFunctor(f))));
    }
  }

  Property(World &world, std::string property_type, Molecule &molecule)
      : operator_vector(molecule.natom() * 3) {
    property = property_type;
    MADNESS_ASSERT(property.compare("nuclear") == 0);

    vecfuncT dv(molecule.natom() * 3);  // default constructor for vector?

    for (size_t atom = 0; atom < molecule.natom(); ++atom) {
      for (int axis = 0; axis < 3; ++axis) {
        // question here....MolecularDerivativeFunctor takes derivative with
        // respect to axis atom and axis
        functorT func(new MolecularDerivativeFunctor(molecule, atom, axis));
        // here we save
        operator_vector[atom * 3 + axis] = functionT(factoryT(world)
                                                         .functor(func)
                                                         .nofence()
                                                         .truncate_on_project()
                                                         .truncate_mode(0));
        // need to project
        //        operator_vector[atom * 3 + axis] = mul_sparse(
        //           world, dv[atom * 3 + axis], Gparams.orbitals,
        //           Rparams.small);

        // project rhs vectors for state

        // core projector contribution
      }
    }
  }
};

Tensor<double> ComputeSecondOrderPropertyTensor(World &world,
                                                FirstOrderDensity rho_b,
                                                Property rho_c);

// namespace madness
#endif  // SRC_APPS_ADRIAN_DENSITY_FREQUENCY_RESPONSE_FUNCTIONS_H_
