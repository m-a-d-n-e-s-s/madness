// Copyright 2021 Adrian Hurtado
#include "molresponse/density.h"

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "../../madness/mra/funcplot.h"
#include "TDDFT.h"
#include "molresponse/global_functions.h"
#include "molresponse/property.h"
#include "molresponse/response_functions.h"

typedef Tensor<double> TensorT;
typedef Function<double, 3> FunctionT;
typedef std::shared_ptr<FunctionFunctorInterface<double, 3>> FunctorT;
typedef FunctionFactory<double, 3> FactoryT;
typedef Vector<double, 3> CoordinateT;
typedef std::vector<real_function_3d> VectorFunction3DT;

// base class for a density
// operator used to create it
// homogeneous sol----x and y functions
// particular sol --- depends on lower order functions used to create it
// it also needs an xc functional
// The Rparams and Gparmas used to create the density
//
FirstOrderDensity::FirstOrderDensity(ResponseParameters Rparams,
                                     GroundParameters Gparams) {
  this->Rparams = Rparams;
  this->Gparams = Gparams;
}
void FirstOrderDensity::ComputeResponse(World &world) {
  // right now everything uses copy
  property = Rparams.response_type;
  // TDHF sets the TDHF with paramaters...
  // sets mask function
  // brodcast serializable
  // sets Function Defaults for calculation
  // sets xcf
  //
  // creating calc should also set up the x and y functions
  //
  x = response_space(world, Rparams.states, Gparams.num_orbitals);
  y = response_space(world, Rparams.states, Gparams.num_orbitals);
  print("Creating Response Functions for X and Y");
  print("X Norms before Computing");
  print(x.norm2());
  print("Y Norms before Computing");
  print(y.norm2());

  TDDFT calc(world, Rparams, Gparams);
  if (calc.Rparams.response_type.compare("excited_state") == 0) {
    print("Entering Excited State Response Runner");
    calc.solve_excited_states(world);
  } else {
    print("Entering Frequency Response Runner");
    calc.ComputeFrequencyResponse(world, property, x, y);
  }
  // omega is determined by the type of calculation
  // property calculation at single frequency
  // excited stat calculation at multipe frequencies
  omega = calc.GetFrequencyOmega();
  property_operator = calc.GetPropertyObject();

  x = calc.GetResponseFunctions("x");
  y = calc.GetResponseFunctions("y");

  print("X Norms before Computing");
  print(x.norm2());
  P = calc.GetPVector();
  Q = calc.GetQVector();

  num_states = x.size();
  num_ground_states = x[0].size();
  // get the response densities for our states
  if (Rparams.omega == 0) {
    rho_omega = calc.transition_density(world, Gparams.orbitals, x, x);
  } else {
    rho_omega = calc.transition_density(world, Gparams.orbitals, x, y);
  }
  if (Rparams.save_density) {
    SaveDensity(world, Rparams.save_density_file);
  }
}

// right now everything uses copy

size_t FirstOrderDensity::GetNumberResponseStates() { return num_states; }
size_t FirstOrderDensity::GetNumberGroundStates() { return num_ground_states; }
VectorFunction3DT FirstOrderDensity::GetDensityVector() { return rho_omega; }
const Molecule FirstOrderDensity::GetMolecule() { return Gparams.molecule; }
TensorT FirstOrderDensity::GetFrequencyOmega() { return omega; }
ResponseParameters FirstOrderDensity::GetResponseParameters() {
  return Rparams;
}

VectorFunction3DT FirstOrderDensity::ComputeDensityVector(World &world,
                                                          bool is_static) {
  std::vector<real_function_3d> densities =
      zero_functions<double, 3>(world, num_states);
  /*
    x.reconstruct_rf();
    y.reconstruct_rf();
    reconstruct(world, Gparams.orbitals);
    */

  if (is_static) {
    for (size_t b = 0; b < num_states; b++) {
      densities[b] = dot(world, x[b], Gparams.orbitals) +
                     dot(world, x[b], Gparams.orbitals);
      /*
        for (size_t    j = 0; j < num_ground_states; j++) {
          densities[b] += mul_sparse(x[b][j], Gparams.orbitals[j],
        Rparams.small); densities[b] += mul_sparse(x[b][j], Gparams.orbitals[j],
        Rparams.small);
        }
        */
    }
  } else {
    for (size_t b = 0; b < num_states; b++) {
      densities[b] = dot(world, x[b], Gparams.orbitals) +
                     dot(world, y[b], Gparams.orbitals);
      /*
        for (size_t    j = 0; j < num_ground_states; j++) {
          densities[b] += mul_sparse(x[b][j], Gparams.orbitals[j],
        Rparams.small); densities[b] += mul_sparse(y[b][j], Gparams.orbitals[j],
        Rparams.small);
        }
        */
    }
  }
  truncate(world, densities);
  return densities;
}
void FirstOrderDensity::PrintDensityInformation() {
  // print
  //
  print("Response Density Information");
  print(property,
        " response at",
        omega(0, 0),
        "frequency using ",
        Rparams.xc,
        " exchange functional");
  print("Number of Response States : ", num_states);
  print("Number of Ground States : ", num_ground_states);
}

void FirstOrderDensity::PlotResponseDensity(World &world) {
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

  for (size_t i = 0; i < num_states; i++) {
    std::snprintf(plotname,
                  sizeof(plotname),
                  "plot_transition_density_%d_%d_x.plt",
                  FunctionDefaults<3>::get_k(),
                  static_cast<int>(i));
    plot_line(plotname, 5001, lo, hi, rho_omega[i]);
  }
}
Tensor<double> FirstOrderDensity::ComputeSecondOrderPropertyTensor(
    World &world) {
  Tensor<double> G(num_states, num_states);
  response_space grp(world, num_states, num_states);

  for (size_t i(0); i < num_states; i++) {
    for (size_t j(0); j < num_states; j++) {
      grp[i][j] = dot(world, P[i], x[j]) + dot(world, Q[i], y[j]);
      G(i, j) = grp[i][j].trace();
      G(i, j) = -2 * G(i, j);
    }
  }

  // Tensor<double> M(num_response_states, num_response_states);
  // do some printing before we compute so we know what we are working with
  //*******************************
  // G algorithim

  print("G");
  print(G);

  // print(M);

  return G;
}

void FirstOrderDensity::PrintSecondOrderAnalysis(
    World &world,
    const Tensor<double> alpha_tensor) {
  Tensor<double> V, epolar;
  syev(alpha_tensor, V, epolar);
  double Dpolar_average = 0.0;
  double Dpolar_iso = 0.0;
  for (size_t i = 0; i < 3; ++i) Dpolar_average = Dpolar_average + epolar[i];
  Dpolar_average = Dpolar_average / 3.0;
  Dpolar_iso =
      sqrt(.5) * sqrt(std::pow(alpha_tensor(0, 0) - alpha_tensor(1, 1), 2) +
                      std::pow(alpha_tensor(1, 1) - alpha_tensor(2, 2), 2) +
                      std::pow(alpha_tensor(2, 2) - alpha_tensor(0, 0), 2));

  size_t num_states = Rparams.states;

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

    for (size_t i = 0; i < num_states; i++) {
      print(epolar[i]);
    }
  }
}
void FirstOrderDensity::SaveDensity(World &world, std::string name) {
  // Archive to write everything to
  archive::ParallelOutputArchive ar(world, name.c_str(), 1);
  // Just going to enforce 1 io server

  ar &property;
  ar &omega;
  ar &num_states;
  ar &num_ground_states;
  // Save response functions x and y
  // x first
  for (size_t i = 0; i < num_states; i++) {
    for (size_t j = 0; j < num_ground_states; j++) {
      ar &x[i][j];
    }
  }

  // y second
  for (size_t i = 0; i < num_states; i++) {
    for (size_t j = 0; j < num_ground_states; j++) {
      ar &y[i][j];
    }
  }
  for (size_t i = 0; i < num_states; i++) {
    ar &rho_omega[i];
  }
  for (size_t i = 0; i < property_operator.num_operators; i++) {
    ar &property_operator.operator_vector[i];
  }

  for (size_t i = 0; i < num_states; i++) {
    for (size_t j = 0; j < num_ground_states; j++) {
      ar &P[i][j];
    }
  }
  for (size_t i = 0; i < num_states; i++) {
    for (size_t j = 0; j < num_ground_states; j++) {
      ar &Q[i][j];
    }
  }
}
// Load a response calculation
void FirstOrderDensity::LoadDensity(World &world,
                                    std::string name,
                                    ResponseParameters Rparams,
                                    GroundParameters Gparams) {
  // create XCF Object
  xcf.initialize(Rparams.xc, false, world, true);

  archive::ParallelInputArchive ar(world, name.c_str());
  // Reading in, in this order;

  ar &property;

  if (property.compare("dipole") == 0) {
    if (world.rank() == 0) print("creating dipole property operator");
    this->property_operator = Property(world, "dipole");
  } else if (property.compare("nuclear") == 0) {
    if (world.rank() == 0) print("creating nuclear property operator");
    this->property_operator = Property(world, "nuclear", Gparams.molecule);
  }
  print("property:", property);

  ar &omega;
  print("omega:", omega);
  ar &num_states;
  print("num_response_states:", num_states);
  ar &num_ground_states;
  print("num_ground_states:", num_ground_states);

  this->x = response_space(world, num_states, num_ground_states);
  this->y = response_space(world, num_states, num_ground_states);

  this->P = response_space(world, num_states, num_ground_states);
  this->Q = response_space(world, num_states, num_ground_states);

  for (size_t i = 0; i < Rparams.states; i++) {
    for (size_t j = 0; j < Gparams.num_orbitals; j++) {
      ar &x[i][j];
      print("norm of x ", x[i][j].norm2());
    }
  }
  world.gop.fence();

  for (size_t i = 0; i < Rparams.states; i++) {
    for (size_t j = 0; j < Gparams.num_orbitals; j++) {
      ar &y[i][j];
      print("norm of y ", y[i][j].norm2());
    }
  }

  world.gop.fence();
  this->rho_omega = zero_functions<double, 3>(world, num_states);
  for (size_t i = 0; i < num_states; i++) {
    ar &rho_omega[i];
    print("norm of rho_omega ", rho_omega[i].norm2());
  }

  for (size_t i = 0; i < property_operator.num_operators; i++) {
    print("norm of operator before ",
          property_operator.operator_vector[i].norm2());
    ar &property_operator.operator_vector[i];
    print("norm of operator after",
          property_operator.operator_vector[i].norm2());
  }

  for (size_t i = 0; i < Rparams.states; i++) {
    for (size_t j = 0; j < Gparams.num_orbitals; j++) {
      ar &P[i][j];
      print("norm of P ", P[i][j].norm2());
    }
  }
  world.gop.fence();

  for (size_t i = 0; i < Rparams.states; i++) {
    for (size_t j = 0; j < Gparams.num_orbitals; j++) {
      ar &Q[i][j];
      print("norm of y ", Q[i][j].norm2());
    }
  }

  world.gop.fence();
}