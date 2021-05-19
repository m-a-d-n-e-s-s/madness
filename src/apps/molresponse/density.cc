// Copyright 2021 Adrian Hurtado
#include "molresponse/density.h"

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "../../madness/mra/funcplot.h"
#include "TDDFT.h"
#include "molresponse/global_functions.h"
#include "molresponse/ground_parameters.h"
#include "molresponse/property.h"
#include "molresponse/response_functions.h"
#include "molresponse/response_parameters.h"

typedef Tensor<double> TensorT;
typedef Function<double, 3> FunctionT;
typedef std::shared_ptr<FunctionFunctorInterface<double, 3>> FunctorT;
typedef FunctionFactory<double, 3> FactoryT;
typedef Vector<double, 3> CoordinateT;
typedef std::vector<real_function_3d> VectorFunction3DT;

density_vector::density_vector(World &world, ResponseParameters r_params, GroundParameters g_params) {
  this->r_params = r_params;
  this->g_params = g_params;
  this->property = r_params.response_type();
  this->num_states = r_params.n_states();
  this->num_ground_states = r_params.num_orbitals();
  this->Chi = X_space(world, num_states, num_ground_states);
  this->PQ = X_space(world, num_states, num_ground_states);
}
void density_vector::compute_response(World &world) {
  // right now everything uses copy
  property = r_params.response_type();
  size_t m_states = r_params.n_states();
  size_t n_orbitals = r_params.num_orbitals();

  Chi = X_space(world, m_states, n_orbitals);
  PQ = X_space(world, m_states, n_orbitals);
  print("Creating Response Functions for X and Y");
  print("X Norms before Computing");
  print(Chi.X.norm2());
  print("Y Norms before Computing");
  print(Chi.Y.norm2());

  TDDFT calc(world, r_params, g_params);
  if (calc.r_params.response_type().compare("excited_state") == 0) {
    print("Entering Excited State Response Runner");
    calc.solve_excited_states(world);
  } else {
    print("Entering Frequency Response Runner");
    calc.compute_freq_response(world, property, Chi, PQ);
  }
  // omega is determined by the type of calculation
  // property calculation at single frequency
  // excited stat calculation at multipe frequencies
  omega = calc.GetFrequencyOmega();
  property_operator = calc.GetPropertyObject();

  Chi = calc.GetXspace();

  print("X Norms before Computing");
  print(Chi.X.norm2());
  PQ = calc.GetPQspace();

  num_states = Chi.X.size();
  num_ground_states = Chi.X[0].size();
  // get the response densities for our states
  if (r_params.omega() == 0) {
    rho_omega = calc.transition_density(world, g_params.orbitals(), Chi.X, Chi.X);
  } else {
    rho_omega = calc.transition_density(world, g_params.orbitals(), Chi.X, Chi.Y);
  }
  if (r_params.save_density()) {
    SaveDensity(world, r_params.save_density_file());
  }
}

// right now everything uses copy

size_t density_vector::GetNumberResponseStates() { return num_states; }
size_t density_vector::GetNumberGroundStates() { return num_ground_states; }
VectorFunction3DT density_vector::GetDensityVector() { return rho_omega; }
const Molecule density_vector::GetMolecule() { return g_params.molecule(); }
TensorT density_vector::GetFrequencyOmega() { return omega; }
ResponseParameters density_vector::GetResponseParameters() { return r_params; }

VectorFunction3DT density_vector::ComputeDensityVector(World &world, bool is_static) {
  std::vector<real_function_3d> densities = zero_functions<double, 3>(world, num_states);
  if (is_static) {
    for (size_t b = 0; b < num_states; b++) {
      densities[b] = dot(world, Chi.X[b], g_params.orbitals()) + dot(world, Chi.X[b], g_params.orbitals());
    }
  } else {
    for (size_t b = 0; b < num_states; b++) {
      densities[b] = dot(world, Chi.X[b], g_params.orbitals()) + dot(world, Chi.Y[b], g_params.orbitals());
    }
  }
  truncate(world, densities);
  return densities;
}
void density_vector::PrintDensityInformation() {
  // print
  //
  print("Response Density Information");
  print(property, " response at", omega(0, 0), "frequency using ", r_params.xc(), " exchange functional");
  print("Number of Response States : ", num_states);
  print("Number of Ground States : ", num_ground_states);
}

void density_vector::PlotResponseDensity(World &world) {
  // Doing line plots along each axis
  if (world.rank() == 0) print("\n\nStarting plots");
  coord_3d lo, hi;
  char plotname[500];
  double Lp = std::min(g_params.get_L(), 24.0);
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
Tensor<double> density_vector::ComputeSecondOrderPropertyTensor(World &world) {
  Tensor<double> H = -2 * inner(Chi, PQ);
  Tensor<double> G(num_states, num_states);
  response_space grp(world, num_states, num_states);

  for (size_t i(0); i < num_states; i++) {
    for (size_t j(0); j < num_states; j++) {
      grp[i][j] = dot(world, PQ.X[i], Chi.X[j]) + dot(world, PQ.Y[i], Chi.Y[j]);
      G(i, j) = grp[i][j].trace();
      G(i, j) = -2 * G(i, j);
    }
  }

  // Tensor<double> M(num_response_states, num_response_states);
  // do some printing before we compute so we know what we are working with
  //*******************************
  // G algorithim
  print("version 1");
  print(H);

  print("version 2");
  print(G);

  // print(M);

  return H;
}

void density_vector::PrintSecondOrderAnalysis(World &world, const Tensor<double> alpha_tensor) {
  Tensor<double> V, epolar;
  syev(alpha_tensor, V, epolar);
  double Dpolar_average = 0.0;
  double Dpolar_iso = 0.0;
  for (size_t i = 0; i < 3; ++i) Dpolar_average = Dpolar_average + epolar[i];
  Dpolar_average = Dpolar_average / 3.0;
  Dpolar_iso = sqrt(.5) * sqrt(std::pow(alpha_tensor(0, 0) - alpha_tensor(1, 1), 2) +
                               std::pow(alpha_tensor(1, 1) - alpha_tensor(2, 2), 2) +
                               std::pow(alpha_tensor(2, 2) - alpha_tensor(0, 0), 2));

  size_t num_states = r_params.n_states();

  if (world.rank() == 0) {
    print("\nTotal Dynamic Polarizability Tensor");
    printf("\nFrequency  = %.6f a.u.\n\n", omega(0, 0));
    // printf("\nWavelength = %.6f a.u.\n\n", r_params.omega() * ???);
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
void density_vector::SaveDensity(World &world, std::string name) {
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
      ar &Chi.X[i][j];
    }
  }

  // y second
  for (size_t i = 0; i < num_states; i++) {
    for (size_t j = 0; j < num_ground_states; j++) {
      ar &Chi.Y[i][j];
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
      ar &PQ.X[i][j];
    }
  }
  for (size_t i = 0; i < num_states; i++) {
    for (size_t j = 0; j < num_ground_states; j++) {
      ar &PQ.Y[i][j];
    }
  }
}
// Load a response calculation
void density_vector::LoadDensity(World &world,
                                 std::string name,
                                 ResponseParameters r_params,
                                 GroundParameters g_params) {
  // create XCF Object
  xcf.initialize(r_params.xc(), false, world, true);

  archive::ParallelInputArchive ar(world, name.c_str());
  // Reading in, in this order;

  ar &property;

  if (property.compare("dipole") == 0) {
    if (world.rank() == 0) print("creating dipole property operator");
    this->property_operator = Property(world, "dipole");
  } else if (property.compare("nuclear") == 0) {
    if (world.rank() == 0) print("creating nuclear property operator");
    this->property_operator = Property(world, "nuclear", g_params.molecule());
  }
  print("property:", property);

  ar &omega;
  print("omega:", omega);
  ar &num_states;
  print("num_response_states:", num_states);
  ar &num_ground_states;
  print("num_ground_states:", num_ground_states);

  this->Chi.X = response_space(world, num_states, num_ground_states);
  this->Chi.Y = response_space(world, num_states, num_ground_states);

  this->PQ.X = response_space(world, num_states, num_ground_states);
  this->PQ.Y = response_space(world, num_states, num_ground_states);

  for (size_t i = 0; i < r_params.n_states(); i++) {
    for (size_t j = 0; j < g_params.n_orbitals(); j++) {
      ar &Chi.X[i][j];
      print("norm of x ", Chi.X[i][j].norm2());
    }
  }
  world.gop.fence();

  for (size_t i = 0; i < r_params.n_states(); i++) {
    for (size_t j = 0; j < g_params.n_orbitals(); j++) {
      ar &Chi.Y[i][j];
      print("norm of y ", Chi.Y[i][j].norm2());
    }
  }

  world.gop.fence();
  this->rho_omega = zero_functions<double, 3>(world, num_states);
  for (size_t i = 0; i < num_states; i++) {
    ar &rho_omega[i];
    print("norm of rho_omega ", rho_omega[i].norm2());
  }

  for (size_t i = 0; i < property_operator.num_operators; i++) {
    print("norm of operator before ", property_operator.operator_vector[i].norm2());
    ar &property_operator.operator_vector[i];
    print("norm of operator after", property_operator.operator_vector[i].norm2());
  }

  for (size_t i = 0; i < r_params.n_states(); i++) {
    for (size_t j = 0; j < g_params.n_orbitals(); j++) {
      ar &PQ.X[i][j];
      print("norm of P ", PQ.X[i][j].norm2());
    }
  }
  world.gop.fence();

  for (size_t i = 0; i < r_params.n_states(); i++) {
    for (size_t j = 0; j < g_params.n_orbitals(); j++) {
      ar &PQ.Y[i][j];
      print("norm of y ", PQ.Y[i][j].norm2());
    }
  }

  world.gop.fence();
}

density_vector set_density_type(World &world, ResponseParameters R, GroundParameters G) {
  if (R.response_type().compare("excited_state") == 0) {
    return excited_state_density_vector(world, R, G);
  } else if (R.response_type().compare("dipole") == 0) {
    return dipole_density_vector(world, R, G);

  } else if (R.response_type().compare("nuclear") == 0) {
    return nuclear_density_vector(world, R, G);
  } else if (R.response_type().compare("order2") == 0) {
    MADNESS_EXCEPTION("not implemented yet", 0);
    return density_vector(world, R, G);
  } else if (R.response_type().compare("order3") == 0) {
    MADNESS_EXCEPTION("not implemented yet", 0);
    return density_vector(world, R, G);

  } else {
    MADNESS_EXCEPTION("what is this????", 0);
    return density_vector(world, R, G);
  }
};
