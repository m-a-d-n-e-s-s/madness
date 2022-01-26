//
// Created by adrianhurtado on 1/24/22.
//

#ifndef MADNESS_SOLVER_HPP
#define MADNESS_SOLVER_HPP

#include "global_functions.h"
#include "load_balance.h"
#include "timer.h"

class Solver {
 public:
  Solver(World& world, const CalcParams& params);

 protected:
  ResponseParameters r_params;
  GroundStateCalculation ground_calc;
  Molecule molecule;
  vector_real_function_3d ground_orbitals;
  Tensor<double> ground_energies;
  real_function_3d mask;
  X_space Chi;
  XCfunctional xcf;

  // Tensors for holding energies
  // residuals, and shifts
  Tensor<double> omega;        // Energies of response functions
  Tensor<double> e_residuals;  // Residuals of energies

  // Information that is inferred from input file
  // Ground state orbitals being used in calculation
  Tensor<double> hamiltonian;  // Ground state hamiltonian tensor
  Tensor<double> ham_no_diag;  // Ground state ham. without diagonal (Used when
  // localized orbitals are given)
  std::vector<int> active;        // The labels of orbitals selected as "active"
  unsigned int act_num_orbitals{};  // Number of ground state orbitals being used
  // in calculation

  std::shared_ptr<PotentialManager> potential_manager;

  // Mask function to handle boundary conditions

  // Functions
  real_function_3d stored_v_nuc;   // Stored nuclear potential from ground state
  real_function_3d stored_v_coul;  // Stored coulomb potential from ground state

  functionT rho0;
  vecfuncT rho_omega;

  response_space stored_potential;  // The ground state potential, stored only
  // if store_potential is true (default is
  poperatorT coulop;
  std::vector<std::shared_ptr<real_derivative_3d>> gradop;
  double vtol{};
  json j_molresponse;

  // Protected Member functions
  void set_protocol(World& world, double thresh);
  void check_k_Xspace(World& world, X_space& X, double thresh, size_t k);
  void check_k(World& world, double thresh, size_t k);
  real_function_3d Coulomb(World& world);
  functionT make_ground_density(World& world, const vecfuncT& v);
  Tensor<double> CreateGroundHamiltonian(World& world, size_t print_level);
  XCOperator<double, 3> create_XCOperator(World& world, std::string xc);
  void save(World& world, const std::string& name);
  void load(World& world, const std::string& name);
  vecfuncT make_density(World& world);
  vector<real_function_3d> transition_density(World& world,
                                              vector<real_function_3d>& orbitals,
                                              response_space& x,
                                              response_space& y);
  vector<real_function_3d> transition_densityTDA(World& world,
                                                 const vector<real_function_3d>& orbitals,
                                                 response_space& x);
  void load_balance(World& world);
  vector<poperatorT> make_bsh_operators_response(World& world, double& shift, double& omega) const;
  X_space Compute_Theta_X(World& world, X_space& Chi, XCOperator<double, 3> xc, std::string calc_type);
  X_space compute_F0X(World& world, X_space& Chi, XCOperator<double, 3> xc, bool compute_Y);
  X_space compute_V0X(World& world, X_space& X, XCOperator<double, 3> xc, bool compute_Y);
  void orbital_load_balance(World& world, vecfuncT& psi0, vecfuncT& psi0_copy, X_space& X, X_space& Chi_copy);
  X_space compute_gamma_tda(World& world, X_space& X, XCOperator<double, 3> xc);
  X_space compute_gamma_full(World& world, X_space& X, const XCOperator<double, 3>& xc);
  X_space compute_gamma_static(World& world, X_space& X, XCOperator<double, 3> xc);
  X_space compute_residual(World& world,
                           X_space& old_Chi,
                           X_space& temp,
                           Tensor<double>& bsh_residualsX,
                           Tensor<double>& bsh_residualsY,
                           std::string calc_type);
  X_space kain_x_space_update(World& world,
                              const X_space& temp,
                              const X_space& res,
                              NonLinearXsolver& kain_x_space,
                              vector<X_vector>& Xvector,
                              vector<X_vector>& Xresidual);
  void x_space_step_restriction(World& world,
                                X_space& old_Chi,
                                X_space& temp,
                                bool restrict_y,
                                Tensor<double>& maxrotn);
  void vector_stats(const vector<double>& v, double& rms, double& maxabsval) const;
  double do_step_restriction(World& world, const vecfuncT& x, vecfuncT& x_new, std::string spin) const;
  double do_step_restriction(World& world,
                             const vecfuncT& x,
                             vecfuncT& x_new,
                             std::string spin,
                             double maxrotn) const;
  double do_step_restriction(World& world,
                             const vecfuncT& x,
                             const vecfuncT& y,
                             vecfuncT& x_new,
                             vecfuncT& y_new,
                             std::string spin) const;
  void PlotGroundandResponseOrbitals(World& world,
                                     size_t iteration,
                                     response_space& x_response,
                                     response_space& y_response,
                                     const ResponseParameters& r_params,
                                     const GroundStateCalculation& g_params);
  void vector_stats_new(const Tensor<double> v, double& rms, double& maxabsval) const;
};

#endif  // MADNESS_SOLVER_HPP
