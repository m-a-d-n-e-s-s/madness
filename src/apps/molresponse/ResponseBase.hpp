//
// Created by adrianhurtado on 1/24/22.
//

#ifndef MADNESS_RESPONSEBASE_HPP
#define MADNESS_RESPONSEBASE_HPP

#include <functional>
#include <numeric>
#include <utility>
#include <vector>

#include "global_functions.h"
#include "madness/chem/SCF.h"
#include "madness/mra/functypedefs.h"
#include "madness/mra/mra.h"
#include "madness/mra/nonlinsol.h"
#include "madness/tensor/tensor.h"
#include "madness/tensor/tensor_json.hpp"
#include "response_macrotask.hpp"
#include "timer.h"
#include "x_space.h"

using namespace madness;

// I need to creat a load X_space function strategy which will be used by the
// FrequencyResponse class and the ExcitedResponse class. This will allow me to
// load the X_space in the FrequencyResponse class and the Quadratic response
// class.  This way I can use the same X_space for both classes. and load
// depending on the type of response I am doing.

typedef std::pair<X_space, Tensor<double>> XData;
typedef std::pair<vector_real_function_3d, vector_real_function_3d>
    double_response_vector;

class LoadXSpaceStrategy {
public:
  virtual ~LoadXSpaceStrategy() = default;
  virtual XData load_x_space(World &world, const std::string &filename,
                             const ResponseParameters &r_params,
                             double omega_state) const = 0;
};
class LoadFrequencyXSpace : public LoadXSpaceStrategy {
public:
  XData load_x_space(World &world, const std::string &filename,
                     const ResponseParameters &r_params,
                     double omega_state) const override {
    if (world.rank() == 0) {
      print("FrequencyResponse::load() -state");
      print("Loading X_space from file: ", filename);
      print("Number of states: ", r_params.num_states());
      print("Number of orbitals: ", r_params.num_orbitals());
      print("Loading omega : ", omega_state);
    }
    //

    X_space chi_new(world, r_params.num_states(), r_params.num_orbitals());

    // The archive to read from
    archive::ParallelInputArchive ar(world, filename.c_str());
    ar & r_params.archive();
    ar & r_params.tda();
    ar & r_params.num_orbitals();
    ar & r_params.num_states();
    for (size_t i = 0; i < r_params.num_states(); i++)
      for (size_t j = 0; j < r_params.num_orbitals(); j++)
        ar & chi_new.x[i][j];
    world.gop.fence();
    if (omega_state == 0.0) {
      for (size_t i = 0; i < r_params.num_states(); i++)
        for (size_t j = 0; j < r_params.num_orbitals(); j++)
          chi_new.y[i][j] = chi_new.x[i][j];
    } else {
      for (size_t i = 0; i < r_params.num_states(); i++)
        for (size_t j = 0; j < r_params.num_orbitals(); j++)
          ar & chi_new.y[i][j];
    }
    world.gop.fence();

    return {chi_new, Tensor<double>()};
  }
};
class LoadExcitedXSpace : public LoadXSpaceStrategy {
public:
  XData load_x_space(World &world, const std::string &filename,
                     const ResponseParameters &r_params,
                     double omega_state) const override {

    Tensor<double> omega;
    auto new_chi =
        X_space(world, r_params.num_states(), r_params.num_orbitals());
    // The archive to read from
    archive::ParallelInputArchive ar(world, filename.c_str());

    // Reading in, in this order;
    //  string           ground-state archive name (garch_name)
    //  bool             TDA flag
    // size_t                number of ground state orbitals (n)
    // size_t                number of excited state orbitals (m)
    //  Tensor<double>   energies of m x-components
    //  for i from 0 to m-1
    //     for j from 0 to n-1
    //        Function<double,3> x_response[i][j]
    //  (If TDA flag == True)
    //  (Tensor<double>  energies of m y-components    )
    //  (for i from 0 to m-1                       )
    //  (   for j from 0 to n-1                    )
    //  (      Function<double,3> y_response[i][j] )

    ar & r_params.archive();
    ar & r_params.tda();
    ar & r_params.num_orbitals();
    ar & r_params.num_states();
    ar & omega;

    for (size_t i = 0; i < r_params.num_states(); i++)
      for (size_t j = 0; j < r_params.num_orbitals(); j++)
        ar & new_chi.x[i][j];
    world.gop.fence();

    if (not r_params.tda()) {
      for (size_t i = 0; i < r_params.num_states(); i++)
        for (size_t j = 0; j < r_params.num_orbitals(); j++)
          ar & new_chi.y[i][j];
      world.gop.fence();
    }
    return {new_chi, omega};
  }
};

class ComputeDensityStrategy {
public:
  virtual ~ComputeDensityStrategy() = default;
  virtual vector_real_function_3d
  compute_density(World &world, const X_space &x,
                  const vector_real_function_3d &phi0,
                  const vector_real_function_3d &rho1, bool update) const = 0;
};

class StaticDensityStrategy : public ComputeDensityStrategy {
public:
  vector_real_function_3d compute_density(World &world, const X_space &x,
                                          const vector_real_function_3d &phi0,
                                          const vector_real_function_3d &rho1,
                                          bool update) const override {

    vector_real_function_3d rho_new;
    if (update) {
      rho_new = copy(world, rho1);
    } else {
      rho_new = zero_functions<double, 3>(world, x.num_states());
    }
    vector_real_function_3d x_phi, y_phi;

    for (const auto &b : x.active) {

      rho_new[b] = 2 * dot(world, x.x[b], phi0);
    }
    world.gop.fence();
    truncate(world, rho_new);
    return rho_new;
  }
};
class FullDensityStrategy : public ComputeDensityStrategy {
public:
  vector_real_function_3d compute_density(World &world, const X_space &x,
                                          const vector_real_function_3d &phi0,
                                          const vector_real_function_3d &rho1,
                                          bool update) const override {
    vector_real_function_3d rho_new;
    if (update) {
      rho_new = copy(world, rho1);
    } else {
      rho_new = zero_functions<double, 3>(world, x.num_states());
    }
    vector_real_function_3d x_phi, y_phi;
    for (const auto &b : x.active) {
      rho_new[b] = dot(world, x.x[b], phi0, true);
    }
    for (const auto &b : x.active) {
      rho_new[b] += dot(world, x.y[b], phi0, true);
    }

    truncate(world, rho_new);

    return rho_new;
  }
};

class VXC1Strategy {
public:
  virtual ~VXC1Strategy() = default;
  virtual X_space compute_VXC1(
      World &world, const X_space &x, const vector_real_function_3d &rho1,
      const std::pair<vector_real_function_3d, vector_real_function_3d> &phi0,
      const XCOperator<double, 3> &xc) const = 0;
};

class VXC1StrategyStandard : public VXC1Strategy {

public:
  X_space compute_VXC1(
      World &world, const X_space &x, const vector_real_function_3d &rho1,
      const std::pair<vector_real_function_3d, vector_real_function_3d> &phi0,
      const XCOperator<double, 3> &xc) const override {

    X_space W =
        X_space::zero_functions(world, x.num_states(), x.num_orbitals());

    auto compute_wx_x = [&, &phi0 = phi0](auto rho_alpha) {
      auto xc_rho = xc.apply_xc_kernel(rho_alpha);
      return mul(world, xc_rho, phi0.first);
    };
    auto compute_wx_y = [&, &phi0 = phi0](auto rho_alpha) {
      auto xc_rho = xc.apply_xc_kernel(rho_alpha);
      return mul(world, xc_rho, phi0.second);
    };
    std::transform(rho1.begin(), rho1.end(), W.x.begin(), compute_wx_x);
    std::transform(rho1.begin(), rho1.end(), W.y.begin(), compute_wx_y);
    return W;
  }
};

class J1Strategy {
public:
  virtual ~J1Strategy() = default;
  virtual X_space compute_J1(
      World &world, const X_space &x, const vector_real_function_3d &rho1,
      const std::pair<vector_real_function_3d, vector_real_function_3d> &phi0,
      const poperatorT &coulomb_ops) const = 0;
};

class J1StrategyFull : public J1Strategy {
public:
  X_space compute_J1(
      World &world, const X_space &x, const vector_real_function_3d &rho1,
      const std::pair<vector_real_function_3d, vector_real_function_3d> &phi0,
      const poperatorT &coulomb_ops) const override {

    X_space J =
        X_space::zero_functions(world, x.num_states(), x.num_orbitals());
    vector_real_function_3d temp_J(3);
    for (const auto &b : x.active) {
      temp_J[b] = apply(*coulomb_ops, rho1[b]);
      J.x[b] = mul(world, temp_J[b], phi0.first, false);
    }
    J.y = J.x.copy();
    return J;
  }
};

class J1StrategyStable : public J1Strategy {
public:
  X_space compute_J1(
      World &world, const X_space &x, const vector_real_function_3d &rho1,
      const std::pair<vector_real_function_3d, vector_real_function_3d> &phi0,
      const poperatorT &coulomb_ops) const override {

    X_space J =
        X_space::zero_functions(world, x.num_states(), x.num_orbitals());
    // if (world.rank() == 0) { print("J1StrategyStable"); }
    vector_real_function_3d temp_J(3);
    for (const auto &b : x.active) {
      temp_J[b] = apply(*coulomb_ops, rho1[b]);
      if (false) {
        auto norm = temp_J[b].norm2();
        if (world.rank() == 0)
          print("norm of temp_J:", norm);
      }
      J.x[b] = mul(world, temp_J[b], phi0.first, false);
      J.y[b] = mul(world, temp_J[b], phi0.second, false);
    }
    world.gop.fence();
    return J;
  }
};

class K1Strategy {
public:
  virtual ~K1Strategy() = default;
  virtual X_space compute_K1(
      World &world, const X_space &x,
      const std::pair<vector_real_function_3d, vector_real_function_3d> &phi_0X,
      const std::pair<vector_real_function_3d, vector_real_function_3d>
          &rhs_vec) const = 0;
  std::string algorithm_;

  [[nodiscard]] auto make_k(const vecfuncT &ket, const vecfuncT &bra) const {
    auto &world = ket[0].world();
    const double lo = 1.e-10;
    Exchange<double, 3> k{world, lo};
    k.set_bra_and_ket(bra, ket);
    if (algorithm_ == "multiworld") {
      k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient);
    } else if (algorithm_ == "multiworld_row") {
      k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient_row);
    } else if (algorithm_ == "largemem") {
      k.set_algorithm(Exchange<double, 3>::Algorithm::large_memory);
    } else if (algorithm_ == "smallmem") {
      k.set_algorithm(Exchange<double, 3>::Algorithm::small_memory);
    }
    return k;
  };
  virtual void set_algorithm(const std::string &algo) = 0;
};

class K1StrategyFull : public K1Strategy {
public:
  X_space compute_K1(
      World &world, const X_space &x,
      const std::pair<vector_real_function_3d, vector_real_function_3d> &phi_0X,
      const std::pair<vector_real_function_3d, vector_real_function_3d>
          &rhs_vec) const override {

    auto K = X_space::zero_functions(world, x.num_states(), x.num_orbitals());

    vector_real_function_3d k1x, k1y, k2x, k2y;
    vector_real_function_3d xb;
    vector_real_function_3d yb;

    auto k1x_temp = create_response_matrix(x.num_states(), x.num_orbitals());
    auto k1y_temp = create_response_matrix(x.num_states(), x.num_orbitals());
    auto k2x_temp = create_response_matrix(x.num_states(), x.num_orbitals());
    auto k2y_temp = create_response_matrix(x.num_states(), x.num_orbitals());

    for (const auto &b : x.active) {
      auto K1X = make_k(x.x[b], phi_0X.first);
      auto K1Y = make_k(phi_0X.first, x.y[b]);

      auto K2X = make_k(x.y[b], phi_0X.first);
      auto K2Y = make_k(phi_0X.second, x.x[b]);
      world.gop.fence();

      k1x_temp[b] = K1X(rhs_vec.first);
      k1y_temp[b] = K1Y(rhs_vec.first);
      k2x_temp[b] = K2X(rhs_vec.second);
      k2y_temp[b] = K2Y(rhs_vec.second);
      world.gop.fence();
      K.x[b] = gaxpy_oop(1.0, k1x_temp[b], 1.0, k1y_temp[b], false);
      K.y[b] = gaxpy_oop(1.0, k2x_temp[b], 1.0, k2y_temp[b], false);
    }
    world.gop.fence();
    return K;
  }
  void set_algorithm(const std::string &algo) override { algorithm_ = algo; }
};

class K1StrategyStatic : public K1Strategy {
public:
  X_space compute_K1(
      World &world, const X_space &x,
      const std::pair<vector_real_function_3d, vector_real_function_3d> &phi_0X,
      const std::pair<vector_real_function_3d, vector_real_function_3d>
          &rhs_vec) const override {
    X_space K =
        X_space::zero_functions(world, x.num_states(), x.num_orbitals());
    vector_real_function_3d k1x, k1y, k2x, k2y;
    const double lo = 1e-10;

    auto k1_temp = create_response_matrix(x.num_states(), x.num_orbitals());
    auto k2_temp = create_response_matrix(x.num_states(), x.num_orbitals());

    for (const auto &b : x.active) {
      auto K1Xs = make_k(x.x[b], phi_0X.first);
      auto K1Ys = make_k(phi_0X.first, x.x[b]);
      k1_temp[b] = K1Xs(rhs_vec.first);
      k2_temp[b] = K1Ys(rhs_vec.first);
      K.x[b] = gaxpy_oop(1.0, k1_temp[b], 1.0, k2_temp[b], false);
    }
    world.gop.fence();
    return K;
  }
  void set_algorithm(const std::string &algo) override { algorithm_ = algo; }
};

class inner_strategy {

public:
  virtual ~inner_strategy() = default;
  [[nodiscard]] virtual Tensor<double>
  compute_inner(const X_space &x, const X_space &y) const = 0;
};

class Context {

private:
  std::unique_ptr<inner_strategy> inner_strategy_;
  std::unique_ptr<J1Strategy> j1_strategy_;
  std::unique_ptr<K1Strategy> k1_strategy_;
  std::unique_ptr<VXC1Strategy> vxc1_strategy_;
  std::unique_ptr<ComputeDensityStrategy> density_strategy_;
  std::unique_ptr<LoadXSpaceStrategy> load_x_space_strategy_;

public:
  explicit Context(
      std::unique_ptr<inner_strategy> &&innerStrategy = {},
      std::unique_ptr<J1Strategy> &&j1Strategy = {},
      std::unique_ptr<K1Strategy> &&k1Strategy = {},
      std::unique_ptr<VXC1Strategy> &&vxc1trategy = {},
      std::unique_ptr<ComputeDensityStrategy> &&densityStrategy = {},
      std::unique_ptr<LoadXSpaceStrategy> &&loadXSpaceStrategy = {})
      : inner_strategy_(std::move(innerStrategy)),
        j1_strategy_(std::move(j1Strategy)),
        k1_strategy_(std::move(k1Strategy)),
        vxc1_strategy_(std::move(vxc1trategy)),
        density_strategy_(std::move(densityStrategy)),
        load_x_space_strategy_(std::move(loadXSpaceStrategy)) {}
  void set_strategy(std::unique_ptr<inner_strategy> &&strategy,
                    std::unique_ptr<J1Strategy> &&j1Strategy,
                    std::unique_ptr<K1Strategy> &&K1Strategy,
                    std::unique_ptr<VXC1Strategy> &&vxc1Strategy,
                    std::unique_ptr<ComputeDensityStrategy> &&densityStrategy,
                    std::unique_ptr<LoadXSpaceStrategy> &&loadXSpaceStrategy,
                    const ResponseParameters &r_params) {
    inner_strategy_ = std::move(strategy);
    j1_strategy_ = std::move(j1Strategy);
    k1_strategy_ = std::move(K1Strategy);
    vxc1_strategy_ = std::move(vxc1Strategy);
    density_strategy_ = std::move(densityStrategy);
    load_x_space_strategy_ = std::move(loadXSpaceStrategy);

    k1_strategy_->set_algorithm(r_params.hfexalg());
  }

  [[nodiscard]] Tensor<double> inner(const X_space &x, const X_space &y) const {
    if (inner_strategy_) {
      return inner_strategy_->compute_inner(x, y);
    } else {
      throw madness::MadnessException("Inner product Strategy isn't set",
                                      "Need to set a strategy", 2, 455, "inner",
                                      "ResponseBase.hpp");
    }
  }

  X_space compute_j1(World &world, const X_space &x,
                     const vector_real_function_3d &rho1,
                     const double_response_vector &phi0,
                     const poperatorT &coulomb_ops) const {
    if (j1_strategy_) {
      return j1_strategy_->compute_J1(world, x, rho1, phi0, coulomb_ops);
    } else {
      throw madness::MadnessException("Compute J1 Strategy isn't set",
                                      "Need to set a strategy", 2, 455, "inner",
                                      "ResponseBase.hpp");
    }
  }

  X_space compute_k1(World &world, const X_space &x,
                     const double_response_vector &phi0X,
                     const double_response_vector &phi0) const {
    if (k1_strategy_) {
      return k1_strategy_->compute_K1(world, x, phi0X, phi0);
    } else {
      throw madness::MadnessException("Compute K1 Strategy isn't set",
                                      "Need to set a strategy", 2, 455, "inner",
                                      "ResponseBase.hpp");
    }
  }
  X_space compute_VXC1(World &world, const X_space &x,
                       const vector_real_function_3d &rho1,
                       const double_response_vector &phi0,
                       const XCOperator<double, 3> &xc) const {
    if (vxc1_strategy_) {
      return vxc1_strategy_->compute_VXC1(world, x, rho1, phi0, xc);
    } else {
      throw madness::MadnessException("Compute VXC1 Strategy isn't set",
                                      "Need to set a strategy", 2, 455, "inner",
                                      "ResponseBase.hpp");
    }
  }

  vector_real_function_3d compute_density(World &world, const X_space &x,
                                          const vector_real_function_3d &phi0,
                                          const vector_real_function_3d &rho1,
                                          bool update) const {
    if (density_strategy_) {
      return density_strategy_->compute_density(world, x, phi0, rho1, update);
    } else {
      throw madness::MadnessException("Compute Density Strategy isn't set",
                                      "Need to set a strategy", 2, 455, "inner",
                                      "ResponseBase.hpp");
    }
  }
  XData load_x_space(World &world, const std::string &filename,
                     const ResponseParameters &r_params,
                     double omega_state) const {
    if (load_x_space_strategy_) {
      return load_x_space_strategy_->load_x_space(world, filename, r_params,
                                                  omega_state);
    } else {
      throw madness::MadnessException("Load X Space Strategy isn't set",
                                      "Need to set a strategy", 2, 455, "inner",
                                      "ResponseBase.hpp");
    }
  }
};

class full_inner_product : public inner_strategy {
public:
  [[nodiscard]] Tensor<double> compute_inner(const X_space &x,
                                             const X_space &y) const override {
    return inner(x, y);
  }
};

class static_inner_product : public inner_strategy {
public:
  [[nodiscard]] Tensor<double> compute_inner(const X_space &x,
                                             const X_space &y) const override {
    return response_space_inner(x.x, y.x);
  }
};
typedef std::vector<XNonlinearSolver<vector_real_function_3d, double,
                                     response_matrix_allocator>>
    response_solver;
typedef std::vector<
    XNonlinearSolver<real_function_3d, double, response_function_allocator>>
    response_function_solver;

class response_timing {
  std::map<std::string, std::vector<double>> wall_time_data;
  std::map<std::string, std::vector<double>> cpu_time_data;
  int iter;

public:
  response_timing();

  void to_json(json &j);

  void print_data();

  void add_data(std::map<std::string, std::pair<double, double>> values);
};
class response_data {
  int iter;
  std::vector<double> thresh;
  std::vector<double> density_target;
  std::vector<double> bsh_target;

public:
  void to_json(json &j);

  void add_data(std::map<std::string, Tensor<double>> values);
  void add_convergence_targets(double p_thresh, double p_density_target,
                               double p_bsh_target);
};

class ResponseTester;

struct residuals {
  X_space residual;
  Tensor<double> residual_norms;
};

using gamma_orbitals =
    std::tuple<X_space, vector_real_function_3d, vector_real_function_3d>;

class ResponseBase {
public:
  friend ResponseTester;

  ResponseBase(World &world, const CalcParams &params);

  void solve(World &world);

  virtual void initialize(World &world) = 0;

  virtual void iterate(World &world) = 0;

  // virtual void iterate();
  auto get_parameter() const -> CalcParams {
    return {ground_calc, molecule, r_params};
  }

  Molecule get_molecule() const;
  auto get_orbitals() const -> vector_real_function_3d {
    return ground_orbitals;
  }

  auto get_chi() const -> X_space { return Chi.copy(); };

  void output_json();

  json j_molresponse{};
  response_timing time_data;
  response_data function_data;
  mutable std::map<std::string, std::pair<double, double>> iter_timing;

  Context response_context;

  static void help() {
    print_header2("help page for MOLRESPONSE ");
    print("The molresponse code computes linear response properties and "
          "excited states");
    print("\nYou can print all available calculation parameters by "
          "running\n");
    print("molresponse --print_parameters\n");
    print("You can perform a simple calculation by running\n");
    print("moldft --geometry=h2o.xyz");
    print("molresponse");
    print("\nprovided you have an xyz file in your directory as well as a "
          "file named 'response.in'.");
    print("with minimal input\n");
    print("response ");
    print("  archive mad.restartdata ");
    print("  excited_state 1 ");
    print("end ");
  }

  static void print_parameters() {
    ResponseParameters rparam;
    print("A molresponse calculation requires a converged moldft "
          "calculations with the ");
    print("corresponding parameters.");
    print("\nDefault parameters for the response part of the molresponse "
          "program are");
    rparam.print("response", "end");
    print("\n\nthe molecular geometry must be specified in a separate "
          "block:");
    Molecule::print_parameters();
  }

  void write_vtk(World &world, int num_points, const double &L,
                 const std::string &filename);

protected:
  // Given molecule returns the nuclear potential of the molecule
  ResponseParameters r_params;
  Molecule molecule;
  GroundStateCalculation ground_calc;
  bool all_done = false;

  XCfunctional xcf;
  real_function_3d mask;

  std::shared_ptr<PotentialManager> potential_manager;
  // shared pointers to Operators
  poperatorT shared_coulomb_operator; // shared pointer to seperated convolution
                                      // operator
  std::vector<std::shared_ptr<real_derivative_3d>> gradop;
  // Stored functions
  mutable real_function_3d
      stored_v_nuc; // Stored nuclear potential from ground state
  mutable real_function_3d
      stored_v_coul; // Stored coulomb potential from ground state

  // Ground state orbitals and energies
  vector_real_function_3d ground_orbitals{};
  Tensor<double> ground_energies;

  // Information that is inferred from input file
  // Ground state orbitals being used in calculation
  Tensor<double> hamiltonian;
  Tensor<double> ham_no_diag;
  // Tensors for holding energies
  // residuals, and shifts

  Tensor<double> e_residuals; // Residuals of energies

  // Mask function to handle boundary conditions

  functionT ground_density; // ground state density

  mutable response_space
      stored_potential; // The ground state potential, stored only
  // if store_potential is true (default is

  double vtol{};

  X_space Chi;
  /// Sets the Function protocol dependent on the truncation threshold.
  /// Sets the polynomial order of basis functions k
  /// Then creates shared coulomb operator, gradient operator, ground density
  /// AS well as the ground state density
  /// \param world
  /// \param thresh
  void set_protocol(World &world, double thresh) {
    int k;
    // Allow for imprecise conversion of threshold
    if (thresh >= 0.9e-2)
      k = 4;
    else if (thresh >= 0.9e-4)
      k = 6;
    else if (thresh >= 0.9e-6)
      k = 8;
    else if (thresh >= 0.9e-8)
      k = 10;
    else
      k = 12;

    // k defaults to make sense with thresh, override by providing k in
    // input file
    if (r_params.k() == -1) {
      FunctionDefaults<3>::set_k(k);
    } else {
      FunctionDefaults<3>::set_k(r_params.k());
    }

    // Set Function Defaults
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_initial_level(2);

    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_apply_randomize(false);

    FunctionDefaults<3>::set_project_randomize(false);
    GaussianConvolution1DCache<double>::map.clear();

    double safety = 0.1;
    vtol = FunctionDefaults<3>::get_thresh() * safety;
    shared_coulomb_operator =
        poperatorT(CoulombOperatorPtr(world, r_params.lo(), thresh));
    gradop = gradient_operator<double, 3>(world);
    potential_manager = std::make_shared<PotentialManager>(molecule, "a");
    potential_manager->make_nuclear_potential(world);
    // Create the masking function
    mask = real_function_3d(
        real_factory_3d(world).f(mask3).initial_level(4).norefine());

    ground_density = make_ground_density(world);
    ground_density.truncate(FunctionDefaults<3>::get_thresh());
    // Basic print
    if (world.rank() == 0) {
      print("\nSolving NDIM=", 3, " with thresh", thresh, "    k",
            FunctionDefaults<3>::get_k(), "  dconv",
            std::max(thresh, r_params.dconv()), "\n");
    }
  }

  virtual void check_k(World &world, double thresh, int k);

  auto make_ground_density(World &world) const -> functionT;

  auto ComputeHamiltonianPair(World &world) const
      -> std::pair<Tensor<double>, Tensor<double>>;

  auto Coulomb(World &world) const -> real_function_3d;

  auto make_xc_operator(World &world) const -> XCOperator<double, 3>;

  void save_x_space(World &world, const std::string &name, const X_space &X) {
    // Archive to write everything to
    archive::ParallelOutputArchive ar(world, name.c_str(), 1);

    ar & X.num_orbitals();
    ar & X.num_states();

    for (size_t i = 0; i < X.num_states(); i++)
      for (size_t j = 0; j < X.num_orbitals(); j++)
        ar & X.x[i][j];
    for (size_t i = 0; i < X.num_states(); i++)
      for (size_t j = 0; j < X.num_orbitals(); j++)
        ar & X.y[i][j];
  }

  // Load a response calculation
  X_space load_x_space(World &world, const std::string &name) {
    archive::ParallelInputArchive ar(world, name.c_str());
    size_t num_orbitals;
    size_t num_states;

    ar & num_orbitals;
    ar & num_states;
    auto result = X_space(world, num_states, num_orbitals);
    for (size_t i = 0; i < num_states; i++)
      for (size_t j = 0; j < num_orbitals; j++)
        ar & result.x[i][j];
    world.gop.fence();
    for (size_t i = 0; i < num_states; i++)
      for (size_t j = 0; j < num_orbitals; j++)
        ar & result.y[i][j];
    world.gop.fence();
    return result;
  }

  virtual void save(World &world, const std::string &name) = 0;

  virtual void load(World &world, const std::string &name) = 0;

  auto make_density(World &world, const X_space &chi) const -> vecfuncT;

  void load_balance_chi(World &world);

  auto
  make_bsh_operators_response(World &world, double &shift,
                              const double omega) const -> vector<poperatorT>;

  auto kain_x_space_update(World &world, const X_space &chi,
                           const X_space &residual_chi,
                           response_solver &kain_x_space) -> X_space;

  void x_space_step_restriction(World &world, const X_space &old_Chi,
                                X_space &temp, bool restrict_y,
                                const double &max_bsh_rotation);

  //    void plotResponseOrbitals(World &world, size_t iteration,
  //                              const response_space &x_response,
  //                              const response_space &y_response,
  //                              const ResponseParameters &responseParameters,
  //                              const GroundStateCalculation &g_params);

  static auto orbital_load_balance(World &world, const gamma_orbitals &,
                                   double load_balance) -> gamma_orbitals;

  auto compute_gamma_tda(World &world, const gamma_orbitals &density,
                         const XCOperator<double, 3> &xc) const -> X_space;

  auto compute_gamma_static(World &world, const gamma_orbitals &,
                            const XCOperator<double, 3> &xc) const -> X_space;

  auto compute_gamma_full(World &world, const gamma_orbitals &,
                          const XCOperator<double, 3> &xc) const -> X_space;
  auto compute_gamma(World &world, const gamma_orbitals &,
                     const XCOperator<double, 3> &xc) const -> X_space;

  auto compute_V0X(World &world, const X_space &X,
                   const XCOperator<double, 3> &xc,
                   bool compute_Y) const -> X_space;

  auto compute_lambda_X(World &world, const X_space &chi,
                        XCOperator<double, 3> &xc,
                        const std::string &calc_type) const -> X_space;

  auto compute_theta_X(World &world, const X_space &chi,
                       const vector_real_function_3d &rho1,
                       const XCOperator<double, 3> &xc,
                       const std::string &calc_type) const -> X_space;

  auto compute_F0X(World &world, const X_space &X,
                   const XCOperator<double, 3> &xc,
                   bool compute_Y) const -> X_space;

  void analyze_vectors(World &world, const vecfuncT &x,
                       const std::string &response_state);

  auto project_ao_basis(World &world,
                        const AtomicBasisSet &aobasis) -> vecfuncT;

  static auto project_ao_basis_only(World &world, const AtomicBasisSet &aobasis,
                                    const Molecule &mol) -> vecfuncT;

  void converged_to_json(json &j);

  auto update_residual(World &world, const X_space &chi, const X_space &g_chi,
                       const std::string &calc_type,
                       const Tensor<double> &old_residuals,
                       const X_space &xres_old) -> residuals;

  auto compute_response_potentials(World &world, const X_space &chi,
                                   XCOperator<double, 3> &xc,
                                   const std::string &calc_type) const
      -> std::tuple<X_space, X_space, X_space>;

  // compute exchange |i><i|J|p>
  auto exchangeHF(const vecfuncT &ket, const vecfuncT &bra,
                  const vecfuncT &vf) const -> vecfuncT {
    World &world = ket[0].world();
    auto n = bra.size();
    auto nf = ket.size();
    double tol = FunctionDefaults<3>::get_thresh(); /// Important this is
    double mul_tol = 0.0;
    const double lo = r_params.lo();

    std::shared_ptr<real_convolution_3d> poisson;
    /// consistent with Coulomb
    vecfuncT Kf = zero_functions_compressed<double, 3>(world, nf, true);

    reconstruct(world, bra);
    reconstruct(world, ket);
    reconstruct(world, vf);

    // i-j sym
    for (int i = 0; i < n; ++i) {
      // for each |i> <i|phi>
      vecfuncT psi_f =
          mul_sparse(world, bra[i], vf, mul_tol, true); /// was vtol
      truncate(world, psi_f, tol, true);
      // apply to vector of products <i|phi>..<i|1> <i|2>...<i|N>
      psi_f = apply(world, *shared_coulomb_operator, psi_f);
      truncate(world, psi_f, tol, true);
      // multiply by ket i  <i|phi>|i>: <i|1>|i> <i|2>|i> <i|2>|i>
      psi_f = mul_sparse(world, ket[i], psi_f, mul_tol, true); /// was vtol
      /// Generalized A*X+y for vectors of functions ---- a[i] = alpha*a[i] +
      // 1*Kf+occ[i]*psi_f
      gaxpy(world, double(1.0), Kf, double(1.0), psi_f);
    }
    truncate(world, Kf, tol, true);
    return Kf;
  }

  void print_inner(World &world, const std::string &name, const X_space &left,
                   const X_space &right) const;

  void function_data_to_json(json &j_mol_in, size_t iter,
                             const Tensor<double> &x_norms,
                             const Tensor<double> &x_abs_norms,
                             const Tensor<double> &rho_norms,
                             const Tensor<double> &rho_res_norms);
  X_space compute_TX(World &world, const X_space &X, bool compute_Y) const;
  vecfuncT update_density(World &world, const X_space &chi,
                          const vecfuncT &old_density) const;
};

// Some helper functions
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Check k given response parameters
/// \param world
/// \param Chi
/// \param thresh
/// \param k
void check_k(World &world, X_space &Chi, double thresh, int k);

auto add_randomness(World &world, const response_space &f,
                    double magnitude) -> response_space;

void normalize(World &world, response_space &f);

void normalize(World &world, X_space &Chi);

static auto kronecker(size_t l, size_t n) -> double {
  if (l == n)
    return 1.0;
  return 0.0;
}

auto solid_harmonics(World &world,
                     int n) -> std::map<std::vector<int>, real_function_3d>;

/***
 * @brief Prints the norms of the functions of a response space
 *
 *
 *
 * @param world
 * @param f
 */
void print_norms(World &world, const response_space &f);

// Returns a list of solid harmonics such that:
// solid_harm.size() * num_ground_orbs > 2 * num. resp. components
auto make_xyz_functions(World &world) -> vector_real_function_3d;

// Selects from a list of functions and energies the k functions with the
// lowest energy
auto select_functions(World &world, response_space f, Tensor<double> &energies,
                      size_t k, size_t print_level) -> response_space;

// Sorts the given tensor of eigenvalues and
// response functions
void sort(World &world, Tensor<double> &vals, response_space &f);

void sort(World &world, Tensor<double> &vals, X_space &f);

// Specialized for response calculations that returns orthonormalized
// functions
auto gram_schmidt(World &world, const response_space &f) -> response_space;

/// Computes the transition density between set of two response functions x and
/// y. Uses std::transform to iterate between x and y vectors \param world
/// \param orbitals
/// \param x
/// \param y
/// \return
auto transition_density(World &world, const vector_real_function_3d &orbitals,
                        const response_space &x,
                        const response_space &y) -> vector_real_function_3d;

auto transition_densityTDA(World &world,
                           const vector_real_function_3d &orbitals,
                           const response_space &x) -> vector_real_function_3d;

auto transform(World &world, const response_space &f,
               const Tensor<double> &U) -> response_space;

auto transform(World &world, const X_space &x,
               const Tensor<double> &U) -> X_space;

// result(i,j) = inner(a[i],b[j]).sum()
auto expectation(World &world, const response_space &A,
                 const response_space &B) -> Tensor<double>;

void inner_to_json(World &world, const std::string &name,
                   const Tensor<double> &m_val,
                   std::map<std::string, Tensor<double>> &data);
class ResponseTester {

public:
  static void load_calc(World &world, ResponseBase *p, double thresh) {
    p->set_protocol(world, thresh);
    p->load(world, p->r_params.restart_file());
    p->check_k(world, thresh, FunctionDefaults<3>::get_k());
  }

  static X_space compute_gamma_full(World &world, ResponseBase *p,
                                    double thresh) {
    XCOperator<double, 3> xc = p->make_xc_operator(world);
    return X_space{};
  }

  X_space compute_lambda_X(World &world, ResponseBase *p, double thresh) {
    XCOperator<double, 3> xc = p->make_xc_operator(world);
    X_space gamma =
        p->compute_lambda_X(world, p->Chi, xc, p->r_params.calc_type());
    return gamma;
  }

  std::pair<X_space, X_space> compute_VFOX(World &world, ResponseBase *p,
                                           bool compute_y) {
    XCOperator<double, 3> xc = p->make_xc_operator(world);
    X_space V = p->compute_V0X(world, p->Chi, xc, compute_y);
    X_space F = p->compute_F0X(world, p->Chi, xc, compute_y);
    return {V, F};
  }
};

#endif // MADNESS_RESPONSEBASE_HPP
