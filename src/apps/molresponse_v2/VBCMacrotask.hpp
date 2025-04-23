#pragma once
#include "FrequencyLoop.hpp"
#include "GroundStateData.hpp"
#include "ResponseIO.hpp"
#include "ResponseState.hpp"
#include "ResponseVector.hpp"
#include "SCF.h"
#include "funcdefaults.h"
#include "functypedefs.h"
#include <madness/mra/macrotaskq.h>
#include <madness/world/cloud.h>

namespace madness {

class VBC_save_load : public MacroTaskOperationBase {

  class Partitioner : public MacroTaskPartitioner {
  public:
    explicit Partitioner() {
      max_batch_size = 1;
      policy = "guided";
    }
  };

public:
  GroundStateData &gs;
  double thresh;
  vector_real_function_3d dipole_vector;
  std::map<char, real_function_3d> dipole_vector_map;
  int num_directions;
  int num_orbitals;
  int num_frequencies;
  vector<int> bfreq_index; // frequency index for state B
  vector<int> cfreq_index; // frequency index for state C
  vector<int> bc_index;    // IJ pair for BC
  std::vector<std::pair<char, char>> BC_pairs;
  std::map<char, ResponseState> state_map;
  std::map<char, vector_real_function_3d> perturbation_map;

  std::string &directions;
  std::vector<ResponseState> active_states;

  std::vector<double> &frequencies;
  bool isSpinRestricted;
  explicit VBC_save_load(GroundStateData &gs, std::string directions,
                         std::vector<double> frequencies, bool isSpinRestricted)
      : gs(gs), directions(directions), frequencies(frequencies),
        isSpinRestricted(isSpinRestricted) {
    partitioner.reset(new Partitioner());
    num_directions = static_cast<int>(directions.size());
    num_orbitals = gs.getNumOrbitals();
    thresh = FunctionDefaults<3>::get_thresh();

    World &world = gs.orbitals[0].world();
    auto num_orbitals = static_cast<int>(gs.orbitals.size());
    auto thresh = FunctionDefaults<3>::get_thresh();

    auto num_directions = static_cast<int>(directions.size());
    std::vector<ResponseState> active_states(num_directions);
    std::vector<vector_real_function_3d> VP(num_directions);
    std::vector<ResponseVector> active_vec(num_directions);

    if (world.rank() == 0) {
      std::cout << "Computing polarizability for frequencies: ";
      for (const auto &freq : frequencies) {
        std::cout << freq << " ";
      }
      std::cout << "\n";
    }
    vector_real_function_3d dipole_vectors(3);
    size_t k = 0;
    // creates a vector of x y z dipole functions
    for (auto &d : dipole_vectors) {
      std::vector<int> f(3, 0);
      f[k++] = 1;
      d = real_factory_3d(world).functor(real_functor_3d(new MomentFunctor(f)));
    }

    for (int i = 0; i < num_directions; i++) {
      auto dir = directions[i];
      dipole_vector_map[dir] = dipole_vectors[i];
    }

    truncate(world, dipole_vectors, FunctionDefaults<3>::get_thresh(), true);

    for (int d = 0; d < directions.size(); d++) {
      if (world.rank() == 0) {
        std::cout << "Computing polarizability for direction: " << directions[d]
                  << "\n";
      }
      DipolePerturbation pert{directions[d]};

      ResponseState state(pert, PerturbationType::Dipole, frequencies,
                          {FunctionDefaults<3>::get_thresh()},
                          isSpinRestricted);
      if (world.rank() == 0) {
        std::cout << "Perturbation: " << state.perturbationDescription()
                  << "\n";
      }
      state_map[directions[d]] = state;
      perturbation_map[directions[d]] = state.perturbation_vector(world, gs);
    }

    // all BC required during one frequency pair
    for (int i = 0; i < num_directions; i++) {
      for (int j = i + 1; j < num_directions; j++) {
        BC_pairs.push_back({directions[i], directions[j]});
      }
    }
    num_frequencies = static_cast<int>(frequencies.size());

    for (int b = 0; b < num_frequencies; b++) {
      for (int c = b; c < num_frequencies; b++) {
        for (int bc = 0; bc < BC_pairs.size(); bc++) {
          bfreq_index.push_back(b);
          cfreq_index.push_back(c);
          bc_index.push_back(bc);
        }
      }
    }

    // setting up the vbc tasks// print

    if (world.rank() == 0) {
      std::cout << "Number of BC pairs: " << BC_pairs.size() << "\n";
      std::cout << "Number of frequency pairs: " << num_frequencies << "\n";
      std::cout << "Number of directions: " << num_directions << "\n";
      std::cout << "Number of orbitals: " << num_orbitals << "\n";
      std::cout << "Number of VBC tasks: " << bfreq_index.size() << "\n";
    }
  }

  typedef std::tuple<const std::vector<int> &> argtupleT;

  using resultT = real_function_3d;

  static auto K(const vector_real_function_3d &bra,
                const vector_real_function_3d &ket) {
    auto &world = bra[0].world();
    Exchange<double, 3> k{world, 1.e-10};
    k.set_bra_and_ket(bra, ket);
    return k;
  }
  static vector_real_function_3d
  make_zeta_bc(World &world, const vector_real_function_3d &by,
               const vector_real_function_3d &cx,
               const vector_real_function_3d &ground_orbitals) {
    auto matrix_bc = matrix_inner(world, by, cx);
    return -1.0 * transform(world, ground_orbitals, matrix_bc, true);
  };

  static auto compute_g(World &world, const vector_real_function_3d &Aleft,
                        const vector_real_function_3d &Aright,
                        const vector_real_function_3d &Bleft,
                        const vector_real_function_3d &Bright,
                        const vector_real_function_3d &phix,
                        const vector_real_function_3d &phiy) {
    auto x_phi = mul(world, Aleft, Aright, true);
    auto y_phi = mul(world, Bleft, Bright, true);
    world.gop.fence();
    auto rho = sum(world, x_phi, true);
    world.gop.fence();
    rho += sum(world, y_phi, true);
    world.gop.fence();
    auto lo = 1.e-10;
    real_convolution_3d op =
        CoulombOperator(world, lo, FunctionDefaults<3>::get_thresh());
    auto temp_J = apply(op, rho);

    auto Jx = mul(world, temp_J, phix, false);
    auto Jy = mul(world, temp_J, phiy, false);

    auto ka = K(Aleft, Aright)(phix); // what happens to k after this?
    auto kb = K(Bleft, Bright)(phix);
    auto ka_conj = K(Aright, Aleft)(phiy);
    auto kb_conj = K(Bright, Bleft)(phiy);
    world.gop.fence();
    // ideally it runs and the Exchange operator is freed

    auto Kx = gaxpy_oop(1.0, ka, 1.0, kb, true);
    auto Ky = gaxpy_oop(1.0, ka_conj, 1.0, kb_conj, true);
    world.gop.fence();

    auto result_x = gaxpy_oop(2.0, Jx, -1.0, Kx, true);
    auto result_y = gaxpy_oop(2.0, Jy, -1.0, Ky, true);

    return std::make_pair(result_x, result_y);
  };
  // This function constructs the J and K operators with A and B and applies
  static auto compute_vbc_i(World &world, const vector_real_function_3d &bx,
                            const vector_real_function_3d &by,
                            const vector_real_function_3d &cx,
                            const vector_real_function_3d &cy,
                            const vector_real_function_3d &zeta_bc,
                            const vector_real_function_3d &phi0,
                            const real_function_3d &v) {
    // left right pairs ka and kb and apply apply
    madness::QProjector<double, 3> Q(phi0);
    auto [gzeta_x, gzeta_y] =
        compute_g(world, bx, cy, phi0, zeta_bc, phi0, phi0);
    gzeta_x = -1.0 * Q(gzeta_x);
    gzeta_y = -1.0 * Q(gzeta_y);

    auto [gbc_x, gbc_y] = compute_g(world, bx, phi0, phi0, by, cx, cy);
    auto [gbc_phi_x, gbc_phi_y] =
        compute_g(world, bx, phi0, phi0, by, phi0, phi0);

    auto vcx = mul(world, v, cx, true);
    auto vcy = mul(world, v, cy, true);

    auto fbx = -1.0 * Q(gbc_x + vcx);
    auto fby = -1.0 * Q(gbc_y + vcy);

    auto vb_phi = mul(world, v, phi0, true);

    auto fb_phi_x = gbc_phi_x + vb_phi;
    auto fb_phi_y = gbc_phi_y + vb_phi;

    auto m_fbx = matrix_inner(world, phi0, fb_phi_x);
    auto m_fby = matrix_inner(world, phi0, fb_phi_y);

    auto fphi_x = transform(world, cx, m_fbx, true);
    auto fphi_y = transform(world, cy, m_fby, true);

    auto result_x = truncate(gzeta_x + fbx + fphi_x,
                             FunctionDefaults<3>::get_thresh(), true);
    auto result_y = truncate(gzeta_y + fby + fphi_y,
                             FunctionDefaults<3>::get_thresh(), true);

    return std::make_pair(result_x, result_y);
  };

  static resultT allocator(World &world, const argtupleT &args) {
    const auto &index = std::get<0>(args);
    return real_factory_3d(world).compressed();
  };

  resultT operator()(const std::vector<int> &index) {

    const int ij = static_cast<int>(batch.result.begin);
    World &world = gs.orbitals[0].world();

    // Here we figure out which BC pair we are going to create and save
    auto b = bfreq_index[ij];
    auto c = cfreq_index[ij];
    auto bc = bc_index[ij];

    auto dirb = BC_pairs[bc].first;
    auto dirc = BC_pairs[bc].second;
    if (world.rank() == 0) {
      print("ij: ", ij, " b: ", b, " c: ", c, " bc: ", bc, " dirb: ", dirb,
            " dirc: ", dirc);
    }

    DipolePerturbation perturbation_b =
        std::get<DipolePerturbation>(state_map[dirb].perturbation);
    DipolePerturbation perturbation_c = std::get<DipolePerturbation>(
        state_map[BC_pairs[bc].second].perturbation);

    auto state_b = state_map[dirb];
    auto state_c = state_map[dirc];
    state_b.set_frequency_index(b);
    state_c.set_frequency_index(c);

    auto omega_b = frequencies[b];
    auto omega_c = frequencies[c];

    SecondOrderResponseState vbc_state(
        perturbation_b, perturbation_c, omega_b, omega_c,
        FunctionDefaults<3>::get_thresh(), isSpinRestricted);

    ResponseVector load_b;
    ResponseVector load_c;
    // load response vectors for freq b
    load_response_vector(world, num_orbitals, state_map[dirb], load_b, 0, b);
    // load response vectors for freq c
    load_response_vector(world, num_orbitals, state_map[dirc], load_c, 0, c);

    // Always work with dynamic, therefore, we need to promote any static
    // Write now it should always be spin restricted
    auto state_vec_b =
        make_response_vector(num_orbitals, false, isSpinRestricted);
    auto state_vec_c =
        make_response_vector(num_orbitals, false, isSpinRestricted);

    promote_response_vector(world, load_b, state_vec_b);
    promote_response_vector(world, load_c, state_vec_c);

    // now we work with state_vec_b and state_vec_c
    // get x_alpha and y_alpha
    if (isSpinRestricted) {
      auto &responseb = std::get<DynamicRestrictedResponse>(state_vec_b);
      auto bx = responseb.x_alpha;
      auto by = responseb.y_alpha;

      auto &responsec = std::get<DynamicRestrictedResponse>(state_vec_c);
      auto cx = responsec.x_alpha;
      auto cy = responsec.y_alpha;

      auto zeta_bc = make_zeta_bc(world, by, cx, gs.orbitals);
      auto zeta_cb = make_zeta_bc(world, cy, bx, gs.orbitals);
      auto vb = perturbation_map[dirb];
      auto vc = perturbation_map[dirc];
      auto compute_result = [&](const int &i) {
        auto results =
            zero_functions_compressed<double, 3>(world, 2 * num_orbitals);

        auto [bcx, bcy] = compute_vbc_i(world, bx, by, cx, cy, zeta_bc,
                                        gs.orbitals, dipole_vector_map[dirb]);
        auto [cbx, cby] = compute_vbc_i(world, cx, cy, bx, by, zeta_cb,
                                        gs.orbitals, dipole_vector_map[dirc]);
        auto result_vector =
            make_response_vector(num_orbitals, false, isSpinRestricted);
        auto result = std::get<DynamicRestrictedResponse>(result_vector);

        auto &result_x =
            std::get<DynamicRestrictedResponse>(result_vector).x_alpha;
        auto &result_y =
            std::get<DynamicRestrictedResponse>(result_vector).y_alpha;

        result_x = gaxpy_oop(1.0, bcx, 1.0, cbx, true);
        result_y = gaxpy_oop(1.0, bcy, 1.0, cby, true);

        for (int j = 0; j < num_orbitals; j++) {

          auto index_x = j;
          auto norm_x = result_x[j].norm2();
          print("i,j = (", i, j, ") index_x: ", index_x, " norm_x: ", norm_x);
        }

        for (int j = 0; j < num_orbitals; j++) {

          auto index_y = j + num_orbitals;
          auto norm_y = result_y[j].norm2();
          print("i,j = (", i, j, ") index_y: ", index_y, " norm_y: ", norm_y);
        }
        result.sync();
        return result_vector;
      };
      auto vbc = compute_result(ij);
      save_response_vector(world, vbc_state, vbc);
    }
    return real_factory_3d(world).compressed();
  }
};

class VBCComputer {
public:
  VBCComputer(World &world, const GroundStateData &gs,
              const std::vector<double> &frequencies,
              const std::string &directions, bool spin_restricted)
      : world_(world), gs_(gs), directions_(directions),
        frequencies_(frequencies), spin_restricted_(spin_restricted),
        num_orbitals_(gs.getNumOrbitals()) {

    for (char dir : directions_) {
      DipolePerturbation pert{dir};
      ResponseState state(pert, PerturbationType::Dipole, frequencies_,
                          {FunctionDefaults<3>::get_thresh()},
                          spin_restricted_);
      if (world.rank() == 0) {
        print("Perturbation: ", state.perturbationDescription());
      }
      state_map_[dir] = state;
      perturbation_vectors_[dir] = state.perturbation_vector(world_, gs_);
    }

    // Generate all BC direction pairs (only unique)
    for (size_t i = 0; i < directions.size(); ++i) {
      for (size_t j = i; j < directions.size(); ++j) {
        BC_pairs_.emplace_back(directions[i], directions[j]);
      }
    }

    if (world.rank() == 0) {
      print("Number of BC pairs: ", BC_pairs_.size());
      print("Number of frequencies: ", frequencies_.size());
      for (const auto &pair : BC_pairs_) {
        print("BC pair: ", pair.first, pair.second);
      }
    }
  }

  // Loopable interface
  size_t num_BC_pairs() const { return BC_pairs_.size(); }
  size_t num_frequencies() const { return frequencies_.size(); }

  // get state from state_map
  [[nodiscard]] ResponseState get_state(char dir) const {
    return state_map_.at(dir);
  }

  std::pair<char, char> get_BC_pair(size_t index) const {
    return BC_pairs_.at(index);
  }

  // Compute and save VBC given BC pair index and frequency indices
  auto compute_and_save(size_t bc_index, size_t b_index,
                        size_t c_index) -> ResponseVector {
    auto [B_dir, C_dir] = BC_pairs_.at(bc_index);

    const auto &omega_b = frequencies_.at(b_index);
    const auto &omega_c = frequencies_.at(c_index);

    auto state_b = state_map_[B_dir];
    auto state_c = state_map_[C_dir];
    state_b.set_frequency_index(b_index);
    state_c.set_frequency_index(c_index);

    auto dipB = std::get<DipolePerturbation>(state_b.perturbation);
    auto dipC = std::get<DipolePerturbation>(state_c.perturbation);

    if (world_.rank() == 0) {
      print("Computing VBC for directions: ", B_dir, C_dir);
      print("Frequencies: ", omega_b, omega_c);
    }

    SecondOrderResponseState vbc_state(dipB, dipC, omega_b, omega_c,
                                       FunctionDefaults<3>::get_thresh(),
                                       spin_restricted_);

    auto vbc_filename = vbc_state.response_filename();
    if (world_.rank() == 0) {
      print("VBC filename: ", vbc_filename);
    }
    if (fs::exists(vbc_filename + ".00000")) {
      if (world_.rank() == 0) {
        print("VBC file already exists, skipping computation.");
      }
      // load and return
      ResponseVector loaded_vbc;
      if (world_.rank() == 0) {
        print("loading ", 2 * num_orbitals_, " orbitals");
      }
      load_response_vector(world_, num_orbitals_, vbc_state, loaded_vbc);

      return loaded_vbc;
    }

    ResponseVector B_loaded, C_loaded;
    load_response_vector(world_, num_orbitals_, state_b, B_loaded, 0, b_index);
    load_response_vector(world_, num_orbitals_, state_c, C_loaded, 0, c_index);

    ResponseVector B_vec =
        make_response_vector(num_orbitals_, false, spin_restricted_);
    ResponseVector C_vec =
        make_response_vector(num_orbitals_, false, spin_restricted_);
    promote_response_vector(world_, B_loaded, B_vec);
    promote_response_vector(world_, C_loaded, C_vec);

    auto result_vec = compute_vbc_response(
        world_, gs_, B_vec, C_vec, perturbation_vectors_.at(B_dir),
        perturbation_vectors_.at(C_dir), dipole_operator(B_dir),
        dipole_operator(C_dir));

    save_response_vector(world_, vbc_state, result_vec);
    return result_vec;
  }

  // output BC pairs
  std::vector<std::pair<char, char>> get_BC_pairs() const { return BC_pairs_; }

private:
  World &world_;
  const GroundStateData &gs_;
  std::string directions_;
  std::vector<double> frequencies_;
  bool spin_restricted_;
  size_t num_orbitals_;

  std::map<char, ResponseState> state_map_;
  std::map<char, vector_real_function_3d> perturbation_vectors_;
  std::vector<std::pair<char, char>> BC_pairs_;

  real_function_3d dipole_operator(char dir) const {
    std::vector<int> f(3, 0);
    f[dir - 'x'] = 1;
    return real_factory_3d(world_).functor(
        real_functor_3d(new MomentFunctor(f)));
  }

  static auto K(const vector_real_function_3d &bra,
                const vector_real_function_3d &ket) {
    auto &world = bra[0].world();
    Exchange<double, 3> k{world, 1.e-10};
    k.set_bra_and_ket(bra, ket);
    return k;
  }

  static auto compute_g(World &world, const vector_real_function_3d &Aleft,
                        const vector_real_function_3d &Aright,
                        const vector_real_function_3d &Bleft,
                        const vector_real_function_3d &Bright,
                        const vector_real_function_3d &phix,
                        const vector_real_function_3d &phiy) {
    auto x_phi = mul(world, Aleft, Aright, true);
    auto y_phi = mul(world, Bleft, Bright, true);
    world.gop.fence();
    auto rho = sum(world, x_phi, true);
    world.gop.fence();
    rho += sum(world, y_phi, true);
    world.gop.fence();
    auto lo = 1.e-10;
    real_convolution_3d op =
        CoulombOperator(world, lo, FunctionDefaults<3>::get_thresh());
    auto temp_J = apply(op, rho);

    auto Jx = mul(world, temp_J, phix, false);
    auto Jy = mul(world, temp_J, phiy, false);

    auto ka = K(Aleft, Aright)(phix); // what happens to k after this?
    auto kb = K(Bleft, Bright)(phix);
    auto ka_conj = K(Aright, Aleft)(phiy);
    auto kb_conj = K(Bright, Bleft)(phiy);
    world.gop.fence();
    // ideally it runs and the Exchange operator is freed

    auto Kx = gaxpy_oop(1.0, ka, 1.0, kb, true);
    auto Ky = gaxpy_oop(1.0, ka_conj, 1.0, kb_conj, true);
    world.gop.fence();

    auto result_x = gaxpy_oop(2.0, Jx, -1.0, Kx, true);
    auto result_y = gaxpy_oop(2.0, Jy, -1.0, Ky, true);

    return std::make_pair(result_x, result_y);
  };
  // This function constructs the J and K operators with A and B and applies
  static auto compute_vbc_i(World &world, const vector_real_function_3d &bx,
                            const vector_real_function_3d &by,
                            const vector_real_function_3d &cx,
                            const vector_real_function_3d &cy,
                            const vector_real_function_3d &zeta_bc,
                            const vector_real_function_3d &phi0,
                            const real_function_3d &v) {
    // left right pairs ka and kb and apply apply
    madness::QProjector<double, 3> Q(phi0);
    auto [gzeta_x, gzeta_y] =
        compute_g(world, bx, cy, phi0, zeta_bc, phi0, phi0);
    gzeta_x = -1.0 * Q(gzeta_x);
    gzeta_y = -1.0 * Q(gzeta_y);

    auto [gbc_x, gbc_y] = compute_g(world, bx, phi0, phi0, by, cx, cy);
    auto [gbc_phi_x, gbc_phi_y] =
        compute_g(world, bx, phi0, phi0, by, phi0, phi0);

    auto vcx = mul(world, v, cx, true);
    auto vcy = mul(world, v, cy, true);

    auto fbx = -1.0 * Q(gbc_x + vcx);
    auto fby = -1.0 * Q(gbc_y + vcy);

    auto vb_phi = mul(world, v, phi0, true);

    auto fb_phi_x = gbc_phi_x + vb_phi;
    auto fb_phi_y = gbc_phi_y + vb_phi;

    auto m_fbx = matrix_inner(world, phi0, fb_phi_x);
    auto m_fby = matrix_inner(world, phi0, fb_phi_y);

    auto fphi_x = transform(world, cx, m_fbx, true);
    auto fphi_y = transform(world, cy, m_fby, true);

    auto result_x = truncate(gzeta_x + fbx + fphi_x,
                             FunctionDefaults<3>::get_thresh(), true);
    auto result_y = truncate(gzeta_y + fby + fphi_y,
                             FunctionDefaults<3>::get_thresh(), true);

    return std::make_pair(result_x, result_y);
  };
  ResponseVector compute_vbc_response(World &world, const GroundStateData &gs,
                                      const ResponseVector &B,
                                      const ResponseVector &C,
                                      const vector_real_function_3d &VB,
                                      const vector_real_function_3d &VC,
                                      const real_function_3d &VB_op,
                                      const real_function_3d &VC_op) {
    const auto &bx = std::get<DynamicRestrictedResponse>(B).x_alpha;
    const auto &by = std::get<DynamicRestrictedResponse>(B).y_alpha;
    const auto &cx = std::get<DynamicRestrictedResponse>(C).x_alpha;
    const auto &cy = std::get<DynamicRestrictedResponse>(C).y_alpha;

    auto zeta_bc = VBCComputer::make_zeta_bc(world, by, cx, gs.orbitals);
    auto zeta_cb = VBCComputer::make_zeta_bc(world, cy, bx, gs.orbitals);

    auto [bcx, bcy] =
        compute_vbc_i(world, bx, by, cx, cy, zeta_bc, gs.orbitals, VB_op);
    auto [cbx, cby] =
        compute_vbc_i(world, cx, cy, bx, by, zeta_cb, gs.orbitals, VC_op);

    auto result_vec = DynamicRestrictedResponse(num_orbitals_);

    result_vec.x_alpha = gaxpy_oop(1.0, bcx, 1.0, cbx, true);
    result_vec.y_alpha = gaxpy_oop(1.0, bcy, 1.0, cby, true);
    result_vec.flatten();

    return result_vec;
  }

  static vector_real_function_3d
  make_zeta_bc(World &world, const vector_real_function_3d &by,
               const vector_real_function_3d &cx,
               const vector_real_function_3d &phi0) {
    auto mat = matrix_inner(world, by, cx);
    return -1.0 * transform(world, phi0, mat, true);
  }
};

} // namespace madness
