
#ifndef MADNESS_APPS_MOLRESPONSE_RESPONSE_MACROTASK_HPP
#define MADNESS_APPS_MOLRESPONSE_RESPONSE_MACROTASK_HPP

// #include "functypedefs.h"
// #include "functypedefs.h"
#include "SCF.h"
#include "functypedefs.h"
#include "vmra.h"
#include "x_space.h"
#include <madness.h>
#include <madness/chem/SCFOperators.h>
#include <madness/chem/projector.h>
#include <madness/mra/macrotaskpartitioner.h>
#include <madness/mra/macrotaskq.h>
#include <madness/world/cloud.h>
#include <tuple>

namespace madness {

// Set the default algorithm for the exchange operator to small memory
// global variable

class VBC_task2 : public MacroTaskOperationBase {

  class Partitioner : public MacroTaskPartitioner {
  public:
    explicit Partitioner(long stride) {
      max_batch_size = 1;
      result_stride = stride;
      policy = "strided";
    }
  };

public:
  long stride;
  explicit VBC_task2(long stride) : stride(stride) {
    partitioner.reset(new Partitioner(stride));
  }
  // index b, index c, response B, response C, zetaBC, zetaCB, phi0,
  // perturbation
  typedef std::tuple<
      const std::vector<int> &, const std::vector<int> &,
      const std::vector<int> &, const vector_real_function_3d &,
      const vector_real_function_3d &, const vector_real_function_3d &,
      const vector_real_function_3d &, const vector_real_function_3d &,
      const vector_real_function_3d &>
      argtupleT;

  using resultT = vector_real_function_3d;

  resultT allocator(World &world, const argtupleT &args) const {

    auto num_states = static_cast<int>(std::get<1>(args).size());
    auto num_orbitals = static_cast<int>(std::get<7>(args).size());
    auto num_functions = 2 * num_states * num_orbitals;
    print("allocator: ", num_states, num_orbitals, num_functions);

    return zero_functions_compressed<double, 3>(world, num_functions);
  };

  resultT
  operator()(const std::vector<int> &dummy_i, const std::vector<int> &b,
             const std::vector<int> &c, const vector_real_function_3d &B,
             const vector_real_function_3d &C,
             const vector_real_function_3d &zeta_BC,
             const vector_real_function_3d &zeta_CB,
             const vector_real_function_3d &phi0,
             const vector_real_function_3d &dipole_perturbations) const {

    World &world = phi0[0].world();
    madness::QProjector<double, 3> Q(phi0);
    auto thresh = FunctionDefaults<3>::get_thresh();

    auto K = [&](const vecfuncT &ket, const vecfuncT &bra) {
      const double lo = 1.e-10;
      auto &world = ket[0].world();
      Exchange<double, 3> k{world, lo};
      k.set_bra_and_ket(bra, ket);

      std::string algorithm_ = "small_memory";
      if (algorithm_ == "multiworld") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient);
      } else if (algorithm_ == "multiworld_row") {
        k.set_algorithm(
            Exchange<double, 3>::Algorithm::multiworld_efficient_row);
      } else if (algorithm_ == "largemem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::large_memory);
      } else if (algorithm_ == "smallmem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::small_memory);
      }

      return k;
    };
    // A and B are the response pairs that make up response density
    // \gamma_{a} = |xa><phi| + |phi><ya|
    // This function constructs the J and K operators with A and B and applies
    // on x
    auto compute_g = [&](const vector_real_function_3d &Aleft,
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
    auto compute_vbc_i = [&](const vector_real_function_3d &bx,
                             const vector_real_function_3d &by,
                             const vector_real_function_3d &cx,
                             const vector_real_function_3d &cy,
                             const vector_real_function_3d &zeta_bc,
                             const vector_real_function_3d &phi,
                             const real_function_3d &v) {
      // left right pairs ka and kb and apply apply
      auto [gzeta_x, gzeta_y] = compute_g(bx, cy, phi0, zeta_bc, phi0, phi0);
      gzeta_x = -1.0 * Q(gzeta_x);
      gzeta_y = -1.0 * Q(gzeta_y);

      auto [gbc_x, gbc_y] = compute_g(bx, phi0, phi0, by, cx, cy);
      auto [gbc_phi_x, gbc_phi_y] = compute_g(bx, phi0, phi0, by, phi0, phi0);

      auto vcx = mul(world, v, cx, true);
      auto vcy = mul(world, v, cy, true);

      auto fbx = -1.0 * Q(gbc_x + vcx);
      auto fby = -1.0 * Q(gbc_y + vcy);

      auto vb_phi = mul(world, v, phi, true);

      auto fb_phi_x = gbc_phi_x + vb_phi;
      auto fb_phi_y = gbc_phi_y + vb_phi;

      auto m_fbx = matrix_inner(world, phi, fb_phi_x);
      auto m_fby = matrix_inner(world, phi, fb_phi_y);

      auto fphi_x = transform(world, cx, m_fbx, true);
      auto fphi_y = transform(world, cy, m_fby, true);

      auto result_x = truncate(gzeta_x + fbx + fphi_x, thresh, true);
      auto result_y = truncate(gzeta_y + fby + fphi_y, thresh, true);

      return std::make_pair(result_x, result_y);
    };

    int num_orbitals = static_cast<int>(phi0.size());
    int num_states = static_cast<int>(b.size());

    auto bc_indexer = x_space_indexer(num_orbitals);
    auto zeta_indexer = response_space_index(num_orbitals);
    auto compute_result = [&](const int &i) {
      auto results =
          zero_functions_compressed<double, 3>(world, 2 * num_orbitals);

      auto bi = b[i];
      auto ci = c[i];
      const auto &bx = bc_indexer.get_x_state(bi, B);
      const auto &by = bc_indexer.get_y_state(bi, B);
      const auto &cx = bc_indexer.get_x_state(ci, C);
      const auto &cy = bc_indexer.get_y_state(ci, C);

      const auto &zeta_bc = zeta_indexer.get_x_state(i, zeta_BC);
      const auto &zeta_cb = zeta_indexer.get_x_state(i, zeta_CB);

      const auto &vb = dipole_perturbations[bi];
      const auto &vc = dipole_perturbations[ci];

      auto [bcx, bcy] = compute_vbc_i(bx, by, cx, cy, zeta_bc, phi0, vb);
      auto [cbx, cby] = compute_vbc_i(cx, cy, bx, by, zeta_cb, phi0, vc);

      auto result_x = gaxpy_oop(1.0, bcx, 1.0, cbx, true);
      auto result_y = gaxpy_oop(1.0, bcy, 1.0, cby, true);

      for (int j = 0; j < num_orbitals; j++) {

        auto index_x = j;
        auto norm_x = result_x[j].norm2();
        print("i,j = (", i, j, ") index_x: ", index_x, " norm_x: ", norm_x);
        results[index_x] = result_x[j];
      }

      for (int j = 0; j < num_orbitals; j++) {

        auto index_y = j + num_orbitals;
        auto norm_y = result_y[j].norm2();
        print("i,j = (", i, j, ") index_y: ", index_y, " norm_y: ", norm_y);
        results[index_y] = result_y[j];
      }

      return results;
    };
    const int ij = static_cast<int>(batch.result.begin);
    return compute_result(ij);
  }
};

class VBC_task_one : public MacroTaskOperationBase {

  class Partitioner : public MacroTaskPartitioner {
  public:
    explicit Partitioner(long stride) {
      max_batch_size = 1;
      result_stride = stride;
      policy = "strided";
    }
  };

public:
  long stride;
  explicit VBC_task_one(long stride) : stride(stride) {
    partitioner.reset(new Partitioner(stride));
  }
  // index b, index c, response B, response C, phi0,
  typedef std::tuple<const std::vector<int> &, const std::vector<int> &,
                     const std::vector<int> &, const vector_real_function_3d &,
                     const vector_real_function_3d &,
                     const vector_real_function_3d &>
      argtupleT;

  using resultT = vector_real_function_3d;

  resultT allocator(World &world, const argtupleT &args) const {
    auto num_states = static_cast<int>(std::get<1>(args).size());
    auto num_orbitals = static_cast<int>(std::get<5>(args).size());
    auto num_functions = 2 * num_states * num_orbitals;
    print("allocator: ", num_states, num_orbitals, num_functions);

    return zero_functions_compressed<double, 3>(world, num_functions);
  };

  resultT operator()(const std::vector<int> &dummy_i, const std::vector<int> &b,
                     const std::vector<int> &c,
                     const vector_real_function_3d &B,
                     const vector_real_function_3d &C,
                     const vector_real_function_3d &phi0) const {

    World &world = phi0[0].world();
    madness::QProjector<double, 3> Q(phi0);
    auto thresh = FunctionDefaults<3>::get_thresh();

    auto K = [&](const vecfuncT &ket, const vecfuncT &bra) {
      const double lo = 1.e-10;
      auto &world = ket[0].world();
      Exchange<double, 3> k{world, lo};
      k.set_bra_and_ket(bra, ket);

      std::string algorithm_ = "small_memory";
      if (algorithm_ == "multiworld") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient);
      } else if (algorithm_ == "multiworld_row") {
        k.set_algorithm(
            Exchange<double, 3>::Algorithm::multiworld_efficient_row);
      } else if (algorithm_ == "largemem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::large_memory);
      } else if (algorithm_ == "smallmem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::small_memory);
      }

      return k;
    };
    // A and B are the response pairs that make up response density
    // \gamma_{a} = |xa><phi| + |phi><ya|
    // This function constructs the J and K operators with A and B and applies
    // on x
    auto compute_g = [&](const vector_real_function_3d &Aleft,
                         const vector_real_function_3d &Aright,
                         const vector_real_function_3d &phi) {
      auto x_phi = mul(world, Aleft, Aright, true);
      auto rho = sum(world, x_phi, true);
      auto lo = 1.e-10;
      real_convolution_3d op =
          CoulombOperator(world, lo, FunctionDefaults<3>::get_thresh());
      auto temp_J = apply(op, rho);

      auto J = mul(world, temp_J, phi, true);
      auto Kx = K(Aleft, Aright)(phi); // what happens to k after this?
      auto Ky = K(Aright, Aleft)(phi);
      world.gop.fence();
      auto result_x = gaxpy_oop(2.0, J, -1.0, Kx, true);
      auto result_y = gaxpy_oop(2.0, J, -1.0, Ky, true);

      return std::make_pair(result_x, result_y);
    };
    auto compute_vbc_i = [&](const vector_real_function_3d &bx,
                             const vector_real_function_3d &cy,
                             const vector_real_function_3d &phi) {
      // left right pairs ka and kb and apply apply
      auto [gzeta_x, gzeta_y] = compute_g(bx, cy, phi0);
      gzeta_x = -1.0 * Q(gzeta_x);
      gzeta_y = -1.0 * Q(gzeta_y);
      return std::make_pair(gzeta_x, gzeta_y);
    };

    int num_orbitals = static_cast<int>(phi0.size());
    int num_states = static_cast<int>(b.size());

    auto bc_indexer = x_space_indexer(num_orbitals);
    auto zeta_indexer = response_space_index(num_orbitals);
    auto compute_result = [&](const int &i) {
      auto results =
          zero_functions_compressed<double, 3>(world, 2 * num_orbitals);

      auto bi = b[i];
      auto ci = c[i];
      const auto &bx = bc_indexer.get_x_state(bi, B);
      const auto &by = bc_indexer.get_y_state(bi, B);

      const auto &cx = bc_indexer.get_x_state(ci, C);
      const auto &cy = bc_indexer.get_y_state(ci, C);

      auto [bcx, bcy] = compute_vbc_i(bx, cy, phi0);
      auto [cbx, cby] = compute_vbc_i(cx, by, phi0);

      auto result_x = gaxpy_oop(1.0, bcx, 1.0, cbx, true);
      auto result_y = gaxpy_oop(1.0, bcy, 1.0, cby, true);

      for (int j = 0; j < num_orbitals; j++) {

        auto index_x = j;
        auto norm_x = result_x[j].norm2();
        print("i,j = (", i, j, ") index_x: ", index_x, " norm_x: ", norm_x);
        results[index_x] = result_x[j];
      }

      for (int j = 0; j < num_orbitals; j++) {

        auto index_y = j + num_orbitals;
        auto norm_y = result_y[j].norm2();
        print("i,j = (", i, j, ") index_y: ", index_y, " norm_y: ", norm_y);
        results[index_y] = result_y[j];
      }

      return results;
    };
    const int ij = static_cast<int>(batch.result.begin);
    return compute_result(ij);
  }
};

class VBC_task_two : public MacroTaskOperationBase {

  class Partitioner : public MacroTaskPartitioner {
  public:
    explicit Partitioner(long stride) {
      max_batch_size = 1;
      result_stride = stride;
      policy = "strided";
    }
  };

public:
  long stride;
  explicit VBC_task_two(long stride) : stride(stride) {
    partitioner.reset(new Partitioner(stride));
  }
  typedef std::tuple<const std::vector<int> &, const vector_real_function_3d &,
                     const vector_real_function_3d &>
      argtupleT;

  using resultT = vector_real_function_3d;

  resultT allocator(World &world, const argtupleT &args) const {

    auto num_states = static_cast<int>(std::get<1>(args).size());
    auto num_orbitals = static_cast<int>(std::get<2>(args).size());
    auto num_functions = 2 * num_states * num_orbitals;
    print("allocator: ", num_states, num_orbitals, num_functions);

    return zero_functions_compressed<double, 3>(world, num_functions);
  };

  resultT operator()(const std::vector<int> &dummy_i,
                     const vector_real_function_3d &zeta_BC,
                     const vector_real_function_3d &phi0) const {

    World &world = phi0[0].world();
    madness::QProjector<double, 3> Q(phi0);
    auto thresh = FunctionDefaults<3>::get_thresh();

    auto K = [&](const vecfuncT &ket, const vecfuncT &bra) {
      const double lo = 1.e-10;
      auto &world = ket[0].world();
      Exchange<double, 3> k{world, lo};
      k.set_bra_and_ket(bra, ket);

      std::string algorithm_ = "small_memory";
      if (algorithm_ == "multiworld") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient);
      } else if (algorithm_ == "multiworld_row") {
        k.set_algorithm(
            Exchange<double, 3>::Algorithm::multiworld_efficient_row);
      } else if (algorithm_ == "largemem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::large_memory);
      } else if (algorithm_ == "smallmem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::small_memory);
      }

      return k;
    };
    // A and B are the response pairs that make up response density
    // \gamma_{a} = |xa><phi| + |phi><ya|
    // This function constructs the J and K operators with A and B and applies
    // on x
    auto compute_g = [&](const vector_real_function_3d &Aleft,
                         const vector_real_function_3d &Aright,
                         const vector_real_function_3d &phi) {
      auto x_phi = mul(world, Aleft, Aright, true);
      auto rho = sum(world, x_phi, true);
      auto lo = 1.e-10;
      real_convolution_3d op =
          CoulombOperator(world, lo, FunctionDefaults<3>::get_thresh());
      auto temp_J = apply(op, rho);

      auto J = mul(world, temp_J, phi, true);
      auto Kx = K(Aleft, Aright)(phi); // what happens to k after this?
      auto Ky = K(Aright, Aleft)(phi);
      world.gop.fence();
      auto result_x = gaxpy_oop(2.0, J, -1.0, Kx, true);
      auto result_y = gaxpy_oop(2.0, J, -1.0, Ky, true);

      return std::make_pair(result_x, result_y);
    };

    auto compute_vbc_i = [&](const vector_real_function_3d &zeta_bc,
                             const vector_real_function_3d &phi) {
      // left right pairs ka and kb and apply apply
      auto [gzeta_x, gzeta_y] = compute_g(zeta_bc, phi, phi);
      gzeta_x = -1.0 * Q(gzeta_x);
      gzeta_y = -1.0 * Q(gzeta_y);
      return std::make_pair(gzeta_x, gzeta_y);
    };

    int num_orbitals = static_cast<int>(phi0.size());
    auto zeta_indexer = response_space_index(num_orbitals);

    auto compute_result = [&](const int &i) {
      auto results =
          zero_functions_compressed<double, 3>(world, 2 * num_orbitals);

      const auto &zeta_bc = zeta_indexer.get_x_state(i, zeta_BC);
      auto [result_x, result_y] = compute_vbc_i(zeta_BC, phi0);
      for (int j = 0; j < num_orbitals; j++) {
        auto index_x = j;
        auto norm_x = result_x[j].norm2();
        print("i,j = (", i, j, ") index_x: ", index_x, " norm_x: ", norm_x);
        results[index_x] = result_x[j];
      }
      for (int j = 0; j < num_orbitals; j++) {

        auto index_y = j + num_orbitals;
        auto norm_y = result_y[j].norm2();
        print("i,j = (", i, j, ") index_y: ", index_y, " norm_y: ", norm_y);
        results[index_y] = result_y[j];
      }

      return results;
    };
    const int ij = static_cast<int>(batch.result.begin);
    return compute_result(ij);
  }
};

class VBC_task_three : public MacroTaskOperationBase {

  class Partitioner : public MacroTaskPartitioner {
  public:
    explicit Partitioner(long stride) {
      max_batch_size = 1;
      result_stride = stride;
      policy = "strided";
    }
  };

public:
  long stride;
  explicit VBC_task_three(long stride) : stride(stride) {
    partitioner.reset(new Partitioner(stride));
  }
  // index b, index c, response B, response C, zetaBC, zetaCB, phi0,
  // perturbation
  typedef std::tuple<const std::vector<int> &, const std::vector<int> &,
                     const std::vector<int> &, const vector_real_function_3d &,
                     const vector_real_function_3d &,
                     const vector_real_function_3d &,
                     const vector_real_function_3d &>
      argtupleT;

  using resultT = vector_real_function_3d;

  resultT allocator(World &world, const argtupleT &args) const {

    auto num_states = static_cast<int>(std::get<1>(args).size());
    auto num_orbitals = static_cast<int>(std::get<5>(args).size());
    auto num_functions = 2 * num_states * num_orbitals;
    print("allocator: ", num_states, num_orbitals, num_functions);

    return zero_functions_compressed<double, 3>(world, num_functions);
  };

  resultT
  operator()(const std::vector<int> &dummy_i, const std::vector<int> &b,
             const std::vector<int> &c, const vector_real_function_3d &B,
             const vector_real_function_3d &C,
             const vector_real_function_3d &phi0,
             const vector_real_function_3d &dipole_perturbations) const {

    World &world = phi0[0].world();
    madness::QProjector<double, 3> Q(phi0);
    auto thresh = FunctionDefaults<3>::get_thresh();

    auto K = [&](const vecfuncT &ket, const vecfuncT &bra) {
      const double lo = 1.e-10;
      auto &world = ket[0].world();
      Exchange<double, 3> k{world, lo};
      k.set_bra_and_ket(bra, ket);

      std::string algorithm_ = "small_memory";
      if (algorithm_ == "multiworld") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient);
      } else if (algorithm_ == "multiworld_row") {
        k.set_algorithm(
            Exchange<double, 3>::Algorithm::multiworld_efficient_row);
      } else if (algorithm_ == "largemem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::large_memory);
      } else if (algorithm_ == "smallmem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::small_memory);
      }

      return k;
    };
    // A and B are the response pairs that make up response density
    // \gamma_{a} = |xa><phi| + |phi><ya|
    // This function constructs the J and K operators with A and B and applies
    // on x
    auto compute_g = [&](const vector_real_function_3d &Bx,
                         const vector_real_function_3d &By,
                         const vector_real_function_3d &Cx,
                         const vector_real_function_3d &Cy,
                         const vector_real_function_3d &phi) {
      auto x_phi = mul(world, Bx, phi, true);
      auto y_phi = mul(world, phi, By, true);
      world.gop.fence();
      auto rho = sum(world, x_phi, true);
      world.gop.fence();
      rho += sum(world, y_phi, true);
      world.gop.fence();
      auto lo = 1.e-10;
      real_convolution_3d op =
          CoulombOperator(world, lo, FunctionDefaults<3>::get_thresh());
      auto temp_JB = apply(op, rho);

      // Two types here, apply on Cx and apply on phi

      auto JBCx = mul(world, temp_JB, Cx, false);
      auto JBCy = mul(world, temp_JB, Cy, false);
      auto JBphi = mul(world, temp_JB, phi, false);

      auto Kapply_x_Bx = K(Bx, phi);
      auto Kapply_x_By = K(phi, By);
      auto Kapply_y_By = K(By, phi);
      auto Kapply_y_Bx = K(phi, Bx);

      auto Kbxcx = Kapply_x_Bx(Cx);
      auto Kbycx = Kapply_x_By(Cx);
      auto Kbycy = Kapply_y_By(Cy);
      auto Kbxcy = Kapply_y_Bx(Cy);

      auto Kbxphi = Kapply_x_Bx(phi);
      auto Kbyphi = Kapply_x_By(phi);
      auto Kbxphi_conj = Kapply_y_Bx(phi);
      auto Kbyphi_conj = Kapply_y_By(phi);

      auto kbcx = gaxpy_oop(1.0, Kbxcx, 1.0, Kbycx, false);
      auto bxcy = gaxpy_oop(1.0, Kbxcy, 1.0, Kbycy, false);
      auto Kbphi = gaxpy_oop(1.0, Kbxphi, 1.0, Kbyphi, false);
      auto Kbphi_conj = gaxpy_oop(1.0, Kbxphi_conj, 1.0, Kbyphi_conj, false);

      world.gop.fence();
      auto result_x = gaxpy_oop(2.0, JBCx, -1.0, kbcx, false);
      auto result_y = gaxpy_oop(2.0, JBCy, -1.0, bxcy, false);
      auto result_phi_x = gaxpy_oop(2.0, JBphi, -1.0, Kbphi, false);
      auto result_phi_y = gaxpy_oop(2.0, JBphi, -1.0, Kbphi_conj, false);
      world.gop.fence();

      return std::make_tuple(result_x, result_y, result_phi_x, result_phi_y);
    };
    auto compute_vbc_i = [&](const vector_real_function_3d &bx,
                             const vector_real_function_3d &by,
                             const vector_real_function_3d &cx,
                             const vector_real_function_3d &cy,
                             const vector_real_function_3d &phi,
                             const real_function_3d &v) {
      // left right pairs ka and kb and apply apply

      auto [gbc_x, gbc_y, gbc_phi_x, gbc_phi_y] =
          compute_g(bx, by, cx, cy, phi);

      auto vcx = mul(world, v, cx, true);
      auto vcy = mul(world, v, cy, true);

      auto fbx = -1.0 * Q(gbc_x + vcx);
      auto fby = -1.0 * Q(gbc_y + vcy);

      auto vb_phi = mul(world, v, phi, true);

      auto fb_phi_x = gbc_phi_x + vb_phi;
      auto fb_phi_y = gbc_phi_y + vb_phi;

      auto m_fbx = matrix_inner(world, phi, fb_phi_x);
      auto m_fby = matrix_inner(world, phi, fb_phi_y);

      auto fphi_x = transform(world, cx, m_fbx, true);
      auto fphi_y = transform(world, cy, m_fby, true);

      auto result_x = fbx + fphi_x;
      auto result_y = fby + fphi_y;

      return std::make_pair(result_x, result_y);
    };

    int num_orbitals = static_cast<int>(phi0.size());
    int num_states = static_cast<int>(b.size());

    auto bc_indexer = x_space_indexer(num_orbitals);
    auto compute_result = [&](const int &i) {
      auto results =
          zero_functions_compressed<double, 3>(world, 2 * num_orbitals);

      auto bi = b[i];
      auto ci = c[i];
      const auto &bx = bc_indexer.get_x_state(bi, B);
      const auto &by = bc_indexer.get_y_state(bi, B);
      const auto &cx = bc_indexer.get_x_state(ci, C);
      const auto &cy = bc_indexer.get_y_state(ci, C);

      const auto &vb = dipole_perturbations[bi];
      const auto &vc = dipole_perturbations[ci];

      auto [bcx, bcy] = compute_vbc_i(bx, by, cx, cy, phi0, vb);
      auto [cbx, cby] = compute_vbc_i(cx, cy, bx, by, phi0, vc);

      auto result_x = gaxpy_oop(1.0, bcx, 1.0, cbx, true);
      auto result_y = gaxpy_oop(1.0, bcy, 1.0, cby, true);

      for (int j = 0; j < num_orbitals; j++) {

        auto index_x = j;
        auto norm_x = result_x[j].norm2();
        print("i,j = (", i, j, ") index_x: ", index_x, " norm_x: ", norm_x);
        results[index_x] = result_x[j];
      }

      for (int j = 0; j < num_orbitals; j++) {

        auto index_y = j + num_orbitals;
        auto norm_y = result_y[j].norm2();
        print("i,j = (", i, j, ") index_y: ", index_y, " norm_y: ", norm_y);
        results[index_y] = result_y[j];
      }

      return results;
    };
    const int ij = static_cast<int>(batch.result.begin);
    return compute_result(ij);
  }
};

class VBC_task_one_i : public MacroTaskOperationBase {

  class Partitioner : public MacroTaskPartitioner {
  public:
    Partitioner() {
      max_batch_size = 1;
      //  result_stride = stride;
      //  policy = "strided";
    }
  };

public:
  VBC_task_one_i() { partitioner.reset(new Partitioner()); }
  // index b, index c, response B, response C, zetaBC, zetaCB, phi0,
  // perturbation
  typedef std::tuple<const std::vector<int> &, const std::vector<int> &,
                     const vector_real_function_3d &,
                     const vector_real_function_3d &,
                     const vector_real_function_3d &>
      argtupleT;

  using resultT = vector_real_function_3d;

  resultT allocator(World &world, const argtupleT &args) const {
    auto num_functions = static_cast<int>(std::get<0>(args).size());
    return zero_functions_compressed<double, 3>(world, num_functions);
  };

  resultT operator()(const std::vector<int> &b, const std::vector<int> &c,
                     const vector_real_function_3d &B,
                     const vector_real_function_3d &C,
                     const vector_real_function_3d &phi0) const {

    World &world = phi0[0].world();
    madness::QProjector<double, 3> Q(phi0);
    auto thresh = FunctionDefaults<3>::get_thresh();

    auto K = [&](const vecfuncT &ket, const vecfuncT &bra) {
      const double lo = 1.e-10;
      auto &world = ket[0].world();
      Exchange<double, 3> k{world, lo};
      k.set_bra_and_ket(bra, ket);
      std::string algorithm_ = "small_memory";
      if (algorithm_ == "multiworld") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient);
      } else if (algorithm_ == "multiworld_row") {
        k.set_algorithm(
            Exchange<double, 3>::Algorithm::multiworld_efficient_row);
      } else if (algorithm_ == "largemem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::large_memory);
      } else if (algorithm_ == "smallmem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::small_memory);
      }

      return k;
    };
    // In this example we only apply once to either x or y orbital.  If x we
    // just reverse left and right in the call
    auto compute_g_i = [&](const vector_real_function_3d &Aleft,
                           const vector_real_function_3d &Aright,
                           const real_function_3d &phi) {
      auto x_phi = mul(world, Aleft, Aright, true);
      world.gop.fence();
      auto rho = sum(world, x_phi, true);
      auto lo = 1.e-10;
      real_convolution_3d op =
          CoulombOperator(world, lo, FunctionDefaults<3>::get_thresh());
      auto temp_J = apply(op, rho);
      auto Jx = temp_J * phi;
      auto kx = K(Aleft, Aright)(phi);
      auto result_x = gaxpy_oop(2.0, Jx.compress(), -1.0, kx.compress(), true);

      return result_x;
    };
    auto compute_vbc_i = [&](const vector_real_function_3d &bx,
                             const vector_real_function_3d &cy,
                             const vector_real_function_3d &phi, int j) {
      auto gzeta_x = compute_g_i(bx, cy, phi0[j]);
      gzeta_x = -1.0 * Q(gzeta_x);
      return gzeta_x;
    };
    int num_orbitals = static_cast<int>(phi0.size());
    auto bc_indexer = x_space_indexer(num_orbitals);
    auto compute_result = [&](const int &i) {
      resultT result = zero_functions_compressed<double, 3>(world, 1);
      auto bi = b[i];
      auto ci = c[i];
      auto response_orbital_index = i % (num_orbitals * 2);
      bool type_x = true;
      if (response_orbital_index < num_orbitals) {
        type_x = true;
      } else {
        type_x = false;
        response_orbital_index -= num_orbitals;
      }
      const auto &bx = bc_indexer.get_x_state(bi, B);
      const auto &by = bc_indexer.get_y_state(bi, B);
      const auto &cx = bc_indexer.get_x_state(ci, C);
      const auto &cy = bc_indexer.get_y_state(ci, C);

      real_function_3d one, two;
      if (type_x) {
        one = compute_vbc_i(bx, cy, phi0, response_orbital_index);
        two = compute_vbc_i(cx, by, phi0, response_orbital_index);
      } else {
        one = compute_vbc_i(cy, bx, phi0, response_orbital_index);
        two = compute_vbc_i(by, cx, phi0, response_orbital_index);
      }
      result[0] = gaxpy_oop(1.0, one, 1.0, two, true);
      return result;
    };
    const long ij = batch.result.begin;
    return compute_result(ij);
  }
};
class VBC_task_two_i_copy : public MacroTaskOperationBase {

  class Partitioner : public MacroTaskPartitioner {
  public:
    Partitioner() {
      max_batch_size = 1;
      //  result_stride = stride;
      //  policy = "strided";
    }
  };

public:
  VBC_task_two_i_copy() { partitioner.reset(new Partitioner()); }
  // index b, index c, response B, response C, zetaBC, zetaCB, phi0,
  // perturbation
  typedef std::tuple<const std::vector<int> &, const vector_real_function_3d &,
                     const vector_real_function_3d &>
      argtupleT;

  using resultT = vector_real_function_3d;

  resultT allocator(World &world, const argtupleT &args) const {

    // In this example the first index is index of all orbitals we are working
    // with
    auto num_functions = static_cast<int>(std::get<0>(args).size());
    // auto num_orbitals = static_cast<int>(std::get<7>(args).size());
    // auto num_functions = 2 * num_states * num_orbitals;
    // print("allocator: ", num_states, num_orbitals, num_functions);

    return zero_functions_compressed<double, 3>(world, num_functions);
  };

  resultT operator()(const std::vector<int> &state_index,
                     const vector_real_function_3d &zeta_BC,
                     const vector_real_function_3d &phi0) const {

    World &world = phi0[0].world();
    madness::QProjector<double, 3> Q(phi0);
    auto thresh = FunctionDefaults<3>::get_thresh();

    auto K = [&](const vecfuncT &ket, const vecfuncT &bra) {
      const double lo = 1.e-10;
      auto &world = ket[0].world();
      Exchange<double, 3> k{world, lo};
      k.set_bra_and_ket(bra, ket);

      std::string algorithm_ = "small_memory";
      if (algorithm_ == "multiworld") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient);
      } else if (algorithm_ == "multiworld_row") {
        k.set_algorithm(
            Exchange<double, 3>::Algorithm::multiworld_efficient_row);
      } else if (algorithm_ == "largemem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::large_memory);
      } else if (algorithm_ == "smallmem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::small_memory);
      }

      return k;
    };
    // In this example we only apply once to either x or y orbital.  If x we
    // just reverse left and right in the call
    auto compute_g_i = [&](const vector_real_function_3d &Aleft,
                           const vector_real_function_3d &Aright,
                           const real_function_3d &phi) {
      auto x_phi = mul(world, Aleft, Aright, true);
      world.gop.fence();
      auto rho = sum(world, x_phi, true);
      auto lo = 1.e-10;
      real_convolution_3d op =
          CoulombOperator(world, lo, FunctionDefaults<3>::get_thresh());
      auto temp_J = apply(op, rho);
      auto Jx = temp_J * phi;
      auto kx = K(Aleft, Aright)(phi);
      auto result_x = gaxpy_oop(2.0, Jx.compress(), -1.0, kx.compress(), true);

      return result_x;
    };

    auto compute_vbc_i = [&](const vector_real_function_3d &phi_bc,
                             const vector_real_function_3d &phi,
                             const real_function_3d phi_apply) {
      // left right pairs ka and kb and apply apply
      // TODO: There maybe an optimization that can tell whether we can skip
      // this orbital or not
      auto gzeta_x = compute_g_i(phi_bc, phi, phi_apply);
      gzeta_x = -1.0 * Q(gzeta_x);
      return gzeta_x;
    };
    int num_orbitals = static_cast<int>(phi0.size());
    auto zeta_indexer = response_space_index(num_orbitals);
    auto compute_result = [&](const int &i) {
      resultT result = zero_functions_compressed<double, 3>(world, 1);
      auto si = state_index[i];
      auto response_orbital_index = i % (num_orbitals * 2);
      bool type_x = true;
      if (response_orbital_index < num_orbitals) {
        type_x = true;
      } else {
        type_x = false;
        response_orbital_index -= num_orbitals;
      }
      const auto &zeta_bc = zeta_indexer.get_x_state(si, zeta_BC);
      real_function_3d one, two;
      if (type_x) {
        one = compute_vbc_i(zeta_bc, phi0, phi0[response_orbital_index]);
      } else {
        one = compute_vbc_i(phi0, zeta_bc, phi0[response_orbital_index]);
      }
      result[0] = one;

      return result;
    };
    const long ij = batch.result.begin;
    return compute_result(ij);
  }
};

class VBC_task_two_i : public MacroTaskOperationBase {

  class Partitioner : public MacroTaskPartitioner {
  public:
    Partitioner() {
      max_batch_size = 1;
      //  result_stride = stride;
      //  policy = "strided";
    }
  };

public:
  VBC_task_two_i() { partitioner.reset(new Partitioner()); }
  // index b, index c, response B, response C, zetaBC, zetaCB, phi0,
  // perturbation

  typedef std::tuple<const std::vector<int> &, const std::vector<int> &,
                     const vector_real_function_3d &,
                     const vector_real_function_3d &,
                     const vector_real_function_3d &>
      argtupleT;

  using resultT = vector_real_function_3d;

  resultT allocator(World &world, const argtupleT &args) const {
    auto num_functions = static_cast<int>(std::get<0>(args).size());
    return zero_functions_compressed<double, 3>(world, num_functions);
  };

  resultT operator()(const std::vector<int> &b, const std::vector<int> &c,
                     const vector_real_function_3d &B,
                     const vector_real_function_3d &C,
                     const vector_real_function_3d &phi0) const {

    World &world = phi0[0].world();
    madness::QProjector<double, 3> Q(phi0);
    auto thresh = FunctionDefaults<3>::get_thresh();

    auto K = [&](const vecfuncT &ket, const vecfuncT &bra) {
      const double lo = 1.e-10;
      auto &world = ket[0].world();
      Exchange<double, 3> k{world, lo};
      k.set_bra_and_ket(bra, ket);
      std::string algorithm_ = "small_memory";
      if (algorithm_ == "multiworld") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient);
      } else if (algorithm_ == "multiworld_row") {
        k.set_algorithm(
            Exchange<double, 3>::Algorithm::multiworld_efficient_row);
      } else if (algorithm_ == "largemem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::large_memory);
      } else if (algorithm_ == "smallmem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::small_memory);
      }

      return k;
    };
    // In this example we only apply once to either x or y orbital.  If x we
    // just reverse left and right in the call
    auto compute_g_i = [&](const vector_real_function_3d &Aleft,
                           const vector_real_function_3d &Aright,
                           const real_function_3d &phi) {
      auto x_phi = mul(world, Aleft, Aright, true);
      world.gop.fence();
      auto rho = sum(world, x_phi, true);
      auto lo = 1.e-10;
      real_convolution_3d op =
          CoulombOperator(world, lo, FunctionDefaults<3>::get_thresh());
      auto temp_J = apply(op, rho);
      auto Jx = temp_J * phi;
      auto kx = K(Aleft, Aright)(phi);
      auto result_x = gaxpy_oop(2.0, Jx.compress(), -1.0, kx.compress(), true);

      return result_x;
    };

    int num_orbitals = static_cast<int>(phi0.size());
    auto bc_indexer = x_space_indexer(num_orbitals);
    auto compute_result = [&](const int &i) {
      resultT result = zero_functions_compressed<double, 3>(world, 1);
      auto bi = b[i];
      auto ci = c[i];
      auto response_orbital_index = i % (num_orbitals * 2);
      bool type_x = true;
      if (response_orbital_index < num_orbitals) {
        type_x = true;
      } else {
        type_x = false;
        response_orbital_index -= num_orbitals;
      }
      const auto &by = bc_indexer.get_y_state(bi, B);
      const auto &cx = bc_indexer.get_x_state(ci, C);

      auto make_zeta_bc = [&](const vector_real_function_3d &by,
                              const vector_real_function_3d &cx,
                              const vector_real_function_3d &ground_orbitals) {
        auto matrix_bc = matrix_inner(world, by, cx);
        return -1.0 * transform(world, ground_orbitals, matrix_bc, true);
      };

      //  auto ztask = ComputeZetaBC(static_cast<int>(ground_orbitals.size()));

      vector_real_function_3d zeta_BC;
      zeta_BC = make_zeta_bc(by, cx, phi0);
      real_function_3d gzeta_x;

      if (type_x) {
        gzeta_x = compute_g_i(phi0, zeta_BC, phi0[response_orbital_index]);
      } else {
        gzeta_x = compute_g_i(zeta_BC, phi0, phi0[response_orbital_index]);
      }
      gzeta_x = -1.0 * Q(gzeta_x);
      result[0] = gzeta_x;
      return result;
    };
    const long ij = batch.result.begin;
    return compute_result(ij);
  }
};

class VBC_task_three_i : public MacroTaskOperationBase {

  class Partitioner : public MacroTaskPartitioner {
  public:
    Partitioner() {
      max_batch_size = 1;
      //  result_stride = stride;
      //  policy = "strided";
    }
  };

public:
  VBC_task_three_i() { partitioner.reset(new Partitioner()); }
  // index b, index c, response B, response C, zetaBC, zetaCB, phi0,
  // perturbation
  typedef std::tuple<
      const std::vector<int> &, const std::vector<int> &,
      const vector_real_function_3d &, const vector_real_function_3d &,
      const vector_real_function_3d &, const vector_real_function_3d &>
      argtupleT;

  using resultT = vector_real_function_3d;

  resultT allocator(World &world, const argtupleT &args) const {

    // In this example the first index is index of all orbitals we are working
    // with
    auto num_functions = static_cast<int>(std::get<0>(args).size());
    // auto num_orbitals = static_cast<int>(std::get<7>(args).size());
    // auto num_functions = 2 * num_states * num_orbitals;
    // print("allocator: ", num_states, num_orbitals, num_functions);

    return zero_functions_compressed<double, 3>(world, num_functions);
  };

  resultT
  operator()(const std::vector<int> &b, const std::vector<int> &c,
             const vector_real_function_3d &B, const vector_real_function_3d &C,
             const vector_real_function_3d &phi0,
             const vector_real_function_3d &dipole_perturbations) const {

    World &world = phi0[0].world();
    madness::QProjector<double, 3> Q(phi0);
    auto thresh = FunctionDefaults<3>::get_thresh();

    auto K = [&](const vecfuncT &ket, const vecfuncT &bra) {
      const double lo = 1.e-10;
      auto &world = ket[0].world();
      Exchange<double, 3> k{world, lo};
      k.set_bra_and_ket(bra, ket);

      std::string algorithm_ = "small_memory";
      if (algorithm_ == "multiworld") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient);
      } else if (algorithm_ == "multiworld_row") {
        k.set_algorithm(
            Exchange<double, 3>::Algorithm::multiworld_efficient_row);
      } else if (algorithm_ == "largemem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::large_memory);
      } else if (algorithm_ == "smallmem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::small_memory);
      }

      return k;
    };
    // In this example we only apply once to either x or y orbital.  If x we
    // just reverse left and right in the call
    auto compute_g_i = [&](const vector_real_function_3d &Bleft,
                           const vector_real_function_3d &Bright,
                           const vector_real_function_3d &phi,
                           const real_function_3d &cx_ii,
                           const real_function_3d &phi_i) {
      auto x_phi = mul(world, Bleft, phi, true);
      auto y_phi = mul(world, phi, Bright, true);
      world.gop.fence();
      auto rho = sum(world, x_phi, true);
      rho += sum(world, y_phi, true);
      auto lo = 1.e-10;
      real_convolution_3d op =
          CoulombOperator(world, lo, FunctionDefaults<3>::get_thresh());
      auto temp_J = apply(op, rho);

      auto Jxc = temp_J * cx_ii;
      auto Jx_phi = temp_J * phi_i;
      // mul(temp_J, phix, true);

      auto kone = K(Bleft, phi);
      auto ktwo = K(phi, Bright);

      auto kone_c = kone(cx_ii);
      auto ktwo_c = ktwo(cx_ii);

      auto kone_phi = kone(phi_i);
      auto ktwo_phi = ktwo(phi_i);

      auto Kcx = kone_c + ktwo_c;
      auto Kphi = kone_phi + ktwo_phi;

      auto result_cx =
          gaxpy_oop(2.0, Jxc.compress(), -1.0, Kcx.compress(), true);
      auto result_phi =
          gaxpy_oop(2.0, Jx_phi.compress(), -1.0, Kphi.compress(), true);

      return std::make_pair(result_cx, result_phi);
    };

    auto compute_vbc_i = [&](const vector_real_function_3d &bx,
                             const vector_real_function_3d &by,
                             const vector_real_function_3d &c,
                             const vector_real_function_3d &phi,
                             const real_function_3d &v, int j) {
      auto [gbc_x, gbc_phi_x] = compute_g_i(bx, by, phi, c[j], phi[j]);
      auto vcx = v * c[j];
      auto fbx = -1.0 * Q(gbc_x + vcx);
      // This is the tricky bit, Here we need to compute the vector of functions
      // and multiply by phi[i] and take the trace of each to get a 1D tensor to
      // multiple by cx auto vb_phi = mul(world, v, phi, true);
      auto vb_phi = v * phi[j];
      auto fb_phi_x = mul(world, gbc_phi_x + vb_phi, phi);
      // this phi*(gbc + vb)_p is a vector of functions

      Tensor<double> m_fbi_row = Tensor<double>(fb_phi_x.size());
      for (int i = 0; i < fb_phi_x.size(); i++) {
        m_fbi_row[i] = fb_phi_x[i].trace();
      }
      auto num_orbitals = static_cast<int>(phi0.size());
      m_fbi_row.reshape({num_orbitals, 1});
      auto fphi_x = transform(world, c, m_fbi_row, true);

      auto result_x = fbx + fphi_x[0];
      return result_x;
    };
    int num_orbitals = static_cast<int>(phi0.size());
    auto bc_indexer = x_space_indexer(num_orbitals);
    auto zeta_indexer = response_space_index(num_orbitals);

    auto compute_result = [&](const int &i) {
      resultT result = zero_functions_compressed<double, 3>(world, 1);
      auto bi = b[i];
      auto ci = c[i];
      auto response_orbital_index = i % (num_orbitals * 2);
      bool type_x = true;
      if (response_orbital_index < num_orbitals) {
        type_x = true;
      } else {
        type_x = false;
        response_orbital_index -= num_orbitals;
      }

      const auto &bx = bc_indexer.get_x_state(bi, B);
      const auto &by = bc_indexer.get_y_state(bi, B);
      const auto &cx = bc_indexer.get_x_state(ci, C);
      const auto &cy = bc_indexer.get_y_state(ci, C);
      const auto &vb = dipole_perturbations[bi];
      const auto &vc = dipole_perturbations[ci];
      real_function_3d one, two;
      if (type_x) {
        one = compute_vbc_i(bx, by, cx, phi0, vb, response_orbital_index);
        two = compute_vbc_i(cx, cy, bx, phi0, vc, response_orbital_index);
      } else {
        one = compute_vbc_i(by, bx, cy, phi0, vb, response_orbital_index);
        two = compute_vbc_i(cy, cx, by, phi0, vc, response_orbital_index);
      }
      result[0] = gaxpy_oop(1.0, one, 1.0, two, true);
      return result;
    };
    const long ij = batch.result.begin;
    return compute_result(ij);
  }
};

class VBC_task_i : public MacroTaskOperationBase {

  class Partitioner : public MacroTaskPartitioner {
  public:
    Partitioner() {
      max_batch_size = 1;
      //  result_stride = stride;
      //  policy = "strided";
    }
  };

public:
  VBC_task_i() { partitioner.reset(new Partitioner()); }
  // index b, index c, response B, response C, zetaBC, zetaCB, phi0,
  // perturbation
  typedef std::tuple<
      const std::vector<int> &, const std::vector<int> &,
      const std::vector<int> &, const std::vector<int> &,
      const vector_real_function_3d &, const vector_real_function_3d &,
      const vector_real_function_3d &, const vector_real_function_3d &,
      const vector_real_function_3d &, const vector_real_function_3d &>
      argtupleT;

  using resultT = vector_real_function_3d;

  resultT allocator(World &world, const argtupleT &args) const {

    // In this example the first index is index of all orbitals we are working
    // with
    auto num_functions = static_cast<int>(std::get<0>(args).size());
    // auto num_orbitals = static_cast<int>(std::get<7>(args).size());
    // auto num_functions = 2 * num_states * num_orbitals;
    // print("allocator: ", num_states, num_orbitals, num_functions);

    return zero_functions_compressed<double, 3>(world, num_functions);
  };

  resultT
  operator()(const std::vector<int> &orb_index,
             const std::vector<int> &state_index, const std::vector<int> &b,
             const std::vector<int> &c, const vector_real_function_3d &B,
             const vector_real_function_3d &C,
             const vector_real_function_3d &zeta_BC,
             const vector_real_function_3d &zeta_CB,
             const vector_real_function_3d &phi0,
             const vector_real_function_3d &dipole_perturbations) const {

    World &world = phi0[0].world();
    madness::QProjector<double, 3> Q(phi0);
    auto thresh = FunctionDefaults<3>::get_thresh();

    auto K = [&](const vecfuncT &ket, const vecfuncT &bra) {
      const double lo = 1.e-10;
      auto &world = ket[0].world();
      Exchange<double, 3> k{world, lo};
      k.set_bra_and_ket(bra, ket);

      std::string algorithm_ = "small_memory";
      if (algorithm_ == "multiworld") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient);
      } else if (algorithm_ == "multiworld_row") {
        k.set_algorithm(
            Exchange<double, 3>::Algorithm::multiworld_efficient_row);
      } else if (algorithm_ == "largemem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::large_memory);
      } else if (algorithm_ == "smallmem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::small_memory);
      }

      return k;
    };
    // In this example we only apply once to either x or y orbital.  If x we
    // just reverse left and right in the call
    auto compute_g_i = [&](const vector_real_function_3d &Aleft,
                           const vector_real_function_3d &Aright,
                           const vector_real_function_3d &Bleft,
                           const vector_real_function_3d &Bright,
                           const real_function_3d &phix) {
      auto x_phi = mul(world, Aleft, Aright, true);
      auto y_phi = mul(world, Bleft, Bright, true);
      world.gop.fence();
      auto rho = sum(world, x_phi, true);
      rho += sum(world, y_phi, true);
      auto lo = 1.e-10;
      real_convolution_3d op =
          CoulombOperator(world, lo, FunctionDefaults<3>::get_thresh());
      auto temp_J = apply(op, rho);

      auto Jx = temp_J * phix;
      // mul(temp_J, phix, true);

      auto ka = K(Aleft, Aright)(phix); // what happens to k after this?
      auto kb = K(Bleft, Bright)(phix);

      // ideally it runs and the Exchange operator is freed
      auto Kx = ka + kb; // gaxpy_oop(1.0, ka, 1.0, kb, true);
      auto result_x = gaxpy_oop(2.0, Jx.compress(), -1.0, Kx.compress(), true);

      return result_x;
    };

    auto compute_vbc_i = [&](const vector_real_function_3d &bx,
                             const vector_real_function_3d &by,
                             const vector_real_function_3d &cx,
                             const vector_real_function_3d &cy,
                             const vector_real_function_3d &bx_zeta,
                             const vector_real_function_3d &cy_zeta,
                             const vector_real_function_3d &phi_bc,
                             const vector_real_function_3d &zeta_bc,
                             const vector_real_function_3d &phi,
                             const real_function_3d &v, int j) {
      // left right pairs ka and kb and apply apply
      // TODO: There maybe an optimization that can tell whether we can skip
      // this orbital or not

      auto gzeta_x = compute_g_i(bx_zeta, cy_zeta, phi_bc, zeta_bc, phi0[j]);
      gzeta_x = -1.0 * Q(gzeta_x);
      auto gbc_x = compute_g_i(bx, phi0, phi0, by, cx[j]);

      auto vcx = v * cx[j];

      auto fbx = -1.0 * Q(gbc_x + vcx);

      // This is the tricky bit, Here we need to compute the vector of functions
      // and multiply by phi[i] and take the trace of each to get a 1D tensor to
      // multiple by cx auto vb_phi = mul(world, v, phi, true);
      auto vb_phi = v * phi[j];
      auto gbc_phi_x = compute_g_i(bx, phi0, phi0, by, phi0[j]);
      auto fb_phi_x = mul(world, gbc_phi_x + vb_phi, phi);
      // this phi*(gbc + vb)_p is a vector of functions
      // We next need trace of each function to form the 1D tensor of
      // coefficients

      Tensor<double> m_fbi_row = Tensor<double>(fb_phi_x.size());
      for (int i = 0; i < fb_phi_x.size(); i++) {
        m_fbi_row[i] = fb_phi_x[i].trace();
      }

      auto num_orbitals = static_cast<int>(phi0.size());
      m_fbi_row.reshape({num_orbitals, 1});

      auto fphi_x = transform(world, cx, m_fbi_row, true);

      auto result_x = gzeta_x + fbx + fphi_x[0];
      return result_x.truncate();
    };
    int num_orbitals = static_cast<int>(phi0.size());
    auto bc_indexer = x_space_indexer(num_orbitals);
    auto zeta_indexer = response_space_index(num_orbitals);
    auto compute_result = [&](const int &i) {
      resultT result = zero_functions_compressed<double, 3>(world, 1);
      // Now I should assume that each b and c will have all of the
      // corresponding [00000,00000,0000,11111,11111,11111,] or I find a way to
      // intelligently map i to bi and ci
      auto bi = b[i];
      auto ci = c[i];
      auto si = state_index[i];

      auto response_orbital_index = i % (num_orbitals * 2);
      bool type_x = true;
      if (response_orbital_index < num_orbitals) {
        type_x = true;
      } else {
        type_x = false;
        response_orbital_index -= num_orbitals;
      }

      const auto &bx = bc_indexer.get_x_state(bi, B);
      const auto &by = bc_indexer.get_y_state(bi, B);
      const auto &cx = bc_indexer.get_x_state(ci, C);
      const auto &cy = bc_indexer.get_y_state(ci, C);

      const auto &zeta_bc = zeta_indexer.get_x_state(si, zeta_BC);
      const auto &zeta_cb = zeta_indexer.get_x_state(si, zeta_CB);

      const auto &vb = dipole_perturbations[bi];
      const auto &vc = dipole_perturbations[ci];
      real_function_3d one, two;
      if (type_x) {
        one = compute_vbc_i(bx, by, cx, cy, bx, cy, phi0, zeta_bc, phi0, vb,
                            response_orbital_index);
        two = compute_vbc_i(cx, cy, bx, by, cx, by, phi0, zeta_cb, phi0, vc,
                            response_orbital_index);
      } else {
        one = compute_vbc_i(by, bx, cy, cx, cy, bx, zeta_bc, phi0, phi0, vb,
                            response_orbital_index);
        two = compute_vbc_i(cy, cx, by, bx, by, cx, zeta_cb, phi0, phi0, vc,
                            response_orbital_index);
      }

      result[0] = gaxpy_oop(1.0, one, 1.0, two, true);
      auto norm_r = result[0].norm2();

      print("VBC_task_i: ", i, "bi: ", bi, "ci: ", ci, "si: ", si,
            "response_orbital_index: ", response_orbital_index,
            "type_x: ", type_x, " norm: ", norm_r);

      return result;
    };
    const long ij = batch.result.begin;
    return compute_result(ij);
  }
};

// class ComputeBetaTask : public MacroTaskOperationBase {

//   int dummy_variable = 24;

//   class Partitioner : public MacroTaskPartitioner {
//   public:
//     Partitioner() {
//       max_batch_size = 1;
//       //  result_stride = stride;
//       //  policy = "strided";
//     }
//   };

// public:
//   ComputeBetaTask() { partitioner.reset(new Partitioner()); }

//   // beta_abc = -2 * <bx|cy * va> + <cx|by * va> + <zeta_bc|phi0 * va> +
//   // <zeta_cb|phi0 * va> + <ax|vbcx> + <ay|vbcy> beta[i], a, b, c, A, B, C,
//   // zeta_BC, zeta_CB, phi0, dipole_perturbations, vbc
//   typedef std::tuple<const std::vector<int> &,        // i
//                      const std::vector<int> &,        // a
//                      const std::vector<int> &,        // b
//                      const std::vector<int> &,        // c
//                      const vector_real_function_3d &, // A
//                      const vector_real_function_3d &, // B
//                      const vector_real_function_3d &, // C
//                      const vector_real_function_3d &, // zeta_BC
//                      const vector_real_function_3d &, // zeta_CB
//                      const vector_real_function_3d &, // phi0
//                      const vector_real_function_3d &, // dipole_perturbations
//                      const vector_real_function_3d &  // vbc
//                      >
//       argtupleT;

//   typedef std::vector<std::shared_ptr<ScalarResult<double>>> resultT;

//   resultT allocator(World &world, const argtupleT &args) const {
//     std::size_t n = std::get<0>(args).size();
//     return scalar_result_shared_ptr_vector<double>(world, n);
//   };

//   resultT operator()(const std::vector<int> &i, const std::vector<int> &a,
//                      const std::vector<int> &b, const std::vector<int> &c,
//                      const vector_real_function_3d &A,
//                      const vector_real_function_3d &B,
//                      const vector_real_function_3d &C,
//                      const vector_real_function_3d &zeta_BC,
//                      const vector_real_function_3d &zeta_CB,
//                      const vector_real_function_3d &phi0,
//                      const vector_real_function_3d &dipole_perturbations,
//                      const vector_real_function_3d &vbc) const

//   {

//     World &world = B[0].world();
//     auto result = scalar_result_shared_ptr_vector<double>(world, 1);

//     auto compute_beta_i = [&](const vector_real_function_3d &ax,
//                               const vector_real_function_3d &ay,
//                               const vector_real_function_3d &bx,
//                               const vector_real_function_3d &by,
//                               const vector_real_function_3d &cx,
//                               const vector_real_function_3d &cy,
//                               const vector_real_function_3d &zeta_bc,
//                               const vector_real_function_3d &zeta_cb,
//                               const vector_real_function_3d &phi0,
//                               const real_function_3d &va,
//                               const vector_real_function_3d &vbcx,
//                               const vector_real_function_3d &vbcy) {
//       auto one = dot(world, bx, cy * va);
//       auto three = dot(world, cx, by * va);

//       auto two = dot(world, zeta_bc, phi0 * va);
//       auto four = dot(world, zeta_cb, phi0 * va);

//       auto five = dot(world, ax, vbcx, true);
//       auto six = dot(world, ay, vbcy, true);

//       auto one_trace = one.trace();
//       auto two_trace = two.trace();
//       auto three_trace = three.trace();
//       auto four_trace = four.trace();
//       auto five_trace = five.trace();
//       auto six_trace = six.trace();

//       auto beta = one_trace + two_trace + three_trace + four_trace +
//                   five_trace + six_trace;
//       *result[0] = -2.0 * beta;
//       return result;
//     };

//     std::cout << "compute_beta: " << dummy_variable << std::endl;

//     const long m = batch.result.begin;

//     auto ak = a[m];
//     auto bk = b[m];
//     auto ck = c[m];
//     auto ik = i[m];

//     auto num_orbitals = static_cast<int>(phi0.size());

//     x_space_indexer x_space_indexer(num_orbitals);
//     response_space_index response_space_index(num_orbitals);
//     auto ax = x_space_indexer.get_x_state(ak, A);
//     auto ay = x_space_indexer.get_y_state(ak, A);
//     auto bx = x_space_indexer.get_x_state(bk, B);
//     auto by = x_space_indexer.get_y_state(bk, B);
//     auto cx = x_space_indexer.get_x_state(ck, C);
//     auto cy = x_space_indexer.get_y_state(ck, C);
//     auto zeta_bc = response_space_index.get_x_state(ik, zeta_BC);
//     auto zeta_cb = response_space_index.get_x_state(ik, zeta_CB);
//     auto va = dipole_perturbations[ak];
//     auto vbcx = x_space_indexer.get_x_state(ik, vbc);
//     auto vbcy = x_space_indexer.get_y_state(ik, vbc);

//     return compute_beta_i(ax, ay, bx, by, cx, cy, zeta_bc, zeta_cb, phi0, va,
//                           vbcx, vbcy);
//   }
// };
class ComputeZetaBC : public MacroTaskOperationBase {

private:
  vector_real_function_3d phi0;
  int num_orbitals;

  class Partitioner : public MacroTaskPartitioner {
  public:
    Partitioner(int stride) {
      max_batch_size = 1;
      result_stride = stride;
      policy = "strided";
    }
  };

public:
  ComputeZetaBC(int num_orbitals) : num_orbitals(num_orbitals) {
    partitioner.reset(new Partitioner(num_orbitals));
  }
  // index b, index c, response B, response C, zetaBC, zetaCB, phi0,
  // perturbation
  typedef std::tuple<const std::vector<int> &, const std::vector<int> &,
                     const vector_real_function_3d &,
                     const vector_real_function_3d &,
                     const vector_real_function_3d &>
      argtupleT;

  using resultT = vector_real_function_3d;

  resultT allocator(World &world, const argtupleT &args) const {

    auto num_states = static_cast<int>(std::get<0>(args).size() * num_orbitals);
    auto num_orbitals = static_cast<int>(std::get<4>(args).size());
    auto num_functions = num_states * num_orbitals;

    return zero_functions_compressed<double, 3>(world, num_functions);
  };

  resultT operator()(const std::vector<int> &b, const std::vector<int> &c,
                     const vector_real_function_3d &B,
                     const vector_real_function_3d &C,
                     const vector_real_function_3d &phi0) const {

    World &world = phi0[0].world();
    auto bc_indexer = x_space_indexer(num_orbitals);

    const long i = batch.result.begin;

    auto bi = b[i];
    auto ci = c[i];

    const auto &by = bc_indexer.get_y_state(bi, B);
    const auto &cx = bc_indexer.get_x_state(ci, C);

    auto mbc = matrix_inner(world, by, cx);
    return -1.0 * transform(world, phi0, mbc, true);
  }
};

class ResponseComputeGammaX : public MacroTaskOperationBase {

  class Partitioner : public MacroTaskPartitioner {
  public:
    Partitioner() {
      max_batch_size = 1;
      //  result_stride = stride;
      //  policy = "strided";
    }
  };

public:
  ResponseComputeGammaX() { partitioner.reset(new Partitioner()); }
  // index b,response B, phi0,
  // perturbation
  typedef std::tuple<const std::vector<int> &,        // orb_index
                     const std::vector<int> &,        // state_index
                     const vector_real_function_3d &, // B
                     const vector_real_function_3d &, // phi0
                     const bool &>
      argtupleT;

  using resultT = vector_real_function_3d;

  resultT allocator(World &world, const argtupleT &args) const {

    // In this example the first index is index of all orbitals we are working
    // with
    auto num_functions = static_cast<int>(std::get<0>(args).size());
    return zero_functions_compressed<double, 3>(world, num_functions);
  };

  resultT operator()(const std::vector<int> &orb_index,
                     const std::vector<int> &state_index,
                     const vector_real_function_3d &B,
                     const vector_real_function_3d &phi0,
                     const bool &static_frequency) const {

    World &world = phi0[0].world();
    madness::QProjector<double, 3> Q(phi0);
    auto thresh = FunctionDefaults<3>::get_thresh();

    auto K = [&](const vecfuncT &ket, const vecfuncT &bra) {
      const double lo = 1.e-10;
      auto &world = ket[0].world();
      Exchange<double, 3> k{world, lo};
      k.set_bra_and_ket(bra, ket);

      std::string algorithm_ = "small_memory";
      if (algorithm_ == "multiworld") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient);
      } else if (algorithm_ == "multiworld_row") {
        k.set_algorithm(
            Exchange<double, 3>::Algorithm::multiworld_efficient_row);
      } else if (algorithm_ == "largemem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::large_memory);
      } else if (algorithm_ == "smallmem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::small_memory);
      }

      return k;
    };
    // In this example we only apply once to either x or y orbital.  If x we
    // just reverse left and right in the call
    auto compute_g1_i = [&](const vector_real_function_3d &Aleft,
                            const vector_real_function_3d &Aright,
                            const vector_real_function_3d &Bleft,
                            const vector_real_function_3d &Bright,
                            const real_function_3d &phix) {
      auto x_phi = mul(world, Aleft, Aright, true);
      auto y_phi = mul(world, Bleft, Bright, true);
      world.gop.fence();
      auto rho = sum(world, x_phi, true);
      rho += sum(world, y_phi, true);
      auto lo = 1.e-10;
      real_convolution_3d op =
          CoulombOperator(world, lo, FunctionDefaults<3>::get_thresh());
      auto temp_J = apply(op, rho);
      auto Jx = temp_J * phix;
      auto ka = K(Aleft, Aright)(phix); // what happens to k after this?
      auto kb = K(Bleft, Bright)(phix);
      auto Kx = ka + kb; // gaxpy_oop(1.0, ka, 1.0, kb, true);
      auto result_x = gaxpy_oop(2.0, Jx.compress(), -1.0, Kx.compress(), true);

      return result_x;
    };

    int num_orbitals = static_cast<int>(phi0.size());
    auto compute_result = [&](const int &i) {
      resultT result = zero_functions_compressed<double, 3>(world, 1);
      // Now I should assume that each b and c will have all of the
      // corresponding [00000,00000,0000,11111,11111,11111,] or I find a way to
      // intelligently map i to bi and ci
      auto si = state_index[i];

      real_function_3d g;
      if (static_frequency) {
        // If static then we pass a response function and only work with x
        auto indexer = response_space_index(num_orbitals);
        const auto &bx = indexer.get_x_state(si, B);
        auto response_orbital_index = i % (num_orbitals);
        // else we do g1
        g = compute_g1_i(bx, phi0, phi0, bx, phi0[response_orbital_index]);
      } else {

        auto indexer = x_space_indexer(num_orbitals);
        const auto &bx = indexer.get_x_state(si, B);
        const auto &by = indexer.get_y_state(si, B);
        auto response_orbital_index = i % (num_orbitals * 2);

        bool type_x = true;
        if (response_orbital_index < num_orbitals) {
          type_x = true;
        } else {
          type_x = false;
          response_orbital_index -= num_orbitals;
        }

        if (type_x) {
          g = compute_g1_i(bx, phi0, phi0, by, phi0[response_orbital_index]);
        } else {
          g = compute_g1_i(by, phi0, phi0, bx, phi0[response_orbital_index]);
        }
      }
      g = Q(g);
      result[0] = g;

      return result;
    };
    const int ij = static_cast<int>(batch.result.begin);
    return compute_result(ij);
  }
};

class ResponseComputeGroundExchange : public MacroTaskOperationBase {

  class Partitioner : public MacroTaskPartitioner {
  public:
    Partitioner() {
      max_batch_size = 1;
      //  result_stride = stride;
      //  policy = "strided";
    }
  };

public:
  ResponseComputeGroundExchange() { partitioner.reset(new Partitioner()); }
  // index b,response B, phi0,
  // perturbation
  typedef std::tuple<const std::vector<int> &,        // orb_index
                     const std::vector<int> &,        // state_index
                     const vector_real_function_3d &, // B
                     const vector_real_function_3d &, // phi0
                     const bool &>
      argtupleT;

  using resultT = vector_real_function_3d;

  resultT allocator(World &world, const argtupleT &args) const {

    // In this example the first index is index of all orbitals we are working
    // with
    auto num_functions = static_cast<int>(std::get<0>(args).size());
    return zero_functions_compressed<double, 3>(world, num_functions);
  };

  resultT operator()(const std::vector<int> &orb_index,
                     const std::vector<int> &state_index,
                     const vector_real_function_3d &B,
                     const vector_real_function_3d &phi0,
                     const bool &static_frequency) const {

    World &world = phi0[0].world();
    madness::QProjector<double, 3> Q(phi0);
    auto thresh = FunctionDefaults<3>::get_thresh();

    auto K = [&](const vecfuncT &ket, const vecfuncT &bra) {
      const double lo = 1.e-10;
      auto &world = ket[0].world();
      Exchange<double, 3> k{world, lo};
      k.set_bra_and_ket(bra, ket);

      std::string algorithm_ = "small_memory";
      if (algorithm_ == "multiworld") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient);
      } else if (algorithm_ == "multiworld_row") {
        k.set_algorithm(
            Exchange<double, 3>::Algorithm::multiworld_efficient_row);
      } else if (algorithm_ == "largemem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::large_memory);
      } else if (algorithm_ == "smallmem") {
        k.set_algorithm(Exchange<double, 3>::Algorithm::small_memory);
      }

      return k;
    };
    auto compute_g0_i = [&](const vector_real_function_3d &Aleft,
                            const vector_real_function_3d &Aright,
                            const real_function_3d &x) {
      auto Kx = K(Aleft, Aright)(x); // what happens to k after this?
      // ideally it runs and the Exchange operator is freed

      return Kx;
    };
    // In this example we only apply once to either x or y orbital.  If x we
    // just reverse left and right in the call

    int num_orbitals = static_cast<int>(phi0.size());
    auto compute_result = [&](const int &i) {
      resultT result = zero_functions_compressed<double, 3>(world, 1);
      // Now I should assume that each b and c will have all of the
      // corresponding [00000,00000,0000,11111,11111,11111,] or I find a way to
      // intelligently map i to bi and ci
      auto si = state_index[i];

      real_function_3d g;
      if (static_frequency) {
        // If static then we pass a response function and only work with x
        auto indexer = response_space_index(num_orbitals);
        const auto &bx = indexer.get_x_state(si, B);
        auto response_orbital_index = i % (num_orbitals);
        // if it's zero then we do g0
        g = compute_g0_i(phi0, phi0, bx[response_orbital_index]);
      } else {

        auto indexer = x_space_indexer(num_orbitals);
        const auto &bx = indexer.get_x_state(si, B);
        const auto &by = indexer.get_y_state(si, B);
        auto response_orbital_index = i % (num_orbitals * 2);

        bool type_x = true;
        if (response_orbital_index < num_orbitals) {
          type_x = true;
        } else {
          type_x = false;
          response_orbital_index -= num_orbitals;
        }

        if (type_x) {
          g = compute_g0_i(phi0, phi0, bx[response_orbital_index]);
        } else {
          g = compute_g0_i(phi0, phi0, by[response_orbital_index]);
        }
      }
      //        g = Q(g);
      result[0] = g;

      return result;
    };
    const int ij = static_cast<int>(batch.result.begin);
    return compute_result(ij);
  }
};

class ResponseApplyBSH : public MacroTaskOperationBase {

  class Partitioner : public MacroTaskPartitioner {
  public:
    Partitioner() {
      max_batch_size = 1;
      //  result_stride = stride;
      //  policy = "strided";
    }
  };

public:
  ResponseApplyBSH() { partitioner.reset(new Partitioner()); }

  typedef SeparatedConvolution<double, 3> operatorT;
  typedef std::shared_ptr<operatorT> poperatorT;
  // index b,response B, phi0,
  // perturbation
  typedef std::tuple<const vector_real_function_3d &, const Tensor<double> &,
                     const double &, const double &, const bool &>
      argtupleT;

  using resultT = vector_real_function_3d;

  resultT allocator(World &world, const argtupleT &args) const {

    // In this example the first index is index of all orbitals we are working
    // with
    auto num_functions = static_cast<int>(std::get<0>(args).size());
    return zero_functions_compressed<double, 3>(world, num_functions);
  };

  resultT operator()(const vector_real_function_3d &x,
                     const Tensor<double> &ground_energies, const double omega,
                     const double e_shift, const bool static_frequency) const {

    World &world = x[0].world();
    auto tol = FunctionDefaults<3>::get_thresh();
    auto num_orbitals = ground_energies.size();
    double lo = 1.e-10;

    auto make_op = [&](const int i, auto omega, auto e_shift) {
      auto mu =
          sqrt(-2.0 * (ground_energies(i % num_orbitals) + omega + e_shift));
      auto bsh = BSHOperatorPtr3D(world, mu, lo, tol);
      return bsh;
    };

    const int ij = static_cast<int>(batch.result.begin);
    resultT result = zero_functions_compressed<double, 3>(world, 1);
    if (static_frequency) {
      auto bsh = make_op(ij, omega, e_shift);
      result[0] = apply(*bsh, x[ij]);
    } else {

      auto response_orbital_index = ij % (num_orbitals * 2);

      if (response_orbital_index >= num_orbitals) {
        auto bsh = make_op(ij, -omega, 0);
        result[0] = apply(*bsh, x[ij]);
      } else {
        auto bsh = make_op(ij, omega, e_shift);
        result[0] = apply(*bsh, x[ij]);
      }
    }

    return result;
  }
};

} // namespace madness
#endif
