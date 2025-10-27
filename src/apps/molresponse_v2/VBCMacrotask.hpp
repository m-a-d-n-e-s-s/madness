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
class SimpleVBCComputer {
public:
  SimpleVBCComputer(World &world, const GroundStateData &gs)
      : world_(world), gs_(gs), num_orbitals_(gs.getNumOrbitals()),
        spin_restricted_(gs.isSpinRestricted()) {}

  std::pair<ResponseVector, ResponseVector>
  get_BC_vecs(const VBCResponseState &vbc_state) {
    auto [stateB, stateC] = vbc_state.get_states();
    // 2) Load their response vectors:
    ResponseVector Bv, Cv;
    load_response_vector(world_, num_orbitals_, stateB, /*th*/ 0, /*f*/ 0, Bv);
    load_response_vector(world_, num_orbitals_, stateC, /*th*/ 0, /*f*/ 0, Cv);

    // 3) Promote to dynamic (x/y) if needed:
    auto Bdyn = make_response_vector(num_orbitals_, /*is_static=*/false,
                                     spin_restricted_);
    auto Cdyn = make_response_vector(num_orbitals_, /*is_static=*/false,
                                     spin_restricted_);
    promote_response_vector(world_, Bv, Bdyn);
    promote_response_vector(world_, Cv, Cdyn);

    return std::make_pair(Bdyn, Cdyn);
  }

  static vector_real_function_3d
  make_zeta_bc(World &world, const vector_real_function_3d &by,
               const vector_real_function_3d &cx,
               const vector_real_function_3d &phi0) {
    auto mat = matrix_inner(world, by, cx);
    return -1.0 * transform(world, phi0, mat, true);
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

    auto zeta_bc = make_zeta_bc(world, by, cx, gs.orbitals);
    auto zeta_cb = make_zeta_bc(world, cy, bx, gs.orbitals);

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

  /// Given a fully–described VBCResponseState, load or compute its vector.
  ResponseVector compute_and_save(const VBCResponseState &vbc_state) {
    auto filename = vbc_state.response_filename();
    // (MADness will append “.00000” etc.)
    if (fs::exists(filename + ".00000")) {
      ResponseVector loaded;
      if (world_.rank() == 0)
        print("📂 Loading existing VBC from", filename);
      load_response_vector(world_, num_orbitals_, vbc_state, 0, 0, loaded);
      return loaded;
    }

    auto [xB, xC] = get_BC_vecs(vbc_state);
    auto [stateB, stateC] = vbc_state.get_states();

    // 2) Load their response vectors:
    auto VB = perturbation_vector(world_, gs_, stateB);
    auto VC = perturbation_vector(world_, gs_, stateC);
    VB.insert(VB.end(), VB.begin(), VB.end());
    VC.insert(VC.end(), VC.begin(), VC.end());
    auto vb = raw_perturbation_operator(world_, gs_, stateB.perturbation);
    auto vc = raw_perturbation_operator(world_, gs_, stateC.perturbation);

    auto result = compute_vbc_response(world_, gs_, xB, xC, VB, VC, vb, vc);

    // 5) Save to disk and return:
    save_response_vector(world_, vbc_state, result);
    if (world_.rank() == 0)
      print("💾 Wrote VBC to", filename);
    return result;
  }

private:
  World &world_;
  const GroundStateData &gs_;
  int num_orbitals_;
  bool spin_restricted_;
};
} // namespace madness
