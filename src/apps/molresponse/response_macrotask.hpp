
#include "vmra.h"
#include <madness.h>
#include <madness/chem/SCFOperators.h>
#include <madness/chem/projector.h>
#include <madness/mra/macrotaskpartitioner.h>
#include <madness/mra/macrotaskq.h>
#include <madness/world/cloud.h>
#include <tuple>

namespace madness {

// In response we have reponse densities
// which are made up 2 pairs of vectors of functions
// rho = x phi + phi y
//
// Here x and phi form a lr product
// and y and phi form a lr product
//
// For the coupled response equations we have to solve for
// the x and y functions which form the x y pairs
//
// For convience we define the following structures

struct response_lr_pair {
  vector_real_function_3d left;
  vector_real_function_3d right;
};

struct response_density {
  response_lr_pair x;
  response_lr_pair y;
};

struct response_xy_pair {
  vector_real_function_3d x;
  vector_real_function_3d y;
};

class compute_vbc : public MacroTaskOperationBase {
 public:
  // response B, response C, response zetaBC, phi0, perturbation
  typedef std::tuple<const response_density&, const response_xy_pair&,
                     const response_density&, const vector_real_function_3d&,
                     const real_function_3d&>
      argtupleT;

  using resultT = response_xy_pair;

  static resultT alloctor(World& world, const argtupleT& args) {
    int n = static_cast<int>(std::get<4>(args).size());
    return {zero_functions_compressed<double, 3>(world, n),
            zero_functions_compressed<double, 3>(world, n)};
  }

  resultT operator()(const response_density& B, const response_xy_pair& C,
                     const response_density& zeta_BC,
                     const vector_real_function_3d& phi0,
                     const real_function_3d& vb) const {

    World& world = phi0[0].world();
    madness::QProjector<double, 3> Q(phi0);
    auto thresh = FunctionDefaults<3>::get_thresh();

    auto K = [&](const vecfuncT& ket, const vecfuncT& bra) {
      const double lo = 1.e-10;
      auto& world = ket[0].world();
      Exchange<double, 3> k{world, lo};
      k.set_bra_and_ket(bra, ket);

      std::string algorithm_ = "smallmem";

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
    // This function constructs the J and K operators with A and B and applies on x
    auto compute_g = [&](const response_lr_pair& A, const response_lr_pair& B,
                         const response_xy_pair& phi) {
      auto x_phi = mul(world, A.left, A.right, true);
      auto y_phi = mul(world, B.left, B.right, true);
      world.gop.fence();
      auto rho = sum(world, x_phi, true);
      world.gop.fence();
      rho += sum(world, y_phi, true);
      world.gop.fence();
      real_convolution_3d op = CoulombOperator(world, 1.e-4, 1.e-5);
      auto temp_J = apply(op, rho);

      response_xy_pair J = {mul(world, temp_J, phi.x, true),
                            mul(world, temp_J, phi.y, true)};
      world.gop.fence();

      auto ka = K(A.left, A.right)(phi.x);  // what happens to k after this?
      world.gop.fence();
      auto kb = K(B.left, B.right)(phi.x);
      world.gop.fence();
      auto ka_conj = K(A.right, A.left)(phi.y);
      world.gop.fence();
      auto kb_conj = K(B.right, B.left)(phi.y);
      world.gop.fence();
      // ideally it runs and the Exchange operator is freed

      response_xy_pair K = {gaxpy_oop(1.0, ka, 1.0, kb, true),
                            gaxpy_oop(1.0, ka_conj, 1.0, kb_conj, true)};
      world.gop.fence();

      response_xy_pair results{gaxpy_oop(2.0, J.x, -1.0, K.x, true),
                               gaxpy_oop(2.0, J.y, -1.0, K.y, true)};
      world.gop.fence();
      return results;
    };
    auto gzeta = compute_g(zeta_BC.x, zeta_BC.y, {phi0, phi0});
    gzeta.x = -1.0 * Q(gzeta.x);
    gzeta.y = -1.0 * Q(gzeta.y);

    auto gBC = compute_g(B.x, B.y, {C.x, C.y});
    auto gBphi = compute_g(B.x, B.y, {phi0, phi0});
    response_xy_pair vbx = {mul(world, vb, C.x, true),
                            mul(world, vb, C.y, true)};
    world.gop.fence();
    response_xy_pair FBX = {-1.0 * Q(gBC.x + vbx.x), -1.0 * Q(gBC.y + vbx.y)};
    auto vb_phi0 = mul(world, vb, phi0, true);
    response_xy_pair FBphi0 = {gBphi.x + vb_phi0, gBphi.y + vb_phi0};

    auto fbx = matrix_inner(world, phi0, FBphi0.x);
    auto fby = matrix_inner(world, phi0, FBphi0.y);

    response_xy_pair FB = {transform(world, C.x, fbx, true),
                           transform(world, C.y, fby, true)};

    response_xy_pair results{truncate(gzeta.x + FBX.x + FB.x, thresh, true),
                             truncate(gzeta.y + FBX.y + FB.y, thresh, true)};
    return results;
  }
};

}  // namespace madness
