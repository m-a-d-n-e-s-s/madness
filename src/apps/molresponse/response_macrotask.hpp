//#include "functypedefs.h"
//#include "functypedefs.h"
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

// In this implmentation we need to represent each x_space as a contigous block of functions.
vector_real_function_3d copyToVector(const X_space& chi) {

  int n = static_cast<int>(chi.num_states());
  int m = static_cast<int>(chi.num_orbitals());

  vector_real_function_3d rf(2 * n * m);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      auto xindex = (2 * i * m) + j;
      auto yindex = (2 * i * m) + j + m;
      rf[xindex] = chi.x[i][j];
      rf[yindex] = chi.y[i][j];
    }
  }
  return rf;
}

void copyToXspace(const vector_real_function_3d& rf, X_space& chi) {

  int n = static_cast<int>(chi.num_states());
  int m = static_cast<int>(chi.num_orbitals());
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      auto xindex = (2 * i * m) + j;
      auto yindex = (2 * i * m) + j + m;
      chi.x[i][j] = rf[xindex];
      chi.y[i][j] = rf[yindex];
    }
  }
}

// In this implmentation we need to represent each x_space as a contigous block of functions.
vector_real_function_3d copyToVector(const response_space& chi) {

  int n = static_cast<int>(chi.num_states);
  int m = static_cast<int>(chi.num_orbitals);

  vector_real_function_3d rf(2 * n * m);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      auto xindex = (2 * i * m) + j;
      rf[xindex] = chi.x[i][j];
    }
  }
  return rf;
}

void copyToResponseSpace(const vector_real_function_3d& rf,
                         response_space& chi) {

  int n = static_cast<int>(chi.num_states);
  int m = static_cast<int>(chi.num_orbitals);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      auto xindex = (2 * i * m) + j;
      chi.x[i][j] = rf[xindex];
    }
  }
}

class x_space_indexer {
  int num_states;
  int num_orbitals;

 public:
  x_space_indexer(int num_states, int num_orbitals)
      : num_states(num_states), num_orbitals(num_orbitals) {}

  [[nodiscard]] vector_real_function_3d
  get_x_state(int i, const vector_real_function_3d& rf) const {

    vector_real_function_3d subset(num_orbitals);
    auto index = 2 * num_states * num_orbitals * i;

    for (int j = 0; j < num_orbitals; j++) {
      subset[j] = rf[index + j];
    }
    return subset;
  }
  [[nodiscard]] vector_real_function_3d
  get_y_state(int i, const vector_real_function_3d& rf) const {

    vector_real_function_3d subset(num_orbitals);
    auto index = 2 * num_states * num_orbitals * i + num_orbitals;

    for (int j = 0; j < num_orbitals; j++) {
      subset[j] = rf[index + j];
    }
    return subset;
  }
};

class response_space_index {
  int num_states;
  int num_orbitals;

 public:
  response_space_index(int num_states, int num_orbitals)
      : num_states(num_states), num_orbitals(num_orbitals) {}

  [[nodiscard]] vector_real_function_3d
  get_x_state(int i, const vector_real_function_3d& rf) const {

    vector_real_function_3d subset(num_orbitals);
    auto index = num_orbitals * i;

    for (int j = 0; j < num_orbitals; j++) {
      subset[j] = rf[index + j];
    }
    return subset;
  }
};

// Next we need to write a function which takes a pair of x and y functions and returns a x_space

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
//

class VBC_task : public MacroTaskOperationBase {
 public:
  // index b, index c, response B, response C, zetaBC, zetaCB, phi0, perturbation
  typedef std::tuple<
      const std::vector<int>&, const std::vector<int>&,
      const vector_real_function_3d&, const vector_real_function_3d&,
      const vector_real_function_3d&, const vector_real_function_3d&,
      const vector_real_function_3d&, const vector_real_function_3d&>
      argtupleT;

  using resultT = vector_real_function_3d;

  resultT allocator(World& world, const argtupleT& args) const {

    auto num_states = static_cast<int>(std::get<0>(args).size());
    auto num_orbitals = static_cast<int>(std::get<6>(args).size());

    return zero_functions_compressed<double, 3>(world,
                                                2 * num_states * num_orbitals);
  };

  resultT
  operator()(const std::vector<int>& b, const std::vector<int>& c,
             const vector_real_function_3d& B, const vector_real_function_3d& C,
             const vector_real_function_3d& zeta_BC,
             const vector_real_function_3d& zeta_CB,
             const vector_real_function_3d& phi0,
             const vector_real_function_3d& dipole_perturbations) const {

    World& world = phi0[0].world();
    madness::QProjector<double, 3> Q(world,phi0);
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
    auto compute_g = [&](const vector_real_function_3d& Aleft,
                         const vector_real_function_3d& Aright,
                         const vector_real_function_3d& Bleft,
                         const vector_real_function_3d& Bright,
                         const vector_real_function_3d& phix,
                         const vector_real_function_3d& phiy) {
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

      auto ka = K(Aleft, Aright)(phix);  // what happens to k after this?
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

    auto compute_vbc_i = [&](const vector_real_function_3d& bx,
                             const vector_real_function_3d& by,
                             const vector_real_function_3d& cx,
                             const vector_real_function_3d& cy,
                             const vector_real_function_3d& zeta_bc,
                             const vector_real_function_3d& phi,
                             const real_function_3d& v) {
      // left right pairs ka and kb and apply apply
      auto [gzeta_x, gzeta_y] = compute_g(cx, by, phi0, zeta_bc, phi0, phi0);
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

    auto bc_indexer = x_space_indexer(3, num_orbitals);
    auto zeta_indexer = response_space_index(10, num_orbitals);
    vector_real_function_3d results(2 * num_states * num_orbitals);

    for (int i = 0; i < num_states; i++) {
      auto bi = b[i];
      auto ci = c[i];

      const auto& bx = bc_indexer.get_x_state(bi, B);
      const auto& by = bc_indexer.get_y_state(bi, B);
      const auto& cx = bc_indexer.get_x_state(ci, C);
      const auto& cy = bc_indexer.get_y_state(ci, C);

      const auto& zeta_bc = zeta_indexer.get_x_state(i, zeta_BC);
      const auto& zeta_cb = zeta_indexer.get_x_state(i, zeta_CB);

      const auto& vb = dipole_perturbations[bi];
      const auto& vc = dipole_perturbations[ci];

      auto [bcx, bcy] = compute_vbc_i(bx, by, cx, cy, zeta_bc, phi0, vb);
      auto [cbx, cby] = compute_vbc_i(cx, cy, bx, by, zeta_cb, phi0, vc);

      auto result_x = gaxpy_oop(1.0, bcx, 1.0, cbx, true);
      auto result_y = gaxpy_oop(1.0, bcy, 1.0, cby, true);

      for (int j = 0; j < num_orbitals; j++) {
        results[(2 * i * num_orbitals) + j] = result_x[j];
        results[(2 * i * num_orbitals) + j + num_orbitals] = result_y[j];
      }
    }
    return results;
  }
};

}  // namespace madness
