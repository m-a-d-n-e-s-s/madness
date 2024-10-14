// #include "functypedefs.h"
// #include "functypedefs.h"
#include "vmra.h"
#include "x_space.h"
#include <madness.h>
#include <madness/chem/SCFOperators.h>
#include <madness/chem/projector.h>
#include <madness/mra/macrotaskpartitioner.h>
#include <madness/mra/macrotaskq.h>
#include <madness/world/cloud.h>
#include <tuple>

namespace madness
{

  // In this implmentation we need to represent each x_space as a contigous block
  // of functions.
  vector_real_function_3d copyToVector(const X_space &chi)
  {

    int n = static_cast<int>(chi.num_states());
    int m = static_cast<int>(chi.num_orbitals());

    vector_real_function_3d rf(2 * n * m);

    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < m; j++)
      {
        auto xindex = (2 * i * m) + j;
        auto yindex = (2 * i * m) + j + m;
        rf[xindex] = chi.x[i][j];
        rf[yindex] = chi.y[i][j];
      }
    }
    return rf;
  }

  void copyToXspace(const vector_real_function_3d &rf, X_space &chi)
  {

    int n = static_cast<int>(chi.num_states());
    int m = static_cast<int>(chi.num_orbitals());
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < m; j++)
      {
        auto xindex = (2 * i * m) + j;
        auto yindex = (2 * i * m) + j + m;
        chi.x[i][j] = rf[xindex];
        chi.y[i][j] = rf[yindex];
      }
    }
  }

  // In this implmentation we need to represent each x_space as a contigous block
  // of functions.
  vector_real_function_3d copyToVector(const response_space &chi)
  {

    int n = static_cast<int>(chi.num_states);
    int m = static_cast<int>(chi.num_orbitals);

    vector_real_function_3d rf(n * m);

    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < m; j++)
      {
        auto xindex = (i * m) + j;
        rf[xindex] = chi.x[i][j];
      }
    }
    return rf;
  }

  void copyToResponseSpace(const vector_real_function_3d &rf,
                           response_space &chi)
  {

    int n = static_cast<int>(chi.num_states);
    int m = static_cast<int>(chi.num_orbitals);
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < m; j++)
      {
        auto xindex = (i * m) + j;
        chi.x[i][j] = rf[xindex];
      }
    }
  }

  class x_space_indexer
  {
    int num_orbitals;

  public:
    x_space_indexer(int num_orbitals)
        : num_orbitals(num_orbitals) {}

    [[nodiscard]] vector_real_function_3d
    get_x_state(int i, const vector_real_function_3d &rf) const
    {

      vector_real_function_3d subset(num_orbitals);
      auto index = 2 * num_orbitals * i;

      for (int j = 0; j < num_orbitals; j++)
      {
        subset[j] = rf[index + j];
      }
      return subset;
    }
    [[nodiscard]] vector_real_function_3d
    get_y_state(int i, const vector_real_function_3d &rf) const
    {

      vector_real_function_3d subset(num_orbitals);
      auto index = 2 * num_orbitals * i + num_orbitals;
      for (int j = 0; j < num_orbitals; j++)
      {
        subset[j] = rf[index + j];
      }
      return subset;
    }
  };

  class response_space_index
  {
    int num_orbitals;

  public:
    response_space_index(int num_orbitals)
        : num_orbitals(num_orbitals) {}

    [[nodiscard]] vector_real_function_3d
    get_x_state(int i, const vector_real_function_3d &rf) const
    {

      vector_real_function_3d subset(num_orbitals);
      auto index = num_orbitals * i;

      for (int j = 0; j < num_orbitals; j++)
      {
        subset[j] = rf[index + j];
      }
      return subset;
    }
  };

  // Next we need to write a function which takes a pair of x and y functions and
  // returns a x_space

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

  class VBC_task : public MacroTaskOperationBase
  {
  public:
    // index b, index c, response B, response C, zetaBC, zetaCB, phi0,
    // perturbation
    typedef std::tuple<
        const std::vector<int> &, const std::vector<int> &,
        const vector_real_function_3d &, const vector_real_function_3d &,
        const vector_real_function_3d &, const vector_real_function_3d &,
        const vector_real_function_3d &, const vector_real_function_3d &>
        argtupleT;

    using resultT = vector_real_function_3d;

    resultT allocator(World &world, const argtupleT &args) const
    {

      auto num_states = static_cast<int>(std::get<0>(args).size());
      auto num_orbitals = static_cast<int>(std::get<6>(args).size());
      print("allocator: ", num_states, num_orbitals);

      return zero_functions_compressed<double, 3>(world,
                                                  2 * num_states * num_orbitals);
    };

    resultT
    operator()(const std::vector<int> &b, const std::vector<int> &c,
               const vector_real_function_3d &B, const vector_real_function_3d &C,
               const vector_real_function_3d &zeta_BC,
               const vector_real_function_3d &zeta_CB,
               const vector_real_function_3d &phi0,
               const vector_real_function_3d &dipole_perturbations) const
    {

      World &world = phi0[0].world();
      madness::QProjector<double, 3> Q(world, phi0);
      auto thresh = FunctionDefaults<3>::get_thresh();

      auto K = [&](const vecfuncT &ket, const vecfuncT &bra)
      {
        const double lo = 1.e-10;
        auto &world = ket[0].world();
        Exchange<double, 3> k{world, lo};
        k.set_bra_and_ket(bra, ket);

        std::string algorithm_ = "smallmem";

        if (algorithm_ == "multiworld")
        {
          k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient);
        }
        else if (algorithm_ == "multiworld_row")
        {
          k.set_algorithm(
              Exchange<double, 3>::Algorithm::multiworld_efficient_row);
        }
        else if (algorithm_ == "largemem")
        {
          k.set_algorithm(Exchange<double, 3>::Algorithm::large_memory);
        }
        else if (algorithm_ == "smallmem")
        {
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
                           const vector_real_function_3d &phiy)
      {
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
                               const real_function_3d &v)
      {
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

      auto bc_indexer = x_space_indexer(num_orbitals);
      auto zeta_indexer = response_space_index(num_orbitals);
      vector_real_function_3d results(2 * num_states * num_orbitals);

      for (int i = 0; i < num_states; i++)
      {
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

        for (int j = 0; j < num_orbitals; j++)
        {
          results[(2 * i * num_orbitals) + j] = result_x[j];
          results[(2 * i * num_orbitals) + j + num_orbitals] = result_y[j];
        }
      }
      return results;
    }
  };
  class VBC_task2 : public MacroTaskOperationBase
  {

    class Partitioner : public MacroTaskPartitioner
    {
    public:
      Partitioner(long stride)
      {
        max_batch_size = 1;
        result_stride = stride;
        policy = "strided";
      }
    };

  public:
    long stride;
    VBC_task2(long stride)
        : stride(stride)
    {
      partitioner.reset(new Partitioner(stride));
    }
    // index b, index c, response B, response C, zetaBC, zetaCB, phi0,
    // perturbation
    typedef std::tuple<
        const std::vector<int> &,
        const std::vector<int> &, const std::vector<int> &,
        const vector_real_function_3d &, const vector_real_function_3d &,
        const vector_real_function_3d &, const vector_real_function_3d &,
        const vector_real_function_3d &, const vector_real_function_3d &>
        argtupleT;

    using resultT = vector_real_function_3d;

    resultT allocator(World &world, const argtupleT &args) const
    {

      auto num_states = static_cast<int>(std::get<1>(args).size());
      auto num_orbitals = static_cast<int>(std::get<7>(args).size());
      auto num_functions = 2 * num_states * num_orbitals;
      print("allocator: ", num_states, num_orbitals, num_functions);

      return zero_functions_compressed<double, 3>(world,
                                                  num_functions);
    };

    resultT
    operator()(
        const std::vector<int> &dummy_i,
        const std::vector<int> &b, const std::vector<int> &c,
        const vector_real_function_3d &B, const vector_real_function_3d &C,
        const vector_real_function_3d &zeta_BC,
        const vector_real_function_3d &zeta_CB,
        const vector_real_function_3d &phi0,
        const vector_real_function_3d &dipole_perturbations) const
    {

      World &world = phi0[0].world();
      madness::QProjector<double, 3> Q(world, phi0);
      auto thresh = FunctionDefaults<3>::get_thresh();

      auto K = [&](const vecfuncT &ket, const vecfuncT &bra)
      {
        const double lo = 1.e-10;
        auto &world = ket[0].world();
        Exchange<double, 3> k{world, lo};
        k.set_bra_and_ket(bra, ket);

        std::string algorithm_ = "smallmem";

        if (algorithm_ == "multiworld")
        {
          k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient);
        }
        else if (algorithm_ == "multiworld_row")
        {
          k.set_algorithm(
              Exchange<double, 3>::Algorithm::multiworld_efficient_row);
        }
        else if (algorithm_ == "largemem")
        {
          k.set_algorithm(Exchange<double, 3>::Algorithm::large_memory);
        }
        else if (algorithm_ == "smallmem")
        {
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
                           const vector_real_function_3d &phiy)
      {
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
      auto compute_vbc_i = [&](
                               const vector_real_function_3d &bx,
                               const vector_real_function_3d &by,
                               const vector_real_function_3d &cx,
                               const vector_real_function_3d &cy,
                               const vector_real_function_3d &zeta_bc,
                               const vector_real_function_3d &phi,
                               const real_function_3d &v)
      {
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
      auto compute_result = [&](const int &i)
      {
        auto results = zero_functions_compressed<double, 3>(world,
                                                            2 * num_orbitals);

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

        for (int j = 0; j < num_orbitals; j++)
        {

          auto index_x = j;
          auto norm_x = result_x[j].norm2();
          print("i,j = (", i, j, ") index_x: ", index_x, " norm_x: ", norm_x);
          results[index_x] = result_x[j];
        }

        for (int j = 0; j < num_orbitals; j++)
        {

          auto index_y = j + num_orbitals;
          auto norm_y = result_y[j].norm2();
          print("i,j = (", i, j, ") index_y: ", index_y, " norm_y: ", norm_y);
          results[index_y] = result_y[j];
        }

        return results;
      };
      const long ij = batch.result.begin;
      return compute_result(ij);
    }
  };

  class molresponse_macrotask_helper
  {

  private:
    int num_states;
    int num_orbitals;
    bool full_response;
    int full_response_factor;

  public:
    molresponse_macrotask_helper(int num_states, int num_orbitals, bool full_response) : num_states(num_states), num_orbitals(num_orbitals), full_response(full_response)
    {
      if (full_response)
      {
        full_response_factor = 2;
      }
      else
      {
        full_response_factor = 1;
      }
    }

    [[nodiscard]] std::tuple<int, int, int, std::string> get_response_indexes(int index) const
    {

      int state_index = index / (num_orbitals * full_response_factor);
      int response_orbital_index = index % (num_orbitals * full_response_factor);
      if (response_orbital_index < num_orbitals)
      {
        return std::make_tuple(state_index, response_orbital_index, response_orbital_index, "x");
      }
      else
      {
        return std::make_tuple(state_index, response_orbital_index, response_orbital_index - index, "y");
      }
    }
  };

  // In this implementation we are going to work with a single index at a time.
  // The dummy vector needs to be the same size as result index and we
  // will need to select the state the orbital i corresponds to.
  // therefore we need to map the index i to the corresponding input and output states which also
  // change between x and y.
  // There will also be an optimization possible which sends to each subworld the necessary data

  class VBC_task_i : public MacroTaskOperationBase
  {

    class Partitioner : public MacroTaskPartitioner
    {
    public:
      Partitioner()
      {
        max_batch_size = 1;
        //  result_stride = stride;
        //  policy = "strided";
      }
    };

  public:
    VBC_task_i()
    {
      partitioner.reset(new Partitioner());
    }
    // index b, index c, response B, response C, zetaBC, zetaCB, phi0,
    // perturbation
    typedef std::tuple<
        const std::vector<int> &,
        const std::vector<int> &,
        const std::vector<int> &, const std::vector<int> &,
        const vector_real_function_3d &, const vector_real_function_3d &,
        const vector_real_function_3d &, const vector_real_function_3d &,
        const vector_real_function_3d &, const vector_real_function_3d &>
        argtupleT;

    using resultT = vector_real_function_3d;

    resultT allocator(World &world, const argtupleT &args) const
    {

      // In this example the first index is index of all orbitals we are working with
      auto num_functions = static_cast<int>(std::get<0>(args).size());
      // auto num_orbitals = static_cast<int>(std::get<7>(args).size());
      // auto num_functions = 2 * num_states * num_orbitals;
      // print("allocator: ", num_states, num_orbitals, num_functions);

      return zero_functions_compressed<double, 3>(world,
                                                  num_functions);
    };

    resultT
    operator()(
        const std::vector<int> &orb_index,
        const std::vector<int> &state_index,
        const std::vector<int> &b, const std::vector<int> &c,
        const vector_real_function_3d &B, const vector_real_function_3d &C,
        const vector_real_function_3d &zeta_BC,
        const vector_real_function_3d &zeta_CB,
        const vector_real_function_3d &phi0,
        const vector_real_function_3d &dipole_perturbations) const
    {

      World &world = phi0[0].world();
      madness::QProjector<double, 3> Q(world, phi0);
      auto thresh = FunctionDefaults<3>::get_thresh();

      auto K = [&](const vecfuncT &ket, const vecfuncT &bra)
      {
        const double lo = 1.e-10;
        auto &world = ket[0].world();
        Exchange<double, 3> k{world, lo};
        k.set_bra_and_ket(bra, ket);

        std::string algorithm_ = "smallmem";

        if (algorithm_ == "multiworld")
        {
          k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient);
        }
        else if (algorithm_ == "multiworld_row")
        {
          k.set_algorithm(
              Exchange<double, 3>::Algorithm::multiworld_efficient_row);
        }
        else if (algorithm_ == "largemem")
        {
          k.set_algorithm(Exchange<double, 3>::Algorithm::large_memory);
        }
        else if (algorithm_ == "smallmem")
        {
          k.set_algorithm(Exchange<double, 3>::Algorithm::small_memory);
        }

        return k;
      };
      // In this example we only apply once to either x or y orbital.  If x we just reverse left and right in the call
      auto compute_g_i = [&](const vector_real_function_3d &Aleft,
                             const vector_real_function_3d &Aright,
                             const vector_real_function_3d &Bleft,
                             const vector_real_function_3d &Bright,
                             const real_function_3d &phix)
      {
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

      auto compute_g_vec = [&](const vector_real_function_3d &Aleft,
                               const vector_real_function_3d &Aright,
                               const vector_real_function_3d &Bleft,
                               const vector_real_function_3d &Bright,
                               const vector_real_function_3d &phix)
      {
        auto x_phi = mul(world, Aleft, Aright, true);
        auto y_phi = mul(world, Bleft, Bright, true);
        world.gop.fence();
        auto rho = sum(world, x_phi, true);
        rho += sum(world, y_phi, true);
        auto lo = 1.e-10;
        real_convolution_3d op =
            CoulombOperator(world, lo, FunctionDefaults<3>::get_thresh());
        auto temp_J = apply(op, rho);

        auto Jx = mul(world, temp_J, phix, false);

        auto ka = K(Aleft, Aright)(phix); // what happens to k after this?
        auto kb = K(Bleft, Bright)(phix);
        // ideally it runs and the Exchange operator is freed
        auto Kx = gaxpy_oop(1.0, ka, 1.0, kb, true);

        auto result_x = gaxpy_oop(2.0, Jx, -1.0, Kx, true);

        return result_x;
      };
      auto compute_vbc_i = [&](
                               const vector_real_function_3d &bx,
                               const vector_real_function_3d &by,
                               const vector_real_function_3d &cx,
                               const vector_real_function_3d &cy,
                               const vector_real_function_3d &bx_zeta,
                               const vector_real_function_3d &cy_zeta,
                               const vector_real_function_3d &phi_bc,
                               const vector_real_function_3d &zeta_bc,
                               const vector_real_function_3d &phi,
                               const real_function_3d &v,
                               int j)
      {
        // left right pairs ka and kb and apply apply

        auto gzeta_x = compute_g_i(bx_zeta, cy_zeta, phi_bc, zeta_bc, phi0[j]);
        gzeta_x = -1.0 * Q(gzeta_x);
        auto gbc_x = compute_g_i(bx, phi0, phi0, by, cx[j]);

        auto vcx = v * cx[j];

        auto fbx = -1.0 * Q(gbc_x + vcx);

        // This is the tricky bit, Here we need to compute the vector of functions and multiply by phi[i] and take the trace of each to get a 1D tensor to multiple by cx
        // auto vb_phi = mul(world, v, phi, true);
        auto vb_phi = v * phi[j];
        auto gbc_phi_x = compute_g_i(bx, phi0, phi0, by, phi0[j]);
        auto fb_phi_x = mul(world, gbc_phi_x + vb_phi, phi);
        // this phi*(gbc + vb)_p is a vector of functions
        // We next need trace of each function to form the 1D tensor of coefficients

        Tensor<double> m_fbi_row = Tensor<double>(fb_phi_x.size());
        for (int i = 0; i < fb_phi_x.size(); i++)
        {
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
      auto compute_result = [&](const int &i)
      {
        resultT result = zero_functions_compressed<double, 3>(world, 1);
        // Now I should assume that each b and c will have all of the corresponding [00000,00000,0000,11111,11111,11111,] or I find a way to intelligently map i to bi and ci
        auto bi = b[i];
        auto ci = c[i];
        auto si = state_index[i];

        auto response_orbital_index = i % (num_orbitals * 2);
        bool type_x = true;
        if (response_orbital_index < num_orbitals)
        {
          type_x = true;
        }
        else
        {
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
        // the difference will be that while i still call this function twice, the order of the arguments if it is x or y

        // i will be the orbital index of bx, by, cx, cy, zeta_bc, zeta_cb, phi0, vb,
        real_function_3d one, two;
        if (type_x)
        {
          one = compute_vbc_i(bx, by, cx, cy, bx, cy, phi0, zeta_bc, phi0, vb, response_orbital_index);
          two = compute_vbc_i(cx, cy, bx, by, cx, by, phi0, zeta_cb, phi0, vc, response_orbital_index);
        }
        else
        {
          // TODO how to flip zeta_bc properly?  I think I need to pass another phi0
          one = compute_vbc_i(by, bx, cy, cx, cy, bx, zeta_bc, phi0, phi0, vb, response_orbital_index);
          two = compute_vbc_i(cy, cx, by, bx, by, cx, zeta_cb, phi0, phi0, vc, response_orbital_index);
        }

        result[0] = gaxpy_oop(1.0, one, 1.0, two, true);
        auto norm_r = result[0].norm2();

        print("VBC_task_i: ", i, "bi: ", bi, "ci: ", ci, "si: ", si, "response_orbital_index: ", response_orbital_index, "type_x: ", type_x,
              " norm: ", norm_r);

        return result;
      };
      const long ij = batch.result.begin;
      return compute_result(ij);
    }
  };

} // namespace madness
