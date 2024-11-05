// Copyright 2021 Adrian Hurtado
#ifndef SRC_APPS_MOLRESPONSE_X_SPACE_H_
#define SRC_APPS_MOLRESPONSE_X_SPACE_H_

#include <madness/mra/mra.h>
#include <madness/tensor/tensor.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <numeric>
#include <vector>

#include "functypedefs.h"
#include "molresponse/response_functions.h"

namespace madness {

typedef std::vector<vector_real_function_3d> response_matrix;

struct X_space;

auto to_response_vector(const vector_real_function_3d &vec)
    -> vector_real_function_3d;
auto create_response_matrix(const size_t &num_state,
                            const size_t &num_orbitals) -> response_matrix;
auto to_response_matrix(const X_space &x) -> response_matrix;
auto to_conjugate_response_matrix(const X_space &x) -> response_matrix;
auto to_flattened_vector(const X_space &x) -> vector_real_function_3d;
auto to_X_space(const response_matrix &x) -> X_space;
auto to_conjugate_X_space(const response_matrix &x) -> X_space;

struct X_space {
private:
  size_t n_states;   // Num. of resp. states
  size_t n_orbitals; // Num. of ground states

public:
  response_space x, y;
  std::list<size_t> active;

public:
  [[nodiscard]] size_t num_states() const { return n_states; }
  [[nodiscard]] size_t num_orbitals() const { return n_orbitals; }
  // default constructor
  X_space() : n_states(0), n_orbitals(0), x(), y(), active(0) {}
  // Copy constructor
  void reset_active() {
    active.resize(n_states);
    x.active.resize(n_states);
    y.active.resize(n_states);
    size_t i{0};
    for (auto &ai : active) {
      ai = i++;
    }
    i = 0;
    for (auto &ai : x.active) {
      ai = i++;
    }
    i = 0;
    for (auto &ai : y.active) {
      ai = i++;
    }
  }
  void set_active(const std::list<size_t> &new_active) {
    active = new_active;
    x.active = new_active;
    y.active = new_active;
  }
  X_space(const X_space &A)
      : n_states(size_states(A)), n_orbitals(size_orbitals(A)), x(A.x), y(A.y),
        active(A.active) {}
  [[nodiscard]] X_space copy() const {
    auto &world = x[0][0].world();
    auto new_x = X_space(*this); // copy
    for (const auto &i : new_x.active) {
      new_x.x[i] = madness::copy(world, x[i], false);
      new_x.y[i] = madness::copy(world, y[i], false);
    }
    world.gop.fence();

    return new_x;
  }
  /// Create a new copy of the function with different distribution and optional
  /// fence

  /// Works in either basis.  Different distributions imply
  /// asynchronous communication and the optional fence is
  /// collective.
  [[nodiscard]] auto
  copy(const std::shared_ptr<WorldDCPmapInterface<Key<3>>> &p_map,
       bool fence = false) const -> X_space {
    auto &world = x[0][0].world();
    auto new_x = X_space(*this); // copy
    for (int i = 0; i < new_x.num_states(); i++) {
      new_x.x[i] = madness::copy(world, x[i], p_map, false);
      new_x.y[i] = madness::copy(world, y[i], p_map, false);
    }
    world.gop.fence();
    return new_x;
  }
  // assignment
  auto operator=(const X_space &B) -> X_space & {
    if (this != &B) { // is it the same object?

      this->n_states = B.num_states();
      this->n_orbitals = B.num_orbitals();
      this->x = B.x;
      this->y = B.y;
      this->active = B.active;
    }
    return *this; // NO SHALLOW COPIES
  }
  X_space(World &world, size_t n_states, size_t n_orbitals)
      : n_states(n_states), n_orbitals(n_orbitals),
        x(world, n_states, n_orbitals), y(world, n_states, n_orbitals),
        active(n_states) {
    reset_active();
  }
  void clear() {
    x.clear();
    y.clear();
    active.clear();
  }
  void push_back(const vector_real_function_3d &vx,
                 const vector_real_function_3d &vy) {
    if (n_orbitals > 0) {
      MADNESS_ASSERT(n_orbitals == vx.size());
      MADNESS_ASSERT(n_orbitals == vy.size());
      MADNESS_ASSERT(vx.size() == vy.size());
    } else { // g_states == 0 (empty vector)
      n_orbitals = vx.size();
    }
    MADNESS_ASSERT(vx.size() == num_orbitals());
    MADNESS_ASSERT(vy.size() == num_orbitals());
    active.push_back(active.back() + 1);
    n_states++;
    x.push_back(vx);
    y.push_back(vy);
    // Be smart with g_states
  }
  void pop_back() {
    x.pop_back();
    y.pop_back();
    active.pop_back();
    n_states--;
    if (n_states == 0) {
      n_orbitals = 0;
    }
  }

  friend auto inplace_apply(
      X_space &A,
      const std::function<void(vector_real_function_3d &)> &func) -> void {
    auto &world = A.x[0][0].world();
    for (auto &i : A.active) {
      func(A.x[i]);
      func(A.y[i]);
    }
    world.gop.fence();
  }

  /**
   * @brief Apply a function to the X_space
   * @param A
   * @param func
   * @return
   */
  friend auto oop_apply(const X_space &A,
                        const std::function<vector_real_function_3d(
                            const vector_real_function_3d &)> &func,
                        bool fence = true) -> X_space {
    auto &world = A.x[0][0].world();
    auto result =
        X_space::zero_functions(world, A.num_states(), A.num_orbitals());
    result.set_active(A.active);
    //    if (world.rank() == 0) { print("oop_apply"); }
    for (auto &i : result.active) {
      //       if (world.rank() == 0) { print("oop_apply", i); }
      result.x[i] = func(A.x[i]);
      result.y[i] = func(A.y[i]);
    }
    if (fence)
      world.gop.fence();
    return result;
  }

  template <typename T>
  friend auto binary_apply(const X_space &A, const X_space &B,
                           T &func) -> X_space {
    MADNESS_ASSERT(same_size(A, B));

    auto &world = A.x[0][0].world();
    X_space result =
        X_space::zero_functions(world, A.num_states(), A.num_orbitals());
    result.set_active(A.active);

    for (const auto &i : result.active) {
      auto ax = A.x[i];
      auto bx = B.x[i];

      auto ay = A.y[i];
      auto by = B.y[i];

      result.x[i] = func(ax, bx);
      result.y[i] = func(ay, by);
    }
    world.gop.fence();
    return result;
  }

  template <class T>
  friend auto binary_inplace(X_space &A, const X_space &B, const T &func) {
    MADNESS_ASSERT(same_size(A, B));
    auto &world = A.x[0][0].world();
    for (const auto &i : A.active) {
      auto ax = A.x[i];
      auto ay = A.y[i];

      auto bx = B.x[i];
      auto by = B.y[i];

      func(ax, bx);
      func(ay, by);
    }
    world.gop.fence();

    return A;
  }

  static X_space zero_functions(World &world, size_t n_states,
                                size_t n_orbitals) {
    auto zeros = X_space(world, n_states, n_orbitals);
    for (int i = 0; i < zeros.num_states(); i++) {
      zeros.x[i] =
          ::madness::zero_functions<double, 3>(world, n_orbitals, true);
      zeros.y[i] =
          ::madness::zero_functions<double, 3>(world, n_orbitals, true);
    }
    world.gop.fence();
    return zeros;
  }

  auto operator+=(const X_space &B) -> X_space & {
    MADNESS_ASSERT(same_size(*this, B));
    auto &world = this->x[B.active.front()][0].world();
    this->active = B.active;
    this->from_vector(this->to_vector() + B.to_vector());

    // auto add_inplace = [&](auto &a, const auto &b)
    // {
    //   gaxpy(world, 1.0, a, 1.0, b, true);
    // };
    // binary_inplace(*this, B, add_inplace);
    return *this;
  }

  friend auto operator+(const X_space &A, const X_space &B) -> X_space {
    MADNESS_ASSERT(same_size(A, B));

    auto result = X_space(A.x[A.active.front()][0].world(), A.num_states(),
                          A.num_orbitals());
    result.set_active(A.active);
    result.from_vector(A.to_vector() + B.to_vector());
    return result;

    auto add_ab = [&](const auto &a, const auto &b) {
      return gaxpy_oop(1.0, a, 1.0, b, true);
    };
    return binary_apply(A, B, add_ab);
  }

  friend X_space operator-(const X_space &A, const X_space &B) {
    MADNESS_ASSERT(same_size(A, B));

    auto result = X_space(A.x[A.active.front()][0].world(), A.num_states(),
                          A.num_orbitals());
    result.set_active(A.active);
    result.from_vector(A.to_vector() - B.to_vector());
    return result;

    auto sub_ab = [&](const auto &a, const auto &b) {
      return gaxpy_oop(1.0, a, -1.0, b, true);
    };
    return binary_apply(A, B, sub_ab);
  }

  friend X_space operator*(const X_space &A, const double &b) {
    World &world = A.x[A.active.front()][0].world();

    auto result = X_space(world, A.num_states(), A.num_orbitals());
    result.set_active(A.active);
    result.from_vector(A.to_vector() * b);
    return result;

    auto scale_a = [&](vector_real_function_3d &vec_ai) {
      scale(world, vec_ai, b, true);
    };
    inplace_apply(result, scale_a);
    return result;
  }
  friend X_space operator*(const double &b, const X_space &A) {
    World &world = A.x[A.active.front()][0].world();
    auto result = X_space(world, A.num_states(), A.num_orbitals());
    result.set_active(A.active);
    result.from_vector(A.to_vector() * b);
    return result;

    auto scale_a = [&](vector_real_function_3d &vec_ai) {
      scale(world, vec_ai, b, true);
    };
    inplace_apply(result, scale_a);
    return result;
  }
  friend X_space operator*(const X_space &B, const X_space &A) {
    World &world = A.x[A.active.front()][0].world();
    auto result = X_space(A.x[0][0].world(), A.num_states(), A.num_orbitals());
    result.set_active(A.active);
    vector_real_function_3d result_vec =
        mul(world, A.to_vector(), B.to_vector());
    result.from_vector(result_vec);
    return result;

    auto mul_ab = [&](vector_real_function_3d &vec_ai,
                      vector_real_function_3d &vec_bi) {
      return mul(world, vec_ai, vec_bi, false);
    };
    return binary_apply(A, B, mul_ab);
  }

  friend X_space operator*(const X_space &A, const Function<double, 3> &f) {
    World &world = A.x[A.active.front()][0].world();

    auto result = X_space(A.x[0][0].world(), A.num_states(), A.num_orbitals());
    result.set_active(A.active);
    result.from_vector(A.to_vector() * f);
    return result;

    auto mul_f = [&](const vector_real_function_3d &vec_ai) {
      return mul(world, f, vec_ai, false);
    };
    return oop_apply(A, mul_f);
  }
  friend auto operator*(const Function<double, 3> &f,
                        const X_space &A) -> X_space {
    World &world = A.x[A.active.front()][0].world();
    auto result = X_space(A.x[0][0].world(), A.num_states(), A.num_orbitals());
    result.set_active(A.active);
    result.from_vector(f * A.to_vector());
    return result;

    auto mul_f = [&](const vector_real_function_3d &vec_ai) {
      return mul(world, f, vec_ai, false);
    };
    return oop_apply(A, mul_f);
  }

  friend auto operator*(const X_space &A, const Tensor<double> &b) -> X_space {
    MADNESS_ASSERT(size_states(A) > 0);
    MADNESS_ASSERT(size_orbitals(A) > 0);

    World &world = A.x[0][0].world();
    auto transform_ai = [&](auto &ai) { return transform(world, ai, b, true); };
    return oop_apply(A, transform_ai);
  }
  /***
   *
   * @param A
   * @param B
   * @return
   */
  friend auto inner(const X_space &A, const X_space &B) -> Tensor<double>;

  void truncate() {

    auto &world = this->x[x.active.front()][0].world();
    this->from_vector(madness::truncate(
        this->to_vector(), FunctionDefaults<3>::get_thresh(), true));
  }

  void compress() const {

    auto x = to_vector();
    auto &world = this->x[active.front()][0].world();
    madness::compress(world, x, true);
  }
  void reconstruct() const {

    auto x = to_vector();
    auto &world = this->x[active.front()][0].world();
    madness::reconstruct(world, x, true);
  }

  void truncate(double thresh) {

    auto &world = this->x[x.active.front()][0].world();
    this->from_vector(madness::truncate(this->to_vector(), thresh, true));
  }

  auto norm2s() const -> Tensor<double> {
    World &world = x[0][0].world();
    Tensor<double> norms(num_states());

    auto x = to_response_matrix(*this);
    int b = 0;
    for (const auto &xb : x) {
      norms[b++] = norm2(world, xb);
    }
    world.gop.fence();
    return norms;
  }

  [[nodiscard]] auto component_norm2s() const -> Tensor<double> {
    World &world = x[0][0].world();
    auto rx = to_flattened_vector(*this);
    auto norms = norm2s_T(world, rx);
    return norms.reshape(n_states, 2 * n_orbitals);
  }

  friend auto size_states(const X_space &x) -> size_t { return x.n_states; }
  friend auto size_orbitals(const X_space &x) -> size_t { return x.n_orbitals; }
  friend auto same_size(const X_space &A, const X_space &B) -> bool {
    return ((size_states(A) == size_states(B) &&
             size_orbitals(A) == size_orbitals(B)));
  }

  [[nodiscard]] auto to_vector() const -> vector_real_function_3d {

    int n = static_cast<int>(active.size());
    int m = static_cast<int>(num_orbitals());

    vector_real_function_3d rf(2 * n * m);

    int i = 0;
    for (const auto &ai : active) {
      for (int j = 0; j < m; j++) {
        auto xindex = (2 * i * m) + j;
        auto yindex = (2 * i * m) + j + m;
        rf[xindex] = x[ai][j];
        rf[yindex] = y[ai][j];
      }
      i++;
    }
    return rf;
  }

  auto from_vector(const vector_real_function_3d &rf) -> void {

    int m = static_cast<int>(num_orbitals());

    int i = 0;
    for (const auto &ai : active) {
      for (int j = 0; j < m; j++) {
        auto xindex = (2 * i * m) + j;
        auto yindex = (2 * i * m) + j + m;

        x[ai][j] = rf[xindex];
        y[ai][j] = rf[yindex];
      }
      i++;
    }
  }
};

// but the solver needs the functions initialized to zero for which we also need
// the world object.

struct X_vector : public X_space {
  X_vector(World &world, size_t n_orbtials) {
    this->X_space::zero_functions(world, size_t(1), n_orbtials);
  }

  X_vector(X_space A, size_t b)
      : X_space(A.x[0][0].world(), size_t(1), A.num_orbitals()) {
    x[0] = A.x[b];
    y[0] = A.y[b];
  }
  friend X_vector operator-(const X_vector &A, const X_vector &B) {
    MADNESS_ASSERT(same_size(A, B));

    World &world = A.x[0][0].world();
    X_vector result(world, size_orbitals(A)); // create zero_functions
    result.x = A.x - B.x;
    result.y = A.y - B.y;
    return result;
  }
  friend X_vector operator*(const X_vector &A, const double &c) {
    World &world = A.x[0][0].world();
    X_vector result(world, size_orbitals(A)); // create zero_functions
    result.x = A.x * c;
    result.y = A.y * c;
    return result;
  }
  X_vector copy() const {
    X_vector copyX(x[0][0].world(), x.num_orbitals);
    copyX.x = x.copy();
    copyX.y = y.copy();
    return copyX;
  }
  auto operator+=(const X_vector &B) -> X_vector & {
    MADNESS_ASSERT(same_size(*this, B));
    this->x += B.x;
    this->y += B.y;
    return *this;
  }
  inline friend auto inner(X_vector &A, X_vector &B) -> double {
    MADNESS_ASSERT(size_states(A) == 1);
    MADNESS_ASSERT(size_orbitals(A) > 0);
    MADNESS_ASSERT(same_size(A, B));

    Tensor<double> G(1, 1);
    Tensor<double> G1(1, 1);
    Tensor<double> G2(1, 1);

    World &world = A.x[0][0].world();

    auto ax = madness::copy(world, A.x[0]);
    auto ay = madness::copy(world, A.y[0]);

    auto bx = madness::copy(world, B.x[0]);
    auto by = madness::copy(world, B.y[0]);

    for (auto &ayi : ay) {
      ax.push_back(madness::copy(ayi));
    }
    for (auto &byi : by) {
      bx.push_back(madness::copy(byi));
    };

    double result = inner(ax, bx);

    return result;
  }
};
// function object with allocator()()
struct response_matrix_allocator {
  World &world;
  const size_t n_orbtials;
  response_matrix_allocator(World &world, size_t n_orbtials)
      : world(world), n_orbtials(n_orbtials) {}
  // overloading the default constructor () operator
  vector_real_function_3d operator()() {
    // print("allocator called with ", int(n_orbtials), " orbitals");
    //  returning constructor of x_vector
    return zero_functions<double, 3>(world, n_orbtials);
  }
  // Copy constructor

  response_matrix_allocator operator=(const response_matrix_allocator &other) {
    return response_matrix_allocator(world, other.n_orbtials);
  }
};

struct response_function_allocator {
  World &world;
  response_function_allocator(World &world) : world(world) {}
  // overloading the default constructor () operator
  real_function_3d operator()() {
    return real_function_3d(real_factory_3d(world).fence(true));
  }
  response_function_allocator
  operator=(const response_function_allocator &other) {
    return response_function_allocator(world);
  }
};

vector_real_function_3d copyToVector(const X_space &chi);

void copyToXspace(const vector_real_function_3d &rf, X_space &chi);
// In this implmentation we need to represent each x_space as a contigous block
// of functions.
vector_real_function_3d copyToVector(const response_space &chi);

void copyToResponseSpace(const vector_real_function_3d &rf,
                         response_space &chi);

class x_space_indexer {
  int num_orbitals;

public:
  x_space_indexer(int num_orbitals) : num_orbitals(num_orbitals) {}

  [[nodiscard]] vector_real_function_3d
  get_x_state(int i, const vector_real_function_3d &rf) const {

    vector_real_function_3d subset(num_orbitals);
    auto index = 2 * num_orbitals * i;

    for (int j = 0; j < num_orbitals; j++) {
      subset[j] = rf[index + j];
    }
    return subset;
  }
  [[nodiscard]] vector_real_function_3d
  get_y_state(int i, const vector_real_function_3d &rf) const {

    vector_real_function_3d subset(num_orbitals);
    auto index = 2 * num_orbitals * i + num_orbitals;
    for (int j = 0; j < num_orbitals; j++) {
      subset[j] = rf[index + j];
    }
    return subset;
  }
};

class response_space_index {
  int num_orbitals;

public:
  response_space_index(int num_orbitals) : num_orbitals(num_orbitals) {}

  [[nodiscard]] vector_real_function_3d
  get_x_state(int i, const vector_real_function_3d &rf) const {

    vector_real_function_3d subset(num_orbitals);
    auto index = num_orbitals * i;

    for (int j = 0; j < num_orbitals; j++) {
      subset[j] = rf[index + j];
    }
    return subset;
  }
};
} // namespace madness

#endif // SRC_APPS_MOLRESPONSE_X_SPACE_H_
