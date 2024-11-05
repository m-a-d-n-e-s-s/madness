/* Copyright 2021 Adrian Hurtado
 *   Small class to hold response functions and to interact with KAIN solver.
 */

#ifndef SRC_APPS_MOLRESPONSE_RESPONSE_FUNCTIONS_H_
#define SRC_APPS_MOLRESPONSE_RESPONSE_FUNCTIONS_H_

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <vector>

// real_function_3d == Function<double,3>
namespace madness {

typedef std::vector<vector_real_function_3d> response_matrix;

/* *
 * @brief response matrix holds response vectors for response state
 *
 */
struct response_space {
  // Member variables
  /**
   * @brief vector of vector of real 3d functions
   *
   */

public:
  size_t num_states;   // Num. of resp. states
  size_t num_orbitals; // Num. of ground states
  response_matrix x;
  std::list<size_t> active;

  // Member functions
public:
  /**
   * @brief default Construct a new response space object
   * num_states(0)
   * num_orbitals(0)
   * x() default constructor of std::vector
   */
  response_space() : num_states(0), num_orbitals(0), x(), active(0) {}

  void reset_active() {
    active.resize(num_states);
    size_t i{0};
    for (auto &ai : active) {
      ai = i++;
    }
  }

  // Copy constructor
  /**
   * @brief copy construct a new response space object
   * we are using copying defined by std:vector
   * we copy madness functions therefore we are copying pointers to function
   * implementations
   * @param y
   */
  response_space(const response_space &y)
      : num_states(y.size()), num_orbitals(y.size_orbitals()), x(y.x),
        active(y.active) {}

  // assignment
  // Copy assignment should copy the members of y and leave y Unchanged
  // The problem is what happens when we copy two functions
  response_space &operator=(const response_space &y) {
    //
    if (this != &y) { // is it the same object?
      this->num_states = y.size();
      this->num_orbitals = y.size_orbitals();
      this->x = y.x;
      this->active = y.active;
      if (x.size() != num_states) {
        x.resize(num_states);
      }
    }
    return *this; //
  }

  response_space &operator=(const response_matrix &y) {
    //
    this->num_states = y.size();
    this->num_orbitals = y[0].size();
    this->x = y;
    return *this; //
  }
  // Initialize functions to zero
  // m = number of response states
  // n = number of ground state orbitals

  // Zero Constructor constructs m vectors
  /**
   * @brief Construct a new response space with zero functions
   *
   * @param world
   * @param num_states
   * @param num_orbitals
   */
  response_space(World &world, size_t num_states, size_t num_orbitals)
      : num_states(num_states), num_orbitals(num_orbitals),
        x(response_matrix(num_states)), active(num_states) {
    for (auto &state : x) {
      state = vector_real_function_3d(num_orbitals);
    }
    reset_active();
    // world.gop.fence();
  }
  // Conversion from respones_matrix
  /**
   * @brief Construct a new response space object from vector of functions
   *
   * @param x
   */
  explicit response_space(const response_matrix &x)
      : num_states(x.size()), num_orbitals(x[0].size()), x(x),
        active(num_states) {
    reset_active();
  }

  // Determines if two ResponseFunctions are the same size
  friend bool same_size(const response_space &a, const response_space &b) {
    return ((a.size() == b.size()) && (a.size_orbitals() == b.size_orbitals()));
  }
  // Returns a deep copy
  [[nodiscard]] response_space copy() const {
    World &world = x[0][0].world();
    response_space result(*this);
    // copy each state
    for (size_t i = 0; i < num_states; i++) {
      result.x[i] = madness::copy(world, x[i], false);
    }
    world.gop.fence();
    return result;
  }
  [[nodiscard]] response_space
  copy(const std::shared_ptr<WorldDCPmapInterface<Key<3>>> &pmap,
       bool fence = false) const {
    auto &world = x[0][0].world();
    response_space result(*this);
    std::transform(x.begin(), x.end(), result.x.begin(), [&](const auto &xi) {
      return madness::copy(world, xi, pmap, fence);
    });
    world.gop.fence();
    return result;
  }

  vector_real_function_3d &operator[](size_t i) { return x.at(i); }
  const vector_real_function_3d &operator[](size_t i) const { return x.at(i); }

  friend auto inplace_unary_apply(
      response_space &A,
      const std::function<void(vector_real_function_3d &)> &func) {
    auto &world = A.x[A.active.front()][0].world();
    for (auto &i : A.active) {
      func(A.x[i]);
    }
    world.gop.fence();
  }

  friend auto oop_unary_apply(const response_space &A,
                              const std::function<vector_real_function_3d(
                                  const vector_real_function_3d &)> &func)
      -> response_space {
    auto result = A.copy();
    auto &world = result.x[A.active.front()][0].world();
    result.active = A.active;
    for (auto &i : result.active) {
      result.x[i] = func(A.x[i]);
    }
    world.gop.fence();
    return result;
  }

  friend auto
  binary_apply(const response_space &A, const response_space &B,
               const std::function<vector_real_function_3d(
                   vector_real_function_3d, vector_real_function_3d)> &func)
      -> response_space {
    MADNESS_ASSERT(same_size(A, B));

    response_space result = A.copy(); // create zero_functions
    auto &world = result.x[A.active.front()][0].world();

    for (const auto &i : result.active) {
      auto ax = result.x[i];
      auto bx = B.x[i];
      result.x[i] = func(ax, bx);
    }
    world.gop.fence();
    return result;
  }

  template <class T>
  friend auto binary_inplace(response_space &A, const response_space &B,
                             const T &func) {
    MADNESS_ASSERT(same_size(A, B));
    auto &world = A.x[A.active.front()][0].world();
    for (const auto &i : A.active) {
      auto ax = A.x[i];
      auto bx = B.x[i];
      func(ax, bx);
    }
    world.gop.fence();

    return A;
  }

  response_space operator+(const response_space &rhs_y) const {

    MADNESS_ASSERT(size() > 0);
    MADNESS_ASSERT(same_size(*this, rhs_y)); // assert that same size
    World &world = x[rhs_y.active.front()][0].world();
    auto result = response_space(world, size(), size_orbitals());
    result.active = rhs_y.active;

    result.from_vector(this->to_vector() + rhs_y.to_vector());
    return result;

    // auto result =
    //     binary_apply(*this, rhs_y, [&](auto xi, auto vi)
    //                  { return gaxpy_oop(1.0, xi, 1.0, vi, false); });
    // return result;
  }

  response_space operator-(const response_space &rhs_y) const {
    MADNESS_ASSERT(size() > 0);
    MADNESS_ASSERT(same_size(*this, rhs_y)); // assert that same size
    World &world = x[rhs_y.active.front()][0].world();
    auto result = this->copy();
    result.active = rhs_y.active;

    result.from_vector(result.to_vector() - rhs_y.to_vector());
    return result;

    // auto result =
    //     binary_apply(*this, rhs_y, [&](auto xi, auto vi)
    //                  { return gaxpy_oop(1.0, xi, -1.0, vi, false); });
    return result;
  }

  friend response_space operator*(const response_space &y, double a) {
    // World &world = y.x.at(0).at(0).world();
    World &world = y.x[y.active.front()][0].world();
    auto multiply_scalar = [&](vector_real_function_3d &vi) {
      madness::scale(world, vi, a, false);
    };

    // auto result = response_space(world, y.size(), y.size_orbitals());
    auto result =
        response_space::zero_functions(world, y.size(), y.size_orbitals());
    result.active = y.active;
    result.from_vector(y.to_vector() * a);
    return result;

    // auto result = y.copy();

    inplace_unary_apply(result, multiply_scalar);
    return result;
  }

  friend response_space operator*(double a, response_space &y) {
    // World &world = y.x.at(0).at(0).world();
    World &world = y.x[y.active.front()][0].world();

    // auto result = response_space(world, y.size(), y.size_orbitals());
    auto result =
        response_space::zero_functions(world, y.size(), y.size_orbitals());
    result.active = y.active;
    result.from_vector(y.to_vector() * a);
    return result;

    // auto multiply_scalar = [&](vector_real_function_3d &vi)
    // { madness::scale(world, vi, a, false); };
    // auto result = y.copy();
    // inplace_unary_apply(result, multiply_scalar);
    // return result;
  }

  response_space &operator*=(double a) {
    // World &world = this->x[0][0].world();
    World &world = x[active.front()][0].world();

    this->from_vector(this->to_vector() * a);
    // auto multiply_scalar = [&](vector_real_function_3d &vi)
    // { madness::scale(world, vi, a, false); };
    // inplace_unary_apply(*this, multiply_scalar);
    return *this;
  }

  // Scaling all internal functions by an external function
  // g[i][j] = x[i][j] * f
  friend response_space operator*(const response_space &a,
                                  const Function<double, 3> &f) {
    // World &world = a.x.at(0).at(0).world();
    World &world = a[a.active.front()][0].world();

    auto result =
        response_space::zero_functions(world, a.size(), a.size_orbitals());
    result.active = a.active;
    result.from_vector(a.to_vector() * f);
    return result;

    auto multiply_scalar_function = [&](const vector_real_function_3d &vi) {
      return mul(world, f, vi, false);
    };
    return oop_unary_apply(a, multiply_scalar_function);
  }

  // Scaling all internal functions by an external function
  // g[i][j] = x[i][j] * f
  friend response_space operator*(const Function<double, 3> &f,
                                  const response_space &a) {

    return a * f;
  }

  response_space operator*(const Function<double, 3> &f) {
    // World &world = x[0][0].world();
    World &world = x[active.front()][0].world();

    // auto result = response_space(world, size(), size_orbitals());
    auto result =
        response_space::zero_functions(world, size(), size_orbitals());
    result.active = active;
    result.from_vector(this->to_vector() * f);

    return result;

    auto multiply_scalar_function = [&](const vector_real_function_3d &vi) {
      return mul(world, f, vi, false);
    };

    return oop_unary_apply(*this, multiply_scalar_function);
  }

  friend response_space operator*(const response_space &a,
                                  const Tensor<double> &b) {
    MADNESS_ASSERT(a.size() > 0);
    MADNESS_ASSERT(!a[0].empty());
    World &world = a[0][0].world();
    auto response_transform = [&](const vector_real_function_3d &vi) {
      return transform(world, vi, b, true);
    };
    return oop_unary_apply(a, response_transform);
  }

  // KAIN must have this
  response_space &operator+=(const response_space &b) {
    MADNESS_ASSERT(same_size(*this, b));
    auto &world = b[b.active.front()][0].world();
    this->active = b.active;

    this->from_vector(this->to_vector() + b.to_vector());
    return *this;

    auto a_plus_equal_b = [&](vector_real_function_3d &a,
                              const vector_real_function_3d &g) {
      gaxpy(world, 1.0, a, 1.0, g, false);
    };
    binary_inplace(*this, b, a_plus_equal_b);
    return *this;

    return *this;
  }

  // Mimicking std::vector with these 4
  void push_back(const vector_real_function_3d &f) {
    x.push_back(f);
    num_states++;
    // print("im calling response_space push back");

    // Be smart with g_states
    if (num_orbitals > 0)
      MADNESS_ASSERT(num_orbitals == f.size());
    else { // g_states == 0 (empty vector)
      num_orbitals = f.size();
    }
  }

  void pop_back() {
    MADNESS_ASSERT(num_states >= 1);
    x.pop_back();
    num_states--;

    // Be smart with g_states
    if (num_states == 0) { // removed last item
      num_orbitals = 0;
    }
  }

  void clear() {
    x.clear();
    num_states = 0;
    num_orbitals = 0;
  }

  auto begin() { return x.begin(); }

  auto end() { return x.end(); }

  const auto begin() const { return x.begin(); }

  [[nodiscard]] const auto end() const { return x.end(); }

  size_t size() const { return num_states; }

  size_t size_orbitals() const { return num_orbitals; }

  // Mimicing standard madness calls with these 3
  void zero() {
    auto &world = x[0][0].world();
    for (int i = 0; i < num_states; i++) {
      x[i] = ::madness::zero_functions<double, 3>(world, num_orbitals, false);
    }
  }

  void compress_rf() const {
    // for (size_t k = 0; k < num_states; k++) { compress(x[0][0].world(), x[k],
    // true); } auto &world = x[0][0].world();
    auto &world = x[active.front()][0].world();
    // compress only active states
    //
    auto xvec = to_vector();
    compress(world, xvec, true);
  }

  void reconstruct_rf() {
    // for (size_t k = 0; k < num_states; k++) { reconstruct(x[0][0].world(),
    // x[k], true); } auto &world = x[0][0].world();
    auto &world = x[active.front()][0].world();
    // reconstruct only active states
    for (auto &i : active) {
      reconstruct(world, x[i], false);
    }
    world.gop.fence();
  }

  void truncate_rf() { truncate_rf(FunctionDefaults<3>::get_thresh()); }

  void truncate_rf(double tol) {
    auto &world = x[active.front()][0].world();
    // truncate only active states
    for (auto &i : active) {
      truncate(world, x[i], tol, false);
    }
    world.gop.fence();
    /*
    for (size_t k = 0; k < num_states; k++) { truncate(x[0][0].world(), x[k],
    tol, true); }
     */
  }

  // Returns norms of each state
  Tensor<double> norm2() {
    auto &world = x[0][0].world();
    Tensor<double> answer(num_states);
    for (size_t i = 0; i < num_states; i++) {
      answer(i) = madness::norm2(world, x[i]);
    }

    world.gop.fence();

    return answer;
  }

  // Scales each state (read: entire row) by corresponding vector element
  //     new[i] = old[i] * mat[i]
  void scale(Tensor<double> &mat) {
    for (size_t i = 0; i < num_states; i++)
      madness::scale(x[0][0].world(), x[i], mat[i], false);
    // x[i] = x[i] * mat[i];
  }

  friend bool operator==(const response_space &a, const response_space &y) {
    if (!same_size(a, y))
      return false;
    for (size_t b = 0; b < a.size(); ++b) {
      for (size_t k = 0; b < a.size_orbitals(); ++k) {
        if ((a[b][k] - y[b][k]).norm2() >
            FunctionDefaults<3>::get_thresh()) // this may be strict
          return false;
      }
    }
    return true;
  }

  static response_space zero_functions(World &world, size_t num_states,
                                       size_t num_orbitals) {
    response_space result(world, num_states, num_orbitals);

    for (int i = 0; i < num_states; i++) {
      result.x[i] =
          ::madness::zero_functions<double, 3>(world, num_orbitals, false);
    }
    return result;
  }

  friend Tensor<double> response_space_inner(const response_space &a,
                                             const response_space &b) {
    MADNESS_ASSERT(a.size() > 0);
    MADNESS_ASSERT(a.size() == b.size());
    MADNESS_ASSERT(!a[0].empty());
    MADNESS_ASSERT(a[0].size() == b[0].size());

    World &world = a[0][0].world();

    size_t dim_1 = a.size();
    size_t dim_2 = b[0].size();
    // Need to take transpose of each input ResponseFunction
    response_space aT(world, dim_2, dim_1);
    response_space bT(world, dim_2, dim_1);
    for (size_t i = 0; i < dim_1; i++) {
      for (size_t j = 0; j < dim_2; j++) {
        aT[j][i] = a[i][j];
        bT[j][i] = b[i][j];
      }
    }
    // Container for result
    Tensor<double> result(dim_1, dim_1);
    for (size_t p = 0; p < dim_2; p++) {
      result += matrix_inner(world, aT[p], bT[p]);
    }
    return result;
  }

  [[nodiscard]] auto to_vector() const -> vector_real_function_3d {

    int n = static_cast<int>(active.size());
    int m = static_cast<int>(num_orbitals);

    vector_real_function_3d rf(n * m);

    int i = 0;
    for (const auto &ai : active) {
      for (int j = 0; j < m; j++) {
        auto xindex = (i * m) + j;
        rf[xindex] = x[ai][j];
      }
      i++;
    }
    return rf;
  }

  auto from_vector(const vector_real_function_3d &rf) -> void {

    int m = static_cast<int>(num_orbitals);

    int i = 0;
    for (const auto &ai : active) {
      for (int j = 0; j < m; j++) {
        auto xindex = (i * m) + j;

        x[ai][j] = rf[xindex];
      }
      i++;
    }
  }
};

// Final piece for KAIN
inline double inner(response_space &a, response_space &b) {
  MADNESS_ASSERT(a.size() > 0);
  MADNESS_ASSERT(a.size() == b.size());
  MADNESS_ASSERT(!a[0].empty());
  MADNESS_ASSERT(a[0].size() == b[0].size());

  double value = 0.0;

  for (unsigned int i = 0; i < a.size(); i++) {
    // vmra.h function
    value += inner(a[i], b[i]);
  }

  return value;
}

auto transposeResponseMatrix(const response_matrix &x) -> response_matrix;
// Final piece for KAIN

} // End namespace madness
#endif // SRC_APPS_MOLRESPONSE_RESPONSE_FUNCTIONS_H_

// Deuces
