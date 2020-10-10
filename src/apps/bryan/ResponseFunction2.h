/*
 *   Small class to hold response functions and to interact with KAIN solver.
 */

#ifndef MADNESS_APPS_TDHF_RESPONSEFUNC_INCLUDE
  #define MADNESS_APPS_TDHF_RESPONSEFUNC_INCLUDE

  #include <madness/mra/mra.h>
  #include <madness/mra/operator.h>

  #include <memory>
  #include <vector>

namespace madness {

class ResponseFunction {
  // Member variables
 public:
  unsigned int r_states;  // Num. of resp. states
  unsigned int g_states;  // Num. of ground states
  std::vector<std::vector<Function<double, 3>>> x;

  // Member functions
 public:
  // Default constructor
  ResponseFunction() : r_states(0), g_states(0) {}

  // Initializes functions to zero
  ResponseFunction(World& world, unsigned int m, unsigned int n)
      : r_states(m), g_states(n) {
    for (unsigned int i = 0; i < m; i++)
      x.push_back(zero_functions<double, 3>(world, n));
    x[0][0].world().gop.fence();
  }

  // Copy constructor
  ResponseFunction(const ResponseFunction& b)
      : r_states(b.r_states), g_states(b.g_states) {
    for (unsigned int i = 0; i < r_states; i++)
      x.push_back(madness::copy(b.x[0][0].world(), b.x[i]));
    x[0][0].world().gop.fence();
  }

  // Determines if two ResponseFunctions are the same size
  bool same_size(const ResponseFunction& b) const {
    return (r_states == b.r_states && g_states == b.g_states);
  }

  // 1D accessor for x
  std::vector<Function<double, 3>>& operator[](long i) { return x[i]; }

  // KAIN must have this
  ResponseFunction operator+(ResponseFunction& b) {
    MADNESS_ASSERT(same_size(b));

    ResponseFunction result(x[0][0].world(), r_states, g_states);

    for (unsigned int i = 0; i < r_states; i++) {
      result[i] = add(x[0][0].world(), x[i], b[i]);
    }

    result[0][0].world().gop.fence();
    return result;
  }

  ResponseFunction operator-(const ResponseFunction& b) const {
    MADNESS_ASSERT(same_size(b));

    ResponseFunction result(x[0][0].world(), r_states, g_states);

    for (unsigned int i = 0; i < r_states; i++) {
      result.x[i] = sub(x[0][0].world(), x[i], b.x[i]);
    }

    return result;
  }

  // KAIN must have this
  // Scaling by a constant
  ResponseFunction operator*(double a) const {
    ResponseFunction result(*this);

    for (unsigned int i = 0; i < r_states; i++) {
      madness::scale(x[0][0].world(), result.x[i], a, false);
    }

    x[0][0].world().gop.fence();
    return result;
  }

  // Scaling all internal functions by an external function
  // g[i][j] = x[i][j] * f
  ResponseFunction operator*(const Function<double, 3>& f) {
    ResponseFunction result;

    for (unsigned int i = 0; i < r_states; i++) {
      // Using vmra.h funciton
      result.push_back(mul(f.world(), f, this->x[i], false));
    }

    f.world().gop.fence();
    return result;
  }

  // KAIN must have this
  ResponseFunction operator+=(const ResponseFunction b) {
    MADNESS_ASSERT(same_size(b));

    for (unsigned int i = 0; i < r_states; i++) {
      x[i] = add(b.x[0][0].world(), x[i], b.x[i]);
    }

    return *this;
  }

  // Returns a deep copy
  ResponseFunction copy() {
    ResponseFunction result(x[0][0].world(), r_states, g_states);

    for (unsigned int i = 0; i < r_states; i++) {
      result.x[i] = madness::copy(x[0][0].world(), x[i], false);
    }
    x[0][0].world().gop.fence();

    return result;
  }

  // Mimicking std::vector with these 4
  void push_back(const std::vector<Function<double, 3>>& f) {
    x.push_back(f);
    r_states++;

    // Be smart with g_states
    if (g_states > 0) {
      MADNESS_ASSERT(g_states = f.size());
    } else {  // g_states == 0 (empty vector)
      g_states = f.size();
    }
  }
  void pop_back() {
    MADNESS_ASSERT(r_states >= 1);
    x.pop_back();
    r_states--;

    // Be smart with g_states
    if (r_states == 0) {  // removed last item
      g_states = 0;
    }
  }
  void clear() {
    x.clear();
    r_states = 0;
    g_states = 0;
  }
  unsigned int size() const { return r_states; }

  // Mimicing standard madness calls with these 3
  void zero() {
    for (unsigned int k = 0; k < r_states; k++) {
      x[k] = zero_functions<double, 3>(x[0][0].world(), g_states);
    }
  }

  void compress_rf() {
    for (unsigned int k = 0; k < r_states; k++) {
      compress(x[0][0].world(), x[k], true);
    }
  }

  void reconstruct_rf() {
    for (unsigned int k = 0; k < r_states; k++) {
      reconstruct(x[0][0].world(), x[k], true);
    }
  }

  void truncate_rf() {
    for (unsigned int k = 0; k < r_states; k++) {
      truncate(x[0][0].world(), x[k], true);
    }
  }

  // Returns norms of each state
  Tensor<double> norm2() {
    Tensor<double> answer(r_states);
    for (unsigned int i = 0; i < r_states; i++)
      answer(i) = sqrt(inner(x[i], x[i]));
    return answer;
  }

  // Scales each state (read: entire row) by corresponding vector element
  //     new[i] = old[i] * mat[i]
  void scale(Tensor<double>& mat) {
    for (unsigned int i = 0; i < r_states; i++)
      madness::scale(x[0][0].world(), x[i], mat[i], false);
    x[0][0].world().gop.fence();
  }
};

// Final piece for KAIN
inline double inner(ResponseFunction& a, ResponseFunction& b) {
  MADNESS_ASSERT(a.size() > 0);
  MADNESS_ASSERT(a.size() == b.size());
  MADNESS_ASSERT(a[0].size() > 0);
  MADNESS_ASSERT(a[0].size() == b[0].size());

  double value = 0.0;

  for (unsigned int i = 0; i < a.size(); i++) {
    // vmra.h function
    value += inner(a[i][0].world(), a[i], b[i]).sum();
  }

  return value;
}

}  // End namespace madness
#endif

// Deuces
