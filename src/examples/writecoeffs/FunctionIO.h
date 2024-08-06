#include "funcdefaults.h"
#include <iostream>
#include <madness/mra/mra.h>

using namespace madness;

constexpr int simple_pow(int a, int b) {
  if (b == 0) {
    return 1;
  } else {
    int result = 1;
    for (int i = 0; i < b; i++) {
      result *= a;
    }
    return result;
  }
}

template <typename T, std::size_t NDIM> class FunctionIO {

private:
  long k = FunctionDefaults<NDIM>::get_k();
  long ndims = NDIM;
  long npts_per_box = simple_pow(k, ndims);

public:
  static size_t count_leaf_nodes(const Function<T, NDIM> &f) {
    const auto &coeffs = f.get_impl()->get_coeffs();
    size_t count = 0;
    for (auto it = coeffs.begin(); it != coeffs.end(); ++it) {
      const auto &key = it->first;
      const auto &node = it->second;
      if (node.has_coeff()) {
        count++;
      }
    }
    f.get_impl()->world.gop.sum(count);
    return count;
  }
  static void write_function_coeffs(const Function<T, NDIM> &f,
                                    std::ostream &out, const Key<NDIM> &key) {
    const auto &coeffs = f.get_impl()->get_coeffs();
    auto it = coeffs.find(key).get();
    if (it == coeffs.end()) {
      for (int i = 0; i < key.level(); ++i)
        out << "  ";
      out << key << "  missing --> " << coeffs.owner(key) << "\n";
    } else {
      const auto &node = it->second;
      if (node.has_coeff()) {
        auto values = f.get_impl()->coeffs2values(key, node.coeff());
        for (int i = 0; i < key.level(); ++i)
          out << "  ";
        out << key.level() << " ";
        for (int i = 0; i < NDIM; ++i)
          out << key.translation()[i] << " ";
        out << std::endl;
        for (size_t i = 0; i < values.size(); i++)
          out << values.ptr()[i] << " ";
        out << std::endl;
      }
      if (node.has_children()) {
        for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
          write_function_coeffs(f, out, kit.key());
        }
      }
    }
  }
  static void write_function(const Function<T, NDIM> &f, std::ostream &out) {
    f.reconstruct();
    std::cout << "NUMBER OF LEAF NODES: " << count_leaf_nodes(f) << std::endl;

    auto flags = out.flags();
    auto precision = out.precision();
    out << std::setprecision(17);
    out << std::scientific;

    if (f.get_impl()->world.rank() == 0) {
      out << NDIM << std::endl;
      const auto &cell = FunctionDefaults<NDIM>::get_cell();
      for (int d = 0; d < NDIM; ++d) {
        for (int i = 0; i < 2; ++i)
          out << cell(d, i) << " ";
        out << std::endl;
      }
      out << f.k() << std::endl;
      out << count_leaf_nodes(f) << std::endl;

      write_function_coeffs(f, out, Key<NDIM>(0));
    }
    f.get_impl()->world.gop.fence();

    out << std::setprecision(precision);
    out.setf(flags);
  }

  static void read_function_coeffs(Function<T, NDIM> &f, std::istream &in,
                                   int num_leaf_nodes) {
    auto &coeffs = f.get_impl()->get_coeffs();

    for (int i = 0; i < num_leaf_nodes; i++) {
      Level n;
      Vector<Translation, NDIM> l;
      long dims[NDIM];
      in >> n;
      if (in.eof())
        break;

      for (int i = 0; i < NDIM; ++i) {
        in >> l[i];
        dims[i] = f.k();
      }
      Key<NDIM> key(n, l);

      Tensor<T> values(NDIM, dims);
      for (size_t i = 0; i < values.size(); i++)
        in >> values.ptr()[i];
      auto t = f.get_impl()->values2coeffs(key, values);

      // f.get_impl()->accumulate2(t, coeffs, key);
      coeffs.task(key, &FunctionNode<T, NDIM>::accumulate2, t, coeffs, key);
    }
  }

  static Function<T, NDIM> read_function(World &world, std::istream &in) {
    size_t ndim;
    in >> ndim;
    MADNESS_CHECK(ndim == NDIM);

    Tensor<double> cell(NDIM, 2);
    for (int d = 0; d < NDIM; ++d) {
      for (int i = 0; i < 2; ++i)
        in >> cell(d, i);
    }
    FunctionDefaults<NDIM>::set_cell(cell);

    int k;
    in >> k;
    int num_leaf_nodes;
    in >> num_leaf_nodes;
    FunctionFactory<T, NDIM> factory(world);
    Function<T, NDIM> f(factory.k(k).empty());
    world.gop.fence();

    read_function_coeffs(f, in, num_leaf_nodes);

    f.verify_tree();

    return f;
  }
};
