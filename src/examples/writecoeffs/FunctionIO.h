#include "ccpairfunction.h"
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
      // const auto &key = it->first;
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
        for (size_t i = 0; i < NDIM; ++i)
          out << key.translation()[i] << " ";
        out << std::endl;
#if HAVE_GENTENSOR
        MADNESS_EXCEPTION("FunctionIO not implemented for GenTensor", 0);
#else
        for (size_t i = 0; i < (size_t)values.size(); i++)
          out << values.ptr()[i] << " ";
#endif

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
      for (size_t d = 0; d < NDIM; ++d) {
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

      for (size_t i = 0; i < NDIM; ++i) {
        in >> l[i];
        dims[i] = f.k();
      }
      Key<NDIM> key(n, l);

      Tensor<T> values(NDIM, dims);
      for (size_t i = 0; i < (size_t)values.size(); i++)
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
    for (size_t d = 0; d < NDIM; ++d) {
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
template <typename T, std::size_t NDIM> struct FunctionIOData {

  long k = 0;
  long npts_per_box = 0;
  std::size_t ndim = NDIM;
  std::array<std::pair<double, double>, NDIM> cell;
  long num_leaf_nodes{};
  std::vector<std::array<long, NDIM + 1>> nl;
  std::vector<std::vector<double>> values;

  FunctionIOData() = default;

  explicit FunctionIOData(const Function<T, NDIM> &f) {

    npts_per_box = simple_pow(f.k(), NDIM);

    f.reconstruct();
    if (f.get_impl()->world.rank() == 0) {
      num_leaf_nodes = FunctionIO<T, NDIM>::count_leaf_nodes(f);
      ndim = NDIM;
      k = f.k();
      const auto &cell_world = FunctionDefaults<NDIM>::get_cell();
      for (size_t d = 0; d < NDIM; ++d) {
        cell[d].first = cell_world(d, 0);
        cell[d].second = cell_world(d, 1);
      }

      initialize_func_coeffs(f, Key<NDIM>(0));
    }
    f.get_impl()->world.gop.fence();
  }

  void initialize_func_coeffs(const Function<T, NDIM> &f,
                              const Key<NDIM> &key) {
    const auto &coeffs = f.get_impl()->get_coeffs();
    auto it = coeffs.find(key).get();
    if (it == coeffs.end()) {
      for (int i = 0; i < key.level(); ++i)
        std::cout << "  ";
      std::cout << key << "  missing --> " << coeffs.owner(key) << "\n";
    } else {
      const auto &node = it->second;
      if (node.has_coeff()) {
        auto node_values = f.get_impl()->coeffs2values(key, node.coeff());
        std::array<long, NDIM + 1> key_i;
        key_i[0] = key.level();
        for (size_t i = 0; i < NDIM; ++i)
          key_i[i + 1] = key.translation()[i];
        nl.push_back(key_i);
        std::vector<double> values_i(npts_per_box);
#if HAVE_GENTENSOR
        MADNESS_EXCEPTION("FunctionIO coeffs not implemented for GenTensor", 0);
#else
        std::copy(node_values.ptr(), node_values.ptr() + npts_per_box,
                  values_i.begin());
#endif
        values.push_back(values_i);
      }
      if (node.has_children()) {
        for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
          initialize_func_coeffs(f, kit.key());
        }
      }
    }
  }
  void set_function_coeffs(Function<T, NDIM> &f, int num_leaf_nodes) {
    auto &coeffs = f.get_impl()->get_coeffs();

    for (int i = 0; i < num_leaf_nodes; i++) {
      Vector<Translation, NDIM> l;
      long dims[NDIM];

      for (size_t i = 0; i < NDIM; ++i) {
        dims[i] = f.k();
      }

      auto n = nl[i][0];
      for (size_t j = 0; j < NDIM; ++j) {
        l[j] = nl[i][j + 1];
      }
      Key<NDIM> key(n, l);

      Tensor<T> values(NDIM, dims);
      std::copy(this->values[i].begin(), this->values[i].end(), values.ptr());
      auto t = f.get_impl()->values2coeffs(key, values);

      // f.get_impl()->accumulate2(t, coeffs, key);
      coeffs.task(key, &FunctionNode<T, NDIM>::accumulate2, t, coeffs, key);
    }
  }

  Function<T, NDIM> create_function(World &world) {

    size_t ndim = this->ndim;
    MADNESS_CHECK(ndim == NDIM);
    Tensor<double> cell_t(NDIM, 2);
    for (size_t d = 0; d < NDIM; ++d) {
      cell_t(d, 0) = cell[d].first;
      cell_t(d, 1) = cell[d].second;
    }

    FunctionDefaults<NDIM>::set_cell(cell_t);

    FunctionFactory<T, NDIM> factory(world);
    Function<T, NDIM> f(factory.k(k).empty());
    world.gop.fence();

    set_function_coeffs(f, num_leaf_nodes);

    f.verify_tree();

    return f;
  }
};

template <typename T, std::size_t NDIM>
void to_json(json &j, const FunctionIOData<T, NDIM> &p) {
  j = json{{"npts_per_box", p.npts_per_box},
           {"k", p.k},
           {"cell", p.cell},
           {"num_leaf_nodes", p.num_leaf_nodes},
           {"nl", p.nl},
           {"ndim", p.ndim},
           {"values", p.values}};
}

template <typename T, std::size_t NDIM>
void from_json(const json &j, FunctionIOData<T, NDIM> &p) {
  j.at("npts_per_box").get_to(p.npts_per_box);
  j.at("k").get_to(p.k);
  j.at("cell").get_to(p.cell);
  j.at("num_leaf_nodes").get_to(p.num_leaf_nodes);
  j.at("nl").get_to(p.nl);
  j.at("values").get_to(p.values);
  j.at("ndim").get_to(p.ndim);
}
