//
// Created by adrianhurtado on 1/10/22.
//

#ifndef MADNESS_TENSOR_JSON_H
#define MADNESS_TENSOR_JSON_H

#include <algorithm>
#include <ostream>

#include "catch.hpp"
#include "json.hpp"
#include "mra.h"
#include "tensor.h"
#include "timers.h"

using namespace madness;
using json = nlohmann::json;

template <typename T>
json tensor_to_json(const Tensor<T> &m) {
  json j = json{};
  // auto dimensions = m.dims();
  long size = m.size();    ///< Number of elements in the tensor
  long n_dims = m.ndim();  ///< Number of dimensions (-1=invalid; 0=no
                           ///< supported; >0=tensor)
  auto dims = m.dims();    // the size of each dimension
  // long id = m.id();       ///< Id from TensorTypeData<T> in type_data.h
  // auto strides = m.strides();
  auto fm = m.flat();
  auto m_vals_vector = std::vector<T>(size);
  auto m_dims_vector = std::vector<long>(n_dims);
  std::copy(&fm[0], &fm[0] + size, m_vals_vector.begin());
  std::copy(dims, dims + n_dims, m_dims_vector.begin());

  // This is everything we need to translate to a numpy vector...
  j["size"] = size;
  j["vals"] = m_vals_vector;
  j["dims"] = m_dims_vector;

  print(j);
  return j;
}

template <typename T>
Tensor<T> tensor_from_json(const json &j) {
  // need to be explicit here about types so we find the proper Tensor
  // constructors
  long size = j["size"];
  std::vector<T> m_vals_vector = j["vals"];
  std::vector<long> m_dims_vector = j["dims"];

  Tensor<T> flat_m(size);
  // copy the values from the vector to the flat tensor
  std::copy(m_vals_vector.begin(), m_vals_vector.end(), &flat_m[0]);
  // reshape the tensor using dimension vector
  Tensor<T> m = flat_m.reshape(m_dims_vector);
  return m;
}

using vec_pair_ints = std::vector<std::pair<std::string, int>>;

template <typename T>
using vec_pair_T = std::vector<std::pair<std::string, T>>;

template <typename T>
using vec_pair_tensor_T = std::vector<std::pair<std::string, Tensor<T>>>;

template <typename T>
void to_json(json &j, std::vector<std::pair<std::string, Tensor<T>>> tensors) {
  for (const auto &val : tensors) {
    j[val.first] = tensor_to_json(val.second);
  }
}

template <typename T>
void to_json(json &j, const vec_pair_T<T> &num_vals) {
  for (const auto &val : num_vals) {
    j[val.first] = val.second;
  }
}
template <typename... Vecs>
void output_schema(std::string schema_name, const json &j) {
  std::ofstream ofs(schema_name + ".json");
  ofs << j;
}

#endif  // MADNESS_TENSOR_JSON_H
