//
// Created by adrianhurtado on 1/10/22.
//

#ifndef MADNESS_TENSOR_JSON_H
#define MADNESS_TENSOR_JSON_H

#include <algorithm>
#include <ostream>

// #include "catch.hpp"
#include <fstream>
#include <madness/external/nlohmann_json/json.hpp>

//#include "mra.h"
#include "tensor.h"

//using namespace madness;
//using json = nlohmann::json;

namespace madness {
    template<typename T>
    nlohmann::json tensor_to_json(const Tensor<T>& m) {
        nlohmann::json j = nlohmann::json{};
        // auto dimensions = m.dims();
        long size = m.size();  ///< Number of elements in the tensor
        long n_dims = m.ndim();///< Number of dimensions (-1=invalid; 0=no
        ///< supported; >0=tensor)
        auto dims = m.dims();// the size of each dimension
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

        //print(j);
        return j;
    }

    template<typename T>
    Tensor<T> tensor_from_json(const nlohmann::json& j) {
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

    template<typename T>
    using vec_pair_T = std::vector<std::pair<std::string, T>>;

    template<typename T>
    using vec_pair_tensor_T = std::vector<std::pair<std::string, Tensor<T>>>;

    template<typename T>
    void to_json(nlohmann::json& j, std::vector<std::pair<std::string, Tensor<T>>> tensors) {
        for (const auto& val: tensors) { j[val.first] = tensor_to_json(val.second); }
    }

    template<typename T>
    void to_json(nlohmann::json& j, const vec_pair_T<T>& num_vals) {
        for (const auto& val: num_vals) { j[val.first] = val.second; }
    }

    template<typename... Vecs>
    nlohmann::json add_time_tag(const nlohmann::json& j1) {

        auto print_time = std::chrono::system_clock::now();
        auto in_time_t = std::chrono::system_clock::to_time_t(print_time);
        std::stringstream ss;
        ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X");

        auto j=j1;
        j["time_tag"] ={};

        j["time_tag"]["time"] = ss.str();
        j["time_tag"]["wall_time"] = wall_time();
        j["time_tag"]["cpu_time"] =cpu_time();
        return j;
    }

    template<typename... Vecs>
    void output_schema(const std::string schema_name, const nlohmann::json& j) {
        auto j1= add_time_tag(j);
        std::ofstream ofs(schema_name + ".json");
        ofs << std::setw(4) << j1;
    }

    template<typename... Vecs>
    nlohmann::json input_schema(const std::string& schema_name) {
        // read old schema
        std::ifstream ifs(schema_name+".json");
        nlohmann::json j;
        if (ifs.is_open()) j=nlohmann::json::parse(ifs);
        return j;
    }

    template<typename... Vecs>
    void update_schema(const std::string schema_name, const nlohmann::json& jnew) {
        auto j=input_schema(schema_name);
        j.update(jnew);
        output_schema(schema_name,j);
    }
}// namespace madness

#endif// MADNESS_TENSOR_JSON_H
