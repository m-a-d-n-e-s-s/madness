//
// Created by adrianhurtado on 3/25/22.
//

#ifndef MADNESS_RESPONSE_DATA_BASE_HPP
#define MADNESS_RESPONSE_DATA_BASE_HPP

#include "apps/external_headers/tensor_json.hpp"
#include <vector>
#include <filesystem>


class ResponseDataBase {
public:
    json j;

    json retrieve_data(const std::string &molecule, const std::string &xc, const std::string &property) const {
        return j.at(molecule).at(xc).at(property);
    }

    void output_data(const std::string &filename, const std::string &molecule, const std::string &xc,
                     const std::string &property) const {

        auto output_json = retrieve_data(molecule, xc, property);
        std::ofstream ofs{filename};
        ofs << output_json;
    }

    explicit ResponseDataBase(const json &j) : j(j) {}

    explicit ResponseDataBase() {}

    void set_data(const json &data) {
        j = json(data);
    }

    void print() {
        for (const auto &[key, value]: j.items()) {
            std::cout << key << " : " << value << "|n" << std::endl;
        }

    }


    std::vector<double>
    get_frequencies(const std::string &molecule, const std::string &xc, const std::string &property) const {
        return retrieve_data(molecule, xc, property).get<std::vector<double>>();
    }
};

json
generate_response_data(const std::filesystem::path &molecule_path, const std::string &xc, const std::string &property,
                       const vector<double> &freq) {
    json data;
    for (const std::filesystem::directory_entry &mol_path:
            std::filesystem::directory_iterator(molecule_path)) {
        if (mol_path.path().extension() == ".mol") {
            auto molecule_name = mol_path.path().stem();
            data[molecule_name][xc][property] = freq;
        }
    }
    std::cout << data << endl;
    return data;
}

#endif //MADNESS_RESPONSE_DATA_BASE_HPP
