//
// Created by adrianhurtado on 3/25/22.
//

#ifndef MADNESS_RESPONSE_DATA_BASE_HPP
#define MADNESS_RESPONSE_DATA_BASE_HPP

#include <filesystem>
#include <utility>
#include <vector>

#include "madness/tensor/tensor_json.hpp"

using path = std::filesystem::path;

class ResponseDataBase {
public:
    json j;

    json retrieve_data(const std::string &molecule, const std::string &xc,
                       const std::string &property) const {
        ::print(j);
        return j.at(molecule).at(xc).at(property);
    }

    void output_data(const std::string &filename, const std::string &molecule,
                     const std::string &xc, const std::string &property) const {

        auto output_json = retrieve_data(molecule, xc, property);
        std::ofstream ofs{filename};
        ofs << output_json;
    }

    explicit ResponseDataBase(json j) : j(std::move(j)) {}

    explicit ResponseDataBase() = default;

    void set_data(const json &data) { j = json(data); }

    void print() {
        for (const auto &[key, value]: j.items()) {
            std::cout << key << " : " << value << "|n" << std::endl;
        }
    }


    size_t get_num_states(const std::string &molecule, const std::string &xc,
                          const std::string &property) const {
        return retrieve_data(molecule, xc, property).get<size_t>();
    }

    std::vector<double> get_frequencies(const std::string &molecule,
                                        const std::string &xc,
                                        const std::string &property) const {
        return retrieve_data(molecule, xc, property).get<std::vector<double>>();
    }

    void add_default_molecule(const json &response_keywords) {
        const std::string molecule_name = response_keywords["molecule"];
        const std::string xc = response_keywords["xc"];
        const std::string op = response_keywords["operator"];
        json j_add;
        if (op == "excited-state") {
            j_add[molecule_name][xc][op] = 8;
        } else if (op == "dipole") {
            if (std::filesystem::exists("molecules/dalton-excited.json")) {
                std::ifstream ifs("molecules/dalton-excited.json");
                try {
                    json dalton_excited;
                    ifs >> dalton_excited;
                    ::print("Read Dalton Excited");
                    ::print(dalton_excited);
                    std::vector<double> freq =
                            dalton_excited[molecule_name][xc]["excited-state"]
                                          ["aug-cc-pVTZ"]["response"]["freq"];
                    auto omega_max = freq.at(0);
                    omega_max = omega_max / 2.0;
                    ::print(omega_max);
                    std::vector<double> omegas = {0, omega_max / 8.0,
                                                  omega_max / 4.0,
                                                  omega_max / 2.0, omega_max};
                    j_add[molecule_name][xc][op] = omegas;
                } catch (const json::out_of_range &e) {
                    std::cout << e.what() << std::endl;
                    // The molecule file exists in the database therefore it is okay to add to frequency.json
                }
            } else {
                std::cout << " did not find dipole-excited.json" << std::endl;
                j_add[molecule_name][xc][op] = {0};
            }
        } else if (op == "nuclear") {
            j_add[molecule_name][xc][op] = {0};
        }
        j.merge_patch(j_add);
        std::ofstream ofs("molecules/frequency.json");
        ofs << std::setw(4) << j << std::endl;
        std::cout << "Added " << j_add << " frequency.json" << std::endl;
    }
};

auto generate_dipole_frequencies(const std::string &molecule_name,
                                 std::string xc) -> vector<double> {

    if (std::filesystem::exists("molecules/dalton-excited.json")) {
        std::ifstream ifs("molecules/dalton-excited.json");
        try {

            json dalton_excited;
            ifs >> dalton_excited;
            ::print("Read Dalton Excited for ", molecule_name);
            ::print(dalton_excited[molecule_name][xc]["excited-state"]
                                  ["aug-cc-pVTZ"]["response"]["freq"]);
            std::vector<double> freq =
                    dalton_excited[molecule_name][xc]["excited-state"]
                                  ["aug-cc-pVTZ"]["response"]["freq"];
            auto omega_max = freq.at(0);
            omega_max = omega_max / 2.0;
            ::print("max frequency at cc-pVTZ", omega_max);

            std::vector<double> omegas = {};
            int Nsteps = 9;
            for (int i = 0; i < Nsteps; i++) {
                omegas.push_back(omega_max * (double) i / 8.0);
            }
            return omegas;

        } catch (const json::out_of_range &e) {
            std::cout << e.what() << std::endl;
            return std::vector<double>(1, 0);
            // The molecule file exists in the database therefore it is okay to add to frequency.json
        } catch (const json::type_error &e) {
            std::cout << e.what() << std::endl;
            return std::vector<double>(1, 0);
        }

    } else {
        std::cout << " did not find dipole-excited.json" << std::endl;
        return {0};
    }
}
json generate_response_data(const std::filesystem::path &molecule_path,
                            const std::string &xc, const std::string &property,
                            const vector<double> &freq) {
    json data;
    for (const std::filesystem::directory_entry &mol_path:
         std::filesystem::directory_iterator(molecule_path)) {
        if (mol_path.path().extension() == ".mol") {
            auto molecule_name = mol_path.path().stem();
            data[molecule_name][xc][property] =
                    generate_dipole_frequencies(molecule_name, xc);
        }
    }
    //std::cout << data << endl;
    return data;
}

json generate_excited_data(const std::filesystem::path &molecule_path,
                           const std::string &xc, int num_states) {
    json data;
    for (const std::filesystem::directory_entry &mol_path:
         std::filesystem::directory_iterator(molecule_path)) {
        if (mol_path.path().extension() == ".mol") {
            auto molecule_name = mol_path.path().stem();
            const std::string property = "excited-state";
            data[molecule_name][xc][property] = num_states;
        }
    }
    return data;
};


#endif//MADNESS_RESPONSE_DATA_BASE_HPP
