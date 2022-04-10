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

    std::vector<double> get_frequencies(const std::string &molecule, const std::string &xc,
                                        const std::string &property) const {
        return retrieve_data(molecule, xc, property).get<std::vector<double>>();
    }

    void add_default_molecule(const json & response_keywords) {

        const std::string molecule_name = response_keywords["molecule"];
        const std::string xc = response_keywords["xc"];
        const std::string op = response_keywords["operator"];
        json j_add;
        if (op == "excited-state") {

            j_add[molecule_name][xc][op] = 4;

        } else if (op == "dipole") {

            if (std::filesystem::exists("molecules/dalton-excited.json")) {
                std::ifstream ifs("molecules/dalton-excited.json");
                try {

                    json dalton_excited;
                    ifs >> dalton_excited;
                    ::print("Read Dalton Excited");
                    ::print(dalton_excited);


                    std::vector<double> freq = dalton_excited[molecule_name][xc]["excited-state"]["aug-cc-pVTZ"]["response"]["freq"];
                    auto omega_max = freq.at(0);
                    omega_max = omega_max / 2.0;
                    ::print(omega_max);

                    std::vector<double> omegas = {0, omega_max / 8.0, omega_max / 4.0, omega_max / 2.0, omega_max};
                    j_add[molecule_name][xc][op] = omegas;


                }
                catch (const json::out_of_range &e) {
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

void addResponseKeyWord(json response_keywords) {
    // Adds response keyword to frequency.json
    // reads in frequency that json and merges

    const std::string molecule_name = response_keywords["molecule"];
    const std::string xc = response_keywords["xc"];
    const std::string op = response_keywords["operator"];

    ResponseDataBase data_base{};

    if (std::filesystem::exists("molecules/frequency.json")) {
        std::ifstream ifs("molecules/frequency.json");
        std::cout << "Trying to read frequency.json\n";
        json j_read;
        ifs >> j_read;
        std::cout << "READ IT\n";
        data_base.set_data(j_read);

            try {
                auto num_states = data_base.retrieve_data(molecule_name, xc, op);
                print(num_states);

            } catch (const json::out_of_range &e) {
                std::cout << e.what() << std::endl;
                if (std::filesystem::exists("molecules/" + molecule_name + ".mol")) {
                    // The molecule file exists in the database therefore it is okay to add to frequency.json
                    data_base.add_default_molecule(response_keywords);
                }

            } catch (const std::exception &e) { print(e.what()); }
            catch (...) {
                std::cout << "uncaught exception" << std::endl;
            }
    } else {
        if (std::filesystem::exists("molecules/" + molecule_name + ".mol")) {
            // The molecule file exists in the database therefore it is okay to add to frequency.json
            data_base.add_default_molecule(response_keywords);
        }

    }
}

#endif//MADNESS_RESPONSE_DATA_BASE_HPP
