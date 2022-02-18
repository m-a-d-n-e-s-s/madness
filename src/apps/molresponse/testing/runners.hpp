//
// Created by adrianhurtado on 2/11/22.
//


#ifndef MADNESS_RUNNERS_HPP
#define MADNESS_RUNNERS_HPP
#include "ExcitedResponse.hpp"
#include "FrequencyResponse.hpp"
#include "ResponseExceptions.hpp"
#include "TDDFT.h"
#include "apps/chem/SCF.h"
#include "apps/external_headers/catch.hpp"
#include "apps/external_headers/tensor_json.hpp"
#include "madness/world/worldmem.h"
#include "response_functions.h"
#include "string"
#include "timer.h"
#include "write_test_input.h"
#include "x_space.h"

// sets the current path to the save path
std::filesystem::path generate_frequency_save_path(std::filesystem::path frequency_run_path) {

    auto save_path = frequency_run_path;
    std::filesystem::current_path(save_path);
    auto run_name = frequency_run_path.filename();
    std::string save_string = "restart_" + run_name.string();
    save_path += "/";
    save_path += save_string;
    save_path += ".00000";
    return save_path;
}

std::filesystem::path generate_response_frequency_run_path(std::filesystem::path moldft_path,
                                                           std::string property, double frequency,
                                                           std::string xc) {
    std::string s_frequency = std::to_string(frequency);
    auto sp = s_frequency.find(".");
    s_frequency = s_frequency.replace(sp, sp, "-");
    std::string run_name = property + "_" + xc + "_" + s_frequency;
    // set r_params to restart true if restart file exist

    auto run_path = moldft_path;
    run_path += "/";
    run_path += std::filesystem::path(run_name);
    std::cout << run_path << endl;
    return run_path;
}
std::filesystem::path generate_response_json_path(std::filesystem::path molecule_path,
                                                  std::string molecule_name, std::string xc,
                                                  std::string property) {
    std::string response_string = "/" + molecule_name + "_" + xc + "_" + property + ".json";
    auto r_p = molecule_path;
    r_p += response_string;
    return r_p;
}
std::filesystem::path generate_moldft_json(std::filesystem::path molecule_path,
                                           std::string molecule_name) {

    auto moldft_results = molecule_path;
    moldft_results += "/";
    moldft_results += molecule_name;
    moldft_results += ".json";
    return moldft_results;
}
std::filesystem::path generate_moldft_path(std::filesystem::path xc_path,
                                           std::string molecule_name) {
    auto moldft = xc_path;
    moldft += "/";
    moldft += molecule_name;
    return moldft;
}
std::filesystem::path generate_moldft_restart_path(std::filesystem::path moldft_path) {
    auto moldft_restart = moldft_path;
    // Now I can check whether restart file exists and calc_info.json exists
    moldft_restart += std::filesystem::path("/restartdata.00000");
    return moldft_restart;
}
std::filesystem::path generate_moldft_calcinfojson_path(std::filesystem::path moldft_path) {
    auto json_path = moldft_path;
    // Now I can check whether restart file exists and calc_info.json exists
    json_path += std::filesystem::path("/calc_info.json");
    return json_path;
}
void runMOLDFT(World &world, std::filesystem::path moldft_path, std::string filename,
               std::string xc) {

    CalculationParameters param1;
    param1.set_user_defined_value("maxiter", 10);
    param1.set_user_defined_value<std::string>("xc", xc);
    param1.set_user_defined_value<double>("l", 200);
    param1.set_user_defined_value<vector<double>>("protocol", {1e-4, 1e-6, 1e-8});
    // write restart file
    // write restart file
    write_test_input test_input(param1, filename, moldft_path);// molecule HF
    commandlineparser parser;
    parser.set_keyval("input", test_input.filename());
    SCF calc(world, parser);
    calc.set_protocol<3>(world, 1e-4);
    MolecularEnergy ME(world, calc);
    // double energy=ME.value(calc.molecule.get_all_coords().flat()); // ugh!
    ME.value(calc.molecule.get_all_coords().flat());// ugh!
    ME.output_calc_info_schema();

    world.gop.fence();
}
void set_default_response_parameters(ResponseParameters &r_params) {

    r_params.set_user_defined_value("maxiter", size_t(10));
    r_params.set_user_defined_value("archive", std::string("../restartdata"));
    r_params.set_user_defined_value("kain", true);
    r_params.set_user_defined_value("maxsub", size_t(10));
    r_params.set_user_defined_value<vector<double>>("protocol", {1e-4, 1e-6, 1e-8});
}

// creates a response input for the test
void initialize_excited_restart(World &world, std::string filename, size_t num_states,
                                std::string xc) {

    // Set the response parameters
    ResponseParameters r_params{};
    set_default_response_parameters(r_params);

    r_params.set_user_defined_value("xc", xc);
    r_params.set_user_defined_value("states", num_states);
    r_params.set_user_defined_value("excited_state", true);


    r_params.set_user_defined_value("restart", true);
    r_params.set_user_defined_value("restart_file", std::string("restart_excited"));

    write_response_input(r_params, filename);
}
void set_and_write(ResponseParameters &r_params, std::string property, std::string xc,
                   double frequency) {
    r_params.set_user_defined_value<vector<double>>("protocol", {1e-4, 1e-6, 1e-8});
    r_params.set_user_defined_value("archive", std::string("../restartdata"));
    r_params.set_user_defined_value("maxiter", size_t(10));
    r_params.set_user_defined_value("maxsub", size_t(10));
    r_params.set_user_defined_value("kain", true);
    r_params.set_user_defined_value("omega", frequency);
    r_params.set_user_defined_value("first_order", true);
    r_params.set_user_defined_value("plot_all_orbitals", true);
    r_params.set_user_defined_value("save", true);
    // set xc, property, frequency,and restart
    r_params.set_user_defined_value("xc", xc);
    // Here
    if (property == "dipole") {
        r_params.set_user_defined_value("dipole", true);
    } else if (property == "nuclear") {
        r_params.set_user_defined_value("nuclear", true);
    }
}
/***
 * sets the run path based on the run type set by r_params
 * creates the run directory and sets current directory to the run data
 * returns the name of parameter file to run from
 *
 * @param parameters
 * @param frequency
 * @param moldft_path
 */
static std::filesystem::path set_frequency_path_and_restart(
        ResponseParameters &parameters, std::string property, double frequency, std::string xc,
        std::filesystem::path moldft_path, std::filesystem::path restart_path, bool restart) {

    // change the logic create save path
    auto frequency_run_path =
            generate_response_frequency_run_path(moldft_path, property, frequency, xc);
    // Crea
    if (std::filesystem::is_directory(frequency_run_path)) {
        cout << "Response directory found " << std::endl;
    } else {// create the file
        std::filesystem::create_directory(frequency_run_path);
        cout << "Creating response_path directory" << std::endl;
    }

    // Calling this function will make the current working directory the frequency save path
    auto save_path = generate_frequency_save_path(frequency_run_path);
    auto save_string = save_path.filename();
    print("save string", save_string);


    parameters.set_user_defined_value("save", true);
    parameters.set_user_defined_value("save_file", save_string);
    // if restart and restartfile exists go ahead and set the restart file
    if (restart) {
        if (std::filesystem::exists(save_path)) {

            parameters.set_user_defined_value("restart", true);
            parameters.set_user_defined_value("restart_file", save_string);
        } else if (std::filesystem::exists(restart_path)) {

            parameters.set_user_defined_value("restart", true);
            parameters.set_user_defined_value("restart_file", restart_path.string());

        } else {
            parameters.set_user_defined_value("restart", false);

            // neither file exists therefore you need to start from fresh
        }
    } else {
        parameters.set_user_defined_value("restart", false);
    }
    std::string filename = "response.in";
    write_response_input(parameters, filename);
    return save_path;
}
std::pair<std::filesystem::path, bool> RunResponse(World &world, std::string filename,
                                                   double frequency, std::string property,
                                                   std::string xc,
                                                   std::filesystem::path moldft_path,
                                                   std::filesystem::path restart_path) {

    // Set the response parameters
    ResponseParameters r_params{};
    set_and_write(r_params, property, xc, frequency);
    auto save_path = set_frequency_path_and_restart(r_params, property, frequency, xc, moldft_path,
                                                    restart_path, true);

    auto calc_params = initialize_calc_params(world, std::string(filename));
    FrequencyResponse calc(world, calc_params, frequency, dipole_generator);
    if (world.rank() == 0) {
        print("\n\n");
        print(" MADNESS Time-Dependent Density Functional Theory Response "
              "Program");
        print(" ----------------------------------------------------------\n");
        print("\n");
        calc_params.molecule.print();
        print("\n");
        calc_params.response_parameters.print("response");
        // put the response parameters in a j_molrespone json object
        calc_params.response_parameters.to_json(calc.j_molresponse);
    }
    // set protocol to the first
    calc.solve(world);
    calc.output_json();
    return {save_path, true};
}

void runFrequencyTests(World &world, std::filesystem::path moldft_path,
                       std::vector<double> frequencies, std::string xc, std::string property) {

    std::filesystem::current_path(moldft_path);
    // add a restart path
    auto restart_path = moldft_path;
    restart_path += "/restart_dipole_hf_0-000000.00000";
    std::pair<std::filesystem::path, bool> success{moldft_path, false};
    bool first = true;
    for (const auto &freq: frequencies) {
        std::filesystem::current_path(moldft_path);
        if (first) {
            first = false;
        } else if (success.second) {
            restart_path += success.first;
        }
        success = RunResponse(world, "response.in", freq, property, xc, moldft_path, restart_path);
    }
}
class ResponseDataBase {


    json j;

public:
    ResponseDataBase() {
        j["10_Be"]["hf"]["dipole"] =
                json{{"freq", vector<double>{0.0, 0.025, 0.050, 0.075, 0.100, 0.125}}};
        j["11_Ne"]["hf"]["dipole"] =
                json{{"freq", vector<double>{0.0, 0.100, 0.200, 0.300, 0.400, 0.500}}};
        j["12_Ar"]["hf"]["dipole"] =
                json{{"freq", vector<double>{0.0, 0.100, 0.200, 0.300, 0.400, 0.500}}};
    }

    json retrieve_data(std::string molecule, std::string xc, std::string property) const {
        return j.at(molecule).at(xc).at(property);
    }
    void output_data(std::string filename, std::string molecule, std::string xc,
                     std::string property) {

        std::ofstream ofs{filename};
        auto output_json = retrieve_data(molecule, xc, property);
        ofs << output_json;
    }
};


#endif//MADNESS_RUNNERS_HPP
