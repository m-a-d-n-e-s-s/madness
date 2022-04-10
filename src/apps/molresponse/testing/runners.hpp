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
#include "response_data_base.hpp"

/**
 * Creates the xc directory in root directory of the
 *
 * Will create the xc directory if it does not already exist. Returns the path of xc directory
 *
 *
 * @param root
 * @param xc
 * @return xc_path
 */
std::filesystem::path create_xc_path_and_directory(
        const std::filesystem::path &root, const std::string &xc) {

    // copy construct the  root path
    auto xc_path = std::filesystem::path(root);
    xc_path += "/";
    xc_path += std::filesystem::path(xc);
    if (std::filesystem::is_directory(xc_path)) {

        cout << "XC directory found " << xc << "\n";

    } else {// create the file
        std::filesystem::create_directory(xc_path);
        cout << "Creating XC directory for " << xc << ":\n";
    }

    return xc_path;
}

// sets the current path to the save path
/**
 * Generates the frequency save path with format
 * /frequency_run_path/restart_[frequency_run_filename].00000
 *
 * @param frequency_run_path
 * @return
 */
std::pair<std::filesystem::path, std::string> generate_frequency_save_path(
        const std::filesystem::path &frequency_run_path) {

    auto save_path = std::filesystem::path(frequency_run_path);
    auto run_name = frequency_run_path.filename();
    std::string save_string = "restart_" + run_name.string();
    save_path += "/";
    save_path += save_string;

    save_path += ".00000";
    return {save_path, save_string};
}

// sets the current path to the save path
/**
 * Generates the frequency save path with format
 * /excited_state/restart_[frequency_run_filename].00000
 *
 * @param excited_state restart path
 * @return
 */
std::pair<std::filesystem::path, std::string> generate_excited_save_path(
        const std::filesystem::path &excited_run_path) {

    auto save_path = std::filesystem::path(excited_run_path);
    std::string save_string = "restart_excited";
    save_path += "/";
    save_path += save_string;

    save_path += ".00000";
    return {save_path, save_string};
}

/**
 * generates the frequency response path using the format
 * [property]_[xc]_[1-100]
 *
 * where 1-100 corresponds a frequency of 1.100
 *
 * @param moldft_path
 * @param property
 * @param frequency
 * @param xc
 * @return
 */
std::filesystem::path generate_response_frequency_run_path(const std::filesystem::path &moldft_path,
                                                           const std::string &property, const double &frequency,
                                                           const std::string &xc) {
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

/**
 * generates the frequency response path using the format
 * [property]_[xc]_[1-100]
 *
 * where 1-100 corresponds a frequency of 1.100
 *
 * @param moldft_path
 * @param property
 * @param frequency
 * @param xc
 * @return
 */
std::filesystem::path generate_excited_run_path(const std::filesystem::path &moldft_path,
                                                const size_t &num_states,
                                                const std::string &xc) {
    std::string s_num_states = std::to_string(num_states);
    std::string run_name = "excited-" + s_num_states;
    // set r_params to restart true if restart file exist

    auto run_path = moldft_path;
    run_path += "/";
    run_path += std::filesystem::path(run_name);
    std::cout << run_path << endl;
    return run_path;
}

/**
 * Generates the response json path with the following format
 * [molecule_name]_[xc]_[property].json
 * @param molecule_path
 * @param molecule_name
 * @param xc
 * @param property
 * @return
 */
std::filesystem::path generate_response_json_path(const std::filesystem::path &molecule_path,
                                                  const std::string &molecule_name, const std::string &xc,
                                                  const std::string &property) {
    std::string response_string = "/" + molecule_name + "_" + xc + "_" + property + ".json";
    auto r_p = std::filesystem::path(molecule_path);
    r_p += response_string;
    return r_p;
}

/**
 * Generates the moldft json path
 * @param molecule_path
 * @param molecule_name
 * @return
 */
std::filesystem::path generate_moldft_json(const std::filesystem::path &molecule_path,
                                           const std::string &molecule_name) {

    auto moldft_results = std::filesystem::path(molecule_path);
    moldft_results += "/";
    moldft_results += molecule_name;
    moldft_results += ".json";
    return moldft_results;
}

/**
 * Generates the moldft path given xc and molecule name
 * moldft_path=/[xc]/[molecule_name]
 * @param xc_path
 * @param molecule_name
 * @return
 */
std::filesystem::path generate_moldft_path(const std::filesystem::path &xc_path,
                                           const std::string &molecule_name) {
    auto moldft = std::filesystem::path(xc_path);;
    moldft += "/";
    moldft += molecule_name;
    return moldft;
}

/**
 * Generates the restart file name for moldft
 * @param moldft_path
 * @return
 */
std::filesystem::path generate_moldft_restart_path(const std::filesystem::path &moldft_path) {
    auto moldft_restart = std::filesystem::path(moldft_path);
    moldft_restart += std::filesystem::path("/restartdata.00000");
    return moldft_restart;
}

// generates the moldft calc info path generated by moldft
std::filesystem::path generate_moldft_calc_info_json_path(const std::filesystem::path &moldft_path) {
    auto json_path = std::filesystem::path(moldft_path);
    json_path += std::filesystem::path("/calc_info.json");
    return json_path;
}

/**
 * Runs moldft in the path provided.  Also generates the moldft input file_name in the directory provided.
 *
 * @param world
 * @param moldft_path
 * @param moldft_filename
 * @param xc
 */
void runMOLDFT(World &world, const std::string &mol_path,
               const std::string &xc) {

    CalculationParameters param1;
    param1.set_user_defined_value("maxiter", 10);
    param1.set_user_defined_value<std::string>("xc", xc);
    param1.set_user_defined_value<double>("l", 200);
    param1.set_user_defined_value<vector<double>>("protocol", {1e-4, 1e-6});
    // write restart file
    // write restart file
    write_test_input test_input(param1, "moldft.in", mol_path);// molecule HF
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

/**
 * Sets the default response properties
 * @param r_params
 */
void set_default_response_parameters(ResponseParameters &r_params) {

    r_params.set_user_defined_value("maxiter", size_t(10));
    r_params.set_user_defined_value("archive", std::string("../restartdata"));
    r_params.set_user_defined_value("kain", true);
    r_params.set_user_defined_value("maxsub", size_t(10));
    r_params.set_user_defined_value<vector<double>>("protocol", {1e-4, 1e-6});
}

/**
 * Sets the response parameters for an excited state restart calculation
 *
 * @param world
 * @param filename
 * @param num_states
 * @param xc
 */
void initialize_excited_restart(World &world, const std::string &filename, const size_t &num_states,
                                const std::string &xc) {

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

/**
 * Sets the response parameters for a frequency response calculation and writes to file
 *
 * @param r_params
 * @param property
 * @param xc
 * @param frequency
 */
void set_and_write_excited_parameters(ResponseParameters &r_params,
                                      const std::string &xc,
                                      const size_t &num_states) {
    r_params.set_user_defined_value<vector<double>>("protocol", {1e-4, 1e-6});
    r_params.set_user_defined_value("archive", std::string("../restartdata"));
    r_params.set_user_defined_value("maxiter", size_t(20));
    r_params.set_user_defined_value("maxsub", size_t(10));
    r_params.set_user_defined_value("kain", true);
    r_params.set_user_defined_value("first_order", true);
    r_params.set_user_defined_value("plot_all_orbitals", true);
    r_params.set_user_defined_value("save", true);
    r_params.set_user_defined_value("guess_xyz", false);
    // set xc, property, num_states,and restart
    r_params.set_user_defined_value("xc", xc);
    r_params.set_user_defined_value("excited_state", true);
    r_params.set_user_defined_value("states", num_states);
    // Here
}

/**
 * Sets the response parameters for a frequency response calculation and writes to file
 *
 * @param r_params
 * @param property
 * @param xc
 * @param frequency
 */
void set_and_write_frequency_response_parameters(ResponseParameters &r_params, const std::string &property,
                                                 const std::string &xc,
                                                 const double &frequency) {
    r_params.set_user_defined_value<vector<double>>("protocol", {1e-4, 1e-6});
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
        ResponseParameters &parameters, const std::string &property, const double &frequency, const std::string &xc,
        const std::filesystem::path &moldft_path, std::filesystem::path restart_path, bool restart) {

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


    std::filesystem::current_path(frequency_run_path);
    // Calling this function will make the current working directory the frequency save path
    auto[save_path, save_string] = generate_frequency_save_path(frequency_run_path);
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

/**
 *
 * @param world
 * @param filename
 * @param frequency
 * @param property
 * @param xc
 * @param moldft_path
 * @param restart_path
 * @return
 */
std::pair<std::filesystem::path, bool> RunResponse(World &world, const std::string &filename,
                                                   double frequency, const std::string &property,
                                                   const std::string &xc,
                                                   const std::filesystem::path &moldft_path,
                                                   std::filesystem::path restart_path) {

    // Set the response parameters
    ResponseParameters r_params{};
    set_and_write_frequency_response_parameters(r_params, property, xc, frequency);
    auto save_path = set_frequency_path_and_restart(r_params, property, frequency, xc, moldft_path,
                                                    restart_path, true);

    auto calc_params = initialize_calc_params(world, std::string(filename));
    RHS_Generator rhs_generator;
    if (property == "dipole") {
        rhs_generator = dipole_generator;
    } else {
        rhs_generator = nuclear_generator;
    }
    FrequencyResponse calc(world, calc_params, frequency, rhs_generator);
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

/***
 * sets the run path based on the run type set by r_params
 * creates the run directory and sets current directory to the run data
 * returns the name of parameter file to run from
 *
 * @param parameters
 * @param frequency
 * @param moldft_path
 */
static std::filesystem::path set_excited_path_and_restart(
        ResponseParameters &parameters, const size_t &num_states, const std::string &xc,
        const std::filesystem::path &moldft_path, std::filesystem::path restart_path, bool restart) {

    // change the logic create save path
    auto excited_state_path =
            generate_excited_run_path(moldft_path, num_states, xc);
    // Crea
    if (std::filesystem::is_directory(excited_state_path)) {
        cout << "Response directory found " << std::endl;
    } else {// create the file
        std::filesystem::create_directory(excited_state_path);
        cout << "Creating response_path directory" << std::endl;
    }


    std::filesystem::current_path(excited_state_path);
    // Calling this function will make the current working directory the frequency save path
    auto[save_path, save_string] = generate_excited_save_path(excited_state_path);
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
/**
 *
 * @param world
 * @param filename
 * @param frequency
 * @param property
 * @param xc
 * @param moldft_path
 * @param restart_path
 * @return
 */
bool RunExcited(World &world, const std::string &filename,
                size_t num_states,
                const std::string &xc,
                const std::filesystem::path &moldft_path,
                std::filesystem::path restart_path) {

    // Set the response parameters
    ResponseParameters r_params{};
    set_and_write_excited_parameters(r_params, xc, num_states);

    // restart is true if possible at the moment
    auto save_path = set_excited_path_and_restart(r_params, num_states, xc, moldft_path,
                                                  restart_path, true);

    auto calc_params = initialize_calc_params(world, std::string(filename));
    ExcitedResponse calc(world, calc_params);
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
    return true;
}

/**
 * Takes in the moldft path where moldft restart file exists
 * runs a response calculations for given property at given frequencies.
 *
 *
 * @param world
 * @param moldft_path
 * @param frequencies
 * @param xc
 * @param property
 */
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
            // if the previous run succeeded then set the restart path
            restart_path += success.first;
        }
        success = RunResponse(world, "response.in", freq, property, xc, moldft_path, restart_path);
    }
}

/**
 * Sets the frequency dependent on the data found in the response_data_base class
 * If the data is not found it will just return an empty vector
 * @param response_data_base
 * @param molecule_path
 * @param molecule_name
 * @param xc
 * @param property
 * @return
 */
std::vector<double>
set_frequencies(const ResponseDataBase &response_data_base, const std::filesystem::path &molecule_path,
                const std::string &molecule_name,
                const std::string &xc, const std::string &property) {

    auto response_json_path =
            generate_response_json_path(molecule_path, molecule_name, xc, property);


    if (std::filesystem::exists(response_json_path)) {
        std::cout << "response_json exists:" << std::endl;
        try {
            response_data_base.output_data(response_json_path.string(),
                                           molecule_name, xc, property);
            return response_data_base.get_frequencies(molecule_name, xc, property);
        } catch (json::exception &e) { std::cout << e.what() << std::endl; }
    } else {
        std::cout << "did not find the frequency data for [" << molecule_name
                  << "][" << xc << "][" << property << "]\n";
        return {0};
    }

}

/**
 * Sets the excited state data found in the response_data_base class
 * If the data is not found it will just return 4
 * @param response_data_base
 * @param molecule_path
 * @param molecule_name
 * @param xc
 * @param property
 * @return
 */
size_t
set_excited_states(const ResponseDataBase &response_data_base, const std::filesystem::path &molecule_path,
                   const std::string &molecule_name,
                   const std::string &xc) {

    const std::string property;
    auto response_json_path =
            generate_response_json_path(molecule_path, molecule_name, xc, property);


    if (std::filesystem::exists(response_json_path)) {
        std::cout << "response_json exists:" << std::endl;
        try {
            response_data_base.output_data(response_json_path.string(),
                                           molecule_name, xc, property);
            return response_data_base.get_num_states(molecule_name, xc, property);
        } catch (json::exception &e) { std::cout << e.what() << std::endl; }
    } else {
        std::cout << "did not find the frequency data for [" << molecule_name
                  << "][" << xc << "][" << property << "]\n";
        return 4;
    }
}

/**
 * Takes in the moldft path where moldft restart file exists
 * runs a response calculations for given property at given frequencies.
 *
 *
 * @param world
 * @param moldft_path
 * @param frequencies
 * @param xc
 * @param property
 */
void runExcitedStates(World &world, std::filesystem::path moldft_path,
                      size_t num_states, std::string xc) {

    std::filesystem::current_path(moldft_path);
    // add a restart path
    auto restart_path = moldft_path;
    restart_path += "/restart_excited.00000";
    std::filesystem::current_path(moldft_path);
    bool success = RunExcited(world, "response.in", num_states, xc, moldft_path, restart_path);
}

/**
 * Generates MOLDFT data in directory if not already generated and returns path to the directory.
 *
 * Tests whether the current data matches with moldft answers in the molecules directory
 * @param world
 * @param xc_path
 * @param xc
 * @param molecule_path
 * @param molecule_name
 * @return
 */
std::filesystem::path
run_moldft_path(World &world, const std::filesystem::path &xc_path, const std::string &xc,
                const std::filesystem::path &mol_path,
                const std::string &molecule_name) {

    auto moldft_json = generate_moldft_json(mol_path.stem(), molecule_name);
    std::cout << moldft_json << std::endl;
    bool moldft_results_exists = false;
    json moldft_answers;

    if (std::filesystem::exists(moldft_json)) {
        moldft_results_exists = true;
        std::ifstream ifs(moldft_json.string());
        //read results into json
        ifs >> moldft_answers;
        // Here are the current answers... check to see if th
        std::cout << "Here are the current answers for" << molecule_name
                  << " check to see if they need to be updated please!"
                  << std::endl;
        cout << moldft_answers;
    } else {
        std::cout << " We do not have moldft answers so please run and save the "
                     "results in the molecule directory"
                  << std::endl;
    }
    auto moldft_path = generate_moldft_path(xc_path, molecule_name);

    if (std::filesystem::is_directory(moldft_path)) {
        cout << "MOLDFT directory found " << molecule_name << "\n";
    } else {// create the file
        std::filesystem::create_directory(moldft_path);
        cout << "Creating MOLDFT directory for " << molecule_name << ":/" << moldft_path
             << ":\n";
    }
    std::filesystem::current_path(moldft_path);
    cout << "Entering : " << moldft_path << " to run MOLDFT \n\n";

    // We are going to look for both a json file and restartdata
    auto calcinfo_json_path = generate_moldft_calc_info_json_path(moldft_path);
    auto moldft_restart = generate_moldft_restart_path(moldft_path);
    // Now I can check whether restart file exists and calc_info.json exists

    // If both the restart file exists and the json file exists then check them against the previous result.
    if (std::filesystem::exists(moldft_restart) &&
        std::filesystem::exists(calcinfo_json_path)) {
        // if both exist, read the calc_info json
        json calc_info_json;
        std::ifstream ifs(calcinfo_json_path);
        ifs >> calc_info_json;

        std::cout << "time: " << calc_info_json["time"] << std::endl;
        std::cout << "MOLDFT return energy: " << calc_info_json["return_energy"]
                  << std::endl;
        std::cout << "MOLDFT return energy answer: "
                  << moldft_answers["return_energy"] << std::endl;
        CHECK(moldft_answers["return_energy"] == calc_info_json["return_energy"]);
        // need to check if they converged somehow
        //
    } else {
        std::cout << "restart file or calc_info.json does not exists for "
                  << molecule_name << " now running MOLDFT";

        runMOLDFT(world, mol_path.string(), xc);
        std::ifstream ifs(calcinfo_json_path);
        nlohmann::json calc_info_json;
        ifs >> calc_info_json;
        // if moldft results does not exist then copy the results into answer
        if (!moldft_results_exists) {

            std::cout << " answers do not exist for comparison.. now saving to "
                      << calc_info_json << "!!!\n";
            std::ofstream ofs(moldft_json.string());
            ofs << calc_info_json;

        } else {
            std::cout << "return energy is equal "
                      << (moldft_answers["return_energy"] ==
                          calc_info_json["return_energy"])
                      << "\n";
            CHECK(moldft_answers["return_energy"] ==
                  calc_info_json["return_energy"]);
        }
    }

    return moldft_path;

}


#endif//MADNESS_RUNNERS_HPP
