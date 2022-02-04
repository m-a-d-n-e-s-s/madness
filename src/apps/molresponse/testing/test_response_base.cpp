//
// Created by adrianhurtado on 1/1/22.
//

#define CATCH_CONFIG_RUNNER
#include <testing/write_test_input.h>

#include <filesystem>

#include "apps/chem/SCF.h"
#include "apps/external_headers/catch.hpp"
#include "madness/world/worldmem.h"
#include "molresponse/ExcitedResponse.hpp"
#include "molresponse/ResponseBase.hpp"
#include "molresponse/ResponseExceptions.hpp"
#include "molresponse/global_functions.h"

int main(int argc, char *argv[]) {
    World &world = madness::initialize(argc, argv);
    int result = 0;
    world.gop.fence();
    startup(world, argc, argv);

    // Here we run all the tests
    result = Catch::Session().run(argc, argv);
    return result;

    // print_meminfo(world.rank(), "startup");
    // std::cout.precision(6);
    // print_stats(world);
}

void runMOLDFT(World &world, std::filesystem::path mol_path, std::string filename) {

    CalculationParameters param1;
    param1.set_user_defined_value("maxiter", 10);
    // write restart file
    // write restart file
    write_test_input test_input(param1, filename, mol_path);// molecule HF
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
}
void runExcitedState(World &world, std::string filename, int num_states,
                     std::filesystem::path runPath) {

    // Set the response parameters
    ResponseParameters r_params{};
    set_default_response_parameters(r_params);

    r_params.set_user_defined_value("xc", std::string("hf"));
    r_params.set_user_defined_value("states", size_t(num_states));
    r_params.set_user_defined_value("excited_state", true);
    r_params.set_user_defined_value("plot_all_orbitals", true);

    r_params.set_user_defined_value("save", true);
    r_params.set_user_defined_value("save_file", std::string("restart_excited"));
    // set r_params to restart true if restart file exist

    auto restartPath = runPath;
    restartPath += "/restart_excited.00000";
    std::cout << restartPath << endl;
    if (std::filesystem::exists(restartPath)) {
        r_params.set_user_defined_value("restart", true);
        r_params.set_user_defined_value("restart_file", std::string("restart_excited"));
    } else {
        std::cout << "Does not exist!!!" << endl;
    }

    write_response_input(r_params, filename);

    auto calc_params = initialize_calc_params(world, std::string(filename));
    auto &[ground_calculation, molecule, r_params_copy] = calc_params;
    vecfuncT ground_orbitals = ground_calculation.orbitals();

    print(norm2s_T(world, ground_orbitals));

    ExcitedResponse calc(world, calc_params);

    calc.solve(world);
    calc.output_json();
}

using json = nlohmann::json;

TEST_CASE("Run MOLDFT and create answers directory") {


    using namespace madness;
    World &world = World::get_default();
    std::cout.precision(6);

    auto root = std::filesystem::current_path();//="/"+molecule_name;
    // first step is to read the molecule directory for molecules... check if it exists else throw error

    auto molecule_path = root;
    molecule_path += "/molecules";

    try {
        if (std::filesystem::is_directory(molecule_path)) {
            for (const std::filesystem::directory_entry &mol_path:
                 std::filesystem::directory_iterator(molecule_path)) {

                std::filesystem::current_path(root);
                //for each molecule in molecules directory
                bool moldft_results_exists = false;
                json moldft_answers;

                if (mol_path.path().extension() == ".mol") {
                    auto molecule_name = mol_path.path().stem();
                    std::cout << "\n\n----------------------------------------------------\n";
                    std::cout << "Beginning Tests for Molecule: " << molecule_name << "\n";

                    // We would like to read the moldft results from corresponding json "molecule_name.json"
                    auto moldft_results = molecule_path;
                    moldft_results += "/";
                    moldft_results += molecule_name;
                    moldft_results += ".json";
                    std::cout << moldft_results << std::endl;

                    // if the results exists then save them into answers
                    if (std::filesystem::exists(moldft_results)) {
                        moldft_results_exists = true;
                        std::ifstream ifs(moldft_results.string());
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
                    //Read the corresponding json results if it exists
                    // TODO Read molecule_name.json file for the results we are testing against
                    auto moldft = root;
                    moldft += "/";
                    moldft += molecule_name;

                    if (std::filesystem::is_directory(moldft)) {

                        cout << "MOLDFT directory found " << molecule_name << "\n";

                    } else {// create the file
                        std::filesystem::create_directory(moldft);
                        cout << "Creating MOLDFT directory for " << molecule_name << ":/" << moldft
                             << ":\n";
                    }


                    std::filesystem::current_path(moldft);
                    cout << "Entering : " << moldft << " to run MOLDFT \n\n";

                    // We are going to look for both a json file and restartdata
                    auto json_path = moldft;
                    auto moldft_restart = moldft;
                    // Now I can check whether restart file exists and calc_info.json exists
                    moldft_restart += std::filesystem::path("/restartdata.00000");
                    json_path += std::filesystem::path("/calc_info.json");

                    // If both the restart file exists and the json file exists then check them against the previous result.
                    if (std::filesystem::exists(moldft_restart) &&
                        std::filesystem::exists(json_path)) {
                        // if both exist, read the calc_info json

                        json calc_info_json;
                        std::ifstream ifs(json_path);
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

                        runMOLDFT(world, mol_path, "moldft.in");
                        std::ifstream ifs(json_path);
                        nlohmann::json calc_info_json;
                        ifs >> calc_info_json;
                        // if moldft results does not exist then copy the results into answer
                        if (!moldft_results_exists) {

                            std::cout << " answers do not exist for comparison.. now saving to "
                                      << json_path << "!!!\n";
                            std::ofstream ofs(moldft_results.string());
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

                    SECTION("Excited Response") {

                        auto response_run_path = moldft;
                        response_run_path += std::filesystem::path("/excited_state");

                        if (std::filesystem::is_directory(response_run_path)) {
                            std::cout << "The Excited State Response Directory Exists" << std::endl;
                            std::cout << response_run_path << ":\n";

                        } else {// create the file
                            bool b = std::filesystem::create_directory(response_run_path);
                        }

                        // Now the current path is the excited state directory
                        std::filesystem::current_path(response_run_path);
                        auto response_filename = "response.in";
                        cout << "response file name:" << response_filename;
                        // make a set of molecule and num state pairs.
                        // so I run for a set of molecules each with different number of response
                        // states.
                        runExcitedState(world, response_filename, 4, response_run_path);
                        // now check if the answers exist.  if the answers do not exist run
                        // response else check the answers
                    }
                }
                // Now check if restart file exists and if calc_info.json exists
            }
        } else {
            std::cout << "did not find molecules" << std::endl;
        }
    } catch (const std::filesystem::filesystem_error &ex) { std::cerr << ex.what() << "\n"; }
}