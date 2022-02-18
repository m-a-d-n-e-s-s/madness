//
// Created by adrianhurtado on 2/11/22.
//
#define CATCH_CONFIG_RUNNER
#include "ExcitedResponse.hpp"
#include "ResponseExceptions.hpp"
#include "TDDFT.h"
#include "apps/chem/SCF.h"
#include "apps/external_headers/catch.hpp"
#include "apps/external_headers/tensor_json.hpp"
#include "madness/world/worldmem.h"
#include "response_functions.h"
#include "runners.hpp"
#include "string"
#include "timer.h"
#include "write_test_input.h"
#include "x_space.h"

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && defined(HAVE_UNISTD_H)

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

static inline int file_exists(const char *input_name) {
    struct stat buffer {};
    size_t rc = stat(input_name, &buffer);
    return (rc == 0);
}

#endif

int main(int argc, char *argv[]) {
    World &world = madness::initialize(argc, argv);
    int result = 0;
    world.gop.fence();
    startup(world, argc, argv);
    { result = Catch::Session().run(argc, argv); }

    return result;

    // print_meminfo(world.rank(), "startup");
    // std::cout.precision(6);
    // print_stats(world);
}

TEST_CASE("Run MOLDFT/RESPONSE") {


    using namespace madness;
    World &world = World::get_default();
    std::cout.precision(6);

    auto root = std::filesystem::current_path();//="/"+molecule_name;
    // first step is to read the molecule directory for molecules... check if it exists else throw error

    auto molecule_path = root;
    molecule_path += "/molecules";
    std::string xc = "hf";
    std::string property = "dipole";
    auto xc_path = root;
    xc_path += "/";
    xc_path += std::filesystem::path(xc);
    if (std::filesystem::is_directory(xc_path)) {

        cout << "XC directory found " << xc << "\n";

    } else {// create the file
        std::filesystem::create_directory(xc_path);
        cout << "Creating XC directory for " << xc << ":\n";
    }

    ResponseDataBase response_data_base = ResponseDataBase();
    try {
        if (std::filesystem::is_directory(molecule_path)) {
            for (const std::filesystem::directory_entry &mol_path:
                 std::filesystem::directory_iterator(molecule_path)) {
                vector<double> frequencies{0};

                std::filesystem::current_path(xc_path);
                //for each molecule in molecules directory
                bool moldft_results_exists = false;
                json moldft_answers;

                if (mol_path.path().extension() == ".mol") {
                    auto molecule_name = mol_path.path().stem();
                    std::cout << "\n\n----------------------------------------------------\n";
                    std::cout << "Beginning Tests for Molecule: " << molecule_name << "\n";
                    auto response_json_path =
                            generate_response_json_path(molecule_path, molecule_name, xc, property);
                    json response_json;
                    if (std::filesystem::exists(response_json_path)) {
                        std::cout << "response_json exists:" << std::endl;
                        try {
                            response_data_base.output_data(response_json_path.string(),
                                                           molecule_name, xc, property);
                            response_json =
                                    response_data_base.retrieve_data(molecule_name, xc, property);
                            frequencies = response_json.at("freq").get<std::vector<double>>();
                        } catch (json::exception &e) { std::cout << e.what() << std::endl; }
                    } else {
                        std::cout << "did not find the frequency data for [" << molecule_name
                                  << "][" << xc << "][" << property << "]\n";
                    }
                    std::cout << response_json << std::endl;

                    // We would like to read the moldft results from corresponding json "molecule_name.json"
                    auto moldft_json = generate_moldft_json(molecule_path,molecule_name);
                    std::cout << moldft_json << std::endl;

                    // if the results exists then save them into answers
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
                    //Read the corresponding json results if it exists
                    // TODO Read molecule_name.json file for the results we are testing against
                    auto moldft_path = generate_moldft_path(xc_path,molecule_name);

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
                    auto calcinfo_json_path = generate_moldft_calcinfojson_path(moldft_path);
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

                        runMOLDFT(world, mol_path, "moldft.in", xc);
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
                    // states.
                    try {
                        runFrequencyTests(world, moldft_path, frequencies, xc,property);
                    } catch (const SafeMPI::Exception &e) {
                        print(e);
                    } catch (const madness::MadnessException &e) {
                        std::cout << e << std::endl;
                    } catch (const madness::TensorException &e) {
                        print(e);
                    } catch (const char *s) { print(s); } catch (const std::string &s) {
                        print(s);
                    } catch (const std::exception &e) { print(e.what()); } catch (...) {
                        error("caught unhandled exception");
                    }
                    // now check if the answers exist.  if the answers do not exist run
                    // response else check the answers
                }
            }
            std::cout << "Please check what happens when I get to this point of the loop"
                      << std::endl;
            // Now check if restart file exists and if calc_info.json exists
        } else {
            std::cout << "did not find molecules" << std::endl;
        }
    } catch (const std::filesystem::filesystem_error &ex) { std::cerr << ex.what() << "\n"; }
}