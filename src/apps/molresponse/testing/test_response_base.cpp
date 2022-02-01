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
    try {


        // Here we run all the tests
        result = Catch::Session().run(argc, argv);
        return result;


    } catch (const SafeMPI::Exception &e) {
        print(e);
        error("caught an MPI exception");
    } catch (const madness::MadnessException &e) {
        print(e);
        error("caught a MADNESS exception");
    } catch (const madness::TensorException &e) {
        print(e);
        error("caught a Tensor exception");
    } catch (const char *s) {
        print(s);
        error("caught a string exception");
    } catch (const std::string &s) {
        print(s);
        error("caught a string (class) exception");
    } catch (const std::exception &e) {
        print(e.what());
        error("caught an STL exception");
    } catch (...) { error("caught unhandled exception"); }

    // print_meminfo(world.rank(), "startup");
    // std::cout.precision(6);
    // print_stats(world);
}



void runMOLDFT(World &world, std::string molecule, std::string filename) {

    CalculationParameters param1;
    param1.set_user_defined_value("maxiter", 2);
    param1.set_user_defined_value("protocol", std::vector<double>({1.e-3}));
    // write restart file
    // write restart file
    write_test_input test_input(param1, filename, molecule);// molecule HF
    commandlineparser parser;
    parser.set_keyval("input", test_input.filename());
    SCF calc(world, parser);
    calc.set_protocol<3>(world, 1e-4);
    MolecularEnergy ME(world, calc);
    //double energy=ME.value(calc.molecule.get_all_coords().flat()); // ugh!
    ME.value(calc.molecule.get_all_coords().flat());// ugh!
    ME.output_calc_info_schema();

    world.gop.fence();
}

void set_default_response_parameters(ResponseParameters & r_params){

    r_params.set_user_defined_value("maxiter", size_t(10));
    r_params.set_user_defined_value("archive",std::string("../restartdata"));
    r_params.set_user_defined_value("kain",true);
    r_params.set_user_defined_value("maxsub",size_t(10));

}
void runExcitedState(World &world, std::string filename,int num_states,std::filesystem::path path) {

    // Set the response parameters
    ResponseParameters r_params{};
    set_default_response_parameters(r_params);
    r_params.set_user_defined_value("xc",std::string("hf"));
    r_params.set_user_defined_value("states",size_t(num_states));
    r_params.set_user_defined_value("excited_state",true);

    r_params.set_user_defined_value("save",true);
    r_params.set_user_defined_value("save_file",std::string("restart_exited"));
    // set r_params to restart true if restart file exist

    write_response_input(r_params,filename);

    auto calc_params = initialize_calc_params(world, std::string(filename));
    auto &[ground_calculation, molecule, r_params_copy] = calc_params;
    vecfuncT ground_orbitals = ground_calculation.orbitals();

    print(norm2s_T(world,ground_orbitals));

    ExcitedResponse calc(world, calc_params);

    calc.solve(world);
    calc.output_json();


}

TEST_CASE("Creating new molecule directory") {

    using namespace madness;
    World &world = World::get_default();
    std::cout.precision(6);

    std::string molecule_name = "hf";


    {
        auto root = std::filesystem::current_path();//="/"+molecule_name;
        std::cout << "root:   " << root << endl;
        auto f_mol = std::filesystem::path("/" + molecule_name);

        auto working_path = root;
        working_path += f_mol;


        if (std::filesystem::is_directory(working_path)) {

            cout << working_path << ":\n";


        } else {// create the file
            bool b = std::filesystem::create_directory(working_path);
        }

        auto restart_path = std::filesystem::path("/restartdata.00000");

        // working path is now the molecule directory
        std::filesystem::current_path(working_path);
        auto molecule_path = std::filesystem::current_path();
        auto moldft_restart = std::filesystem::current_path();
        moldft_restart += restart_path;


        std::string response_type = "excited_state";


        working_path += restart_path;

        // now check if a restart file exist for that molecule
        // Check if the restartdata exists check if...if it doesn't then run moldft with settings
        if (std::filesystem::exists(moldft_restart)) {
            // need to check if they converged somehow


            //
        } else {
            std::cout << "restart file does not exists for " << f_mol << std::endl;
            //run moldft to generate restartdata file
            runMOLDFT(world, molecule_name, "moldft.in");
        }
        // now we know that restart data exists


        auto response_type_path = std::filesystem::path("/excited_state");
        auto response_run_path = std::filesystem::current_path();
        response_run_path += response_type_path;

        if (std::filesystem::is_directory(response_run_path)) {
            cout << response_run_path << ":\n";
        } else {// create the file
            bool b = std::filesystem::create_directory(response_run_path);
        }
        std::filesystem::current_path(response_run_path);
        auto response_filename = response_run_path.string() + "/response.in";
        cout << "response file name:" << response_filename;

        // make a set of molecule and num state pairs.
        // so I run for a set of molecules each with different number of response states.
        runExcitedState(world,response_filename,4);


        // now check if the answers exist.  if the answers do not exist run response else check the answers
    }

    /*

    CalculationParameters param1;
    param1.set_user_defined_value("maxiter", 2);
    param1.set_user_defined_value("protocol", std::vector<double>({1.e-3}));
    // write restart file
    // write restart file
    write_test_input test_input(param1, "hf");// molecule HF
    commandlineparser parser;
    parser.set_keyval("input", test_input.filename());
    SCF calc(world, parser);
    calc.set_protocol<3>(world, 1e-4);
    MolecularEnergy ME(world, calc);
    //double energy=ME.value(calc.molecule.get_all_coords().flat()); // ugh!
    ME.value(calc.molecule.get_all_coords().flat());// ugh!
    ME.output_calc_info_schema();

    world.gop.fence();

    // Set the response parameters
    ResponseParameters r_params{};
    r_params.set_user_defined_value("maxiter", 10);
    r_params.set_user_defined_value("xc", 10);


    /*
    auto calc_params = initialize_calc_params(world, std::string(input_file));
    auto &[ground_calculation, molecule, r_params] = calc_params;
    vecfuncT ground_orbitals = ground_calculation.orbitals();
    print(norm2s_T(world,ground_orbitals));

    r_params.print("ResponseParameters", "Not sure footer");

    ExcitedResponse calc(world, calc_params);

    calc.solve(world);
    calc.output_json();
     */
}
