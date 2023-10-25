//
// Created by adrianhurtado on 2/11/22.
//

#ifndef MADNESS_RUNNERS_HPP
#define MADNESS_RUNNERS_HPP

#include <utility>

#include "ExcitedResponse.hpp"
#include "FrequencyResponse.hpp"
#include "ResponseExceptions.hpp"
#include "madness/chem/SCF.h"
#include "madness/tensor/tensor_json.hpp"
#include "madness/world/worldmem.h"
#include "response_data_base.hpp"
#include "response_functions.h"
#include "sstream"
#include "string"
#include "timer.h"
#include "write_test_input.h"
#include "x_space.h"

using path = std::filesystem::path;

auto split(const std::string &s, char delim) -> vector<std::string> {
    vector<std::string> result;
    std::stringstream ss(s);
    std::string item;

    while (getline(ss, item, delim)) { result.push_back(item); }

    return result;
}

auto addPath(const path &root, const std::string &branch) -> path {
    path p_branch = root;
    p_branch += branch;
    return p_branch;
}

struct runSchema {
    path root;               // root directory
    path json_database;      // json database
    path molecules;          // molecule directory
    path xc_path;            // create xc path
    path freq_json;          // path to freq_json
    path dalton_dipole_json; // path to dalton to dipole json
    path dalton_excited_json;// path to dalton excited json
    ResponseDataBase freq_json_data;

    explicit runSchema(World &world, const std::string &xc) {
        root = std::filesystem::current_path();//="/"+molecule_name;
        molecules = root / "molecules";
        xc_path = root / xc;
        json_database = root / "json_data";

        freq_json = json_database / "frequency.json";
        dalton_dipole_json = json_database / "dalton-dipole.json";
        dalton_excited_json = json_database / "dalton-excited.json";

        world.gop.fence();
        if (std::filesystem::exists(xc_path)) {
            if (world.rank() == 0) { ::print(xc, " Directory Exists"); }
        } else {
            if (world.rank() == 0) {
                ::print("Creating ", xc, " directory");
                std::filesystem::create_directory(xc_path);
            }
        }
        // Get the database where the calculation will be run from

        freq_json_data = ResponseDataBase();
        if (std::filesystem::exists(freq_json)) {
            std::ifstream ifs(freq_json);
            json j_read;
            ifs >> j_read;
            freq_json_data.set_data(j_read);
        }
        if (world.rank() == 0) { print(); }
    }

    void print() const {
        ::print("------------Database Runner---------------");
        ::print("Root: ", root);
        ::print("Molecule Directory: ", molecules);
        ::print("XC Path: ", xc_path);
        ::print("Freq Json Path: ", freq_json);
        ::print("Dalton Dipole Json Path: ", dalton_dipole_json);
        ::print("Dalton Excited Json Path: ", dalton_excited_json);
    }
};

struct moldftSchema {

    path moldft_path;

    path moldft_restart;
    path calc_info_json_path;
    json calc_info_json;
    path mol_path;
    std::string mol_name;
    std::string xc;

    moldftSchema(World &world, std::string molecule_name, std::string m_xc, const runSchema &schema)
        : mol_name(std::move(molecule_name)), xc(std::move(m_xc)) {
        moldft_path = addPath(schema.xc_path, '/' + mol_name);
        moldft_restart = addPath(moldft_path, "/moldft.restartdata.00000");
        calc_info_json_path = addPath(moldft_path, "/moldft.calc_info.json");
        mol_path = addPath(schema.molecules, "/" + mol_name + ".mol");
        if (std::filesystem::exists(moldft_restart) && std::filesystem::exists(calc_info_json_path)) {
            // if both exist, read the calc_info json
            std::ifstream ifs(calc_info_json_path);
            ifs >> calc_info_json;
            if (world.rank() == 0) {

                std::cout << "time: " << calc_info_json["time"] << std::endl;
                std::cout << "MOLDFT return energy: " << calc_info_json["return_energy"] << std::endl;
            }
        }
        if (world.rank() == 0) { print(); }
    }

    void print() const {
        ::print("----------------- Moldft Paths --------------------");
        ::print("moldft path :", moldft_path);
        ::print("moldft restart path :", moldft_restart);
        ::print("molecule path  path :", mol_path);
        ::print("calc_info json path :", calc_info_json_path);
    }
};

struct frequencySchema {

    const std::string mol_name;
    const std::string xc;
    const std::string op;

    const path moldft_path;
    vector<double> freq;

    frequencySchema(World &world, const runSchema &run_schema, const moldftSchema &m_schema, std::string r_operator,
                    bool static_calc);

    void print_schema() {
        print("Frequency Calculation");
        print("Molecule Name: ", mol_name);
        print("Functional: ", xc);
        print("Operator: ", op);
        print("MOLDFT PATH: ", moldft_path);
        print("Frequencies : ", freq);
    }
};
frequencySchema::frequencySchema(World &world, const runSchema &run_schema, const moldftSchema &m_schema,
                                 std::string r_operator, bool static_calc)
    : mol_name(m_schema.mol_name), xc(m_schema.xc), op(std::move(r_operator)), moldft_path(m_schema.moldft_path) {
    if (!static_calc) {
        if (world.rank() == 0) {
            print_schema();
            freq = run_schema.freq_json_data.get_frequencies(mol_name, xc, op);
        }
    } else {
        if (world.rank() == 0) { freq = {0.0}; }
    }

    world.gop.broadcast_serializable(freq, 0);
    return;
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
size_t set_excited_states(const ResponseDataBase &response_data_base, const std::string &molecule_name,
                          const std::string &xc) {

    const std::string property = "excited-state";

    try {
        return response_data_base.get_num_states(molecule_name, xc, property);
    } catch (json::exception &e) {
        std::cout << e.what() << std::endl;
        std::cout << "did not find the frequency data for [" << molecule_name << "][" << xc << "][" << property
                  << "]\n";
        return 4;
    }
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
std::filesystem::path generate_excited_run_path(const std::filesystem::path &moldft_path, const size_t &num_states) {
    std::string s_num_states = std::to_string(num_states);
    std::string run_name = "excited-" + s_num_states;
    // set r_params to restart true if restart file exist

    auto run_path = moldft_path;
    run_path += "/";
    run_path += std::filesystem::path(run_name);
    std::cout << run_path << endl;
    return run_path;
}
// sets the current path to the save path
/**
 * Generates the frequency save path with format
 * /excited_state/restart_[frequency_run_filename].00000
 *
 * @param excited_state restart path
 * @return
 */
std::pair<std::filesystem::path, std::string>
generate_excited_save_path(const std::filesystem::path &excited_run_path) {

    auto save_path = std::filesystem::path(excited_run_path);
    std::string save_string = "restart_excited";
    save_path += "/";
    save_path += save_string;

    save_path += ".00000";
    return {save_path, save_string};
}

struct excitedSchema {
    std::string xc;
    size_t num_states;
    path excited_state_run_path;
    path save_path;
    std::string save_string;

    path rb_json;


    excitedSchema(const runSchema &run_schema, const moldftSchema &m_schema) : xc(m_schema.xc) {
        num_states = set_excited_states(run_schema.freq_json_data, m_schema.mol_name, xc);
        excited_state_run_path = generate_excited_run_path(m_schema.moldft_path, num_states);
        auto [sp, s] = generate_excited_save_path(excited_state_run_path);
        save_path = sp;
        save_string = s;
        rb_json = addPath(excited_state_run_path, "/response_base.json");
    }

    void print() const {

        ::print("xc: ", xc);
        ::print("num states: ", num_states);
        ::print("excited_state run_path: ", excited_state_run_path);
        ::print("save_path: ", save_path);
        ::print("save_string: ", save_string);
    }
};

/**
 * Creates the xc directory in root directory of the
 *
 * Will create the xc directory if it does not already exist. Returns the path
 * of xc directory
 *
 *
 * @param root
 * @param xc
 * @return xc_path
 */
std::filesystem::path create_xc_path_and_directory(const std::filesystem::path &root, const std::string &xc) {

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
auto generate_frequency_save_path(const std::filesystem::path &frequency_run_path)
        -> std::pair<std::filesystem::path, std::string> {

    auto save_path = std::filesystem::path(frequency_run_path);
    auto run_name = frequency_run_path.filename();
    std::string save_string = "restart_" + run_name.string();
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
auto generate_response_frequency_run_path(const std::filesystem::path &moldft_path, const std::string &property,
                                          const double &frequency, const std::string &xc) -> std::filesystem::path {
    std::string s_frequency = std::to_string(frequency);
    auto sp = s_frequency.find('.');
    s_frequency = s_frequency.replace(sp, sp, "-");
    std::string run_name = property + "_" + xc + "_" + s_frequency;
    // set r_params to restart true if restart file exist

    auto run_path = moldft_path;
    run_path += "/";
    run_path += std::filesystem::path(run_name);
    return run_path;
}

/**
 * Runs moldft in the path provided.  Also generates the moldft input file_name
 * in the directory provided.
 *
 * @param world
 * @param moldft_path
 * @param moldft_filename
 * @param xc
 */
void runMOLDFT(World &world, const moldftSchema &moldftSchema, bool try_run, bool restart,
               const std::string &precision) {

    CalculationParameters param1;
    json calcInfo;

    if (world.rank() == 0) {
        param1.set_user_defined_value("maxiter", 30);
        //param1.set_user_defined_value("Kain", true);
        param1.set_user_defined_value<std::string>("xc", moldftSchema.xc);
        param1.set_user_defined_value<double>("l", 50);

        if (precision == "low") {
            param1.set_user_defined_value<vector<double>>("protocol", {1e-4, 1e-6, 1e-7});
            param1.set_user_defined_value<double>("dconv", 1e-4);
        } else if (precision == "high") {
            param1.set_user_defined_value<vector<double>>("protocol", {1e-4, 1e-6, 1e-8});
            param1.set_user_defined_value<double>("dconv", 1e-6);
        } else {
            param1.set_user_defined_value<vector<double>>("protocol", {1e-8});
            param1.set_user_defined_value<double>("dconv", 1e-6);
        }
        param1.set_user_defined_value<std::string>("localize", "new");
        CalculationParameters param_calc;
        if (std::filesystem::exists(moldftSchema.calc_info_json_path)) {
            std::cout << "Reading Calc Info JSON" << std::endl;
            std::ifstream ifs(moldftSchema.calc_info_json_path);
            ifs >> calcInfo;
            param_calc.from_json(calcInfo["parameters"]);
            print(param1.print_to_string());
            print(param_calc.print_to_string());
            print("param1 != param_calc = ", param1 != param_calc);
        }
    }
    world.gop.broadcast_serializable(param1, 0);

    //If the parameters are exactly equal do not run
    // If calc info doesn't exist the param_calc will definitely be different

    // if parameters are different or if I want to run no matter what run
    // if I want to restart and if I can. restart
    if (try_run) {
        if (world.rank() == 0) { print("-------------Running moldft------------"); }
        // if params are different run and if restart exists and if im asking to restar
        if (std::filesystem::exists(moldftSchema.moldft_restart) && restart) {
            param1.set_user_defined_value<bool>("restart", true);
        }
        world.gop.fence();
        if (world.rank() == 0) {
            molresponse::write_test_input test_input(param1, "moldft.in",
                                                     moldftSchema.mol_path);// molecule HF
        }
        world.gop.fence();
        commandlineparser parser;
        parser.set_keyval("input", "moldft.in");

        if (world.rank() == 0) print("input filename: ", parser.value("input"));


        print_meminfo(world.rank(), "startup");
        FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));

        std::cout.precision(6);
        SCF calc(world, parser);
        if (world.rank() == 0) {
            print("\n\n");
            print(" MADNESS Hartree-Fock and Density Functional Theory "
                  "Program");
            print(" ----------------------------------------------------------"
                  "\n");
            print("\n");
            calc.molecule.print();
            print("\n");
            calc.param.print("dft");
        }
        if (world.size() > 1) {
            calc.set_protocol<3>(world, 1e-4);
            calc.make_nuclear_potential(world);
            calc.initial_load_bal(world);
        }
        //vama
        calc.set_protocol<3>(world, calc.param.protocol()[0]);
        //calc.set_protocol<3>(world, 1e-4);
        world.gop.fence();
        MolecularEnergy ME(world, calc);
        world.gop.fence();
        // double energy=ME.value(calc.molecule.get_all_coords().flat()); // ugh!
        ME.value(calc.molecule.get_all_coords().flat());// ugh!
        world.gop.fence();
        const real_function_3d rho = 2.0 * calc.make_density(world, calc.aocc, calc.amo);
        auto dipole_t = calc.dipole(world, rho);
        std::map<std::string, double> results;
        results["scf_energy"] = calc.current_energy;
        world.gop.fence();
        if (world.rank() == 0) {
            calc.output_scf_info_schema(results, dipole_t);
            ME.output_calc_info_schema();
        }
    } else {
        if (world.rank() == 0) {
            print("Skipping Calculation and printing CALC INFO");
            std::cout << calcInfo;
        }
    }
}

/**
 * Sets the response parameters for a frequency response calculation and writes
 * to file
 *
 * @param r_params
 * @param property
 * @param xc
 * @param frequency
 */
void set_excited_parameters(World &world, ResponseParameters &r_params, const std::string &xc, const size_t &num_states,
                            const std::string &precision) {


    if (world.rank() == 0) {
        if (precision == "high") {
            r_params.set_user_defined_value<vector<double>>("protocol", {1e-4, 1e-6, 1e-8});
            r_params.set_user_defined_value<double>("dconv", 1e-6);
        } else if (precision == "low") {
            r_params.set_user_defined_value<vector<double>>("protocol", {1e-4, 1e-6});
            r_params.set_user_defined_value<double>("dconv", 1e-4);
        } else {
            r_params.set_user_defined_value<vector<double>>("protocol", {1e-9});
            r_params.set_user_defined_value<double>("dconv", 1e-6);
        }
        //r_params.set_user_defined_value("archive", std::string("../restartdata"));
        r_params.set_user_defined_value("maxiter", size_t(15));
        r_params.set_user_defined_value("maxsub", size_t(10));
        // if it's too large then bad guess is very strong
        r_params.set_user_defined_value("kain", true);
        r_params.set_user_defined_value("plot_all_orbitals", false);
        r_params.set_user_defined_value("save", true);
        r_params.set_user_defined_value("guess_xyz", true);
        r_params.set_user_defined_value("print_level", 1);
        // set xc, property, num_states,and restart
        r_params.set_user_defined_value("xc", xc);
        r_params.set_user_defined_value("excited_state", true);
        r_params.set_user_defined_value("states", num_states);
    }
    // Here
}

/**
 * Sets the response parameters for a frequency response calculation and writes
 * to file
 *
 * @param r_params
 * @param property
 * @param xc
 * @param frequency
 */
void setHyperpolarizabilityParameters(World &world, ResponseParameters &r_params, const std::string &property,
                                      const std::string &xc, const std::vector<double> &frequency,
                                      const std::string &precision) {
    if (world.rank() == 0) {
        r_params.set_user_defined_value("quadratic", true);
        r_params.set_user_defined_value("freq_range", frequency);
        r_params.set_user_defined_value("xc", xc);

        if (property == "dipole") {
            r_params.set_user_defined_value("dipole", true);
            r_params.set_derived_value<size_t>("states", 3);
        }

        // r_params.set_user_defined_value("property", property);

        if (precision == "high") {
            r_params.set_user_defined_value<vector<double>>("protocol", {1e-8});
            r_params.set_user_defined_value<double>("dconv", 1e-6);
        } else if (precision == "low") {
            r_params.set_user_defined_value<vector<double>>("protocol", {1e-6});
            r_params.set_user_defined_value<double>("dconv", 1e-4);
        }
    }
    world.gop.broadcast_serializable(r_params, 0);
}
/**
 * Sets the response parameters for a frequency response calculation and writes
 * to file
 *
 * @param r_params
 * @param property
 * @param xc
 * @param frequency
 */
void set_frequency_response_parameters(World &world, ResponseParameters &r_params, const std::string &property,
                                       const std::string &xc, const double &frequency, const std::string &precision) {
    if (world.rank() == 0) {
        if (precision == "high") {
            r_params.set_user_defined_value<vector<double>>("protocol", {1e-4, 1e-6, 1e-8});
            r_params.set_user_defined_value<double>("dconv", 1e-6);
        } else if (precision == "low") {
            r_params.set_user_defined_value<vector<double>>("protocol", {1e-4, 1e-6});
            r_params.set_user_defined_value<double>("dconv", 1e-4);
        } else {
            r_params.set_user_defined_value<vector<double>>("protocol", {1e-9});
            r_params.set_user_defined_value<double>("dconv", 1e-6);
        }
        //r_params.set_user_defined_value("archive", std::string("../restartdata"));
        r_params.set_user_defined_value("maxiter", size_t(35));
        r_params.set_user_defined_value("maxsub", size_t(3));
        r_params.set_user_defined_value("kain", true);
        r_params.set_user_defined_value("omega", frequency);
        r_params.set_user_defined_value("first_order", true);
        r_params.set_user_defined_value("plot_all_orbitals", true);
        r_params.set_user_defined_value("print_level", 1);
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
    world.gop.broadcast_serializable(r_params, 0);
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
static auto set_frequency_path_and_restart(World &world, ResponseParameters &parameters, const std::string &property,
                                           const double &frequency, const std::string &xc,
                                           const std::filesystem::path &moldft_path,
                                           std::filesystem::path &restart_path, bool restart) -> std::filesystem::path {

    if (world.rank() == 0) { print("restart path", restart_path); }


    // change the logic create save path
    auto frequency_run_path = generate_response_frequency_run_path(moldft_path, property, frequency, xc);
    world.gop.fence();
    if (world.rank() == 0) { print("frequency run path", frequency_run_path); }
    if (std::filesystem::is_directory(frequency_run_path)) {
        if (world.rank() == 0) { cout << "Response directory found " << std::endl; }
    } else {// create the file
        std::filesystem::create_directory(frequency_run_path);
        if (world.rank() == 0) { cout << "Creating response_path directory" << std::endl; }
    }
    std::filesystem::current_path(frequency_run_path);
    // frequency save path
    auto [save_path, save_string] = generate_frequency_save_path(frequency_run_path);
    if (world.rank() == 0) { print("save path", save_path); }
    if (world.rank() == 0) { print("save string", save_string); }


    if (world.rank() == 0) {
        parameters.set_user_defined_value("save", true);
        parameters.set_user_defined_value("save_file", save_string);
        if (restart) {//if we are trying a restart calculation
            if (std::filesystem::exists(save_path)) {
                //if the save path exists then we know we can
                // restart from the previous save
                parameters.set_user_defined_value("restart", true);
                parameters.set_user_defined_value("restart_file", save_string);
            } else if (std::filesystem::exists(restart_path)) {
                parameters.set_user_defined_value("restart", true);
                auto split_restart_path = split(restart_path.replace_extension("").string(), '/');
                std::string restart_file_short = "../" + split_restart_path[split_restart_path.size() - 2] + "/" +
                                                 split_restart_path[split_restart_path.size() - 1];
                parameters.set_user_defined_value("restart_file", restart_file_short);
                // Then we restart from the previous file instead
            } else {
                parameters.set_user_defined_value("restart", false);
            }
            // neither file exists therefore you need to start from fresh
        }
    }
    world.gop.broadcast_serializable(parameters, 0);
    // if restart and restartfile exists go ahead and set the restart file
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
auto RunResponse(World &world, const std::string &filename, double frequency, const std::string &property,
                 const std::string &xc, const std::filesystem::path &moldft_path, std::filesystem::path restart_path,
                 const std::string &precision) -> std::pair<std::filesystem::path, bool> {
    // Set the response parameters
    ResponseParameters r_params{};
    set_frequency_response_parameters(world, r_params, property, xc, frequency, precision);
    auto save_path =
            set_frequency_path_and_restart(world, r_params, property, frequency, xc, moldft_path, restart_path, true);
    if (world.rank() == 0) { molresponse::write_response_input(r_params, filename); }
    // if rbase exists and converged I just return save path and true
    if (std::filesystem::exists("response_base.json")) {
        std::ifstream ifs("response_base.json");
        json response_base;
        ifs >> response_base;
        if (response_base["converged"] && response_base["precision"]["dconv"] == r_params.dconv()) {
            return {save_path, true};
        }
    }
    auto calc_params = initialize_calc_params(world, std::string(filename));
    RHS_Generator rhs_generator;
    if (property == "dipole") {
        rhs_generator = dipole_generator;
    } else {
        rhs_generator = nuclear_generator;
    }
    FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));
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
    calc.solve(world);
    world.gop.fence();
    // set protocol to the first
    if (world.rank() == 0) {
        //calc.time_data.to_json(calc.j_molresponse);
        calc.output_json();
    }
    //calc.time_data.print_data();
    return {save_path, calc.j_molresponse["converged"]};
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
auto WriteVTKOutputs(World &world, const std::string &filename, double frequency, const std::string &property,
                     const std::string &xc, const std::filesystem::path &moldft_path,
                     std::filesystem::path restart_path, const std::string &precision)
        -> std::pair<std::filesystem::path, bool> {
    // Set the response parameters
    ResponseParameters r_params{};
    set_frequency_response_parameters(world, r_params, property, xc, frequency, precision);
    auto save_path =
            set_frequency_path_and_restart(world, r_params, property, frequency, xc, moldft_path, restart_path, true);
    if (world.rank() == 0) { molresponse::write_response_input(r_params, filename); }
    // if rbase exists and converged I just return save path and true
    auto calc_params = initialize_calc_params(world, std::string(filename));
    RHS_Generator rhs_generator;
    if (property == "dipole") {
        rhs_generator = dipole_generator;
    } else {
        rhs_generator = nuclear_generator;
    }
    FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));
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
    calc.write_vtk(world);
    world.gop.fence();
    // set protocol to the first
    //calc.time_data.print_data();
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
static void set_and_write_restart_excited_parameters(ResponseParameters &parameters, excitedSchema &schema,
                                                     bool restart) {

    parameters.set_user_defined_value("save", true);
    parameters.set_user_defined_value("save_file", schema.save_string);
    // if restart and restartfile exists go ahead and set the restart file
    if (restart && std::filesystem::exists(schema.save_path)) {
        print("setting restart");
        parameters.set_user_defined_value("restart", true);
        parameters.set_user_defined_value("restart_file", schema.save_string);
    } else {
        parameters.set_user_defined_value("restart", false);
    }
    std::string filename = "response.in";
    molresponse::write_response_input(parameters, filename);
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
static void create_excited_paths(excitedSchema &schema) {

    if (std::filesystem::is_directory(schema.excited_state_run_path)) {
        cout << "Response directory found " << std::endl;
    } else {// create the file
        std::filesystem::create_directory(schema.excited_state_run_path);
        cout << "Creating response_path directory" << std::endl;
    }
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
auto runExcited(World &world, excitedSchema schema, bool restart, const std::string &precision) -> bool {


    // Set the response parameters
    ResponseParameters r_params{};

    set_excited_parameters(world, r_params, schema.xc, schema.num_states, precision);
    create_excited_paths(schema);
    std::filesystem::current_path(schema.excited_state_run_path);
    set_and_write_restart_excited_parameters(r_params, schema, restart);

    auto calc_params = initialize_calc_params(world, "response.in");
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
void runFrequencyTests(World &world, const frequencySchema &schema, const std::string &high_prec) {
    std::filesystem::current_path(schema.moldft_path);
    // add a restart path
    auto restart_path =
            addPath(schema.moldft_path, "/" + schema.op + "_0-000000.00000/restart_" + schema.op + "_0-000000.00000");
    std::pair<std::filesystem::path, bool> success{schema.moldft_path, false};
    bool first = true;
    for (const auto &freq: schema.freq) {
        if (world.rank() == 0) { print(success.second); }
        std::filesystem::current_path(schema.moldft_path);
        if (first) {
            first = false;
        } else if (success.second) {
            // if the previous run succeeded then set the restart path
            restart_path = success.first;
            if (world.rank() == 0) { print("restart_path", restart_path); }
        } else {
            throw Response_Convergence_Error{};
        }
        success = RunResponse(world, "response.in", freq, schema.op, schema.xc, schema.moldft_path, restart_path,
                              high_prec);
        if (world.rank() == 0) { print("Frequency ", freq, " completed"); }
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
void write_VTK_outputs(World &world, const frequencySchema &schema, const std::string &high_prec) {
    std::filesystem::current_path(schema.moldft_path);
    // add a restart path
    auto restart_path =
            addPath(schema.moldft_path, "/" + schema.op + "_0-000000.00000/restart_" + schema.op + "_0-000000.00000");
    std::pair<std::filesystem::path, bool> success{schema.moldft_path, false};
    bool first = true;
    for (const auto &freq: schema.freq) {
        if (world.rank() == 0) { print(success.second); }
        std::filesystem::current_path(schema.moldft_path);
        if (first) {
            first = false;
        } else if (success.second) {
            // if the previous run succeeded then set the restart path
            restart_path = success.first;
            if (world.rank() == 0) { print("restart_path", restart_path); }
        } else {
            throw Response_Convergence_Error{};
        }
        success = WriteVTKOutputs(world, "response.in", freq, schema.op, schema.xc, schema.moldft_path, restart_path,
                                  high_prec);
        if (world.rank() == 0) { print("Frequency ", freq, " completed"); }
    }
}

// set up a function that creates a beta_json with the fields defined  below.  in each field there will
// a vector of values.

nlohmann::ordered_json create_beta_json() {
    // i need A B C to hold char values and A-freq, B-freq, C-freq to hold double values


    nlohmann::ordered_json beta_json = {{"A-freq", json::array()}, {"B-freq", json::array()}, {"C-freq", json::array()},
                                        {"A", json::array()},      {"B", json ::array()},     {"C", json::array()},
                                        {"Beta", json::array()}};
    return beta_json;
}


// for a set of frequencies create a table from the beta values
void append_to_beta_json(const std::array<double, 3> &freq, const std::array<double, 18> &beta,
                         nlohmann::ordered_json &beta_json) {


    // create 3 columns of directions for each A,B,C
    std::array<char, 18> direction_A{'X', 'X', 'X', 'X', 'X', 'X', 'Y', 'Y', 'Y',
                                     'Y', 'Y', 'Y', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z'};
    std::array<char, 18> direction_B{'X', 'X', 'X', 'Y', 'Y', 'Z', 'X', 'X', 'X',
                                     'Y', 'Y', 'Z', 'X', 'X', 'X', 'Y', 'Y', 'Z'};
    std::array<char, 18> direction_C{'X', 'Y', 'Z', 'Y', 'Z', 'Z', 'X', 'Y', 'Z',
                                     'Y', 'Z', 'Z', 'X', 'Y', 'Z', 'Y', 'Z', 'Z'};

    // append each value of the columns to the beta json
    // for each value of beta
    // capitalize the direction


    for (int i = 0; i < 18; i++) {
        beta_json["A-freq"].push_back(freq[0]);
        beta_json["B-freq"].push_back(freq[1]);
        beta_json["C-freq"].push_back(freq[2]);


        beta_json["A"].push_back(std::string(1, direction_A[i]));
        beta_json["B"].push_back(std::string(1, direction_B[i]));
        beta_json["C"].push_back(std::string(1, direction_C[i]));
        beta_json["Beta"].push_back(beta[i]);
    }
}
/**
 * Takes in a series of frequencies and runs a quadratic response calculations
 * for given property at given frequencies.
 *
 * The frequencies are given in a vector of doubles
 * For quadratic response the set of frequencies that a the first order response calculations have been run.
 *
 * The assumption is that response vectors are store in restart path.  If they aren't already then we need to
 * run the first order responnse calculation with the mad-freq executable.  This would be equivalent to running
 * the runFrequencyTests function.
 *
 * So the first step of this function is to check if the restart path exists for all frequencies.
 * If it does, then we can run the quadratic response calculation. which does a double for loop through the frequencies
 *
 * Also, this function will save the quadratic response output to a json file formatted
 *
 * A-freq, B-freq, C-freq, A, B, C, , Beta(value)
 *
 * where A-freq, B-freq, C-freq are the frequencies used in the first order response calculations
 * and A, B, C are the perturbation operator used in the quadratic response calculation
 * and Beta is the value of the quadratic response calculation
 *
 * As an example, if A,B,C are dipole operators then this calculations gives the hyperpolarizability
 *
 *
 *
 *
 * @param world
 * @param moldft_path
 * @param frequencies
 * @param xc
 * @param property
 */
void runQuadraticResponse(World &world, const frequencySchema &schema, const std::string precision) {
    std::filesystem::current_path(schema.moldft_path);

    bool run_first_order = false;

    // add a restart path
    auto restart_path =
            addPath(schema.moldft_path, "/" + schema.op + "_0-000000.00000/restart_" + schema.op + "_0-000000.00000");
    try {
        for (const auto &freq: schema.freq) {
            std::filesystem::current_path(schema.moldft_path);
            ResponseParameters r_params{};
            auto save_path = set_frequency_path_and_restart(world, r_params, schema.op, freq, schema.xc,
                                                            schema.moldft_path, restart_path, true);


            if (std::filesystem::exists("response_base.json")) {
                std::ifstream ifs("response_base.json");
                json response_base;
                ifs >> response_base;
                if (response_base["converged"] && response_base["precision"]["dconv"] == r_params.dconv()) {
                    { print("Response calculation already converged"); }
                    continue;
                } else {
                    if (world.rank() == 0) { print("Response calculation not converged"); }
                    run_first_order = true;
                    break;
                }
            }
            if (!std::filesystem::exists(save_path)) { throw Response_Convergence_Error{}; }
        }
        world.gop.fence();

        if (world.rank() == 0) { print("Running quadratic response calculations"); }
        std::filesystem::current_path(schema.moldft_path);
        RHS_Generator rhs_generator;
        if (schema.op == "dipole") {
            rhs_generator = dipole_generator;
        } else {
            rhs_generator = nuclear_generator;
        }
        ResponseParameters quad_parameters{};
        if (world.rank() == 0) { print("Set up rhs generator"); }

        setHyperpolarizabilityParameters(world, quad_parameters, schema.op, schema.xc, schema.freq, std::string());
        if (world.rank() == 0) { molresponse::write_response_input(quad_parameters, "quad.in"); }

        //auto calc_params = initialize_calc_params(world, std::string("quad.in"));
        commandlineparser parser;
        std::string moldft_archive = "moldft.restartdata";
        GroundStateCalculation ground_calculation{world, moldft_archive};
        if (world.rank() == 0) { ground_calculation.print_params(); }
        Molecule molecule = ground_calculation.molecule();
        quad_parameters.set_ground_state_calculation_data(ground_calculation);
        if (world.rank() == 0) { quad_parameters.print(); }

        world.gop.fence();
        FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));

        nlohmann::ordered_json beta_json = create_beta_json();

        QuadraticResponse quad_calculation{
                world,
                {ground_calculation, molecule, quad_parameters},
                rhs_generator,
        };
        //if beta.json exists remove it
        if (std::filesystem::exists("beta.json")) { std::filesystem::remove("beta.json"); }

        for (const auto &omega_b: schema.freq) {
            for (const auto &omega_c: schema.freq) {


                auto generate_omega_restart_path = [&](double frequency) {
                    auto linear_response_calc_path =
                            generate_response_frequency_run_path(schema.moldft_path, schema.op, frequency, schema.xc);
                    auto restart_file_and_path = generate_frequency_save_path(linear_response_calc_path);
                    return restart_file_and_path;
                };


                auto omega_a = omega_b + omega_c;
                if (world.rank() == 0) {
                    print("New combination of frequencies ", omega_a, " ", omega_b, " ", omega_c);
                    print("-------------------------------------------");
                }
                if (omega_a <= schema.freq.back()) {

                    auto restartA = generate_omega_restart_path(omega_a);
                    auto restartB = generate_omega_restart_path(omega_b);
                    auto restartC = generate_omega_restart_path(omega_c);
                    if (world.rank() == 0) {
                        print("Restart file for A", restartA.first);
                        print("Restart file for B", restartB.first);
                        print("Restart file for C", restartC.first);
                    }

                    // check if restartA exists
                    if (!std::filesystem::exists(restartA.first)) {
                        if (world.rank() == 0) { print("Restart file for omega_a = ", omega_a, " doesn't exist"); }
                        throw Response_Convergence_Error{};
                    } else {
                        if (world.rank() == 0) { print("Restart file for omega_a = ", omega_a, " exists"); }

                        std::array<double, 3> omegas{omega_a, omega_b, omega_c};

                        // remove .00000 from restartA.first

                        std::array<path, 3> restarts{restartA.first.replace_extension(""),
                                                     restartB.first.replace_extension(""),
                                                     restartC.first.replace_extension("")};

                        quad_calculation.set_x_data(world, omegas, restarts);
                        auto beta_abc = quad_calculation.compute_beta(world);
                        nlohmann::ordered_json beta_entry;
                        //beta_entry["omega_a"] = omega_a;
                        //beta_entry["omega_b"] = omega_b;
                        //beta_entry["omega_c"] = omega_c;


                        std::array<double, 18> beta_vector{};
                        std::copy(beta_abc.ptr(), beta_abc.ptr() + 3 * 6, beta_vector.begin());
                        append_to_beta_json({omega_a, omega_b, omega_c}, beta_vector, beta_json);

                        std::ofstream outfile("beta.json");
                        if (outfile.is_open()) {
                            outfile << beta_json.dump(4);
                            outfile.close();
                        } else {
                            std::cerr << "Failed to open file for writing: "
                                      << "beta.json" << std::endl;
                        }
                    }


                } else {
                    if (world.rank() == 0) {
                        print("Skipping omega_a = ", omega_a, " because it is greater than the max frequency");
                    }
                    continue;
                }
            }
        }

        // print the beta table
        if (world.rank() == 0) { print(beta_json.dump(4)); }

        // write the beta json to file
        std::ofstream ofs("beta.json");
        ofs << std::setw(4) << beta_json << std::endl;
        ofs.close();


    } catch (Response_Convergence_Error &e) {
        if (true) {
            // if the first order response calculations haven't been run then run them
            if (world.rank() == 0) { print("Running first order response calculations"); }

            runFrequencyTests(world, schema, precision);
        } else {
            if (world.rank() == 0) {
                print("First order response calculations haven't been run and can't be run");
                print("Quadratic response calculations can't be run");
            }
        }
    }
}

/**
 *
 * @param world
 * @param m_schema
 * @param try_moldft do we try moldft or not... if we try we still may restart
 * @param restart  do we force a restart or not
 * @param precision high precision or no?
 */
void moldft(World &world, moldftSchema &m_schema, bool try_moldft, bool restart, const std::string &precision) {

    if (std::filesystem::is_directory(m_schema.moldft_path)) {
        if (world.rank() == 0) { cout << "MOLDFT directory found " << m_schema.mol_path << "\n"; }
    } else {// create the file
        if (world.rank() == 0) {
            std::filesystem::create_directory(m_schema.moldft_path);
            cout << "Creating MOLDFT directory for " << m_schema.mol_name << ":/" << m_schema.moldft_path << ":\n";
        }
        world.gop.fence();
    }
    std::filesystem::current_path(m_schema.moldft_path);
    if (world.rank() == 0) { cout << "Entering : " << m_schema.moldft_path << " to run MOLDFT \n\n"; }
    runMOLDFT(world, m_schema, try_moldft, restart, precision);
}

#endif// MADNESS_RUNNERS_HPP
