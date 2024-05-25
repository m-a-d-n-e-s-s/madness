#ifndef MADNESS_RUNNERS_HPP
#define MADNESS_RUNNERS_HPP

#include <algorithm>
#include <filesystem>
#include <utility>

#include "CalculationParameters.h"
#include "maddft/write_test_input.h"
#include "madness/chem/SCF.h"
#include "madness/world/world.h"
#include "madness/world/worldmem.h"
#include "molresponse/FrequencyResponse.hpp"
#include "molresponse/ResponseExceptions.hpp"
#include "molresponse/response_parameters.h"
#include "sstream"
#include "string"

using path = std::filesystem::path;

auto split(const std::string &s, char delim) -> vector<std::string>
{
    vector<std::string> result;
    std::stringstream ss(s);
    std::string item;

    while (getline(ss, item, delim))
    {
        result.push_back(item);
    }

    return result;
}

auto addPath(const path &root, const std::string &branch) -> path
{
    path p_branch = root;
    p_branch += branch;
    return p_branch;
}

void write_molecule_json_to_input_file(const json &molecule_json, std::ostream &output_stream)
{
    output_stream << "geometry" << std::endl;
    // get parameters seperately
    auto parameters = molecule_json["parameters"];

    for (auto &[key, value] : parameters.items())
    {
        output_stream << "    " << key << " " << value << std::endl;
    }

    // grab the geometry and symbols as two separate arrays
    std::vector<std::string> symbols;
    std::vector<std::vector<double>> geometry;

    auto coords = molecule_json["geometry"];
    auto symbols_json = molecule_json["symbols"];

    auto n = symbols_json.size();
    std::vector<std::pair<std::string, std::vector<double>>> lines(n);

    for (int i = 0; i < n; i++)
    {
        symbols.push_back(symbols_json[i]);
        geometry.push_back(coords[i]);
    }

    for (int i = 0; i < n; i++)
    {
        output_stream << "    " << symbols[i] << " " << geometry[i][0] << " " << geometry[i][1] << " " << geometry[i][2] << std::endl;
    }
    output_stream << "end" << std::endl << std::endl;
}

void write_json_to_input_file(const json &input_json, const std::vector<std::string> &keys, std::ostream &output_stream)
{
    for (const auto &key : keys)
    {
        if (input_json.find(key) == input_json.end())
        {
            throw std::runtime_error("Key not found in input json");
        }

        output_stream << key << std::endl;
        // for each key within the block write the key value pair to the file line
        // by line
        for (auto &[key, value] : input_json[key].items())
        {
            output_stream << "    " << key << " " << value << std::endl;
        }
        output_stream << "end" << std::endl << std::endl;
    }
}

void write_moldft_input(const json &input_json, std::ostream &out)
{
    write_json_to_input_file(input_json, {"dft"}, out);
    write_molecule_json_to_input_file(input_json["molecule"], out);
}

void write_response_input(const json &input_json, std::ostream &out) { write_json_to_input_file(input_json, {"response"}, out); }

// Three ways to construct the ParameterManagerClass.
// The first is a simple input file which contains
// moldft and response parameters and geometry
//
// The second is a json file which contains the same information
//
// The last is a mix of input file and molecule file
// which the input file contains the moldft and response
// parameters and the molecule file contains the geometry
// and the input file can be either a json or a text file

class ParameterManager
{
private:
    path input_file_path;
    path input_file_json_path;

    path final_input_file_path = "final_input";
    path final_input_json_path = "final_input.json";

    json all_input_json;
    commandlineparser parser;

    Molecule molecule;
    CalculationParameters moldft_params{};
    ResponseParameters molresponse_params{};

    bool run_moldft = false;
    bool run_response = false;

public:
    ParameterManager() = default;

    void print_file_paths() const
    {
        ::print("------------Parameter Manager---------------");
        ::print("Input File Path: ", input_file_path);
        ::print("Final Input File Path: ", final_input_file_path);
        ::print("Input File Json Path: ", input_file_json_path);
        ::print("-------------------------------------------");
    }

    static void help()
    {
        print_header2("help page for MADNESS DFT and Response Properties Code ");
        print("This code is designed to run DFT and Response Property calculations");
        print("Within the input one defines both the ground and response "
              "calculations in the input file by specifiying the dft and response "
              "blocks");
        print("By defining the quadratic block one can compute quadratic response "
              "properties such as the hyperpolarizability");
    }
    void print_params() const
    {
        ::print("------------Parameter Manager---------------");
        ::print("Molecule: ");
        molecule.print();
        ::print("Moldft Parameters: ");
        moldft_params.print();
        ::print("Molresponse Parameters: ");
        molresponse_params.print();
        ::print("-------------------------------------------");
    }

    /** Reads the molecule from the json file and orients it
     * if paramater is set to true.  The molecule is then copied
     * back to the json file to ensure that the geometry is
     * consistent
     *
     *
     * @param j
     */
    void read_molecule_from_json_and_orient(json &j)
    {
        if (j.contains("molecule"))
        {
            molecule.from_json(j["molecule"]);
            j["molecule"] = molecule.to_json();
        }
        else
        {
            throw std::runtime_error("Molecule not found in input file");
        }
    }

    /**
     * Reads the chemical parameters from the json file
     * Available parameters are for
     *  - dft
     *  - response
     *
     *
     * @param j
     */
    void read_chem_params_from_json(const json &j)
    {
        if (j.contains("dft"))
        {
            all_input_json["dft"] = j["dft"];
            moldft_params.from_json(j["dft"]);
        }
        if (j.contains("response"))
        {
            all_input_json["response"] = j["response"];
            molresponse_params.from_json(j["response"]);
        }
    }

    /**
     * Reads the chemical parameters from standard madness input file
     * and writes the input to json if parameters are defined
     * Available parameters are for
     *  - dft
     *  - response
     *
     *
     * @param world : madness world
     * @param input_file : path to the input file
     */
    void read_chem_params_and_write_json(World &world, const path &input_file)
    {
        parser.set_keyval("input", input_file.string());
        moldft_params = CalculationParameters(world, this->parser);
        auto moldft_json = moldft_params.to_json_if_precedence("defined");
        if (!moldft_json.is_null())
        {
            all_input_json["dft"] = moldft_params.to_json_if_precedence("defined");
        }
        molresponse_params = ResponseParameters(world, this->parser);
        auto molresponse_json = molresponse_params.to_json_if_precedence("defined");
        if (!molresponse_json.is_null())
        {
            all_input_json["response"] = molresponse_json;
        }
    }

    /**
     * Sets the calculation flags for moldft and response
     * based on the input json
     *
     */
    void set_calc_flags()
    {
        if (all_input_json.contains("dft"))
        {
            moldft_params.from_json(all_input_json["dft"]);
            run_moldft = true;
        }
        if (all_input_json.contains("response"))
        {
            molresponse_params.from_json(all_input_json["response"]);
            run_response = true;
        }
    }

    /**
     * Broadcasts the parameters to all the nodes from 0
     *
     * @param world
     */

    explicit ParameterManager(World &world, const path &input_file)
    {
        // First read the molecule file because we have one

        std::ifstream input_file_stream(input_file);
        bool is_json = json::accept(input_file_stream);
        input_file_stream.close();
        if (world.rank() == 0)
        {
            print("Input file path: ", input_file);
            print("Input file is json: ", is_json);
        }

        if (is_json)
        {
            input_file_stream.open(input_file);
            all_input_json = json::parse(input_file_stream);
            input_file_stream.close();
            read_chem_params_from_json(all_input_json);
            read_molecule_from_json_and_orient(all_input_json);
        }
        else
        {
            read_molecule_file(input_file);
            all_input_json["molecule"] = molecule.to_json();
            read_chem_params_and_write_json(world, input_file);
        }

        set_calc_flags();

        if (world.rank() == 0)
        {
            print(all_input_json.dump(4));
        }
    }

    void read_molecule_file(const path &mol_file)
    {
        bool is_json = is_json_file(mol_file);
        if (is_json)
        {
            print("Molecule file is json");

            std::ifstream mol_stream(mol_file);
            auto mol_json = json::parse(mol_stream);
            molecule.from_json(mol_json["molecule"]);
            mol_stream.close();
        }
        else
        {
            std::ifstream mol_stream(mol_file);
            molecule.read(mol_stream);
            mol_stream.close();
        }
        all_input_json["molecule"] = molecule.to_json();
    }

    [[nodiscard]] bool is_json_file(const path &input_file) const
    {
        std::ifstream input_file_stream(input_file);
        bool is_json = json::accept(input_file_stream);
        input_file_stream.close();
        return is_json;
    }

    json read_json_file(const path &input_file)
    {
        std::ifstream input_file_stream(input_file);
        auto j = json::parse(input_file_stream);
        input_file_stream.close();
        return j;
    }

    // The intention of this second constructor is to allow the user to
    // specify the molecule file and the input file separately.  The specific
    // use cases are when one wants to create a database of molecule calculations
    // all with the same basic input.
    explicit ParameterManager(World &world, const std::pair<path, path> &input_files)
    {
        auto [input_file, mol_file] = input_files;

        read_molecule_file(mol_file);
        if (world.rank() == 0)
        {
            print("Input file path: ", input_file);
        }
        auto is_json = is_json_file(input_file);

        if (is_json)
        {
            if (world.rank() == 0)
            {
                print(input_file, " is a json file");
            }
            auto input_j = read_json_file(input_file);
            read_chem_params_from_json(input_j);
        }
        else
        {
            if (world.rank() == 0)
            {
                print(input_file, " is not json file");
            }
            read_chem_params_and_write_json(world, input_file);
        }
        set_calc_flags();
        if (world.rank() == 0)
        {
            print("Input files: ", input_file, mol_file);
            print(all_input_json.dump(4));
        }
    }

    [[nodiscard]] json get_input_json() const { return all_input_json; }

    void write_input_file(std::ostream &os) const
    {
        write_json_to_input_file(all_input_json, {"dft", "response"}, os);
        write_molecule_json_to_input_file(all_input_json["molecule"], os);
    }

    void write_json_input(std::ostream &os) const { os << all_input_json.dump(4); }

    void write_moldft_input_file(std::ostream &os) const { write_moldft_input(all_input_json, os); }

    void write_response_input_file(std::ostream &os) const { write_response_input(all_input_json, os); }

    [[nodiscard]] auto get_moldft_params() const -> const CalculationParameters & { return moldft_params; }
    [[nodiscard]] auto get_molresponse_params() const -> const ResponseParameters & { return molresponse_params; }
    [[nodiscard]] auto get_molecule() const -> const Molecule & { return molecule; }

    void write_moldft_json(std::ostream &os)
    {
        os << std::setw(4) << all_input_json["dft"];
        os << std::setw(4) << all_input_json["molecule"];
    }

    void write_response_json(std::ostream &os) { os << std::setw(4) << all_input_json["response"]; }
    [[nodiscard]] bool get_run_moldft() const { return run_moldft; }
    [[nodiscard]] bool get_run_response() const { return run_response; }
};

// create a helper class for checking equivalence of two parameter class
//
bool operator==(const ParameterManager &lhs, const ParameterManager &rhs) { return lhs.get_input_json() == rhs.get_input_json(); }

class ResponseCalcManager
{
    path root; // root directory
    ParameterManager parameter_manager;

    path moldft_path; // base moldft path
    path moldft_restart; // path to moldft restart file
    path moldft_calc_info_path;
    path moldft_scf_calc_info_path;
    //
    path calc_path_json_path; // path to the calc path json file
    json calc_path_json = {}; // json object to store the calc paths
    //
    vector<path> response_paths; // paths to response calculations
    vector<path> response_restart_paths;
    vector<path> response_outfiles_paths;

    path quadratic_json_path; // path to the quadratic json file
                              // paths to the moldft and response calculations

    std::string op;
    std::string xc;
    std::vector<double> freq;
    Molecule molecule;

public:
    [[nodiscard]] auto get_calc_path_json() const -> const json & { return calc_path_json; }
    [[nodiscard]] auto get_root() const -> const path & { return root; }
    [[nodiscard]] auto get_moldft_path() const -> const path & { return moldft_path; }
    [[nodiscard]] auto get_moldft_restart() const -> const path & { return moldft_restart; }
    [[nodiscard]] auto get_moldft_calc_info_path() const -> const path & { return moldft_calc_info_path; }
    [[nodiscard]] auto get_moldft_scf_calc_info() const -> json { return moldft_scf_calc_info_path; }
    [[nodiscard]] auto get_op() const -> const std::string & { return op; }
    [[nodiscard]] auto get_xc() const -> const std::string & { return xc; }
    [[nodiscard]] auto get_freq() const -> const std::vector<double> & { return freq; }
    [[nodiscard]] auto get_molecule() const -> const Molecule & { return molecule; }
    [[nodiscard]] auto get_quadratic_path() const -> const path & { return quadratic_json_path; }

    // TODO: I should be allowed to run multiple xc and ops and freqs as well as excited states
    explicit ResponseCalcManager(World &world, ParameterManager pm) : parameter_manager(std::move(pm))
    {
        xc = parameter_manager.get_moldft_params().xc();
        op = parameter_manager.get_molresponse_params().perturbation();
        freq = parameter_manager.get_molresponse_params().freq_range();

        // if we are doing a quadratic response calculation then the freq range is considered as input frequencies and we compute
        // the linear response at all possible sums of the input frequencies
        //
        //


        if (parameter_manager.get_molresponse_params().quadratic())
        {

            vector<double> freqs_copy = freq;
            auto num_freqs = freq.size();

            auto compare_freqs = [](double x, double y) { return std::abs(x - y) < 1e-3; };

            for (int i = 0; i < num_freqs; i++)
            {
                for (int j = i; j < num_freqs; j++)
                {
                    auto omega_1 = freq[i];
                    auto omega_2 = freq[j];
                    auto omega_3 = omega_1 + omega_2;

                    if (std::find_if(freqs_copy.begin(), freqs_copy.end(), [&](double x) { return compare_freqs(x, omega_3); }) != freqs_copy.end())
                    {
                        continue;
                    }
                    if (omega_2 == 0.0)
                        continue;
                    freqs_copy.push_back(omega_3);
                }
            }


            freq = freqs_copy;
            std::sort(freq.begin(), freq.end());
            // only unique frequencies
            freq.erase(std::unique(freq.begin(), freq.end()), freq.end());
        }


        root = std::filesystem::current_path();

        json input_json = parameter_manager.get_input_json();
        input_json.erase("response");

        std::hash<json> json_hash; // hash the json object
        // create a hash of the json object using only molecule and dft parameters
        auto hash = json_hash(input_json);
        std::string output_directory = "output"; // + std::to_string(hash);

        moldft_path = root / output_directory;
        moldft_restart = moldft_path / "moldft.restartdata.00000";
        moldft_calc_info_path = moldft_path / "moldft.calc_info.json";
        moldft_scf_calc_info_path = moldft_path / "moldft.scf_info.json";

        auto define_response_paths = [&]()
        {
            for (const auto &freq : freq)
            {
                response_paths.push_back(generate_response_frequency_run_path(freq));
                auto [restart_path, restart_string] = generate_frequency_save_path(response_paths.back());
                response_restart_paths.push_back(response_paths.back() / restart_path);
            }
        };

        if (parameter_manager.get_run_response())
        {
            define_response_paths();
        }

        if (parameter_manager.get_molresponse_params().quadratic())
        {
            quadratic_json_path = moldft_path / "beta.json";
        }

        auto create_calc_path_json = [&]()
        {
            calc_path_json["moldft"] = {};
            calc_path_json["moldft"]["calc_dirs"] = getAbsolutePath(moldft_path, root);
            calc_path_json["moldft"]["restart"] = getAbsolutePath(moldft_restart, root);
            calc_path_json["moldft"]["outfiles"] = {};

            calc_path_json["moldft"]["outfiles"]["calc_info"] = getAbsolutePath(moldft_calc_info_path, root);
            calc_path_json["moldft"]["outfiles"]["scf_info"] = getAbsolutePath(moldft_scf_calc_info_path, root);

            calc_path_json["response"] = {};
            calc_path_json["response"]["calc_dirs"] = {};
            for (const auto &path : response_paths)
            {
                calc_path_json["response"]["calc_dirs"].push_back(getAbsolutePath(path, root));
            }
            calc_path_json["response"]["restarts"] = {};
            for (const auto &path : response_restart_paths)
            {
                calc_path_json["response"]["restarts"].push_back(getAbsolutePath(path, root));
            }
            calc_path_json["response"]["outfiles"] = {};
            for (const auto &path : response_paths)
            {
                calc_path_json["response"]["outfiles"].push_back(getAbsolutePath(path / "response_base.json", root));
            }

            calc_path_json["quadratic_json_path"] = getAbsolutePath(quadratic_json_path, root);
        };

        create_calc_path_json();
        output_calc_path_json();

        if (world.rank() == 0)
        {
            print();
        }
    }

    void print() const
    {
        ::print("------------Database Runner---------------");
        ::print("Root: ", root);
        ::print("MOLDFT Directory: ", moldft_path);
        ::print("----------------- Moldft Paths --------------------");
        ::print("moldft path :", moldft_path);
        ::print("moldft restart path :", moldft_restart);
        ::print("moldft calc_info json path :", moldft_calc_info_path);
        ::print("moldft scf_info json path :", moldft_calc_info_path);
        ::print("----------------- Parameters --------------------");
        ::print("xc: ", xc);
        ::print("op: ", op);
        ::print("freq: ", freq);
        ::print("------------------------------------------------");
        ::print("Response Paths: ");
        for (const auto &path : response_paths)
        {
            ::print(path);
        }
        ::print("Quadratic Path: ", quadratic_json_path);
        ::print("------------------------------------------------");
        ::print("Calc Path Json: ");
        ::print(calc_path_json.dump(4));
    }

    [[nodiscard]] path relative_to_root(const path &p) const { return std::filesystem::relative(p, root); }

    path getAbsolutePath(const path &relativePath, const path &rootDirectory) const
    {
        return std::filesystem::relative(relativePath, rootDirectory);
        return std::filesystem::weakly_canonical(std::filesystem::relative(relativePath, rootDirectory));
    }

    void output_calc_path_json() const
    {
        std::ofstream ofs("calc_path.json");
        ofs << std::setw(4) << calc_path_json;
        ofs.close();
    }

    /**
     * Runs moldft in the path provided.  Also generates the moldft input
     * file_name in the directory provided.
     *
     * @param world
     * @param moldft_path
     * @param moldft_filename
     * @param xc
     */
    void run_moldft(World &world, bool restart) const
    {
        // first thing to do is change the current directory to the moldft path
        // create the moldft path if it does not exist

        if (!std::filesystem::exists(moldft_path))
        {
            std::filesystem::create_directory(moldft_path);
        }

        std::filesystem::current_path(moldft_path);

        json calcInfo;
        auto param1 = parameter_manager.get_moldft_params();
        world.gop.broadcast_serializable(param1, 0);

        if (world.rank() == 0)
        {
            ::print("-------------Running moldft------------");
        }

        if (std::filesystem::exists(moldft_restart) && std::filesystem::exists(moldft_calc_info_path))
        {
            // if both exist, read the calc_info json
            std::ifstream ifs(moldft_calc_info_path);

            auto moldft_calc_info = json::parse(ifs);
            if (world.rank() == 0)
            {
                std::cout << "time: " << moldft_calc_info["time"] << std::endl;
                std::cout << "MOLDFT return energy: " << moldft_calc_info["return_energy"] << std::endl;
            }
        }
        else
        {
            // if params are different run and if restart exists and if im asking to
            if (std::filesystem::exists(moldft_restart) && restart)
            {
                param1.set_user_defined_value<bool>("restart", true);
            }
            world.gop.fence();
            if (world.rank() == 0)
            {
                parameter_manager.write_moldft_input_file(std::cout);

                std::ofstream ofs("moldft.in");
                parameter_manager.write_moldft_input_file(ofs);
                ofs.close();
            }
            world.gop.fence();
            commandlineparser parser;
            parser.set_keyval("input", "moldft.in");

            if (world.rank() == 0)
                ::print("input filename: ", parser.value("input"));

            print_meminfo(world.rank(), "startup");
            FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));

            std::cout.precision(6);
            SCF calc(world, parser);
            if (world.rank() == 0)
            {
                ::print("\n\n");
                ::print(" MADNESS Hartree-Fock and Density Functional Theory "
                        "Program");
                ::print(" ----------------------------------------------------------"
                        "\n");
                calc.param.print("dft");
            }
            if (world.size() > 1)
            {
                calc.set_protocol<3>(world, 1e-4);
                calc.make_nuclear_potential(world);
                calc.initial_load_bal(world);
            }
            // vama
            calc.set_protocol<3>(world, calc.param.protocol()[0]);
            // calc.set_protocol<3>(world, 1e-4);
            world.gop.fence();
            MolecularEnergy ME(world, calc);
            world.gop.fence();
            // double energy=ME.value(calc.molecule.get_all_coords().flat()); // ugh!
            ME.value(calc.molecule.get_all_coords().flat()); // ugh!
            world.gop.fence();
            const real_function_3d rho = 2.0 * calc.make_density(world, calc.aocc, calc.amo);
            auto dipole_t = calc.dipole(world, rho);
            std::map<std::string, double> results;
            results["scf_energy"] = calc.current_energy;
            world.gop.fence();
            if (world.rank() == 0)
            {
                calc.output_scf_info_schema(results, dipole_t);
                ME.output_calc_info_schema();
            }
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
    [[nodiscard]] auto generate_response_frequency_run_path(const double &frequency) const -> std::filesystem::path
    {
        std::string s_frequency = std::to_string(frequency);
        auto sp = s_frequency.find('.');
        s_frequency = s_frequency.replace(sp, sp, "-");
        std::string run_name = op + "_" + xc + "_" + s_frequency;
        return moldft_path / std::filesystem::path(run_name);
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
    auto set_frequency_path_and_restart(World &world, const double &frequency, std::filesystem::path &restart_path, bool restart, ResponseParameters &parameters) const -> std::filesystem::path
    {
        if (world.rank() == 0)
        {
            ::print("restart path", restart_path);
        }
        // change the logic create save path
        auto frequency_run_path = generate_response_frequency_run_path(frequency);
        world.gop.fence();
        if (world.rank() == 0)
        {
            ::print("frequency run path", frequency_run_path);
        }
        if (std::filesystem::is_directory(frequency_run_path))
        {
            if (world.rank() == 0)
            {
                cout << "Response directory found " << std::endl;
            }
        }
        else
        { // create the file
            std::filesystem::create_directory(frequency_run_path);
            if (world.rank() == 0)
            {
                cout << "Creating response_path directory" << std::endl;
            }
        }
        std::filesystem::current_path(frequency_run_path);
        // frequency save path
        auto [save_path, save_string] = generate_frequency_save_path(frequency_run_path);
        if (world.rank() == 0)
        {
            ::print("save path", save_path);
        }
        if (world.rank() == 0)
        {
            ::print("save string", save_string);
        }

        if (world.rank() == 0)
        {
            parameters.set_user_defined_value("save", true);
            parameters.set_user_defined_value("save_file", save_string);
            if (restart)
            { // if we are trying a restart calculation
                if (std::filesystem::exists(save_path))
                {
                    // if the save path exists then we know we can
                    //  restart from the previous save
                    parameters.set_user_defined_value("restart", true);
                    parameters.set_user_defined_value("restart_file", save_string);
                }
                else if (std::filesystem::exists(restart_path))
                {
                    ::print("restart path exists", restart_path);
                    parameters.set_user_defined_value("restart", true);
                    ::print(restart_path.parent_path().stem());
                    ::print(restart_path.filename().stem());

                    // get the directory of the restart path
                    auto new_restart_path = path("../") / restart_path.parent_path().stem() / restart_path.filename().stem();

                    // format restart path to be ../restart_path/restart_path

                    //
                    ::print("new restart file: ", restart_path);
                    parameters.set_user_defined_value("restart_file", new_restart_path.string());
                    // Then we restart from the previous file instead
                }
                else
                {
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
    auto RunResponse(World &world, const std::string &filename, double frequency, std::filesystem::path restart_path) -> std::pair<std::filesystem::path, bool>
    {
        // Set the response parameters
        ResponseParameters r_params = parameter_manager.get_molresponse_params();
        r_params.set_user_defined_value("omega", frequency);

        auto save_path = set_frequency_path_and_restart(world, frequency, restart_path, true, r_params);

        if (world.rank() == 0)
        {
            molresponse::write_response_input(r_params, filename);
        }
        // if rbase exists and converged I just return save path and true
        if (std::filesystem::exists("response_base.json"))
        {
            std::ifstream ifs("response_base.json");
            json response_base;
            ifs >> response_base;
            if (response_base["converged"] && response_base["precision"]["dconv"] == r_params.dconv())
            {
                return {save_path, true};
            }
        }
        auto calc_params = initialize_calc_params(world, std::string(filename));
        RHS_Generator rhs_generator;
        if (op == "dipole")
        {
            rhs_generator = dipole_generator;
        }
        else
        {
            rhs_generator = nuclear_generator;
        }
        FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));
        FrequencyResponse calc(world, calc_params, frequency, rhs_generator);
        if (world.rank() == 0)
        {
            ::print("\n\n");
            ::print(" MADNESS Time-Dependent Density Functional Theory Response "
                    "Program");
            ::print(" ----------------------------------------------------------\n");
            ::print("\n");
            calc_params.molecule.print();
            ::print("\n");
            calc_params.response_parameters.print("response");
            // put the response parameters in a j_molrespone json object
            calc_params.response_parameters.to_json(calc.j_molresponse);
        }
        calc.solve(world);
        world.gop.fence();
        // set protocol to the first
        if (world.rank() == 0)
        {
            // calc.time_data.to_json(calc.j_molresponse);
            calc.output_json();
        }
        // calc.time_data.print_data();
        return {save_path, calc.j_molresponse["converged"]};
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
    void run_molresponse(World &world)
    {
        std::filesystem::current_path(moldft_path);
        // add a restart path
        auto restart_path = addPath(moldft_path, "/" + op + "_0-000000.00000/restart_" + op + "_0-000000.00000");
        std::pair<std::filesystem::path, bool> success{moldft_path, false};
        bool first = true;


        for (const auto &freq : freq)
        {
            if (world.rank() == 0)
            {
                ::print(success.second);
            }
            std::filesystem::current_path(moldft_path);
            if (first)
            {
                first = false;
            }
            else if (success.second)
            {
                // if the previous run succeeded then set the restart path
                restart_path = success.first;
                if (world.rank() == 0)
                {
                    ::print("restart_path", restart_path);
                }
            }
            else
            {
                throw Response_Convergence_Error{};
            }
            success = RunResponse(world, "response.in", freq, restart_path);
            //                      high_prec);
            if (world.rank() == 0)
            {
                ::print("Frequency ", freq, " completed");
            }
        }
    }

    void run_quadratic_response(World &world)
    {
        std::filesystem::current_path(moldft_path);

        bool run_first_order = false;

        // add a restart path
        auto restart_path = addPath(moldft_path, "/" + op + "_0-000000.00000/restart_" + op + "_0-000000.00000");
        try
        {
            for (const auto &freq : freq)
            {
                std::filesystem::current_path(moldft_path);
                ResponseParameters r_params{};
                auto save_path = set_frequency_path_and_restart(world, freq, restart_path, true, r_params);

                if (std::filesystem::exists("response_base.json"))
                {
                    std::ifstream ifs("response_base.json");
                    json response_base;
                    ifs >> response_base;
                    if (response_base["converged"] && response_base["precision"]["dconv"] == r_params.dconv())
                    {
                        {
                            ::print("Response calculation already converged");
                        }
                        continue;
                    }
                    else
                    {
                        if (world.rank() == 0)
                        {
                            ::print("Response calculation not converged");
                        }
                        run_first_order = true;
                        break;
                    }
                }
                if (!std::filesystem::exists(save_path))
                {
                    throw Response_Convergence_Error{};
                }
            }
            world.gop.fence();

            if (world.rank() == 0)
            {
                ::print("Running quadratic response calculations");
            }
            std::filesystem::current_path(moldft_path);
            RHS_Generator rhs_generator;
            if (op == "dipole")
            {
                rhs_generator = dipole_generator;
            }
            else
            {
                rhs_generator = nuclear_generator;
            }
            if (world.rank() == 0)
            {
                ::print("Set up rhs generator");
            }

            auto set_hyperpolarizability_parameters = [&]()
            {
                ResponseParameters quad_parameters{};

                auto moldft_params = parameter_manager.get_moldft_params();
                auto molresponse_params = parameter_manager.get_molresponse_params();

                quad_parameters.set_user_defined_value("quadratic", true);
                quad_parameters.set_user_defined_value("freq_range", molresponse_params.freq_range());
                quad_parameters.set_user_defined_value("xc", moldft_params.xc());

                if (op == "dipole")
                {
                    quad_parameters.set_user_defined_value("dipole", true);
                    quad_parameters.set_derived_value<size_t>("states", 3);
                }

                auto final_protocol = molresponse_params.protocol().back();
                quad_parameters.set_user_defined_value<vector<double>>("protocol", {final_protocol});

                return quad_parameters;
            };

            auto quad_parameters = set_hyperpolarizability_parameters();

            // auto calc_params = initialize_calc_params(world,
            // std::string("quad.in"));
            commandlineparser parser;
            std::string moldft_archive = "moldft.restartdata";
            GroundStateCalculation ground_calculation{world, moldft_archive};
            if (world.rank() == 0)
            {
                ground_calculation.print_params();
            }
            Molecule molecule = ground_calculation.molecule();
            quad_parameters.set_ground_state_calculation_data(ground_calculation);
            if (world.rank() == 0)
            {
                quad_parameters.print();
            }

            world.gop.fence();
            FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));

            nlohmann::ordered_json beta_json = create_beta_json();

            if (world.rank() == 0)
            {
                molresponse::write_response_input(quad_parameters, "quad.in");
            }

            QuadraticResponse quad_calculation{
                world,
                {ground_calculation, molecule, quad_parameters},
                rhs_generator,
            };
            // if beta.json exists remove it
            if (std::filesystem::exists("beta.json"))
            {
                std::filesystem::remove("beta.json");
            }

            auto max_freq = freq.back();
            auto num_freqs = (freq.size() / 2) + 1;
            if (world.rank() == 0)
            {
                ::print("Number of frequencies: ", num_freqs);
            }


            for (int b = 0; b < num_freqs; b++)
            {
                for (int c = b; c < num_freqs; c++)
                {

                    auto generate_omega_restart_path = [&](double frequency)
                    {
                        auto linear_response_calc_path = generate_response_frequency_run_path(frequency);
                        auto restart_file_and_path = generate_frequency_save_path(linear_response_calc_path);
                        return restart_file_and_path;
                    };

                    auto omega_b = freq[b];
                    auto omega_c = freq[c];
                    auto omega_a = freq[b + c];

                    if (world.rank() == 0)
                    {
                        ::print("New combination of frequencies ", omega_a, " ", omega_b, " ", omega_c);
                        ::print("-------------------------------------------");
                    }
                    if (omega_a <= max_freq)
                    {

                        auto restartA = generate_omega_restart_path(omega_a);
                        auto restartB = generate_omega_restart_path(omega_b);
                        auto restartC = generate_omega_restart_path(omega_c);
                        if (world.rank() == 0)
                        {
                            ::print("Restart file for A", restartA.first);
                            ::print("Restart file for B", restartB.first);
                            ::print("Restart file for C", restartC.first);
                        }

                        // check if restartA exists
                        if (!std::filesystem::exists(restartA.first))
                        {
                            if (world.rank() == 0)
                            {
                                ::print("Restart file for omega_a = ", omega_a, " doesn't exist");
                            }
                            throw Response_Convergence_Error{};
                        }
                        else
                        {
                            if (world.rank() == 0)
                            {
                                ::print("Restart file for omega_a = ", omega_a, " exists");
                            }

                            std::array<double, 3> omegas{omega_a, omega_b, omega_c};

                            // remove .00000 from restartA.first

                            std::array<path, 3> restarts{restartA.first.replace_extension(""), restartB.first.replace_extension(""), restartC.first.replace_extension("")};

                            quad_calculation.set_x_data(world, omegas, restarts);
                            auto beta_abc = quad_calculation.compute_beta_v2(world);

                            // Make a table printing the beta value for each direction
                            //[XXX,XYY,XZZ,YXX,YYY,YZZ,ZXX,ZYY,ZZZ,XYZ]

                            vector<std::string> directions{"XXX", "XYY", "XZZ", "YXX", "YYY", "YZZ", "ZXX", "ZYY", "ZZZ", "XYZ"};

                            if (world.rank() == 0)
                            {
                                for (int i = 0; i < 10; i++)
                                {
                                    if (world.rank() == 0)
                                    {
                                        ::print(directions[i], " : ", beta_abc[i]);
                                    }
                                }
                            }


                            nlohmann::ordered_json beta_entry;

                            std::array<double, 10> beta_vector{};
                            std::copy(beta_abc.ptr(), beta_abc.ptr() + 10, beta_vector.begin());
                            append_to_beta_json({-1.0 * omega_a, omega_b, omega_c}, beta_vector, beta_json);

                            std::ofstream outfile("beta.json");
                            if (outfile.is_open())
                            {
                                outfile << beta_json.dump(4);
                                outfile.close();
                            }
                            else
                            {
                                std::cerr << "Failed to open file for writing: " << "beta.json" << std::endl;
                            }
                        }
                    }
                    else
                    {
                        if (world.rank() == 0)
                        {
                            ::print("Skipping omega_a = ", omega_a, " because it is greater than the max frequency");
                        }
                        continue;
                    }
                }
            }

            // write the beta json to file
            std::ofstream ofs("beta.json");
            ofs << std::setw(4) << beta_json << std::endl;
            ofs.close();
        }
        catch (Response_Convergence_Error &e)
        {
            if (true)
            {
                // if the first order response calculations haven't been run then run
                // them
                if (world.rank() == 0)
                {
                    ::print("Running first order response calculations");
                }
                run_molresponse(world);
            }
            else
            {
                if (world.rank() == 0)
                {
                    ::print("First order response calculations haven't been run and "
                            "can't be run");
                    ::print("Quadratic response calculations can't be run");
                }
            }
        }
    }
    // sets the current path to the save path
    /**
     * Generates the frequency save path with format
     * /frequency_run_path/restart_[frequency_run_filename].00000
     *
     * @param frequency_run_path
     * @return
     */
    static auto generate_frequency_save_path(const std::filesystem::path &frequency_run_path) -> std::pair<std::filesystem::path, std::string>
    {
        auto save_path = std::filesystem::path(frequency_run_path);
        auto run_name = frequency_run_path.filename();
        std::string save_string = "restart_" + run_name.string();
        save_path += "/";
        save_path += save_string;
        save_path += ".00000";

        return {save_path, save_string};
    }
    // set up a function that creates a beta_json with the fields defined  below.
    // in each field there will a vector of values.

    static nlohmann::ordered_json create_beta_json()
    {
        // i need A B C to hold char values and A-freq, B-freq, C-freq to hold
        // double values

        nlohmann::ordered_json beta_json = {{"A-freq", json::array()}, {"B-freq", json::array()}, {"C-freq", json::array()}, {"A", json::array()},
                                            {"B", json::array()},      {"C", json::array()},      {"Beta", json::array()}};
        return beta_json;
    }

    // for a set of frequencies create a table from the beta values
    static void append_to_beta_json(const std::array<double, 3> &freq, const std::array<double, 10> &beta, nlohmann::ordered_json &beta_json)
    {
        // create 3 columns of directions for each A,B,C
        std::array<char, 10> direction_A{'X', 'X', 'X', 'Y', 'Y', 'Y', 'Z', 'Z', 'Z', 'X'};
        std::array<char, 10> direction_B{'X', 'Y', 'Z', 'X', 'Y', 'Z', 'X', 'Y', 'Z', 'Y'};
        std::array<char, 10> direction_C{'X', 'Y', 'Z', 'X', 'Y', 'Z', 'X', 'Y', 'Z', 'Z'};

        // append each value of the columns to the beta json
        // for each value of beta
        // capitalize the direction

        for (int i = 0; i < 10; i++)
        {
            beta_json["A-freq"].push_back(freq[0]);
            beta_json["B-freq"].push_back(freq[1]);
            beta_json["C-freq"].push_back(freq[2]);

            beta_json["A"].push_back(std::string(1, direction_A[i]));
            beta_json["B"].push_back(std::string(1, direction_B[i]));
            beta_json["C"].push_back(std::string(1, direction_C[i]));
            beta_json["Beta"].push_back(beta[i]);
        }
    }
};

#endif // MADNESS_RUNNERS_HPP
