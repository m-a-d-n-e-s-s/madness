#include <algorithm>
#include <filesystem>
#include <madchem.h>
#include <madness/chem/CalculationParameters.h>
#include <madness/chem/SCF.h>
#include <madness/chem/molecule.h>
#include <madness/external/nlohmann_json/json.hpp>
#include <apps/molresponse/response_parameters.h>
#include <utility>

using path = std::filesystem::path;
using json = nlohmann::json;
using commandlineparser = madness::commandlineparser;

using namespace madness;

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
    output_stream << "end" << std::endl
                  << std::endl;
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
        output_stream << "end" << std::endl
                      << std::endl;
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


