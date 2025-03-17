/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHO
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/
#ifndef SRC_APPS_MADQC_CALC_MANAGER_HPP_
#define SRC_APPS_MADQC_CALC_MANAGER_HPP_

// #include <madchem.h>
#include <madness/chem/CalculationParameters.h>
// #include <madness/chem/SCF.h>
#include "madqc/utils.hpp"
#include "parameter_manager.hpp"
#include <apps/molresponse/FrequencyResponse.hpp>
#include <apps/molresponse/ResponseExceptions.hpp>
#include <filesystem>
#include <madness/chem/molecule.h>
#include <madness/external/nlohmann_json/json.hpp>
#include <madness/tensor/solvers.h>
#include <madness/tensor/tensor_json.hpp>
#include <map>
#include <memory>
#include <mpi.h>
#include <utility>
#include <vector>

using json = nlohmann::json;
using path = std::filesystem::path;
using ResponseInput = std::tuple<std::string, std::string, std::vector<double>>;
namespace fs = std::filesystem;

// I define this path because this way I can be more explicit what I expect out
// of the json objects Each Calculation needs to define,
//
// 1. The calculation paths, where the calculation will be run
// 2. The restart paths, where the calculation will be restarted
// 3. The output paths, where the output of the calculation will be stored,
//
//
// For consistency, each of theses are vectors of paths
struct CalculationTemplate
{
public:
  vector<path> calc_paths = {};
  vector<path> restarts = {};
  std::map<std::string, vector<path>> outputs = {};

  void print();
};

inline void to_json(json &j, const CalculationTemplate &p)
{
  j = json{{"calculations", p.calc_paths},
           {"restarts", p.restarts},
           {"outputs", p.outputs}};
}

inline void from_json(const json &j, CalculationTemplate &p)
{
  j.at("calculations").get_to(p.calc_paths);
  j.at("restarts").get_to(p.restarts);
  j.at("outputs").get_to(p.outputs);
}

inline void CalculationTemplate::print()
{
  json j = *this;
  std::cout << j.dump(4);
}

template <typename T>
void to_json(json &j, const Tensor<T> &m)
{
  // auto dimensions = m.dims();
  long size = m.size(); ///< Number of elements in the tensor
  if (size == 0)
  {
    return;
  }
  long n_dims = m.ndim(); ///< Number of dimensions (-1=invalid; 0=no
  ///< supported; >0=tensor)
  auto dims = m.dims(); // the size of each dimension
  // auto strides = m.strides();
  auto m_vals_vector = std::vector<T>(size);
  auto m_dims_vector = std::vector<long>(n_dims);
  std::copy(&m[0], &m[0] + size, m_vals_vector.begin());
  std::copy(dims, dims + n_dims, m_dims_vector.begin());
  // This is everything we need to translate to a numpy vector...
  j["size"] = size;
  j["vals"] = m_vals_vector;
  j["dims"] = m_dims_vector;
}

template <typename T>
void from_json(const nlohmann::json &j, Tensor<T> &m)
{
  // need to be explicit here about types so we find the proper Tensor
  // constructors
  long size = j["size"];
  std::vector<T> m_vals_vector = j["vals"];
  std::vector<long> m_dims_vector = j["dims"];

  Tensor<T> flat_m(size);
  // copy the values from the vector to the flat tensor
  std::copy(m_vals_vector.begin(), m_vals_vector.end(), &flat_m[0]);
  // reshape the tensor using dimension vector
  m = flat_m.reshape(m_dims_vector);
}

// Defines the type of each possible output
template <typename T>
struct OutputTemplate
{
  double energy{};
  Tensor<T> dipole{};
  Tensor<T> gradient{};
  Tensor<T> hessian{};
  nlohmann::ordered_json alpha;
  nlohmann::ordered_json beta;
  OutputTemplate() = default;
};

template <typename T>
void to_json(json &j, const OutputTemplate<T> &p)
{

  j = json();

  j["energy"] = p.energy;
  to_json(j["dipole"], p.dipole);
  to_json(j["gradient"], p.gradient);
}

template <typename T>
void from_json(const json &j, OutputTemplate<T> &p)
{
  p.energy = j.at("energy");
  p.dipole = j.at("dipole");
  p.gradient = j.at("gradient");
  p.hessian = j.at("hessian");
  p.alpha = j.at("alpha");
  p.beta = j.at("beta");
}

class PathManager
{
private:
  json path_data;

public:
  void setPath(const std::string &key, const CalculationTemplate &value)
  {
    path_data[key] = value;
  }

  CalculationTemplate getPath(const std::string &key)
  {
    return path_data[key].get<CalculationTemplate>();
  }

  [[nodiscard]] path get_output_path() const { return path_data["outputs"]; }

  [[nodiscard]] json getAllPaths() const { return path_data; }

  void createDirectories()
  {
    json temp = path_data;
    temp.erase("outputs");
    auto cwd = fs::current_path();
    for (const auto &[key, calc_paths] : temp.items())
    {
      auto calc_path = calc_paths.get<CalculationTemplate>();

      for (const auto &calc_dir : calc_path.calc_paths)
      {
        print("Creating directory: ", calc_dir);
        if (!std::filesystem::exists(cwd / calc_dir.parent_path()))
        {
          std::filesystem::create_directory(cwd / calc_dir.parent_path());
        }
        // if base direcotry does not exist create it
        if (!std::filesystem::exists(calc_dir))
        {
          std::filesystem::create_directory(cwd / calc_dir);
        }
      }
    }
  }
  // This function is used to set the paths for a specific calculation
  void createDirectories(const std::string &key)
  {
    auto calc_paths = path_data[key].get<CalculationTemplate>();
    for (const auto &calc_dir : calc_paths.calc_paths)
    {
      if (!std::filesystem::exists(calc_dir))
      {
        std::filesystem::create_directory(fs::current_path() / calc_dir);
      }
    }
  }

  explicit PathManager(const path &root)
  {
    path_data["outputs"] = root / "outputs.json";
  }
};

// A word for a class which takes input from another class and
// therefore has to define the input interface
class InputInterface
{

protected:
  vector<std::string> input_names;

public:
  explicit InputInterface(const std::vector<std::string> &names)
      : input_names(names) {}
};

class CalculationStrategy
{
public:
  virtual void setPaths(PathManager &path_manager, const path &root) = 0;
  virtual void compute(World &world, PathManager &path_manager,
                       const path &root, const Tensor<double> &coords) = 0;
  virtual ~CalculationStrategy() = default;

  CalculationStrategy(std::string calc_name,
                      std::map<std::string, bool> properties)
      : name(std::move(calc_name)),
        requested_properties(std::move(properties)) {}
  CalculationStrategy() = default;
  [[nodiscard]] virtual std::unique_ptr<CalculationStrategy> clone() const = 0;

protected:
  std::string name;
  std::map<std::string, bool> requested_properties;
};

class CompositeCalculationStrategy : public CalculationStrategy
{
private:
  std::vector<std::unique_ptr<CalculationStrategy>> strategies;

public:
  void addStrategy(std::unique_ptr<CalculationStrategy> strategy)
  {
    strategies.push_back(std::move(strategy));
  }

  std::vector<std::unique_ptr<CalculationStrategy>> &getStrategies()
  {
    return strategies;
  }

  void setPaths(PathManager &path_manager, const path &root) override
  {
    for (const auto &strategy : strategies)
    {
      strategy->setPaths(path_manager, root);
    }
  }

  [[nodiscard]] std::unique_ptr<CalculationStrategy> clone() const override
  {
    auto new_strategy = std::make_unique<CompositeCalculationStrategy>();

    for (const auto &strategy : strategies)
    {
      new_strategy->addStrategy(strategy->clone());
    }
    return new_strategy;
  }

  void compute(World &world, PathManager &path_manager, const path &root,
               const Tensor<double> &coords) override
  {
    for (const auto &strategy : strategies)
    {
      strategy->compute(world, path_manager, root, coords);
    }
  }
};

class MoldftCalculationStrategy : public CalculationStrategy
{
  CalculationParameters parameters;
  Molecule molecule;
  json paths;
  path output_path;

  std::vector<std::string> available_properties = {"energy", "gradient",
                                                   "dipole"};

public:
  [[nodiscard]] std::unique_ptr<CalculationStrategy> clone() const override
  {
    return std::make_unique<MoldftCalculationStrategy>(*this);
  }
  MoldftCalculationStrategy(const CalculationParameters &params, Molecule mol,
                            std::string calc_name = "moldft",
                            std::map<std::string, bool> properties = {{"energy",
                                                                       true}})
      : parameters(params), molecule(std::move(mol)),
        CalculationStrategy(std::move(calc_name), std::move(properties)) {}

  void setPaths(PathManager &path_manager, const path &root) override
  {
    // Build the paths for the moldf calculation
    //

    CalculationTemplate moldft = {};
    auto base_path = root / name;
    base_path = fs::relative(base_path, root);

    moldft.calc_paths.push_back(base_path);
    moldft.restarts.push_back(base_path / "moldft.restartdata.00000");
    moldft.outputs["calc_info"] = {base_path / "moldft.calc_info.json"};
    moldft.outputs["scf_info"] = {base_path / "moldft.calc_info.json"};

    path_manager.setPath(base_path, moldft);
  }

  void compute(World &world, PathManager &path_manager, const path &root,
               const Tensor<double> &coords) override
  {
    // Get the paths for the moldft calculation by name
    /*auto base_path = root / name;*/
    /*base_path = fs::relative(base_path, root);*/
    auto moldft_paths = path_manager.getPath(name);
    auto moldft_path = moldft_paths.calc_paths[0];
    auto restart_path = moldft_paths.restarts[0]; //
    path calc_info_path = moldft_paths.outputs["calc_info"][0];
    // Set the current path to the moldft path
    std::filesystem::current_path(moldft_path);

    json calcInfo;
    auto param1 = parameters;
    world.gop.broadcast_serializable(param1, 0);
    world.gop.fence();

    if (world.rank() == 0)
    {
      ::print("-------------Running moldft------------");
    }

    bool run_moldft = true;

    OutputTemplate<double> otemp;
    // Adjust the parameters for the calculation
    // if restart and calc_info exists the read the calc_info json
    if (std::filesystem::exists(root / restart_path) &&
        std::filesystem::exists(root / calc_info_path))
    {
      // if both exist, read the calc_info json
      //

      std::ifstream ifs(root / calc_info_path);
      auto moldft_calc_info = json::parse(ifs);
      if (world.rank() == 0)
      {
        std::cout << "time: " << moldft_calc_info["time_tag"] << std::endl;
        std::cout << "MOLDFT return energy: "
                  << moldft_calc_info["return_energy"] << std::endl;
      }
      // read in the molecule from the moldft_calc_info
      Molecule read_molecule;
      if (moldft_calc_info.contains("molecule"))
      {

        read_molecule.from_json(moldft_calc_info["molecule"]);
        auto last_coords = read_molecule.get_all_coords().flat();
        auto molecule_changed = (last_coords - coords).normf() > 1e-6;
        if (world.rank() == 0)
        {
          print("Last molecule: ", last_coords);
          print("Current molecule: ", coords);
          print("Molecule changed: ", molecule_changed);
        }
        if (molecule_changed)
        {
          if (world.rank() == 0)
          {
            ::print("Molecule has changed, restarting calculation");
          }
          param1.set_user_defined_value<bool>("restart", true);
          run_moldft = true;
        }
        else
        {
          if (world.rank() == 0)
          {
            ::print("Molecule has not changed, skipping calculation");
          }
          run_moldft = false;

          for (const auto &[property, requested] : requested_properties)
          {
            if (requested)
            {
              if (property == "energy")
              {
                otemp.energy = moldft_calc_info["return_energy"];
              }
              if (property == "gradient")
              {
                otemp.gradient = moldft_calc_info["gradient"];
              }
              if (property == "dipole")
              {
                otemp.dipole = moldft_calc_info["dipole"];
              }
            }
          }
        }
      }
      else
      {
        run_moldft = false;
      }
    }
    if (run_moldft)
    {
      {
        // if params are different run and if restart exists and if im asking to
        if (std::filesystem::exists(root / restart_path))
        {
          param1.set_user_defined_value<bool>("restart", true);
        }
        world.gop.fence();
        if (world.rank() == 0)
        {
          json moldft_input_json = {};
          moldft_input_json["dft"] =
              parameters.to_json_if_precedence("defined");
          moldft_input_json["molecule"] = molecule.to_json();
          print("moldft_input_json: ", moldft_input_json.dump(4));
          std::ofstream ofs("moldft.in");
          write_moldft_input(moldft_input_json, ofs);
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
        calc.set_protocol<3>(world, calc.param.protocol()[0]);
        MolecularEnergy E(world, calc);
        double energy = E.value(coords); // ugh!

        world.gop.fence();
        if (world.rank() == 0)
        {
          calc.output_calc_info_schema();
        }
        world.gop.fence();
        properties(world, calc, energy, path_manager);
      }
    }
    else
    {
      if (world.rank() == 0)
      {
        ::print("Skipping moldft calculation");
      }
      if (world.rank() == 0)
      {
        json persistent_output = {};
        json output = {};
        to_json<double>(output, otemp);
        // print("output: ", output.dump(4));
        if (std::filesystem::exists(path_manager.get_output_path()))
        {
          std::ifstream ifs(path_manager.get_output_path());
          ifs >> persistent_output;
          ifs.close();
        }
        output["molecule"] = molecule.to_json();
        persistent_output[name] = output;

        // print("output: ", output.dump(4));
        std::ofstream ofs(path_manager.get_output_path());
        ofs << persistent_output.dump(4);
        ofs.close();
      }

      // still output properites into output.json
    }

    // Add actual calculation logic here
  }
  template <typename T>
  static void compute_property(World &world, const std::string &property,
                               SCF &calc, double energy,
                               OutputTemplate<T> &result)
  {
    if (property == "energy")
    {
      result.energy = energy;
    }
    else if (property == "gradient")
    {
      auto gradient = calc.derivatives(
          world, calc.make_density(world, calc.aocc, calc.amo));
      result.gradient = gradient;
      // we need to make json representation of a tensor....
    }
    else if (property == "dipole")
    {
      auto dipole =
          calc.dipole(world, calc.make_density(world, calc.aocc, calc.amo));
      result.dipole = dipole;
    }
    else
    {
      throw std::runtime_error("Property not available");
    }
  }

  void properties(World &world, SCF &calc, double energy, PathManager &pm)
  {
    // This is where a property interface would be useful. for properties use
    // the properties interface
    OutputTemplate<double> otemp;
    // to get the properties of the calculation
    json results;
    paths[name]["output"]["properties"] = {};

    for (const auto &[property, compute] : requested_properties)
    {
      if (compute)
      {
        compute_property(world, property, calc, energy, otemp);
      }
    }

    json output = {};
    to_json<double>(output, otemp);

    // Read the output json and write it to the output path

    if (world.rank() == 0)
    {
      json persistent_output = {};
      // print("output: ", output.dump(4));
      if (std::filesystem::exists(pm.get_output_path()))
      {
        std::ifstream ifs(pm.get_output_path());
        ifs >> persistent_output;
        ifs.close();
      }
      output["molecule"] = molecule.to_json();
      persistent_output[name] = output;

      // print("output: ", output.dump(4));
      std::ofstream ofs(pm.get_output_path());
      ofs << persistent_output.dump(4);
      ofs.close();
    }
  }
};

class ResponseConfig
{
public:
  std::string calc_name = "response";
  std::string perturbation;
  std::string xc;
  std::vector<double> frequencies;
  static path restart_path(const std::filesystem::path &calc_path)
  {
    auto save_path = std::filesystem::path(calc_path);
    auto run_name = calc_path.filename();
    std::string save_string = "restart_" + run_name.string();
    save_path += "/";
    save_path += save_string;
    save_path += ".00000";

    return save_path;
  }

  [[nodiscard]] auto calc_path(const path &root, const double &frequency) const
      -> std::filesystem::path
  {
    std::string s_frequency = std::to_string(frequency);
    auto sp = s_frequency.find('.');
    s_frequency = s_frequency.replace(sp, sp, "-");
    std::string run_name =
        this->perturbation + "_" + this->xc + "_" + s_frequency;
    return root / std::filesystem::path(run_name);
  }
  explicit ResponseConfig(ResponseInput input,
                          std::string calc_name = "response")
      : calc_name(std::move(calc_name)), perturbation(std::get<0>(input)),
        xc(std::get<1>(input)), frequencies(std::get<2>(input)) {}
};

class LinearResponseStrategy : public CalculationStrategy, InputInterface
{
  ResponseParameters parameters;
  std::string op;
  std::string xc;
  std::vector<double> freqs;
  ResponseConfig config;
  path output_path;
  std::vector<std::string> available_properties = {"alpha"};

public:
  explicit LinearResponseStrategy(
      const ResponseParameters &params, const ResponseInput &r_input,
      std::string name = "response",
      const std::vector<std::string> &input_names = {"moldft"})
      : parameters(params), config(r_input, std::move(name)),
        CalculationStrategy(name, {{"alpha", true}}),
        InputInterface(input_names) {}
  [[nodiscard]] std::unique_ptr<CalculationStrategy> clone() const override
  {
    return std::make_unique<LinearResponseStrategy>(*this);
  }

  void setPaths(PathManager &path_manager, const path &root) override
  {
    // We always start at root+calc_name
    auto base_path = root / name;
    base_path = fs::relative(base_path, root);
    CalculationTemplate response_paths;

    response_paths.outputs["response_base"] = {};

    for (const auto &frequency : config.frequencies)
    {
      auto frequency_run_path = config.calc_path(base_path, frequency);
      auto restart_path = ResponseConfig::restart_path(frequency_run_path);
      response_paths.calc_paths.push_back(frequency_run_path);
      response_paths.restarts.push_back(restart_path);
      response_paths.outputs["response_base"].push_back(frequency_run_path /
                                                        "response_base.json");
    }
    response_paths.outputs["alpha"].push_back(base_path / "alpha.json");
    path_manager.setPath(base_path, response_paths);
  }

  static void append_to_alpha_json(const double &omega,
                                   const std::vector<std::string> &ij,
                                   const Tensor<double> &alpha,
                                   nlohmann::ordered_json &alpha_json)
  {
    auto num_unique_elements = ij.size();
    for (int i = 0; i < num_unique_elements; i++)
    {
      alpha_json["omega"].push_back(omega);
      alpha_json["ij"].push_back(ij[i]);
      alpha_json["alpha"].push_back(alpha[i]);
    }
  }
  static void add_alpha_i_to_json(nlohmann::ordered_json &alpha_i,
                                  nlohmann::ordered_json &alpha_json)
  {
    // print(alpha_json.dump(4));

    alpha_json["omega"].insert(alpha_json["omega"].end(),
                               alpha_i["omega"].begin(),
                               alpha_i["omega"].end());
    alpha_json["ij"].insert(alpha_json["ij"].end(), alpha_i["ij"].begin(),
                            alpha_i["ij"].end());
    alpha_json["alpha"].insert(alpha_json["alpha"].end(),
                               alpha_i["alpha"].begin(),
                               alpha_i["alpha"].end());
  }

  bool runFrequency(World &world, ResponseParameters &r_params,
                    GroundStateCalculation &ground_calculation,
                    const Molecule &molecule)
  {
    auto op = r_params.perturbation();
    if (world.rank() == 0)
    {
      json input_json = {};
      input_json["response"] = r_params.to_json_if_precedence("defined");
      std::ofstream out("response.in");
      write_json_to_input_file(input_json, {"response"}, out);
    }
    bool converged = false;
    // if rbase exists and converged I just return save path and true
    if (world.rank() == 0)
    {
      if (std::filesystem::exists("response_base.json"))
      {
        {
          std::ifstream ifs("response_base.json");
          json response_base;
          ifs >> response_base;
          if (response_base["converged"] &&
              response_base["precision"]["dconv"] == r_params.dconv())
            converged = true;
        }
      }
    }
    world.gop.fence();
    world.gop.broadcast(converged, 0);
    if (converged)
    {
      return true;
    }
    else
    {

      r_params.set_ground_state_calculation_data(ground_calculation);
      r_params.set_derived_values(world, molecule);
      if (r_params.omega() == 0)
      {
        r_params.set_derived_value<std::string>("calc_type", "static");
      }
      else
      {
        r_params.set_derived_value<std::string>("calc_type", "full");
      }
      CalcParams calc_params = {ground_calculation, molecule, r_params};
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
      FrequencyResponse calc(world, calc_params, rhs_generator);
      if (world.rank() == 0)
      {
        ::print("\n\n");
        ::print(" MADNESS Time-Dependent Density Functional Theory Response "
                "Program");
        ::print(
            " ----------------------------------------------------------\n");
        calc_params.response_parameters.to_json(calc.j_molresponse);
      }
      calc.solve(world);
      world.gop.fence();

      if (world.rank() == 0)
      {
        auto [omega, polar_omega] = calc.get_response_data();
        // flatten polar_omega

        auto alpha = polar_omega.flat();

        std::vector<std::string> ij{"XX", "XY", "XZ", "YX", "YY",
                                    "YZ", "ZX", "ZY", "ZZ"};
        nlohmann::ordered_json alpha_json;
        alpha_json["omega"] = {};
        alpha_json["ij"] = {};
        alpha_json["alpha"] = {};

        append_to_alpha_json(omega, ij, alpha, alpha_json);
        calc.j_molresponse["properties"] = {};
        calc.j_molresponse["properties"]["alpha"] = alpha_json;
        calc.output_json();

        converged = calc.j_molresponse["converged"];
      }
      world.gop.broadcast(converged, 0);
      world.gop.fence();
      // calc.time_data.print_data();
      return converged;
    }
  }

  void compute(World &world, PathManager &path_manager, const path &root,
               const Tensor<double> &coords) override
  {

    auto path_key = name; // access from path manager
    auto moldft_paths = path_manager.getPath(input_names[0]);
    auto response_paths = path_manager.getPath(path_key);

    std::string moldft_restart =
        root / moldft_paths.restarts[0].replace_extension("").string();
    auto alpha_outpath = response_paths.outputs["alpha"][0];

    // I need to analyze alpha.json to see which frequencies are missing.
    // Or I just write them to seperate files and then combine them later?
    auto freqs = config.frequencies;
    auto num_freqs = freqs.size();

    bool last_converged = false;
    nlohmann::ordered_json alpha_json;
    alpha_json["omega"] = json::array();
    alpha_json["ij"] = json::array();
    alpha_json["alpha"] = json::array();

    auto set_restart_params = [&](ResponseParameters &r_params,
                                  const vector<double> &freqs, const size_t &i,
                                  const CalculationTemplate &response_paths)
    {
      auto freq = freqs[i];
      auto restart_path_i = root / response_paths.restarts[i];
      bool restart = true;
      const path &save_path =
          restart_path_i; // current restart path aka save path
      path restart_path = (i > 0) ? root / response_paths.restarts[i - 1]
                                  : root / response_paths.restarts[i];
      std::string save_string = save_path.filename().stem();
      r_params.set_user_defined_value("omega", freq);
      r_params.set_user_defined_value("archive", moldft_restart);
      if (last_converged || i == 0)
      {
        r_params.set_user_defined_value("save", true);
        r_params.set_user_defined_value("save_file", save_string);
        if (restart)
        { // if we are trying a restart calculation
          if (std::filesystem::exists(save_path))
          {
            r_params.set_user_defined_value("restart", true);
            r_params.set_user_defined_value("restart_file", save_string);
          }
          else if (std::filesystem::exists(restart_path))
          {
            r_params.set_user_defined_value("restart", true);
            auto new_restart_path = path("../") /
                                    restart_path.parent_path().stem() /
                                    restart_path.filename().stem();
            r_params.set_user_defined_value("restart_file",
                                            new_restart_path.string());
          }
          else
          {
            r_params.set_user_defined_value("restart", false);
          }
        }
      }
      else
      {
        throw Response_Convergence_Error{};
      }
    };

    GroundStateCalculation ground_calculation{world, moldft_restart};
    Molecule molecule = ground_calculation.molecule();

    for (size_t i = 0; i < num_freqs; i++)
    {
      auto freq_i = freqs[i];
      auto calc_path_i = root / response_paths.calc_paths[i];
      auto restart_path_i = root / response_paths.restarts[i];
      auto response_base_i = root / response_paths.outputs["response_base"][i];

      std::filesystem::current_path(calc_path_i);

      ResponseParameters r_params = parameters;
      if (world.rank() == 0)
      {
        set_restart_params(r_params, freqs, i, response_paths);
      }
      world.gop.broadcast_serializable(r_params, 0);

      last_converged =
          runFrequency(world, r_params, ground_calculation, molecule);

      nlohmann::ordered_json alpha_i;

      if (last_converged)
      {
        // read output_paths[i]
        std::ifstream ifs(response_base_i);
        json response_base_i;
        ifs >> response_base_i;
        alpha_i = response_base_i["properties"]["alpha"];
      }
      if (world.rank() == 0)
      {
        /*print("last converged: ", last_converged);*/
        /*//print(alpha_i.dump(4));*/
      }
      add_alpha_i_to_json(alpha_i, alpha_json);
    }

    std::ofstream out_file(alpha_outpath);
    if (world.rank() == 0)
    {
      out_file << alpha_json.dump(4);
    }
    // Read in output json and write alpha.json to output path
    //
    if (world.rank() == 0)
    {
      json persistent_output = {};
      // read in the output json
      if (std::filesystem::exists(path_manager.get_output_path()))
      {
        std::ifstream ifs(path_manager.get_output_path());
        ifs >> persistent_output;
        ifs.close();
      }
      json output = {};
      output["alpha"] = alpha_json;

      persistent_output[name] = output;
      std::ofstream ofs(path_manager.get_output_path());
      ofs << persistent_output.dump(4);
      ofs.close();
    }
  }
};

struct BetaData
{
  std::array<std::string, 27> ijk;
  std::map<std::tuple<int, int, int>,
           std::pair<std::array<double, 3>, std::array<double, 27>>>
      beta_data;

  BetaData()
  {

    // This is what we expect
    auto index_B = {0, 0, 0, 1, 1, 1, 2, 2, 2};
    auto index_C = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    std::vector<std::string> xyz = {"X", "Y", "Z"};
    std::vector<std::string> jk = {"XX", "XY", "XZ", "YX", "YY",
                                   "YZ", "ZX", "ZY", "ZZ"};
    int index = 0;
    for (const auto &a : xyz)
    {
      for (const auto &bc : jk)
      {
        ijk[index++] = a + bc;
      }
    }
  };

  void add_data(
      const std::tuple<int, int, int> &abc,
      const std::pair<std::array<double, 3>, std::array<double, 27>> &data)
  {
    beta_data[abc] = data;
  }

  nlohmann::ordered_json to_json_table()
  {

    nlohmann::ordered_json beta_json;

    beta_json["Afreq"] = json::array();
    beta_json["Bfreq"] = json::array();
    beta_json["Cfreq"] = json::array();
    beta_json["A"] = json::array();
    beta_json["B"] = json::array();
    beta_json["C"] = json::array();
    beta_json["Beta"] = json::array();

    auto num_unique_elements = ijk.size();
    // print(num_unique_elements);

    // loop through beta data and add to json
    for (const auto &[index, data] : beta_data)
    {
      auto freqs = data.first;
      auto beta_i = data.second;
      auto Afreq = freqs[0];
      auto Bfreq = freqs[1];
      auto Cfreq = freqs[2];

      for (int i = 0; i < num_unique_elements; i++)
      {

        const std::string dir = ijk[i];
        // print("Dir: ", dir);
        const double beta_ijk = beta_i[i];

        std::string A{};
        std::string B{};
        std::string C{};

        A += dir[0];
        B += dir[1];
        C += dir[2];

        /*print("A: ", A);*/
        /*print("B: ", B);*/
        /*print("C: ", C);*/
        /*print("Beta: ", beta_ijk);*/
        /*print("Afreq: ", Afreq);*/
        /*print("Bfreq: ", Bfreq);*/
        /*print("Cfreq: ", Cfreq);*/

        beta_json["Afreq"].push_back(Afreq);
        beta_json["Bfreq"].push_back(Bfreq);
        beta_json["Cfreq"].push_back(Cfreq);
        beta_json["A"].push_back(A);
        beta_json["B"].push_back(B);
        beta_json["C"].push_back(C);
        beta_json["Beta"].push_back(beta_ijk);
      }
    }
    return beta_json;
  }
  [[nodiscard]] size_t length() const { return beta_data.size(); }
};

void to_json(json &j, const BetaData &p) { j = json{{"beta", p.beta_data}}; }
void from_json(const json &j, BetaData &p) { j.at("beta").get_to(p.beta_data); }

using beta_indexes = std::vector<
    std::pair<std::tuple<int, int, int>, std::tuple<double, double, double>>>;

class ResponseHyper : public CalculationStrategy, InputInterface
{
  ResponseParameters parameters;
  std::string op;
  std::string xc;
  ResponseConfig config;
  json paths;
  path output_path;
  std::vector<std::string> available_properties = {"beta"};

  beta_indexes abc_freqs;

public:
  [[nodiscard]] std::unique_ptr<CalculationStrategy> clone() const override
  {
    return std::make_unique<ResponseHyper>(*this);
  }
  explicit ResponseHyper(
      const ResponseParameters &params, const ResponseInput &r_input,
      const beta_indexes &abc_freqs = {}, const std::string &name = "hyper",
      const std::vector<std::string> &input_names = {"moldft", "response"})
      : parameters(params), abc_freqs(abc_freqs), config(r_input, name),
        InputInterface(input_names),
        CalculationStrategy(name, {{"beta", true}})
  {
    if (input_names.size() != 2)
    {
      throw std::runtime_error("ResponseHyper requires two input names");
    }
  }

  void setPaths(PathManager &path_manager, const path &root) override
  {

    auto base_path = root / name;
    base_path = fs::relative(base_path, root);
    CalculationTemplate hyper_paths;
    hyper_paths.calc_paths.push_back(base_path);
    hyper_paths.outputs["beta"] = {base_path / "beta.json",
                                   base_path / "beta2.json"};
    path_manager.setPath(base_path, hyper_paths);
  }

  void compute(World &world, PathManager &path_manager, const path &root,
               const Tensor<double> &coords) override
  {
    auto moldft_name = input_names[0];
    auto response_name = input_names[1];
    if (world.rank() == 0)
    {
      print("moldft_name: ", moldft_name);
      print("response_name: ", response_name);
    }

    auto hyper_paths = path_manager.getPath(name);
    auto calc_path = root / hyper_paths.calc_paths[0];

    auto moldft_paths = path_manager.getPath(moldft_name);
    auto response_paths = path_manager.getPath(response_name);

    auto moldft_restart = root / moldft_paths.restarts[0];
    auto beta_outpath = root / hyper_paths.outputs["beta"][0];
    auto beta_2_path = root / hyper_paths.outputs["beta"][1];
    auto response_restarts = response_paths.restarts;

    // print("freqs: ", config.frequencies);
    //  Run the calculations
    bool last_converged = false;

    BetaData beta_data;

    try
    {
      auto num_freqs = config.frequencies.size();

      if (world.rank() == 0)
      {
        ::print("Running quadratic response calculations");
      }
      std::filesystem::current_path(calc_path);
      RHS_Generator rhs_generator;
      if (config.perturbation == "dipole")
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

        auto molresponse_params = parameters;

        quad_parameters.set_user_defined_value("quadratic", true);
        quad_parameters.set_user_defined_value("freq_range",
                                               molresponse_params.freq_range());
        quad_parameters.set_user_defined_value("hfexalg",
                                               molresponse_params.hfexalg());

        if (config.perturbation == "dipole")
        {
          quad_parameters.set_user_defined_value("dipole", true);
          quad_parameters.set_derived_value<size_t>("states", 3);
        }

        auto final_protocol = molresponse_params.protocol().back();
        quad_parameters.set_user_defined_value<vector<double>>(
            "protocol", {final_protocol});

        quad_parameters.set_user_defined_value("dconv",
                                               molresponse_params.dconv());

        return quad_parameters;
      };

      auto quad_parameters = set_hyperpolarizability_parameters();

      // auto calc_params = initialize_calc_params(world,
      // std::string("quad.in"));
      commandlineparser parser;
      auto moldft_archive = moldft_restart.replace_extension().string();

      GroundStateCalculation ground_calculation{world, moldft_archive};
      /*if (world.rank() == 0) {*/
      /*  ground_calculation.print_params();*/
      /*}*/
      Molecule molecule = ground_calculation.molecule();
      quad_parameters.set_ground_state_calculation_data(ground_calculation);
      /*if (world.rank() == 0) {*/
      /*  quad_parameters.print();*/
      /*}*/

      world.gop.fence();
      FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));

      QuadraticResponse quad_calculation{
          world,
          {ground_calculation, molecule, quad_parameters},
          rhs_generator,
      };

      json beta_json_2;

      if (std::filesystem::exists(beta_2_path))
      {
        if (world.rank() == 0)
        {
          std::ifstream ifs(beta_2_path);
          ifs >> beta_json_2;
          // print("beta_json_2: ", beta_json_2.dump(4));

          beta_data = beta_json_2.get<BetaData>();
        }
        world.gop.broadcast_serializable(beta_data.beta_data, 0);
      }

      int startb = 0;
      int startc = 0;

      for (const auto &[index, omegas] : abc_freqs)
      {

        auto [a, b, c] = index;
        auto [omega_a, omega_b, omega_c] = omegas;
        if (beta_data.beta_data.find({a, b, c}) == beta_data.beta_data.end())
        {

          std::array<int, 2> indexs = {b, c};

          auto restartA = root / response_restarts[a];
          auto restartB = root / response_restarts[b];
          auto restartC = root / response_restarts[c];

          std::array<double, 3> omegas{omega_a, omega_b, omega_c};
          std::array<path, 3> restarts{restartA.replace_extension(""),
                                       restartB.replace_extension(""),
                                       restartC.replace_extension("")};

          quad_calculation.set_x_data(world, omegas, restarts);
          auto [beta, beta_directions] =
              quad_calculation.compute_beta_v2(world, omega_b, omega_c);
          world.gop.fence();
          // Need to call copy since beta_ptr belongs to all?

          /*if (world.rank() == 0) {*/
          /*  print("Beta outside of compute_beta_v2: \n", beta);*/
          /*}*/
          std::array<double, 27> beta_i{};
          auto beta_1D = beta.reshape(beta_i.size());

          int i = 0;
          for (auto &beta_ij : beta_i)
          {
            beta_ij = beta_1D[i++];
          }
          world.gop.fence(); // this fence is essential

          if (world.rank() == 0)
          {
            ::print("Beta values for omega_A", " = -(omega_", b, " + omega_", c,
                    ") = -", omega_a, " = (", omega_b, " + ", omega_c, ")");
            {
              for (int i = 0; i < beta_directions.size(); i++)
              {
                std::cout << std::fixed << std::setprecision(5)
                          << "i = " << i + 1 << ", beta[" << beta_directions[i]
                          << "]" << " = " << beta_i[i] << std::endl;
              }
            }
            beta_data.add_data({a, b, c}, {omegas, beta_i});

            std::ofstream out_file(beta_2_path);
            if (out_file.is_open())
            {
              json beta2_json = beta_data;
              out_file << beta2_json.dump(4);
              out_file.close();
            }
          }
        }
        else
        {
          if (world.rank() == 0)
          {
            ::print("Beta data for omega_", b, " + omega_", c,
                    " already exists");
          }
        }
      }

      // add beta data to json
      if (world.rank() == 0)
      {
        json persistent_output = {};
        // read in the output json
        if (std::filesystem::exists(path_manager.get_output_path()))
        {
          std::ifstream ifs(path_manager.get_output_path());
          ifs >> persistent_output;
          ifs.close();
        }
        json out_json = {};
        out_json["beta"] = beta_data.to_json_table();
        // print("beta_data: ", out_json.dump(4));
        persistent_output[name] = out_json;
        std::ofstream ofs(path_manager.get_output_path());
        ofs << persistent_output.dump(4);
        ofs.close();
      }
    }
    catch (Response_Convergence_Error &e)
    {
      if (world.rank() == 0)
      {
        ::print("First order response calculations haven't been run and "
                "can't be run");
        ::print("Quadratic response calculations can't be run");
      }
    }
  }

  void add_beta_i_to_json(nlohmann::ordered_json &beta_i,
                          nlohmann::ordered_json &beta_json)
  {
    // print(beta_json.dump(4));

    for (auto &[key, value] : beta_i.items())
    {
      beta_json[key].insert(beta_json[key].end(), value.begin(), value.end());
    }
  }
};

/**/
/*class WriteResponseVTKOutputStrategy : public CalculationStrategy {*/
/**/
/*  ResponseParameters parameters;*/
/*  std::string op;*/
/*  std::string name;*/
/**/
/* public:*/
/*  explicit WriteResponseVTKOutputStrategy(const ResponseParameters& params) :
 * parameters(params){};*/
/*  json calcPaths(const path& root) override { return {}; }*/
/*  void runCalculation(World& world, const json& paths) override {*/
/**/
/*    auto& moldft_paths = paths["moldft"];*/
/*    auto moldft_restart = moldft_paths["restart"].get<std::string>();*/
/*    moldft_restart =
 * std::string(path(moldft_restart).replace_extension(""));*/
/**/
/*    auto& response_paths = paths["response"];*/
/*    auto& calc_paths = response_paths["calculation"];*/
/**/
/*    auto restart_paths =
 * response_paths["restart"].get<std::vector<std::string>>();*/
/*    auto output_paths =
 * response_paths["output"].get<std::vector<std::string>>();*/
/*    auto alpha_path =
 * response_paths["properties"]["alpha"].get<std::string>();*/
/*    auto freqs = response_paths["frequencies"].get<std::vector<double>>();*/
/*    auto num_freqs = freqs.size();*/
/**/
/*    if (world.rank() == 0) {*/
/*      print("Running VTK output for response calculations");*/
/*      print("Number of frequencies: ", num_freqs);*/
/*      print("Frequencies: ", freqs);*/
/*      print("Output paths: ", output_paths);*/
/*      print("Restart paths: ", restart_paths);*/
/*    }*/
/**/
/*    for (size_t i = 0; i < num_freqs; i++) {*/
/*      auto freq_i = freqs[i];*/
/**/
/*      if (!std::filesystem::exists(calc_paths[i])) {*/
/*        std::filesystem::create_directory(calc_paths[i]);*/
/*      }*/
/*      std::filesystem::current_path(calc_paths[i]);*/
/*      print("current path: ", std::filesystem::current_path());*/
/*      print("calc path: ", calc_paths[i]);*/
/*      print("freq: ", freq_i);*/
/**/
/*      path restart_file_i = restart_paths[i];*/
/*      path response_base_i = output_paths[i];*/
/**/
/*      double last_protocol;*/
/*      bool converged = false;*/
/**/
/*      if (world.rank() == 0) {*/
/*        std::ifstream ifs(response_base_i);*/
/*        json response_base;*/
/*        ifs >> response_base;*/
/*        last_protocol =
 * *response_base["parameters"]["protocol"].get<std::vector<double>>().end();*/
/*        converged = response_base["converged"].get<bool>();*/
/*        print("Thresh of restart data: ", last_protocol);*/
/*        print("Converged: ", converged);*/
/*      }*/
/**/
/*      world.gop.broadcast(converged, 0);*/
/*      world.gop.broadcast(last_protocol, 0);*/
/**/
/*      ResponseParameters r_params = parameters;*/
/*      r_params.set_user_defined_value("omega", freq_i);*/
/*      r_params.set_user_defined_value("archive", moldft_restart);*/
/*      r_params.set_user_defined_value("restart", true);*/
/*      std::string restart_file_string = restart_file_i.filename().stem();*/
/*      r_params.set_user_defined_value("restart_file", restart_file_string);*/
/*      if (converged) {*/
/**/
/*        GroundStateCalculation ground_calculation{world,
 * r_params.archive()};*/
/*        Molecule molecule = ground_calculation.molecule();*/
/*        r_params.set_ground_state_calculation_data(ground_calculation);*/
/*        r_params.set_derived_values(world, molecule);*/
/*        CalcParams calc_params = {ground_calculation, molecule, r_params};*/
/**/
/*        RHS_Generator rhs_generator;*/
/*        if (op == "dipole") {*/
/*          rhs_generator = dipole_generator;*/
/*        } else {*/
/*          rhs_generator = nuclear_generator;*/
/*        }*/
/**/
/*        FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));*/
/*        FrequencyResponse calc(world, calc_params, freq_i, rhs_generator);*/
/*        if (world.rank() == 0) {*/
/*          print("\n\n");*/
/*          print(*/
/*              " MADNESS Time-Dependent Density Functional Theory Response "*/
/*              "Program");*/
/*          print("
 * ----------------------------------------------------------\n");*/
/*          // put the response parameters in a j_molrespone json object*/
/*        }*/
/*        std::string plot_name = "response_" + std::to_string(freq_i);*/
/*        calc.write_vtk(world, 100, parameters.L() / 20, plot_name);*/
/**/
/*        world.gop.fence();*/
/*        // Then we restart from the previous file instead*/
/*      } else {*/
/*      }*/
/*    }*/
/*  }*/
/*};*/
/** /*/
/*class OptimizationCalculationStrategy : public CalculationStrategy {*/
/*  json paths;*/
/*  path output_path;*/
/*  ParameterManager parameters;*/
/*  std::unique_ptr<OptimizationStrategy> optimization_strategy;*/
/**/
/*  std::vector<std::string> available_properties = {};*/
/**/
/* public:*/
/*  OptimizationCalculationStrategy(ParameterManager params,*/
/*                                  std::unique_ptr<OptimizationStrategy>&
 * target,*/
/*                                  std::string calc_name = "optimization")*/
/*      : parameters(std::move(params)),
 * optimization_strategy(std::move(target)),*/
/*        CalculationStrategy(std::move(calc_name),*/
/*                            std::move(std::map<std::string, bool>{})) {}*/
/**/
/*  void setPaths(PathManager& path_manager, const path& root) override {*/
/*    CalculationTemplate opt;*/
/*    auto opt_path = root / name;*/
/**/
/*    opt.calc_paths.push_back(opt_path);*/
/*    opt.restarts.push_back(opt_path / "optimization_restart.00000");*/
/*    opt.outputs["input_molecule"] = {opt_path / "input.mol"};*/
/*    opt.outputs["output_molecule"] = {opt_path / "output.mol"};*/
/*    opt.outputs["properties"] = {opt_path / "output.json"};*/
/**/
/*    path_manager.setPath(name, opt);*/
/*  }*/
/**/
/*  void compute(World& world, PathManager& path_manager) override {*/
/*    auto opt_paths = path_manager.getPath(name);*/
/*    auto calc_path = opt_paths.calc_paths[0];*/
/*    auto restart_path = opt_paths.restarts[0];*/
/*    auto output_path = opt_paths.outputs["properties"][0];*/
/*    auto input_molecule = opt_paths.outputs["input_molecule"][0];*/
/*    auto output_molecule = opt_paths.outputs["output_molecule"][0];*/
/**/
/*    if (world.rank() == 0) {*/
/*      print("Running optimization calculations");*/
/*      print("Calculation path: ", calc_path);*/
/*      print("Restart path: ", restart_path);*/
/*      print("Output path: ", output_path);*/
/*    }*/
/**/
/*    path_manager.createDirectories(name);*/
/**/
/*    Molecule molecule = parameters.get_molecule();*/
/**/
/*    // write the input molecule to the input_molecule file*/
/*    if (world.rank() == 0) {*/
/*      std::ofstream ofs(input_molecule);*/
/*      json j_molecule = molecule.to_json();*/
/*      ofs << j_molecule.dump(4);*/
/*      ofs.close();*/
/*    }*/
/*    auto opt_mol = optimization_strategy->optimize(world, parameters);*/
/*    if (world.rank() == 0) {*/
/*      std::ofstream ofs(output_molecule);*/
/*      json j_molecule = opt_mol.to_json();*/
/*      ofs << j_molecule.dump(4);*/
/*      ofs.close();*/
/*    }*/
/*  }*/
/*};*/

class CalculationDriver
{
private:
  std::string method;
  std::unique_ptr<CompositeCalculationStrategy> strategies;
  PathManager path_manager;
  path root;
  World &world;
  int calculation_number = 0;

public:
  explicit CalculationDriver(World &world,
                             path base_root = std::filesystem::current_path(),
                             std::string method = "moldft")
      : world(world), root(std::move(base_root)), method(std::move(method)),
        path_manager(base_root)
  {
    strategies = std::make_unique<CompositeCalculationStrategy>();
  };
  // Copy constructor with new root
  [[nodiscard]] std::unique_ptr<CalculationDriver> clone() const
  {
    auto clone = std::make_unique<CalculationDriver>(world);
    clone->method = method;
    clone->setRoot(root);
    for (const auto &strategy : strategies->getStrategies())
    {
      clone->addStrategy(strategy->clone());
    }
    return clone;
  }
  void setRoot(const path &new_root)
  {
    root = new_root;
    path_manager = PathManager(new_root);
  }
  [[nodiscard]] path getRoot() const { return root; }
  void setMethod(const std::string &new_method) { method = new_method; }
  [[nodiscard]] std::string getMethod() const { return method; }

  void addStrategy(std::unique_ptr<CalculationStrategy> newStrategy)
  {
    strategies->addStrategy(std::move(newStrategy));
  }
  [[nodiscard]] int getCalcNumber() const { return calculation_number; }
  void
  setStrategies(std::unique_ptr<CompositeCalculationStrategy> newStrategies)
  {
    strategies = std::move(newStrategies);
  }
  path get_output_path() { return path_manager.get_output_path(); }

  void runCalculations(const Tensor<double> &coords)
  {

    strategies->setPaths(path_manager, root);

    auto paths = path_manager.getAllPaths();

    // First step is to look for the paths.json file and read it in if it exists

    if (world.rank() == 0)
    {
      if (std::filesystem::exists("paths.json"))
      {
        std::ifstream ifs("paths.json");
        json read_paths;
        ifs >> read_paths;

        paths.merge_patch(read_paths);
      }
    }

    // TODO: world.gop.broadcast(paths, 0); // broadcast the paths to all the
    // nodes isn't possible because of the json object Might not be a problem
    // but I should look into it

    if (world.rank() == 0)
    {

      // print(paths.dump(4));
      std::ofstream ofs("paths.json");
      ofs << paths.dump(4);
      ofs.close();

      path_manager.createDirectories();
    }
    world.gop.fence();

    strategies->compute(world, path_manager, root, coords);
  }

  // Returns the value at a molecular position
  double value(const Tensor<double> &x)
  {

    // if opt driver setting is set to save then we copy if not we set the
    // current directory to a new name and run test_optimize/optimize
    // test_optimize/optimize/value_0
    // test_optimize/optimize/value_1
    // test_optimize/optimize/value_2
    path root_i;
    if (calculation_number == 0)
    {
      root_i = root / ("value_" + std::to_string(calculation_number));
    }
    else
    {
      root_i =
          root.parent_path() / ("value_" + std::to_string(calculation_number));
    }

    this->setRoot(root_i);
    this->runCalculations(x);
    this->calculation_number++;

    double energy = 0.0;
    if (world.rank() == 0)
    {
      std::ifstream ifs(path_manager.get_output_path());
      json output;
      ifs >> output;
      energy = output[method]["energy"].get<double>();
    }
    world.gop.broadcast(energy, 0);
    return energy;
  }
  void energy_and_gradient(const Molecule &molecule, double &energy,
                           Tensor<double> &gradient)
  {

    energy = this->value(molecule.get_all_coords().flat());
    if (world.rank() == 0)
    {
      std::ifstream ifs(path_manager.get_output_path());
      json output;
      ifs >> output;
      energy = output[method]["energy"].get<double>();
      from_json(output[method]["gradient"], gradient);
    }
    world.gop.broadcast_serializable(gradient, 0);
  }
};

#endif // SRC_APPS_MADQC_CALCMANAGER_HPP_
