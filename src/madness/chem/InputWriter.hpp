#pragma once
#include <madchem.h>
#include <madness/chem/CalculationParameters.h>
#include <madness/chem/SCF.h>
#include <madness/chem/molecule.h>
#include <apps/molresponse_v2/ResponseParameters.hpp>
#include <madness/external/nlohmann_json/json.hpp>

/// Write out the “molecule” block from a JSON object, in MADNESS .inp format
void write_molecule_json_to_input_file(nlohmann::json const& molecule_json,
                                       std::ostream& output_stream);

/// Write out named blocks (by key) from a JSON object, in MADNESS .inp format
void write_json_to_input_file(nlohmann::json const& input_json,
                              std::vector<std::string> const& keys,
                              std::ostream& output_stream);

/// Top‐level writer of the DFT (.dft) section + molecule block
void write_moldft_input(nlohmann::json const& input_json, std::ostream& out);

/// Top‐level writer of the response section
void write_response_input(nlohmann::json const& input_json, std::ostream& out);
