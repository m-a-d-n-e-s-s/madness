#include <apps/molresponse/response_parameters.h>
#include <madchem.h>
#include <madness/chem/CalculationParameters.h>
#include <madness/chem/SCF.h>
#include <madness/chem/molecule.h>
#include <algorithm>
#include <filesystem>
#include <madness/external/nlohmann_json/json.hpp>
#include <utility>
#include "tasks.hpp"
#include "utils.hpp"

using path = std::filesystem::path;
using json = nlohmann::json;
using commandlineparser = madness::commandlineparser;

using namespace madness;




