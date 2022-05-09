//
// Created by adrianhurtado on 1/1/22.
//
#define CATCH_CONFIG_RUNNER

#include "ExcitedResponse.hpp"
#include <fstream>
#include "response_parameters.h"
#include "FrequencyResponse.hpp"
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


using path = std::filesystem::path;

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

TEST_CASE("Run ground and excited-state") {
    // Set up the run directories
    using namespace madness;

    World &world = World::get_default();

    std::cout.precision(6);

    auto root = std::filesystem::current_path();//="/"+molecule_name;
    auto molecule_path = root;
    molecule_path += "/molecules";

    const std::string molecule_name = "Be";
    const std::string xc = "lda";
    const std::string op = "excited-state";
    // A calculation is defined by a molecule, functional, and operator
    // xc include (hf/lda)
    // operators include (excited-state)
    json response_keyword = {{"molecule", molecule_name}, {"xc", xc}, {"operator", op}};

    auto xc_path = create_xc_path_and_directory(root, xc);
    addResponseKeyWord(response_keyword);
}

TEST_CASE("Create Excited Json") {
    // Set up the run directories
    using namespace madness;

    World &world = World::get_default();

    std::cout.precision(6);

    auto root = std::filesystem::current_path();//="/"+molecule_name;
    auto molecule_path = root;
    molecule_path += "/molecules";

    const std::string xc = "hf";
    const auto ops = vector<std::string>{"excited-state"};
    // A calculation is defined by a molecule, functional, and operator
    // xc include (hf/lda)
    // operators include (excited-state)
    // A calculation is defined by a molecule, functional, and operator
    // xc include (hf/lda)
    for (const auto &op: ops) {
        for (const std::filesystem::directory_entry &mol_path:
             std::filesystem::directory_iterator(molecule_path)) {
            auto molecule_name = mol_path.path().stem();
            json response_keyword = {{"molecule", molecule_name}, {"xc", xc}, {"operator", op}};
            addResponseKeyWord(response_keyword);
        }
    }
}

TEST_CASE("Create Dipole Json") {
    // Set up the run directories
    using namespace madness;

    World &world = World::get_default();

    std::cout.precision(6);

    auto root = std::filesystem::current_path();//="/"+molecule_name;
    auto molecule_path = root;
    molecule_path += "/molecules";

    const std::string xc = "hf";
    const auto ops = vector<std::string>{"dipole"};
    // A calculation is defined by a molecule, functional, and operator
    // xc include (hf/lda)
    // operators include (excited-state)
    // A calculation is defined by a molecule, functional, and operator
    // xc include (hf/lda)
    for (const auto &op: ops) {
        for (const std::filesystem::directory_entry &mol_path:
             std::filesystem::directory_iterator(molecule_path)) {
            auto molecule_name = mol_path.path().stem();
            json response_keyword = {{"molecule", molecule_name}, {"xc", xc}, {"operator", op}};
            addResponseKeyWord(response_keyword);
        }
    }
}

TEST_CASE("response parameters json") {
    // Set up the run directories

    World &world = World::get_default();
    int result = 0;
    world.gop.fence();

    std::cout.precision(6);

    const std::string molecule_name{"Be"};
    const std::string xc{"hf"};
    const std::string op = "excited-state";

    auto schema = runSchema(xc);
    auto m_schema = moldftSchema(molecule_name, xc, schema);
    auto e_schema = excitedSchema(schema, m_schema);

    ResponseParameters original{};
    ResponseParameters params{};
    print("--------------------------default parameters----------------------\n", params.print_to_string());

    std::ifstream ifs(e_schema.rb_json);
    json rb;
    ifs >> rb;

    from_json(rb["response_parameters"], params);

    CHECK(original!=params);
    CHECK(original==original);

}


TEST_CASE("Run if moldft json ==") {
    // Set up the run directories

    World &world = World::get_default();
    int result = 0;
    world.gop.fence();

    std::cout.precision(6);

    const std::string molecule_name{"Be"};
    const std::string xc{"hf"};
    const std::string op = "excited-state";

    auto schema = runSchema(xc);
    auto m_schema = moldftSchema(molecule_name, xc, schema);

    moldft(world, m_schema, false, false, 0);

}


TEST_CASE("Symmetry Tests") {
    // Set up the run directories

    World &world = World::get_default();
    int result = 0;
    world.gop.fence();

    std::cout.precision(6);

    const std::string molecule_name{"Be"};
    const std::string xc{"hf"};
    const std::string op = "excited-state";

    projector_irrep p{};

    print(p.get_pointgroup());


}
