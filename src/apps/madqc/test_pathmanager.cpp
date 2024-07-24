// Copyright 2024 ahurta92"
#define CATCH_CONFIG_RUNNER

#include <filesystem>
#include "madness/external/catch/catch.hpp"
#include "path_manager.hpp"

using path = std::filesystem::path;

int main(int argc, char* argv[]) {
  int result = 0;
  {
    World& world = madness::initialize(argc, argv);

    world.gop.fence();
    startup(world, argc, argv);
    { result = Catch::Session().run(argc, argv); }
    world.gop.fence();

    finalize();
  }
  return result;
}

TEST_CASE("Path Manager", "MoldftPathStrategy") {
  World& world = World::get_default();

  PathManager path_manager;

  // a default path strategy is added
  path_manager.addStrategy(std::make_unique<MoldftPathStrategy>());
  path_manager.addStrategy(std::make_unique<MoldftPathStrategy>("moldft_1"));
  path_manager.addStrategy(std::make_unique<MoldftPathStrategy>("moldft_2"));

  json paths = path_manager.generateCalcPaths("root");
  // output the paths
  std::cout << paths.dump(4) << std::endl;
};
