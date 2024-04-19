//
// Created by adrianhurtado on 1/1/22.
#include "coordinator.hpp"

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && defined(HAVE_UNISTD_H)

#include <sys/stat.h>
#include <unistd.h>

static inline int file_exists(const char *input_name) {
    struct stat buffer {};
    size_t rc = stat(input_name, &buffer);
    return (rc == 0);
}

#endif

using path = std::filesystem::path;

using namespace madness;


auto main(int argc, char *argv[]) -> int {

    madness::initialize(argc, argv);
    std::cout.precision(6);

    {
        World world(SafeMPI::COMM_WORLD);
        startup(world, argc, argv, true);

        try {


            path input_json("resources/inputs/freq_input.json");
            path mol_input("resources/molecules/H2.mol");

            ParameterManager params(world, input_json, mol_input);
            auto response_manager = ResponseCalcManager(world, params);

            if(world.rank() == 0) {
                print("Running MOLDFT");
                print("Calc Info Path: ", response_manager.calc_info_json_path);
                print("Moldft Path: ", response_manager.moldft_path);
            }

            if (std::filesystem::exists(response_manager.calc_info_json_path) &&
                std::filesystem::exists(response_manager.moldft_path)) {
                response_manager.run_molresponse(world);
            } else {
                response_manager.run_moldft(world,true);
                world.gop.fence();
                response_manager.run_molresponse(world);
                world.gop.fence();
            }
        } catch (const SafeMPI::Exception &e) {
            print(e.what());
            error("caught an MPI exception");
        } catch (const madness::MadnessException &e) {
            print(e);
            error("caught a MADNESS exception");
        } catch (const madness::TensorException &e) {
            print(e.what());
            error("caught a Tensor exception");
        } catch (const nlohmann::detail::exception &e) {
            print(e.what());
            error("Caught JSON exception");
        } catch (const std::filesystem::filesystem_error &ex) {
            std::cerr << ex.what() << "\n";
        } catch (const std::exception &e) {
            print(e.what());
            error("caught an STL exception");
        } catch (...) { error("caught unhandled exception"); }
        // Nearly all memory will be freed at this point
        print_stats(world);
        if (world.rank() == 0) { print("Finished All Frequencies"); }
    }
    finalize();
    return 0;
}
