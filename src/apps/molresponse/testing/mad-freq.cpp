//
// Created by adrianhurtado on 1/1/22.
//
#include "ResponseExceptions.hpp"
#include "madness/tensor/tensor_json.hpp"
#include "madness/world/worldmem.h"
#include "response_functions.h"
#include "runners.hpp"

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

using path = std::filesystem::path;

using namespace madness;


auto main(int argc, char *argv[]) -> int {
    World &world = madness::initialize(argc, argv);
    int result = 0;
    world.gop.fence();
    startup(world, argc, argv);

    std::cout.precision(6);

    if (argc != 5) {

        std::cout << "Wrong number of inputs" << std::endl;
        return 1;
    }

    const std::string molecule_name{argv[1]};
    const std::string xc{argv[2]};
    const std::string op{argv[3]};
    const std::string is_high_prec{argv[4]};

    bool high_prec;

    if (is_high_prec == "high") {
        high_prec = true;
    } else {
        high_prec = false;
    }
    try {

        auto schema = runSchema(xc);
        auto m_schema = moldftSchema(molecule_name, xc, schema);
        m_schema.print();
        auto f_schema = frequencySchema(schema, m_schema, op);
        try {

            moldft(world, m_schema, false, false, high_prec);
            runFrequencyTests(world, f_schema, high_prec);
        } catch (MadnessException &madnessException) {
            print(madnessException);
            moldft(world, m_schema, true, false, high_prec);
            runFrequencyTests(world, f_schema, high_prec);
        }

    } catch (const SafeMPI::Exception &e) { print(e); } catch (const madness::MadnessException &e) {
        std::cout << e << std::endl;
    } catch (const madness::TensorException &e) { print(e); } catch (const char *s) {
        print(s);
    } catch (const std::string &s) { print(s); } catch (const std::exception &e) {
        print(e.what());
    } catch (const std::filesystem::filesystem_error &ex) {
        std::cerr << ex.what() << "\n";
    } catch (...) { error("caught unhandled exception"); }

    if (world.rank() == 0) { print("Finished All Frequencies"); }

    return result;

    // print_meminfo(world.rank(), "startup");
    // std::cout.precision(6);
    // print_stats(world);
}
