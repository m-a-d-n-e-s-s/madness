/*
  pymadness - Python bindings for MADNESS
  World initialization, finalization, and World class bindings.
*/

#include <pybind11/pybind11.h>
#include <madness/world/MADworld.h>
#include <madness/mra/mra.h>

namespace py = pybind11;
using namespace madness;

namespace {

// Static storage for argc/argv needed by MPI_Init
struct InitData {
    static int argc;
    static const char* argv_str;
    static char* argv[2];
    static bool initialized;
};

int InitData::argc = 1;
const char* InitData::argv_str = "pymadness";
char* InitData::argv[2] = {const_cast<char*>(InitData::argv_str), nullptr};
bool InitData::initialized = false;

}  // namespace

void bind_world(py::module_& m) {

    py::class_<World>(m, "World",
        "MADNESS parallel runtime world (MPI communicator wrapper)")
        .def_property_readonly("rank", &World::rank,
            "MPI rank of this process")
        .def_property_readonly("size", &World::size,
            "Total number of MPI processes")
        .def("fence", [](World& w) { w.gop.fence(); },
            "Global synchronization barrier");

    m.def("initialize", [](bool quiet) -> World& {
        if (InitData::initialized) {
            // Already initialized - return the default world
            return *World::world_from_id(0);
        }
        char** argv_ptr = InitData::argv;
        World& world = madness::initialize(InitData::argc, argv_ptr, quiet);
        madness::startup(world, InitData::argc, argv_ptr);
        InitData::initialized = true;
        return world;
    }, py::arg("quiet") = true,
       py::return_value_policy::reference,
       "Initialize the MADNESS runtime. Returns the default World.");

    m.def("finalize", []() {
        if (InitData::initialized) {
            madness::finalize();
            InitData::initialized = false;
        }
    }, "Finalize the MADNESS runtime and release all resources.");

    m.def("is_initialized", []() {
        return InitData::initialized;
    }, "Check whether the MADNESS runtime has been initialized.");
}
