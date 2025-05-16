#include "core.hpp"

hdf5::file::File create_file(const std::string &name) { 
   return hdf5::file::create(name);
}
