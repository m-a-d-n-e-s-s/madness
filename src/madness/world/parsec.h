#ifndef MADNESS_PARSEC_INCLUED
#define MADNESS_PARSEC_INCLUED

#include <madness/madness_config.h>

#ifdef HAVE_PARSEC

#include <dague.h>
#include <dague/dague_internal.h>
#include <dague/devices/device.h>
#include <dague/execution_unit.h>
#include <dague/scheduling.h>  
#include <iostream>

namespace madness{
  extern "C"{

    extern const dague_function_t madness_function;
    extern dague_handle_t madness_handle;
  }

}

#endif // HAVE_PARSEC

#endif // MADNESS_PARSEC_INCLUED
