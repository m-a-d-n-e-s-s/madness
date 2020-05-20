#ifndef MADNESS_PARSEC_INCLUED
#define MADNESS_PARSEC_INCLUED

#include <madness/madness_config.h>

#ifdef HAVE_PARSEC

#include <parsec.h>
#include <parsec/parsec_config.h>
#include <parsec/parsec_internal.h>
#include <parsec/mca/device/device.h>
#include <parsec/execution_stream.h>
#include <parsec/scheduling.h>

#include <iostream>

namespace madness{
  extern "C"{

    extern const dague_function_t madness_function;
    extern dague_handle_t madness_handle;
  }

}

#endif // HAVE_PARSEC

#endif // MADNESS_PARSEC_INCLUED
