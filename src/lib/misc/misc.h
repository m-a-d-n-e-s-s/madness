#ifndef MISC_H
#define MISC_H

/// \file misc.h
/// \brief Header to declare stuff which has not yet found a home

#include "communicator.h"
namespace madness {
    double walltime();
    double cputime();
    unsigned long checksum_file(const char* filename);
    void redirectio(const Communicator& comm);
}

#endif



