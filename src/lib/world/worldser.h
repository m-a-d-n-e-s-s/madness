#ifndef WORLDSER_H
#define WORLDSER_H

/// \file worldser.h
/// \brief Includes archive headers, and also the archive names into madness namespace 

#include <world/bufar.h>
#include <world/vecar.h>
#include <world/binfsar.h>
#include <world/textfsar.h>

namespace madness {
    using archive::BufferInputArchive;
    using archive::BufferOutputArchive;
    using archive::VectorInputArchive;
    using archive::VectorOutputArchive;
    using archive::BinaryFstreamInputArchive;
    using archive::BinaryFstreamOutputArchive;
    using archive::TextFstreamInputArchive;
    using archive::TextFstreamOutputArchive;
}    
#endif
