#ifndef MAD_TYPES_H
#define MAD_TYPES_H

/// \file mad_types.h
/// \brief Defines types widely used by madness (translation, level, processid, ...)

/// Translation in 1d ... more than 31 levels of refinement will require wide integers
typedef unsigned long Translation;

/// Level
typedef int Level;

/// ID (rank) of a parallel process
typedef int ProcessID;

/// Cost (compute cost) of a node (can be changed but set to integer for now)
typedef unsigned long Cost;

#endif
