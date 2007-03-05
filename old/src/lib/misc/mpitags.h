#ifndef MPITAGS_H_
#define MPITAGS_H_

/// \file mpitags.h
/// \brief Constants defining internal MPI tags 

/// MADNESS does not use tags between 0 and 1023 inclusive
/// Non-dynamic MADNESS tags are between 1024 and 2047 inclusive
/// Dynamically generated tags are between 2047 and TAG_UB   

namespace madness {
    static const int BASE_TAG = 1024;
    static const int AM_TAG = BASE_TAG + 1;
    static const int AUTOREF_TAG1 = BASE_TAG + 2;
    static const int AUTOREF_TAG2 = BASE_TAG + 3;
    static const int NORM2SQ_TAG = BASE_TAG + 4;
    static const int COMPRESS_TAG = BASE_TAG + 5;
    static const int RECONSTRUCT_TAG = BASE_TAG + 6;
    static const int REFINE_TAG = BASE_TAG + 7;
    static const int TRUNCATE_TAG1 = BASE_TAG + 8;
    static const int TRUNCATE_TAG2 = BASE_TAG + 9;
    static const int INNER_TAG = BASE_TAG + 10;
    static const int SAVE_TAG = BASE_TAG + 11;
    static const int TASK_GENERIC_DATA_TAG = BASE_TAG + 12;
    static const int AM_LEFT_TAG = BASE_TAG + 13;
    static const int AM_RIGHT_TAG = BASE_TAG + 14; // Must be LEFT_TAG+1
    static const int MPIAR_TAG = BASE_TAG + 15;
    static const int AM_BVALUE_TAG = BASE_TAG + 16;
    
    static const ProcessID INVALID_RANK = -5551212;
    static const int INVALID_TAG = -5551212;

//    static const int _TAG = BASE_TAG + ;
    
}

#endif /*MPITAGS_H_*/
