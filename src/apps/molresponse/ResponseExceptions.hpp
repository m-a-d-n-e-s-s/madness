//
// Created by adrianhurtado on 1/21/22.
//

#ifndef MADNESS_RESPONSEEXCEPTIONS_HPP
#define MADNESS_RESPONSEEXCEPTIONS_HPP
#include "madness_exception.h"
/// Capturing the line/function/filename info is best done with the
/// macros listed below.
/// \param[in] m The error message.
/// \param[in] a String describing the exception.
/// \param[in] v Value associated with the exception.
/// \param[in] l Line number where the exception occurred.
/// \param[in] fn Function where the exception occurred.
/// \param[in] f File where the exception occurred.
using MadException = madness::MadnessException;

class Input_Error : public MadException {
public:
    explicit Input_Error()
        : MadException("input file not found", nullptr, 25, __LINE__, __FUNCTION__, __FILE__) {}
};
class Response_Input_Error : public MadException {
public:
    explicit Response_Input_Error()
        : MadException("Response input not correct", nullptr, 25, __LINE__, __FUNCTION__,
                       __FILE__) {}
};

class Response_Convergence_Error : public MadException {
public:
    explicit Response_Convergence_Error()
        : MadException("Previous frequency response calculation did not converge ", nullptr, 25,
                       __LINE__, __FUNCTION__, __FILE__) {}
};

class Quad_Response_Error : public MadException {
public:
    explicit Quad_Response_Error()
            : MadException("Previous frequency response calculation did not converge ", nullptr, 25,
                           __LINE__, __FUNCTION__, __FILE__) {}
};

#endif// MADNESS_RESPONSEEXCEPTIONS_HPP
