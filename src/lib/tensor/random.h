#ifndef MAD_RANDOM_H
#define MAD_RANDOM_H

/// \file random.h
/// \brief Defines and implements RandomNumber

#include <complex>
typedef std::complex<float> float_complex;
typedef std::complex<double> double_complex;

#include <tensor/mtrand.h>

namespace madness {

    /// A simple random number class that currently wraps the Mersenne twister
    template <class T>
    T RandomNumber() {
        //return (T) std::rand();
        return (T) genrand_int31();
    }
    template <> double RandomNumber<double> ();
    template <> float RandomNumber<float> ();
    template <> double_complex RandomNumber<double_complex> ();
    template <> float_complex RandomNumber<float_complex> ();
}

#endif
