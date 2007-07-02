#ifndef MAD_RANDOM_H
#define MAD_RANDOM_H

/// \file random.h
/// \brief Defines and implements RandomNumber

#include <tensor/mtrand.h>

namespace madness {

    /// A simple random number class that currently wraps the Mersenne twister
    template <class T>
    T RandomNumber() {
        //return (T) std::rand();
        return (T) genrand_int31();
    }
    template <> double RandomNumber<double> () {
        //return std::rand()/(RAND_MAX+1.0);
        return genrand_res53();
    }
    template <> float RandomNumber<float> () {
        //return std::rand()/(RAND_MAX+1.0);
        return float(genrand_real2());
    }
    template <> double_complex RandomNumber<double_complex> () {
        return double_complex(RandomNumber<double>(),RandomNumber<double>());
    }
    template <> float_complex RandomNumber<float_complex> () {
        return float_complex(RandomNumber<float>(),RandomNumber<float>());
    }
}

#endif
