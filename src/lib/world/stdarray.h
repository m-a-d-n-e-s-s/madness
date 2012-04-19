#ifndef MADNESS_STDARRAY_H
#define MADNESS_STDARRAY_H

#include <madness_config.h>

#if defined(MADNESS_USE_ARRAY)
#  include <array>
#elif defined(MADNESS_USE_TR1_ARRAY)
#  include <tr1/array>
#elif defined(MADNESS_USE_BOOST_TR1_ARRAY_HPP)
#  include <boost/tr1/array.hpp>
#else
#  define MADNESS_HAS_STD_ARRAY 1
#  include <world/stdarray_bits.h>
   namespace std {
      using namespace madness::tr1::array;
   }
#endif

namespace std {
#if defined(MADNESS_HAS_STD_TR1_ARRAY) && !defined(MADNESS_HAS_STD_ARRAY)
#   define MADNESS_HAS_STD_ARRAY 1
    using ::std::tr1::array;
    using ::std::tr1::swap;
    using ::std::tr1::tuple_size;
    using ::std::tr1::tuple_element;
    using ::std::tr1::tuple_size;
    using ::std::tr1::tuple_element;
    using ::std::tr1::get;
#endif
}

#endif
