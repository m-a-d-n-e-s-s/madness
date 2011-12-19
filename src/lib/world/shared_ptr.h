#ifndef MADNESS_SHARED_PTR_H_
#define MADNESS_SHARED_PTR_H_

#include <madness_config.h>

#include <world/boost_checked_delete_bits.h>

// Select header that contains shared_ptr

#if defined(MADNESS_USE_MEMORY)
#  include <memory>
#elif defined(MADNESS_USE_TR1_MEMORY)
#  include <tr1/memory>
#elif defined(MADNESS_USE_BOOST_TR1_MEMORY_HPP)
#  include <boost/tr1/memory.hpp>
#else
#  define MADNESS_HAS_STD_SHARED_PTR 1
#  include <world/shared_ptr_bits.h>
   namespace std {
       using namespace madness::tr1::shptr;
   }
#endif // MEMORY

#if defined(MADNESS_HAS_STD_TR1_SHARED_PTR) && !defined(MADNESS_HAS_STD_SHARED_PTR)
#define MADNESS_HAS_STD_SHARED_PTR 1
// shared_ptr is in std::tr1 but we want it in std namespace
namespace std {
    using ::std::tr1::bad_weak_ptr;
    using ::std::tr1::shared_ptr;
    using ::std::tr1::swap;
    using ::std::tr1::static_pointer_cast;
    using ::std::tr1::dynamic_pointer_cast;
    using ::std::tr1::const_pointer_cast;
    using ::std::tr1::get_deleter;
    using ::std::tr1::weak_ptr;
    using ::std::tr1::enable_shared_from_this;
}

#endif

#endif // MADNESS_SHARED_PTR_H_
