/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

#ifndef MADNESS_WORLD_TYPE_TRAITS_H__INCLUDED
#define MADNESS_WORLD_TYPE_TRAITS_H__INCLUDED

#include <madness/madness_config.h>
#  include <type_traits>

/// \file typestuff.h
/// \brief Grossly simplified Boost-like type traits and templates

#include <cstddef>
#include <stdint.h>
#include <madness/madness_config.h>
#include <madness/world/function_traits.h>

namespace madness {

    template <typename> class Future;
    template <typename> struct add_future;
    template <typename> struct remove_future;

    // Remove Future, const, volatile, and reference qualifiers from the type
    template <typename T>
    struct remove_fcvr {
        typedef typename remove_future<typename std::remove_cv<
                   typename std::remove_reference<T>::type>::type>::type type;
    };
    template <typename T>
    using remove_fcvr_t = typename remove_fcvr<T>::type;

    /// This defines stuff that is serialiable by default rules ... basically anything contiguous
    template <typename T>
    struct is_serializable {
        static const bool value = std::is_fundamental<T>::value || std::is_member_function_pointer<T>::value || std::is_function<T>::value  || std::is_function<typename std::remove_pointer<T>::type>::value || std::is_pod<T>::value;
    };


    /* Macros to make some of this stuff more readable */

    /**
       \def REMCONST(TYPE)
       \brief Macro to make remove_const<T> easier to use

       \def MEMFUN_RETURNT(TYPE)
       \brief Macro to make member function type traits easier to use
    */

#define REMCONST(TYPE)  typename std::remove_const< TYPE >::type
#define MEMFUN_RETURNT(MEMFUN) typename madness::detail::memfunc_traits< MEMFUN >::result_type

} // namespace madness

#endif // MADNESS_WORLD_TYPE_TRAITS_H__INCLUDED
