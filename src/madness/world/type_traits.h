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
#include <type_traits>

/// \file typestuff.h
/// \brief type traits and templates

#include <cstddef>
#include <cstdint>
#include <madness/madness_config.h>
#include <madness/world/function_traits.h>

#include <type_traits>

namespace std {

// C++17 features
#if __cplusplus <= 201402L

// GNU stdlibc++ provides void_t if -gnu++11 or -gnu++14 are given
#if __GNUC__ && defined(__GLIBCXX__) && !__STRICT_ANSI__ && __cplusplus >= 201103L
#define HAVE_VOID_T
#endif

#ifndef HAVE_VOID_T
template <typename... Ts>
struct make_void {
  using type = void;
};
template <typename... Ts>
using void_t = typename make_void<Ts...>::type;
#endif

#endif  // C++17 features

}  // namespace std

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

    /// This defines stuff that is serialiable by bitwise copy N.B. This reports true
    /// for \c T that is an aggregate type (struct or array) that includes pointers.
    template <typename T>
    struct is_trivially_serializable {
      static const bool value = \
        std::is_arithmetic<T>::value || \
        std::is_member_function_pointer<T>::value || \
        std::is_function<T>::value  || \
        std::is_function<typename std::remove_pointer<T>::type>::value || \
//        (std::is_pod<T>::value && !std::is_pointer<T>::value);
        ((std::is_class<T>::value || std::is_array<T>::value) && std::is_trivially_copyable<T>::value);
    };

    /// True for types that are "serialiable" to a std::ostream
    template <typename T, typename = void>
    struct is_ostreammable : std::false_type {};
    template <typename T>
    struct is_ostreammable<T, std::void_t<decltype(std::declval<std::ostream&>() << std::declval<const T&>())>> : std::true_type {};
    /// True for types that are "deserialiable" from an std::istream
    template <typename T, typename = void>
    struct is_istreammable : std::false_type {};
    template <typename T>
    struct is_istreammable<T, std::void_t<decltype(std::declval<std::istream&>() >> std::declval<T&>())>> : std::true_type {};

    template <typename T> constexpr bool is_always_serializable =
    std::is_arithmetic<T>::value || \
    std::is_same<std::nullptr_t, typename std::remove_cv<T>::type>::value || \
    std::is_member_function_pointer<T>::value || \
    std::is_function<T>::value  || \
    std::is_function<typename std::remove_pointer<T>::type>::value;

    template <typename Archive, typename T, typename = void>
    struct is_serializable : std::false_type {};

    // forward declare archives to provide archive-specific overloads
    namespace archive {
    class BinaryFstreamOutputArchive;
    class BinaryFstreamInputArchive;
    class BufferOutputArchive;
    class BufferInputArchive;
    class VectorOutputArchive;
    class VectorInputArchive;
    class TextFstreamOutputArchive;
    class TextFstreamInputArchive;
    class MPIRawOutputArchive;
    class MPIRawInputArchive;
    class MPIOutputArchive;
    class MPIInputArchive;
    }
    template <typename T>
    struct is_serializable<archive::BinaryFstreamOutputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_serializable<archive::BinaryFstreamInputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_serializable<archive::BufferOutputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_serializable<archive::BufferInputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_serializable<archive::VectorOutputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_serializable<archive::VectorInputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_serializable<archive::TextFstreamOutputArchive, T, std::enable_if_t<is_ostreammable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_serializable<archive::TextFstreamInputArchive, T, std::enable_if_t<is_istreammable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_serializable<archive::MPIRawOutputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_serializable<archive::MPIRawInputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_serializable<archive::MPIOutputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_serializable<archive::MPIInputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};

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
