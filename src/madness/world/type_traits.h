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

/// \file typestuff.h
/// \brief type traits and templates

/*
 * N.B. this must be pure c++, usable without any context other than
 *      the current compiler + library + C++ standard.
 *      DO NOT include non-standard headers here!
 */

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <iosfwd>
#include <madness/world/meta.h>

namespace madness {

   namespace operators {
   class __x {};
   std::ostream& operator<<(std::ostream&, const __x&);
   std::ostream& operator>>(std::ostream&, __x&);
   }

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

    /// is true type if \p T is a pointer to a free function
    template <typename T, typename Enabler = void> struct is_function_pointer : public std::false_type {};
    template <typename T> struct is_function_pointer<T, std::enable_if_t<std::is_function<typename std::remove_pointer<T>::type>::value>> : public std::true_type {};
    template <typename T> constexpr bool is_function_pointer_v = is_function_pointer<T>::value;

    // use std::is_member_function_pointer<T> if looking for is_member_function_pointer

    /// is true type if \p T is a pointer to free or member function
    template <typename T, typename Enabler = void> struct is_any_function_pointer : public std::false_type {};
    template <typename T> struct is_any_function_pointer<T, std::enable_if_t<std::is_member_function_pointer<T>::value || is_function_pointer_v<T>>> : public std::true_type {};
    template <typename T> constexpr bool is_any_function_pointer_v = is_any_function_pointer<T>::value;

    /// This defines stuff that is serialiable by bitwise copy.
    /// \warning This reports true for \c T that is an aggregate type
    ///          (struct or array) that includes pointers.
    template <typename T>
    struct is_trivially_serializable {
      static const bool value = \
        std::is_arithmetic<T>::value || \
        std::is_function<T>::value  || \
        is_any_function_pointer_v<T> || \
        (std::is_pod<T>::value && !std::is_pointer<T>::value);
//        ((std::is_class<T>::value || std::is_array<T>::value) && std::is_trivially_copyable<T>::value);
    };

    // namespace hiding implementation details of is_ostreammable ... by ensuring that the detector lives in a different namespace branch than the operators we do not accidentally pick them up
    namespace is_ostreammable_ns {

    template <typename To, typename From> using left_shift = decltype(std::declval<To>() << std::declval<From>());
    template <typename To, typename From> using left_shift_in_ns_madness_operators = decltype(madness::operators::operator<<(std::declval<To>(), std::declval<From>()));

    template <typename T> struct impl : public meta::disjunction<meta::is_detected_exact<std::ostream&, left_shift, std::ostream&, const T&>,
                                                                 meta::is_detected_exact<std::ostream&, left_shift_in_ns_madness_operators, std::ostream&, const T&>> {};
    }  // namespace is_ostreammable_ns

    /// True for types that are "serialiable" to a std::ostream
    /// \note \c operator<<(std::ostream&,const T&) must be visible via ADL or defined in namespace madness::operators
    template <typename T>
    struct is_ostreammable : public is_ostreammable_ns::impl<T> {};

    /// Shortcut for \c is_ostreammable<T>::value
    template <typename T> constexpr bool is_ostreammable_v = is_ostreammable<T>::value;

    // namespace hiding implementation details of is_istreammable ... by ensuring that the detector lives in a different namespace branch than the operators we do not accidentally pick them up
    namespace is_istreammable_ns {

    template <typename From, typename To> using right_shift = decltype(std::declval<From>() >> std::declval<To>());
    template <typename From, typename To> using right_shift_in_ns_madness_operators = decltype(madness::operators::operator<<(std::declval<From>(), std::declval<To>()));

    template <typename T> struct impl : public meta::disjunction<meta::is_detected_exact<std::istream&, right_shift, std::istream&, T&>,
                                                                 meta::is_detected_exact<std::istream&, right_shift_in_ns_madness_operators, std::istream&, T&>> {};

    }  // namespace is_istreammable_ns

    /// True for types that are "deserialiable" from an std::istream
    /// \note \c operator>>(std::ostream&,T&) must be visible via ADL or defined in namespace madness::operators
    template <typename T>
    struct is_istreammable : public is_istreammable_ns::impl<T> {};

    /// Shortcut for \c is_istreammable<T>::value
    template <typename T> constexpr bool is_istreammable_v = is_istreammable<T>::value;

    /// providing automatic support for serializing to/from std streams requires bidirectional streammability
    template <typename T> constexpr bool is_iostreammable_v = is_istreammable_v<T> && is_ostreammable_v<T>;

    template <typename T> constexpr bool is_always_serializable =
    std::is_arithmetic<T>::value || \
    std::is_same<std::nullptr_t, typename std::remove_cv<T>::type>::value || \
    is_any_function_pointer_v<T> || \
    std::is_function<T>::value;

    /// \brief is \c std::true_type if \c T can be serialized to \c Archive
    ///        without specialized \c serialize() method
    ///
    /// For text stream-based \c Archive this is \c std::true_type if \c is_iostreammable<T>::value is true.
    /// For other \c Archive types this is \c std::true_type if \c is_trivially_serializable<T>::value is true.
    /// \tparam Archive an Archive type
    /// \tparam T a type
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
    struct is_serializable<archive::TextFstreamOutputArchive, T, std::enable_if_t<is_iostreammable_v<T>>> : std::true_type {};
    template <typename T>
    struct is_serializable<archive::TextFstreamInputArchive, T, std::enable_if_t<is_iostreammable_v<T>>> : std::true_type {};
    template <typename T>
    struct is_serializable<archive::MPIRawOutputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_serializable<archive::MPIRawInputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_serializable<archive::MPIOutputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_serializable<archive::MPIInputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};

    /// \brief This trait types tests if \c Archive is a text archive
    /// \tparam Archive an archive type
    /// \note much be specialized for each archive
    template <typename Archive, typename Enabler = void>
    struct is_text_archive : std::false_type {};

    template <>
    struct is_text_archive<archive::TextFstreamOutputArchive> : std::true_type {};
    template <>
    struct is_text_archive<archive::TextFstreamInputArchive> : std::true_type {};

    /// \brief \c is_text_archive_v<A> is a shorthand for \c is_text_archive<A>::value
    /// \tparam Archive an archive type
    template <typename Archive>
    constexpr const bool is_text_archive_v = is_text_archive<Archive>::value;

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
