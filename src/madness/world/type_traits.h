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
#include <vector>

namespace madness {

   namespace operators {
   class __x {};
   std::ostream& operator<<(std::ostream&, const __x&);
   std::ostream& operator>>(std::ostream&, __x&);
   }  // namespace operators

   // fwd decls
   namespace archive {
     template <typename Archive, typename T, typename Enabler = void>
     struct ArchiveSerializeImpl;

     template <class Archive, class T, typename Enabler = void>
     struct ArchiveLoadImpl;

     template <class Archive, class T, typename Enabler = void>
     struct ArchiveStoreImpl;

     template <class Archive, class T, typename Enabler = void>
     struct ArchiveImpl;
   }

    template <typename> class Future;

    /// test if a type is a future.

    /// \tparam T The type to test.
    template <typename T>
    struct is_future : public std::false_type { };

    template <typename T>
    struct is_future< Future<T> > : public std::true_type { };

    /// maps type \c T to \c Future<T>.

    /// \tparam T The type to have future added.
    template <typename T>
    struct add_future {
        /// Type with \c Future added.
        typedef Future<T> type;
    };

    /// maps \c Future<T> to \c Future<T>.

    /// Specialization of \c add_future<T> that properly forbids the type
    /// \c Future< Future<T> >.
    /// \tparam T The underlying data type.
    template <typename T>
    struct add_future< Future<T> > {
       /// Type with \c Future added.
       typedef Future<T> type;
    };

    /// maps \c Future<T> to \c T.

    /// \tparam T The type to have future removed; in this case, do nothing.
    template <typename T>
    struct remove_future {
        /// Type with \c Future removed.
        typedef T type;
    };

    /// This metafunction maps \c Future<T> to \c T.

    /// \internal Future is a wrapper for T (it acts like an Identity monad), so this
    /// unwraps T. It makes sense that the result should preserve the access traits
    /// of the Future, i.e. const Future<T> should map to const T, etc.

    /// Specialization of \c remove_future for \c Future<T>
    /// \tparam T The type to have future removed.
    template <typename T>
    struct remove_future< Future<T> > {
        /// Type with \c Future removed.
        typedef T type;
    };

    /// Specialization of \c remove_future for \c Future<T>
    /// \tparam T The type to have future removed.
    template <typename T>
    struct remove_future< const Future<T> > {
        /// Type with \c Future removed.
        typedef const T type;
    };

    /// Specialization of \c remove_future for \c Future<T>&
    /// \tparam T The type to have future removed.
    template <typename T>
    struct remove_future< Future<T>& > {
        /// Type with \c Future removed.
        typedef T& type;
    };

    /// Specialization of \c remove_future for \c Future<T>&&
    /// \tparam T The type to have future removed.
    template <typename T>
    struct remove_future< Future<T>&& > {
        /// Type with \c Future removed.
        typedef T&& type;
    };

    /// Specialization of \c remove_future for \c const \c Future<T>&
    /// \tparam T The type to have future removed.
    template <typename T>
    struct remove_future< const Future<T>& > {
        /// Type with \c Future removed.
        typedef const T& type;
    };

    /// Macro to determine type of future (by removing wrapping \c Future template).

    /// \param T The type (possibly with \c Future).
    #define REMFUTURE(T) typename remove_future< T >::type

    /// C++11 version of REMFUTURE
    template <typename T>
    using remove_future_t = typename remove_future< T >::type;

    /// Similar to remove_future , but future_to_ref<Future<T>> evaluates to T& ,whereas
    /// remove_future<Future<T>> evaluates to T .
    /// \tparam T The type to have future removed; in this case, do nothing.
    template <typename T>
    struct future_to_ref {
        typedef T type;
    };
    template <typename T>
    struct future_to_ref<Future<T>> {
        typedef T& type;
    };
    template <typename T>
    struct future_to_ref<Future<T>*> {
        typedef T& type;
    };
    template <typename T>
    struct future_to_ref<Future<T>&> {
        typedef T& type;
    };
    template <typename T>
    struct future_to_ref<const Future<T>&> {
         typedef T& type;
    };
    template <typename T>
    using future_to_ref_t = typename future_to_ref< T >::type;


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

    /// trait for trivial (=bitwise) copyability of T, defaults to std::is_trivially_copyable<T> but can be specialized as needed
    template <typename T>
    struct is_trivially_copyable : std::is_trivially_copyable<T> {};

    template <typename T>
    inline constexpr bool is_trivially_copyable_v = is_trivially_copyable<T>::value;

    /// This defines stuff that is serializable by bitwise copy.
    /// \warning This reports true for \c T that is an aggregate type
    ///          (struct or array) that includes pointers.
    template <typename T>
    struct is_trivially_serializable {
      static const bool value = \
        std::is_arithmetic<T>::value || \
        std::is_function<T>::value  || \
        is_any_function_pointer_v<T> || \
        (std::is_standard_layout<T>::value && std::is_trivial<T>::value && !std::is_pointer<T>::value);
//        ((std::is_class<T>::value || std::is_array<T>::value) && std::is_trivially_copyable<T>::value);
    };

    template <typename T>
    inline constexpr bool is_trivially_serializable_v = is_trivially_serializable<T>::value;

    // namespace hiding implementation details of is_ostreammable ... by ensuring that the detector lives in a different namespace branch than the operators we do not accidentally pick them up
    namespace is_ostreammable_ns {

    template <typename To, typename From> using left_shift = decltype(std::declval<To>() << std::declval<From>());
    template <typename To, typename From> using left_shift_in_ns_madness_operators = decltype(madness::operators::operator<<(std::declval<To>(), std::declval<From>()));

    template <typename T> struct impl : public meta::disjunction<meta::is_detected_exact<std::ostream&, left_shift, std::ostream&, std::add_const_t<std::add_lvalue_reference_t<T>>>,
                                                                 meta::is_detected_exact<std::ostream&, left_shift_in_ns_madness_operators, std::ostream&, std::add_const_t<std::add_lvalue_reference_t<T>>>> {};
    }  // namespace is_ostreammable_ns

    /// True for types that are "serializable" to a std::ostream
    /// \note \c operator<<(std::ostream&,const T&) must be visible via ADL or defined in namespace madness::operators
    template <typename T>
    struct is_ostreammable : public is_ostreammable_ns::impl<T> {};

    /// Shortcut for \c is_ostreammable<T>::value
    template <typename T> constexpr bool is_ostreammable_v = is_ostreammable<T>::value;

    // namespace hiding implementation details of is_istreammable ... by ensuring that the detector lives in a different namespace branch than the operators we do not accidentally pick them up
    namespace is_istreammable_ns {

    template <typename From, typename To> using right_shift = decltype(std::declval<From>() >> std::declval<To>());
    template <typename From, typename To> using right_shift_in_ns_madness_operators = decltype(madness::operators::operator<<(std::declval<From>(), std::declval<To>()));

    template <typename T> struct impl : public meta::disjunction<meta::is_detected_exact<std::istream&, right_shift, std::istream&, std::add_lvalue_reference_t<T>>,
                                                                 meta::is_detected_exact<std::istream&, right_shift_in_ns_madness_operators, std::istream&, std::add_lvalue_reference_t<T>>> {};

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

    /// helps to detect that `T` has a member serialization method that
    /// accepts single argument of type `Archive`
    /// @note use in combination with madness::meta::is_detected_v
    template<typename T, typename Archive>
    using has_member_serialize_t = decltype(std::declval<T&>().serialize(std::declval<Archive&>()));

    /// helps to detect that `T` has a member serialization method that
    /// accepts one argument of type `Archive` and an unsigned version
    /// @note use in combination with madness::meta::is_detected_v
    template<typename T, typename Archive>
    using has_member_serialize_with_version_t = decltype(std::declval<T&>().serialize(std::declval<Archive&>(),0u));

    /// helps to detect that `T` supports nonintrusive symmetric serialization
    /// @note use in combination with madness::meta::is_detected_v
    template<typename T, typename Archive>
    using has_nonmember_serialize_t = decltype(madness::archive::ArchiveSerializeImpl<Archive, T>::serialize(std::declval<Archive&>(), std::declval<T&>()));

    /// helps to detect that `T` supports nonintrusive asymmetric serialization via load
    /// @note use in combination with madness::meta::is_detected_v
    template<typename T, typename Archive>
    using has_nonmember_load_t = decltype(madness::archive::ArchiveLoadImpl<Archive, T>::load(std::declval<Archive&>(), std::declval<T&>()));

    /// helps to detect that `T` supports nonintrusive asymmetric serialization via store
    /// @note use in combination with madness::meta::is_detected_v
    template<typename T, typename Archive>
    using has_nonmember_store_t = decltype(madness::archive::ArchiveStoreImpl<Archive, T>::store(std::declval<Archive&>(), std::declval<T&>()));

    /// helps to detect that `T` supports nonintrusive asymmetric serialization via wrap_load
    /// @note use in combination with madness::meta::is_detected_v
    template<typename T, typename Archive>
    using has_nonmember_wrap_load_t = decltype(madness::archive::ArchiveImpl<Archive, T>::wrap_load(std::declval<Archive&>(), std::declval<T&>()));

    /// helps to detect that `T` supports nonintrusive asymmetric serialization via wrap_store
    /// @note use in combination with madness::meta::is_detected_v
    template<typename T, typename Archive>
    using has_nonmember_wrap_store_t = decltype(madness::archive::ArchiveImpl<Archive, T>::wrap_store(std::declval<Archive&>(), std::declval<T&>()));

    /// helps to detect that `T` supports freestanding `serialize` function
    /// @note use in combination with madness::meta::is_detected_v
    template<typename T, typename Archive>
    using has_freestanding_serialize_t = decltype(serialize(std::declval<Archive&>(), std::declval<T&>()));

    /// helps to detect that `T=U*` supports freestanding `serialize` function
    /// @note use in combination with madness::meta::is_detected_v
    template<typename T, typename Archive, typename = std::enable_if_t<std::is_pointer_v<T>>>
    using has_freestanding_serialize_with_size_t = decltype(serialize(std::declval<Archive&>(), std::declval<T&>(), 1u));

    /// helps to detect that `T` supports freestanding `serialize` function that accepts version
    /// @note use in combination with madness::meta::is_detected_v
    template<typename T, typename Archive, typename = std::enable_if_t<!std::is_pointer_v<T>>>
    using has_freestanding_serialize_with_version_t = decltype(serialize(std::declval<Archive&>(), std::declval<T&>(), 0u));

    /// helps to detect that `T` supports freestanding `default_serialize` function
    /// @note use in combination with madness::meta::is_detected_v
    template<typename T, typename Archive>
    using has_freestanding_default_serialize_t = decltype(default_serialize(std::declval<Archive&>(), std::declval<T&>()));

    /// helps to detect that `T=U*` supports freestanding `default_serialize` function
    /// @note use in combination with madness::meta::is_detected_v
    template<typename T, typename Archive, typename = std::enable_if_t<std::is_pointer_v<T>>>
    using has_freestanding_default_serialize_with_size_t = decltype(default_serialize(std::declval<Archive&>(), std::declval<const T&>(), 1u));

    /// true if this is well-formed:
    /// \code
    ///   // T t; Archive ar;
    ///   t.serialize(ar);
    /// \endcode
    template <typename T, typename Archive>
    inline constexpr bool has_member_serialize_v = madness::meta::is_detected_v<madness::has_member_serialize_t,T,Archive>;

    /// true if this is well-formed:
    /// \code
    ///   // T t; Archive ar;
    ///   t.serialize(ar, 0u);
    /// \endcode
    template <typename T, typename Archive>
    inline constexpr bool has_member_serialize_with_version_v = madness::meta::is_detected_v<madness::has_member_serialize_with_version_t,T,Archive>;

    /// true if this is well-formed:
    /// \code
    ///   // T t; Archive ar;
    ///   madness::archive::ArchiveSerializeImpl<Archive, T>::serialize(ar, t);
    /// \endcode
    template <typename T, typename Archive>
    inline constexpr bool has_nonmember_serialize_v = madness::meta::is_detected_v<madness::has_nonmember_serialize_t,T,Archive>;

    /// true if this is well-formed:
    /// \code
    ///   // T t; Archive ar;
    ///   madness::archive::ArchiveLoadImpl<Archive, T>::load(ar, t);
    /// \endcode
    template <typename T, typename Archive>
    inline constexpr bool has_nonmember_load_v = madness::meta::is_detected_v<madness::has_nonmember_load_t,T,Archive>;

    /// true if this is well-formed:
    /// \code
    ///   // T t; Archive ar;
    ///   madness::archive::ArchiveStoreImpl<Archive, T>::store(ar, t);
    /// \endcode
    template <typename T, typename Archive>
    inline constexpr bool has_nonmember_store_v = madness::meta::is_detected_v<madness::has_nonmember_store_t,T,Archive>;

    template <typename T, typename Archive>
    inline constexpr bool has_nonmember_load_and_store_v = has_nonmember_load_v<T, Archive> && has_nonmember_store_v<T, Archive>;

    /// true if this is well-formed:
    /// \code
    ///   // T t; Archive ar;
    ///   madness::archive::ArchiveImpl<Archive, T>::wrap_load(ar, t);
   /// \endcode
    template <typename T, typename Archive>
    inline constexpr bool has_nonmember_wrap_load_v = madness::meta::is_detected_v<madness::has_nonmember_wrap_load_t,T,Archive>;

    /// true if this is well-formed:
    /// \code
    ///   // T t; Archive ar;
    ///   madness::archive::ArchiveImpl<Archive, T>::wrap_store(ar, t);
    /// \endcode
    template <typename T, typename Archive>
    inline constexpr bool has_nonmember_wrap_store_v = madness::meta::is_detected_v<madness::has_nonmember_wrap_store_t,T,Archive>;

    template <typename T, typename Archive>
    inline constexpr bool has_nonmember_wrap_load_and_store_v = has_nonmember_wrap_load_v<T, Archive> && has_nonmember_wrap_store_v<T, Archive>;

    /// true if this is well-formed:
    /// \code
    ///   // T t; Archive ar;
    ///   serialize(ar, t);
    /// \endcode
    template <typename T, typename Archive>
    inline constexpr bool has_freestanding_serialize_v = madness::meta::is_detected_v<madness::has_freestanding_serialize_t,T,Archive>;

    /// true if this is well-formed:
    /// \code
    ///   // T t; Archive ar;
    ///   serialize(ar, &t, 1u);
    /// \endcode
    template <typename T, typename Archive>
    inline constexpr bool has_freestanding_serialize_with_size_v = madness::meta::is_detected_v<madness::has_freestanding_serialize_with_size_t,T,Archive>;

    /// true if this is well-formed:
    /// \code
    ///   // T t; Archive ar;
    ///   serialize(ar, t, 0u);
    /// \endcode
    template <typename T, typename Archive>
    inline constexpr bool has_freestanding_serialize_with_version_v = madness::meta::is_detected_v<madness::has_freestanding_serialize_with_version_t,T,Archive>;

    /// true if this is well-formed:
    /// \code
    ///   // T t; Archive ar;
    ///   default_serialize(ar, t);
    /// \endcode
    template <typename T, typename Archive>
    inline constexpr bool has_freestanding_default_serialize_v = madness::meta::is_detected_v<madness::has_freestanding_default_serialize_t,T,Archive>;

    /// true if this is well-formed:
    /// \code
    ///   // T t; Archive ar;
    ///   default_serialize(ar, &t, 1u);
    /// \endcode
    template <typename T, typename Archive>
    inline constexpr bool has_freestanding_default_serialize_with_size_v = madness::meta::is_detected_v<madness::has_freestanding_default_serialize_with_size_t,T,Archive>;

    template <typename Archive, typename T, typename Enabler = void>
    struct is_default_serializable_helper : public std::false_type {};

    /// \brief is \c std::true_type if \c T can be serialized to \c Archive
    ///        without specialized \c serialize() method
    ///
    /// For text stream-based \c Archive this is \c std::true_type if \c is_iostreammable<T>::value is true.
    /// For other \c Archive types this is \c std::true_type if \c is_trivially_serializable<T>::value is true.
    /// \tparam Archive an Archive type
    /// \tparam T a type
    template <typename Archive, typename T>
    struct is_default_serializable {
      static constexpr bool value = is_default_serializable_helper<std::remove_cv_t<std::remove_reference_t<Archive>>,std::remove_cv_t<std::remove_reference_t<T>>>::value;
    };

    template <typename Archive, typename T>
    inline constexpr bool is_default_serializable_v = is_default_serializable<Archive, T>::value;

    template <typename Archive, typename T>
    inline constexpr bool is_default_serializable_v<Archive, const T> = is_default_serializable_v<Archive, T>;

    // forward declare archives to provide archive-specific overloads
    namespace archive {
    class BaseArchive;
    class BaseInputArchive;
    class BaseOutputArchive;
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
    class ContainerRecordInputArchive;
    class ContainerRecordOutputArchive;
    template <class localarchiveT>
    class ParallelOutputArchive;
    template <class localarchiveT>
    class ParallelInputArchive;
    template <typename T>
    class archive_array;
    }  // namespace archive

    /// Checks if \c T is an archive type.

    /// If \c T is an archive type, then \c is_archive will be inherited
    /// from \c std::true_type, otherwise it is inherited from
    /// \c std::false_type.
    /// \tparam T The type to check.
    /// \note define for your custom MADNESS archive type
    template <typename T, typename Enabler = void>
    struct is_archive;

    template <typename T>
    using is_archive_defined_t = typename is_archive<std::remove_reference_t<std::remove_cv_t<T>>>::type;

    template <typename T>
    inline constexpr bool is_archive_v = meta::is_detected_v<is_archive_defined_t,T>;


    /// Checks if \c T is an input archive type.

    /// If \c T is an input archive type, then \c is_input_archive will be
    /// inherited from \c std::true_type, otherwise it is inherited from
    /// \c std::false_type.
    /// \tparam T The type to check.
    /// \note define for your custom MADNESS input archive type
    template <typename T, typename Enabler = void>
    struct is_input_archive;

    template <typename T>
    using is_input_archive_defined_t = typename is_input_archive<std::remove_reference_t<std::remove_cv_t<T>>>::type;

    template <typename T>
    inline constexpr bool is_input_archive_v = meta::is_detected_v<is_input_archive_defined_t, T>;

    /// Checks if \c T is an output archive type.

    /// If \c T is an output archive type, then \c is_output_archive will
    /// be inherited from \c std::true_type, otherwise it is inherited from
    /// \c std::false_type.
    /// \tparam T The type to check.
    /// \note define for your custom MADNESS output archive type
    template <typename T, typename Enabler = void>
    struct is_output_archive;

    template <typename T>
    using is_output_archive_defined_t = typename is_output_archive<std::remove_reference_t<std::remove_cv_t<T>>>::type;

    template <typename T>
    inline constexpr bool is_output_archive_v = meta::is_detected_v<is_output_archive_defined_t, T>;

    template <typename T>
    struct is_default_serializable_helper<archive::BinaryFstreamOutputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_default_serializable_helper<archive::BinaryFstreamInputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_default_serializable_helper<archive::BufferOutputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_default_serializable_helper<archive::BufferInputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_default_serializable_helper<archive::VectorOutputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_default_serializable_helper<archive::VectorInputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    // N.B. if type can be printed but can't be read it's not serializable
    // N.N.B. functions and function pointers will be converted to integers, hence will be always serializable
    template <typename T>
    struct is_default_serializable_helper<archive::TextFstreamOutputArchive, T, std::enable_if_t<is_iostreammable_v<T> || std::is_function_v<T> || is_any_function_pointer_v<T>>> : std::true_type {};
    template <typename T>
    struct is_default_serializable_helper<archive::TextFstreamInputArchive, T, std::enable_if_t<is_iostreammable_v<T> || is_any_function_pointer_v<T>>> : std::true_type {};
    template <typename T>
    struct is_default_serializable_helper<archive::MPIRawOutputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_default_serializable_helper<archive::MPIRawInputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_default_serializable_helper<archive::MPIOutputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_default_serializable_helper<archive::MPIInputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_default_serializable_helper<archive::ContainerRecordOutputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T>
    struct is_default_serializable_helper<archive::ContainerRecordInputArchive, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T, class localarchiveT>
    struct is_default_serializable_helper<archive::ParallelOutputArchive<localarchiveT>, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename T, class localarchiveT>
    struct is_default_serializable_helper<archive::ParallelInputArchive<localarchiveT>, T, std::enable_if_t<is_trivially_serializable<T>::value>> : std::true_type {};
    template <typename Archive, typename T>
    struct is_default_serializable_helper<Archive, archive::archive_array<T>, std::enable_if_t<is_default_serializable_helper<Archive,T>::value>> : std::true_type {};

    template <>
    struct is_archive<archive::BinaryFstreamOutputArchive> : std::true_type {};
    template <>
    struct is_archive<archive::BinaryFstreamInputArchive> : std::true_type {};
    template <>
    struct is_archive<archive::BufferOutputArchive> : std::true_type {};
    template <>
    struct is_archive<archive::BufferInputArchive> : std::true_type {};
    template <>
    struct is_archive<archive::VectorOutputArchive> : std::true_type {};
    template <>
    struct is_archive<archive::VectorInputArchive> : std::true_type {};
    template <>
    struct is_archive<archive::TextFstreamOutputArchive> : std::true_type {};
    template <>
    struct is_archive<archive::TextFstreamInputArchive> : std::true_type {};
    template <>
    struct is_archive<archive::MPIRawOutputArchive> : std::true_type {};
    template <>
    struct is_archive<archive::MPIRawInputArchive> : std::true_type {};
    template <>
    struct is_archive<archive::MPIOutputArchive> : std::true_type {};
    template <>
    struct is_archive<archive::MPIInputArchive> : std::true_type {};
    template <>
    struct is_archive<archive::ContainerRecordOutputArchive> : std::true_type {};
    template <>
    struct is_archive<archive::ContainerRecordInputArchive> : std::true_type {};
    template <class localarchiveT>
    struct is_archive<archive::ParallelOutputArchive<localarchiveT> > : std::true_type {};
    template <class localarchiveT>
    struct is_archive<archive::ParallelInputArchive<localarchiveT> > : std::true_type {};

    template <>
    struct is_output_archive<archive::BinaryFstreamOutputArchive> : std::true_type {};
    template <>
    struct is_output_archive<archive::BufferOutputArchive> : std::true_type {};
    template <>
    struct is_output_archive<archive::VectorOutputArchive> : std::true_type {};
    template <>
    struct is_output_archive<archive::TextFstreamOutputArchive> : std::true_type {};
    template <>
    struct is_output_archive<archive::MPIRawOutputArchive> : std::true_type {};
    template <>
    struct is_output_archive<archive::MPIOutputArchive> : std::true_type {};
    template <>
    struct is_output_archive<archive::ContainerRecordOutputArchive> : std::true_type {};
    template <class localarchiveT>
    struct is_output_archive<archive::ParallelOutputArchive<localarchiveT> > : std::true_type {};

    template <>
    struct is_input_archive<archive::BinaryFstreamInputArchive> : std::true_type {};
    template <>
    struct is_input_archive<archive::BufferInputArchive> : std::true_type {};
    template <>
    struct is_input_archive<archive::VectorInputArchive> : std::true_type {};
    template <>
    struct is_input_archive<archive::TextFstreamInputArchive> : std::true_type {};
    template <>
    struct is_input_archive<archive::MPIRawInputArchive> : std::true_type {};
    template <>
    struct is_input_archive<archive::MPIInputArchive> : std::true_type {};
    template <>
    struct is_input_archive<archive::ContainerRecordInputArchive> : std::true_type {};
    template <class localarchiveT>
    struct is_input_archive<archive::ParallelInputArchive<localarchiveT> > : std::true_type {};

    /// Evaluates to true if can serialize an object of type `T` to an object of type `Archive` using user-provided methods
    /// \tparam Archive
    /// \tparam T
    template <typename Archive, typename T>
    inline constexpr bool is_user_serializable_v = is_archive_v<Archive> && (has_member_serialize_v<T, Archive> ||
                                                                    has_nonmember_serialize_v<T, Archive> ||
                                                                    ((has_nonmember_load_v<T, Archive> || has_nonmember_wrap_load_v<T, Archive>) && is_input_archive_v<Archive> && !has_freestanding_default_serialize_v<T, Archive>) ||
                                                                    ((has_nonmember_store_v<T, Archive> || has_nonmember_wrap_store_v<T, Archive>) && is_output_archive_v<Archive> && !has_freestanding_default_serialize_v<T, Archive>));

    template <typename Archive, typename T>
    inline constexpr bool is_user_serializable_v<Archive, const T> = is_user_serializable_v<Archive, T>;

    /// Evaluates to true if can serialize an object of type `T` to an object of type `Archive`,
    /// using either user-provided methods or, if `T` is default-serializable to `Archive`,
    /// using default method for this `Archive`
    /// \tparam Archive
    /// \tparam T
    template <typename Archive, typename T>
    inline constexpr bool is_serializable_v = is_archive_v<Archive> && (is_default_serializable_v<Archive, T> ||
        is_user_serializable_v<Archive,T>);

    template <typename Archive, typename T>
    inline constexpr bool is_serializable_v<Archive, const T> = is_serializable_v<Archive, T>;

    template <typename Archive, typename T>
    struct is_serializable : std::bool_constant<is_serializable_v<Archive,T>> {};

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

    /*
     * Some traits for Functions and alike
     */

    /// next two structs loop over type and dimension
    /// loop over N=1..6 and apply Functor<T,N> to functor_args..., then call the resulting functor with call_args...
    /// returns array of results
    /// Functor must be a template with two parameters: type and integer
    /// functor_args is a tuple of arguments to be forwarded to Functor<T,N>
    /// call_args are arguments to be forwarded to the resulting functor
    /// Example:
    /// \code
    ///   template<typename T, std::size_t N>
    ///   struct MyFunctor {
    ///       MyFunctor(int a, double b) : a_(a), b_(b) {}
    ///       T operator()(const std::string& s) { return T(a_ * N, b_ * N, s); }
    ///       int a_;
    ///       double b_;
    ///   };
    ///   ...
    ///   loop_types<MyFunctor, double, float, double_complex, float_complex>(std::tuple<int,double>(1,2.0),std::string("hello"));
    ///   results is std::array
    /// \endcode
    template<template<typename, std::size_t> class Functor, typename T, std::size_t... Is, typename... FunctorArgs, typename... CallArgs>
    auto loop_N(std::index_sequence<Is...>, std::tuple<FunctorArgs...>&& functor_args, CallArgs&&... call_args)
        -> std::array<decltype(Functor<T, 1>(std::forward<FunctorArgs>(std::get<FunctorArgs>(functor_args))...)(std::forward<CallArgs>(call_args)...)), sizeof...(Is)>
    {
        return { Functor<T, Is + 1>(std::forward<FunctorArgs>(std::get<FunctorArgs>(functor_args))...)(std::forward<CallArgs>(call_args)...)... };
    }

    template<template<typename, std::size_t> class Functor, typename... Ts, typename... FunctorArgs, typename... CallArgs>
    auto loop_types(std::tuple<FunctorArgs...>&& functor_args, CallArgs&&... call_args)
    {
        return std::make_tuple(loop_N<Functor, Ts>(std::make_index_sequence<6>{}, std::move(functor_args), std::forward<CallArgs>(call_args)...)...);
    }

    /// loop over a tuple and apply unary operator op to each element
    template<typename tupleT, typename opT, std::size_t I=0>
    static void unary_tuple_loop(tupleT& tuple, opT& op) {
        if constexpr(I < std::tuple_size_v<tupleT>) {
            auto& element1=std::get<I>(tuple);
            op(element1);
            unary_tuple_loop<tupleT,opT, I+1>(tuple,op);
        }
    }

    /// loop over the tuple elements of both tuples and execute the operation op on each element pair
    template<typename tupleT, typename tupleR, typename opT, std::size_t I=0>
    static void binary_tuple_loop(tupleT& tuple1, tupleR& tuple2, opT& op) {
        if constexpr(I < std::tuple_size_v<tupleT>) {
            auto& element1=std::get<I>(tuple1);
            auto& element2=std::get<I>(tuple2);
            op(element1,element2);
            binary_tuple_loop<tupleT, tupleR, opT, I+1>(tuple1,tuple2,op);
        }
    }

    /// check if objT is a std::vector of Function<T,NDIM>
    /// forward declaration of Function
    /// usage: is_madness_function_vector<objT>::value
    template<typename T, std::size_t NDIM> class Function;
    template<typename>
    struct is_madness_function_vector : std::false_type { };
    template<typename T, std::size_t NDIM>
    struct is_madness_function_vector<std::vector<Function<T, NDIM>>> : std::true_type { };



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
