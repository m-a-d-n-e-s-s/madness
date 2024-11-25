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

#ifndef MADNESS_WORLD_ARCHIVE_H__INCLUDED
#define MADNESS_WORLD_ARCHIVE_H__INCLUDED

/**
 \file archive.h
 \brief Interface templates for the archives (serialization).
 \ingroup serialization
*/

#include <algorithm>
#include <type_traits>
#include <complex>
#include <iostream>
#include <cassert>
#include <cstdio>
#include <cstddef>
#include <cstring>
#include <array>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <optional>
#include <tuple>
#include <madness/config.h>
//#include <madness/world/worldprofile.h>
#include <madness/world/type_traits.h>
#include <madness/world/madness_exception.h>

/// \todo Brief description needed.
#define ARCHIVE_COOKIE "archive"

/// Major version number for archive.
#define ARCHIVE_MAJOR_VERSION 0
/// Minor version number for archive.
#define ARCHIVE_MINOR_VERSION 1

//#define MAD_ARCHIVE_DEBUG_ENABLE

/// Macro for helping debug archive tools.
#ifdef MAD_ARCHIVE_DEBUG_ENABLE
#define MAD_ARCHIVE_DEBUG(s) s
//using std::endl;
#else
#define MAD_ARCHIVE_DEBUG(s)
#endif

namespace madness {

    // Forward declarations
    template <typename T> class Tensor;

    namespace archive {

        // Forward declarations
        template <class>
        class archive_array;
        template <class T>
        inline archive_array<T> wrap(const T*, unsigned int);
        template <class T>
        inline archive_array<unsigned char> wrap_opaque(const T*, unsigned int);
        template <class T>
        inline archive_array<unsigned char> wrap_opaque(const T&);

        /// \addtogroup serialization
        /// @{

        // There are 64 empty slots for user types. Free space for
        // registering user types begins at cookie=128.



        /// The list of type names for use in archives.
        extern const char *archive_type_names[256];

        /// Initializes the type names for the archives.
        // Implemented in archive_type_names.cc
        void archive_initialize_type_names();

        /// Used to enable type checking inside archives.

        /// \tparam T The data type.
        template <typename T>
        struct archive_typeinfo {
            static const unsigned char cookie = 255; ///< Numeric ID for the type; 255 indicates unknown type.
        };

        /// Returns the name of the type, or unknown if not registered.

        /// \tparam T The data type.
        /// \return The name of the type.
        template <typename T>
        const char* get_type_name() {
            return archive_type_names[archive_typeinfo<T>::cookie];
        }

        /// Used to associate a type with a cookie value inside archive.

        /// Makes a specialization of \c archive_typeinfo for type \c T that
        /// specifies the correct cookie value.
        /// \param[in] T The type.
        /// \param[in] cooky The cookie value.
#define ARCHIVE_REGISTER_TYPE(T, cooky) \
        template <> \
        struct archive_typeinfo< T > { \
            static const unsigned char cookie = cooky; \
        }


        /// Used to associate a type and a pointer to the type with a cookie value inside archive.

        /// \param[in] T The type.
        /// \param[in] cooky The cookie value.
#define ARCHIVE_REGISTER_TYPE_AND_PTR(T, cooky) \
        ARCHIVE_REGISTER_TYPE(T, cooky); \
        ARCHIVE_REGISTER_TYPE(T*, cooky+64)


        /// Alias for \c archive_type_names.
#define ATN ::madness::archive::archive_type_names
        /// Alias for \c archive_typeinfo.
#define ATI ::madness::archive::archive_typeinfo


        /// Used to associate names with types.

        /// \param[in] T The type.
#define ARCHIVE_REGISTER_TYPE_NAME(T) \
        if (strcmp( ATN[ATI< T >::cookie], "invalid") ) { \
            std::cout << "archive_register_type_name: slot/cookie already in use! " << #T << " " << ATN[ATI< T >::cookie] << std::endl; \
            MADNESS_EXCEPTION("archive_register_type_name: slot/cookie already in use!", 0); \
         } \
         ATN[ATI< T >::cookie] = #T        


        /// Used to associate names with types and pointers to that type.

        /// \param[in] T The type.
#define ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(T) \
        ARCHIVE_REGISTER_TYPE_NAME(T); \
        ARCHIVE_REGISTER_TYPE_NAME(T*)



#ifndef DOXYGEN_SHOULD_SKIP_THIS
        // **********
        // Register standard types and "common" MADNESS types.
        //
        // doxygen interprets these all as functions, not as calls to a macro.
        // Thus, we force doxygen to skip this block.

        ARCHIVE_REGISTER_TYPE_AND_PTR(unsigned char,0);
        ARCHIVE_REGISTER_TYPE_AND_PTR(unsigned short,1);
        ARCHIVE_REGISTER_TYPE_AND_PTR(unsigned int,2);
        ARCHIVE_REGISTER_TYPE_AND_PTR(unsigned long,3);
        ARCHIVE_REGISTER_TYPE_AND_PTR(unsigned long long,4);
        ARCHIVE_REGISTER_TYPE_AND_PTR(signed char,5);
        ARCHIVE_REGISTER_TYPE_AND_PTR(char,5);	// Needed, but why?
        ARCHIVE_REGISTER_TYPE_AND_PTR(signed short,6);
        ARCHIVE_REGISTER_TYPE_AND_PTR(signed int,7);
        ARCHIVE_REGISTER_TYPE_AND_PTR(signed long,8);
        ARCHIVE_REGISTER_TYPE_AND_PTR(signed long long,9);
        ARCHIVE_REGISTER_TYPE_AND_PTR(bool,10);
        ARCHIVE_REGISTER_TYPE_AND_PTR(float,11);
        ARCHIVE_REGISTER_TYPE_AND_PTR(double,12);
        ARCHIVE_REGISTER_TYPE_AND_PTR(long double,13);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::complex<float>,14);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::complex<double>,15);

        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<char>,20);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<unsigned char>,21);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<short>,22);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<unsigned short>,23);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<int>,24);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<unsigned int>,25);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<long>,26);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<unsigned long>,27);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<bool>,28);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<float>,29);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<double>,30);

        ARCHIVE_REGISTER_TYPE_AND_PTR(std::string,31);

        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor<int>,32);
        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor<long>,33);
        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor<float>,34);
        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor<double>,35);
        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor< std::complex<float> >,36);
        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor< std::complex<double> >,37);
        // **********
#endif

        /// \name function pointer serialization
        /// \note relative function pointers are represented by std::ptrdiff_t , with
        ///       member function pointers represented by std::array<std::ptrdiff_t, N> (with
        ///       N=2 on most common platforms, and a type-dependent constant on some (Microsoft))
        /// @{
        /// \return function pointer to serve as the reference for computing relative pointers
        /// \note the value returned by this function is a pointer to a non-virtual member function,
        ///       this helps on the platforms that use the parity to distinguish non-virtual and virtual pointers (e.g. Itanium ABI)
        std::ptrdiff_t fn_ptr_origin();

        /// \brief converts function or (free or static member) function pointer to the relative function pointer
        /// \param[in] fn a function or function pointer
        template <typename T, typename = std::enable_if_t<std::is_function<T>::value || is_function_pointer<T>::value>>
        std::ptrdiff_t to_rel_fn_ptr(const T& fn) {
          std::ptrdiff_t retval;
          if
#if __cplusplus >= 201703L
            constexpr
#endif
            (std::is_function<T>::value) {
            static_assert(sizeof(std::ptrdiff_t) == sizeof(T*));
            retval = reinterpret_cast<std::ptrdiff_t>(&fn) - fn_ptr_origin();
          }
          else {
#if __cplusplus >= 201703L
            static_assert(sizeof(std::ptrdiff_t) == sizeof(T));
#endif
            retval = reinterpret_cast<std::ptrdiff_t>(fn) - fn_ptr_origin();
          }
          return retval; // To silence stupid Intel 19* compiler
        }

        /// \brief converts nonstatic member function pointer to the relative equivalent
        /// \param[in] fn a member function pointer
        template <typename T, typename = std::enable_if_t<std::is_member_function_pointer<T>::value>>
        auto to_rel_memfn_ptr(const T& fn) {

#if MADNESS_CXX_ABI == MADNESS_CXX_ABI_GenericItanium || MADNESS_CXX_ABI == MADNESS_CXX_ABI_GenericARM
          static_assert(sizeof(T) % sizeof(ptrdiff_t) == 0);
          using result_t = std::array<std::ptrdiff_t, sizeof(T) / sizeof(ptrdiff_t)>;
#elif MADNESS_CXX_ABI == MADNESS_CXX_ABI_Microsoft
          // T will consist of std::ptrdiff_t and up to 3 ints:
          // static_assert((sizeof(T) - sizeof(ptrdiff_t)) % sizeof(int) == 0);
          // using result_t = std::pair<std::ptrdiff_t, std::array<int, (sizeof(T) - sizeof(ptrdiff_t))/sizeof(int)>>;
          //
          //  ... but due to padding it's sufficient to just represent as array of std::ptrdiff_t
          static_assert(sizeof(T) % sizeof(ptrdiff_t) == 0);
          using result_t = std::array<std::ptrdiff_t, sizeof(T) / sizeof(ptrdiff_t)>;
#endif

          result_t result;
          std::memcpy(&result, &fn, sizeof(result_t));
//          typedef union {T abs_ptr; result_t rel_ptr; } cast_t;
//          cast_t caster = {fn};
//          result_t result = caster.rel_ptr;

          // ABI-dependent code:
#if MADNESS_CXX_ABI == MADNESS_CXX_ABI_GenericItanium
          // http://itanium-cxx-abi.github.io/cxx-abi/abi.html#member-pointers
          if (result[0] == 0) { // 0 = null ptr, set adjustment to std::numeric_limits<std::ptrdiff_t>::min() to indicate this
            result[1] = std::numeric_limits<std::ptrdiff_t>::min();
          }
          else if ((result[0] & 1) == 0) { // even pointer = real ptr; odd pointer -> virtual function (does not need translation)
            result[0] -= fn_ptr_origin();
          }
#elif MADNESS_CXX_ABI == MADNESS_CXX_ABI_GenericARM
          // https://developer.arm.com/docs/ihi0041/latest
          const auto even_adj = (result[1] & 1) == 0;
          if (even_adj) {
            if (result[0] == 0) { // {0,even} = null ptr, set ptr to std::numeric_limits<std::ptrdiff_t>::min() to indicate this
              result[0] = std::numeric_limits<std::ptrdiff_t>::min();
            } else { // {nonzero,even} = real ptr
              result[0] -= fn_ptr_origin();
            }
          } else { // odd adjustment -> virtual function (does not need translation)
          }
#elif MADNESS_CXX_ABI == MADNESS_CXX_ABI_Microsoft
          // details are https://www.codeproject.com/Articles/7150/Member-Function-Pointers-and-the-Fastest-Possible and http://www.openrce.org/articles/files/jangrayhood.pdf
          // but I still don't understand them completely, so assume everything is kosher
#warning "Support for serializing member function pointers in Microsoft ABI is not complete"
#endif

          return result;
        }

        /// \brief converts relative free or static member function pointer to the absolute function pointer
        /// \param[in] fn a relative function pointer
        /// \return absolute function pointer
        /// \sa to_rel_fn_ptr
        template <typename T, typename = std::enable_if_t<is_function_pointer_v<T>>>
        T to_abs_fn_ptr(std::ptrdiff_t rel_fn_ptr) {
          static_assert(sizeof(std::ptrdiff_t) == sizeof(T));
          return reinterpret_cast<T>(rel_fn_ptr + fn_ptr_origin());
        }

        /// \brief converts relative (nonstatic) member function pointer to the absolute function pointer
        /// \param[in] fn a relative function pointer
        /// \return absolute function pointer
        /// \sa to_rel_memfn_ptr
        /// \note see notes in to_rel_memfn_ptr re: Microsoft ABI
        template <typename T, std::size_t N, typename = std::enable_if_t<std::is_member_function_pointer<std::remove_reference_t<T>>::value>>
        auto to_abs_memfn_ptr(std::array<std::ptrdiff_t, N> rel_fn_ptr) {
          // ABI-dependent code
#if MADNESS_CXX_ABI == MADNESS_CXX_ABI_GenericItanium
          if (rel_fn_ptr[0] == 0 && rel_fn_ptr[1] == std::numeric_limits<std::ptrdiff_t>::min()) {  // null ptr
            rel_fn_ptr[1] = 0;
          }
          else if ((rel_fn_ptr[0] & 1) == 0) {  // nonvirtual relative ptr is even since fn_ptr_origin is guaranteed to be even also
            rel_fn_ptr[0] += fn_ptr_origin();
          }
#elif MADNESS_CXX_ABI == MADNESS_CXX_ABI_GenericARM
          const auto even_adj = (rel_fn_ptr[1] & 1) == 0;
          if (even_adj) {
            if (rel_fn_ptr[0] == std::numeric_limits<std::ptrdiff_t>::min()) { // null ptr
              rel_fn_ptr[0] = 0;
            } else { // real ptr
              rel_fn_ptr[0] += fn_ptr_origin();
            }
          } else { // odd adjustment -> virtual function (does not need translation)
          }
#elif MADNESS_CXX_ABI == MADNESS_CXX_ABI_Microsoft  // nothing yet done here
#endif
          using result_t = std::remove_reference_t<T>;
          result_t result;
          std::memcpy(&result, &rel_fn_ptr, sizeof(result_t));
          return result;
//          typedef union {result_t abs_ptr; std::array<std::ptrdiff_t, N> rel_ptr; } cast_t;
//          cast_t caster = { .rel_ptr = rel_fn_ptr};
//          return caster.abs_ptr;
        }

        ///@}

        /// Base class for all archive classes.
        class BaseArchive {
        public:
            static constexpr bool is_archive = true; ///< Flag to determine if this object is an archive.
            static constexpr bool is_input_archive = false; ///< Flag to determine if this object is an input archive.
            static constexpr bool is_output_archive = false; ///< Flag to determine if this object is an output archive.
            using is_loading = std::false_type; ///< Type used by Boost.Serialization to determine if this object is an input archive.
            using is_saving = std::false_type; ///< Type used by Boost.Serialization to determine if this object is an output archive.
            static constexpr bool is_parallel_archive = false; ///< Flag to determine if this object is a parallel archive.

            BaseArchive() {
                archive_initialize_type_names();
            }
        }; // class BaseArchive


        /// Base class for input archive classes.
        class BaseInputArchive : public BaseArchive {
        public:
            static constexpr bool is_input_archive = true; ///< Flag to determine if this object is an input archive.
            using is_loading = std::true_type; ///< Type used by Boost.Serialization to determine if this object is an input archive.
        }; // class BaseInputArchive


        /// Base class for output archive classes.
        class BaseOutputArchive : public BaseArchive {
        public:
            static constexpr bool is_output_archive = true; ///< Flag to determine if this object is an output archive.
            using is_saving = std::true_type; ///< Type used by Boost.Serialization to determine if this object is an output archive.
        }; // class BaseOutputArchive


        /// Serialize an array of fundamental stuff.

        /// The function only appears (via \c enable_if) if \c T is
        /// serializable and \c Archive is an output archive.
        /// \tparam Archive The archive type.
        /// \tparam T The type of data in the array.
        /// \param[in] ar The archive.
        /// \param[in] t Pointer to the start of the array.
        /// \param[in] n Number of data items to be serialized.
        template <class Archive, class T>
        typename std::enable_if_t< is_output_archive<Archive>::value &&
                                   is_default_serializable<Archive,T>::value &&
                                   is_function_pointer_v<T>, void>
        default_serialize(const Archive& ar, const T* t, unsigned int n) {
            MAD_ARCHIVE_DEBUG(std::cout << "serialize fund array" << std::endl);
            // transform function pointers to relative pointers
            std::vector<std::ptrdiff_t> t_rel(n);
            std::transform(t, t+n, begin(t_rel), [](auto& fn_ptr) {
              return to_rel_fn_ptr(fn_ptr);
            });
            ar.store(t_rel.data(), n);
        }

        template <class Archive, class T>
        typename std::enable_if_t< is_output_archive<Archive>::value &&
                                   is_default_serializable<Archive,T>::value &&
                                   std::is_member_function_pointer<T>::value, void>
        default_serialize(const Archive& ar, const T* t, unsigned int n) {
            MAD_ARCHIVE_DEBUG(std::cout << "serialize fund array" << std::endl);
            using rel_memfn_ptr_t = decltype(to_rel_memfn_ptr(t[0]));
            static_assert(sizeof(rel_memfn_ptr_t) % sizeof(std::ptrdiff_t) == 0);
            constexpr std::size_t memfn_ptr_width = sizeof(rel_memfn_ptr_t) / sizeof(std::ptrdiff_t);
            std::vector<rel_memfn_ptr_t> t_rel(n);
            std::transform(t, t+n, begin(t_rel), [](auto& fn_ptr) {
              return to_rel_memfn_ptr(fn_ptr);
            });
            ar.store(reinterpret_cast<std::ptrdiff_t*>(t_rel.data()), n * memfn_ptr_width);
        }

        template <class Archive, class T>
        typename std::enable_if_t< is_output_archive<Archive>::value &&
                                      !(is_function_pointer_v<T> || std::is_member_function_pointer<T>::value) &&
                                      is_default_serializable<Archive,T>::value
                                   , void>
        default_serialize(const Archive& ar, const T* t, unsigned int n) {
            MAD_ARCHIVE_DEBUG(std::cout << "serialize fund array" << std::endl);
            ar.store(t, n);
        }


        /// Deserialize an array of fundamental stuff.

        /// The function only appears (via \c enable_if) if \c T is
        /// serializable and \c Archive is an input archive.
        /// \tparam Archive The archive type.
        /// \tparam T The type of data in the array.
        /// \param[in] ar The archive.
        /// \param[in] t Pointer to the start of the array.
        /// \param[in] n Number of data items to be deserialized.
        template <class Archive, class T>
        typename std::enable_if_t< is_input_archive<Archive>::value &&
                                   is_default_serializable<Archive,T>::value &&
                                   is_function_pointer_v<T>, void>
        default_serialize(const Archive& ar, const T* t, unsigned int n) {
            MAD_ARCHIVE_DEBUG(std::cout << "deserialize fund array" << std::endl);
            // transform function pointers to relative offsets
            std::vector<std::ptrdiff_t> t_rel(n);
            ar.load(t_rel.data(),n);
            std::transform(begin(t_rel), end(t_rel), (T*)t, [](auto& fn_ptr_rel) {
              return to_abs_fn_ptr<T>(fn_ptr_rel);
            });
        }

        template <class Archive, class T>
        typename std::enable_if_t< is_input_archive<Archive>::value &&
                                   is_default_serializable<Archive,T>::value &&
                                   std::is_member_function_pointer<T>::value, void>
        default_serialize(const Archive& ar, const T* t, unsigned int n) {
            MAD_ARCHIVE_DEBUG(std::cout << "deserialize fund array" << std::endl);
            using rel_memfn_ptr_t = decltype(to_rel_memfn_ptr(t[0]));
            static_assert(sizeof(rel_memfn_ptr_t) % sizeof(std::ptrdiff_t) == 0);
            constexpr std::size_t memfn_ptr_width = sizeof(rel_memfn_ptr_t) / sizeof(std::ptrdiff_t);
            std::vector<rel_memfn_ptr_t> t_rel(n);
            ar.load(reinterpret_cast<std::ptrdiff_t*>(t_rel.data()), n * memfn_ptr_width);
            std::transform(begin(t_rel), end(t_rel), (T*)t, [](auto& memfn_ptr_rel) {
              return to_abs_memfn_ptr<T>(memfn_ptr_rel);
            });
        }

        template <class Archive, class T>
        typename std::enable_if_t< is_input_archive<Archive>::value && is_default_serializable<Archive,T>::value &&
                                   !(is_function_pointer_v<T> || std::is_member_function_pointer<T>::value), void>
        default_serialize(const Archive& ar, const T* t, unsigned int n) {
            MAD_ARCHIVE_DEBUG(std::cout << "deserialize fund array" << std::endl);
            ar.load((T*) t,n);
        }

        /// Serialize (or deserialize) an array of non-fundamental stuff.

        /// The function only appears (via \c enable_if) if \c T is
        /// not serializable and \c Archive is an archive.
        /// \tparam Archive The archive type.
        /// \tparam T The type of data in the array.
        /// \param[in] ar The archive.
        /// \param[in] t Pointer to the start of the array.
        /// \param[in] n Number of data items to be serialized.
        template <class Archive, class T>
        typename std::enable_if_t< ! is_default_serializable<Archive,T>::value && is_archive<Archive>::value, void>
        serialize(const Archive& ar, const T* t, unsigned int n) {
            MAD_ARCHIVE_DEBUG(std::cout << "(de)serialize non-fund array" << std::endl);
            for (unsigned int i=0; i<n; ++i)
                ar & t[i];
        }


        /// Default implementation of the pre/postamble for type checking.

        /// \tparam Archive The archive class.
        /// \tparam T The type to serialized or to expect upon deserialization.
        template <class Archive, class T>
        struct ArchivePrePostImpl {
            /// Deserialize a cookie and check the type.

            /// \param[in] ar The archive.
            static inline void preamble_load(const Archive& ar) {
                unsigned char ck = archive_typeinfo<T>::cookie;
                unsigned char cookie;
                ar.load(&cookie, 1); // cannot use >>
                if (cookie != ck) {
                    const std::size_t bufsize=255;
                    char msg[bufsize];
                    std::snprintf(msg,bufsize,"InputArchive type mismatch: expected cookie "
                                 "%u (%s) but got %u (%s) instead",
                                 ck, archive_type_names[ck],
                                 cookie,archive_type_names[cookie]);
                    std::cerr << msg << std::endl;
                    MADNESS_EXCEPTION(msg, static_cast<int>(cookie));
                }
                else {
                    MAD_ARCHIVE_DEBUG(std::cout << "read cookie " << archive_type_names[cookie] << std::endl);
                }
            }

            /// Serialize a cookie for type checking.

            /// \param[in] ar The archive.
            static inline void preamble_store(const Archive& ar) {
                unsigned char ck = archive_typeinfo<T>::cookie;
                ar.store(&ck, 1); // cannot use <<
                MAD_ARCHIVE_DEBUG(std::cout << "wrote cookie " << archive_type_names[ck] << std::endl);
            }

            /// By default there is no postamble.
            static inline void postamble_load(const Archive& /*ar*/) {}

            /// By default there is no postamble.
            static inline void postamble_store(const Archive& /*ar*/) {}
        };


        /// Default symmetric serialization of a non-fundamental type that has `serialize` method

        /// \tparam Archive The archive type.
        /// \tparam T The type to symmetrically serialize.
        template <class Archive, class T, typename Enabler>
        struct ArchiveSerializeImpl {
            /// Serializes the type.

            /// \param[in] ar The archive.
            /// \param[in,out] t The data.
            template <typename A = Archive, typename U = T, typename = std::enable_if_t<has_member_serialize_v<U,A>>>
            static inline void serialize(const Archive& ar, T& t) {
              constexpr auto has_serialize = has_member_serialize_v<T,Archive>;
              if constexpr (has_serialize)
                t.serialize(ar);
              else {  // this is never instantiated
                constexpr auto T_has_serialize_method = !has_serialize;
                static_assert(T_has_serialize_method);
              }
            }
        };


        /// Redirect `serialize(ar, t)` to `serialize(ar, &t, 1)` for fundamental types.

        /// The function only appears (due to \c enable_if) if \c T is
        /// serializable and \c Archive is an archive.
        /// \tparam Archive The archive type.
        /// \tparam T The data type.
        /// \param[in] ar The archive.
        /// \param[in] t The data to be serialized.
        template <class Archive, class T>
        inline
        std::enable_if_t< is_default_serializable<Archive,T>::value && is_archive<Archive>::value, void>
        default_serialize(const Archive& ar, const T& t) {
            MAD_ARCHIVE_DEBUG(std::cout << "default_serialize(ar,t) -> default_serialize(ar,&t,1)" << std::endl);
            default_serialize(ar, &t, 1);
        }


        /// Redirect `serialize(ar,t)` to \c ArchiveSerializeImpl for non-fundamental types.

        /// The function only appears (due to \c enable_if) if \c T is not
        /// serializable and \c Archive is an archive.
        /// \tparam Archive The archive type.
        /// \tparam T The data type.
        /// \param[in] ar The archive.
        /// \param[in] t The data to be serialized.
        template <class Archive, class T>
        inline
        std::enable_if_t< (!is_default_serializable<Archive,T>::value && has_nonmember_serialize_v<T, Archive>) && is_archive<Archive>::value, void >
        serialize(const Archive& ar, const T& t) {
            MAD_ARCHIVE_DEBUG(std::cout << "serialize(ar,t) -> ArchiveSerializeImpl" << std::endl);
            ArchiveSerializeImpl<Archive, T>::serialize(ar, const_cast<T&>(t));
        }


        /// Default store of an object via `serialize(ar, t)`.

        /// \tparam Archive The archive type.
        /// \tparam T The data type.
        template <class Archive, class T, typename Enabler>
        struct ArchiveStoreImpl {

            template <typename A = Archive, typename U = T>
            static inline std::enable_if_t<is_output_archive_v<A> &&
                !std::is_function<U>::value &&
                (has_member_serialize_v<U,A> ||
                 has_nonmember_serialize_v<U,A> ||
                 has_freestanding_serialize_v<U, A> ||
                 has_freestanding_default_serialize_v<U, A>),
                                           void>
            store(const A& ar, const U& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "store(ar,t) default" << std::endl);
                if constexpr (has_member_serialize_v<U,A>) {
                  const_cast<U&>(t).serialize(ar);
                }
                else if constexpr (has_nonmember_serialize_v<U,A>) {
                  ArchiveSerializeImpl<A, U>::serialize(ar, const_cast<U&>(t));
                }
                else if constexpr (has_freestanding_serialize_v<U, A>) {
                  serialize(ar, const_cast<U&>(t));
                }
                else if constexpr (has_freestanding_default_serialize_v<U, A>) {
                  default_serialize(ar, const_cast<U&>(t));
                }
                else  // unreachable
                  assert(false);
             }

          /// Store function reference as a function pointer

          /// \param[in] ar The archive.
          /// \param[in] t The data.
          template <typename A = Archive, typename U = T,
              typename = std::enable_if_t<is_output_archive_v<A> &&
                                          std::is_function<U>::value>>
          static inline void store(const Archive& ar, const U& t) {
            MAD_ARCHIVE_DEBUG(std::cout << "store(ar,t) default" << std::endl);
            // convert function to function ptr
            std::add_pointer_t<U> fn_ptr = &t;
            if constexpr (has_freestanding_default_serialize_v<std::add_pointer_t<U>,A>)
              default_serialize(ar,fn_ptr);
            else if constexpr (has_freestanding_serialize_v<std::add_pointer_t<U>,A>)
              serialize(ar,fn_ptr);
            else
              assert(false);
          }

        };


        /// Default load of an object via `serialize(ar, t)`.

        /// \tparam Archive The archive type.
        /// \tparam T The data type.
        template <class Archive, class T, typename Enabler>
        struct ArchiveLoadImpl {
            /// Load an object.

            /// \param[in] ar The archive.
            /// \param[in] t The data.
            template <typename A = Archive,
                      typename U = T,
                      typename = std::enable_if_t<is_input_archive_v<A> &&
                          (has_member_serialize_v<U,A>  ||
                           has_nonmember_serialize_v<U,A> ||
                           has_freestanding_serialize_v<U, A> ||
                           has_freestanding_default_serialize_v<U, A>)>>
            static inline void load(const A& ar, const U& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "load(ar,t) default" << std::endl);
              if constexpr (has_member_serialize_v<U,A>) {
                const_cast<U&>(t).serialize(ar);
              }
              else if constexpr (has_nonmember_serialize_v<U,A>) {
                ArchiveSerializeImpl<A, U>::serialize(ar, const_cast<U&>(t));
              }
              else if constexpr (has_freestanding_serialize_v<U, A>) {
                serialize(ar, const_cast<U&>(t));
              }
              else if constexpr (has_freestanding_default_serialize_v<U, A>) {
                default_serialize(ar, const_cast<U&>(t));
              }
              else  // unreachable
                assert(false);
            }

          /// Load function reference stored as function pointer

          /// \param[in] ar The archive.
          /// \param[in] t The data.
          template <typename A = Archive, typename U = T,
              typename = std::enable_if_t<is_input_archive_v<A> &&
                                          std::is_function<U>::value>>
          static inline void load(const Archive& ar, U& t) {
            MAD_ARCHIVE_DEBUG(std::cout << "load(ar,t) default" << std::endl);
            // convert function ptr to function ref
            std::add_pointer_t<U> fn_ptr;
            if constexpr (has_freestanding_default_serialize_v<std::add_pointer_t<U>,A>)
              default_serialize(ar,fn_ptr);
            else if constexpr (has_freestanding_serialize_v<std::add_pointer_t<U>,A>)
              serialize(ar,fn_ptr);
            else
              assert(false);
            t = *fn_ptr;
          }

        };


        /// Default implementations of \c wrap_store and \c wrap_load.

        /// "Wrapping" refers to the addition of the type's preamble and
        /// postamble around the data to provide runtime type-checking.
        /// \tparam Archive The archive type.
        /// \tparam T The data type.
        template <class Archive, class T, typename Enabler>
        struct ArchiveImpl {
            /// Store an object sandwiched between its preamble and postamble.

            /// \param[in] ar The archive.
            /// \param[in] t The data.
            /// \return The archive.
            template <typename A = Archive, typename = std::enable_if_t<is_output_archive_v<A> && has_nonmember_store_v<T,Archive>>>
            static inline const Archive& wrap_store(const Archive& ar, const T& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_store for default" << std::endl);
                ArchivePrePostImpl<Archive,T>::preamble_store(ar);
                ArchiveStoreImpl<Archive,T>::store(ar,t);
                ArchivePrePostImpl<Archive,T>::postamble_store(ar);
                return ar;
            }

            /// Load an object sandwiched between its preamble and postamble.

            /// \param[in] ar The archive.
            /// \param[in] t The data.
            /// \return The archive.
            template <typename A = Archive, typename = std::enable_if_t<is_input_archive_v<A> && has_nonmember_load_v<T,Archive>>>
            static inline const Archive& wrap_load(const Archive& ar, const T& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_load for default" << std::endl);
                ArchivePrePostImpl<Archive,T>::preamble_load(ar);
                ArchiveLoadImpl<Archive,T>::load(ar,const_cast<T&>(t));  // Loses constness here!
                ArchivePrePostImpl<Archive,T>::postamble_load(ar);
                return ar;
            }
        };


        /// Redirect \c << to \c ArchiveImpl::wrap_store for output archives.

        /// The function only appears (due to \c enable_if) if \c Archive
        /// is an output archive.
        /// \tparam Archive The archive type.
        /// \tparam T The data type.
        /// \param[in] ar The archive.
        /// \param[in] t The data.
        template <class Archive, class T>
        inline
        std::enable_if_t<is_output_archive_v<Archive>, const Archive&>
        operator<<(const Archive& ar, const T& t) {
            // assert that T is serializable
            if constexpr (!is_serializable_v<Archive,T>) {
              constexpr bool T_is_serializable = is_serializable_v<Archive,T>;
              static_assert(T_is_serializable, "T is neither trivially serializable not provided serialization methods");
            }
            //PROFILE_FUNC;
            return ArchiveImpl<Archive,T>::wrap_store(ar,t);
        }

        /// Redirect \c >> to `ArchiveImpl::wrap_load` for input archives.

        /// The function only appears (due to \c enable_if) if \c Archive
        /// is an input archive.
        /// \tparam Archive The archive type.
        /// \tparam T The data type.
        /// \param[in] ar The archive.
        /// \param[in] t The data.
        template <class Archive, class T>
        inline
        std::enable_if_t<is_input_archive_v<Archive>, const Archive&>
        operator>>(const Archive& ar, const T& t) {
            // assert that T is serializable
            if constexpr (!is_serializable_v<Archive,T>) {
              constexpr bool T_is_serializable = is_serializable_v<Archive,T>;
              static_assert(T_is_serializable, "T is neither trivially serializable not provided serialization methods");
            }
            //PROFILE_FUNC;
            return ArchiveImpl<Archive,T>::wrap_load(ar,t);
        }

        /// Redirect \c & to `ArchiveImpl::wrap_store` for output archives.

        /// The function only appears (due to \c enable_if) if \c Archive
        /// is an output archive.
        /// \tparam Archive The archive type.
        /// \tparam T The data type.
        /// \param[in] ar The archive.
        /// \param[in] t The data.
        template <class Archive, class T>
        inline
        std::enable_if_t<is_output_archive_v<Archive>, const Archive&>
        operator&(const Archive& ar, const T& t) {
            // assert that T is serializable
            if constexpr (!is_serializable_v<Archive,T>) {
              constexpr bool T_is_serializable = is_serializable_v<Archive,T>;
              static_assert(T_is_serializable, "T is neither trivially serializable not provided serialization methods");
            }
            //PROFILE_FUNC;
            return ArchiveImpl<Archive,T>::wrap_store(ar,t);
        }

        /// Redirect \c & to `ArchiveImpl::wrap_load` for input archives.

        /// The function only appears (due to \c enable_if) if \c Archive
        /// is an input archive.
        /// \tparam Archive The archive type.
        /// \tparam T The data type.
        /// \param[in] ar The archive.
        /// \param[in] t The data.
        template <class Archive, class T>
        inline
        std::enable_if_t<is_input_archive_v<Archive>, const Archive&>
        operator&(const Archive& ar, const T& t) {
            // assert that T is serializable
            if constexpr (!is_serializable_v<Archive,T>) {
              constexpr bool T_is_serializable = is_serializable_v<Archive,T>;
              static_assert(T_is_serializable, "T is neither trivially serializable not provided serialization methods");
            }
            //PROFILE_FUNC;
            return ArchiveImpl<Archive,T>::wrap_load(ar,t);
        }


        // -----------------------------------------------------------------

        /// Wrapper for an opaque pointer for serialization purposes.

        /// Performs a bitwise copy of the pointer without any remapping.
        /// \tparam T The type of object being pointed to.
        /// \todo Verify this documentation.
        template <class T>
        class archive_ptr {
        public:
            T* ptr; ///< The pointer.

            /// Constructor specifying \c nullptr by default.

            /// \param[in] t The pointer.
            archive_ptr(T* t = nullptr)
                : ptr(t) {}

            /// Dereference the pointer.

            /// \return The dereferenced pointer.
            T& operator*() {
                return *ptr;
            }

            /// Serialize the pointer.

            /// \tparam Archive The archive type.
            /// \param[in] ar The archive.
            template <class Archive>
            void serialize(const Archive& ar) {ar & wrap_opaque(&ptr, 1);}
        };


        /// Wrapper for pointers.

        /// \tparam T The type of object being pointed to.
        /// \param[in] p The pointer.
        /// \return The wrapped pointer.
      	template <class T>
      	inline archive_ptr<T> wrap_ptr(T* p) {
            return archive_ptr<T>(p);
        }

        /// Wrapper for dynamic arrays and pointers.

        /// \tparam T The type of object being pointed to.
        template <class T>
        class archive_array {
        public:
            const T* ptr; ///< The pointer.
            unsigned int n; ///< The number of objects in the array.

            /// Constructor specifying a memory location and size.

            /// \param[in] ptr The pointer.
            /// \param[in] n The number of objects in the array.
            archive_array(const T *ptr, unsigned int n) : ptr(ptr), n(n) {}

            /// Constructor specifying no array and of 0 length.
            archive_array() : ptr(nullptr), n(0) {}
        };


        /// Factory function to wrap a dynamically allocated pointer as a typed \c archive_array.

        /// \tparam T The data type.
        /// \param[in] ptr The pointer.
        /// \param[in] n The number of data elements in the array.
        /// \return The wrapped pointer.
        template <class T>
        inline archive_array<T> wrap(const T* ptr, unsigned int n) {
            return archive_array<T>(ptr,n);
        }


        /// Factory function to wrap a pointer to contiguous data as an opaque (\c uchar) \c archive_array.

        /// \tparam T The data type.
        /// \param[in] ptr The pointer.
        /// \param[in] n The number of data elements in the array.
        /// \return The wrapped pointer, as an opaque \c archive_array.
        template <class T>
        inline archive_array<unsigned char> wrap_opaque(const T* ptr, unsigned int n) {
            return archive_array<unsigned char>((unsigned char*) ptr, n*sizeof(T));
        }

        /// Factory function to wrap a contiguous scalar as an opaque (\c uchar) \c archive_array.

        /// \tparam T The data type.
        /// \param[in] t The data.
        /// \return The wrapped data.
        template <class T>
        inline archive_array<unsigned char> wrap_opaque(const T& t) {
            return archive_array<unsigned char>((unsigned char*) &t,sizeof(t));
        }


        /// Serialize a function pointer.

        /// \tparam Archive The archive type.
        /// \tparam resT The function's return type.
        /// \tparam paramT Parameter pack for the function's arguments.
        template <class Archive, typename resT, typename... paramT>
        struct ArchiveSerializeImpl<Archive, resT(*)(paramT...), std::enable_if_t<!is_default_serializable_v<Archive, resT(*)(paramT...)> >> {
            /// Serialize the function pointer.

            /// \param[in] ar The archive.
            /// \param[in] fn The function pointer.
            template <typename A = Archive, typename = std::enable_if_t<is_output_archive<A>::value>>
            static inline void serialize(const A& ar, resT(*(&fn))(paramT...)) {
                ar &wrap_opaque(to_rel_fn_ptr(fn));
            }

            template <typename A = Archive>
            static inline std::enable_if_t<!is_output_archive<A>::value, void> serialize(const A& ar, resT(*(&fn))(paramT...)) {
                std::ptrdiff_t rel_fn_ptr{};
                ar & wrap_opaque(rel_fn_ptr);
                fn = to_abs_fn_ptr<std::remove_reference_t<decltype(fn)>>(rel_fn_ptr);
            }
        };


        /// Serialize a member function pointer.

        /// \tparam Archive The archive type.
        /// \tparam resT The member function's return type.
        /// \tparam objT The object type.
        /// \tparam paramT Parameter pack for the member function's arguments.
        template <class Archive, typename resT, typename objT, typename... paramT>
        struct ArchiveSerializeImpl<Archive, resT(objT::*)(paramT...),
                                    std::enable_if_t<!is_default_serializable_v<Archive, resT(objT::*)(paramT...)>>> {
            /// Serialize the member function pointer.

            /// \param[in] ar The archive.
            /// \param[in] memfn The member function pointer.
            template <typename A = Archive, typename = std::enable_if_t<is_output_archive<A>::value>>
            static inline void serialize(const A& ar, resT(objT::*(&memfn))(paramT...)) {
                ar &wrap_opaque(to_rel_memfn_ptr(memfn));
            }

            template <typename A = Archive>
            static inline std::enable_if_t<!is_output_archive<A>::value, void> serialize(const A& ar, resT(objT::*(&memfn))(paramT...)) {
                using rel_memfn_ptr_t = decltype(to_rel_memfn_ptr(memfn));
                rel_memfn_ptr_t rel_fn_ptr{};
                ar & wrap_opaque(rel_fn_ptr);
                memfn = to_abs_memfn_ptr<decltype(memfn)>(rel_fn_ptr);
            }
        };

        /// Serialize a const member function pointer.

        /// \tparam Archive The archive type.
        /// \tparam resT The const member function's return type.
        /// \tparam objT The object type.
        /// \tparam paramT Parameter pack for the const member function's arguments.
        template <class Archive, typename resT, typename objT, typename... paramT>
        struct ArchiveSerializeImpl<Archive, resT(objT::*)(paramT...) const,
            std::enable_if_t<!is_default_serializable_v<Archive, resT(objT::*)(paramT...) const>>
            > {
            /// Serialize the const member function pointer.

            /// \param[in] ar The archive.
            /// \param[in] memfn The const member function pointer.
            static inline void serialize(const Archive& ar, resT(objT::*memfn)(paramT...) const) {
                ar & wrap_opaque(memfn);
            }
        };


        /// Partial specialization of \c ArchiveImpl for \c archive_array.

        /// \tparam Archive The archive type.
        /// \tparam T The data type in the \c archive_array.
        template <class Archive, class T>
        struct ArchiveImpl< Archive, archive_array<T> > {
            /// Store the \c archive_array, wrapped by the preamble/postamble.

            /// \param[in] ar The archive.
            /// \param[in] t The \c archive_array.
            /// \return The archive.
            static inline const Archive& wrap_store(const Archive& ar, const archive_array<T>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_store for archive_array" << std::endl);
                ArchivePrePostImpl<Archive,T*>::preamble_store(ar);
                //ar << t.n;
                //ArchivePrePostImpl<Archive,T>::preamble_store(ar);
                if constexpr (has_freestanding_serialize_with_size_v<std::add_pointer_t<T>, Archive>)
                  serialize(ar,(T *) t.ptr,t.n);
                else if constexpr (has_freestanding_default_serialize_with_size_v<std::add_pointer_t<T>, Archive>)
                  default_serialize(ar,(T *) t.ptr,t.n);
                else
                  assert(false);
                //ArchivePrePostImpl<Archive,T>::postamble_store(ar);
                ArchivePrePostImpl<Archive,T*>::postamble_store(ar);
                return ar;
            }

            /// Load the \c archive_array, using the preamble and postamble to perform runtime type-checking.

            /// \param[in] ar The archive.
            /// \param[out] t The \c archive_array.
            /// \return The archive.
            static inline const Archive& wrap_load(const Archive& ar, const archive_array<T>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_load for archive_array" << std::endl);
                ArchivePrePostImpl<Archive,T*>::preamble_load(ar);
                //unsigned int n;
                //ar >> n;
                //if (n != t.n)
                //    MADNESS_EXCEPTION("deserializing archive_array: dimension mismatch", n);
                //ArchivePrePostImpl<Archive,T>::preamble_load(ar);
                if constexpr (has_freestanding_serialize_with_size_v<std::add_pointer_t<T>, Archive>)
                  serialize(ar,(T *) t.ptr,t.n);
                else if constexpr (has_freestanding_default_serialize_with_size_v<std::add_pointer_t<T>, Archive>)
                  default_serialize(ar,(T *) t.ptr,t.n);
                else
                  assert(false);
                //ArchivePrePostImpl<Archive,T>::postamble_load(ar);
                ArchivePrePostImpl<Archive,T*>::postamble_load(ar);
                return ar;
            }
        };


        /// Partial specialization of \c ArchiveImpl for fixed-dimension arrays that redirects to \c archive_array.

        /// \tparam Archive The archive type.
        /// \tparam T The data type.
        /// \tparam n The array size.
        template <class Archive, class T, std::size_t n>
        struct ArchiveImpl<Archive, T[n], std::enable_if_t<!std::is_same_v<T,char> && is_serializable_v<Archive, T>>> {
            /// Store the array, wrapped by the preamble/postamble.

            /// \param[in] ar The archive.
            /// \param[in] t The array.
            /// \return The archive.
            static inline const Archive& wrap_store(const Archive& ar, const T(&t)[n]) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_store for array" << std::endl);
                ar << wrap(&t[0],n);
                return ar;
            }

            /// Load the array, using the preamble and postamble to perform runtime type-checking.

            /// \param[in] ar The archive.
            /// \param[out] t The array.
            /// \return The archive.
            static inline const Archive& wrap_load(const Archive& ar, const T(&t)[n]) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_load for array" << std::endl);
                ar >> wrap(&t[0],n);
                return ar;
            }
        };


        /// Serialize a complex number.

        /// \tparam Archive The archive type.
        /// \tparam T The data type underlying the complex number.
        template <class Archive, typename T>
        struct ArchiveStoreImpl< Archive, std::complex<T>, std::enable_if_t<is_serializable_v<Archive, T>> > {
            /// Store a complex number.

            /// \param[in] ar The archive.
            /// \param[in] c The complex number.
            static inline void store(const Archive& ar, const std::complex<T>& c) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize complex number" << std::endl);
                ar & c.real() & c.imag();
            }
        };


        /// Deserialize a complex number.

        /// \tparam Archive the archive type.
        /// \tparam T The data type underlying the complex number.
        template <class Archive, typename T>
        struct ArchiveLoadImpl< Archive, std::complex<T>, std::enable_if_t<is_serializable_v<Archive, T>> > {
            /// Load a complex number.

            /// \param[in] ar The archive.
            /// \param[out] c The complex number.
            static inline void load(const Archive& ar, std::complex<T>& c) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize complex number" << std::endl);
                T r = 0, i = 0;
                ar & r & i;
                c = std::complex<T>(r,i);
            }
        };

        /// Serialize a \c std::allocator.

        /// This is a no-op.
        /// \tparam Archive the archive type.
        /// \tparam T The data type allocated by the \c allocator.
        template <class Archive, typename T>
        struct ArchiveStoreImpl< Archive, std::allocator<T>, std::enable_if_t<!is_future<T>::value && is_serializable_v<Archive, T>> > {

            /// Storing  a \c std::allocator is a no-op

            /// \param[in] ar The archive.
            /// \param[in] v The \c allocator.
            static inline void store(const Archive& ar, const std::allocator<T>& v) {
            }
        };


        /// Deserialize a \c std::allocator.

        /// This is a no-op.
        /// \tparam Archive the archive type.
        /// \tparam T The data type alllocated by in the \c allocator.
        template <class Archive, typename T>
        struct ArchiveLoadImpl< Archive, std::allocator<T>, std::enable_if_t<!is_future<T>::value && is_serializable_v<Archive, T>> > {

            /// Loading  a \c std::allocator is a no-op

            /// \param[in] ar The archive.
            /// \param[out] v The \c allocator.
            static void load(const Archive& ar, std::allocator<T>& v) {
            }
        };

        /// Serialize a \c std::vector.

        /// \tparam Archive the archive type.
        /// \tparam T The data type stored in the \c vector.
        /// \tparam Alloc The allocator type.
        template <class Archive, typename T, typename Alloc>
        struct ArchiveStoreImpl< Archive, std::vector<T, Alloc>, std::enable_if_t<!is_future<T>::value && is_serializable_v<Archive, T>> > {

            /// Store a \c std::vector of plain data.

            /// \param[in] ar The archive.
            /// \param[in] v The \c vector.
            static inline void store(const Archive& ar, const std::vector<T, Alloc>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize std::vector of plain data" << std::endl);
                if constexpr (!std::allocator_traits<Alloc>::is_always_equal::value) {
                    ar & v.get_allocator();
                }
                ar & v.size();
                ar & wrap(v.data(),v.size());
            }
        };


        /// Deserialize a \c std::vector. Clears and resizes as necessary.

        /// \tparam Archive the archive type.
        /// \tparam T The data type stored in the \c vector.
        /// \tparam Alloc The allocator type.
        template <class Archive, typename T, typename Alloc>
        struct ArchiveLoadImpl< Archive, std::vector<T, Alloc>, std::enable_if_t<!is_future<T>::value && is_serializable_v<Archive, T>> > {

            /// Load a \c std::vector.

            /// Clears and resizes the \c vector as necessary.
            /// \param[in] ar The archive.
            /// \param[out] v The \c vector.
            static void load(const Archive& ar, std::vector<T, Alloc>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize std::vector of plain data" << std::endl);
                if constexpr (!std::allocator_traits<Alloc>::is_always_equal::value) {
                    Alloc allocator;
                    ar & allocator;
                    v = std::vector<T, Alloc>(allocator);
                }
                std::size_t n = 0ul;
                ar & n;
                if (n != v.size()) {
                    v.clear();
                    v.resize(n);
                }
                ar & wrap((T *) v.data(),n);
            }

        };


        /// Serialize a \c std::vector<bool> (as a plain array of bool).

        /// \tparam Archive The archive type.
        /// \tparam Alloc The allocator type.
        template <class Archive, typename Alloc>
        struct ArchiveStoreImpl< Archive, std::vector<bool, Alloc> > {
            /// Store a \c vector<bool>.

            /// \param[in] ar The archive.
            /// \param[in] v The \c vector.
            static inline void store(const Archive& ar, const std::vector<bool, Alloc>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize std::vector<bool>" << std::endl);
                if constexpr (!std::allocator_traits<Alloc>::is_always_equal::value) {
                    ar & v.get_allocator();
                }
                std::size_t n = v.size();
                bool* b = new bool[n];
                for (std::size_t i=0; i<n; ++i) b[i] = v[i];
                ar & n & wrap(b,v.size());
                delete [] b;
            }
        };


        /// Deserialize a std::vector<bool>. Clears and resizes as necessary.

        /// \tparam Archive The archive type.
        /// \tparam Alloc The allocator type.
        template <class Archive, typename Alloc>
        struct ArchiveLoadImpl< Archive, std::vector<bool, Alloc> > {
            /// Load a \c vector<bool>.

            /// Clears and resizes the \c vector as necessary.
            /// \param[in] ar The archive.
            /// \param[out] v The \c vector.
            static void load(const Archive& ar, std::vector<bool, Alloc>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize std::vector<bool>" << std::endl);
                if constexpr (!std::allocator_traits<Alloc>::is_always_equal::value) {
                    Alloc allocator;
                    ar & allocator;
                    v = std::vector<bool, Alloc>(allocator);
                }
                std::size_t n = 0ul;
                ar & n;
                if (n != v.size()) {
                    v.clear();
                    v.resize(n);
                }
                bool* b = new bool[n];
                ar & wrap(b,v.size());
                for (std::size_t i=0; i<n; ++i) v[i] = b[i];
                delete [] b;
            }
        };

        /// Serialize a \c std::array.

        /// \tparam Archive the archive type.
        /// \tparam T The data type stored in the \c std::array.
        /// \tparam N The size of the \c std::array.
        template <class Archive, typename T, std::size_t N>
        struct ArchiveStoreImpl< Archive, std::array<T, N>, std::enable_if_t<is_serializable_v<Archive, T>> > {

            /// Store a \c std::array.

            /// \param[in] ar The archive.
            /// \param[in] v The array object to be serialized.
            static inline void store(const Archive& ar, const std::array<T, N>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize std::array<T," << N << ">, with T plain data" << std::endl);
                ar & v.size();
                ar & wrap(v.data(),v.size());
            }

        };

        /// Deserialize a \c std::array. \c MADNESS_ASSERT 's that the size matches.

        /// \tparam Archive the archive type.
        /// \tparam T The data type stored in the \c std::array.
        /// \tparam N The size of the \c std::array.
        template <class Archive, typename T, std::size_t N>
        struct ArchiveLoadImpl< Archive, std::array<T, N>, std::enable_if_t<is_serializable_v<Archive, T>> > {

            /// Load a \c std::array.

            /// \param[in] ar The archive.
            /// \param[out] v The array to be deserialized.
            static void load(const Archive& ar, std::array<T, N>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize std::array<T," << N << ">, with T plain data" << std::endl);
                std::size_t n = 0ul;
                ar & n;
                MADNESS_ASSERT(n == v.size());
                ar & wrap((T *) v.data(),n);
            }

        };

        /// Serialize a 'std::string'.

        /// \tparam Archive The archive type.
        template <class Archive>
        struct ArchiveStoreImpl< Archive, std::string > {
            /// Store a string.

            /// \param[in] ar The archive.
            /// \param[in] v The string.
            static void store(const Archive& ar, const std::string& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize std::string" << std::endl);
                ar & v.size();
                ar & wrap((const char*) v.data(),v.size());
            }
        };


        /// Deserialize a std::string. Clears and resizes as necessary.

        /// \tparam Archive The archive type.
        template <class Archive>
        struct ArchiveLoadImpl< Archive, std::string > {
            /// Load a string.

            /// Clears and resizes the string as necessary.
            /// \param[in] ar The archive.
            /// \param[out] v The string.
            static void load(const Archive& ar, std::string& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize std::string" << std::endl);
                std::size_t n = 0ul;
                ar & n;
                if (n != v.size()) {
                    v.clear();
                    v.resize(n);
                }
                ar & wrap((char*) v.data(),n);
            }
        };


        /// Serialize (deserialize) an std::pair.

        /// \tparam Archive The archive type.
        /// \tparam T The first data type in the pair.
        /// \tparam Q The second data type in the pair.
        template <class Archive, typename T, typename Q>
        struct ArchiveSerializeImpl< Archive, std::pair<T, Q>, std::enable_if_t<is_serializable_v<Archive, T> && is_serializable_v<Archive, Q>> > {
            /// Serialize the \c pair.

            /// \param[in] ar The archive.
            /// \param[in,out] t The \c pair.
            static inline void serialize(const Archive& ar, std::pair<T,Q>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "(de)serialize std::pair" << std::endl);
                ar & t.first & t.second;
            }
        };

        /// Serialize (deserialize) an std::optional.

        /// \tparam Archive The archive type.
        /// \tparam T The data type stored in the optional object
        template <class Archive, typename T>
        struct ArchiveSerializeImpl< Archive, std::optional<T>, std::enable_if_t<is_serializable_v<Archive, T>> > {
            /// Serialize the \c std::optional.

            /// \param[in] ar The archive.
            /// \param[in,out] t The \c optional.
            static inline void serialize(const Archive& ar, std::optional<T>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "(de)serialize std::optional" << std::endl);
                if constexpr (is_output_archive_v<Archive>) {  // serialize
                    ar & t.has_value();
                    if (t.has_value())
                      ar & t.value();
                } else {
                    bool has_value;
                    ar & has_value;
                    if (has_value) {
                      T value;
                      ar & value;
                      t = std::move(value);
                    }
                }
            }
        };

        namespace {

          template <size_t idx, class Archive, typename... Types>
          struct tuple_serialize_helper;

          template <class Archive, typename... Types>
          struct tuple_serialize_helper<0,Archive,Types...> {
        	  static void exec(const Archive& ar, std::tuple<Types...>& t) {
        		  ar & std::get<0>(t);
        	  }
          };
          template <size_t idx, class Archive, typename... Types>
          struct tuple_serialize_helper {
        	  static void exec(const Archive& ar, std::tuple<Types...>& t) {
        		  ar & std::get<idx>(t);
        		  tuple_serialize_helper<idx-1,Archive,Types...>::exec(ar,t);
        	  }
          };

        };

        /// Serialize (deserialize) a std::tuple

        /// \tparam Archive The archive type.
        /// \tparam Types The tuple payload
        template <class Archive, typename... Types>
        struct ArchiveSerializeImpl< Archive, std::tuple<Types...>, std::enable_if_t<(is_serializable_v<Archive, Types> && ... ) >> {
            /// Serialize the \c std::tuple.

            /// \param[in] ar The archive.
            /// \param[in,out] t The \c tuple.
            static inline void serialize(const Archive& ar, std::tuple<Types...>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "(de)serialize std::tuple" << std::endl);
                constexpr auto size = std::tuple_size<std::tuple<Types...>>::value;
                tuple_serialize_helper<size-1,Archive,Types...>::exec(ar, t);
            }
        };

        /// Serialize an \c std::map.

        /// \tparam Archive The archive type.
        /// \tparam T The map's key type.
        /// \tparam Q The map's data type.
        /// \tparam Compare The map's comparer type.
        /// \tparam Alloc The map's allocator type.
        template <class Archive, typename T, typename Q, typename Compare, typename Alloc>
        struct ArchiveStoreImpl< Archive, std::map<T,Q,Compare,Alloc>, std::enable_if_t<is_serializable_v<Archive, T> && is_serializable_v<Archive, Q>>  > {
            /// Store a \c map.

            /// \param[in] ar The archive.
            /// \param[in] t The \c map.
            static void store(const Archive& ar, const std::map<T,Q,Compare,Alloc>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize std::map" << std::endl);
                if constexpr (!std::allocator_traits<Alloc>::is_always_equal::value) {
                    ar & t.get_allocator();
                }
                ar << t.size();
                for (auto p = t.begin();
                        p != t.end(); ++p) {
                    // Fun and games here since IBM's iterator (const or
                    // otherwise) gives us a const qualified key
                    // (p->first) which buggers up the type matching
                    // unless the user defines pair(T,Q) and pair(const
                    // T,Q) to have cookie (which is tedious).
                    std::pair<T,Q> pp = *p;
                    ar & pp;
                }
            }
        };


        /// Deserialize an \c std::map. The \c map is \em not cleared; duplicate elements are replaced.

        /// \tparam Archive The archive type.
        /// \tparam T The map's key type.
        /// \tparam Q The map's data type.
        /// \tparam Compare The map's comparer type.
        /// \tparam Alloc The map's allocator type.
        template <class Archive, typename T, typename Q, typename Compare, typename Alloc>
        struct ArchiveLoadImpl< Archive, std::map<T,Q,Compare,Alloc>, std::enable_if_t<is_serializable_v<Archive, T> && is_serializable_v<Archive, Q>>  > {
            /// Load a \c map.

            /// The \c map is \em not cleared; duplicate elements are replaced.
            /// \param[in] ar The archive.
            /// \param[out] t The \c map.
            static void load(const Archive& ar, std::map<T,Q,Compare,Alloc>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize std::map" << std::endl);
                if constexpr (!std::allocator_traits<Alloc>::is_always_equal::value) {
                    Alloc allocator;
                    ar & allocator;
                    t = std::map<T,Q,Compare,Alloc>(allocator);
                }
                else
                    t.clear();
                std::size_t n = 0;
                ar & n;
                while (n--) {
                    std::pair<T,Q> p;
                    ar & p;
                    t.emplace(std::move(p.first), std::move(p.second));
                }
            }
        };


        /// Serialize a \c std::set.

        /// \tparam Archive the archive type.
        /// \tparam T The data type stored in the \c set.
        /// \tparam Compare The comparison operator.
        /// \tparam Alloc The allocator.
        template <class Archive, typename T, typename Compare, typename Alloc>
        struct ArchiveStoreImpl< Archive, std::set<T, Compare, Alloc>, std::enable_if_t<!is_future<T>::value && is_serializable_v<Archive, T>> > {

            /// Store a \c std::set.

            /// \param[in] ar The archive.
            /// \param[in] s The \c set.
            static inline void store(const Archive& ar, const std::set<T, Compare, Alloc>& s) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize std::set" << std::endl);
                if constexpr (!std::allocator_traits<Alloc>::is_always_equal::value) {
                    ar & s.get_allocator();
                }
                ar << s.size();
                for (const auto &i : s)
                  ar << i;
            }
        };


        /// Deserialize a \c std::set. Clears and resizes as necessary.

        /// \tparam Archive the archive type.
        /// \tparam T The data type stored in the \c set.
        /// \tparam Compare The comparison operator.
        /// \tparam Alloc The allocator.
        template <class Archive, typename T, typename Compare, typename Alloc>
        struct ArchiveLoadImpl< Archive, std::set<T, Compare, Alloc>, std::enable_if_t<!is_future<T>::value && is_serializable_v<Archive, T>> > {

            /// Load a \c std::set.
            /// \param[in] ar The archive.
            /// \param[out] s The \c set.
            static void load(const Archive& ar, std::set<T, Compare, Alloc>& s) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize std::set" << std::endl);
                if constexpr (!std::allocator_traits<Alloc>::is_always_equal::value) {
                  Alloc allocator;
                  ar & allocator;
                  s = std::set<T, Compare, Alloc>(allocator);
                }
                std::size_t size=0;
                ar >> size;
                s.clear();
                auto hint = s.begin();
                for (std::size_t i = 0; i < size; ++i)
                {
                  typename std::set<T, Compare, Alloc>::key_type key=0;
                  ar >> key;
                  hint = s.emplace_hint(hint, std::move(key));
                }
            }

        };


        /// Serialize a \c std::list.

        /// \tparam Archive the archive type.
        /// \tparam T The data type stored in the \c list.
        /// \tparam Alloc The allocator.
        template <class Archive, typename T, typename Alloc>
        struct ArchiveStoreImpl< Archive, std::list<T,Alloc>, std::enable_if_t<!is_future<T>::value && is_serializable_v<Archive, T>> > {

            /// Store a \c std::list.

            /// \param[in] ar The archive.
            /// \param[in] s The \c list.
            static inline void store(const Archive& ar, const std::list<T, Alloc>& s) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize std::list" << std::endl);
                if constexpr (!std::allocator_traits<Alloc>::is_always_equal::value) {
                  ar & s.get_allocator();
                }
                ar << s.size();
                for (const auto &i : s)
                    ar << i;
            }
        };


        /// Deserialize a \c std::list. Clears and resizes as necessary.

        /// \tparam Archive the archive type.
        /// \tparam T The data type stored in the \c list.
        /// \tparam Alloc The allocator.
        template <class Archive, typename T, typename Alloc>
        struct ArchiveLoadImpl< Archive, std::list<T, Alloc>, std::enable_if_t<!is_future<T>::value && is_serializable_v<Archive, T>> > {

        /// Load a \c std::list.
        /// \param[in] ar The archive.
        /// \param[out] s The \c list.
        static void load(const Archive& ar, std::list<T, Alloc>& s) {
            MAD_ARCHIVE_DEBUG(std::cout << "deserialize std::list" << std::endl);
            if constexpr (!std::allocator_traits<Alloc>::is_always_equal::value) {
                    Alloc allocator;
                    ar & allocator;
                    s = std::list<T, Alloc>(allocator);
            }
            std::size_t size=0;
            ar >> size;
            s.clear();
            for (std::size_t i = 0; i < size; ++i)
            {
                T elem;
                ar >> elem;
                s.emplace_back(std::move(elem));
            }
        }

        };

        /// @}

}  // namespace madness::archive
}  // namespace madness

#endif // MADNESS_WORLD_ARCHIVE_H__INCLUDED
