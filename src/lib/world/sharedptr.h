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


  $Id$
*/


#ifndef MADNESS_WORLD_SHAREDPTR_H__INCLUDED
#define MADNESS_WORLD_SHAREDPTR_H__INCLUDED

/// \file sharedptr.h
/// \brief Includes TR1 shared_ptr. If shared_ptr is in std::tr1 namespace, it
/// is imported into the std namespace. It also includes make_shared and
/// allocate_shared helper functions which are a part of the current C++0x
/// draft.

#include <madness_config.h>

// Select header that contains shared_ptr
#if defined(MADNESS_USE_MEMORY)
#include <memory>
#elif defined(MADNESS_USE_TR1_MEMORY)
#include <tr1/memory>
#elif defined(MADNESS_USE_BOOST_TR1_MEMORY_HPP)
#include <boost/tr1/memory.hpp>
#else
#error No acceptable memory include directive was found.
#endif // MEMORY

#if defined(MADNESS_HAS_STD_TR1_SHARED_PTR) && !defined(MADNESS_HAS_STD_SHARED_PTR)
#define MADNESS_HAS_STD_SHARED_PTR 1
// shard_ptr is in std::tr1 but we want it in std namespace
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

#endif // defined(MADNESS_HAS_STD_TR1_SHARED_PTR) && !defined(MADNESS_HAS_STD_SHARED_PTR)

namespace madness {
    namespace detail {

        // These checked delete and deleters are copied from Boost.
        // They ensure that compilers issue warnings if T is an incomplete type.

        /// Checked pointer delete function

        /// This function ensures that the pointer is a complete type.
        /// \tparam T The pointer type (must be a complete type).
        /// \param p The pointer to be deleted.
        template<typename T>
        inline void checked_delete(T* p) {
            // intentionally complex - simplification causes regressions
            typedef char type_must_be_complete[ sizeof(T)? 1: -1 ];
            (void) sizeof(type_must_be_complete);
            delete p;
        }


        /// Checked array pointer delete function

        /// This function ensures that the pointer is a complete type.
        /// \tparam T The pointer type (must be a complete type).
        /// \param a The array pointer to be deleted.
        template<typename T>
        inline void checked_array_delete(T* a) {
            typedef char type_must_be_complete[ sizeof(T)? 1: -1 ];
            (void) sizeof(type_must_be_complete);
            delete [] a;
        }

        /// Function to free memory for a shared_ptr using free()

        /// Checks the pointer to make sure it is a complete type, you will get
        /// a compiler error if it is not.
        template <typename T>
        inline void checked_free(T* t) {
            typedef char type_must_be_complete[ sizeof(T)? 1: -1 ];
            (void) sizeof(type_must_be_complete);
            free(t);
        }

        /// Use this function with shared_ptr to do nothing for the pointer cleanup
        template <typename T>
        inline void no_delete(T*) { }

        /// Checked pointer delete functor

        /// This functor is used to delete a pointer. It ensures that the
        /// pointer is a complete type.
        /// \tparam T The pointer type (must be a complete type).
        template<typename T>
        struct CheckedDeleter {
            typedef void result_type;
            typedef T * argument_type;

            void operator()(T* p) const { checked_delete(p); }
        };

        /// Checked array pointer delete functor

        /// This functor is used to delete an array pointer. It ensures that the
        /// pointer is a complete type.
        /// \tparam T The pointer type (must be a complete type).
        template<typename T>
        struct CheckedArrayDeleter {
            typedef void result_type;
            typedef T * argument_type;

            void operator()(T* a) const { checked_array_delete(a); }
        };

        /// Deleter to free memory for a shared_ptr using free()

        /// Checks the pointer to make sure it is a complete type, you will get
        /// a compiler error if it is not.
        template<typename T>
        struct CheckedFree {
            typedef void result_type;
            typedef T * argument_type;

            void operator()(T* p) const { checked_fr   (p); }
        };

        /// Use this deleter with shared_ptr to do nothing for the pointer cleanup
        template<typename T>
        struct NoDeleter {
            typedef void result_type;
            typedef T * argument_type;

            void operator()(T*) const { no_delete(static_cast<T*>(NULL)); }
        };

    } // namespace detail
} // namespace madness

// make_shared / allocate_shared
#ifndef MADNESS_HAS_STD_MAKE_SHARED
#define MADNESS_HAS_STD_MAKE_SHARED 1

#if defined(MADNESS_USE_BOOST_TR1_MEMORY_HPP) || (defined(__Intel) && defined(MADNESS_HAS_BOOST_TR1))

// We are using Boost tr1 so we can use the Boost make_shared function
#include <boost/make_shared.hpp>
namespace std {
    using ::boost::make_shared;
    using ::boost::allocate_shared;
} // namespace std

#else

// We do not have make_shared/allocate_shared so we need to implement it here.
namespace std {

    /// Construct an shared pointer

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T());
    /// \endcode
    /// \tparam T The shared_ptr type
    /// \return A shared_ptr constructed with the given arguments
    template <class T>
    std::shared_ptr<T> make_shared() {
        return std::shared_ptr<T>(new T());
    }

    /// Construct an shared pointer with an allocator

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(), deleter(), a);
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A the allocator type
    /// \param a The allocator object used to allocate the shared_ptr
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A> std::shared_ptr<T> allocate_shared(A const & a) {
        return std::shared_ptr<T>(new T(), &madness::detail::checked_delete, a);
    }

    /// Construct an shared pointer

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(a1));
    /// \endcode
    /// \tparam T The shared_ptr type
    /// \tparam A1 pointer constructor argument 1 type
    /// \param a1 pointer constructor argument 1
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A1>
    std::shared_ptr<T> make_shared(A1 const & a1) {
        return std::shared_ptr<T>(new T(a1));
    }

    /// Construct an shared pointer with an allocator

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(a1), deleter(), a);
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A the allocator type
    /// \tparam A1 pointer constructor argument 1 type
    /// \param a The allocator object used to allocate the shared_ptr
    /// \param a1 pointer constructor argument 1
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A, class A1>
    std::shared_ptr<T> allocate_shared(A const & a, A1 const & a1) {
        return std::shared_ptr<T>(new T(a1), &madness::detail::checked_delete, a);
    }

    /// Construct an shared pointer

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(a1, a2));
    /// \endcode
    /// \tparam T The shared_ptr type
    /// \tparam A1 pointer constructor argument 1 type
    /// \tparam A2 pointer constructor argument 2 type
    /// \param a1 pointer constructor argument 1
    /// \param a2 pointer constructor argument 2
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A1, class A2>
    std::shared_ptr<T> make_shared(A1 const & a1, A2 const & a2) {
        return std::shared_ptr<T>(new T(a1, a2));
    }

    /// Construct an shared pointer with an allocator

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(a1, a2), deleter(), a);
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A the allocator type
    /// \tparam A1 pointer constructor argument 1 type
    /// \tparam A2 pointer constructor argument 2 type
    /// \param a The allocator object used to allocate the shared_ptr
    /// \param a1 pointer constructor argument 1
    /// \param a2 pointer constructor argument 2
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A, class A1, class A2>
    std::shared_ptr<T> allocate_shared(A const & a, A1 const & a1, A2 const & a2) {
        return std::shared_ptr<T>(new T(a1, a2), &madness::detail::checked_delete, a);
    }

    /// Construct an shared pointer

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(a1, a2, a3));
    /// \endcode
    /// \tparam T The shared_ptr type
    /// \tparam A1 pointer constructor argument 1 type
    /// \tparam A2 pointer constructor argument 2 type
    /// \tparam A3 pointer constructor argument 3 type
    /// \param a1 pointer constructor argument 1
    /// \param a2 pointer constructor argument 2
    /// \param a3 pointer constructor argument 3
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A1, class A2, class A3>
    std::shared_ptr<T> make_shared(A1 const & a1, A2 const & a2, A3 const & a3) {
        return std::shared_ptr<T>(new T(a1, a2, a3));
    }

    /// Construct an shared pointer with an allocator

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(a1, a2, a3), deleter(), a);
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A the allocator type
    /// \tparam A1 pointer constructor argument 1 type
    /// \tparam A2 pointer constructor argument 2 type
    /// \tparam A3 pointer constructor argument 3 type
    /// \param a The allocator object used to allocate the shared_ptr
    /// \param a1 pointer constructor argument 1
    /// \param a2 pointer constructor argument 2
    /// \param a3 pointer constructor argument 3
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A, class A1, class A2, class A3>
    std::shared_ptr<T> allocate_shared(A const & a, A1 const & a1, A2 const & a2, A3 const & a3) {
        return std::shared_ptr<T>(new T(a1, a2, a3), &madness::detail::checked_delete, a);
    }

    /// Construct an shared pointer

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(a1, a2, a3, a4));
    /// \endcode
    /// \tparam T The shared_ptr type
    /// \tparam A1 pointer constructor argument 1 type
    /// \tparam A2 pointer constructor argument 2 type
    /// \tparam A3 pointer constructor argument 3 type
    /// \tparam A4 pointer constructor argument 4 type
    /// \param a1 pointer constructor argument 1
    /// \param a2 pointer constructor argument 2
    /// \param a3 pointer constructor argument 3
    /// \param a4 pointer constructor argument 4
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A1, class A2, class A3, class A4>
    std::shared_ptr<T> make_shared(A1 const & a1, A2 const & a2, A3 const & a3, A4 const & a4) {
        return std::shared_ptr<T>(new T(a1, a2, a3, a4));
    }

    /// Construct an shared pointer with an allocator

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(a1, a2, a3, a4), deleter(), a);
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A the allocator type
    /// \tparam A1 pointer constructor argument 1 type
    /// \tparam A2 pointer constructor argument 2 type
    /// \tparam A3 pointer constructor argument 3 type
    /// \tparam A4 pointer constructor argument 4 type
    /// \param a The allocator object used to allocate the shared_ptr
    /// \param a1 pointer constructor argument 1
    /// \param a2 pointer constructor argument 2
    /// \param a3 pointer constructor argument 3
    /// \param a4 pointer constructor argument 4
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A, class A1, class A2, class A3, class A4>
    std::shared_ptr<T> allocate_shared(A const & a, A1 const & a1, A2 const & a2, A3 const & a3,
            A4 const & a4) {
        return std::shared_ptr<T>(new T(a1, a2, a3, a4), &madness::detail::checked_delete, a);
    }

    /// Construct an shared pointer

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(a1, a2, a3, a4, a5));
    /// \endcode
    /// \tparam T The shared_ptr type
    /// \tparam A1 pointer constructor argument 1 type
    /// \tparam A2 pointer constructor argument 2 type
    /// \tparam A3 pointer constructor argument 3 type
    /// \tparam A4 pointer constructor argument 4 type
    /// \tparam A5 pointer constructor argument 5 type
    /// \param a1 pointer constructor argument 1
    /// \param a2 pointer constructor argument 2
    /// \param a3 pointer constructor argument 3
    /// \param a4 pointer constructor argument 4
    /// \param a5 pointer constructor argument 5
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A1, class A2, class A3, class A4, class A5>
    std::shared_ptr<T> make_shared(A1 const & a1, A2 const & a2, A3 const & a3, A4 const & a4,
            A5 const & a5) {
        return std::shared_ptr<T>(new T(a1, a2, a3, a4, a5));
    }

    /// Construct an shared pointer with an allocator

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(a1, a2, a3, a4, a5, a6), deleter(), a);
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A the allocator type
    /// \tparam A1 pointer constructor argument 1 type
    /// \tparam A2 pointer constructor argument 2 type
    /// \tparam A3 pointer constructor argument 3 type
    /// \tparam A4 pointer constructor argument 4 type
    /// \tparam A5 pointer constructor argument 5 type
    /// \param a The allocator object used to allocate the shared_ptr
    /// \param a1 pointer constructor argument 1
    /// \param a2 pointer constructor argument 2
    /// \param a3 pointer constructor argument 3
    /// \param a4 pointer constructor argument 4
    /// \param a5 pointer constructor argument 5
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A, class A1, class A2, class A3, class A4, class A5>
    std::shared_ptr<T> allocate_shared(A const & a, A1 const & a1, A2 const & a2, A3 const & a3,
            A4 const & a4, A5 const & a5) {
        return std::shared_ptr<T>(new T(a1, a2, a3, a4, a5), &madness::detail::checked_delete, a);
    }

    /// Construct an shared pointer

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(a1, a2, a3, a4, a5, a6));
    /// \endcode
    /// \tparam T The shared_ptr type
    /// \tparam A1 pointer constructor argument 1 type
    /// \tparam A2 pointer constructor argument 2 type
    /// \tparam A3 pointer constructor argument 3 type
    /// \tparam A4 pointer constructor argument 4 type
    /// \tparam A5 pointer constructor argument 5 type
    /// \tparam A6 pointer constructor argument 6 type
    /// \param a1 pointer constructor argument 1
    /// \param a2 pointer constructor argument 2
    /// \param a3 pointer constructor argument 3
    /// \param a4 pointer constructor argument 4
    /// \param a5 pointer constructor argument 5
    /// \param a6 pointer constructor argument 6
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A1, class A2, class A3, class A4, class A5, class A6>
    std::shared_ptr<T> make_shared(A1 const & a1, A2 const & a2, A3 const & a3, A4 const & a4,
            A5 const & a5, A6 const & a6) {
        return std::shared_ptr<T>(new T(a1, a2, a3, a4, a5, a6));
    }

    /// Construct an shared pointer with an allocator

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(a1, a2, a3, a4, a5, a6), deleter(), a);
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A the allocator type
    /// \tparam A1 pointer constructor argument 1 type
    /// \tparam A2 pointer constructor argument 2 type
    /// \tparam A3 pointer constructor argument 3 type
    /// \tparam A4 pointer constructor argument 4 type
    /// \tparam A5 pointer constructor argument 5 type
    /// \tparam A6 pointer constructor argument 6 type
    /// \param a The allocator object used to allocate the shared_ptr
    /// \param a1 pointer constructor argument 1
    /// \param a2 pointer constructor argument 2
    /// \param a3 pointer constructor argument 3
    /// \param a4 pointer constructor argument 4
    /// \param a5 pointer constructor argument 5
    /// \param a6 pointer constructor argument 6
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A, class A1, class A2, class A3, class A4, class A5, class A6>
    std::shared_ptr<T> allocate_shared(A const & a, A1 const & a1, A2 const & a2, A3 const & a3,
            A4 const & a4, A5 const & a5, A6 const & a6) {
        return std::shared_ptr<T>(new T(a1, a2, a3, a4, a5, a6), &madness::detail::checked_delete,
                a);
    }

    /// Construct an shared pointer

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(a1, a2, a3, a4, a5, a6, a7));
    /// \endcode
    /// \tparam T The shared_ptr type
    /// \tparam A1 pointer constructor argument 1 type
    /// \tparam A2 pointer constructor argument 2 type
    /// \tparam A3 pointer constructor argument 3 type
    /// \tparam A4 pointer constructor argument 4 type
    /// \tparam A5 pointer constructor argument 5 type
    /// \tparam A6 pointer constructor argument 6 type
    /// \tparam A7 pointer constructor argument 7 type
    /// \param a1 pointer constructor argument 1
    /// \param a2 pointer constructor argument 2
    /// \param a3 pointer constructor argument 3
    /// \param a4 pointer constructor argument 4
    /// \param a5 pointer constructor argument 5
    /// \param a6 pointer constructor argument 6
    /// \param a7 pointer constructor argument 7
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A1, class A2, class A3, class A4, class A5, class A6, class A7>
    std::shared_ptr<T> make_shared(A1 const & a1, A2 const & a2, A3 const & a3, A4 const & a4,
            A5 const & a5, A6 const & a6, A7 const & a7) {
        return std::shared_ptr<T>(new T(a1, a2, a3, a4, a5, a6, a7));
    }

    /// Construct an shared pointer with an allocator

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(a1, a2, a3, a4, a5, a6, a7), deleter(), a);
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A the allocator type
    /// \tparam A1 pointer constructor argument 1 type
    /// \tparam A2 pointer constructor argument 2 type
    /// \tparam A3 pointer constructor argument 3 type
    /// \tparam A4 pointer constructor argument 4 type
    /// \tparam A5 pointer constructor argument 5 type
    /// \tparam A6 pointer constructor argument 6 type
    /// \tparam A7 pointer constructor argument 7 type
    /// \param a The allocator object used to allocate the shared_ptr
    /// \param a1 pointer constructor argument 1
    /// \param a2 pointer constructor argument 2
    /// \param a3 pointer constructor argument 3
    /// \param a4 pointer constructor argument 4
    /// \param a5 pointer constructor argument 5
    /// \param a6 pointer constructor argument 6
    /// \param a7 pointer constructor argument 7
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A, class A1, class A2, class A3, class A4, class A5, class A6,
            class A7>
    std::shared_ptr<T> allocate_shared(A const & a, A1 const & a1, A2 const & a2, A3 const & a3,
            A4 const & a4, A5 const & a5, A6 const & a6, A7 const & a7) {
        return std::shared_ptr<T>(new T(a1, a2, a3, a4, a5, a6, a7),
                &madness::detail::checked_delete, a);
    }

    /// Construct an shared pointer

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(a1, a2, a3, a4, a5, a6, a7, a8));
    /// \endcode
    /// \tparam T The shared_ptr type
    /// \tparam A1 pointer constructor argument 1 type
    /// \tparam A2 pointer constructor argument 2 type
    /// \tparam A3 pointer constructor argument 3 type
    /// \tparam A4 pointer constructor argument 4 type
    /// \tparam A5 pointer constructor argument 5 type
    /// \tparam A6 pointer constructor argument 6 type
    /// \tparam A7 pointer constructor argument 7 type
    /// \tparam A8 pointer constructor argument 8 type
    /// \param a1 pointer constructor argument 1
    /// \param a2 pointer constructor argument 2
    /// \param a3 pointer constructor argument 3
    /// \param a4 pointer constructor argument 4
    /// \param a5 pointer constructor argument 5
    /// \param a6 pointer constructor argument 6
    /// \param a7 pointer constructor argument 7
    /// \param a8 pointer constructor argument 8
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A1, class A2, class A3, class A4, class A5, class A6, class A7,
            class A8>
    std::shared_ptr<T> make_shared(A1 const & a1, A2 const & a2, A3 const & a3, A4 const & a4,
            A5 const & a5, A6 const & a6, A7 const & a7, A8 const & a8) {
        return std::shared_ptr<T>(new T(a1, a2, a3, a4, a5, a6, a7, a8));
    }

    /// Construct an shared pointer with an allocator

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(a1, a2, a3, a4, a5, a6, a7, a8), deleter(), a);
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A the allocator type
    /// \tparam A1 pointer constructor argument 1 type
    /// \tparam A2 pointer constructor argument 2 type
    /// \tparam A3 pointer constructor argument 3 type
    /// \tparam A4 pointer constructor argument 4 type
    /// \tparam A5 pointer constructor argument 5 type
    /// \tparam A6 pointer constructor argument 6 type
    /// \tparam A7 pointer constructor argument 7 type
    /// \tparam A8 pointer constructor argument 8 type
    /// \param a The allocator object used to allocate the shared_ptr
    /// \param a1 pointer constructor argument 1
    /// \param a2 pointer constructor argument 2
    /// \param a3 pointer constructor argument 3
    /// \param a4 pointer constructor argument 4
    /// \param a5 pointer constructor argument 5
    /// \param a6 pointer constructor argument 6
    /// \param a7 pointer constructor argument 7
    /// \param a8 pointer constructor argument 8
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A, class A1, class A2, class A3, class A4, class A5, class A6,
            class A7, class A8>
    std::shared_ptr<T> allocate_shared(A const & a, A1 const & a1, A2 const & a2, A3 const & a3,
            A4 const & a4, A5 const & a5, A6 const & a6, A7 const & a7, A8 const & a8) {
        return std::shared_ptr<T>(new T(a1, a2, a3, a4, a5, a6, a7, a8),
                &madness::detail::checked_delete, a);
    }

    /// Construct an shared pointer

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(a1, a2, a3, a4, a5, a6, a7, a8, a9));
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A1 pointer constructor argument 1 type
    /// \tparam A2 pointer constructor argument 2 type
    /// \tparam A3 pointer constructor argument 3 type
    /// \tparam A4 pointer constructor argument 4 type
    /// \tparam A5 pointer constructor argument 5 type
    /// \tparam A6 pointer constructor argument 6 type
    /// \tparam A7 pointer constructor argument 7 type
    /// \tparam A8 pointer constructor argument 8 type
    /// \tparam A9 pointer constructor argument 9 type
    /// \param a1 pointer constructor argument 1
    /// \param a2 pointer constructor argument 2
    /// \param a3 pointer constructor argument 3
    /// \param a4 pointer constructor argument 4
    /// \param a5 pointer constructor argument 5
    /// \param a6 pointer constructor argument 6
    /// \param a7 pointer constructor argument 7
    /// \param a8 pointer constructor argument 8
    /// \param a9 pointer constructor argument 9
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A1, class A2, class A3, class A4, class A5, class A6, class A7,
            class A8, class A9>
    std::shared_ptr<T> make_shared(A1 const & a1, A2 const & a2, A3 const & a3, A4 const & a4,
            A5 const & a5, A6 const & a6, A7 const & a7, A8 const & a8, A9 const & a9) {
        return std::shared_ptr<T>(new T(a1, a2, a3, a4, a5, a6, a7, a8, a9));
    }

    /// Construct an shared pointer with an allocator

    /// This function has the same effect as:
    /// \code
    /// return std::shared_ptr<T>(new T(a1, a2, a3, a4, a5, a6, a7, a8, a9), deleter(), a);
    /// \endcode
    /// where deleter is a deleter object supplied by this function.
    /// \tparam T The shared_ptr type
    /// \tparam A the allocator type
    /// \tparam A1 pointer constructor argument 1 type
    /// \tparam A2 pointer constructor argument 2 type
    /// \tparam A3 pointer constructor argument 3 type
    /// \tparam A4 pointer constructor argument 4 type
    /// \tparam A5 pointer constructor argument 5 type
    /// \tparam A6 pointer constructor argument 6 type
    /// \tparam A7 pointer constructor argument 7 type
    /// \tparam A8 pointer constructor argument 8 type
    /// \tparam A9 pointer constructor argument 9 type
    /// \param a The allocator object used to allocate the shared_ptr
    /// \param a1 pointer constructor argument 1
    /// \param a2 pointer constructor argument 2
    /// \param a3 pointer constructor argument 3
    /// \param a4 pointer constructor argument 4
    /// \param a5 pointer constructor argument 5
    /// \param a6 pointer constructor argument 6
    /// \param a7 pointer constructor argument 7
    /// \param a8 pointer constructor argument 8
    /// \param a9 pointer constructor argument 9
    /// \return A shared_ptr constructed with the given arguments
    template <class T, class A, class A1, class A2, class A3, class A4, class A5, class A6,
            class A7, class A8, class A9>
    std::shared_ptr<T> allocate_shared(A const & a, A1 const & a1, A2 const & a2, A3 const & a3,
            A4 const & a4, A5 const & a5, A6 const & a6, A7 const & a7, A8 const & a8,
            A9 const & a9) {
        return std::shared_ptr<T>(new T(a1, a2, a3, a4, a5, a6, a7, a8, a9),
                &madness::detail::checked_delete, a);
    }

} // namespace std

#endif // MADNESS_USE_BOOST_TR1_MEMORY_HPP
#endif // MADNESS_HAS_STD_MAKE_SHARED

#endif // MADNESS_WORLD_SHAREDPTR_H__INCLUDED
