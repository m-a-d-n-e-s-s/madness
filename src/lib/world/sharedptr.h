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
/// \brief Minimal, thread safe, modified (and renamed) Boost-like SharedPtr & SharedArray

#include <madness_config.h>
#include <memory> // for shared_ptr

#if defined(MADNESS_HAS_STD_TR1_SHARED_PTR) && !defined(MADNESS_HAS_STD_SHARED_PTR)
#define MADNESS_HAS_STD_SHARED_PTR
// shard_ptr is in std::tr1 but we want it in std namespace
namespace std {
    using ::std::tr1::bad_weak_ptr;
    using ::std::tr1::shared_ptr;
#if !BOOST_WORKAROUND(__BORLANDC__, < 0x0582)
    using ::std::tr1::swap;
#endif
    using ::std::tr1::static_pointer_cast;
    using ::std::tr1::dynamic_pointer_cast;
    using ::std::tr1::const_pointer_cast;
    using ::std::tr1::get_deleter;
    using ::std::tr1::weak_ptr;
    using ::std::tr1::enable_shared_from_this;
}

#endif

namespace madness {

    template <typename>
    class RemoteReference;

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

    }
}
#endif // MADNESS_WORLD_SHAREDPTR_H__INCLUDED
