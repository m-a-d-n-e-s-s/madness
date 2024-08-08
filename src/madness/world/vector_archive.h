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

#ifndef MADNESS_WORLD_VECTOR_ARCHIVE_H__INCLUDED
#define MADNESS_WORLD_VECTOR_ARCHIVE_H__INCLUDED

/**
 \file vector_archive.h
 \brief Implements an archive wrapping an STL \c vector.
 \ingroup serialization

 \todo With a bit of thought this could be generalized to several STL containers.
*/

#include <type_traits>
#include <vector>
#include <cstring>
#include <madness/world/archive.h>

namespace madness {
    namespace archive {

        /// \addtogroup serialization
        /// @{

        /// Wraps an archive around an STL \c vector for output.
        class VectorOutputArchive : public BaseOutputArchive {
        public:
            mutable std::vector<unsigned char>* v; ///< The STL vector being wrapped.

        public:
            /// Create a buffer to wrap the specified \c vector.

            /// \param[in] v The \c vector.
            /// \param[in] hint The minimum capacity of the vector.
            VectorOutputArchive(std::vector<unsigned char>& v, std::size_t hint=262144)
                    : v(&v) {
                open(hint);
            };

            /// Appends data to the end of the vector.

            /// \todo Verify/complete the documentation.
            /// \tparam T The type of data to be appended.
            /// \param[in] t Pointer to the data to be appended.
            /// \param[in] n The number of data items to be appended.
            /// \return Description needed.
            template <class T>
            inline
            typename std::enable_if< madness::is_trivially_serializable<T>::value, void >::type
            store(const T* t, long n) const {
                const unsigned char* ptr = (unsigned char*) t;
                v->insert(v->end(),ptr,ptr+n*sizeof(T));
            }

            /// Clear any data in the vector and ensure its capacity is at least \c hint.

            /// \param[in] hint The minimum capacity for the vector.
            void open(std::size_t hint=262144) {
                v->clear();
                v->reserve(hint);
            };

            /// Close the archive.
            void close() {};

            /// Flush the archive.
            void flush() {};
        };


        /// Wraps an archive around an STL \c vector for input.
        class VectorInputArchive : public BaseInputArchive {
            mutable std::vector<unsigned char>* v; ///< The STL vector being wrapped.
            mutable std::size_t i; ///< Current input location.

        public:
            /// Create a buffer to wrap the specified \c vector.

            /// \param[in] v The \c vector.
            VectorInputArchive(std::vector<unsigned char>& v) : v(&v) , i(0) {}

            /// Load data from the vector.

            /// The function only appears (due to \c enable_if) if \c T is
            /// serializable.
            /// \tparam T The type of data to be loaded.
            /// \param[out] t Where to store the loaded data.
            /// \param[in] n The number of data items to be loaded.
            template <class T>
            inline
            typename std::enable_if< madness::is_trivially_serializable<T>::value, void >::type
            load(T* t, long n) const {
                std::size_t m = n*sizeof(T);
                if (m+i >  v->size()) MADNESS_EXCEPTION("VectorInputArchive: reading past end", m+1);
                memcpy((unsigned char*) t, &((*v)[i]), m);
                i += m;
            }

            /// Open the archive.
            void open() {};

            /// Reset the read location to the beginning of the \c vector.
            void rewind() const {
                i=0;
            };

            /// Get the amount of space left to be read from the \c vector.

            /// \return The amount of space left to be read from the \c vector.
            std::size_t nbyte_avail() const {
                return v->size()-i;
            };

            /// Close the archive.
            void close() {}
        };

        /// Implementation of functions for storing the pre/postamble in Vector archives.

        /// \attention No type checking over Vector buffers, for efficiency.
        /// \tparam T The data type.
        template <class T>
        struct ArchivePrePostImpl<VectorOutputArchive,T> {
            /// Store the preamble.

            /// \param[in] ar The archive.
            static void preamble_store(const VectorOutputArchive& ar) {};

            /// Store the postamble.

            /// \param[in] ar The archive.
            static inline void postamble_store(const VectorOutputArchive& ar) {};
        };

        /// Implementation of functions for loading the pre/postamble in Vector archives.

        /// \attention No type checking over Vector buffers, for efficiency.
        /// \tparam T The data type.
        template <class T>
        struct ArchivePrePostImpl<VectorInputArchive,T> {
            /// Load the preamble.

            /// \param[in] ar The archive.
            static inline void preamble_load(const VectorInputArchive& ar) {};

            /// Load the postamble.

            /// \param[in] ar The archive.
            static inline void postamble_load(const VectorInputArchive& ar) {};
        };

        /// @}
    }
}

#endif // MADNESS_WORLD_VECTOR_ARCHIVE_H__INCLUDED
