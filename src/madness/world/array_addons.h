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

#ifndef MADNESS_WORLD_ARRAY_ADDONS_H__INCLUDED
#define MADNESS_WORLD_ARRAY_ADDONS_H__INCLUDED

/**
 \file array_addons.h
 \brief Supplements to the \c std::array class, such as I/O operations,
    for convenience.
 \ingroup containers
*/

#include <madness/madness_config.h>
#include <madness/world/madness_exception.h>
#include <madness/world/worldhash.h>
#include <array>
#include <iostream>

namespace madness {

    /// Output \c std::array to stream for human consumption.

    /// \tparam T The type of data stored in the array.
    /// \tparam N The size of the array.
    /// \param[in,out] s The output stream.
    /// \param[in] a The array to be output.
    /// \return The output stream.
    template <typename T, std::size_t N>
    std::ostream& operator<<(std::ostream& s, const std::array<T,N>& a) {
        s << "[";
        for(std::size_t i=0; i<N; ++i) {
            s << a[i];
            if (i != (N-1)) s << ",";
        }
        s << "]";
        return s;
    }

    /// Hash std::array with madness hash.

    /// \tparam T The type of data stored in the array.
    /// \tparam N The size of the array.
    /// \param[in] a The array.
    /// \return The hash.
    template <typename T, std::size_t N>
    madness::hashT hash_value(const std::array<T,N>& a) {
        // Use this version of range for potential optimization.
        return madness::hash_range(a.data(), N);
    }

    // Serialize std::array objects
    namespace archive {

        template <class Archive, class T>
        struct ArchiveStoreImpl;
        template <class Archive, class T>
        struct ArchiveLoadImpl;

        /// Serialize a \c std::array.

        /// \tparam Archive the archive type.
        /// \tparam T The data type stored in the \c std::array.
        /// \tparam N The size of the \c std::array.
        template <class Archive, typename T, std::size_t N>
        struct ArchiveStoreImpl< Archive, std::array<T, N> > {

            /// Store a \c std::array of plain data.

            /// \param[in] ar The archive.
            /// \param[in] v The array object to be serialized.
            template <typename U = T, typename = std::enable_if_t<is_serializable<Archive,U>::value>>
            static inline void store(const Archive& ar, const std::array<U, N>& v) {
                //MAD_ARCHIVE_DEBUG(std::cout << "serialize std::array<T," << N << ">, with T plain data" << std::endl);
                ar & v.size();
                ar & wrap(v.data(),v.size());
            }

            /// Store a \c std::array of non-plain data.

            /// \param[in] ar The archive.
            /// \param[in] v The array object to be serialized.
            template <typename U = T>
            static inline void store(const Archive& ar, const std::array<U, N>& v, std::enable_if_t<!is_serializable<Archive,U>::value>* = nullptr) {
                //MAD_ARCHIVE_DEBUG(std::cout << "serialize std::array<T," << N << ">, with T non-plain data" << std::endl);
                ar & v.size();
                for(auto&& elem: v) {
                    ar & elem;
                }
            }

        };

        /// Deserialize a \c std::array. \c MADNESS_ASSERT 's that the size matches.

        /// \tparam Archive the archive type.
        /// \tparam T The data type stored in the \c std::array.
        /// \tparam N The size of the \c std::array.
        template <class Archive, typename T, std::size_t N>
        struct ArchiveLoadImpl< Archive, std::array<T, N> > {

            /// Load a \c std::array of plain data.

            /// \param[in] ar The archive.
            /// \param[out] v The array to be deserialized.
            template <typename U = T, typename = std::enable_if_t<is_serializable<Archive,U>::value>>
            static void load(const Archive& ar, std::array<U, N>& v) {
                //MAD_ARCHIVE_DEBUG(std::cout << "deserialize std::array<T," << N << ">, with T plain data" << std::endl);
                std::size_t n = 0ul;
                ar & n;
                MADNESS_ASSERT(n == v.size());
                ar & wrap((T *) v.data(),n);
            }

            /// Load a \c std::array of non-plain data.

            /// \param[in] ar The archive.
            /// \param[out] v The array to be deserialized.
            template <typename U = T>
            static void load(const Archive& ar, std::array<U, N>& v, std::enable_if_t<!is_serializable<Archive,U>::value>* = nullptr) {
                //MAD_ARCHIVE_DEBUG(std::cout << "deserialize std::array<T," << N << ">, with T non-plain data" << std::endl);
                std::size_t n = 0ul;
                ar & n;
                MADNESS_ASSERT(n == v.size());
                for(auto&& elem: v) {
                    ar & elem;
                }
            }

        };


    } // namespace archive

} // namespace madness

#endif // MADNESS_WORLD_ARRAY_ADDONS_H__INCLUDED
