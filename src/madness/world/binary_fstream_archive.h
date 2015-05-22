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

#ifndef MADNESS_WORLD_BINARY_FSTREAM_ARCHIVE_H__INCLUDED
#define MADNESS_WORLD_BINARY_FSTREAM_ARCHIVE_H__INCLUDED

/**
 \file binary_fstream_archive.h
 \brief Implements an archive wrapping a binary filestream.
 \ingroup serialization
*/

#include <type_traits>
#include <fstream>
#include <memory>
#include <madness/world/archive.h>

namespace madness {
    namespace archive {

        /// \addtogroup serialization
        /// @{

        /// Wraps an archive around a binary filestream for output.
        class BinaryFstreamOutputArchive : public BaseOutputArchive {
            static const std::size_t IOBUFSIZE = 4*1024*1024; ///< Buffer size.
            std::shared_ptr<char> iobuf; ///< Buffer.
            mutable std::ofstream os; ///< The filestream.

        public:
            /// Default constructor.

            /// The filename and open modes are optional here; they can be
            /// specified later by calling \c open().
            /// \param[in] filename Name of the file to write to.
            /// \param[in] mode I/O attributes for opening the file.
            BinaryFstreamOutputArchive(const char* filename = nullptr,
                                       std::ios_base::openmode mode = std::ios_base::binary | \
                                                                      std::ios_base::out | std::ios_base::trunc);

            /// Write to the filestream.

            /// The function only appears (due to \c enable_if) if \c T is
            /// serializable.
            /// \tparam T The type of data to be written.
            /// \param[in] t Location of the data to be written.
            /// \param[in] n The number of data items to be written.
            template <class T>
            inline
            typename std::enable_if< madness::is_serializable<T>::value, void >::type
            store(const T* t, long n) const {
                os.write((const char *) t, n*sizeof(T));
            }

            /// Open the filestream.

            /// \param[in] filename The name of the file.
            /// \param[in] mode I/O attributes for opening the file.
            void open(const char* filename,
                      std::ios_base::openmode mode = std::ios_base::binary | \
                                                     std::ios_base::out | std::ios_base::trunc);

            /// Close the filestream.
            void close();

            /// Flush the filestream.
            void flush();
        };

        /// Wraps an archive around a binary filestream for input.
        class BinaryFstreamInputArchive : public BaseInputArchive {
            static const std::size_t IOBUFSIZE = 4*1024*1024; ///< Buffer size.
            std::shared_ptr<char> iobuf; ///< Buffer.
            mutable std::ifstream is; ///< The filestream.

        public:
            /// Default constructor.

            /// The filename and open modes are optional here; they can be
            /// specified later by calling \c open().
            /// \param[in] filename Name of the file to read from.
            /// \param[in] mode I/O attributes for opening the file.
            BinaryFstreamInputArchive(const char* filename = nullptr, std::ios_base::openmode mode = std::ios_base::binary | std::ios_base::in);

            /// Load from the filestream.

            /// The function only appears (due to \c enable_if) if \c T is
            /// serializable.
            /// \tparam T The type of data to be read.
            /// \param[out] t Where to put the loaded data.
            /// \param[in] n The number of data items to be loaded.
            template <class T>
            inline
            typename std::enable_if< madness::is_serializable<T>::value, void >::type
            load(T* t, long n) const {
                is.read((char *) t, n*sizeof(T));
            }

            /// Open the filestream.

            /// \param[in] filename Name of the file to read from.
            /// \param[in] mode I/O attributes for opening the file.
            void open(const char* filename, std::ios_base::openmode mode = std::ios_base::binary | std::ios_base::in);

            /// Close the filestream.
            void close();
        };

        /// @}
    }
}
#endif // MADNESS_WORLD_BINARY_FSTREAM_ARCHIVE_H__INCLUDED
