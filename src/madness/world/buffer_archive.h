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

#ifndef MADNESS_WORLD_BUFFER_ARCHIVE_H__INCLUDED
#define MADNESS_WORLD_BUFFER_ARCHIVE_H__INCLUDED

/**
 \file buffer_archive.h
 \brief Implements an archive wrapping a memory buffer.
 \ingroup serialization
*/

#include <type_traits>
#include <madness/world/archive.h>
#include <madness/world/print.h>
#include <cstring>

namespace madness {
    namespace archive {

        /// \addtogroup serialization
        /// @{

        /// Wraps an archive around a memory buffer for output.

        /// \note Type checking is disabled for efficiency.
        ///
        /// \throw madness::MadnessException in case of buffer overflow.
        ///
        /// The default constructor can also be used to count stuff.
        class BufferOutputArchive : public BaseOutputArchive {
        private:
            unsigned char * const ptr; ///< The memory buffer.
            const std::size_t nbyte; ///< Buffer size.
            mutable std::size_t i; /// Current output location.
            bool countonly; ///< If true just count, don't copy.

        public:
            /// Default constructor; the buffer will only count data.
            BufferOutputArchive()
                    : ptr(nullptr), nbyte(0), i(0), countonly(true) {}

            /// Constructor that assigns a buffer.

            /// \param[in] ptr Pointer to the buffer.
            /// \param[in] nbyte Size of the buffer.
            BufferOutputArchive(void* ptr, std::size_t nbyte)
                    : ptr((unsigned char *) ptr), nbyte(nbyte), i(0), countonly(false) {}

            /// Stores (counts) data into the memory buffer.

            /// The function only appears (due to \c enable_if) if \c T is
            /// serializable.
            /// \tparam T Type of the data to be stored (counted).
            /// \param[in] t Pointer to the data to be stored (counted).
            /// \param[in] n Size of data to be stored (counted).
            template <typename T>
            inline
            typename std::enable_if< madness::is_trivially_serializable<T>::value, void >::type
            store(const T* t, long n) const {
                std::size_t m = n*sizeof(T);
                if (countonly) {
                    i += m;
                }
                else if (i+m > nbyte) {
                    madness::print("BufferOutputArchive:ptr,nbyte,i,n,m,i+m:",(void *)ptr,nbyte,i,n,m,i+m);
                    MADNESS_ASSERT(i+m<=nbyte);
                }
                else {
MADNESS_PRAGMA_GCC(diagnostic push)
MADNESS_PRAGMA_GCC(diagnostic ignored "-Wmaybe-uninitialized")
		  memcpy(ptr+i, t, m);
MADNESS_PRAGMA_GCC(diagnostic pop)
                    i += m;
                }
            }

            /// Open a buffer with a specific size.
            void open(std::size_t /*hint*/) {}

            /// Close the archive.
            void close() {}

            /// Flush the archive.
            void flush() {}

            /// Determine if this buffer is used for counting.

            /// \return True if this buffer is only used for counting.
            bool count_only() const { return countonly; }

            /// Return the amount of data stored (counted) in the buffer.

            /// \return The amount of data stored (counted) in the buffer.
            inline std::size_t size() const {
                return i;
            };
        };


        /// Wraps an archive around a memory buffer for input.

        /// \note Type checking is disabled for efficiency.
        ///
        /// \throw madness::MadnessException in case of buffer overrun.
        class BufferInputArchive : public BaseInputArchive {
        private:
            const unsigned char* const ptr; ///< The memory buffer.
            const std::size_t nbyte; ///< Buffer size.
            mutable std::size_t i; ///< Current input location.

        public:
            /// Constructor that assigns a buffer.

            /// \param[in] ptr Pointer to the buffer.
            /// \param[in] nbyte Size of the buffer.
            BufferInputArchive(const void* ptr, std::size_t nbyte)
                    : ptr((const unsigned char *) ptr), nbyte(nbyte), i(0) {};

            /// Reads data from the memory buffer.

            /// The function only appears (due to \c enable_if) if \c T is
            /// serializable.
            /// \tparam T Type of the data to be read.
            /// \param[out] t Where to store the read data.
            /// \param[in] n Size of data to be read.
            template <class T>
            inline
            typename std::enable_if< madness::is_trivially_serializable<T>::value, void >::type
            load(T* t, long n) const {
                std::size_t m = n*sizeof(T);
                MADNESS_ASSERT(m+i <=  nbyte);
                memcpy((unsigned char*) t, ptr+i, m);
                i += m;
            }

            /// Open the archive.
            void open() {};

            /// Reset the read location to the beginning of the buffer.
            void rewind() const {
                i=0;
            };

            /// Get the amount of space yet to be read from the buffer.

            /// \return The amount of space yet to be read from the buffer.
            std::size_t nbyte_avail() const {
                return nbyte-i;
            };

            /// Close the archive.
            void close() {}
        };

        /// Implement pre/postamble storage routines for a \c BufferOutputArchive.

        /// \note No type checking over the buffer stream, for efficiency.
        /// \tparam T The type to be stored.
        template <class T>
        struct ArchivePrePostImpl<BufferOutputArchive, T> {
            /// Write the preamble to the archive.
            static inline void preamble_store(const BufferOutputArchive& /*ar*/) {}

            /// Write the postamble to the archive.
            static inline void postamble_store(const BufferOutputArchive& /*ar*/) {}
        };

        /// Implement pre/postamble load routines for a \c BufferInputArchive.

        /// \note No type checking over \c Buffer stream, for efficiency.
        /// \tparam T The type to be loaded.
        template <class T>
        struct ArchivePrePostImpl<BufferInputArchive, T> {
            /// Load the preamble.
            static inline void preamble_load(const BufferInputArchive& /*ar*/) {}

            /// Load the postamble.
            static inline void postamble_load(const BufferInputArchive& /*ar*/) {}
        };

        /// @}
    }
}
#endif // MADNESS_WORLD_BUFFER_ARCHIVE_H__INCLUDED
