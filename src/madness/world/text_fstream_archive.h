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

#ifndef MADNESS_WORLD_TEXT_FSTREAM_ARCHIVE_H__INCLUDED
#define MADNESS_WORLD_TEXT_FSTREAM_ARCHIVE_H__INCLUDED

/**
 \file text_fstream_archive.h
 \brief Implements an archive wrapping text filestream.
 \ingroup serialization
*/

#include <type_traits>
#include <fstream>
#include <cstring>
#include <madness/world/archive.h>
#include <madness/world/print.h>  // this injects operator<<(std::ostream,T) for common Ts

namespace madness {
    namespace archive {

        /// \addtogroup serialization
        /// @{

        /// Wraps an archive around a text filestream for output.
        class TextFstreamOutputArchive : public BaseOutputArchive {
            mutable std::ofstream os; ///< The filestream.

        public:
            /// Default constructor.

            /// The filename and open modes are optional here; they can be
            /// specified later by calling \c open().
            /// \param[in] filename Name of the file to write to.
            /// \param[in] mode I/O attributes for opening the file.
            TextFstreamOutputArchive(const char* filename = nullptr,
                    std::ios_base::openmode mode=std::ios_base::binary | std::ios_base::out | std::ios_base::trunc)
            {
                if (filename)
                    open(filename, mode);
            }

            /// Stores the cookie for runtime type-checking into the archive as a tag.

            /// \tparam T The type of data to be stored between the tags.
            template <class T>
            void store_start_tag() const {
                const std::size_t bufsize=256;
                char tag[bufsize];
                unsigned char cookie = archive_typeinfo<T>::cookie;
                snprintf(tag,bufsize,"<t%d>", cookie);
                os << tag << std::endl;
                MAD_ARCHIVE_DEBUG(std::cout << "textarchive: tag = " << tag << std::endl);
            }

            /// Closes the "cookie" tag for runtime type-checking.

            /// \tparam T The type of data to be stored between the tags.
            template <class T>
            void store_end_tag() const {
                const std::size_t bufsize=256;
                char tag[bufsize];
                unsigned char cookie = archive_typeinfo<T>::cookie;
                snprintf(tag,bufsize,"</t%d>",cookie);
                os << tag << std::endl;
            }

            /// Store data to the filestream.

            /// The function only appears (due to \c enable_if) if \c T is
            /// serializable.
            /// \tparam T The type of data to be written.
            /// \param[in] t Location of the data to be written.
            /// \param[in] n The number of data items to be written.
            template <class T>
            typename std::enable_if< madness::is_ostreammable_v<T> >::type
            store(const T* t, long n) const {
                using madness::operators::operator<<;
                for (long i=0; i<n; ++i)
                    os << t[i] << std::endl;
            }

            /// Store a character string, escaping '&', '<' and '> along the way.

            /// \param[in] t The character string.
            void store(const char* t, long /*n*/) const;

            /// Store a character string, without escaping characters.

            /// \param[in] t The character string.
            /// \param[in] n The number of characters to store.
            void store(const unsigned char* t, long n) const {
                for (long i=0; i<n; ++i)
                    os << (unsigned int) t[i] << std::endl;
            }

            /// Open the filestream.

            /// \param[in] filename The name of the file.
            /// \param[in] mode I/O attributes for opening the file.
            void open(const char* filename, std::ios_base::openmode mode = std::ios_base::out | std::ios_base::trunc);

            /// Close the filestream.
            void close();

            /// Flush the filestream.
            void flush() {
                if (os.is_open())
                    os.flush();
            }

            /// Destructor.
            ~TextFstreamOutputArchive() {
                close();
            }
        }; // class TextFstreamOutputArchive


        /// Wraps an archive around a text filestream for input.
        class TextFstreamInputArchive : public BaseInputArchive {
        private:
            mutable std::ifstream is; ///< The filestream.

            /// Eat the EOL after each entry to enable a `char`-by-`char` read of strings.
            void eat_eol() const;

        public:
            /// Default constructor.

            /// The filename and open modes are optional here; they can be
            /// specified later by calling \c open().
            /// \param[in] filename Name of the file to read from.
            /// \param[in] mode I/O attributes for opening the file.
            TextFstreamInputArchive(const char* filename = nullptr,
                    std::ios_base::openmode mode = std::ios_base::in)
            {
                if (filename)
                    open(filename, mode);
            }

            /// Check the "cookie" tag in the archive for runtime type-checking.

            /// \tparam T The expected data type.
            /// \throw MadnessException if the tag does not match that of the
            ///     expected type.
            template <class T>
            void check_start_tag(bool end=false) const {
                const std::size_t bufsize=256;
                char tag[bufsize], ftag[bufsize];
                is.getline(ftag,bufsize);
                unsigned char cookie = archive_typeinfo<T>::cookie;
                if (end)
                    snprintf(tag,bufsize,"</t%d>",cookie);
                else
                    snprintf(tag,bufsize,"<t%d>",cookie);

                if (strcmp(tag,ftag) != 0) {
                    std::cout << "TextFstreamInputArchive: type mismatch: expected=" << tag
                              << " "
                              << archive_type_names[cookie]
                              << " "
                              << " got=" << ftag << std::endl;
                    MADNESS_EXCEPTION("TextFstreamInputArchive: check_tag: types do not match/corrupt file", 1);
                }
            }

            /// Read the closing "cookie" tag.

            /// \tparam T The expected data type between the tags.
            template <class T>
            inline void check_end_tag() const {
                check_start_tag<T>(true);
            }

            /// Load from the filestream.

            /// The function only appears (due to \c enable_if) if \c T is
            /// serializable.
            /// \tparam T The type of data to be read.
            /// \param[out] t Where to put the loaded data.
            /// \param[in] n The number of data items to be loaded.
            template <class T>
            typename std::enable_if< madness::is_istreammable_v<T> >::type
            load(T* t, long n) const {
                using madness::operators::operator>>;
                for (long i=0; i<n; ++i) is >> t[i];
                eat_eol();
            }

            /// Load characters from the filestream, interpreting escaped characters along the way.

            /// \param[out] t Where to put the loaded characters.
            /// \param[in] n The number of characters to be loaded.
            void load(unsigned char* t, long n) const;

            /// Load characters from the filestream, without converting escaped characters.

            /// \param[out] t Where to put the loaded characters.
            /// \param[in] n The number of characters to be loaded.
            void load(char* t, long n) const;

            /// Open the filestream.

            /// \param[in] filename The name of the file.
            /// \param[in] mode I/O attributes for opening the file.
            void open(const char* filename,
                      std::ios_base::openmode mode = std::ios_base::in);

            /// Close the filestream.
            void close() {
                is.close();
            }
        }; // class TextFstreamInputArchive

        /// Implement pre/postamble storage routines for a \c TextFstreamOutputArchive.

        /// \tparam T The type to be stored.
        template <class T>
        struct ArchivePrePostImpl<TextFstreamOutputArchive, T> {
            /// Write the preamble to the archive.

            /// \param[in] ar The archive.
            static void preamble_store(const TextFstreamOutputArchive& ar) {
                ar.store_start_tag<T>();
            }

            /// Write the postamble to the archive.

            /// \param[in] ar The archive.
            static inline void postamble_store(const TextFstreamOutputArchive& ar) {
                ar.store_end_tag<T>();
            }
        }; // struct ArchivePrePostImpl<TextFstreamOutputArchive,T>

        /// Implement pre/postamble load routines for a \c TextFstreamInputArchive.

        /// \tparam T The expected type to be loaded.
        template <class T>
        struct ArchivePrePostImpl<TextFstreamInputArchive,T> {
            /// Load the preamble and perform type-checking.

            /// \param[in] ar The archive.
            static inline void preamble_load(const TextFstreamInputArchive& ar) {
                ar.check_start_tag<T>();
            }

            /// Load the postamble and perform type-checking.

            /// \param[in] ar The archive.
            static inline void postamble_load(const TextFstreamInputArchive& ar) {
                ar.check_end_tag<T>();
            }
        }; // struct ArchivePrePostImpl<TextFstreamInputArchive,T>

        /// @}
    }
}

#endif // MADNESS_WORLD_TEXT_FSTREAM_ARCHIVE_H__INCLUDED
