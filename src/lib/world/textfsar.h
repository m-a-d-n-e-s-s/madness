/*
  This file is part of MADNESS.

  Copyright (C) <2007> <Oak Ridge National Laboratory>

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


#ifndef MAD_TEXTFSAR_H
#define MAD_TEXTFSAR_H

/// \file textfsar.h
/// \brief Implements archive wrapping text file stream

#include <fstream>
#include <world/archive.h>

namespace madness {
    namespace archive {

        /// Wraps an archive around a text file stream for output
        class TextFstreamOutputArchive : public BaseOutputArchive {
            mutable std::ofstream os;

        public:
            TextFstreamOutputArchive(const char* filename = 0,
                                     std::ios_base::openmode mode=std::ios_base::binary |
                                                                  std::ios_base::out | std::ios_base::trunc) {
                if (filename) open(filename, mode);
            }

            template <class T>
            void store_start_tag() const {
                char tag[256];
                unsigned char cookie = archive_typeinfo<T>::cookie;
                sprintf(tag,"<t%d>", cookie);
                os << tag << std::endl;
                MAD_ARCHIVE_DEBUG(std::cout << "textarchive: tag = " << tag << std::endl);
            }

            template <class T>
            void store_end_tag() const {
                char tag[256];
                unsigned char cookie = archive_typeinfo<T>::cookie;
                sprintf(tag,"</t%d>",cookie);
                os << tag << std::endl;
            }

            template <class T>
            typename madness::enable_if< madness::is_serializable<T>, void >::type
            store(const T* t, long n) const {
                for (long i=0; i<n; i++) os << t[i] << std::endl;
            }

            void store(const char* t, long n) const {
                // store character string, escaping &, < and > along the way
                while (*t) {
                    char c = *t++;
                    if (c == '\\') {
                        os.put('\\');
                        os.put('\\');
                    }
                    else if (c == '<') {
                        os.put('\\');
                        os.put('l');
                    }
                    else if (c == '>') {
                        os.put('\\');
                        os.put('r');
                    }
                    else {
                        os.put(c);
                    }
                }
                os << std::endl;
            };

            void store(const unsigned char* t, long n) const {
                for (long i=0; i<n; i++) os << (unsigned int) t[i] << std::endl;
            }

            void open(const char* filename,
                      std::ios_base::openmode mode = std::ios_base::out |
                                                     std::ios_base::trunc) {
                os.open(filename, mode);
                os.setf(std::ios::scientific);
                os.precision(17);
                char tag[256];
                os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << std::endl;
                sprintf(tag,"<archive major_version=\"%d\" minor_version=\"%d\">",
                        ARCHIVE_MAJOR_VERSION, ARCHIVE_MINOR_VERSION);
                os << tag << std::endl;
                os << "<typemap>" << std::endl;
                for (int i=0; i<256; i++) {
                    sprintf(tag,"%d \"%s\"",i,archive_type_names[i]);
                    store(tag,strlen(tag)); // Must use store to escape characters
                }
                os << "</typemap>" << std::endl;
            };


            void close() {
                if (os.is_open()) {
                    os << "</archive>" << std::endl;
                    os.close();
                }
            };


            void flush() {
                if (os.is_open()) os.flush();
            };

            ~TextFstreamOutputArchive() {
                close();
            };
        };


        /// Wraps an archive around a text file stream for input
        class TextFstreamInputArchive : public BaseInputArchive {
        private:
            mutable std::ifstream is;

            // Eat EOL after each entry to enable char-by-char read of strings
            void eat_eol() const {
                char eol;
                is.get(eol);
                if (eol != '\n')
                    MADNESS_EXCEPTION("TextFstreamInputArchive: eat_eol: indigestion", static_cast<int>(eol));
            };

        public:
            TextFstreamInputArchive(const char* filename = 0,
                                    std::ios_base::openmode mode = std::ios_base::in) {
                if (filename) open(filename, mode);
            }

            template <class T>
            void check_start_tag(bool end=false) const {
                char tag[256], ftag[256];
                is.getline(ftag,256);
                unsigned char cookie = archive_typeinfo<T>::cookie;
                if (end)
                    sprintf(tag,"</t%d>",cookie);
                else
                    sprintf(tag,"<t%d>",cookie);

                if (strcmp(tag,ftag) != 0) {
                    std::cout << "TextFstreamInputArchive: type mismatch: expected=" << tag
                              << " "
                              << archive_type_names[cookie]
                              << " "
                              << " got=" << ftag << std::endl;
                    MADNESS_EXCEPTION("TextFstreamInputArchive: check_tag: types do not match/corrupt file", 1);
                }
            }

            template <class T>
            inline
            void check_end_tag() const {
                check_start_tag<T>(true);
            }

            template <class T>
            typename madness::enable_if< madness::is_serializable<T>, void >::type
            load(T* t, long n) const {
                for (long i=0; i<n; i++) is >> t[i];
                eat_eol();
            }

            void load(unsigned char* t, long n) const {
                for (long i=0; i<n; i++) {
                    unsigned int x;
                    is >> x;
                    t[i] = (unsigned char) x;
                }
                eat_eol();
            };

            void load(char* t, long n) const {
                for (long i=0; i<n; i++) {
                    char c0;
                    is.get(c0);
                    if (c0 == '\\') {
                        char c1;
                        is.get(c1);
                        if (c1 == '\\') {
                            t[i] = '\\';
                        }
                        else if (c1 == 'l') {
                            t[i] = '<';
                        }
                        else if (c1 == 'r') {
                            t[i] = '>';
                        }
                        else {
                          MADNESS_EXCEPTION("TextFstreamInputArchive: malformed string?", static_cast<int>(c1));
                        }
                    }
                    else {
                        t[i] = c0;
                    }
                }
                eat_eol();
            };

            void open(const char* filename,
                      std::ios_base::openmode mode = std::ios_base::in) {
                is.open(filename, mode);
                char buf[256];
                is.getline(buf,256);		// skip xml header
                is.getline(buf,256);

                char tag[256];
                sprintf(tag,"<archive major_version=\"%d\" minor_version=\"%d\">",
                        ARCHIVE_MAJOR_VERSION, ARCHIVE_MINOR_VERSION);
                if (strcmp(buf,tag)) {
                    std::cout << "TextFstreamInputArchive: not an archive/bad version?" << std::endl;
                    std::cout << "Found this: " << buf;
                    std::cout << "Expected  : " << tag;
                    MADNESS_EXCEPTION("TextFstreamInputArchive: not an archive/bad version?", 1);
                }

                // For now just skip over typemap
                for (int i=0; i<258; i++) is.getline(buf,256);
            };

            void close() {
                is.close();
            };
        };

        template <class T>
        struct ArchivePrePostImpl<TextFstreamOutputArchive,T> {
            static void preamble_store(const TextFstreamOutputArchive& ar) {
                ar.store_start_tag<T>();
            };
            static inline void postamble_store(const TextFstreamOutputArchive& ar) {
                ar.store_end_tag<T>();
            };
        };

        template <class T>
        struct ArchivePrePostImpl<TextFstreamInputArchive,T> {
            static inline void preamble_load(const TextFstreamInputArchive& ar) {
                ar.check_start_tag<T>();
            };
            static inline void postamble_load(const TextFstreamInputArchive& ar) {
                ar.check_end_tag<T>();
            };
        };
    }
}

#endif
