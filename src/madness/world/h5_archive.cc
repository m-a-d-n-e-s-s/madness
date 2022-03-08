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

/**
 \file text_fstream_archive.cc
 \brief Implements an archive wrapping text filestream.
 \ingroup serialization
*/

#include <madness/world/h5_archive.h>


namespace madness {
    namespace archive {

//vama        void TextFstreamOutputArchive::store(const char* t, long /*n*/) const {
//vama            // store character string, escaping &, < and > along the way
//vama            while (*t) {
//vama                char c = *t++;
//vama                if (c == '\\') {
//vama                    os.put('\\');
//vama                    os.put('\\');
//vama                }
//vama                else if (c == '<') {
//vama                    os.put('\\');
//vama                    os.put('l');
//vama                }
//vama                else if (c == '>') {
//vama                    os.put('\\');
//vama                    os.put('r');
//vama                }
//vama                else {
//vama                    os.put(c);
//vama                }
//vama            }
//vama            os << std::endl;
//vama        }
//vama
//vama        void TextFstreamOutputArchive::open(const char* filename,
//vama                  std::ios_base::openmode mode)
//vama        {
//vama            os.open(filename, mode);
//vama            os.setf(std::ios::scientific);
//vama            os.precision(17);
//vama            char tag[256];
//vama            os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << std::endl;
//vama            sprintf(tag,"<archive major_version=\"%d\" minor_version=\"%d\">",
//vama                    ARCHIVE_MAJOR_VERSION, ARCHIVE_MINOR_VERSION);
//vama            os << tag << std::endl;
//vama            os << "<typemap>" << std::endl;
//vama            for (int i=0; i<256; ++i) {
//vama                sprintf(tag,"%d \"%s\"",i,archive_type_names[i]);
//vama                store(tag,strlen(tag)); // Must use store to escape characters
//vama            }
//vama            os << "</typemap>" << std::endl;
//vama        }
//vama
//vama
//vama        void TextFstreamOutputArchive::close() {
//vama            if (os.is_open()) {
//vama                os << "</archive>" << std::endl;
//vama                os.close();
//vama            }
//vama        }
//vama
//vama
//vama        // Eat EOL after each entry to enable char-by-char read of strings
//vama        void TextFstreamInputArchive::eat_eol() const {
//vama            char eol;
//vama            is.get(eol);
//vama            if (eol != '\n')
//vama                MADNESS_EXCEPTION("TextFstreamInputArchive: eat_eol: indigestion", static_cast<int>(eol));
//vama        }
//vama
//vama
//vama        void TextFstreamInputArchive::load(unsigned char* t, long n) const {
//vama            for (long i=0; i<n; ++i) {
//vama                unsigned int x;
//vama                is >> x;
//vama                t[i] = (unsigned char) x;
//vama            }
//vama            eat_eol();
//vama        }
//vama
//vama        void TextFstreamInputArchive::load(char* t, long n) const {
//vama            for (long i=0; i<n; ++i) {
//vama                char c0;
//vama                is.get(c0);
//vama                if (c0 == '\\') {
//vama                    char c1;
//vama                    is.get(c1);
//vama                    if (c1 == '\\') {
//vama                        t[i] = '\\';
//vama                    }
//vama                    else if (c1 == 'l') {
//vama                        t[i] = '<';
//vama                    }
//vama                    else if (c1 == 'r') {
//vama                        t[i] = '>';
//vama                    }
//vama                    else {
//vama                      MADNESS_EXCEPTION("TextFstreamInputArchive: malformed string?",
//vama                              static_cast<int>(c1));
//vama                    }
//vama                }
//vama                else {
//vama                    t[i] = c0;
//vama                }
//vama            }
//vama            eat_eol();
//vama        }
//vama
//vama        void TextFstreamInputArchive::open(const char* filename, std::ios_base::openmode mode) {
//vama            is.open(filename, mode);
//vama            char buf[256];
//vama            is.getline(buf,256);        // skip xml header
//vama            is.getline(buf,256);
//vama
//vama            char tag[256];
//vama            sprintf(tag,"<archive major_version=\"%d\" minor_version=\"%d\">",
//vama                    ARCHIVE_MAJOR_VERSION, ARCHIVE_MINOR_VERSION);
//vama            if (strcmp(buf,tag)) {
//vama                std::cout << "TextFstreamInputArchive: not an archive/bad version?" << std::endl;
//vama                std::cout << "Found this: " << buf;
//vama                std::cout << "Expected  : " << tag;
//vama                MADNESS_EXCEPTION("TextFstreamInputArchive: not an archive/bad version?", 1);
//vama            }
//vama
//vama            // For now just skip over typemap
//vama            for (int i=0; i<258; ++i) is.getline(buf,256);
//vama        }
//vama
    }
}
