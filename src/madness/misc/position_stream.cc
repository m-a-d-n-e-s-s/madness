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
#include <madness/misc/misc.h>
#include <sstream>
#include <stdexcept>

namespace madness {

    /// message thrown when a tag cannot be found

    /// Callers classify the exception by its message prefix
    /// (QCCalculationParametersBase::is_missing_datagroup_exception, to tell a
    /// missing data group from a malformed one), so this text is load-bearing --
    /// keep it in sync with the prefix there. The tag itself is deliberately not
    /// interpolated: the caller knows the tag it asked for anyway.
    static constexpr const char* not_found_msg = "position_stream: failed to locate the requested tag";

    std::istream& position_stream(std::istream& f, const std::string& tag, bool rewind) {
        if (rewind) f.seekg(0);
        std::string s;
        while (std::getline(f,s)) {
            std::string::size_type loc = s.find(tag, 0);
            if(loc != std::string::npos) return f;
        }
        throw std::runtime_error(not_found_msg);
    }

    /// position the input stream to tag, which must be a word (not part of a word)

    /// \param f        input stream
    /// \param tag      the word to look for
    /// \param comment  a comment character for (parts of) a line
    /// \param rewind   rewind to the beginning of the stream
    /// \param silent   throws if not successful, but doesn't print error message
    /// \return         a stream
    std::istream& position_stream_to_word(std::istream& f, const std::string& tag, const char comment, bool rewind,
                                          bool silent) {
        if (rewind) f.seekg(0);
        std::string line, word;
        while (std::getline(f,line)) {
            // remove comments from line
            std::size_t last = line.find_first_of(comment);
            std::string line1=line.substr(0,last);

            // check for tag in line
            std::stringstream sline(line1);
            while (sline >> word) {
                if (tag==word) {        // found tag as a full word in line
                    std::string::size_type loc = line.find(tag, 0);
                    if (loc != std::string::npos) return f;
                }
            }
        }

        if (silent) {
            throw std::runtime_error(not_found_msg);
        } else {
            printf("position_stream: failed to locate %s\n",tag.c_str());
            throw std::runtime_error(not_found_msg);
        }
        return f;
    }


    std::string lowercase(const std::string& s) {
        std::string r(s);
        for (unsigned int i=0; i<r.size(); ++i) r[i] = tolower(r[i]);
        return r;
    }

}

