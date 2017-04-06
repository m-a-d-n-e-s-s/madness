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
 \file madness_exception.h
 \brief Defines \c madness::MadnessException for exception handling.
 \ingroup libraries

 By default, the \c MADNESS_ASSERT macro throws a \c madness::MadnessException.
 Configure options can specify other behaviors: \c MADNESS_ASSERT can be
    - THROW: (default) throw a \c madness::MadnessException.
    - ASSERT: defer to the standard \c assert function.
    - DISABLE: do nothing (ignore the assertion).
    - ABORT: abort execution.
*/

#ifndef MADNESS_WORLD_MADNESS_EXCEPTION_H__INCLUDED
#define MADNESS_WORLD_MADNESS_EXCEPTION_H__INCLUDED

#include <iosfwd>
#include <exception>
#include <madness/madness_config.h>
#ifdef MADNESS_ASSERTIONS_ASSERT
#  include <cassert>
#endif

#ifndef MADNESS_DISPLAY_EXCEPTION_BREAK_MESSAGE
/// Display the exception break message unless otherwise specified.
#define MADNESS_DISPLAY_EXCEPTION_BREAK_MESSAGE 1
#endif

namespace madness {

    /// Base class for exceptions thrown in MADNESS.

    /// Most exceptions thrown in MADNESS should be derived from this.
    class MadnessException : public std::exception {
    public:
        const char* msg; ///< The error message.
        const char* assertion; ///< String describing the assertion.
        const int value; ///< Value associated with the exception.
        const int line; ///< Line number where the exception occurred.
        const char *function; ///< Function where the exception occurred.
        const char *filename; ///< File where the exception occurred.

        /// Constructor that processes the requisite information.

        /// Capturing the line/function/filename info is best done with the
        /// macros listed below.
        /// \param[in] m The error message.
        /// \param[in] a String describing the exception.
        /// \param[in] v Value associated with the exception.
        /// \param[in] l Line number where the exception occurred.
        /// \param[in] fn Function where the exception occurred.
        /// \param[in] f File where the exception occurred.
        MadnessException(const char* m, const char *a, int v,
                         int l, const char *fn, const char *f)
                : msg(m)
                , assertion(a)
                , value(v)
                , line(l)
                , function(fn)
                , filename(f) {}

        /// Returns the error message, as specified by `std::exception`.

        /// \return The error message.
        virtual const char* what() const throw() {
            return msg;
        }
    };

    /// Enables easy printing of a \c MadnessException.

    /// \param[in,out] out Output stream.
    /// \param[in] e The \c MadnessException.
    /// \return The output stream.
    std::ostream& operator <<(std::ostream& out, const MadnessException& e);

    /// This function is executed just before a \c MadnessException is thrown.

    /// \param[in] message True to print an error message to \c cerr; false otherwise.
    void exception_break(bool message);

/// Macro for throwing a MADNESS exception.

/// \throws A \c madness::MadnessException.
/// \param[in] msg The error message.
/// \param[in] value The value associated with the exception.
#define MADNESS_EXCEPTION(msg,value)  { \
    madness::exception_break(MADNESS_DISPLAY_EXCEPTION_BREAK_MESSAGE); \
    throw madness::MadnessException(msg,0,value,__LINE__,__FUNCTION__,__FILE__); \
}

// the following define/undef are for documentation purposes only.
/// Assert a condition.

/// Depending on the configuration, one of the following happens if
/// \c condition is false:
/// - a \c madness::MadnessException is thrown.
/// - `assert(condition)` is called.
/// - execution is aborted.
/// - nothing.
/// \param[in] condition The condition to be asserted.
#define MADNESS_ASSERT(condition)
#undef MADNESS_ASSERT

#ifdef MADNESS_ASSERTIONS_ABORT
#  define MADNESS_ASSERT(condition) \
     do {if (!(condition)) ((void (*)())0)();} while(0)
#endif

#ifdef MADNESS_ASSERTIONS_DISABLE
#  define MADNESS_ASSERT(condition)
#endif

#ifdef MADNESS_ASSERTIONS_ASSERT
#  define MADNESS_ASSERT(condition) assert(condition)
#endif

#ifdef MADNESS_ASSERTIONS_THROW
#  define MADNESS_STRINGIZE(X) #X
#  define MADNESS_EXCEPTION_AT(F, L) MADNESS_STRINGIZE(F) "(" MADNESS_STRINGIZE(L) ")"
#  define MADNESS_ASSERT(condition) \
    do { \
        if (!(condition)) { \
            madness::exception_break(MADNESS_DISPLAY_EXCEPTION_BREAK_MESSAGE); \
            throw madness::MadnessException("MADNESS ASSERTION FAILED: " MADNESS_EXCEPTION_AT( __FILE__, __LINE__ ), \
                #condition,0,__LINE__,__FUNCTION__,__FILE__); \
        } \
    } while (0)
#endif

} // namespace madness

#endif // MADNESS_WORLD_MADNESS_EXCEPTION_H__INCLUDED
