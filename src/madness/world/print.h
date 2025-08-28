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

#ifndef MADNESS_WORLD_PRINT_H__INCLUDED
#define MADNESS_WORLD_PRINT_H__INCLUDED

/**
 \file print.h
 \brief Defines simple templates for printing to \c std::cout "a la Python".
 \ingroup libraries
*/

#include <type_traits>
#include <iostream>
#include <fstream>
#include <complex>
#include <list>
#include <vector>
#include <madness/world/worldmutex.h>
#include <madness/world/array_addons.h>

#ifdef BRAINDEAD
// Cray XT nonsense
#define ENDL "\n"

#else

#define ENDL std::endl

#endif


namespace madness {

namespace operators {

/// \addtogroup libraries
/// @{

/// Easy printing of complex numbers.

/// \tparam T The "real" type of the complex number.
/// \param[in,out] s The output stream.
/// \param[in] c The complex number.
/// \return The output stream (for chaining).
template <typename T>
std::ostream &operator<<(std::ostream &s, const std::complex<T> &c) {
  s << c.real() << "+" << c.imag() << "j";
  return s;
}

/// Easy printing of pairs.

/// \tparam T Type 1 of the pair.
/// \tparam U Type 2 of the pair.
/// \param[in,out] s The output stream.
/// \param[in] p The pair.
/// \return The output stream (for chaining).
template <typename T, typename U>
std::ostream &operator<<(std::ostream &s, const std::pair<T, U> &p) {
  s << "(" << p.first << "," << p.second << ")";
  return s;
}

///// Easy printing of std::arrays.
//
///// \tparam U Type 2 of the pair.
///// \param[in,out] s The output stream.
///// \param[in] p The pair.
///// \return The output stream (for chaining).
//template <typename T, std::size_t NDIM>
//std::ostream &operator<<(std::ostream &s, const std::array<T, NDIM> &p) {
//    s << "[";
//    for (int i=0; i<p.size(); ++i) s << p[i];
//    s << "]";
//    return s;
//}

/// Easy printing of lists.

/// \tparam T Type stored in the list.
/// \param[in,out] s The output stream.
/// \param[in] c The list.
/// \return The output stream (for chaining).
template <typename T>
std::ostream &operator<<(std::ostream &s, const std::list<T> &c) {
  s << "[";
  typename std::list<T>::const_iterator it = c.begin();
  while (it != c.end()) {
    s << *it;
    ++it;
    if (it != c.end())
      s << ", ";
  };
  s << "]";
  return s;
}

/// Easy printing of vectors.

/// \tparam T Type stored in the vector.
/// \param[in,out] s The output stream.
/// \param[in] c The vector.
/// \return The output stream (for chaining).
template <typename T>
std::ostream &operator<<(std::ostream &s, const std::vector<T> &c) {
  s << "[";
  typename std::vector<T>::const_iterator it = c.begin();
  while (it != c.end()) {
    s << *it;
    ++it;
    if (it != c.end())
      s << ", ";
  };
  s << "]";
  return s;
}

/// Easy printing of fixed dimension arrays.

/// STL I/O already does char (thus the \c enable_if business).
/// \tparam T Type of data in the array.
/// \tparam N Size of the array.
/// \param[in,out] s The output stream.
/// \param[in] v The array.
/// \return The output stream (for chaining).
template <typename T, std::size_t N>
typename std::enable_if<!std::is_same<T, char>::value, std::ostream &>::type
operator<<(std::ostream &s, const T (&v)[N]) {
  s << "[";
  for (std::size_t i = 0; i < N; ++i) {
    s << v[i];
    if (i != (N - 1))
      s << ",";
  }
  s << "]";
  return s;
}

}  // namespace operators

    /// big section heading
    void print_header1(const std::string& s);

    /// medium section heading
    void print_header2(const std::string& s);

    /// small section heading
    void print_header3(const std::string& s);

    /// Print a string justified on the left to start at the given column with optional underlining.
    void print_justified(const char* s, int column=0, bool underline=true);

    /// Print a string centered at the given column with optional underlining.
    void print_centered(const char* s, int column=40, bool underline=true);

    ///
    void printf_msg_energy_time(const std::string msg, const double energy, const double time);

    // the "print" function and functions to help handle the variadic templates

    /// Helper function for \c print. Base case.

    /// This gets called recursively when there are no items left to print.
    /// \param[in,out] out Output stream.
    /// \return The output stream (for chaining).
    inline std::ostream& print_helper(std::ostream& out) {
        return out;
    }

    /// \brief Helper function for \c print. Prints the first item (\c t) and
    ///    recursively passes on the other items.

    /// \tparam T Type of the item to print in this call.
    /// \tparam Ts Argument pack type for the remaining items.
    /// \param[in,out] out Output stream.
    /// \param[in] t The item to print in this call.
    /// \param[in] ts The remaining items in the argument pack (they get
    ///    recursively passed on).
    /// \return The output stream (for chaining).
    template <typename T, typename... Ts>
    inline std::ostream& print_helper(std::ostream& out,
                                      const T& t, const Ts&... ts) {
        using madness::operators::operator<<;
        out << ' ' << t;
        return print_helper(out, ts...);
    }

    /// \brief Print items to \c std::cout (items separated by spaces) and
    ///    terminate with a new line

    /// The first item is printed here so that it isn't preceded by a space.
    /// \tparam T Type of the first item to be printed.
    /// \tparam Ts Argument pack type for the items to be printed.
    /// \param[in] t The first item to be printed.
    /// \param[in] ts The remaining items to be printed in the argument pack.
    template<typename T, typename... Ts>
    void print(const T& t, const Ts&... ts) {
        using madness::operators::operator<<;
        ScopedMutex<Mutex> safe(detail::printmutex);
        std::cout << t;
        print_helper(std::cout, ts...) << ENDL;
    }

    /// \brief Print items to \c std::cerr (items separated by spaces) and
    ///    terminate with a new line

    /// The first item is printed here so that it isn't preceded by a space.
    /// \tparam T Type of the first item to be printed.
    /// \tparam Ts Argument pack type for the items to be printed.
    /// \param[in] t The first item to be printed.
    /// \param[in] ts The remaining items to be printed in the argument pack.
    template<typename T, typename... Ts>
    void print_error(const T& t, const Ts&... ts) {
        using madness::operators::operator<<;
        ScopedMutex<Mutex> safe(detail::printmutex);
        std::cerr << t;
        print_helper(std::cerr, ts...) << ENDL;
    }


    /// RAII class to redirect cout to a file
    struct io_redirect {
        std::streambuf* stream_buffer_cout;
        static std::streambuf* stream_buffer_cout_default; ///< default stream buffer for cout, used to restore cout
        std::ofstream ofile;
        bool debug = false;

        io_redirect(const long task_number, std::string filename, bool debug = false) : debug(debug) {
            stream_buffer_cout_default = std::cout.rdbuf();
            constexpr std::size_t bufsize = 256;
            char cfilename[bufsize];
            std::snprintf(cfilename, bufsize, "%s.%5.5ld", filename.c_str(), task_number);
            ofile = std::ofstream(cfilename);
            if (debug) std::cout << "redirecting to file " << cfilename << std::endl;
            stream_buffer_cout = std::cout.rdbuf(ofile.rdbuf());
            std::cout.sync_with_stdio(true);
        }

        ~io_redirect() {
            std::cout.rdbuf(stream_buffer_cout);
            ofile.close();
            std::cout.sync_with_stdio(true);
            if (debug) std::cout << "redirecting back to cout" << std::endl;
        }
    };

    /// class to temporarily redirect output to cout
    struct io_redirect_cout {
        std::streambuf* stream_buffer_cout;

        io_redirect_cout() {
            stream_buffer_cout = std::cout.rdbuf(io_redirect::stream_buffer_cout_default);
            std::cout.sync_with_stdio(true);
        }

        ~io_redirect_cout() {
            std::cout.rdbuf(stream_buffer_cout);
            std::cout.sync_with_stdio(true);
        }
    };

}
#endif // MADNESS_WORLD_PRINT_H__INCLUDED
