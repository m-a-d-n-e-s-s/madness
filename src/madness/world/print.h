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
#include <complex>
#include <list>
#include <vector>
#include <madness/world/worldmutex.h>

#ifdef BRAINDEAD
// Cray XT nonsense
#define ENDL "\n"

#else

#define ENDL std::endl

#endif


namespace madness {

    /// \addtogroup libraries
    /// @{

    /// Easy printing of complex numbers.

    /// \tparam T The "real" type of the complex number.
    /// \param[in,out] s The output stream.
    /// \param[in] c The complex number.
    /// \return The output stream (for chaining).
    template <typename T>
    std::ostream& operator<<(std::ostream& s, const std::complex<T>& c) {
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
    std::ostream& operator<<(std::ostream& s, const std::pair<T,U>& p) {
        s << "(" << p.first << "," << p.second << ")";
        return s;
    }

    /// Easy printing of lists.

    /// \tparam T Type stored in the list.
    /// \param[in,out] s The output stream.
    /// \param[in] c The list.
    /// \return The output stream (for chaining).
    template <typename T>
    std::ostream& operator<<(std::ostream& s, const std::list<T>& c) {
        s << "[";
        typename std::list<T>::const_iterator it = c.begin();
        while (it != c.end()) {
            s << *it;
            ++it;
            if (it != c.end()) s << ", ";
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
    std::ostream& operator<<(std::ostream& s, const std::vector<T>& c) {
        s << "[";
        typename std::vector<T>::const_iterator it = c.begin();
        while (it != c.end()) {
            s << *it;
            ++it;
            if (it != c.end()) s << ", ";
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
    typename std::enable_if<!std::is_same<T,char>::value, std::ostream&>::type
    operator<<(std::ostream& s, const T(&v)[N]) {
        s << "[";
        for (std::size_t i=0; i<N; ++i) {
            s << v[i];
            if (i != (N-1)) s << ",";
        }
        s << "]";
        return s;
    }

    /// Print a string justified on the left to start at the given column with optional underlining.
    void print_justified(const char* s, int column=0, bool underline=true);

    /// Print a string centered at the given column with optional underlining.
    void print_centered(const char* s, int column=40, bool underline=true);


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
        ScopedMutex<Mutex> safe(detail::printmutex);
        std::cerr << t;
        print_helper(std::cerr, ts...) << ENDL;
    }

    /// @}

}
#endif // MADNESS_WORLD_PRINT_H__INCLUDED
