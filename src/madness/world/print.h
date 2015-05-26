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

/// \file print.h
/// \brief Defines simple templates for printing to \c std::cout "a la Python".

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

    /// Easy printing of complex numbers
    template <typename T>
    std::ostream& operator<<(std::ostream& s, const std::complex<T>& c) {
        s << c.real() << "+" << c.imag() << "j";
        return s;
    }

    /// Easy printing of pairs
    template <typename T, typename U>
    std::ostream& operator<<(std::ostream& s, const std::pair<T,U>& p) {
        s << "(" << p.first << "," << p.second << ")";
        return s;
    }

    /// Easy printing of lists
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

    /// Easy printing of vectors
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

    /// Easy printing of fixed dimension arrays

    /// STL I/O already does char.
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

    /// Print a string justified on the left to start at the given column with optional underlining
    void print_justified(const char* s, int column=0, bool underline=true);

    /// Print a string centered at the given column with optional underlining
    void print_centered(const char* s, int column=40, bool underline=true);



    // the "print" function and structs to help handle the variadic templates

    /// Helper class to print an argument pack.

    /// This is just a forward declaration for the compiler.
    /// \tparam Ts The argument pack type of items to print.
    template<typename... Ts>
    struct print_helper;

    /// Helper class to print an argument pack.

    /// Recursively calls itself (with one less argument in the pack)
    /// after printing the first argument in the pack.
    /// \tparam T Type of the first item (the one that gets printed).
    /// \tparam Ts Argument pack type for the remaining items.
    template<typename T, typename... Ts>
    struct print_helper<T, Ts...> {
        /// \brief Print the first item from the argument pack; recursively
        ///    pass the remaining items on.

        /// \param[in] t The first item in the argument pack.
        /// \param[in] ts The remaining items in the argument pack (they get
        ///    recursively passed on).
        static inline
        void print_nomutex(const T& t, const Ts&... ts) {
            std::cout << ' ' << t; // print the first item
            print_helper<Ts...>::print_nomutex(ts...); // print the remaining items
        }
    };

    /// Helper class to print an argument pack.

    /// This is the base case (no arguments remaining); do nothing.
    template<>
    struct print_helper<> {
        /// Print the... there aren't any arguments left. Do nothing.
        static inline
        void print_nomutex() {
        }
    };

    /// \brief Print items to \c std::cout (items separated by spaces) and
    ///    terminate with a new line

    /// The first item is printed here so that it isn't preceded by a space.
    /// \tparam T The type of the first item to print.
    /// \tparam Ts Argument pack type for the remaining items.
    template<typename T, typename... Ts>
    void print(const T &t, const Ts&... ts) {
        ScopedMutex<Mutex> safe(detail::printmutex);
        std::cout << t;
        print_helper<Ts...>::print_nomutex(ts...);
        std::cout << ENDL;
    }

}
#endif // MADNESS_WORLD_PRINT_H__INCLUDED
