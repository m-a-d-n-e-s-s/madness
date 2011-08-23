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

#ifndef MADNESS_WORLD_TR1_ARRAY_H__INCLUDED
#define MADNESS_WORLD_TR1_ARRAY_H__INCLUDED

#include <madness_config.h>

// Select header that contains array
#if defined(MADNESS_USE_ARRAY)
#include <array>
#elif defined(MADNESS_USE_TR1_ARRAY)
#include <tr1/array>
#elif defined(MADNESS_USE_BOOST_TR1_ARRAY_HPP)
#include <boost/tr1/array.hpp>
#else
#error No acceptable include directive for TR1 array was found.
#endif // ARRAY

#ifndef MADNESS_BEGIN_NAMESPACE_TR1

#if defined(BOOST_TR1_ARRAY_INCLUDED) || defined(BOOST_TR1_ARRAY_HPP_INCLUDED)

// We are using boost
#define MADNESS_BEGIN_NAMESPACE_TR1 namespace boost {
#define MADNESS_END_NAMESPACE_TR1 } // namespace boost

#elif defined(MADNESS_HAS_STD_TR1_ARRAY)

// We are using TR1
#define MADNESS_BEGIN_NAMESPACE_TR1 namespace std { namespace tr1 {
#define MADNESS_END_NAMESPACE_TR1 } } // namespace std namespace tr1

#elif defined(MADNESS_HAS_STD_ARRAY)

// We are using C++0x
#define MADNESS_BEGIN_NAMESPACE_TR1 namespace std {
#define MADNESS_END_NAMESPACE_TR1 } // namespace std

#else
// We do not know.
#error Unable to determine the correct namespace for TR1 fuctional.

#endif

#endif // MADNESS_BEGIN_NAMESPACE_TR1

// Insert the tr1 array class into the std namespace.
namespace std {

#if defined(MADNESS_HAS_STD_TR1_ARRAY) && !defined(MADNESS_HAS_STD_ARRAY)
#define MADNESS_HAS_STD_ARRAY 1

    using ::std::tr1::array;
    using ::std::tr1::swap;
    using ::std::tr1::tuple_size;
    using ::std::tr1::tuple_element;
    using ::std::tr1::tuple_size;
    using ::std::tr1::tuple_element;
    using ::std::tr1::get;

#endif
} // namespace std

#endif // MADNESS_WORLD_TR1_ARRAY_H__INCLUDED
