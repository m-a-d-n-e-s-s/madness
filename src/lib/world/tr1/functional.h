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


#ifndef MADNESS_WORLD_TR1_FUNCTIONAL_H__INCLUDED
#define MADNESS_WORLD_TR1_FUNCTIONAL_H__INCLUDED

/// \file worldhash.h
/// \brief Include C++ TR1 hash function objects.It imports all TR1 type traits into the
/// std namespace.

#include <madness_config.h>

// Select header that contains hash
#if defined(MADNESS_USE_FUNCTIONAL)
#include <functional>

#elif defined(MADNESS_USE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(MADNESS_USE_BOOST_TR1_FUNCTIONAL_HPP)
#include <boost/tr1/functional.hpp>
#else
#error No acceptable functional include directive was found.
#endif // FUNCTIONAL

#ifndef MADNESS_BEGIN_NAMESPACE_TR1

#if defined(BOOST_TR1_FUNCTIONAL_INCLUDED) || defined(BOOST_TR1_FUNCTIONAL_HPP_INCLUDED)

// We are using boost
#define MADNESS_BEGIN_NAMESPACE_TR1 namespace boost {
#define MADNESS_END_NAMESPACE_TR1 } // namespace std

#elif defined(MADNESS_HAS_STD_TR1_HASH)

// We are using TR1
#define MADNESS_BEGIN_NAMESPACE_TR1 namespace std { namespace tr1 {
#define MADNESS_END_NAMESPACE_TR1 } } // namespace std namespace tr1

#elif defined(MADNESS_HAS_STD_HASH)

// We are using C++0x
#define MADNESS_BEGIN_NAMESPACE_TR1 namespace std {
#define MADNESS_END_NAMESPACE_TR1 } // namespace std

#else
// We do not know.
#error Unable to determine the correct namespace for TR1 fuctional.

#endif

#endif // MADNESS_BEGIN_NAMESPACE_TR1


#if defined(MADNESS_HAS_STD_TR1_HASH) && !defined(MADNESS_HAS_STD_HASH)
#define MADNESS_HAS_STD_HASH 1
// hash is in std::tr1 but we want it in std namespace
namespace std {
    using ::std::tr1::hash;
}
#endif


#endif // MADNESS_WORLD_TR1_FUNCTIONAL_H__INCLUDED
