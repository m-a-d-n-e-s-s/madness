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


  $Id: typestuff.h 2446 2011-07-22 15:21:22Z justus.c79@gmail.com $
*/

#include <madness_config.h>

#if defined(MADNESS_USE_TYPE_TRAITS)
#  include <type_traits>
#elif defined(MADNESS_USE_TR1_TYPE_TRAITS)
#  include <tr1/type_traits>
#elif defined(MADNESS_USE_BOOST_TR1_TYPE_TRAITS_HPP)
#  include <boost/tr1/type_traits.hpp>
#else
#  define MADNESS_HAS_STD_TYPE_TRAITS 1
#  include <world/type_traits_bits.h>
   namespace std {
      using madness::tr1::is_void;
      using madness::tr1::is_integral;
      using madness::tr1::is_floating_point;
      using madness::tr1::is_array;
      using madness::tr1::is_pointer;
      using madness::tr1::is_reference;
      using madness::tr1::is_function;
      using madness::tr1::is_member_function_pointer;
      using madness::tr1::is_arithmetic;
      using madness::tr1::is_fundamental;
      using madness::tr1::is_object;
      using madness::tr1::is_compound;
      using madness::tr1::is_const;

      using madness::tr1::is_same;
      using madness::tr1::is_base_of;

      using madness::tr1::remove_const;
      using madness::tr1::remove_volatile;
      using madness::tr1::remove_cv;
      using madness::tr1::add_const;
      using madness::tr1::remove_reference;
      using madness::tr1::remove_pointer;
      using madness::tr1::enable_if;
   }
#endif

#if defined(MADNESS_HAS_STD_TR1_TYPE_TRAITS) && !defined(MADNESS_HAS_STD_TYPE_TRAITS)
#define MADNESS_HAS_STD_TYPE_TRAITS 1

// Insert the tr1 type traits into the std namespace.
namespace std {
    using ::std::tr1::integral_constant;
    using ::std::tr1::true_type;
    using ::std::tr1::false_type;
    using ::std::tr1::is_void;
    using ::std::tr1::is_integral;
    using ::std::tr1::is_floating_point;
    using ::std::tr1::is_array;
    using ::std::tr1::is_pointer;
    using ::std::tr1::is_reference;
    using ::std::tr1::is_member_object_pointer;
    using ::std::tr1::is_member_function_pointer;
    using ::std::tr1::is_enum;
    using ::std::tr1::is_union;
    using ::std::tr1::is_class;
    using ::std::tr1::is_function;
    using ::std::tr1::is_arithmetic;
    using ::std::tr1::is_fundamental;
    using ::std::tr1::is_object;
    using ::std::tr1::is_scalar;
    using ::std::tr1::is_compound;
    using ::std::tr1::is_member_pointer;
    using ::std::tr1::is_const;
    using ::std::tr1::is_volatile;
    using ::std::tr1::is_pod;
    using ::std::tr1::is_empty;
    using ::std::tr1::is_polymorphic;
    using ::std::tr1::is_abstract;
    using ::std::tr1::has_trivial_constructor;
    using ::std::tr1::has_trivial_copy;
    using ::std::tr1::has_trivial_assign;
    using ::std::tr1::has_trivial_destructor;
    using ::std::tr1::has_nothrow_constructor;
    using ::std::tr1::has_nothrow_copy;
    using ::std::tr1::has_nothrow_assign;
    using ::std::tr1::has_virtual_destructor;
    using ::std::tr1::is_signed;
    using ::std::tr1::is_unsigned;
    using ::std::tr1::alignment_of;
    using ::std::tr1::rank;
    using ::std::tr1::extent;
    using ::std::tr1::is_same;
    using ::std::tr1::is_base_of;
    using ::std::tr1::is_convertible;
    using ::std::tr1::remove_const;
    using ::std::tr1::remove_volatile;
    using ::std::tr1::remove_cv;
    using ::std::tr1::add_const;
    using ::std::tr1::add_volatile;
    using ::std::tr1::add_cv;
    using ::std::tr1::remove_reference;
    using ::std::tr1::add_reference;
    using ::std::tr1::remove_extent;
    using ::std::tr1::remove_all_extents;
    using ::std::tr1::remove_pointer;
    using ::std::tr1::add_pointer;
    using ::std::tr1::aligned_storage;
} // namespace std
#endif
