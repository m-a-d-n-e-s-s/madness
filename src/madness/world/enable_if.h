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

#ifndef MADNESS_WORLD_ENABLE_IF_H__INCLUDED
#define MADNESS_WORLD_ENABLE_IF_H__INCLUDED

namespace madness {

    /// Mirror of \c std::enable_if for conditionally instantiating (disabling) templates based on type.

    /// Evaluates to \c returnT if \c B is false, otherwise to an invalid type expression,
    /// which causes the template expression in which it is used to not be considered for
    /// overload resolution.
    /// \tparam B The bool value.
    /// \tparam returnT The type.
    template <bool B, class returnT = void>
    struct disable_if {
        typedef returnT type;
    };

    /// Mirror of \c std::enable_if for conditionally instantiating (disabling) templates based on type.

    /// Specialization that disables \c type when \c B is true.
    /// \tparam returnT The type.
    template <class returnT>
    struct disable_if<true, returnT> {};

    /// Mirror of \c std::enable_if f for conditionally instantiating templates based on type, when the type \c T may only be meaningful if \c B is true.

    /// Evaluates to \c returnT if \c B is true, otherwise to an invalid type expression
    /// which causes the template expression in which it is used to not be considered for
    /// overload resolution. This "lazy" version is used if \c T is only valid when
    /// B is true. Note: typename T::type is the return type and must be well formed.
    /// \tparam B The bool value.
    /// \tparam returnT The type.
    template <bool B, class returnT>
    struct lazy_enable_if {
      typedef typename returnT::type type;
    };

    /// Mirror of \c std::enable_if for conditionally instantiating templates based on type, when the type \c T may only be meaningful if \c B is true.

    /// Specialization that disables \c type when \c B is false.
    /// \tparam returnT The type.
    template <class returnT>
    struct lazy_enable_if<false, returnT> { };

    /// Mirror of \c madness:lazy_enable_if for conditionally instantiating (disabling) templates based on type, when the type \c T may only be meaningful if \c B is false.

    /// Evaluates to \c returnT if \c B is false, otherwise to an invalid type expression
    /// which causes the template expression in which it is used to not be considered for
    /// overload resolution. This "lazy" version is used if \c returnT is only valid
    /// when B is false. Note: typename T::type is the return type and must be well formed.
    /// \tparam B The bool value.
    /// \tparam returnT The type.
    template <bool B, class returnT>
    struct lazy_disable_if {
      typedef typename returnT::type type;
    };

    /// Mirror of madness::lazy_enable_if for conditionally instantiating (disabling) templates based on type, when the type \c T may only be meaningful if \c B is false.

    /// Specialization that disables \c type when \c B is true.
    /// \tparam returnT The type.
    template <class returnT>
    struct lazy_disable_if<true, returnT> {};

    /// enable_if_same (from Boost?) for conditionally instantiating templates if two types are equal

    /// Use example
    /// \code
    ///     template <class T> A(T& other, typename enable_if_same<A const,T>::type = 0) {
    /// \endcode
//    template <class T, class U, class returnT = void>
//    struct enable_if_same : public enable_if<madness::is_same<T,U>, returnT> {};

        template <bool Cond, typename T1, typename T2>
        struct if_c {
            typedef T1 type;
        };
        
        template <typename T1, typename T2>
        struct if_c<false, T1, T2> {
            typedef T2 type;
        };
        
        template <typename Cond, typename T1, typename T2>
        struct if_ : public if_c<Cond::value, T1, T2> {};



} // namespace madness


/* Macros to make some of this stuff more readable */

/**

   \def ENABLE_IF(CONDITION,TYPEIFTRUE)
   \brief Macro to make enable_if<> template easier to use

   \def DISABLE_IF(CONDITION,TYPEIFTRUE)
   \brief Macro to make disable_if<> template easier to use


   \def DISABLE_IF(CONDITION,TYPEIFTRUE)
   \brief Macro to make enable_if<madness::is_same< A , B > > template easier to use

*/

#define ENABLE_IF(CONDITION,TYPEIFTRUE)  typename madness::enable_if< CONDITION, TYPEIFTRUE >::type
#define DISABLE_IF(CONDITION,TYPEIFTRUE) typename madness::disable_if< CONDITION, TYPEIFTRUE >::type
#define ENABLE_IF_SAME(A,B,TYPEIFTRUE) typename madness::enable_if<madness::is_same< A , B >, TYPEIFTRUE >::type
#define DISABLE_IF_SAME(A,B,TYPEIFTRUE) typename madness::disable_if<madness::is_same< A , B >, TYPEIFTRUE >::type

#endif // MADNESS_WORLD_ENABLE_IF_H__INCLUDED
