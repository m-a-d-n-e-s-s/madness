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
        typedef returnT type; ///< The type.
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
      typedef typename returnT::type type; ///< The type.
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
      typedef typename returnT::type type; ///< The type.
    };

    /// Mirror of madness::lazy_enable_if for conditionally instantiating (disabling) templates based on type, when the type \c T may only be meaningful if \c B is false.

    /// Specialization that disables \c type when \c B is true.
    /// \tparam returnT The type.
    template <class returnT>
    struct lazy_disable_if<true, returnT> {};

    /// Use \c Cond to determine the type, \c T1 or \c T2.

    /// \c type will have type \c T1 if \c Cond is true; otherwise it will
    /// have type \c T2.
    /// \tparam Cond The bool value.
    /// \tparam T1 Type of \c type if \c Cond is true.
    /// \tparam T2 Type of \c type if \c Cond is false.
    template <bool Cond, typename T1, typename T2>
    struct switch_type {
        typedef T1 type; ///< The type.
    };
    
    /// Specialization of \c switch_type for when \c Cond is false.

    /// \c type will have type \c T2.
    /// \tparam T1 Type of \c type if \c Cond is true (not used).
    /// \tparam T2 Type of \c type if \c Cond is false.
    template <typename T1, typename T2>
    struct switch_type<false, T1, T2> {
        typedef T2 type; ///< The type.
    };

} // namespace madness

#endif // MADNESS_WORLD_ENABLE_IF_H__INCLUDED
