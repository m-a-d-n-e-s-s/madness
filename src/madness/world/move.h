/*
  This file is part of MADNESS.

  Copyright (C) 2012 Justus Calvin

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

#ifndef MADNESS_WORLD_MOVE_H__INCLUDED
#define MADNESS_WORLD_MOVE_H__INCLUDED

namespace madness {
    namespace detail {

        /// Wrapper for movable objects

        /// This object wraps movable objects. It may be used to implement move
        /// constructor and assignment operators. It may be used to set futures.
        /// \tparam T The object type
        template <typename T>
        class MoveWrapper {
            T* t_; ///< A pointer to the wrapped object

        public:
            /// Constructor

            /// \param t The object to be wrapped
            MoveWrapper(T& t) : t_(&t) { }

            /// Copy constructor

            /// \param other The object to be copied
            MoveWrapper(const MoveWrapper<T>& other) : t_(other.t_) { }

            /// Assignment operator

            /// \param other The object to be copied
            /// \return A reference to this object
            MoveWrapper<T>& operator=(const MoveWrapper<T>& other) {
                t_ = other.t_;
                return *this;
            }

            /// Get the wrapped object reference

            /// \return A reference to the wrapped object
            T& get() const { return *t_; }

            /// Get the wrapped object pointer

            /// \return A pointer to the wrapped object
            T* get_pointer() const { return t_; }
        }; // class MoveWrapper

        /// Type trait for movable objects

        /// \tparam T The type to test
        template <typename T>
        struct is_movable : public std::false_type { };

        /// Type trait for movable objects

        /// \tparam T The type to test
        template <typename T>
        struct is_movable<MoveWrapper<T> > : public std::true_type { };

    } // namspace detail

    /// Move wrapper factory function

    /// Construct a move wrapper for a movable object.
    /// \tparam T The wrapped object type
    /// \param t The obect to be wrapped
    /// \return MoveWrapper for t.
    template <typename T>
    detail::MoveWrapper<T> move(T& t) { return detail::MoveWrapper<T>(t); }

    /// Move wrapper factory function

    /// Const objects are not movable so they are not placed in a move wrapper.
    /// \tparam T The object type
    /// \param t The obect
    /// \return A const reference to \c t .
    template <typename T>
    const T& move(const T& t) { return t; }

    /// Remove move wrapper from a movable object

    /// \tparam T The wrapped object type
    /// \param t The object move warpper
    /// \return A reference to the wrapped object
    template <typename T>
    T& unwrap_move(const detail::MoveWrapper<T>& t) { return t.get(); }

    /// Passthrough function for a non-movable object

    /// \tparam T The object type
    /// \param t The objecmove warpper
    /// \return A reference to the object
    template <typename T>
    typename disable_if<detail::is_movable<T>, T&>::type unwrap_move(T& t) { return t; }

} // namespace madness

#endif // MADNESS_WORLD_MOVE_H__INCLUDED
