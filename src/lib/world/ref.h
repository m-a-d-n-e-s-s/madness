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

#ifndef MADNESS_WORLD_REF_H__INCLUDED
#define MADNESS_WORLD_REF_H__INCLUDED

#include <world/worldexc.h>
#include <world/type_traits.h>

//  ref.hpp - ref/cref, useful helper functions
//
//  Copyright (C) 1999, 2000 Jaakko Jarvi (jaakko.jarvi@cs.utu.fi)
//  Copyright (C) 2001, 2002 Peter Dimov
//  Copyright (C) 2002 David Abrahams
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
//  See http://www.boost.org/libs/bind/ref.html for documentation.
//
//  Ported by JAC into madness 6/27/2012

namespace madness {
    namespace detail {

        /// Reference wrapper class

        /// Used to hold references for task functions where a copy would
        /// otherwise be used. This wrapper object is default constructable,
        /// copy constructable,  assignable, and serializable, while a normal
        /// reference would not be. This object can also be used to pass
        /// arguments as a reference to a function that would otherwise be
        /// passed by value. Use the factory functions \c ref and \c cref to
        /// create \c ReferenceWrapper objects.
        /// \tparam T The reference type
        template<class T>
        class ReferenceWrapper {
        private:
            T* t_; ///< A pointer to the object that is being referenced

        public:
            typedef T type; ///< The reference type

            /// Default constructor

            /// The reference references nothing
            ReferenceWrapper() : t_(NULL) { }

            /// Constructor

            /// \param t The object to reference
            explicit ReferenceWrapper(T& t) : t_(&t) { }

            // Compiler generated copy constructor and assignment operator are OK here

            /// Reference accessor

            /// \return A reference to the referenced object
            /// \throw madness::MadnessException Reference has not been set
            T& get() const {
                MADNESS_ASSERT(t_);
                return *t_;
            }

            /// Type conversion operator

            /// This has the same effect as \c get() .
            /// \return A reference to the reference object
            /// \throw madness::MadnessException Reference has not been set
            operator T& () const { return get(); }

            /// Obect pointer accessor

            /// \return A pointer to the referenced object
            /// \throw nothing
            T* get_pointer() const { return t_; }

            /// Serialization

            /// This function is here for compatibility with task functions.
            /// Since serializing a reference to a local object is inheirently
            /// wrong, this function simply throws an exception.
            /// \throw madness::MadnessException Always
            template <typename Archive>
            void serialize(const Archive&) {
                MADNESS_EXCEPTION("ReferenceWrapper serialization not supported.", 0);
            }

        }; // class ReferenceWrapper

    }  // namespace detail


    /// Reference wrapper factory function

    /// \tparam T The reference type (may be const or non-const)
    /// \param t The object to be wrapped
    /// \return A reference wrapper object
    template<class T>
    inline detail::ReferenceWrapper<T> const ref(T& t) {  return detail::ReferenceWrapper<T>(t); }

    /// Constant reference wrapper factory function

    /// \tparam T The reference type (without const)
    /// \param t The object to be wrapped
    /// \return Constant reference wrapper object
    template<class T>
    inline detail::ReferenceWrapper<const T> const cref(const T& t) { return detail::ReferenceWrapper<const T>(t); }


    /// Type trait for reference wrapper

    /// \tparam T The test type
    template<typename T>
    class is_reference_wrapper : public std::false_type { };

    /// \c ReferenceWrapper type trait accessor
    template<typename T>
    class UnwrapReference {
    public:
        typedef T type; ///< The reference type
    }; // class UnwrapReference

    template<typename T>
    class is_reference_wrapper<detail::ReferenceWrapper<T> > : public std::true_type { };

    template<typename T>
    class is_reference_wrapper<detail::ReferenceWrapper<T> const> : public std::true_type { };

    template<typename T>
    class is_reference_wrapper<detail::ReferenceWrapper<T> volatile> : public std::true_type { };

    template<typename T>
    class is_reference_wrapper<detail::ReferenceWrapper<T> const volatile> : public std::true_type { };

    template<typename T>
    class UnwrapReference<detail::ReferenceWrapper<T> > {
    public:
        typedef T type;
    }; // class UnwrapReference<ReferenceWrapper<T> >

    template<typename T>
    class UnwrapReference<detail::ReferenceWrapper<T> const> {
    public:
        typedef T type;
    }; // class UnwrapReference<ReferenceWrapper<T> const>

    template<typename T>
    class UnwrapReference<detail::ReferenceWrapper<T> volatile> {
    public:
        typedef T type;
    }; // class UnwrapReference<ReferenceWrapper<T> volatile>

    template<typename T>
    class UnwrapReference<detail::ReferenceWrapper<T> const volatile> {
    public:
        typedef T type;
    }; // class UnwrapReference<ReferenceWrapper<T> const volatile>


    /// Function for retreaving the referenced object

    /// \tparam The reference type
    /// \param t The reference being unwrapped
    /// \return A reference to the original objects
    template <class T>
    inline typename UnwrapReference<T>::type& unwrap_ref(T& t) { return t; }

    /// Function for retreaving a pointer to the referenced object

    /// \tparam The reference type
    /// \param t The ReferenceWrapper object
    /// \return A reference to the original objects
    template<class T>
    inline T* get_pointer(const detail::ReferenceWrapper<T>& r ) { return r.get_pointer(); }

    // Forward declairation
    template <typename> class Future;

    /// Future for holding \c ReferenceWrapper objects

    /// This is not a standard future. It is intended to hold a \c ReferenceWapper
    /// in a task function. This avoids the coping data into the task.
    /// \tparam T The type held by the reference
    template <typename T>
    class Future<detail::ReferenceWrapper<T> > {
    private:
        detail::ReferenceWrapper<T> ref_;

        // Not needed
        Future(const Future<detail::ReferenceWrapper<T> >&);
        Future<detail::ReferenceWrapper<T> >&
        operator=(const Future<detail::ReferenceWrapper<T> >&);

    public:
        /// Constructor

        /// Construct a future with a reference wrapper object
        Future(const detail::ReferenceWrapper<T>& ref) :
            ref_(ref)
        { }

        // Not needed right now. Uncomment these lines if that changes.
//        Future() : ref_() { }
//        Future(const Future<detail::ReferenceWrapper<T> >& other) :
//            ref_(other.ref_)
//        { }
//        Future<detail::ReferenceWrapper<T> >&
//        operator=(const Future<detail::ReferenceWrapper<T> >& other) {
//            ref_ = other.ref_;
//            return *this;
//        }

        /// Probe set state of the future

        /// \return \c true
        bool probe() const { return true; }

        /// Register a callback

        /// Since this future will never be unset, the \c callback notify function
        /// will always be immediately invocked and return.
        /// \param callback The callback object pointer
        void register_callback(CallbackInterface* callback) {
            callback->notify();
        }

        /// Type convertion to the type being referenced.
        operator T&() const { return ref_.get(); }

    }; // class Future<detail::ReferenceWrapper<T>  >

} // namespace madness

#endif // MADNESS_WORLD_REF_H__INCLUDED
