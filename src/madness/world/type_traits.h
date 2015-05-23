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

#ifndef MADNESS_WORLD_TYPE_TRAITS_H__INCLUDED
#define MADNESS_WORLD_TYPE_TRAITS_H__INCLUDED

#include <madness/madness_config.h>
#  include <type_traits>

/// \file typestuff.h
/// \brief Grossly simplified Boost-like type traits and templates

#include <cstddef>
#include <stdint.h>
#include <madness/madness_config.h>
#include <madness/world/type_traits.h>
#include <madness/world/function_traits.h>

namespace madness {

    template <typename> class Future;
    template <typename> struct add_future;
    template <typename> struct remove_future;

    // Remove Future, const, volatile, and reference qualifiers from the type
    template <typename T>
    struct remove_fcvr {
        typedef typename remove_future<typename std::remove_cv<
                                           typename std::remove_reference<T>::type>::type>::type type;
    };

    /// This defines stuff that is serialiable by default rules ... basically anything contiguous
    template <typename T>
    struct is_serializable {
        static const bool value = std::is_fundamental<T>::value || std::is_member_function_pointer<T>::value || std::is_function<T>::value  || std::is_function<typename std::remove_pointer<T>::type>::value;
    };

    /// Simple binder for member functions with no arguments
    template <class T, typename resultT>
    class BindNullaryMemFun {
    private:
        T* t_;
        resultT(T::*op_)();
    public:
        BindNullaryMemFun(T* t, resultT(T::*op)()) : t_(t), op_(op) {}
        resultT operator()() {
            return (t_->*op_)();
        }
    };

    /// Specialization of BindNullaryMemFun for void return
    template <class T>
    class BindNullaryMemFun<T,void> {
    private:
        T* t_;
        void (T::*op_)();
    public:
        BindNullaryMemFun(T* t, void (T::*op)()) : t_(t), op_(op) {}
        void operator()() {
            (t_->*op_)();
        }
    };


    /// Simple binder for const member functions with no arguments
    template <class T, typename resultT>
    class BindNullaryConstMemFun {
    private:
        const T* t_;
        resultT(T::*op_)() const;
    public:
        BindNullaryConstMemFun(const T* t, resultT(T::*op)() const) : t_(t), op_(op) {}
        resultT operator()() const {
            return (t_->*op_)();
        }
    };

    /// Specialization of BindNullaryConstMemFun for void return
    template <class T>
    class BindNullaryConstMemFun<T,void> {
    private:
        const T* t_;
        void (T::*op_)() const;
    public:
        BindNullaryConstMemFun(const T* t, void (T::*op)() const) : t_(t), op_(op) {}
        void operator()() const {
            (t_->*op_)();
        }
    };

    /// Factory function for BindNullaryMemFun
    template <class T, typename resultT>
    inline BindNullaryMemFun<T,resultT> bind_nullary_mem_fun(T* t, resultT(T::*op)()) {
        return BindNullaryMemFun<T,resultT>(t,op);
    }

    /// Factory function for BindNullaryConstMemFun
    template <class T, typename resultT>
    inline BindNullaryConstMemFun<T,resultT> bind_nullary_mem_fun(const T* t, resultT(T::*op)() const) {
        return BindNullaryConstMemFun<T,resultT>(t,op);
    }

    /// A type you can return when you want to return void ... use "return None"
    struct Void {};

    /// None, a la Python
    static const Void None = Void();

    /// Wrapper so that can return something even if returning void
    template <typename T>
    struct ReturnWrapper {
        typedef T type;
    };

    /// Wrapper so that can return something even if returning void
    template <>
    struct ReturnWrapper<void> {
        typedef Void type;
    };


    template <typename T>
    struct has_result_type {
        // yes and no are guaranteed to have different sizes,
        // specifically sizeof(yes) == 1 and sizeof(no) == 2
        typedef char yes[1];
        typedef char no[2];

        template <typename C>
        static yes& test(typename C::result_type*);

        template <typename>
        static no& test(...);

        // if the sizeof the result of calling test<T>(0) is equal to the sizeof(yes),
        // the first overload worked and T has a nested type named type.
        static const bool value = sizeof(test<T>(0)) == sizeof(yes);
    };


    /* Macros to make some of this stuff more readable */

    /**
       \def REMCONST(TYPE)
       \brief Macro to make remove_const<T> easier to use

       \def MEMFUN_RETURNT(TYPE)
       \brief Macro to make member function type traits easier to use
    */

#define REMCONST(TYPE)  typename std::remove_const< TYPE >::type
#define MEMFUN_RETURNT(MEMFUN) typename madness::detail::memfunc_traits< MEMFUN >::result_type

} // namespace madness

#endif // MADNESS_WORLD_TYPE_TRAITS_H__INCLUDED
