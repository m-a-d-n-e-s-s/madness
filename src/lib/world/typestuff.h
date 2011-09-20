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


#ifndef MADNESS_WORLD_TYPESTUFF_H__INCLUDED
#define MADNESS_WORLD_TYPESTUFF_H__INCLUDED

/// \file typestuff.h
/// \brief Grossly simplified Boost-like type traits and templates

#include <cstddef>
#include <stdint.h>
#include <madness_config.h>

// Next 3 use TR1 or BOOST if available, otherwise a minimal
// version we support ... in the std:: namespace
#include <world/type_traits.h>

// Boost-like enable if stuff in the madness:: namespace
#include <world/enable_if.h>

// Boost-like stuff in the madnes::tr1::details namespace
#include <world/function_traits_bits.h>

namespace madness {

    /// True if A is derived from B and is not B
    template <class A, class B>
    struct is_derived_from {
        typedef char yes;
        typedef int no;
        static no f(...);
        static yes f(B*);
        static const bool value = (sizeof(f((A*)0)) == sizeof(yes));
    };

    /// True if A is derived from B and is not B
    template <class A>
    struct is_derived_from<A,A> : public std::false_type {};

    template <typename> class Future;
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
        static const bool value = std::is_fundamental<T>::value || std::is_member_function_pointer<T>::value || std::is_function<T>::value;
    };

    /// Simple binder for member functions with no arguments
    template <class T, typename resultT>
    class BindNullaryMemFun {
    private:
        T* t;
        resultT(T::*op)();
    public:
        BindNullaryMemFun(T* t, resultT(T::*op)()) : t(t), op(op) {};
        resultT operator()() {
            return (t->*op)();
        };
    };

    /// Specialization of BindNullaryMemFun for void return
    template <class T>
    class BindNullaryMemFun<T,void> {
    private:
        T* t;
        void (T::*op)();
    public:
        BindNullaryMemFun(T* t, void (T::*op)()) : t(t), op(op) {};
        void operator()() {
            (t->*op)();
        };
    };


    /// Simple binder for const member functions with no arguments
    template <class T, typename resultT>
    class BindNullaryConstMemFun {
    private:
        const T* t;
        resultT(T::*op)() const;
    public:
        BindNullaryConstMemFun(const T* t, resultT(T::*op)() const) : t(t), op(op) {};
        resultT operator()() const {
            return (t->*op)();
        };
    };

    /// Specialization of BindNullaryConstMemFun for void return
    template <class T>
    class BindNullaryConstMemFun<T,void> {
    private:
        const T* t;
        void (T::*op)() const;
    public:
        BindNullaryConstMemFun(const T* t, void (T::*op)() const) : t(t), op(op) {};
        void operator()() const {
            (t->*op)();
        };
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


    /// Used to provide rvalue references to support move semantics
    template <typename T> class Reference {
        T* p;
    public:
        Reference(T* x) : p(x) {}
        T& operator*() const {return *p;}
        T* operator->() const {return p;}
    };


    /* Macros to make some of this stuff more readable */

    /**
       \def REMREF(TYPE)
       \brief Macro to make remove_reference<T> easier to use

       \def REMCONST(TYPE)
       \brief Macro to make remove_const<T> easier to use

       \def MEMFUN_RETURNT(TYPE)
       \brief Macro to make member function type traits easier to use

       \def MEMFUN_OBJT(TYPE)
       \brief Macro to make member function type traits easier to use

       \def MEMFUN_ARITY(TYPE)
       \brief Macro to make member function type traits easier to use

       \def MEMFUN_ARG1T(TYPE)
       \brief Macro to make member function type traits easier to use

       \def MEMFUN_ARG2T(TYPE)
       \brief Macro to make member function type traits easier to use

       \def MEMFUN_ARG3T(TYPE)
       \brief Macro to make member function type traits easier to use

       \def MEMFUN_ARG4T(TYPE)
       \brief Macro to make member function type traits easier to use

       \def FUNCTION_RETURNT(TYPE)
       \brief Macro to make function type traits easier to use

       \def FUNCTION_ARITY(TYPE)
       \brief Macro to make function type traits easier to use

       \def FUNCTION_ARG1T(TYPE)
       \brief Macro to make function type traits easier to use

       \def FUNCTION_ARG2T(TYPE)
       \brief Macro to make function type traits easier to use

       \def FUNCTION_ARG3T(TYPE)
       \brief Macro to make function type traits easier to use

       \def FUNCTION_ARG4T(TYPE)
       \brief Macro to make function type traits easier to use

       \def IS_SAME(TYPEA, TYPEB)
       \brief Macro to make is_same<T> template easier to use

       \def RETURN_WRAPPERT(TYPE)
       \brief Returns TYPE unless TYPE is void when returns ReturnWrapper<void>

    */

#define REMREF(TYPE)    typename std::remove_reference< TYPE >::type
#define REMCONST(TYPE)  typename std::remove_const< TYPE >::type
#define REMCONSTX(TYPE) std::remove_const< TYPE >::type
#define RETURN_WRAPPERT(TYPE) typename madness::ReturnWrapper< TYPE >::type

#define MEMFUN_RETURNT(MEMFUN) typename madness::detail::memfunc_traits< MEMFUN >::result_type
#define MEMFUN_CONSTNESS(MEMFUN) madness::detail::memfunc_traits< MEMFUN >::constness
#define MEMFUN_OBJT(MEMFUN)    typename madness::detail::memfunc_traits< MEMFUN >::obj_type
#define MEMFUN_ARITY(MEMFUN)   madness::detail::memfunc_traits< MEMFUN >::arity
#define MEMFUN_ARG1T(MEMFUN)   typename madness::detail::memfunc_traits< MEMFUN >::arg1_type
#define MEMFUN_ARG2T(MEMFUN)   typename madness::detail::memfunc_traits< MEMFUN >::arg2_type
#define MEMFUN_ARG3T(MEMFUN)   typename madness::detail::memfunc_traits< MEMFUN >::arg3_type
#define MEMFUN_ARG4T(MEMFUN)   typename madness::detail::memfunc_traits< MEMFUN >::arg4_type

#define FUNCTION_RETURNT(FUNCTION) typename madness::detail::function_traits< FUNCTION >::result_type
#define FUNCTION_ARITY(FUNCTION)   madness::detail::function_traits< FUNCTION >::arity
#define FUNCTION_ARG1T(FUNCTION)   typename madness::detail::function_traits< FUNCTION >::arg1_type
#define FUNCTION_ARG2T(FUNCTION)   typename madness::detail::function_traits< FUNCTION >::arg2_type
#define FUNCTION_ARG3T(FUNCTION)   typename madness::detail::function_traits< FUNCTION >::arg3_type
#define FUNCTION_ARG4T(FUNCTION)   typename madness::detail::function_traits< FUNCTION >::arg4_type

#define RESULT_OF(FUNCTION) typename madness::detail::result_of< FUNCTION >::type

#define IS_SAME(A, B) std::is_same< A, B >
#define IS_EQ(A, B) std::is_eq< A, B >

} // end of namespace madness
#endif // MADNESS_WORLD_TYPESTUFF_H__INCLUDED
