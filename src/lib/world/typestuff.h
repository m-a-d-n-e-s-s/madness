/*
  This file is part of MADNESS.
  
  Copyright (C) <2007> <Oak Ridge National Laboratory>
  
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

  
#ifndef TYPESTUFF_H
#define TYPESTUFF_H

/// \file typestuff.h
/// \brief Grossly simplified Boost-like type traits and templates

// Aims to be a compatible subset of Boost ... wish we could rely on
// Boost being there but between bjam (theirs), bloat (theirs), brain
// deficiences (mine, not theirs), bleeding edge massively parallel
// systems, and compiler deficiences, we have to be independent for the
// time being.

#include <stdint.h>

namespace madness {
    
    /// type_or_c<bool A, bool B>::value will be true if (A || B)
    template <bool A, bool B>
    class type_or_c {
    public: 
        static const bool value = true;
    };
    
    /// type_or_c<bool A, bool B>::value will be true if (A || B)
    template<> class type_or_c<false,false> {
    public: 
        static const bool value = false; 
    };

    /// type_or<CondA,CondB>::value  will be true if (CondA::value || CondB::value)
    template <class CondA, class CondB>
    class type_or : public type_or_c<CondA::value, CondB::value> {};


    /// type_and_c<bool A, bool B>::value will be true if (A && B)
    template <bool A, bool B>
    class type_and_c {
    public: 
        static const bool value = false;
    };
    
    /// type_and_c<bool A, bool B>::value will be true if (A && B)
    template<> class type_and_c<true, true> {
    public: 
        static const bool value = true; 
    };

    /// type_and<CondA,CondB>::value  will be true if (CondA::value && CondB::value)
    template <class CondA, class CondB>
    class type_and: public type_and_c<CondA::value, CondB::value> {};
    
    /// is_integral<T>::value will be true if T is a fundamental integer type
    template <class T>
    class is_integral {
    public: 
        static const bool value = false;
    };
    
    /// is_float<T>::value will be true if T is a fundamental floating-point type
    template <class T>
    class is_float {
    public: 
        static const bool value = false;
    };
    
    /// is_fundamental<T>::value will be true if T is a fundamental type 
    template <class T>
    class is_fundamental {
    public: 
        static const bool value = type_or< is_integral<T>, is_float<T> >::value;  
    };

    
    /// is_pointer<T>::value will be true if T is a pointer type
    template <typename T> 
    struct is_pointer 
    { static const bool value = false; };
    
    /// is_pointer<T>::value will be true if T is a pointer type
    template <typename T> 
    struct is_pointer<T*> 
    { static const bool value = true; };
    
    /// is_array<T>::value will be true if T is an array type
    template <typename T> 
    struct is_array
    { static const bool value = false; };
    
    /// is_array<T>::value will be true if T is an array type
    template <typename T, std::size_t n> 
    struct is_array<T[n]>
    { static const bool value = true; };


    /// is_reference<T>::value will be true if T is a reference type
    template <typename T>
    struct is_reference {
        static const bool value = false;
    };
    
    /// is_reference<T>::value will be true if T is a reference type
    template <typename T>
    struct is_reference<T&> {
        static const bool value = true;
    };

    /// is_const<T>::value will be true if T is const
    template <typename T>
    struct is_const {
        static const bool value = false;
    };
    
    /// is_const<T>::value will be true if T is const
    template <typename T>
    struct is_const<const T> {
        static const bool value = true;
    };


    /// is_same<A,B> returns true if A and B are the same type
    template <typename A, typename B>
    struct is_same {
        static const bool value = false;
    };
    
    /// is_same<A,B> returns true if A and B are the same type
    template <typename A>
    struct is_same<A,A> {
        static const bool value = true;
    };


    /// remove_reference<&T>::type will be T
    template <typename T>
    struct remove_reference {
        typedef T type;
    };
    
    /// remove_reference<&T>::type will be T
    template <typename T>
    struct remove_reference<T&> {
        typedef T type;
    };


    /// remove_const<const T>::type will be T
    template <typename T>
    struct remove_const {
        typedef T type;
    };
    
    /// remove_const<const T>::type will be T
    template <typename T>
    struct remove_const<const T> {
        typedef T type;
    };
    
    /// enable_if_c from Boost for conditionally instantiating templates based on type

    /// Evaluates to \c returnT if \c B is true, otherwise to an invalid type expression
    /// which causes the template expression in which it is used to not be considered for
    /// overload resolution.
    template <bool B, class returnT>
    struct enable_if_c {
        typedef returnT type;
    };

    /// enable_if_c from Boost for conditionally instantiating templates based on type
    template <class returnT>
    struct enable_if_c<false, returnT> {};

    /// enable_if from Boost for conditionally instantiating templates based on type

    /// Evaluates to \c returnT if \c Cond::value is true, otherwise to an invalid type expression
    /// which causes the template expression in which it is used to not be considered for
    /// overload resolution.
    template <class Cond, class returnT>
    struct enable_if : public enable_if_c<Cond::value, returnT> {};

    /// disable_if from Boost for conditionally instantiating templates based on type

    /// Evaluates to \c returnT if \c Cond::value is false, otherwise to an invalid type expression
    /// which causes the template expression in which it is used to not be considered for
    /// overload resolution.
    template <bool B, class returnT>
    struct disable_if_c {
        typedef returnT type;
    };

    /// disable_if from Boost for conditionally instantiating templates based on type
    template <class returnT>
    struct disable_if_c<true, returnT> {};

    /// disable_if from Boost for conditionally instantiating templates based on type
    template <class Cond, class returnT>
    struct disable_if : public disable_if_c<Cond::value, returnT> {};

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
    struct is_derived_from<A,A> {
        static const bool value = false;
    };


    /// Function traits in the spirt of boost function traits
    template <typename functionT> 
    struct function_traits {
        static const bool value = false;
    };
    
    /// Function traits in the spirt of boost function traits
    template <typename returnT>
    struct function_traits<returnT (*)()> {
        static const bool value = true;
        static const int arity = 0;
        typedef returnT result_type;
    };
    
    /// Function traits in the spirt of boost function traits
    template <typename returnT, typename arg1T>
    struct function_traits<returnT (*)(arg1T)> {
        static const bool value = true;
        static const int arity = 1;
        typedef returnT result_type;
        typedef arg1T arg1_type;
    };
    
    /// Function traits in the spirt of boost function traits
    template <typename returnT, typename arg1T, typename arg2T>
    struct function_traits<returnT (*)(arg1T,arg2T)> {
        static const bool value = true;
        static const int arity = 2;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
    };
    
    /// Function traits in the spirt of boost function traits
    template <typename returnT, typename arg1T, typename arg2T, typename arg3T>
    struct function_traits<returnT (*)(arg1T,arg2T,arg3T)> {
        static const bool value = true;
        static const int arity = 3;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
    };
    
    /// Function traits in the spirt of boost function traits
    template <typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
    struct function_traits<returnT (*)(arg1T,arg2T,arg3T,arg4T)> {
        static const bool value = true;
        static const int arity = 4;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
    };

    /// Function traits in the spirt of boost function traits
    template <typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
    struct function_traits<returnT (*)(arg1T,arg2T,arg3T,arg4T,arg5T)> {
        static const bool value = true;
        static const int arity = 5;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
    };


    /// Function traits in the spirt of boost function traits
    template <typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
    struct function_traits<returnT (*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T)> {
        static const bool value = true;
        static const int arity = 6;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
    };


    /// Function traits in the spirt of boost function traits
    template <typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
    struct function_traits<returnT (*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T)> {
        static const bool value = true;
        static const int arity = 7;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
        typedef arg7T arg7_type;
    };

    /// Function traits in the spirt of boost function traits
    template <typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T, typename arg8T>
    struct function_traits<returnT (*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T)> {
        static const bool value = true;
        static const int arity = 8;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
        typedef arg7T arg7_type;
        typedef arg8T arg8_type;
    };

    /// Function traits in the spirt of boost function traits
    template <typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T, typename arg8T, typename arg9T>
    struct function_traits<returnT (*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T,arg9T)> {
        static const bool value = true;
        static const int arity = 9;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
        typedef arg7T arg7_type;
        typedef arg8T arg8_type;
        typedef arg9T arg9_type;
    };


    /// Member function traits in the spirt of boost function traits
    template <typename memfuncT> 
    struct memfunc_traits {
        static const bool value = false;
    };
    
    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT>
    struct memfunc_traits<returnT (objT::*)()> {
        static const bool value = true;
        static const int arity = 0;
        static const bool constness = false;
        typedef objT obj_type;
        typedef returnT result_type;
    };
    
    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T>
    struct memfunc_traits<returnT (objT::*)(arg1T)> {
        static const bool value = true;
        static const int arity = 1;
        static const bool constness = false;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
    };
    
    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T>
    struct memfunc_traits<returnT (objT::*)(arg1T,arg2T)> {
        static const bool value = true;
        static const int arity = 2;
        static const bool constness = false;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
    };
    
    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T>
    struct memfunc_traits<returnT (objT::*)(arg1T,arg2T,arg3T)> {
        static const bool value = true;
        static const int arity = 3;
        static const bool constness = false;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
    };
    
    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
    struct memfunc_traits<returnT (objT::*)(arg1T,arg2T,arg3T,arg4T)> {
        static const bool value = true;
        static const int arity = 4;
        static const bool constness = false;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
    };
    

    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
    struct memfunc_traits<returnT (objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T)> {
        static const bool value = true;
        static const int arity = 5;
        static const bool constness = false;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
    };
    

    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
    struct memfunc_traits<returnT (objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T)> {
        static const bool value = true;
        static const int arity = 6;
        static const bool constness = false;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
    };
    

    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
    struct memfunc_traits<returnT (objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T)> {
        static const bool value = true;
        static const int arity = 7;
        static const bool constness = false;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
        typedef arg7T arg7_type;
    };
    

    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T, typename arg8T>
    struct memfunc_traits<returnT (objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T)> {
        static const bool value = true;
        static const int arity = 8;
        static const bool constness = false;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
        typedef arg7T arg7_type;
        typedef arg8T arg8_type;
    };
    

    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T, typename arg8T, typename arg9T>
    struct memfunc_traits<returnT (objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T,arg9T)> {
        static const bool value = true;
        static const int arity = 9;
        static const bool constness = false;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
        typedef arg7T arg7_type;
        typedef arg8T arg8_type;
        typedef arg9T arg9_type;
    };
    

    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT>
    struct memfunc_traits<returnT (objT::*)() const> {
        static const bool value = true;
        static const int arity = 0;
        static const bool constness = true;
        typedef objT obj_type;
        typedef returnT result_type;
    };
    
    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T>
    struct memfunc_traits<returnT (objT::*)(arg1T) const> {
        static const bool value = true;
        static const int arity = 1;
        static const bool constness = true;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
    };
    
    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T>
    struct memfunc_traits<returnT (objT::*)(arg1T,arg2T) const> {
        static const bool value = true;
        static const int arity = 2;
        static const bool constness = true;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
    };
    
    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T>
    struct memfunc_traits<returnT (objT::*)(arg1T,arg2T,arg3T) const> {
        static const bool value = true;
        static const int arity = 3;
        static const bool constness = true;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
    };
    

    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
    struct memfunc_traits<returnT (objT::*)(arg1T,arg2T,arg3T,arg4T) const> {
        static const bool value = true;
        static const int arity = 4;
        static const bool constness = true;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
    };
    


    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
    struct memfunc_traits<returnT (objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T) const> {
        static const bool value = true;
        static const int arity = 5;
        static const bool constness = true;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
    };


    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T>
    struct memfunc_traits<returnT (objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T) const> {
        static const bool value = true;
        static const int arity = 6;
        static const bool constness = true;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
    };


    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T>
    struct memfunc_traits<returnT (objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T) const> {
        static const bool value = true;
        static const int arity = 7;
        static const bool constness = true;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
        typedef arg7T arg7_type;
    };



    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T, typename arg8T>
    struct memfunc_traits<returnT (objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T) const> {
        static const bool value = true;
        static const int arity = 8;
        static const bool constness = true;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
        typedef arg7T arg7_type;
        typedef arg8T arg8_type;
    };


    /// Member function traits in the spirt of boost function traits
    template <typename objT, typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T, typename arg6T, typename arg7T, typename arg8T, typename arg9T>
    struct memfunc_traits<returnT (objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T,arg9T) const> {
        static const bool value = true;
        static const int arity = 9;
        static const bool constness = true;
        typedef objT obj_type;
        typedef returnT result_type;
        typedef arg1T arg1_type;
        typedef arg2T arg2_type;
        typedef arg3T arg3_type;
        typedef arg4T arg4_type;
        typedef arg5T arg5_type;
        typedef arg6T arg6_type;
        typedef arg7T arg7_type;
        typedef arg8T arg8_type;
        typedef arg9T arg9_type;
    };


    /// is_function_ptr<T>::value is true if T is a function pointer
    template <typename T>
    struct is_function_ptr {
        static const bool value = function_traits<T>::value;
    };


    /// is_memfunc_ptr<T>::value is true if T is a member function pointer
    template <typename T>
    struct is_memfunc_ptr {
        static const bool value = memfunc_traits<T>::value;
    };


    /// This defines stuff that is serialiable by default rules ... basically anything contiguous
    template <typename T>
    struct is_serializable {
        static const bool value = madness::is_fundamental<T>::value || madness::is_memfunc_ptr<T>::value || madness::is_function_ptr<T>::value;
    };
        


    

#define SET_TYPE_TRAIT(trait, T, val) \
 template<> class trait < T > {public: static const bool value = val; }
    
    SET_TYPE_TRAIT(is_integral,unsigned char,true);
    SET_TYPE_TRAIT(is_integral,unsigned short,true);
    SET_TYPE_TRAIT(is_integral,unsigned int,true);
    SET_TYPE_TRAIT(is_integral,unsigned long,true);
    SET_TYPE_TRAIT(is_integral,unsigned long long,true);
    SET_TYPE_TRAIT(is_integral,signed char,true);
    SET_TYPE_TRAIT(is_integral,signed short,true);
    SET_TYPE_TRAIT(is_integral,signed int,true);
    SET_TYPE_TRAIT(is_integral,signed long,true);
    SET_TYPE_TRAIT(is_integral,signed long long,true);
    SET_TYPE_TRAIT(is_integral,bool,true);
    SET_TYPE_TRAIT(is_integral,char,true);
    
    SET_TYPE_TRAIT(is_float,float,true);
    SET_TYPE_TRAIT(is_float,double,true);
    SET_TYPE_TRAIT(is_float,long double,true);
    
#undef SET_TYPE_TRAIT

    /// Simple binder for member functions with no arguments
    template <class T, typename resultT> 
    class BindNullaryMemFun {
    private:
        T* t;
        resultT (T::*op)();
    public:
        BindNullaryMemFun(T* t, resultT (T::*op)()) : t(t), op(op) {};
        resultT operator()() {return (t->*op)();};
    };

    /// Specialization of BindNullaryMemFun for void return
    template <class T> 
    class BindNullaryMemFun<T,void> {
    private:
        T* t;
        void (T::*op)();
    public:
        BindNullaryMemFun(T* t, void (T::*op)()) : t(t), op(op) {};
        void operator()() {(t->*op)();};
    };
    

    /// Simple binder for const member functions with no arguments
    template <class T, typename resultT> 
    class BindNullaryConstMemFun {
    private:
        const T* t;
        resultT (T::*op)() const;
    public:
        BindNullaryConstMemFun(const T* t, resultT (T::*op)() const) : t(t), op(op) {};
        resultT operator()() const {return (t->*op)();};
    };

    /// Specialization of BindNullaryConstMemFun for void return
    template <class T> 
    class BindNullaryConstMemFun<T,void> {
    private:
        const T* t;
        void (T::*op)() const;
    public:
        BindNullaryConstMemFun(const T* t, void (T::*op)() const) : t(t), op(op) {};
        void operator()() const {(t->*op)();};
    };
    
    /// Factory function for BindNullaryMemFun
    template <class T, typename resultT>
    inline BindNullaryMemFun<T,resultT> bind_nullary_mem_fun(T* t, resultT (T::*op)()) {
        return BindNullaryMemFun<T,resultT>(t,op);
    }

    /// Factory function for BindNullaryConstMemFun
    template <class T, typename resultT>
    inline BindNullaryConstMemFun<T,resultT> bind_nullary_mem_fun(const T* t, resultT (T::*op)() const) {
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
    struct ReturnWrapper<void> 
    {
        typedef Void type;
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

       \def ENABLE_IF(CONDITION,TYPEIFTRUE)
       \brief Macro to make enable_if<> template easier to use

       \def DISABLE_IF(CONDITION,TYPEIFTRUE)
       \brief Macro to make disable_if<> template easier to use

       \def IS_SAME(TYPEA, TYPEB)
       \brief Macro to make is_same<T> template easier to use

       \def RETURN_WRAPPERT(TYPE)
       \brief Returns TYPE unless TYPE is void when returns ReturnWrapper<void>

    */


#define REMREF(TYPE)    typename remove_reference< TYPE >::type
#define REMCONST(TYPE)  typename remove_const< TYPE >::type
#define RETURN_WRAPPERT(TYPE) typename ReturnWrapper< TYPE >::type

#define MEMFUN_RETURNT(MEMFUN) typename memfunc_traits< MEMFUN >::result_type
#define MEMFUN_OBJT(MEMFUN)    typename memfunc_traits< MEMFUN >::obj_type
#define MEMFUN_ARITY(MEMFUN)   memfunc_traits< MEMFUN >::arity
#define MEMFUN_ARG1T(MEMFUN)   typename memfunc_traits< MEMFUN >::arg1_type
#define MEMFUN_ARG2T(MEMFUN)   typename memfunc_traits< MEMFUN >::arg2_type
#define MEMFUN_ARG3T(MEMFUN)   typename memfunc_traits< MEMFUN >::arg3_type
#define MEMFUN_ARG4T(MEMFUN)   typename memfunc_traits< MEMFUN >::arg4_type

#define FUNCTION_RETURNT(FUNCTION) typename function_traits< FUNCTION >::result_type
#define FUNCTION_ARITY(FUNCTION)   function_traits< FUNCTION >::arity
#define FUNCTION_ARG1T(FUNCTION)   typename function_traits< FUNCTION >::arg1_type
#define FUNCTION_ARG2T(FUNCTION)   typename function_traits< FUNCTION >::arg2_type
#define FUNCTION_ARG3T(FUNCTION)   typename function_traits< FUNCTION >::arg3_type
#define FUNCTION_ARG4T(FUNCTION)   typename function_traits< FUNCTION >::arg4_type

#define DISABLE_IF(CONDITION,TYPEIFTRUE) typename disable_if< CONDITION, TYPEIFTRUE >::type
#define ENABLE_IF(CONDITION,TYPEIFTRUE)  typename  enable_if< CONDITION, TYPEIFTRUE >::type
#define IS_SAME(A, B) is_same< A, B >
    
} // end of namespace madness
#endif
