#ifndef MADNESS_TYPE_TRAITS_BITS_H
#define MADNESS_TYPE_TRAITS_BITS_H

// Minimal functionality

#include <cstddef>
#include <stdint.h>
#include <world/enable_if.h>
#include <world/integral_constant.h>

#include <world/function_traits_bits.h>

namespace madness {
    namespace tr1 {

        // for simplicity in the same order as
        // http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2003/n1424.htm
        /// true if T is void

        template <typename T>
        struct is_void : public false_type {} ;

        /// true if T is void
        template <>
        struct is_void <void> : public true_type {};

        /// is_integral<T>::value will be true if T is a fundamental integer type
        template <class T>
        class is_integral {
        public:
            static const bool value = false;
        };

        /// is_floating_point<T>::value will be true if T is a fundamental floating-point type
        template <class T>
        struct is_floating_point : public false_type {};

    #define SET_TYPE_TRAIT(trait, T) template<> struct trait < T > : public true_type {}

        SET_TYPE_TRAIT(is_integral,unsigned char);
        SET_TYPE_TRAIT(is_integral,unsigned short);
        SET_TYPE_TRAIT(is_integral,unsigned int);
        SET_TYPE_TRAIT(is_integral,unsigned long);
        SET_TYPE_TRAIT(is_integral,unsigned long long);
        SET_TYPE_TRAIT(is_integral,signed char);
        SET_TYPE_TRAIT(is_integral,signed short);
        SET_TYPE_TRAIT(is_integral,signed int);
        SET_TYPE_TRAIT(is_integral,signed long);
        SET_TYPE_TRAIT(is_integral,signed long long);
        SET_TYPE_TRAIT(is_integral,bool);
        SET_TYPE_TRAIT(is_integral,char);

        SET_TYPE_TRAIT(is_floating_point,float);
        SET_TYPE_TRAIT(is_floating_point,double);
        SET_TYPE_TRAIT(is_floating_point,long double);

    #undef SET_TYPE_TRAIT

        /// is_array<T>::value will be true if T is an array type
        template <typename T>
        struct is_array : public false_type {};

        /// is_array<T>::value will be true if T is an array type
        template <typename T, std::size_t n>
        struct is_array<T[n]> : public true_type {};

        /// is_pointer<T>::value will be true if T is a pointer type
        template <typename T>
        struct is_pointer : public false_type {};

        /// is_pointer<T>::value will be true if T is a pointer type
        template <typename T>
        struct is_pointer<T*> : public true_type {};

        /// is_reference<T>::value will be true if T is a reference type
        template <typename T>
        struct is_reference : public false_type {};

        /// is_reference<T>::value will be true if T is a reference type
        template <typename T>
        struct is_reference<T&> : public true_type {};
        
        /// is_function<T>::value is true if T is a function pointer
        template <typename T>
        struct is_function {
            static const bool value = detail::function_traits<T>::value;
        };

        /// is_member_function_pointer<T>::value is true if T is a member function pointer
        template <typename T>
        struct is_member_function_pointer {
            static const bool value = detail::memfunc_traits<T>::value;
        };

        // Not using these
        // is_member_object_pointer
        // is_enum
        // is_union
        // is_class

        template <class T> struct is_arithmetic { 
            static const bool value = is_integral<T>::value || is_floating_point<T>::value; 
            typedef bool value_type;
            typedef integral_constant<value_type,value> type;
            //operator type()const;
        };

        /// is_fundamental<T>::value will be true if T is a fundamental type
        template <class T>
        struct is_fundamental {
            typedef bool value_type;
            static const bool value = is_integral<T>::value  || is_floating_point<T>::value || is_void<T>::value; 
            //operator type()const;
        };

        template <class T> struct is_object{
            static const bool value = 
                !(is_function<T>::value 
                  || is_reference<T>::value 
                  || is_void<T>::value); 
            typedef bool value_type;
            typedef integral_constant<value_type,value> type;
            //operator type()const;
        };

        // is_scalar

        template <class T> struct is_compound{
            static const bool value = !is_fundamental<T>::value; 
            typedef bool value_type;
            typedef integral_constant<value_type,value> type;
            //operator type()const;
        };

        // is_member_pointer

        /// is_const<T>::value will be true if T is const
        template <typename T>
        struct is_const : public false_type {} ;

        /// is_const<T>::value will be true if T is const
        template <typename T>
        struct is_const<const T> : public true_type {};

        // is_volatile
        // is_pod
        // is_empty
        // is_polymorphic
        // is_abstract;
        // has_trivial_constructor;
        // has_trivial_copy;
        // has_trivial_assign;
        // has_trivial_destructor;
        // has_nothrow_constructor;
        // has_nothrow_copy;
        // has_nothrow_assign;
        // is_signed;
        // is_unsigned;

        /// is_same<A,B> returns true if A and B are the same type
        template <typename A, typename B> 
        struct is_same : public false_type {};

        /// is_same<A,B> returns true if A and B are the same type
        template <typename A> 
        struct is_same<A,A> : public true_type {};

        // is_convertible

        namespace detail {
            typedef char (&yes)[1];
            typedef char (&no)[2];

            template <typename B, typename D>
            struct Host
            {
                operator B*() const;
                operator D*();
            };
        }

        template <typename B, typename D>
        struct is_base_of
        {
            template <typename T>  
            static detail::yes check(D*, T);
            static detail::no check(B*, int);

            static const bool value = sizeof(check(detail::Host<B,D>(), int())) == sizeof(detail::yes);
        };

        // is_convertible
        // is_explicitly_convertible

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

        /// remove_volatile<volatile T>::type will be T
        template <typename T>
        struct remove_volatile {
            typedef T type;
        };

        /// remove_volatile<volatile T>::type will be T
        template <typename T>
        struct remove_volatile<volatile T> {
            typedef T type;
        };

        /// remove_cv<const T>::type or remove_cv<volatile T>::type will be T
        template <typename T>
        struct remove_cv {
            typedef typename remove_volatile<typename remove_const<T>::type>::type type;
        };

        /// add_const<T>::type will be const T
        template <typename T>
        struct add_const {
            typedef const T type;
        };

        /// add_const<T>::type will be const T
        template <typename T>
        struct add_const<const T> {
            typedef const T type;
        };


        // template <class T> struct add_volatile;
        // template <class T> struct add_cv;

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

        // template <class T> struct add_lvalue_reference;
        // template <class T> struct add_rvalue_reference;

        // // 20.6.6.3, sign modifications:
        // template <class T> struct make_signed;
        // template <class T> struct make_unsigned;

        // // 20.6.6.4, array modifications:
        // template <class T> struct remove_extent;
        // template <class T> struct remove_all_extents;
        // // 20.6.6.5, pointer modifications:

        /// remove_pointer<T*>::type will be T
        template <typename T>
        struct remove_pointer {
            typedef T type;
        };

        /// remove_pointer<T>::type will be T
        template <typename T>
        struct remove_pointer<T*> {
            typedef T type;
        };

        // template <class T> struct add_pointer;

        // // 20.6.7, other transformations:
        // template <std::size_t Len, std::size_t Align> struct aligned_storage;
        // template <std::size_t Len, class... Types> struct aligned_union;
        // template <class T> struct decay;

        /// std::enable_if is similar to boost::enable_if_c for conditionally instantiating templates based on type

        /// Evaluates to \c returnT if \c B is true, otherwise to an invalid type expression
        /// which causes the template expression in which it is used to not be considered for
        /// overload resolution.
        template <bool B, class returnT>
        struct enable_if {
            typedef returnT type;
        };

        /// std::enable_if is similar to boost::enable_if_c for conditionally instantiating templates based on type
        template <class returnT>
        struct enable_if<false, returnT> {};

        // template <bool, class T, class F> struct conditional;
        // template <class... T> struct common_type;
    }
}

#endif
