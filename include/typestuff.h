#ifndef TYPESTUFF_H
#define TYPESTUFF_H

/// \file typestuff.h
/// \brief Grossly simplified Boost-like type traits and templates

/// \full Aim to be a compatible subset of Boost ... wish we could rely on
/// Boost being there but between bloat, bjam, and compiler
/// deficiences we have to be independent for the time being.

namespace madness {
    
    /// type_or_c<bool A, bool B>::value will be true if (A || B)
    template <bool A, bool B>
    class type_or_c {
    public: 
        static const bool value = true;
    };
    
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

    /// enable_if_c from Boost for conditionally instantiating templates based on type

    /// Evaluates to \c returnT if \c B is true, otherwise to an invalid type expression
    /// which causes the template expression in which it is used to not be considered for
    /// overload resolution.
    template <bool B, class returnT>
    struct enable_if_c {
        typedef returnT type;
    };

    template <class returnT>
    struct enable_if_c<false, returnT> {};

    /// enable_if from Boost for conditionally instantiating templates based on type

    /// Evaluates to \c returnT if \c Cond::value is true, otherwise to an invalid type expression
    /// which causes the template expression in which it is used to not be considered for
    /// overload resolution.
    template <class Cond, class returnT>
    struct enable_if : public enable_if_c<Cond::value, returnT> {};

    template <bool B, class returnT>
    struct disable_if_c {
        typedef returnT type;
    };

    template <class returnT>
    struct disable_if_c<true, returnT> {};

    template <class Cond, class returnT>
    struct disable_if : public disable_if_c<Cond::value, returnT> {};

    
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
    
} // end of namespace madness
#endif
