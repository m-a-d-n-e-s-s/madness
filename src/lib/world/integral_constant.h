
#ifndef INTEGRAL_CONSTANT_H
#define INTEGRAL_CONSTANT_H
namespace madness {
    namespace tr1 {
        template <class T, T val>
        struct integral_constant
        {
            typedef integral_constant<T, val>  type;
            typedef T                          value_type;
            static const T value = val;
        };
        
        typedef integral_constant<bool, true>  true_type;
        typedef integral_constant<bool, false> false_type;
    }
}

namespace std {
    using madness::tr1::integral_constant;
    using madness::tr1::true_type;
    using madness::tr1::false_type;
}

#endif
