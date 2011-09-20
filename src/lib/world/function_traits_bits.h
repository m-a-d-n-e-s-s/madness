#ifndef MADNESS_FUNCTION_TRAITS
#define MADNESS_FUNCTION_TRAITS
// Used only if we don't have boost or C++0x

namespace madness {
    namespace detail {
        /// Function traits in the spirt of boost function traits
        template <typename functionT>
        struct function_traits : public std::false_type {};
        
        /// Function traits in the spirt of boost function traits
        template <typename returnT>
        struct function_traits<returnT(*)()> {
            static const bool value = true;
            static const int arity = 0;
            typedef returnT result_type;
        };
        
        /// Function traits in the spirt of boost function traits
        template <typename returnT, typename arg1T>
        struct function_traits<returnT(*)(arg1T)> {
            static const bool value = true;
            static const int arity = 1;
            typedef returnT result_type;
            typedef arg1T arg1_type;
        };
        
        /// Function traits in the spirt of boost function traits
        template <typename returnT, typename arg1T, typename arg2T>
        struct function_traits<returnT(*)(arg1T,arg2T)> {
            static const bool value = true;
            static const int arity = 2;
            typedef returnT result_type;
            typedef arg1T arg1_type;
            typedef arg2T arg2_type;
        };
        
        /// Function traits in the spirt of boost function traits
        template <typename returnT, typename arg1T, typename arg2T, typename arg3T>
        struct function_traits<returnT(*)(arg1T,arg2T,arg3T)> {
            static const bool value = true;
            static const int arity = 3;
            typedef returnT result_type;
            typedef arg1T arg1_type;
            typedef arg2T arg2_type;
            typedef arg3T arg3_type;
        };
        
        /// Function traits in the spirt of boost function traits
        template <typename returnT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        struct function_traits<returnT(*)(arg1T,arg2T,arg3T,arg4T)> {
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
        struct function_traits<returnT(*)(arg1T,arg2T,arg3T,arg4T,arg5T)> {
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
        struct function_traits<returnT(*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T)> {
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
        struct function_traits<returnT(*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T)> {
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
        struct function_traits<returnT(*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T)> {
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
        struct function_traits<returnT(*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T,arg9T)> {
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
        struct memfunc_traits<returnT(objT::*)()> {
            static const bool value = true;
            static const int arity = 0;
            static const bool constness = false;
            typedef objT obj_type;
            typedef returnT result_type;
        };
        
        /// Member function traits in the spirt of boost function traits
        template <typename objT, typename returnT, typename arg1T>
        struct memfunc_traits<returnT(objT::*)(arg1T)> {
            static const bool value = true;
            static const int arity = 1;
            static const bool constness = false;
            typedef objT obj_type;
            typedef returnT result_type;
            typedef arg1T arg1_type;
        };
        
        /// Member function traits in the spirt of boost function traits
        template <typename objT, typename returnT, typename arg1T, typename arg2T>
        struct memfunc_traits<returnT(objT::*)(arg1T,arg2T)> {
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
        struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T)> {
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
        struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T)> {
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
        struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T)> {
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
        struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T)> {
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
        struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T)> {
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
        struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T)> {
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
        struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T,arg9T)> {
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
        struct memfunc_traits<returnT(objT::*)() const> {
            static const bool value = true;
            static const int arity = 0;
            static const bool constness = true;
            typedef objT obj_type;
            typedef returnT result_type;
        };
        
        /// Member function traits in the spirt of boost function traits
        template <typename objT, typename returnT, typename arg1T>
        struct memfunc_traits<returnT(objT::*)(arg1T) const> {
            static const bool value = true;
            static const int arity = 1;
            static const bool constness = true;
            typedef objT obj_type;
            typedef returnT result_type;
            typedef arg1T arg1_type;
        };
        
        /// Member function traits in the spirt of boost function traits
        template <typename objT, typename returnT, typename arg1T, typename arg2T>
        struct memfunc_traits<returnT(objT::*)(arg1T,arg2T) const> {
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
        struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T) const> {
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
        struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T) const> {
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
        struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T) const> {
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
        struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T) const> {
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
        struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T) const> {
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
        struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T) const> {
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
        struct memfunc_traits<returnT(objT::*)(arg1T,arg2T,arg3T,arg4T,arg5T,arg6T,arg7T,arg8T,arg9T) const> {
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
               
        template <typename fnT, typename Enabler = void>
        struct result_of {
            typedef typename fnT::result_type type;
        };
        
        template <typename fnT>
        struct result_of<fnT, typename enable_if_c<function_traits<fnT>::value>::type> {
            typedef typename function_traits<fnT>::result_type type;
        };
        
        template <typename fnT>
        struct result_of<fnT, typename enable_if_c<memfunc_traits<fnT>::value>::type> {
            typedef typename memfunc_traits<fnT>::result_type type;
        };
    }
}
#endif
