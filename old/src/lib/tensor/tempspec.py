
typelist = ["int","long","double","float","double_complex","float_complex"]

complex_typelist = ["double_complex","float_complex"]

gemm_typelist = ["double","float","double_complex","float_complex"]

f = open("tensor_spec.h","w")
for t in typelist:
    f.write("\n// Instantiations for %s\n" % t)
    f.write("template class Tensor<%s>;\n" % t)
    f.write("template class SliceTensor<%s>;\n" % t)
    f.write("template std::ostream& operator << (std::ostream& s, const Tensor<%s>& t);\n" % t)
    f.write("template Tensor<%s> copy(const Tensor<%s>& t);\n" % (t,t))
    f.write("template Tensor<%s> outer(const Tensor<%s>& left, const Tensor<%s>& right);\n" % (t,t,t))
    f.write("template void inner_result(const Tensor<%s>& left, const Tensor<%s>& right,\n" % (t,t))
    f.write("                           long k0, long k1, Tensor<%s>& result);\n" % t)
    f.write("template Tensor<%s> inner(const Tensor<%s>& left, const Tensor<%s>& right,\n" % (t,t,t))
    f.write("                         long k0, long k1);\n")
    f.write("template Tensor<%s> transform(const Tensor<%s>& t, const Tensor<%s>& c);\n" % (t,t,t))
    f.write("template void fast_transform(const Tensor<%s>& t, const Tensor<%s>& c, " \
        "Tensor<%s>& result, Tensor<%s>& workspace);\n"% (t,t,t,t))
    f.write("template Tensor< Tensor<%s>::scalar_type > abs(const Tensor<%s>& t);\n" % (t,t))
    f.write("template Tensor<%s> transpose(const Tensor<%s>& t);\n" % (t,t))


f.write("\n// Instantiations only for complex types\n")
for t in complex_typelist:
    f.write("\n// Instantiations for %s" % t)
    f.write("template Tensor< Tensor<%s>::scalar_type > arg(const Tensor<%s>& t);\n" % (t,t))
    f.write("template Tensor< Tensor<%s>::scalar_type > real(const Tensor<%s>& t);\n" % (t,t))
    f.write("template Tensor< Tensor<%s>::scalar_type > imag(const Tensor<%s>& t);\n" % (t,t))
    f.write("template Tensor<%s> conj(const Tensor<%s>& t);\n" % (t,t))
    f.write("template Tensor<%s> conj_transpose(const Tensor<%s>& t);\n" % (t,t))

f.close()
    
## f.write("\n\nPUT THESE AT THE BOTTOM OF mxm.cc\n"
## for t in typelist:
##     f.write("template void mTxm<%s>(long dimi, long dimj, long dimk, %s * c, const %s *a, const %s *b);" % (t,t,t,t)
##     f.write("template void mxmT<%s>(long dimi, long dimj, long dimk, %s * c, const %s *a, const %s *b);" % (t,t,t,t)
##     f.write("template void mTxmT<%s>(long dimi, long dimj, long dimk, %s * c, const %s *a, const %s *b);" % (t,t,t,t)
##     f.write("template void mxm<%s>(long dimi, long dimj, long dimk, %s * c, const %s *a, const %s *b);" % (t,t,t,t)
    

##f.write("\n\nPUT THESE AT THE BOTTOM OF tensoriter.cc\n"
f = open("tensoriter_spec.h","w")

a = []
for t in typelist:
    for q in typelist:
        for r in typelist:
            a.append("template class TensorIterator<%s,%s,%s>;\n" % (t,q,r))
a.sort()

prev = a[0]
f.write(prev)
for cur in a[1:]:
    if cur != prev:
        f.write(cur)
        prev = cur
        
f.close()
