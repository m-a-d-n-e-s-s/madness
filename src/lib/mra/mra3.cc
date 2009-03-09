#include <mra/mraimpl.h>
#ifdef FUNCTION_INSTANTIATE_3

namespace madness {
    template class FunctionDefaults<3>;
    template class Function<double, 3>;
    template class Function<std::complex<double>, 3>;
    template class FunctionImpl<double, 3>;
    template class FunctionImpl<std::complex<double>, 3>;
    template class FunctionCommonData<double, 3>;
    template class FunctionCommonData<double_complex, 3>;
    template class Displacements<3>;

    template void plotdx<double,3>(const Function<double,3>&, const char*, const Tensor<double>&,
                                   const std::vector<long>&, bool binary);
    template void plotdx<double_complex,3>(const Function<double_complex,3>&, const char*, const Tensor<double>&,
                                   const std::vector<long>&, bool binary);
}

#endif

