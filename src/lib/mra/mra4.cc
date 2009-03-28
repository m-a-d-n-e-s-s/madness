
#include <mra/mraimpl.h>
#ifdef FUNCTION_INSTANTIATE_4

namespace madness {
    template class FunctionDefaults<4>;
    template class Function<double, 4>;
    template class Function<std::complex<double>, 4>;
    template class FunctionImpl<double, 4>;
    template class FunctionImpl<std::complex<double>, 4>;
    template class FunctionCommonData<double, 4>;
    template class FunctionCommonData<double_complex, 4>;
    template class Displacements<4>;

    template void plotdx<double,4>(const Function<double,4>&, const char*, const Tensor<double>&,
                                   const std::vector<long>&, bool binary);
    template void plotdx<double_complex,4>(const Function<double_complex,4>&, const char*, const Tensor<double>&,
                                           const std::vector<long>&, bool binary);
}
#endif

