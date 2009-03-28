#include <mra/mra.h>
#ifdef FUNCTION_INSTANTIATE_6
#include <mra/mraimpl.h>

namespace madness {
    template class FunctionDefaults<6>;
    template class Function<double, 6>;
    template class Function<std::complex<double>, 6>;
    template class FunctionImpl<double, 6>;
    template class FunctionImpl<std::complex<double>, 6>;
    template class FunctionCommonData<double, 6>;
    template class FunctionCommonData<double_complex, 6>;
    template class Displacements<6>;

    template void plotdx<double,6>(const Function<double,6>&, const char*, const Tensor<double>&,
                                   const std::vector<long>&, bool binary);
    template void plotdx<double_complex,6>(const Function<double_complex,6>&, const char*, const Tensor<double>&,
                                           const std::vector<long>&, bool binary);
}
#endif

