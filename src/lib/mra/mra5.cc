#include <mra/mra.h>
#ifdef FUNCTION_INSTANTIATE_5
#include <mra/mraimpl.h>

namespace madness {
    template class FunctionDefaults<5>;
    template class Function<double, 5>;
    template class Function<std::complex<double>, 5>;
    template class FunctionImpl<double, 5>;
    template class FunctionImpl<std::complex<double>, 5>;
    template class FunctionCommonData<double, 5>;
    template class FunctionCommonData<double_complex, 5>;
    template class Displacements<5>;

    template void plotdx<double,5>(const Function<double,5>&, const char*, const Tensor<double>&,
                                   const std::vector<long>&, bool binary);
    template void plotdx<double_complex,5>(const Function<double_complex,5>&, const char*, const Tensor<double>&,
                                   const std::vector<long>&, bool binary);
}
#endif

