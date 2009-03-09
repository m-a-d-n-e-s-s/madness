#include <mra/mraimpl.h>

#ifdef FUNCTION_INSTANTIATE_2
namespace madness {
    template class FunctionDefaults<2>;
    template class Function<double, 2>;
    template class Function<std::complex<double>, 2>;
    template class FunctionImpl<double, 2>;
    template class FunctionImpl<std::complex<double>, 2>;
    template class FunctionCommonData<double, 2>;
    template class FunctionCommonData<double_complex, 2>;
    template class Displacements<2>;

    template void plotdx<double,2>(const Function<double,2>&, const char*, const Tensor<double>&,
                                   const std::vector<long>&, bool binary);
    template void plotdx<double_complex,2>(const Function<double_complex,2>&, const char*, const Tensor<double>&,
                                   const std::vector<long>&, bool binary);
}
#endif
