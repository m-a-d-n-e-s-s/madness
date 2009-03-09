#include <mra/mraimpl.h>

namespace madness {

    template <>
    ConcurrentHashMap< double, SharedPtr< GaussianConvolution1D<double> > >  GaussianConvolution1DCache<double>::map = ConcurrentHashMap< double, SharedPtr< GaussianConvolution1D<double> > >();

    template <>
    ConcurrentHashMap< double, SharedPtr< GaussianConvolution1D<double_complex> > > GaussianConvolution1DCache<double_complex>::map = ConcurrentHashMap< double, SharedPtr< GaussianConvolution1D<double_complex> > >();

#ifdef FUNCTION_INSTANTIATE_1
    template class FunctionDefaults<1>;
    template class Function<double, 1>;
    template class Function<std::complex<double>, 1>;
    template class FunctionImpl<double, 1>;
    template class FunctionImpl<std::complex<double>, 1>;
    template class FunctionCommonData<double, 1>;
    template class FunctionCommonData<double_complex, 1>;
    template class Displacements<1>;

    template void plotdx<double,1>(const Function<double,1>&, const char*, const Tensor<double>&,
                                   const std::vector<long>&, bool binary);
    template void plotdx<double_complex,1>(const Function<double_complex,1>&, const char*, const Tensor<double>&,
                                           const std::vector<long>&, bool binary);
#endif

}

/// Quietly used as a global lock when looking for bugs with multiple threads
madness::Mutex THELOCK;

