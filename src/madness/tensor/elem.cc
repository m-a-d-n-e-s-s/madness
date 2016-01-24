#include <madness/madness_config.h>

#ifdef MADNESS_HAS_ELEMENTAL

#include <madness/tensor/elem.h>

namespace madness {
    template
    void gesvp(World& world, const Tensor<double>& a, const Tensor<double>& b, Tensor<double>& x);

    template
    void sygvp(World& world, const Tensor<double>& A, const Tensor<double>& B, int itype,
              Tensor<double>& V, Tensor<double>& e);
}

#else 

/////int this_is_not_used_junk_junk_junk;

#endif //MADNESS_HAS_ELEMENTAL
