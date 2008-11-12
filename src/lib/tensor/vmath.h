#ifndef MAD_VMATH_H
#define MAD_VMATH_H

#include <madness_config.h>

#ifdef HAVE_MKL
#include <mkl.h>

#elif defined(HAVE_ACML)
void vzExp(int n, const double_complex* x, double_complex* y);
#endif

#endif

