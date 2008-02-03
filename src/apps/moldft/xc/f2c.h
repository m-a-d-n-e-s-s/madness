#ifndef MADNESS_F2C_H
#define MADNESS_F2C_H

#include <math.h>
#include <madness_config.h>

typedef double doublereal;
typedef MADNESS_FORINT integer;
typedef int logical;

#define max(a,b) ((a)>=(b)?(a):(b))

static double pow_dd(doublereal* a, doublereal* b) {
    return pow(*a,*b);
}

#endif
