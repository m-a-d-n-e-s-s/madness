//
// Created by Jonathon Misiewicz on 2/11/25.
//

#ifndef MPQC_SAP_INTERPOLATORS_H
#define MPQC_SAP_INTERPOLATORS_H

#include <madness/misc/interpolation_1d.h>

namespace madness {
    extern std::vector<CubicInterpolationTable<double>> SAPCharges;
}

#endif  // MPQC_SAP_INTERPOLATORS_H
