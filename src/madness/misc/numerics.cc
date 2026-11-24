//
// Created by Jonathon Misiewicz on 5/16/23.
//

#include <madness/misc/misc.h>

#include <cmath>
#include <limits>

namespace madness {
    // Soft equality check
    bool nearlyEqual(double a, double b, double epsilon) {
        if (madness::is_nan(a) || madness::is_nan(b)) return false;
        if (!madness::is_finite(a) || !madness::is_finite(b)) {
            return !madness::is_finite(a) && !madness::is_finite(b) &&
                   (std::signbit(a) == std::signbit(b));
        }
        if (a == b) return true;

        double diff = std::abs(a - b);
        double norm = std::min((std::abs(a) + std::abs(b)), std::numeric_limits<double>::max());

        // Use relative epsilon for large numbers, absolute for small numbers
        return diff < std::max(epsilon, epsilon * norm);
    }
}
