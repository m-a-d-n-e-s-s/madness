//
// Created by Jonathon Misiewicz on 5/16/23.
//

#include <madness/misc/misc.h>

#include <limits>

namespace madness {
    // Soft equality check
    bool nearlyEqual(double a, double b, double epsilon) {
        // Handle equality (also covers infinity)
        if (a == b) return true;

        double diff = std::abs(a - b);
        double norm = std::min((std::abs(a) + std::abs(b)), std::numeric_limits<double>::max());

        // Use relative epsilon for large numbers, absolute for small numbers
        return diff < std::max(epsilon, epsilon * norm);
    }
}
