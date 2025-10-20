#pragma once
#include <madness/mra/vmra.h>
#include <madness/world/world.h>
#include "ResponseState.hpp"
#include <variant> // C++ 17

using namespace madness;

struct StaticRestrictedResponse;
struct DynamicRestrictedResponse;
struct StaticUnrestrictedResponse;
struct DynamicUnrestrictedResponse;


using ResponseVector =
std::variant<StaticRestrictedResponse, DynamicRestrictedResponse,
    StaticUnrestrictedResponse, DynamicUnrestrictedResponse>;
/**
 * @brief Holds all response functions describing
 * a response state, one, two, or four vectors of real functions.
 *
 * Static/Restricted shell:    x_alpha only
 * Dynamic/Restricted shell:   x_alpha and y_alpha
 * Static/Unrestricted shell:  x_alpha and x_beta
 * Dynamic/Unrestricted shell: x_alpha, y_alpha, x_beta and y_beta
 *
 * Each type has an allocator depending on the number of orbitals for the system.
 *
 * If you want to access all the response functions, use the `flat` vector.
 *
 * If you want to access individual vectors, use the `x_alpha`, `y_alpha`,
 *
 * std::get<StaticRestrictedResponse>(response).x_alpha
 * std::get<DynamicRestrictedResponse>(response).y_alpha
 *
 *
 *
 *
 *
 *
 *
 */

// Definition of StaticRestrictedResponse.
struct StaticRestrictedResponse {
    vector_real_function_3d x_alpha;
    vector_real_function_3d flat;

    StaticRestrictedResponse() = default;

    explicit StaticRestrictedResponse(const size_t &num_orbitals)
        : x_alpha(num_orbitals) {
        flatten();
    }

    void sync() {
        for (size_t i = 0; i < x_alpha.size(); ++i) {
            x_alpha[i] = flat[i];
        }
    }

    void flatten() { flat = x_alpha; }
};

struct DynamicRestrictedResponse {
    vector_real_function_3d x_alpha;
    vector_real_function_3d y_alpha;
    vector_real_function_3d flat;

    DynamicRestrictedResponse() = default;

    explicit DynamicRestrictedResponse(const size_t &num_orbitals)
        : x_alpha(num_orbitals), y_alpha(num_orbitals) {
        flatten();
    }

    // x_alpha + y_alpha in a single vector
    void sync() {
        for (size_t i = 0; i < x_alpha.size(); ++i) {
            x_alpha[i] = flat[i];
            y_alpha[i] = flat[i + x_alpha.size()];
        }
    }

    void flatten() {
        flat = x_alpha;
        flat.insert(flat.end(), y_alpha.begin(), y_alpha.end());
    }
};

struct StaticUnrestrictedResponse {
    vector_real_function_3d x_alpha;
    vector_real_function_3d x_beta;
    vector_real_function_3d flat;

    StaticUnrestrictedResponse() = default;

    explicit StaticUnrestrictedResponse(const size_t &num_orbitals)
        : x_alpha(num_orbitals), x_beta(num_orbitals) {
        flatten();
    }

    void sync() {
        for (size_t i = 0; i < x_alpha.size(); ++i) {
            x_alpha[i] = flat[i];
            x_beta[i] = flat[i + x_alpha.size()];
        }
    }

    void flatten() {
        flat = x_alpha;
        flat.insert(flat.end(), x_beta.begin(), x_beta.end());
    }

    // x_alpha + x_beta in a single vector
};

struct DynamicUnrestrictedResponse {
    vector_real_function_3d x_alpha;
    vector_real_function_3d y_alpha;
    vector_real_function_3d x_beta;
    vector_real_function_3d y_beta;
    vector_real_function_3d flat;

    DynamicUnrestrictedResponse() = default;

    explicit DynamicUnrestrictedResponse(const size_t &num_orbitals)
        : x_alpha(num_orbitals), y_alpha(num_orbitals), x_beta(num_orbitals),
          y_beta(num_orbitals) {
        flatten();
        // flat.insert(flat.end(), y_beta.begin(), y_beta.end());
    }

    void flatten() {
        flat = x_alpha;
        flat.insert(flat.end(), y_alpha.begin(), y_alpha.end());
        flat.insert(flat.end(), x_beta.begin(), x_beta.end());
        flat.insert(flat.end(), y_beta.begin(), y_beta.end());
    }

    void sync() {
        for (size_t i = 0; i < x_alpha.size(); ++i) {
            x_alpha[i] = flat[i];
            y_alpha[i] = flat[i + x_alpha.size()];
            x_beta[i] = flat[ (2 * x_alpha.size()) + i];
            y_beta[i] = flat[i + (3 * x_alpha.size())];
        }
    }

    // x_alpha + y_alpha + x_beta + y_beta in a single vector
};


/**
 * @brief Create a ResponseVector based on the number of orbitals and
 * response type.
 *
 * @param num_orbitals The number of orbitals.
 * @param is_static Whether the response is static.
 * @param is_unrestricted Whether the response is unrestricted.
 * @return ResponseVector The created ResponseVector.
 */
inline ResponseVector make_response_vector(size_t num_orbitals, bool is_static,
                                           bool is_unrestricted) {
    if (!is_unrestricted && is_static) {
        return StaticRestrictedResponse(num_orbitals);
    } else if (!is_unrestricted && !is_static) {
        return DynamicRestrictedResponse(num_orbitals);
    } else if (is_unrestricted && is_static) {
        return StaticUnrestrictedResponse(num_orbitals);
    } else {
        // unrestricted && dynamic
        return DynamicUnrestrictedResponse(num_orbitals);
    }
}


/**
 * @brief Get the flat vector from a ResponseVector.
 *
 * @param vec The ResponseVector to extract the flat vector from.
 * @return const vector_real_function_3d& The flat vector.
 */
inline const vector_real_function_3d &get_flat(const ResponseVector &vec) {
    return std::visit(
        [](const auto &v) -> const vector_real_function_3d & { return v.flat; },
        vec);
}
