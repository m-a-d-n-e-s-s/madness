// Copyright 2021 Adrian Hurtado


#include "x_space.h"

#include "response_functions.h"

namespace madness {
    auto to_response_matrix(const X_space& x) -> response_matrix {
        World& world = x.X[0][0].world();
        auto mX = response_matrix(x.num_states());
        int b = 0;
        for (auto& mi: mX) {
            mi = copy(world, x.X[b]);
            std::for_each(x.Y[b].begin(), x.Y[b].end(),
                          [&](const auto& yi) { mi.push_back(copy(yi)); });
            b++;
        }
        return mX;
    }

    /**
     * Flattens all response functions into a single vector of functions
     * @param x
     * @return
     */
    auto to_flattened_vector(const X_space& x) -> vector_real_function_3d {

        World& world = x.X[0][0].world();
        auto num_orbitals = 2 * x.num_orbitals();
        auto vij = vector_real_function_3d(x.num_states()*num_orbitals);
        auto mx = to_response_matrix(x);

        int b = 0;
        for (auto& mi: mx) {
            std::copy(mi.begin(), mi.end(), vij.begin() + b * num_orbitals);
            b++;
        }
        return vij;
    }

    auto to_X_space(const response_matrix& x) -> X_space {

        World& world = x[0][0].world();

        auto num_states = x.size();
        auto num_orbitals = size_t(x[0].size() / 2);
        auto x_space = X_space(world, num_states, num_orbitals);

        int b = 0;
        std::for_each(x.begin(), x.end(), [&](auto x_vec) {
            //auto norm_vi = norm2(world, x_vec);
            //if (world.rank() == 0) { print("norm in xvec i", norm_vi); }
            std::copy(x_vec.begin(), x_vec.begin() + num_orbitals, x_space.X[b].begin());
            std::copy(x_vec.begin() + num_orbitals, x_vec.end(), x_space.Y[b].begin());
            b++;
        });

        //  auto norms = x_space.norm2s();
        // if (world.rank() == 0) { print("norms after copy ", norms); }


        return x_space;
    }

    auto transposeResponseMatrix(const response_matrix& x) -> response_matrix {

        auto XT = response_matrix(x[0].size());

        auto b = 0;
        for (auto& xi: XT) {
            auto j = 0;
            xi = vector_real_function_3d(x.size());
            for (auto& xji: xi) { xji = copy(x[j++][b]); }
            b++;
        }
        return XT;
    }
    auto inner(const X_space& A, const X_space& B) -> Tensor<double> {
        MADNESS_ASSERT(size_states(A) > 0);
        MADNESS_ASSERT(size_orbitals(A) > 0);
        MADNESS_ASSERT(same_size(A, B));

        long size = static_cast<long>(A.n_states);

        auto a = to_response_matrix(A);
        auto b = to_response_matrix(B);

        World& world = a[0][0].world();

        auto a_transpose = transposeResponseMatrix(a);
        auto b_transpose = transposeResponseMatrix(b);
        /*
        std::for_each(a_transpose.begin(), a_transpose.end(),
                      [&](const auto& ai) { print(norm2(world, ai)); });
        std::for_each(b_transpose.begin(), b_transpose.end(),
                      [&](const auto& ai) { print(norm2(world, ai)); });
                      */
        // Container for result
        Tensor<double> result(a.size(), a.size());
        int p = 0;
        std::for_each(a_transpose.begin(), a_transpose.end(), [&](const auto& ati) {
            result += matrix_inner(world, ati, b_transpose[p++]);
        });

        /*
         * TODO Figure out why this won't work
    std::inner_product(
            a_transpose.begin(), a_transpose.end(), b_transpose.begin(), result,
            [](Tensor<double> a, const Tensor<double>& b) {
                print("a", a);
                print("b", b);
                auto ab = a + b;
                print("a + b = ", ab);
                return ab;
            },
            [&](const auto& ai, const auto& bi) {
                auto m = matrix_inner(world, ai, bi);
                print("m: ", m);
                return m;
            });
            */

        //print("results: ", result);
        return result;
    }
}// namespace madness
