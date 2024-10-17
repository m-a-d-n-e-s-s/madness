// Copyright 2021 Adrian Hurtado

#include "x_space.h"

#include "response_functions.h"

namespace madness
{
    /**
     * @phi -> [phi conjugate(phi)]
     *
     * @param vec
     * @return vector_real_function_3d
     */
    auto to_response_vector(const vector_real_function_3d &vec)
        -> vector_real_function_3d
    {
        auto &world = vec[0].world();
        // copy the vector
        auto response_vector = vector_real_function_3d(2 * vec.size());
        auto n = vec.size();

        int i = 0;
        for (auto &r0i : response_vector)
        {
            r0i = vec[i++ % n];
        }
        auto copy_vec = copy(world, response_vector);
        world.gop.fence();
        return copy_vec;
    }

    /**
     * @brief Create a response matrix object
     *
     * @param num_states
     * @param num_orbitals
     * @return response_matrix
     */
    auto create_response_matrix(const size_t &num_states,
                                const size_t &num_orbitals) -> response_matrix
    {

        auto matrix = response_matrix(num_states);
        std::for_each(matrix.begin(), matrix.end(), [&](auto &xi)
                      { xi = vector_real_function_3d(num_orbitals); });
        return matrix;
    }

    /**
     * @ Converts Xspace object to response_matrix object
     *
     * @param x
     * @return response_matrix
     */
    auto to_response_matrix(const X_space &x) -> response_matrix
    {
        auto mX = response_matrix(x.num_states());
        auto num_orbitals = x.num_orbitals();
        int b = 0;
        std::for_each(mX.begin(), mX.end(), [&](vector_real_function_3d &mi)
                      {
            mi = vector_real_function_3d(2 * num_orbitals);
            std::copy(x.x[b].begin(), x.x[b].end(), mi.begin());// shallow copy
            std::copy(x.y[b].begin(), x.y[b].end(),
                      mi.begin() + num_orbitals);// shallow copy
            b++; });
        return mX;
    }

    auto to_conjugate_response_matrix(const X_space &x) -> response_matrix
    {
        World &world = x.x[0][0].world();
        auto mX = response_matrix(x.num_states());
        int b = 0;
        auto num_orbitals = x.num_orbitals();
        std::for_each(mX.begin(), mX.end(), [&](auto &mi)
                      {
            mi = vector_real_function_3d(2 * num_orbitals);
            std::copy(x.y[b].begin(), x.y[b].end(), mi.begin());// shallow copy
            std::copy(x.x[b].begin(), x.x[b].end(),
                      mi.begin() + num_orbitals);// shallow copy
            b++; });
        return mX;
    }

    /**
     * @brief Flattens all response functions into a single vector of functions
     * @param x
     * @return
     */
    auto to_flattened_vector(const X_space &x) -> vector_real_function_3d
    {

        World &world = x.x[0][0].world();
        auto num_orbitals = 2 * x.num_orbitals();
        auto vij = vector_real_function_3d(x.num_states() * num_orbitals);
        auto mx = to_response_matrix(x);
        int b = 0;
        for (const auto &mi : mx)
        {
            std::copy(mi.begin(), mi.end(), vij.begin() + b * num_orbitals);
            b++;
        }
        return vij;
    }

    auto to_X_space(const response_matrix &x) -> X_space
    {

        World &world = x[0][0].world();
        auto num_states = x.size();
        auto num_orbitals = size_t(x[0].size() / 2);
        auto x_space = X_space(world, num_states, num_orbitals);
        int b = 0;
        for (const auto &x_vec : x)
        {
            std::copy(x_vec.begin(), x_vec.begin() + num_orbitals,
                      x_space.x[b].begin());
            std::copy(x_vec.begin() + num_orbitals, x_vec.end(),
                      x_space.y[b].begin());
            b++;
        };
        return x_space;
    }

    /**
     * @brief response_matrix [x,y] -> Xspace X.x=y X.y=conjugate(x)
     *
     * @param x
     * @return X_space
     */
    auto to_conjugate_X_space(const response_matrix &x) -> X_space
    {

        World &world = x[0][0].world();

        auto num_states = x.size();
        auto num_orbitals = size_t(x[0].size() / 2);
        auto x_space = X_space(world, num_states, num_orbitals);

        int b = 0;
        std::for_each(x.begin(), x.end(), [&](auto x_vec)
                      {
            std::transform(x_vec.begin(), x_vec.begin() + num_orbitals,
                           x_space.y[b].begin(),
                           [&](const auto &xi) { return copy(xi, false); });
            std::transform(x_vec.begin() + num_orbitals, x_vec.end(),
                           x_space.x[b].begin(),
                           [&](const auto &xi) { return copy(xi, false); });
            b++; });
        world.gop.fence();

        return x_space;
    }

    auto transposeResponseMatrix(const response_matrix &x) -> response_matrix
    {

        auto XT = response_matrix(x[0].size());

        auto b = 0;
        for (auto &xi : XT)
        {
            auto j = 0;
            xi = vector_real_function_3d(x.size());
            for (auto &xji : xi)
            {
                xji = x[j++][b];
            }
            b++;
        }
        return XT;
    }

    /**
     * @brief Computes the matrix elements between two response spaces
     *
     * cij=inner(ai,bj)
     *
     * @param A
     * @param B
     * @return Tensor<double>
     */
    auto inner(const X_space &A, const X_space &B) -> Tensor<double>
    {
        MADNESS_ASSERT(size_states(A) > 0);
        MADNESS_ASSERT(size_orbitals(A) > 0);
        MADNESS_ASSERT(A.num_orbitals() == B.num_orbitals());

        // return response_space_inner(A.x, B.x) + response_space_inner(A.y, B.y);

        auto a = to_response_matrix(A);
        auto b = to_response_matrix(B);
        World &world = a[0][0].world();
        world.gop.fence();

        auto a_transpose = transposeResponseMatrix(a);
        auto b_transpose = transposeResponseMatrix(b);

        world.gop.fence();
        Tensor<double> result(A.num_states(), B.num_states());
        int p = 0;
        std::for_each(a_transpose.begin(), a_transpose.end(),
                      [&](const auto &ati)
                      {
                          result += matrix_inner(world, ati, b_transpose[p++]);
                          world.gop.fence();
                      });
        return result;
    }
    vector_real_function_3d copyToVector(const X_space &chi)
    {

        int n = static_cast<int>(chi.num_states());
        int m = static_cast<int>(chi.num_orbitals());

        vector_real_function_3d rf(2 * n * m);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                auto xindex = (2 * i * m) + j;
                auto yindex = (2 * i * m) + j + m;
                rf[xindex] = chi.x[i][j];
                rf[yindex] = chi.y[i][j];
            }
        }
        return rf;
    }

    void copyToXspace(const vector_real_function_3d &rf, X_space &chi)
    {

        int n = static_cast<int>(chi.num_states());
        int m = static_cast<int>(chi.num_orbitals());
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                auto xindex = (2 * i * m) + j;
                auto yindex = (2 * i * m) + j + m;
                chi.x[i][j] = rf[xindex];
                chi.y[i][j] = rf[yindex];
            }
        }
    }

    // In this implmentation we need to represent each x_space as a contigous block
    // of functions.
    vector_real_function_3d copyToVector(const response_space &chi)
    {

        int n = static_cast<int>(chi.num_states);
        int m = static_cast<int>(chi.num_orbitals);

        vector_real_function_3d rf(n * m);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                auto xindex = (i * m) + j;
                rf[xindex] = chi.x[i][j];
            }
        }
        return rf;
    }

    void copyToResponseSpace(const vector_real_function_3d &rf,
                             response_space &chi)
    {

        int n = static_cast<int>(chi.num_states);
        int m = static_cast<int>(chi.num_orbitals);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                auto xindex = (i * m) + j;
                chi.x[i][j] = rf[xindex];
            }
        }
    }

} // namespace madness
