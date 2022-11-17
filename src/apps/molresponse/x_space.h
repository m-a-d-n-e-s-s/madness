// Copyright 2021 Adrian Hurtado
#ifndef SRC_APPS_MOLRESPONSE_X_SPACE_H_
#define SRC_APPS_MOLRESPONSE_X_SPACE_H_

#include <madness/mra/mra.h>
#include <madness/tensor/tensor.h>

#include <algorithm>
#include <cstdint>
#include <memory>
#include <numeric>
#include <vector>

#include "molresponse/response_functions.h"

namespace madness {
    struct X_space;

    auto to_response_vector(const vector_real_function_3d &vec) -> vector_real_function_3d;
    auto create_response_matrix(const size_t &num_state, const size_t &num_orbitals)
            -> response_matrix;
    auto to_response_matrix(const X_space &x) -> response_matrix;
    auto to_conjugate_response_matrix(const X_space &x) -> response_matrix;
    auto to_flattened_vector(const X_space &x) -> vector_real_function_3d;
    auto to_X_space(const response_matrix &x) -> X_space;
    auto to_conjugate_X_space(const response_matrix &x) -> X_space;
    struct X_space {
    private:
        size_t n_states;  // Num. of resp. states
        size_t n_orbitals;// Num. of ground states

    public:
        response_space X, Y;

    public:
        size_t num_states() const { return n_states; }
        size_t num_orbitals() const { return n_orbitals; }
        // default constructor
        X_space() : n_states(0), n_orbitals(0), X(), Y() {}
        // Copy constructor
        X_space(const X_space &A)
            : n_states(size_states(A)), n_orbitals(size_orbitals(A)), X(A.X), Y(A.Y) {}
        X_space copy() const {
            auto &world = X[0][0].world();
            auto m = to_response_matrix(*this);
            auto copy_m = create_response_matrix(num_states(), num_orbitals());
            std::transform(m.begin(), m.end(), copy_m.begin(),
                           [&](const auto &mi) { return madness::copy(world, mi, true); });
            return to_X_space(copy_m);
        }
        /// Create a new copy of the function with different distribution and optional
        /// fence

        /// Works in either basis.  Different distributions imply
        /// asynchronous communication and the optional fence is
        /// collective.
        auto copy(const std::shared_ptr<WorldDCPmapInterface<Key<3>>> &pmap,
                  bool fence = false) const -> X_space {
            auto &world = X[0][0].world();
            auto m = to_response_matrix(*this);
            auto copy_m = create_response_matrix(num_states(), num_orbitals());
            std::transform(m.begin(), m.end(), copy_m.begin(),
                           [&](const auto &mi) { return madness::copy(world, mi, pmap, true); });
            return to_X_space(copy_m);
        }
        // assignment
        auto operator=(const X_space &B) -> X_space & {
            if (this != &B) {// is it the same object?
                this->n_states = B.num_states();
                this->n_orbitals = B.num_orbitals();

                this->X = B.X;
                this->Y = B.Y;
            }
            return *this;// NO SHALLOW COPIES
        }
        // Zero Constructor
        X_space(World &world, size_t n_states, size_t n_orbitals)
            : n_states(n_states), n_orbitals(n_orbitals), X(world, n_states, n_orbitals),
              Y(world, n_states, n_orbitals) {}
        // explicit constructor from 2 resonse_space
        explicit X_space(response_space &X, response_space &Y) {
            MADNESS_ASSERT(X.size() == Y.size());
            MADNESS_ASSERT(X[0].size() == Y[0].size());
            this->n_states = X.size();
            this->n_orbitals = X[0].size();
            this->X = X.copy();
            this->Y = Y.copy();
        }
        void clear() {
            X.clear();
            Y.clear();
        }
        auto operator+(const X_space &B) -> X_space {
            MADNESS_ASSERT(same_size(*this, B));
            World &world = this->X[0][0].world();

            auto ax = to_response_matrix(*this);
            world.gop.fence();
            auto bx = to_response_matrix(B);
            world.gop.fence();

            response_matrix add_x(num_states());

            std::transform(ax.begin(), ax.end(), bx.begin(), add_x.begin(),
                           [&](const auto &a, const auto &b) { return add(world, a, b, false); });
            world.gop.fence();

            return to_X_space(add_x);
        }

        auto operator+=(const X_space &B) -> X_space & {
            MADNESS_ASSERT(same_size(*this, B));
            this->X += B.X;
            this->Y += B.Y;
            return *this;
        }

        void push_back(vector_real_function_3d x, vector_real_function_3d y) {
            if (n_orbitals > 0) {
                MADNESS_ASSERT(n_orbitals == x.size());
                MADNESS_ASSERT(n_orbitals == y.size());
                MADNESS_ASSERT(x.size() == y.size());
            } else {// g_states == 0 (empty vector)
                n_orbitals = x.size();
            }
            MADNESS_ASSERT(x.size() == num_orbitals());
            MADNESS_ASSERT(y.size() == num_orbitals());

            n_states++;
            X.push_back(x);
            Y.push_back(y);
            // Be smart with g_states
        }
        void pop_back() {
            X.pop_back();
            Y.pop_back();
            n_states--;
            if (n_states == 0) { n_orbitals = 0; }
        }

        static X_space zero_functions(World &world, size_t n_states, size_t n_orbitals) {
            auto zeros = X_space(world, n_states, n_orbitals);

            for (auto &xi: zeros.X) {
                xi = ::madness::zero_functions<double, 3>(world, n_orbitals, false);
            }
            for (auto &yi: zeros.Y) {
                yi = ::madness::zero_functions<double, 3>(world, n_orbitals, false);
            }
            world.gop.fence();
            return zeros;
        }

        friend auto operator+(const X_space &A, const X_space &B) -> X_space {
            MADNESS_ASSERT(same_size(A, B));

            World &world = A.X[0][0].world();
            X_space result(world, A.n_states, A.n_orbitals);// create zero_functions

            auto ax = to_response_matrix(A);
            auto bx = to_response_matrix(B);

            response_matrix add_x(A.num_states());

            std::transform(ax.begin(), ax.end(), bx.begin(), add_x.begin(),
                           [&](auto a, auto b) { return add(world, a, b); });

            return to_X_space(add_x);

            result.X = A.X + B.X;
            result.Y = A.Y + B.Y;
            return result;
        }

        X_space operator-(const X_space B) {
            MADNESS_ASSERT(same_size(*this, B));
            World &world = this->X[0][0].world();
            X_space result(world, n_states, n_orbitals);
            result.X = X - B.X;
            result.Y = Y - B.Y;
            return result;
        }

        friend X_space operator-(const X_space &A, const X_space &B) {
            MADNESS_ASSERT(same_size(A, B));

            World &world = A.X[0][0].world();
            X_space result(world, A.n_states, A.n_orbitals);// create zero_functions

            result.X = A.X - B.X;
            result.Y = A.Y - B.Y;
            return result;
        }

        friend X_space operator*(const X_space &A, const double &b) {
            World &world = A.X[0][0].world();
            X_space result(world, A.n_states, A.n_orbitals);// create zero_functions

            result.X = A.X * b;
            result.Y = A.Y * b;
            return result;
        }
        friend X_space operator*(const double &b, const X_space &A) {
            World &world = A.X[0][0].world();
            X_space result(world, A.n_states, A.n_orbitals);// create zero_functions

            result.X = A.X * b;
            result.Y = A.Y * b;
            return result;
        }
        X_space operator*(const double &b) {
            this->X *= b;
            this->Y *= b;
            return *this;
        }

        friend X_space operator*(const X_space &A, const Function<double, 3> &f) {
            World &world = A.X[0][0].world();
            X_space result(world, A.n_states, A.n_orbitals);// create zero_functions

            result.X = A.X * f;
            result.Y = A.Y * f;
            return result;
        }
        friend auto operator*(const Function<double, 3> &f, const X_space &A) -> X_space {
            World &world = A.X[0][0].world();
            X_space result(world, A.n_states, A.n_orbitals);// create zero_functions

            result.X = A.X * f;
            result.Y = A.Y * f;
            return result;
        }

        friend auto operator*(const X_space &A, const Tensor<double> &b) -> X_space {
            MADNESS_ASSERT(size_states(A) > 0);
            MADNESS_ASSERT(size_orbitals(A) > 0);

            World &world = A.X[0][0].world();
            X_space result(world, A.n_states, A.n_orbitals);
            result.X = A.X * b;
            result.Y = A.Y * b;

            return result;
        }
        /***
         *
         * @param A
         * @param B
         * @return
         */
        friend auto inner(const X_space &A, const X_space &B) -> Tensor<double>;

        void truncate() {
            X.truncate_rf();
            Y.truncate_rf();
        }

        auto norm2s() -> Tensor<double> {
            World &world = X[0][0].world();

            Tensor<double> norms(num_states());
            for (size_t b = 0; b < num_states(); b++) {
                auto xb = madness::copy(world, X[b]);
                for (auto &yb: Y[b]) { xb.push_back(madness::copy(yb, true)); }
                norms[b] = sqrt(inner(xb, xb));
            }
            return norms;
        }

        auto component_norm2s() const -> Tensor<double> {
            World &world = X[0][0].world();

            auto rx = to_flattened_vector(*this);
            auto norms = norm2s_T(world, rx);
            return norms.reshape(n_states, 2 * n_orbitals);
        }

        friend auto size_states(const X_space &x) -> size_t { return x.n_states; }
        friend auto size_orbitals(const X_space &x) -> size_t { return x.n_orbitals; }
        friend auto same_size(const X_space &A, const X_space &B) -> bool {
            return ((size_states(A) == size_states(B) && size_orbitals(A) == size_orbitals(B)));
        }
    };

    // but the solver needs the functions initialized to zero for which we also need
    // the world object.

    struct X_vector : public X_space {
        X_vector(World &world, size_t n_orbtials) {
            this->X_space::zero_functions(world, size_t(1), n_orbtials);
        }

        X_vector(X_space A, size_t b) : X_space(A.X[0][0].world(), size_t(1), A.num_orbitals()) {
            X[0] = A.X[b];
            Y[0] = A.Y[b];
        }
        friend X_vector operator-(const X_vector &A, const X_vector &B) {
            MADNESS_ASSERT(same_size(A, B));

            World &world = A.X[0][0].world();
            X_vector result(world, size_orbitals(A));// create zero_functions
            result.X = A.X - B.X;
            result.Y = A.Y - B.Y;
            return result;
        }
        friend X_vector operator*(const X_vector &A, const double &c) {
            World &world = A.X[0][0].world();
            X_vector result(world, size_orbitals(A));// create zero_functions
            result.X = A.X * c;
            result.Y = A.Y * c;
            return result;
        }
        X_vector copy() const {
            X_vector copyX(X[0][0].world(), X.num_orbitals);
            copyX.X = X.copy();
            copyX.Y = Y.copy();
            return copyX;
        }
        auto operator+=(const X_vector &B) -> X_vector & {
            MADNESS_ASSERT(same_size(*this, B));


            this->X += B.X;
            this->Y += B.Y;

            return *this;
        }
        inline friend auto inner(X_vector &A, X_vector &B) -> double {
            MADNESS_ASSERT(size_states(A) == 1);
            MADNESS_ASSERT(size_orbitals(A) > 0);
            MADNESS_ASSERT(same_size(A, B));

            Tensor<double> G(1, 1);
            Tensor<double> G1(1, 1);
            Tensor<double> G2(1, 1);

            World &world = A.X[0][0].world();

            auto ax = madness::copy(world, A.X[0]);
            auto ay = madness::copy(world, A.Y[0]);

            auto bx = madness::copy(world, B.X[0]);
            auto by = madness::copy(world, B.Y[0]);

            for (auto &ayi: ay) { ax.push_back(madness::copy(ayi)); }
            for (auto &byi: by) { bx.push_back(madness::copy(byi)); };

            double result = inner(ax, bx);

            return result;
        }
    };
    // function object with allocator()()
    struct response_matrix_allocator {
        World &world;
        const size_t n_orbtials;
        response_matrix_allocator(World &world, size_t n_orbtials)
            : world(world), n_orbtials(n_orbtials) {}
        // overloading the default constructor () operator
        vector_real_function_3d operator()() {
            //print("allocator called with ", int(n_orbtials), " orbitals");
            // returning constructor of x_vector
            return zero_functions<double, 3>(world, n_orbtials);
        }
        // Copy constructor

        response_matrix_allocator operator=(const response_matrix_allocator &other) {
            return response_matrix_allocator(world, other.n_orbtials);
        }
    };
}// namespace madness

#endif// SRC_APPS_MOLRESPONSE_X_SPACE_H_
