#ifndef MOLRESPONSE_V3_RESPONSEFUNCTIONS_HPP
#define MOLRESPONSE_V3_RESPONSEFUNCTIONS_HPP

#include <madness/mra/mra.h>
#include <madness/world/world.h>
#include <string>

namespace molresponse_v3 {

using namespace madness;

/// Response function storage for one perturbation point.
///
/// Templated on the value type T (double for electric response,
/// std::complex<double> for magnetic response) and dimension NDIM.
/// The base storage is std::vector<Function<T,NDIM>> — the same
/// pattern SCF uses for amo/bmo. All vmra.h operations work directly
/// on the members since they are templated on T.
///
/// The three response types map to storage as:
///   Static (ω=0): x only, y empty
///   Full (ω≠0):   x and y populated
///   TDA:          x only, y empty
///
/// No arithmetic operators — the solver works directly on the
/// vector<Function<T,NDIM>> members using vmra routines.
template <typename T = double, int NDIM = 3>
struct ResponseState {

    using functionT = Function<T, NDIM>;
    using vecfuncT = std::vector<functionT>;

    // ---- Storage ----

    vecfuncT x_alpha;   // always present
    vecfuncT y_alpha;   // empty if static/TDA
    vecfuncT x_beta;    // empty if spin-restricted
    vecfuncT y_beta;    // empty if spin-restricted or static/TDA

    // ---- Queries ----

    bool is_restricted() const { return x_beta.empty(); }
    bool has_y() const { return !y_alpha.empty(); }

    long num_alpha() const { return static_cast<long>(x_alpha.size()); }
    long num_beta() const { return static_cast<long>(x_beta.size()); }
    long num_orbitals() const { return num_alpha(); }

    long total_size() const {
        return static_cast<long>(
            x_alpha.size() + y_alpha.size() +
            x_beta.size() + y_beta.size());
    }

    // ---- Flat representation ----

    /// Concatenate all vectors: [x_alpha | y_alpha | x_beta | y_beta]
    /// Function<T,NDIM> has shared_ptr semantics — in-place vmra
    /// operations on the returned vector mutate the originals.
    vecfuncT flat() const {
        vecfuncT f;
        f.reserve(total_size());
        f.insert(f.end(), x_alpha.begin(), x_alpha.end());
        f.insert(f.end(), y_alpha.begin(), y_alpha.end());
        f.insert(f.end(), x_beta.begin(), x_beta.end());
        f.insert(f.end(), y_beta.begin(), y_beta.end());
        return f;
    }

    /// Split a flat vector back into the component vectors.
    void from_flat(const vecfuncT& f) {
        MADNESS_CHECK(static_cast<long>(f.size()) == total_size());

        auto it = f.begin();

        x_alpha.assign(it, it + num_alpha());
        it += num_alpha();

        if (has_y()) {
            y_alpha.assign(it, it + num_alpha());
            it += num_alpha();
        }

        if (!is_restricted()) {
            x_beta.assign(it, it + num_beta());
            it += num_beta();

            if (has_y()) {
                y_beta.assign(it, it + num_beta());
                it += num_beta();
            }
        }
    }

    // ---- Archive I/O ----

    void save(World& world, const std::string& filename) const {
        archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive>
            ar(world, filename.c_str());

        int k = FunctionDefaults<NDIM>::get_k();
        bool _has_y = has_y();
        bool _has_beta = !is_restricted();
        long na = num_alpha();

        ar & k & _has_y & _has_beta & na;

        for (const auto& f : x_alpha) ar & f;
        if (_has_y) {
            for (const auto& f : y_alpha) ar & f;
        }

        if (_has_beta) {
            long nb = num_beta();
            ar & nb;
            for (const auto& f : x_beta) ar & f;
            if (_has_y) {
                for (const auto& f : y_beta) ar & f;
            }
        }
    }

    static ResponseState load(World& world, const std::string& filename) {
        archive::ParallelInputArchive<archive::BinaryFstreamInputArchive>
            ar(world, filename.c_str());

        int loaded_k;
        bool _has_y, _has_beta;
        long na;

        ar & loaded_k & _has_y & _has_beta & na;

        ResponseState state;

        state.x_alpha.resize(na);
        for (long i = 0; i < na; i++) ar & state.x_alpha[i];

        if (_has_y) {
            state.y_alpha.resize(na);
            for (long i = 0; i < na; i++) ar & state.y_alpha[i];
        }

        if (_has_beta) {
            long nb;
            ar & nb;

            state.x_beta.resize(nb);
            for (long i = 0; i < nb; i++) ar & state.x_beta[i];

            if (_has_y) {
                state.y_beta.resize(nb);
                for (long i = 0; i < nb; i++) ar & state.y_beta[i];
            }
        }

        // Reconcile wavelet order if needed
        int current_k = FunctionDefaults<NDIM>::get_k();
        double thresh = FunctionDefaults<NDIM>::get_thresh();
        if (loaded_k != current_k) {
            auto f = state.flat();
            reconstruct(world, f);
            for (auto& fn : f) {
                fn = project(fn, current_k, thresh, true);
            }
            truncate(world, f, thresh);
            state.from_flat(f);
        }

        return state;
    }

    // ---- Factory ----

    static ResponseState allocate(World& world, long n_alpha, long n_beta,
                                   bool include_y) {
        ResponseState state;
        state.x_alpha = zero_functions<T, NDIM>(world, n_alpha);
        if (include_y) {
            state.y_alpha = zero_functions<T, NDIM>(world, n_alpha);
        }
        if (n_beta > 0) {
            state.x_beta = zero_functions<T, NDIM>(world, n_beta);
            if (include_y) {
                state.y_beta = zero_functions<T, NDIM>(world, n_beta);
            }
        }
        return state;
    }

    // ---- Info ----

    void print_info(const std::string& label = "") const {
        if (!label.empty()) print("ResponseState:", label);
        print("  x_alpha:", num_alpha(), "functions");
        if (has_y()) print("  y_alpha:", num_alpha(), "functions");
        if (!is_restricted()) {
            print("  x_beta: ", num_beta(), "functions");
            if (has_y()) print("  y_beta: ", num_beta(), "functions");
        }
        print("  restricted:", is_restricted(), " has_y:", has_y(),
              " total:", total_size());
    }
};

// Convenience aliases
using RealResponseState = ResponseState<double, 3>;
using ComplexResponseState = ResponseState<std::complex<double>, 3>;

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_RESPONSEFUNCTIONS_HPP
