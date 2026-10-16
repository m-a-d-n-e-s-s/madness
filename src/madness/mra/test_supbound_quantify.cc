/// \file test_supbound_quantify.cc
/// \brief The sup bound of supbound.h on model functions, and the stopping error of the screened multiplication
/// against today's estimate and the Hoelder form: the measurement before the Hoelder screen is implemented
/// (~/exchange_screening/HANDOFF_mulsparse_hoelder.md section 2; it decides the leaf value, D1, and the second term, D3).
///
/// Model functions, normalized to ||f|| = 1 like orbitals: Gaussians exp(-r^2/w^2), w = 0.3, 1, 3 bohr; cusps
/// exp(-Z r), Z = 1, 8, 35 (H, O, Br); an orbital-like exp(-8 r) + 0.25 exp(-r^2/4); a p-like x exp(-r^2). Each is
/// centred on a box corner (the origin, a corner at every level), at the centre of the box at the level of the finest
/// leaf at the origin in the corner placement, and off the grid. Settings (k, thresh) = (6, 1e-4), (8, 1e-6),
/// (10, 1e-8) in the cell [-20, 20]^3 with truncate mode 1, projected and truncated as moldft's orbitals are.
///
/// (i)   Per leaf: M, the largest |f| at the k^3 quadrature points (coeffs2values); S, the coefficient sum
///       (SupBoundMap::leaf_bound); T, the largest |f| of the leaf's polynomial on a 24^3 tensor grid with the faces
///       and corners (T24), refined three times around its argmax on 9^3 grids of a quarter of the previous spacing
///       (Tref; T is the largest of T24, Tref and M); A, the largest |f| of the analytic function on the dense grid, and
///       dev, the largest pointwise deviation of the polynomial from it (the representation check).
/// (ii)  Per interior box of the redundant tree: S_own, M_own and T_own of the box's own polynomial (its scaling
///       coefficients) against sup_tree, the max of S over the leaves below (sup_bound_map()), and against the max of
///       M, T and A over the leaves below.
/// (iii) Per box of the union tree of a pair: the stopping error E(B) = ||f g - s_f s_g||_B, by Gauss quadrature with
///       2k points per dimension on every union leaf below B (exact for the represented functions: E_rep; against the
///       analytic product: E_ana), against today's estimate of mulXXveca and four Hoelder forms (f on the left, g on
///       the right, d = dnorm_tree, sup = sup_tree): sup_M(f) d_g + d_f sup_M(g) with M at the leaves; sup(f) d_g +
///       d_f sup(g) (D3 as decided); sup(f) d_g + d_f S(own g), S of g's coefficients on the box (D3's fallback); and
///       S(own f) d_g + d_f sup(g), from the mirrored split f g - s_f s_g = d_f g + s_f d_g. The last two are bounds.
///       The boxes where each criterion stops the descent at the development budget truncate_tol(0.01 thresh, B). For
///       the record, the collocation error of do_mul, C(B) = ||s_f s_g - q_B||, and the kernel's error if it formed
///       the product on B, K(B) = ||f g - q_B||, with q_B the collocated product, are outside every estimate.
///
/// Writes functions.csv, leaves.csv, interior.csv and union.csv into --out=DIR and summary tables to stdout. Checked
/// (exit code 1 if one fails beyond roundoff): no value of a polynomial exceeds its S (leaves and own polynomials);
/// both (S, S_own) estimates are at least E_rep on every box; and the kernel itself agrees with the union tree:
/// mul_sparse at moldft's tolerance forms its products exactly at today's stops, at tol = 0 exactly at the union
/// leaves, and the norms of their difference and of the exact product are what the quadrature gives.
/// Options: --out=DIR (write the CSV files there; none without it), --k=6,8,10 (the settings to run), --quick (two
/// functions, two placements and one pair).

#include <madness/mra/mra.h>
#include <madness/mra/legendre.h>
#include <madness/mra/supbound.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace madness;

namespace {

const std::size_t D = 3;
typedef Key<D> keyT;
typedef Vector<double, D> coordT;
typedef Function<double, D> functionT;
typedef FunctionImpl<double, D> implT;
typedef FunctionNode<double, D> nodeT;
typedef Range<int> rangeT;
typedef std::unordered_map<keyT, double, Hash<keyT>> keymapT;

/// the cell is [-HALF, HALF]^3
const double HALF = 20.0;
/// points per dimension of the dense grid on which T, T_own and A are read: j/23, faces and corners included
const int NDENSE = 24;
/// the kernel's development budget: the callers pass 0.1 thresh, mulXXvec multiplies by MUL_SCREENING_SAFETY = 0.1
const double BUDGET = 0.01;
/// the orbital-like model: exp(-Z r) + ORB_C exp(-r^2 / ORB_W^2)
const double ORB_C = 0.25, ORB_W = 2.0;
/// relative slack of the inequalities checked as theorems (roundoff)
const double SLACK = 1.e-9;

coordT vec(const double x, const double y, const double z) {
    coordT r;
    r[0] = x;
    r[1] = y;
    r[2] = z;
    return r;
}

/// the off-grid centre: no coordinate is a dyadic fraction of the cell
const coordT OFFGRID = vec(0.3183098861837907, -0.2718281828459045, 0.1414213562373095);
/// the direction from the first to the second function of a pair
const coordT PAIRDIR = vec(0.6, 0.64, 0.48);

enum class Family { gauss, cusp, orbital, p };

struct ModelSpec {
    std::string name;
    Family family;
    double a;   ///< the width w (gauss), Z (cusp, orbital); the p-like function has width 1
};

/// a model function centred at c, normalized to ||f|| = 1
class Model : public FunctionFunctorInterface<double, D> {
    Family family;
    double a, norm;
    coordT c;
    Level lev;

public:
    Model(const ModelSpec& s, const coordT& c) : family(s.family), a(s.a), norm(1.0), c(c), lev(2) {
        const double pi = constants::pi;
        double n2 = 1.0, scale = 1.0;   // the squared norm; the smallest length scale, for the special level
        switch (family) {
            case Family::gauss: n2 = std::pow(0.5 * pi * a * a, 1.5); scale = a; break;
            case Family::cusp: n2 = pi / (a * a * a); scale = 1.0 / a; break;
            case Family::orbital: {
                // the cross term 2 C 4 pi int r^2 exp(-Z r - r^2/W^2) dr by Simpson's rule on [0, 40] bohr
                const int n = 400000;
                const double h = 40.0 / n;
                double s = 0.0;
                for (int i = 0; i <= n; ++i) {
                    const double r = i * h;
                    const double w = (i == 0 or i == n) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
                    s += w * r * r * std::exp(-a * r - r * r / (ORB_W * ORB_W));
                }
                s *= h / 3.0;
                n2 = pi / (a * a * a) + ORB_C * ORB_C * std::pow(0.5 * pi * ORB_W * ORB_W, 1.5) + 8.0 * pi * ORB_C * s;
                scale = 1.0 / a;
                break;
            }
            case Family::p: n2 = 0.25 * std::pow(0.5 * pi, 1.5); scale = 1.0; break;
        }
        norm = 1.0 / std::sqrt(n2);
        // the level whose boxes are about as wide as the smallest length scale
        const double l = std::ceil(std::log2(2.0 * HALF / scale));
        lev = Level(std::min(double(FunctionDefaults<D>::get_max_refine_level()), std::max(2.0, l)));
    }

    double operator()(const coordT& r) const override {
        const double x = r[0] - c[0], y = r[1] - c[1], z = r[2] - c[2];
        const double rr2 = x * x + y * y + z * z;
        switch (family) {
            case Family::gauss: return norm * std::exp(-rr2 / (a * a));
            case Family::cusp: return norm * std::exp(-a * std::sqrt(rr2));
            case Family::orbital:
                return norm * (std::exp(-a * std::sqrt(rr2)) + ORB_C * std::exp(-rr2 / (ORB_W * ORB_W)));
            case Family::p: return norm * x * std::exp(-rr2);
        }
        return 0.0;
    }

    std::vector<coordT> special_points() const override { return {c}; }
    Level special_level() const override { return lev; }
    const coordT& centre() const { return c; }
};

std::vector<ModelSpec> all_models() {
    return {{"gauss_w0.3", Family::gauss, 0.3}, {"gauss_w1", Family::gauss, 1.0}, {"gauss_w3", Family::gauss, 3.0},
            {"cusp_Z1", Family::cusp, 1.0},     {"cusp_Z8", Family::cusp, 8.0},   {"cusp_Z35", Family::cusp, 35.0},
            {"orbital", Family::orbital, 8.0},  {"p", Family::p, 1.0}};
}

ModelSpec model_named(const std::string& name) {
    for (const auto& m : all_models())
        if (m.name == name) return m;
    MADNESS_EXCEPTION("no such model", 0);
}

struct PairSpec {
    std::string name;
    ModelSpec f, g;
    double distance;    ///< g's centre is f's (OFFGRID) plus distance * PAIRDIR
};

std::vector<PairSpec> all_pairs() {
    return {{"gauss_w1*cusp_Z35", model_named("gauss_w1"), model_named("cusp_Z35"), 1.5},
            {"cusp_Z8*cusp_Z8", model_named("cusp_Z8"), model_named("cusp_Z8"), 2.0},
            {"orbital*orbital", model_named("orbital"), model_named("orbital"), 2.5}};
}

/// the width of a box at `level` (the cell is cubic)
double box_width(const int level) { return 2.0 * HALF * std::pow(0.5, level); }

/// the lower corner of the box `key` in user coordinates
coordT box_lo(const keyT& key) {
    const double h = box_width(key.level());
    coordT lo;
    for (std::size_t d = 0; d < D; ++d) lo[d] = -HALF + h * double(key.translation()[d]);
    return lo;
}

/// whether the closed box `key` contains the point r
bool box_contains(const keyT& key, const coordT& r) {
    const double h = box_width(key.level()), eps = 1.e-12 * h;
    const coordT lo = box_lo(key);
    for (std::size_t d = 0; d < D; ++d)
        if (r[d] < lo[d] - eps or r[d] > lo[d] + h + eps) return false;
    return true;
}

/// phi(i, j) = phi_i(u_j), the normalized Legendre scaling functions at the box coordinates u
Tensor<double> scaling_matrix(const int k, const std::vector<double>& u) {
    Tensor<double> phi(long(k), long(u.size()));
    std::vector<double> p(k);
    for (std::size_t j = 0; j < u.size(); ++j) {
        legendre_scaling_functions(u[j], k, p.data());
        for (int i = 0; i < k; ++i) phi(long(i), long(j)) = p[i];
    }
    return phi;
}

/// the polynomial with scaling coefficients c on a box at `level` on the tensor grid of phi (one matrix per dimension)
Tensor<double> poly_values(const Tensor<double>& c, const int level, const std::array<Tensor<double>, D>& phi) {
    Tensor<double> v = general_transform(c, phi.data());
    v.scale(std::pow(2.0, 0.5 * D * level) / std::sqrt(FunctionDefaults<D>::get_cell_volume()));
    return v;
}

/// the box coordinates of the dense grid: NDENSE points, faces included
std::vector<double> dense_points() {
    std::vector<double> u(NDENSE);
    for (int j = 0; j < NDENSE; ++j) u[j] = double(j) / double(NDENSE - 1);
    return u;
}

/// the largest |p| of a polynomial on the dense grid (phi24 from dense_points())
double dense_max(const Tensor<double>& c, const int level, const Tensor<double>& phi24) {
    return poly_values(c, level, {phi24, phi24, phi24}).absmax();
}

/// median, 10th and 90th percentile (linear interpolation), min and max
struct Dist {
    long n = 0;
    double med = 0, p10 = 0, p90 = 0, min = 0, max = 0;
};

Dist dist(std::vector<double> v) {
    Dist r;
    r.n = long(v.size());
    if (v.empty()) return r;
    std::sort(v.begin(), v.end());
    auto q = [&v](const double x) {
        const double pos = x * double(v.size() - 1);
        const std::size_t i = std::size_t(pos);
        if (i + 1 >= v.size()) return v.back();
        return v[i] + (pos - double(i)) * (v[i + 1] - v[i]);
    };
    r.med = q(0.5);
    r.p10 = q(0.1);
    r.p90 = q(0.9);
    r.min = v.front();
    r.max = v.back();
    return r;
}

/// the share of v below x
double share_below(const std::vector<double>& v, const double x) {
    if (v.empty()) return 0.0;
    long n = 0;
    for (double y : v) n += (y < x);
    return double(n) / double(v.size());
}

struct Setting {
    int k;
    double thresh;
};

Setting setting_for(const int k) {
    if (k == 6) return {6, 1.e-4};
    if (k == 8) return {8, 1.e-6};
    if (k == 10) return {10, 1.e-8};
    MADNESS_EXCEPTION("k must be 6, 8 or 10", k);
}

/// the protocol of moldft (SCF::set_protocol, truncate mode 1 as SCF sets it)
void set_defaults(const Setting& s) {
    FunctionDefaults<D>::set_k(s.k);
    FunctionDefaults<D>::set_thresh(s.thresh);
    FunctionDefaults<D>::set_refine(true);
    FunctionDefaults<D>::set_initial_level(2);
    FunctionDefaults<D>::set_truncate_mode(1);
    FunctionDefaults<D>::set_autorefine(false);
    FunctionDefaults<D>::set_apply_randomize(false);
    FunctionDefaults<D>::set_project_randomize(false);
    FunctionDefaults<D>::set_cubic_cell(-HALF, HALF);
}

functionT project(World& world, const std::shared_ptr<Model>& m) {
    const std::shared_ptr<FunctionFunctorInterface<double, D>> functor = m;
    functionT f = FunctionFactory<double, D>(world).functor(functor);
    f.truncate();
    f.reconstruct();
    return f;
}

/// the theorems checked: violations and the worst ratio seen
struct Theorems {
    long leaf = 0, own = 0, stop = 0, kernel = 0;
    double worst_leaf = 0, worst_own = 0, worst_stop = 0;
};

/// one row of functions.csv
struct FunctionRecord {
    Setting s;
    std::string name, place;
    coordT c;
    long leaves = 0, interior = 0;
    int lmin = 99, lmax = 0;
    int level_at_c = -1;    ///< the deepest leaf level whose box contains c
    std::array<double, D> frac_c{{-1, -1, -1}};    ///< c's position in the deepest leaf containing it
    double norm = 0, t_project = 0, t_S = 0, t_M = 0, t_map = 0;
};

struct LeafRow {
    keyT key;
    double M = 0, S = 0, T24 = 0, Tref = 0, A = 0, dev = 0, norm = 0;
    bool at_c = false;
    double T() const { return std::max(std::max(T24, Tref), M); }
};

/// (i) the leaves of f
std::vector<LeafRow> measure_leaves(World& world, const functionT& f, const Model& m, FunctionRecord& rec,
                                    Theorems& th) {
    const implT& impl = *f.get_impl();
    const int k = impl.get_k();
    const Tensor<double> weights = SupBoundMap<D>::leaf_weights(k);
    std::vector<keyT> keys;
    std::vector<Tensor<double>> coeffs;
    for (auto it = impl.get_coeffs().begin(); it != impl.get_coeffs().end(); ++it) {
        if (it->second.has_children()) {
            ++rec.interior;
            continue;
        }
        keys.push_back(it->first);
        coeffs.push_back(it->second.has_coeff() ? copy(it->second.coeff().full_tensor()) : Tensor<double>());
    }
    const int n = int(keys.size());
    rec.leaves = n;

    // the cost of S and M per leaf, one thread
    double sink = 0.0, t0 = wall_time();
    for (int i = 0; i < n; ++i)
        if (coeffs[i].size()) sink += SupBoundMap<D>::leaf_bound(keys[i], coeffs[i], weights);
    rec.t_S = 1.e9 * (wall_time() - t0) / std::max(1, n);
    t0 = wall_time();
    for (int i = 0; i < n; ++i)
        if (coeffs[i].size()) sink += impl.coeffs2values(keys[i], coeffs[i]).absmax();
    rec.t_M = 1.e9 * (wall_time() - t0) / std::max(1, n);
    if (sink < 0.0) print("never");

    std::vector<LeafRow> rows(n);
    const std::vector<double> u = dense_points();
    const Tensor<double> phi24 = scaling_matrix(k, u);
    auto leaf_op = [&](const rangeT::iterator& it) -> bool {
        const int i = it;
        LeafRow& r = rows[i];
        const keyT& key = keys[i];
        const Tensor<double>& c = coeffs[i];
        r.key = key;
        r.at_c = box_contains(key, m.centre());
        const coordT lo = box_lo(key);
        const double h = box_width(key.level());
        Tensor<double> v;
        if (c.size()) {
            r.S = SupBoundMap<D>::leaf_bound(key, c, weights);
            r.M = impl.coeffs2values(key, c).absmax();
            r.norm = c.normf();
            v = poly_values(c, key.level(), {phi24, phi24, phi24});
            long qmax = 0;
            double vmax = -1.0;
            const double* pv = v.ptr();
            for (long q = 0; q < v.size(); ++q)
                if (std::abs(pv[q]) > vmax) {
                    vmax = std::abs(pv[q]);
                    qmax = q;
                }
            r.T24 = vmax;
            // refine around the argmax, three times: 9 points per dimension spanning one step either side, at a
            // quarter of the step each time (the last spacing is 1/64 of the dense grid's)
            std::array<double, D> x0{{u[qmax / (NDENSE * NDENSE)], u[(qmax / NDENSE) % NDENSE], u[qmax % NDENSE]}};
            double step = 1.0 / double(NDENSE - 1);
            r.Tref = vmax;
            for (int stage = 0; stage < 3; ++stage) {
                std::array<std::vector<double>, D> ul;
                std::array<Tensor<double>, D> phl;
                for (std::size_t d = 0; d < D; ++d) {
                    for (int s = -4; s <= 4; ++s) ul[d].push_back(std::min(1.0, std::max(0.0, x0[d] + 0.25 * s * step)));
                    phl[d] = scaling_matrix(k, ul[d]);
                }
                const Tensor<double> vl = poly_values(c, key.level(), phl);
                const double* pl = vl.ptr();
                for (long q = 0; q < vl.size(); ++q)
                    if (std::abs(pl[q]) > r.Tref) {
                        r.Tref = std::abs(pl[q]);
                        x0 = {{ul[0][q / 81], ul[1][(q / 9) % 9], ul[2][q % 9]}};
                    }
                step *= 0.25;
            }
        }
        double amax = 0.0, dev = 0.0;
        for (int a = 0; a < NDENSE; ++a)
            for (int b = 0; b < NDENSE; ++b)
                for (int e = 0; e < NDENSE; ++e) {
                    const double fa = m(vec(lo[0] + h * u[a], lo[1] + h * u[b], lo[2] + h * u[e]));
                    amax = std::max(amax, std::abs(fa));
                    dev = std::max(dev, std::abs((c.size() ? v(a, b, e) : 0.0) - fa));
                }
        r.A = amax;
        r.dev = dev;
        return true;
    };
    world.taskq.for_each(rangeT(0, n), leaf_op).get();

    for (const auto& r : rows) {
        rec.lmin = std::min(rec.lmin, int(r.key.level()));
        rec.lmax = std::max(rec.lmax, int(r.key.level()));
        if (r.at_c and int(r.key.level()) > rec.level_at_c) {
            rec.level_at_c = r.key.level();
            const coordT lo = box_lo(r.key);
            for (std::size_t d = 0; d < D; ++d) rec.frac_c[d] = (m.centre()[d] - lo[d]) / box_width(r.key.level());
        }
        const double worst = std::max(std::max(r.T24, r.Tref), r.M);
        if (worst > r.S * (1.0 + 1.e-12)) ++th.leaf;
        if (r.S > 0.0) th.worst_leaf = std::max(th.worst_leaf, worst / r.S);
    }
    return rows;
}

struct InteriorRow {
    keyT key;
    int deepest = 0;
    double Sown = 0, Mown = 0, Town = 0, supS = 0, supM = 0, Ttree = 0, Atree = 0, nt = 0, dnt = 0;
};

/// (ii) the interior boxes of the redundant tree of f
std::vector<InteriorRow> measure_interior(World& world, const functionT& f, const std::vector<LeafRow>& leaves,
                                          FunctionRecord& rec, Theorems& th) {
    const implT& impl = *f.get_impl();
    MADNESS_CHECK(impl.is_redundant());
    const int k = impl.get_k();
    const Tensor<double> weights = SupBoundMap<D>::leaf_weights(k);
    double t0 = wall_time();
    const SupBoundMap<D> supS = impl.sup_bound_map();
    rec.t_map = 1.e3 * (wall_time() - t0);

    // the max of M, T, A and the deepest level over the leaves below, propagated from the leaves up
    keymapT Mt, Tt, At, deep;
    for (const auto& r : leaves) {
        Mt[r.key] = r.M;
        Tt[r.key] = r.T();
        At[r.key] = r.A;
        deep[r.key] = r.key.level();
    }
    std::vector<InteriorRow> rows;
    std::vector<Tensor<double>> coeffs;
    for (auto it = impl.get_coeffs().begin(); it != impl.get_coeffs().end(); ++it) {
        if (not it->second.has_children()) continue;
        InteriorRow r;
        r.key = it->first;
        r.nt = it->second.get_norm_tree();
        r.dnt = it->second.get_dnorm_tree();
        rows.push_back(r);
        MADNESS_CHECK(it->second.has_coeff());
        coeffs.push_back(copy(it->second.coeff().full_tensor()));
    }
    std::vector<int> order(rows.size());
    for (std::size_t i = 0; i < rows.size(); ++i) order[i] = int(i);
    std::sort(order.begin(), order.end(), [&rows](int a, int b) { return rows[a].key.level() > rows[b].key.level(); });
    for (int i : order) {
        InteriorRow& r = rows[i];
        double m = 0.0, t = 0.0, a = 0.0, dl = 0.0;
        for (KeyChildIterator<D> kit(r.key); kit; ++kit) {
            m = std::max(m, Mt.at(kit.key()));
            t = std::max(t, Tt.at(kit.key()));
            a = std::max(a, At.at(kit.key()));
            dl = std::max(dl, deep.at(kit.key()));
        }
        Mt[r.key] = r.supM = m;
        Tt[r.key] = r.Ttree = t;
        At[r.key] = r.Atree = a;
        deep[r.key] = dl;
        r.deepest = int(dl);
        r.supS = supS(r.key);
    }

    const Tensor<double> phi24 = scaling_matrix(k, dense_points());
    auto own_op = [&](const rangeT::iterator& it) -> bool {
        const int i = it;
        InteriorRow& r = rows[i];
        r.Sown = SupBoundMap<D>::leaf_bound(r.key, coeffs[i], weights);
        r.Mown = impl.coeffs2values(r.key, coeffs[i]).absmax();
        r.Town = dense_max(coeffs[i], r.key.level(), phi24);
        return true;
    };
    world.taskq.for_each(rangeT(0, int(rows.size())), own_op).get();
    for (const auto& r : rows) {
        const double worst = std::max(r.Mown, r.Town);
        if (worst > r.Sown * (1.0 + 1.e-12)) ++th.own;
        if (r.Sown > 0.0) th.worst_own = std::max(th.worst_own, worst / r.Sown);
    }
    return rows;
}

/// the criteria of the union-tree measurement, in the order of UBox::est: today's, the Hoelder form with M leaves,
/// with sup_tree in both terms, with S of the right operand's own coefficients, with S of the left operand's
const int NCRIT = 5;
const char* crit_name[NCRIT] = {"cs", "hM", "hS", "hSown", "hSown_gf"};

/// a box of the union tree of f and g
struct UBox {
    keyT key;
    int parent = -1;
    bool in_f = false, in_g = false, leaf = false;
    Tensor<double> sf, sg;  ///< scaling coefficients on the box; restricted from the leaf above if the box is below it
    Tensor<double> qc;      ///< the product do_mul forms on the box from sf and sg (collocation at the k^3 points)
    double nf = 0, ndf = 0, ng = 0, ndg = 0;            ///< norm_tree, dnorm_tree (0 below a leaf)
    double supSf = 0, supMf = 0, supSg = 0, supMg = 0;  ///< max of S, M over the leaves below; of the restriction below a leaf
    double Sf = 0, Mf = 0, Tf = 0, Sg = 0, Mg = 0, Tg = 0;  ///< S, M and the dense-grid max of the box's own polynomials
    double theta = 0;
    /// squared, summed over the union leaves below: the stopping error ||f g - s_f s_g|| of the represented and of
    /// the analytic functions, the projection floor (E_ana of the union leaves), ||f g||, the collocation error
    /// ||s_f s_g - qc|| and the kernel's error ||f g - qc|| if it formed the product here (represented, analytic)
    double e2rep = 0, e2ana = 0, floor2 = 0, n2fg = 0, c2 = 0, k2rep = 0, k2ana = 0;
    /// squared, summed over the union leaves below: the product formed here against the products formed on the
    /// union leaves (the kernel at tol = 0); over the stops of today's criterion it is ||mul_sparse(tol) - mul_sparse(0)||^2
    double d2q = 0;
    std::array<double, NCRIT> est{};
    std::array<char, NCRIT> stop{}, above{};
};

/// the max of M over the leaves below each box of f's tree (reconstructed or redundant)
keymapT M_tree(const implT& impl) {
    keymapT t;
    std::vector<keyT> interior;
    for (auto it = impl.get_coeffs().begin(); it != impl.get_coeffs().end(); ++it) {
        if (it->second.has_children())
            interior.push_back(it->first);
        else
            t[it->first] =
                    it->second.has_coeff() ? impl.coeffs2values(it->first, it->second.coeff().full_tensor()).absmax() : 0.0;
    }
    std::sort(interior.begin(), interior.end(), [](const keyT& a, const keyT& b) { return a.level() > b.level(); });
    for (const keyT& key : interior) {
        double m = 0.0;
        for (KeyChildIterator<D> kit(key); kit; ++kit) m = std::max(m, t.at(kit.key()));
        t[key] = m;
    }
    return t;
}

/// (iii) the union tree of a pair
std::vector<UBox> measure_union(World& world, functionT& f, functionT& g, const Model& fm, const Model& gm,
                                const Setting& s, Theorems& th) {
    const implT& fi = *f.get_impl();
    const implT& gi = *g.get_impl();
    const int k = fi.get_k();
    const keymapT Mf = M_tree(fi), Mg = M_tree(gi);
    const SupBoundMap<D> Sf = fi.sup_bound_map(), Sg = gi.sup_bound_map();
    f.make_redundant(true);
    g.make_redundant(true);
    const auto& cdata = fi.get_cdata();

    auto node = [](const implT& impl, const keyT& key) -> const nodeT* {
        const auto it = impl.get_coeffs().find(key).get();
        return (it == impl.get_coeffs().end()) ? nullptr : &it->second;
    };
    // breadth first: every parent precedes its children
    std::vector<UBox> U(1);
    {
        UBox& r = U[0];
        r.key = cdata.key0;
        const nodeT* nf = node(fi, r.key);
        const nodeT* ng = node(gi, r.key);
        MADNESS_CHECK(nf and ng);
        r.in_f = r.in_g = true;
        r.sf = copy(nf->coeff().full_tensor());
        r.sg = copy(ng->coeff().full_tensor());
        r.nf = nf->get_norm_tree();
        r.ndf = nf->get_dnorm_tree();
        r.ng = ng->get_norm_tree();
        r.ndg = ng->get_dnorm_tree();
    }
    for (std::size_t i = 0; i < U.size(); ++i) {
        const keyT key = U[i].key;
        const bool fch = U[i].in_f and node(fi, key)->has_children();
        const bool gch = U[i].in_g and node(gi, key)->has_children();
        if (not fch and not gch) {
            U[i].leaf = true;
            continue;
        }
        // the children's coefficients of an operand that has no node below this box: one level of unfilter
        Tensor<double> fss, gss;
        if (not fch) {
            Tensor<double> d(cdata.v2k);
            d(cdata.s0) = U[i].sf;
            fss = fi.unfilter(d);
        }
        if (not gch) {
            Tensor<double> d(cdata.v2k);
            d(cdata.s0) = U[i].sg;
            gss = gi.unfilter(d);
        }
        for (KeyChildIterator<D> kit(key); kit; ++kit) {
            UBox b;
            b.key = kit.key();
            b.parent = int(i);
            if (fch) {
                const nodeT* n = node(fi, b.key);
                MADNESS_CHECK(n and n->has_coeff() and n->get_dnorm_tree() < NORM_TREE_UNCOMPUTED);
                b.in_f = true;
                b.sf = copy(n->coeff().full_tensor());
                b.nf = n->get_norm_tree();
                b.ndf = n->get_dnorm_tree();
            } else {
                b.sf = copy(fss(fi.child_patch(b.key)));
                b.nf = b.sf.normf();
                b.ndf = 0.0;
            }
            if (gch) {
                const nodeT* n = node(gi, b.key);
                MADNESS_CHECK(n and n->has_coeff() and n->get_dnorm_tree() < NORM_TREE_UNCOMPUTED);
                b.in_g = true;
                b.sg = copy(n->coeff().full_tensor());
                b.ng = n->get_norm_tree();
                b.ndg = n->get_dnorm_tree();
            } else {
                b.sg = copy(gss(gi.child_patch(b.key)));
                b.ng = b.sg.normf();
                b.ndg = 0.0;
            }
            U.push_back(std::move(b));
        }
    }

    // the bounds and the estimates per box
    const Tensor<double> weights = SupBoundMap<D>::leaf_weights(k);
    const Tensor<double> phi24 = scaling_matrix(k, dense_points());
    auto box_op = [&](const rangeT::iterator& it) -> bool {
        UBox& b = U[int(it)];
        const int n = b.key.level();
        b.Sf = SupBoundMap<D>::leaf_bound(b.key, b.sf, weights);
        b.Mf = fi.coeffs2values(b.key, b.sf).absmax();
        b.Tf = dense_max(b.sf, n, phi24);
        b.Sg = SupBoundMap<D>::leaf_bound(b.key, b.sg, weights);
        b.Mg = gi.coeffs2values(b.key, b.sg).absmax();
        b.Tg = dense_max(b.sg, n, phi24);
        // do_mul: the values at the quadrature points multiplied and projected back
        Tensor<double> vfg = fi.coeffs2values(b.key, b.sf);
        vfg.emul(gi.coeffs2values(b.key, b.sg));
        b.qc = transform(vfg, cdata.quad_phiw).scale(std::pow(0.5, 0.5 * D * n) *
                                                     std::sqrt(FunctionDefaults<D>::get_cell_volume()));
        b.supSf = b.in_f ? Sf(b.key) : b.Sf;
        b.supMf = b.in_f ? Mf.at(b.key) : b.Mf;
        b.supSg = b.in_g ? Sg(b.key) : b.Sg;
        b.supMg = b.in_g ? Mg.at(b.key) : b.Mg;
        // in the kernel's own floating-point order: the caller's 0.1 thresh, times MUL_SCREENING_SAFETY in mulXXvec
        b.theta = fi.truncate_tol((0.1 * s.thresh) * MUL_SCREENING_SAFETY, b.key);
        b.est[0] = b.ng * b.ndf + b.nf * b.ndg + b.ndf * b.ndg;     // mulXXveca today, f on the left
        b.est[1] = b.supMf * b.ndg + b.ndf * b.supMg;
        b.est[2] = b.supSf * b.ndg + b.ndf * b.supSg;
        b.est[3] = b.supSf * b.ndg + b.ndf * b.Sg;                  // f g - s_f s_g = f d_g + d_f s_g
        b.est[4] = b.supSg * b.ndf + b.ndg * b.Sf;                  // f g - s_f s_g = d_f g + s_f d_g
        return true;
    };
    world.taskq.for_each(rangeT(0, int(U.size())), box_op).get();

    // E(B): Gauss quadrature with 2k points per dimension on every union leaf, carried up to each ancestor
    const int npt = 2 * k;
    std::vector<double> tq(npt), wq(npt);
    MADNESS_CHECK(gauss_legendre(npt, 0.0, 1.0, tq.data(), wq.data()));
    const Tensor<double> phiq = scaling_matrix(k, tq);
    std::vector<int> leaves;
    for (std::size_t i = 0; i < U.size(); ++i)
        if (U[i].leaf) leaves.push_back(int(i));
    std::vector<std::vector<std::array<double, 7>>> contrib(leaves.size());
    auto leaf_op = [&](const rangeT::iterator& it) -> bool {
        const int li = it;
        const UBox& lb = U[leaves[li]];
        const int m = lb.key.level();
        const double h = box_width(m);
        const coordT lo = box_lo(lb.key);
        Tensor<double> W(npt, npt, npt), FGa(npt, npt, npt);
        for (int a = 0; a < npt; ++a)
            for (int b = 0; b < npt; ++b)
                for (int e = 0; e < npt; ++e) {
                    W(a, b, e) = wq[a] * wq[b] * wq[e] * h * h * h;
                    const coordT x = vec(lo[0] + h * tq[a], lo[1] + h * tq[b], lo[2] + h * tq[e]);
                    FGa(a, b, e) = fm(x) * gm(x);
                }
        Tensor<double> FG = poly_values(lb.sf, m, {phiq, phiq, phiq});
        FG.emul(poly_values(lb.sg, m, {phiq, phiq, phiq}));
        // the product the kernel forms on the leaf itself (its tol = 0 result there)
        const Tensor<double> QL = poly_values(lb.qc, m, {phiq, phiq, phiq});
        auto& out = contrib[li];
        for (int bi = leaves[li]; bi >= 0; bi = U[bi].parent) {
            const UBox& B = U[bi];
            const int n = B.key.level();
            // the leaf's points in the frame of the box B
            const double scale = std::pow(0.5, m - n);
            std::array<Tensor<double>, D> ph;
            for (std::size_t d = 0; d < D; ++d) {
                const double off = double(lb.key.translation()[d]) * scale - double(B.key.translation()[d]);
                std::vector<double> uu(npt);
                for (int j = 0; j < npt; ++j) uu[j] = off + scale * tq[j];
                ph[d] = scaling_matrix(k, uu);
            }
            Tensor<double> SS = poly_values(B.sf, n, ph);
            SS.emul(poly_values(B.sg, n, ph));
            const Tensor<double> QQ = poly_values(B.qc, n, ph);
            double er = 0.0, ea = 0.0, nfg = 0.0, cc = 0.0, kr = 0.0, ka = 0.0, dq = 0.0;
            const double *pw = W.ptr(), *pfg = FG.ptr(), *pfa = FGa.ptr(), *pss = SS.ptr(), *pq = QQ.ptr(),
                         *pl = QL.ptr();
            for (long q = 0; q < W.size(); ++q) {
                er += pw[q] * (pfg[q] - pss[q]) * (pfg[q] - pss[q]);
                ea += pw[q] * (pfa[q] - pss[q]) * (pfa[q] - pss[q]);
                nfg += pw[q] * pfg[q] * pfg[q];
                cc += pw[q] * (pss[q] - pq[q]) * (pss[q] - pq[q]);
                kr += pw[q] * (pfg[q] - pq[q]) * (pfg[q] - pq[q]);
                ka += pw[q] * (pfa[q] - pq[q]) * (pfa[q] - pq[q]);
                dq += pw[q] * (pq[q] - pl[q]) * (pq[q] - pl[q]);
            }
            out.push_back({er, ea, nfg, cc, kr, ka, dq});
        }
        return true;
    };
    world.taskq.for_each(rangeT(0, int(leaves.size())), leaf_op).get();
    for (std::size_t li = 0; li < leaves.size(); ++li) {
        std::size_t j = 0;
        for (int bi = leaves[li]; bi >= 0; bi = U[bi].parent, ++j) {
            U[bi].e2rep += contrib[li][j][0];
            U[bi].e2ana += contrib[li][j][1];
            U[bi].n2fg += contrib[li][j][2];
            U[bi].c2 += contrib[li][j][3];
            U[bi].k2rep += contrib[li][j][4];
            U[bi].k2ana += contrib[li][j][5];
            U[bi].d2q += contrib[li][j][6];
            U[bi].floor2 += contrib[li][0][1];
        }
    }

    // where each criterion stops the descent: the first box on the path down whose estimate is within the budget
    for (auto& b : U) {
        for (int c = 0; c < NCRIT; ++c) {
            b.above[c] = (b.parent >= 0) and (U[b.parent].stop[c] or U[b.parent].above[c]);
            b.stop[c] = (not b.above[c]) and (b.est[c] <= b.theta);
        }
        // the theorem: the (S, S_own) estimate bounds the stopping error of the represented functions
        for (int c : {3, 4}) {
            const double e = std::sqrt(b.e2rep);
            if (e > b.est[c] * (1.0 + SLACK) + 1.e-12 * std::sqrt(b.n2fg)) ++th.stop;
            if (b.est[c] > 0.0) th.worst_stop = std::max(th.worst_stop, e / b.est[c]);
        }
    }
    return U;
}

/// The kernel itself, against the union tree: mul_sparse at moldft's tolerance (today's criterion) must form its
/// products exactly at the stops of "cs", mul_sparse at tol = 0 exactly at the union leaves, and the difference of the
/// two products, and the norm of the exact one, must be what the quadrature says (all of it exact arithmetic up to
/// roundoff: the polynomials' squares are integrated exactly)
void check_kernel(const functionT& f, const functionT& g, const std::vector<UBox>& U, const Setting& s, Theorems& th) {
    functionT q = mul_sparse(f, g, 0.1 * s.thresh);
    functionT q0 = mul_sparse(f, g, 0.0);
    auto leaves_of = [](const functionT& h) {
        std::unordered_set<keyT, Hash<keyT>> l;
        for (auto it = h.get_impl()->get_coeffs().begin(); it != h.get_impl()->get_coeffs().end(); ++it)
            if (not it->second.has_children()) l.insert(it->first);
        return l;
    };
    const auto lq = leaves_of(q), l0 = leaves_of(q0);
    std::unordered_set<keyT, Hash<keyT>> stops, uleaves;
    double d2 = 0.0, k2 = 0.0;
    for (const auto& b : U) {
        if (b.stop[0]) {
            stops.insert(b.key);
            d2 += b.d2q;
        }
        if (b.leaf) {
            uleaves.insert(b.key);
            k2 += b.k2rep;
        }
    }
    auto differ = [](const std::unordered_set<keyT, Hash<keyT>>& a, const std::unordered_set<keyT, Hash<keyT>>& b) {
        long n = 0;
        for (const auto& x : a) n += (b.count(x) == 0);
        for (const auto& x : b) n += (a.count(x) == 0);
        return n;
    };
    const long miss = differ(lq, stops), miss0 = differ(l0, uleaves);
    const double n0 = q0.norm2(), nfg = std::sqrt(U[0].n2fg), kex = std::sqrt(k2);
    const double dk = (q - q0).norm2(), dme = std::sqrt(d2);
    const bool ok = miss == 0 and miss0 == 0 and std::abs(dk - dme) <= 1.e-6 * dk + 1.e-13 * n0 and
                    std::abs(n0 - nfg) <= kex * (1.0 + 1.e-6) + 1.e-13 * n0;
    std::printf("    kernel: mul_sparse(0.1 thresh) has %zu leaves, today's criterion %zu stops, %ld differ; mul_sparse(0) "
                "%zu leaves, union leaves %zu, %ld differ; ||q - q0|| %.6e, quadrature %.6e; ||q0|| %.10e against ||f g|| "
                "%.10e, within the collocation error %.2e: %s\n", lq.size(), stops.size(), miss, l0.size(),
                uleaves.size(), miss0, dk, dme, n0, nfg, kex, ok ? "agree" : "DISAGREE");
    if (not ok) ++th.kernel;
}

// ------------------------------------------------------------------ output

/// the CSV file `name` in `dir` with its header line; none without a directory
FILE* open_csv(const std::string& dir, const std::string& name, const char* header) {
    if (dir.empty()) return nullptr;
    const std::string path = dir + "/" + name;
    FILE* fp = std::fopen(path.c_str(), "w");
    MADNESS_CHECK_THROW(fp, ("cannot write " + path).c_str());
    std::fprintf(fp, "%s\n", header);
    return fp;
}

std::string key_cols(const keyT& key) {
    char buf[96];
    std::snprintf(buf, sizeof(buf), "%d,%ld,%ld,%ld", int(key.level()), long(key.translation()[0]),
                  long(key.translation()[1]), long(key.translation()[2]));
    return buf;
}

void print_dist_row(const char* label, const Dist& m, const Dist& s, const double f95, const double f80) {
    std::printf("    %5s %6ld   M/T %6.3f %6.3f %6.3f   S/T %6.3f %6.3f %6.3f   M<0.95T %5.3f  M<0.8T %5.3f\n", label, m.n,
                m.med, m.p10, m.min, s.med, s.p90, s.max, f95, f80);
}

/// the summary of (i) for one function: per level, the leaves with T > thresh
void summarize_leaves(const FunctionRecord& rec, const std::vector<LeafRow>& rows) {
    std::printf("  (i) %s %s: %ld leaves (levels %d-%d), %ld interior; leaf level at the centre %d at (%.3f %.3f %.3f);"
                " ||f|| %.6f; projected in %.2f s; S %.0f ns, M %.0f ns per leaf\n",
                rec.name.c_str(), rec.place.c_str(), rec.leaves, rec.lmin, rec.lmax, rec.interior, rec.level_at_c,
                rec.frac_c[0], rec.frac_c[1], rec.frac_c[2], rec.norm, rec.t_project, rec.t_S, rec.t_M);
    std::printf("    level  leaves(T>thresh)  M/T median p10 min      S/T median p90 max     shares of M below 0.95 T, 0.8 T\n");
    for (int l = rec.lmin - 1; l <= rec.lmax; ++l) {
        std::vector<double> mt, st;
        for (const auto& r : rows) {
            if (l >= rec.lmin and int(r.key.level()) != l) continue;
            if (r.T() <= rec.s.thresh) continue;
            mt.push_back(r.M / r.T());
            st.push_back(r.S / r.T());
        }
        if (mt.empty()) continue;
        const std::string label = (l < rec.lmin) ? "all" : std::to_string(l);
        print_dist_row(label.c_str(), dist(mt), dist(st), share_below(mt, 0.95), share_below(mt, 0.8));
    }
    // the leaf holding the function's largest value
    const LeafRow* peak = nullptr;
    double ta = 0.0, worst_dev = 0.0, tref = 0.0;
    for (const auto& r : rows) {
        if (not peak or r.T() > peak->T()) peak = &r;
        if (r.T() > rec.s.thresh) ta = std::max(ta, std::abs(r.T() / std::max(r.A, 1.e-300) - 1.0));
        worst_dev = std::max(worst_dev, r.dev);
        if (r.T24 > 0.0) tref = std::max(tref, r.Tref / r.T24 - 1.0);
    }
    if (peak)
        std::printf("    peak leaf: level %d, T %.4e, M/T %.3f, S/T %.3f; largest |T/A - 1| (T > thresh) %.2e, largest "
                    "|p - f| %.2f thresh, refinement raised T by at most %.2e\n",
                    int(peak->key.level()), peak->T(), peak->M / peak->T(), peak->S / peak->T(), ta,
                    worst_dev / rec.s.thresh, tref);
}

/// the summary of (ii) for one function: per level, interior boxes with T_tree > thresh, near (within two levels of
/// their deepest leaf) and far
void summarize_interior(const FunctionRecord& rec, const std::vector<InteriorRow>& rows) {
    std::printf("    (ii) interior boxes with T_tree > thresh: S_own/sup_tree median p90 max, T_own/sup_tree max, "
                "boxes with T_own > sup_tree; near = within two levels of the deepest leaf below\n");
    int lmax = 0;
    for (const auto& r : rows) lmax = std::max(lmax, int(r.key.level()));
    for (int l = -1; l <= lmax; ++l) {
        for (int near = 1; near >= 0; --near) {
            std::vector<double> so, to;
            for (const auto& r : rows) {
                if (l >= 0 and int(r.key.level()) != l) continue;
                if (r.Ttree <= rec.s.thresh) continue;
                if ((r.deepest - int(r.key.level()) <= 2) != bool(near)) continue;
                so.push_back(r.Sown / r.supS);
                to.push_back(r.Town / r.supS);
            }
            if (so.empty()) continue;
            const Dist a = dist(so), b = dist(to);
            long over = 0;
            for (double x : to) over += (x > 1.0 + 1.e-12);
            std::printf("    %5s %4s %6ld   S_own/sup_tree %6.3f %6.3f %6.3f   T_own/sup_tree max %6.3f  over %ld\n",
                        l < 0 ? "all" : std::to_string(l).c_str(), near ? "near" : "far", a.n, a.med, a.p90, a.max,
                        b.max, over);
        }
    }
}

/// the summary of (iii) for one pair
void summarize_union(const std::string& name, const Setting& s, const std::vector<UBox>& U) {
    int lmax = 0;
    long nleaf = 0;
    for (const auto& b : U) {
        lmax = std::max(lmax, int(b.key.level()));
        nleaf += b.leaf;
    }
    std::printf("  (iii) %s, k %d thresh %.0e: %zu union boxes, %ld union leaves, levels 0-%d; f (the first) on the "
                "left; hSown bounds g's own polynomial, hSown_gf f's\n", name.c_str(), s.k, s.thresh, U.size(), nleaf,
                lmax);
    std::printf("    criterion  stops  total E_rep  total E_ana  total K_rep  budget(stops)  max E_rep/theta  min est/E_rep at stops\n");
    for (int c = 0; c < NCRIT; ++c) {
        double e2 = 0.0, a2 = 0.0, k2 = 0.0, t2 = 0.0, worst = 0.0, minr = 1.e300;
        long n = 0;
        for (const auto& b : U) {
            if (not b.stop[c]) continue;
            ++n;
            e2 += b.e2rep;
            a2 += b.e2ana;
            k2 += b.k2rep;
            t2 += b.theta * b.theta;
            worst = std::max(worst, std::sqrt(b.e2rep) / b.theta);
            if (b.e2rep > 0.0) minr = std::min(minr, b.est[c] / std::sqrt(b.e2rep));
        }
        std::printf("    %-9s %6ld   %.3e    %.3e    %.3e    %.3e      %8.3f         %s\n", crit_name[c], n,
                    std::sqrt(e2), std::sqrt(a2), std::sqrt(k2), std::sqrt(t2), worst,
                    minr < 1.e300 ? std::to_string(minr).c_str() : "(all exact)");
    }
    // D3: at the boxes where (S, sup_tree) stops, the right operand's own polynomial against its sup_tree
    long n_s = 0, n_t = 0, m_s = 0, m_t = 0, nstop = 0;
    for (const auto& b : U) {
        if (not b.stop[2] or b.leaf) continue;
        ++nstop;
        n_s += (b.Sg > b.supSg * (1.0 + 1.e-12));
        n_t += (b.Tg > b.supSg * (1.0 + 1.e-12));
        m_s += (b.Sf > b.supSf * (1.0 + 1.e-12));
        m_t += (b.Tf > b.supSf * (1.0 + 1.e-12));
    }
    std::printf("    D3, at the %ld stops of hS off the union leaves: own S above sup_tree for g %ld (f %ld), own T above "
                "sup_tree for g %ld (f %ld)\n", nstop, n_s, m_s, n_t, m_t);
    std::printf("    level  boxes   est/E_rep median|min per criterion (boxes off the union leaves)        stops per criterion\n");
    std::printf("                ");
    for (int c = 0; c < NCRIT; ++c) std::printf("  %-17s", crit_name[c]);
    std::printf("\n");
    for (int l = 0; l <= lmax; ++l) {
        std::array<std::vector<double>, NCRIT> r;
        std::array<long, NCRIT> stops{};
        long nb = 0;
        for (const auto& b : U) {
            if (int(b.key.level()) != l) continue;
            for (int c = 0; c < NCRIT; ++c) stops[c] += b.stop[c];
            if (b.leaf or b.e2rep <= 0.0) continue;
            ++nb;
            for (int c = 0; c < NCRIT; ++c) r[c].push_back(b.est[c] / std::sqrt(b.e2rep));
        }
        std::printf("    %5d %6ld ", l, nb);
        for (int c = 0; c < NCRIT; ++c) {
            const Dist d = dist(r[c]);
            std::printf("  %8.2e|%8.2e", d.med, d.min);
        }
        std::printf("   ");
        for (int c = 0; c < NCRIT; ++c) std::printf(" %5ld", stops[c]);
        std::printf("\n");
    }
}

struct Options {
    std::string out;    ///< no CSV files without it
    std::vector<int> ks{6, 8, 10};
    bool quick = false;
};

Options parse(int argc, char** argv) {
    Options o;
    bool kset = false;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a.rfind("--out=", 0) == 0) o.out = a.substr(6);
        else if (a == "--quick") o.quick = true;
        else if (a.rfind("--k=", 0) == 0) {
            o.ks.clear();
            kset = true;
            std::string list = a.substr(4);
            for (std::size_t p = 0; p < list.size();) {
                const std::size_t q = std::min(list.find(',', p), list.size());
                o.ks.push_back(std::stoi(list.substr(p, q - p)));
                p = q + 1;
            }
        }
    }
    if (o.quick and not kset) o.ks = {6};
    return o;
}

}   // namespace

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    startup(world, argc, argv, true);
    const Options opt = parse(argc, argv);
    MADNESS_CHECK_THROW(world.size() == 1, "test_supbound_quantify runs on one rank");

    Theorems th;
    FILE* ffun = open_csv(opt.out, "functions.csv",
                          "k,thresh,function,placement,cx,cy,cz,leaves,interior,lmin,lmax,level_at_c,fx,fy,fz,norm,"
                          "t_project_s,t_S_ns,t_M_ns,t_map_ms");
    FILE* fleaf = open_csv(opt.out, "leaves.csv",
                           "k,thresh,function,placement,level,l0,l1,l2,at_c,M,S,T24,Tref,A,dev,norm");
    FILE* fint = open_csv(opt.out, "interior.csv",
                          "k,thresh,function,placement,level,l0,l1,l2,deepest,Sown,Mown,Town,supS,supM,Ttree,Atree,"
                          "norm_tree,dnorm_tree");
    FILE* fu = open_csv(opt.out, "union.csv",
                        "k,thresh,pair,level,l0,l1,l2,leaf,in_f,in_g,nf,ndf,ng,ndg,supSf,supMf,supSg,supMg,Sf,Mf,Tf,"
                        "Sg,Mg,Tg,theta,E_rep,E_ana,floor,n_fg,C,K_rep,K_ana,D_q,est_cs,est_hM,est_hS,est_hSown,"
                        "est_hSown_gf,stop_cs,stop_hM,stop_hS,stop_hSown,stop_hSown_gf");

    std::vector<ModelSpec> models = all_models();
    std::vector<std::string> places{"corner", "centre", "offgrid"};
    std::vector<PairSpec> pairs = all_pairs();
    if (opt.quick) {
        models = {model_named("gauss_w1"), model_named("cusp_Z8")};
        places = {"corner", "offgrid"};
        pairs = {pairs[1]};
    }

    for (const int kk : opt.ks) {
        const Setting s = setting_for(kk);
        set_defaults(s);
        std::printf("\n==== k %d, thresh %.0e, cell [%g, %g]^3, truncate mode 1, budget truncate_tol(%.2f thresh, B)\n",
                    s.k, s.thresh, -HALF, HALF, BUDGET);
        for (const auto& spec : models) {
            int level_corner = -1;  // the deepest leaf level at the origin in the corner placement
            for (const auto& place : places) {
                coordT c(0.0);
                if (place == "centre") {
                    MADNESS_CHECK(level_corner >= 0);
                    const double h = box_width(level_corner);
                    c = vec(0.5 * h, 0.5 * h, 0.5 * h);
                } else if (place == "offgrid") {
                    c = OFFGRID;
                }
                FunctionRecord rec;
                rec.s = s;
                rec.name = spec.name;
                rec.place = place;
                rec.c = c;
                const auto m = std::make_shared<Model>(spec, c);
                double t0 = wall_time();
                functionT f = project(world, m);
                rec.t_project = wall_time() - t0;
                rec.norm = f.norm2();
                const std::vector<LeafRow> leaves = measure_leaves(world, f, *m, rec, th);
                if (place == "corner") level_corner = rec.level_at_c;
                f.make_redundant(true);
                const std::vector<InteriorRow> interior = measure_interior(world, f, leaves, rec, th);

                summarize_leaves(rec, leaves);
                summarize_interior(rec, interior);
                std::fflush(stdout);

                if (not ffun) continue;
                std::fprintf(ffun, "%d,%.0e,%s,%s,%.16g,%.16g,%.16g,%ld,%ld,%d,%d,%d,%.6f,%.6f,%.6f,%.12f,%.4f,%.1f,%.1f,%.3f\n",
                             s.k, s.thresh, rec.name.c_str(), place.c_str(), c[0], c[1], c[2], rec.leaves, rec.interior,
                             rec.lmin, rec.lmax, rec.level_at_c, rec.frac_c[0], rec.frac_c[1], rec.frac_c[2], rec.norm,
                             rec.t_project, rec.t_S, rec.t_M, rec.t_map);
                for (const auto& r : leaves)
                    std::fprintf(fleaf, "%d,%.0e,%s,%s,%s,%d,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e\n", s.k, s.thresh,
                                 rec.name.c_str(), place.c_str(), key_cols(r.key).c_str(), int(r.at_c), r.M, r.S, r.T24,
                                 r.Tref, r.A, r.dev, r.norm);
                for (const auto& r : interior)
                    std::fprintf(fint, "%d,%.0e,%s,%s,%s,%d,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e\n", s.k,
                                 s.thresh, rec.name.c_str(), place.c_str(), key_cols(r.key).c_str(), r.deepest, r.Sown,
                                 r.Mown, r.Town, r.supS, r.supM, r.Ttree, r.Atree, r.nt, r.dnt);
            }
        }

        for (const auto& p : pairs) {
            double t0 = wall_time();
            const auto fm = std::make_shared<Model>(p.f, OFFGRID);
            const auto gm = std::make_shared<Model>(p.g, OFFGRID + p.distance * PAIRDIR);
            functionT f = project(world, fm);
            functionT g = project(world, gm);
            const std::vector<UBox> U = measure_union(world, f, g, *fm, *gm, s, th);
            summarize_union(p.name, s, U);
            check_kernel(f, g, U, s, th);
            std::printf("    measured in %.1f s\n", wall_time() - t0);
            std::fflush(stdout);
            if (not fu) continue;
            for (const auto& b : U) {
                std::fprintf(fu, "%d,%.0e,%s,%s,%d,%d,%d", s.k, s.thresh, p.name.c_str(), key_cols(b.key).c_str(),
                             int(b.leaf), int(b.in_f), int(b.in_g));
                for (double x : {b.nf, b.ndf, b.ng, b.ndg, b.supSf, b.supMf, b.supSg, b.supMg, b.Sf, b.Mf, b.Tf, b.Sg,
                                 b.Mg, b.Tg, b.theta, std::sqrt(b.e2rep), std::sqrt(b.e2ana), std::sqrt(b.floor2),
                                 std::sqrt(b.n2fg), std::sqrt(b.c2), std::sqrt(b.k2rep), std::sqrt(b.k2ana), std::sqrt(b.d2q)})
                    std::fprintf(fu, ",%.7e", x);
                for (double x : b.est) std::fprintf(fu, ",%.7e", x);
                for (char x : b.stop) std::fprintf(fu, ",%d", int(x));
                std::fprintf(fu, "\n");
            }
        }
    }
    for (FILE* fp : {ffun, fleaf, fint, fu})
        if (fp) std::fclose(fp);

    std::printf("\ntheorems: polynomial values above S on %ld leaves (worst value/S %.12f), on %ld own polynomials "
                "(worst %.12f); E_rep above the (S, S_own) estimates on %ld boxes (worst E_rep/est %.6f); the kernel "
                "disagrees with the union tree on %ld pairs\n",
                th.leaf, th.worst_leaf, th.own, th.worst_own, th.stop, th.worst_stop, th.kernel);
    const bool ok = (th.leaf == 0 and th.own == 0 and th.stop == 0 and th.kernel == 0);
    std::printf("test_supbound_quantify %s\n", ok ? "passed" : "FAILED");
    world.gop.fence();
    finalize();
    return ok ? 0 : 1;
}
