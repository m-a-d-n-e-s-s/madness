


#ifndef SRC_APPS_CHEM_EXCHANGEOPERATOR_H_
#define SRC_APPS_CHEM_EXCHANGEOPERATOR_H_

#include<madness.h>
#include<madness/world/cloud.h>
#include<madness/mra/macrotaskq.h>

using namespace madness;

namespace madness {

// forward declaration
class SCF;

class Nemo;


template<typename T, std::size_t NDIM>
class Exchange {
    typedef Function<T, NDIM> functionT;
    typedef std::vector<functionT> vecfuncT;

    static inline std::atomic<long> apply_timer;
    static inline std::atomic<long> mul2_timer;
    static inline std::atomic<long> mul1_timer; ///< timing

    static void reset_timer() {
        mul1_timer = 0l;
        mul2_timer = 0l;
        apply_timer = 0l;
    }

    static void print_timer(World& world) {
        double t1 = double(mul1_timer) * 0.001;
        double t2 = double(apply_timer) * 0.001;
        double t3 = double(mul2_timer) * 0.001;
        world.gop.sum(t1);
        world.gop.sum(t2);
        world.gop.sum(t3);
        if (world.rank() == 0) {
            printf(" cpu time spent in multiply1   %8.2fs\n", t1);
            printf(" cpu time spent in apply       %8.2fs\n", t2);
            printf(" cpu time spent in multiply2   %8.2fs\n", t3);
        }
    }

public:

public:
//    class MacroTaskExchange : public MacroTaskIntermediate<MacroTaskExchange> {
//
//        typedef std::vector<std::shared_ptr<MacroTaskBase> > taskqT;
//
//        static inline std::shared_ptr<vecfuncT> mo_ket, mo_bra, vf;
//        static inline std::shared_ptr<real_convolution_3d> poisson;
//
//        std::pair<long,long> row_range, column_range;
//        long nocc=0;
//        long nf=0;
//        double lo=1.e-4;
//        double mul_tol=1.e-7;
//        bool symmetric=true;
//
//
//    public:
//        MacroTaskExchange(const std::pair<long,long>& row_range, const std::pair<long,long>& column_range,
//                          const long nocc, const long nf, const double lo, const double mul_tol, const bool symmetric)
//                : row_range(row_range)
//                , column_range(column_range)
//                , nocc(nocc)
//                , nf(nf)
//                , lo(lo)
//                , mul_tol(mul_tol)
//                , symmetric(symmetric) {
//            this->priority=compute_priority();
//        }
//
//        /// compute the priority of this task for non-dumb scheduling
//
//        /// \return the priority as double number (no limits)
//        double compute_priority() const {
//            long nrow=row_range.second-row_range.first;
//            long ncol=column_range.second-column_range.first;
//            return double(nrow*ncol);
//        }
//
//        /// run a macrotask
//
//        /// input and output is done through the cloud
//        /// \param subworld     the world this task is executed in
//        /// \param cloud        a storage class for input functions
//        /// \param taskq        the taskq (for later reference..)
//        void run(World& subworld, Cloud& cloud, taskqT& taskq);
//
//        /// compute a batch of the exchange matrix, with identical ranges, exploiting the matrix symmetry
//
//        /// \param subworld     the world we're computing in
//        /// \param cloud        where to store the results
//        /// \param bra_batch    the bra batch of orbitals (including the nuclear correlation factor square)
//        /// \param ket_batch    the ket batch of orbitals, i.e. the orbitals to premultiply with
//        /// \param vf_batch     the argument of the exchange operator
//        vecfuncT compute_diagonal_batch_in_symmetric_matrix(World& subworld,
//                                       const vecfuncT& bra_batch, const vecfuncT& ket_batch, const vecfuncT& vf_batch) const {
//            double mul_tol=0.0;
//            double symmetric=true;
//            return Exchange<T,NDIM>::compute_K_tile(subworld,bra_batch,ket_batch,vf_batch,poisson,symmetric,mul_tol);
//        }
//
//        /// compute a batch of the exchange matrix, with non-identical ranges
//
//        /// \param subworld     the world we're computing in
//        /// \param cloud        where to store the results
//        /// \param bra_batch    the bra batch of orbitals (including the nuclear correlation factor square)
//        /// \param ket_batch    the ket batch of orbitals, i.e. the orbitals to premultiply with
//        /// \param vf_batch     the argument of the exchange operator
//        vecfuncT compute_batch_in_asymmetric_matrix(World& subworld,
//                                       const vecfuncT& bra_batch, const vecfuncT& ket_batch, const vecfuncT& vf_batch) const {
//            double mul_tol=0.0;
//            double symmetric=false;
//            return Exchange<T,NDIM>::compute_K_tile(subworld,bra_batch,ket_batch,vf_batch,poisson,symmetric,mul_tol);
//        }
//
//        /// compute a batch of the exchange matrix, with non-identical ranges
//
//        /// \param subworld     the world we're computing in
//        /// \param cloud        where to store the results
//        /// \param bra_batch    the bra batch of orbitals (including the nuclear correlation factor square)
//        /// \param ket_batch    the ket batch of orbitals, i.e. the orbitals to premultiply with
//        /// \param vf_batch     the argument of the exchange operator
//        std::pair<vecfuncT, vecfuncT> compute_offdiagonal_batch_in_symmetric_matrix(World &subworld,
//                                       const vecfuncT& bra_batch, const vecfuncT& ket_batch, const vecfuncT& vf_batch) const;
//
//        void cleanup() {
//            if (mo_ket) mo_ket.reset();
//            if (mo_bra) mo_bra.reset();
//            if (vf) vf.reset();
//            if (poisson) poisson.reset();
//        }
//
//        template <typename Archive>
//        void serialize(const Archive& ar) {
//            ar & row_range & column_range & nocc & lo &  mul_tol ;
//        }
//
//        void print_me(std::string s="") const {
//            print("K apply task", s, this->stat, "priority",this->get_priority());
//        }
//
//    };
public:

    enum Algorithm {
        small_memory, large_memory, multiworld_efficient
    };
    Algorithm algorithm_ = multiworld_efficient;

    /// default ctor
    Exchange(World& world) : world(world), symmetric_(false) {};

    /// ctor with a conventional calculation
    Exchange(World& world, const SCF *calc, const int ispin);

    /// ctor with a nemo calculation
    Exchange(World& world, const Nemo *nemo, const int ispin);

    /// set the bra and ket orbital spaces, and the occupation

    /// @param[in]	bra		bra space, must be provided as complex conjugate
    /// @param[in]	ket		ket space
    void set_parameters(const vecfuncT& bra, const vecfuncT& ket, const double lo1) {
        mo_bra = copy(world, bra);
        mo_ket = copy(world, ket);
        lo = lo1;
    }

    static auto set_poisson(World& world, const double lo, const double econv = FunctionDefaults<3>::get_thresh()) {
        return std::shared_ptr<real_convolution_3d>(CoulombOperatorPtr(world, lo, econv));
    }

    Function<T, NDIM> operator()(const Function<T, NDIM>& ket) const {
        vecfuncT vket(1, ket);
        vecfuncT vKket = this->operator()(vket);
        return vKket[0];
    }

    /// apply the exchange operator on a vector of functions

    /// note that only one spin is used (either alpha or beta orbitals)
    /// @param[in]  vket       the orbitals |i> that the operator is applied on
    /// @return     a vector of orbitals  K| i>
    vecfuncT operator()(const vecfuncT& vket, const double mul_tol = 0.0) const;

    /// compute the matrix element <bra | K | ket>

    /// @param[in]  bra    real_function_3d, the bra state
    /// @param[in]  ket    real_function_3d, the ket state
    T operator()(const Function<T, NDIM>& bra, const Function<T, NDIM>& ket) const {
        return inner(bra, this->operator()(ket));
    }

    /// compute the matrix < vbra | K | vket >

    /// @param[in]  vbra    vector of real_function_3d, the set of bra states
    /// @param[in]  vket    vector of real_function_3d, the set of ket states
    /// @return K_ij
    Tensor<T> operator()(const vecfuncT& vbra, const vecfuncT& vket) const {
        const auto bra_equiv_ket = &vbra == &vket;
        vecfuncT vKket = this->operator()(vket);
        auto n = norm2s(world, vKket);
        auto result = matrix_inner(world, vbra, vKket, bra_equiv_ket);
        return result;
    }

    bool is_symmetric() const { return symmetric_; }

    Exchange& symmetric(const bool flag) {
        symmetric_ = flag;
        return *this;
    }

    Exchange& set_algorithm(const Algorithm& alg) {
        algorithm_ = alg;
        return *this;
    }

private:

    /// exchange using macrotasks, i.e. apply K on a function in individual worlds
    vecfuncT K_macrotask_efficient(const vecfuncT& vket, const double mul_tol = 0.0) const;

    /// computing the full square of the double sum (over vket and the K orbitals)
    vecfuncT K_small_memory(const vecfuncT& vket, const double mul_tol = 0.0) const;

    /// computing the upper triangle of the double sum (over vket and the K orbitals)
    vecfuncT K_large_memory(const vecfuncT& vket, const double mul_tol = 0.0) const;

    /// computing the upper triangle of the double sum (over vket and the K orbitals)
    static vecfuncT compute_K_tile(World& world, const vecfuncT& mo_bra, const vecfuncT& mo_ket,
                                   const vecfuncT& vket, std::shared_ptr<real_convolution_3d> poisson,
                                   const bool symmetric, const double mul_tol = 0.0);

    inline bool do_print_timings() const { return (world.rank() == 0) and (printlevel >= 3); }

    inline bool printdebug() const { return (world.rank() == 0) and (printlevel >= 10); }

    World& world;
    bool symmetric_ = false;      /// is the exchange matrix symmetric? K phi_i = \sum_k \phi_k \int \phi_k \phi_i
    vecfuncT mo_bra, mo_ket;    ///< MOs for bra and ket
    double lo = 1.e-4;
    long printlevel = 0;

    class MacroTaskExchangeSimple : public MicroTaskBase {

        long nresult;
        double lo = 1.e-4;
        double mul_tol = 1.e-7;
        bool symmetric = true;

        /// custom partitioning for the exchange operator in exchangeoperator.h

        /// arguments are: result[i] += sum_k vket[k] \int 1/r vbra[k] f[i]
        /// with f and vbra being batched, result and vket being passed on as a whole
        class MacroTaskPartitionerExchange : public MacroTaskPartitioner {
        public:
            MacroTaskPartitionerExchange(const bool symmetric) : symmetric(symmetric) {}
            bool symmetric=false;
            partitionT do_partitioning(const std::size_t &vsize1, const std::size_t &vsize2,
                                       const std::string policy) const override {

                partitionT partition1=do_1d_partition(vsize1,policy);
                partitionT partition2=do_1d_partition(vsize2,policy);
                partitionT result;
                for (auto i=partition1.begin(); i!=partition1.end(); ++i) {
                    if (symmetric) {
                        for (auto j=i; j!=partition1.end(); ++j) {
                            result.push_back(Batch(i->input[0],j->input[0],_));
                        }
                    } else {
                        for (auto j=partition2.begin(); j!=partition2.end(); ++j) {
                            result.push_back(Batch(i->input[0],j->input[0],_));
                        }
                    }
                }
                return result;
            }
        };

    public:
        MacroTaskExchangeSimple(const long nresult, const double lo, const double mul_tol, const bool symmetric)
                : nresult(nresult), lo(lo), mul_tol(mul_tol), symmetric(symmetric) {
            partitioner.reset(new MacroTaskPartitionerExchange(symmetric));
        }

        /// compute the priority of this task for non-dumb scheduling

        /// \return the priority as double number (no limits)
        double compute_priority() const {
            MADNESS_CHECK(batch.input.size() == 2);   // must be quadratic batches
            long nrow = batch.input[0].size();
            long ncol = batch.input[1].size();
            return double(nrow * ncol);
        }

        // you need to define the exact argument(s) of operator() as tuple
        typedef std::tuple<const std::vector<Function<T,NDIM>>&,
                           const std::vector<Function<T,NDIM>>&,
                           const std::vector<Function<T,NDIM>>&> argtupleT;

        using resultT = std::vector<Function<T,NDIM>>;

        // you need to define an empty constructor for the result
        // resultT must implement operator+=(const resultT&)
        resultT allocator(World& world, const argtupleT& argtuple) const {
            std::size_t n = std::get<0>(argtuple).size();
            resultT result = zero_functions_compressed<T,NDIM>(world, n);
            return result;
        }

        std::vector<Function<T,NDIM>> operator()(const std::vector<Function<T,NDIM>>& vf_batch,     // will be batched (column)
                                                 const std::vector<Function<T,NDIM>>& bra_batch,    // will be batched (row)
                                                 const std::vector<Function<T,NDIM>>& vket) {       // will not be batched

            World& world = vf_batch.front().world();
            resultT Kf = zero_functions_compressed<T,NDIM>(world, nresult);

            bool diagonal_block = batch.input[0] == batch.input[1];
            auto& bra_range = batch.input[1];    // corresponds to vbra
            auto& vf_range = batch.input[0];       // corresponds to vf_batch

            if (vf_range.is_full_size()) vf_range.end=vf_batch.size();
            if (bra_range.is_full_size()) bra_range.end=bra_batch.size();

            MADNESS_CHECK(vf_range.end <= nresult);
            if (symmetric) MADNESS_CHECK(bra_range.end <= nresult);

            if (symmetric and diagonal_block) {
                auto ket_batch=bra_range.copy_batch(vket);
                vecfuncT resultcolumn = compute_diagonal_batch_in_symmetric_matrix(world, ket_batch, bra_batch, vf_batch);

                for (int i = vf_range.begin; i < vf_range.end; ++i)
                    Kf[i] += resultcolumn[i - vf_range.begin];

            } else if (symmetric and not diagonal_block) {
                auto[resultcolumn, resultrow]=compute_offdiagonal_batch_in_symmetric_matrix(world, vket, bra_batch, vf_batch);

                for (int i = bra_range.begin; i < bra_range.end; ++i)
                    Kf[i] += resultcolumn[i - bra_range.begin];
                for (int i = vf_range.begin; i < vf_range.end; ++i)
                    Kf[i] += resultrow[i - vf_range.begin];
            } else {
                auto ket_batch=vf_range.copy_batch(vket);
                vecfuncT resultcolumn = compute_batch_in_asymmetric_matrix(world, ket_batch, bra_batch, vf_batch);
                for (int i = bra_range.begin; i < bra_range.end; ++i)
                    Kf[i] += resultcolumn[i - bra_range.begin];
            }
            return Kf;
        }

        /// compute a batch of the exchange matrix, with identical ranges, exploiting the matrix symmetry

        /// \param subworld     the world we're computing in
        /// \param cloud        where to store the results
        /// \param bra_batch    the bra batch of orbitals (including the nuclear correlation factor square)
        /// \param ket_batch    the ket batch of orbitals, i.e. the orbitals to premultiply with
        /// \param vf_batch     the argument of the exchange operator
        vecfuncT compute_diagonal_batch_in_symmetric_matrix(World& subworld,
                                                            const vecfuncT& ket_batch,      // is batched
                                                            const vecfuncT& bra_batch,      // is batched
                                                            const vecfuncT& vf_batch        // is batched
                                                            ) const {
            double mul_tol = 0.0;
            double symmetric = true;
            auto poisson = Exchange<double, 3>::set_poisson(subworld, lo);
            return Exchange<T, NDIM>::compute_K_tile(subworld, bra_batch, ket_batch, vf_batch, poisson, symmetric,
                                                     mul_tol);
        }

        /// compute a batch of the exchange matrix, with non-identical ranges

        /// \param subworld     the world we're computing in
        /// \param cloud        where to store the results
        /// \param bra_batch    the bra batch of orbitals (including the nuclear correlation factor square)
        /// \param ket_batch    the ket batch of orbitals, i.e. the orbitals to premultiply with
        /// \param vf_batch     the argument of the exchange operator
        vecfuncT compute_batch_in_asymmetric_matrix(World& subworld,
                                                    const vecfuncT& ket_batch,
                                                    const vecfuncT& bra_batch,
                                                    const vecfuncT& vf_batch) const {
            double mul_tol = 0.0;
            double symmetric = false;
            auto poisson = Exchange<double, 3>::set_poisson(subworld, lo);
            return Exchange<T, NDIM>::compute_K_tile(subworld, bra_batch, ket_batch, vf_batch, poisson, symmetric,
                                                     mul_tol);
        }

        /// compute a batch of the exchange matrix, with non-identical ranges

        /// \param subworld     the world we're computing in
        /// \param cloud        where to store the results
        /// \param bra_batch    the bra batch of orbitals (including the nuclear correlation factor square)
        /// \param ket_batch    the ket batch of orbitals, i.e. the orbitals to premultiply with
        /// \param vf_batch     the argument of the exchange operator
        std::pair<vecfuncT, vecfuncT> compute_offdiagonal_batch_in_symmetric_matrix(World& subworld,
                                                                                    const vecfuncT& ket, // not batched
                                                                                    const vecfuncT& bra_batch, // batched
                                                                                    const vecfuncT& vf_batch) const; // batched

    };

};

} /* namespace madness */

#endif /* SRC_APPS_CHEM_EXCHANGEOPERATOR_H_ */
