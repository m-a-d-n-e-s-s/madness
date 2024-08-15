#ifndef SRC_APPS_CHEM_EXCHANGEOPERATOR_H_
#define SRC_APPS_CHEM_EXCHANGEOPERATOR_H_

#include<madness.h>
#include<madness/world/cloud.h>
#include<madness/mra/macrotaskq.h>
#include<madness/chem/SCFOperators.h>

namespace madness {

// forward declaration
class SCF;
class Nemo;


template<typename T, std::size_t NDIM>
class Exchange<T,NDIM>::ExchangeImpl {
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

    typedef Exchange<T,NDIM>::Algorithm Algorithm;
    Algorithm algorithm_ = multiworld_efficient;

    /// default ctor
    ExchangeImpl(World& world, const double lo, const double thresh) : world(world), lo(lo), thresh(thresh) {}

    /// ctor with a conventional calculation
    ExchangeImpl(World& world, const SCF *calc, const int ispin) ;

    /// ctor with a nemo calculation
    ExchangeImpl(World& world, const Nemo *nemo, const int ispin);

    /// set the bra and ket orbital spaces, and the occupation

    /// @param[in]	bra		bra space, must be provided as complex conjugate
    /// @param[in]	ket		ket space
    void set_bra_and_ket(const vecfuncT& bra, const vecfuncT& ket) {
        mo_bra = copy(world, bra);
        mo_ket = copy(world, ket);
    }

    std::string info() const {return "K";}

    static auto set_poisson(World& world, const double lo, const double econv = FunctionDefaults<3>::get_thresh()) {
        return std::shared_ptr<real_convolution_3d>(CoulombOperatorPtr(world, lo, econv));
    }

    /// apply the exchange operator on a vector of functions

    /// note that only one spin is used (either alpha or beta orbitals)
    /// @param[in]  vket       the orbitals |i> that the operator is applied on
    /// @return     a vector of orbitals  K| i>
    vecfuncT operator()(const vecfuncT& vket) const;

    bool is_symmetric() const { return symmetric_; }

    ExchangeImpl& set_taskq(std::shared_ptr<MacroTaskQ> taskq1) {
        this->taskq=taskq1;
        return *this;
    }

    ExchangeImpl& symmetric(const bool flag) {
        symmetric_ = flag;
        return *this;
    }

    ExchangeImpl& set_algorithm(const Algorithm& alg) {
        algorithm_ = alg;
        return *this;
    }

    ExchangeImpl& set_printlevel(const long& level) {
        printlevel=level;
        return *this;
    }

private:

    /// exchange using macrotasks, i.e. apply K on a function in individual worlds via tiles
    vecfuncT K_macrotask_efficient(const vecfuncT& vket, const double mul_tol = 0.0) const;

    /// exchange using macrotasks, i.e. apply K on a function in individual worlds row-wise
    vecfuncT K_macrotask_efficient_row(const vecfuncT& vket, const double mul_tol = 0.0) const;

    /// computing the full square of the double sum (over vket and the K orbitals)
    vecfuncT K_small_memory(const vecfuncT& vket, const double mul_tol = 0.0) const;

    /// computing the upper triangle of the double sum (over vket and the K orbitals)
    vecfuncT K_large_memory(const vecfuncT& vket, const double mul_tol = 0.0) const;

    /// computing the upper triangle of the double sum (over vket and the K orbitals)
    static vecfuncT compute_K_tile(World& world, const vecfuncT& mo_bra, const vecfuncT& mo_ket,
                                   const vecfuncT& vket, std::shared_ptr<real_convolution_3d> poisson,
                                   const bool symmetric, const double mul_tol = 0.0);

    inline bool do_print_timings() const { return (world.rank() == 0) and (printlevel >= 3); }

    inline bool printdebug() const {return printlevel >= 10; }

    World& world;
    std::shared_ptr<MacroTaskQ> taskq;
    bool symmetric_ = false;      /// is the exchange matrix symmetric? K phi_i = \sum_k \phi_k \int \phi_k \phi_i
    vecfuncT mo_bra, mo_ket;    ///< MOs for bra and ket
    double lo = 1.e-4;
    double thresh = FunctionDefaults<NDIM>::get_thresh();
    long printlevel = 0;
    double mul_tol = FunctionDefaults<NDIM>::get_thresh()*0.1;

    class MacroTaskExchangeSimple : public MacroTaskOperationBase {

        long nresult;
        double lo = 1.e-4;
        double mul_tol = 1.e-7;
        bool symmetric = false;

        /// custom partitioning for the exchange operator in exchangeoperator.h

        /// arguments are: result[i] += sum_k vket[k] \int 1/r vbra[k] f[i]
        /// with f and vbra being batched, result and vket being passed on as a whole
        class MacroTaskPartitionerExchange : public MacroTaskPartitioner {
        public:
            MacroTaskPartitionerExchange(const bool symmetric) : symmetric(symmetric) {
                max_batch_size=30;
            }

            bool symmetric = false;

            partitionT do_partitioning(const std::size_t& vsize1, const std::size_t& vsize2,
                                       const std::string policy) const override {

                partitionT partition1 = do_1d_partition(vsize1, policy);
                partitionT partition2 = do_1d_partition(vsize2, policy);
                partitionT result;
                for (auto i = partition1.begin(); i != partition1.end(); ++i) {
                    if (symmetric) {
                        for (auto j = i; j != partition1.end(); ++j) {
                            Batch batch(i->first.input[0], j->first.input[0], _);
                            double priority=compute_priority(batch);
                            result.push_back(std::make_pair(batch,priority));
                        }
                    } else {
                        for (auto j = partition2.begin(); j != partition2.end(); ++j) {
                            Batch batch(i->first.input[0], j->first.input[0], _);
                            double priority=compute_priority(batch);
                            result.push_back(std::make_pair(batch,priority));
                        }
                    }
                }
                return result;
            }

            /// compute the priority of this task for non-dumb scheduling

            /// \return the priority as double number (no limits)
            double compute_priority(const Batch& batch) const override {
                MADNESS_CHECK(batch.input.size() == 2);   // must be quadratic batches
                long nrow = batch.input[0].size();
                long ncol = batch.input[1].size();
                return double(nrow * ncol);
            }
        };

    public:
        MacroTaskExchangeSimple(const long nresult, const double lo, const double mul_tol, const bool symmetric)
                : nresult(nresult), lo(lo), mul_tol(mul_tol), symmetric(symmetric) {
            partitioner.reset(new MacroTaskPartitionerExchange(symmetric));
        }


        // you need to define the exact argument(s) of operator() as tuple
        typedef std::tuple<const std::vector<Function<T, NDIM>>&,
                const std::vector<Function<T, NDIM>>&,
                const std::vector<Function<T, NDIM>>&> argtupleT;

        using resultT = std::vector<Function<T, NDIM>>;

        // you need to define an empty constructor for the result
        // resultT must implement operator+=(const resultT&)
        resultT allocator(World& world, const argtupleT& argtuple) const {
            std::size_t n = std::get<0>(argtuple).size();
            resultT result = zero_functions_compressed<T, NDIM>(world, n);
            return result;
        }

        std::vector<Function<T, NDIM>>
        operator()(const std::vector<Function<T, NDIM>>& vf_batch,     // will be batched (column)
                   const std::vector<Function<T, NDIM>>& bra_batch,    // will be batched (row)
                   const std::vector<Function<T, NDIM>>& vket) {       // will not be batched

            World& world = vf_batch.front().world();
            resultT Kf = zero_functions_compressed<T, NDIM>(world, nresult);

            bool diagonal_block = batch.input[0] == batch.input[1];
            auto& bra_range = batch.input[1];    // corresponds to vbra
            auto& vf_range = batch.input[0];       // corresponds to vf_batch

            if (vf_range.is_full_size()) vf_range.end = vf_batch.size();
            if (bra_range.is_full_size()) bra_range.end = bra_batch.size();

            MADNESS_CHECK(vf_range.end <= nresult);
            if (symmetric) MADNESS_CHECK(bra_range.end <= nresult);

            if (symmetric and diagonal_block) {
                auto ket_batch = bra_range.copy_batch(vket);
                vecfuncT resultcolumn = compute_diagonal_batch_in_symmetric_matrix(world, ket_batch, bra_batch,
                                                                                   vf_batch);

                for (int i = vf_range.begin; i < vf_range.end; ++i){
                    Kf[i] += resultcolumn[i - vf_range.begin];}

            } else if (symmetric and not diagonal_block) {
                auto[resultcolumn, resultrow]=compute_offdiagonal_batch_in_symmetric_matrix(world, vket, bra_batch,
                                                                                            vf_batch);

                for (int i = bra_range.begin; i < bra_range.end; ++i){
                    Kf[i] += resultcolumn[i - bra_range.begin];}
                for (int i = vf_range.begin; i < vf_range.end; ++i){
                    Kf[i] += resultrow[i - vf_range.begin];}
            } else {
                auto ket_batch = bra_range.copy_batch(vket);
                vecfuncT resultcolumn = compute_batch_in_asymmetric_matrix(world, ket_batch, bra_batch, vf_batch);
                for (int i = vf_range.begin; i < vf_range.end; ++i)
                    Kf[i] += resultcolumn[i - vf_range.begin];
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
            auto poisson = Exchange<double, 3>::ExchangeImpl::set_poisson(subworld, lo);
            return Exchange<T, NDIM>::ExchangeImpl::compute_K_tile(subworld, bra_batch, ket_batch, vf_batch, poisson, symmetric,
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
            auto poisson = Exchange<double, 3>::ExchangeImpl::set_poisson(subworld, lo);
            return Exchange<T, NDIM>::ExchangeImpl::compute_K_tile(subworld, bra_batch, ket_batch, vf_batch, poisson, symmetric,
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

    class MacroTaskExchangeRow : public MacroTaskOperationBase {

        long nresult;
        double lo = 1.e-4;
        double mul_tol = 1.e-7;
        bool symmetric = false;

        /// custom partitioning for the exchange operator in exchangeoperator.h
        class MacroTaskPartitionerRow : public MacroTaskPartitioner {
        public:
            MacroTaskPartitionerRow() {
              max_batch_size=1;
            }       
        };

    public:
        MacroTaskExchangeRow(const long nresult, const double lo, const double mul_tol)
                : nresult(nresult), lo(lo), mul_tol(mul_tol) {
            partitioner.reset(new MacroTaskPartitionerRow());
        }

        // you need to define the exact argument(s) of operator() as tuple
        typedef std::tuple<const std::vector<Function<T, NDIM>>&,
                           const std::vector<Function<T, NDIM>>&,
                           const std::vector<Function<T, NDIM>>&> argtupleT;

        using resultT = std::vector<Function<T, NDIM>>;

        // you need to define an empty constructor for the result
        // resultT must implement operator+=(const resultT&)
        resultT allocator(World& world, const argtupleT& argtuple) const {
            std::size_t n = std::get<0>(argtuple).size();
            resultT result = zero_functions_compressed<T, NDIM>(world, n);
            return result;
        }

        /// compute exchange row-wise for a fixed orbital phi_i of vket
        std::vector<Function<T, NDIM>>
        operator()(const std::vector<Function<T, NDIM>>& vket,
                   const std::vector<Function<T, NDIM>>& mo_bra, 
                   const std::vector<Function<T, NDIM>>& mo_ket) {       

            World& world = vket.front().world();
            mul_tol = 0.0;
            print("mul_tol ", mul_tol);
            

            resultT Kf = zero_functions_compressed<T, NDIM>(world, 1);
            auto poisson = Exchange<double, 3>::ExchangeImpl::set_poisson(world, lo);

            auto& i = batch.input[0].begin;  

            size_t min_tile = 10;
            size_t ntile = std::min(mo_bra.size(), min_tile);
            vecfuncT psif = zero_functions_compressed<T,NDIM>(world, mo_bra.size()); 

            // compute psif
            for (size_t ilo=0; ilo<mo_bra.size(); ilo+=ntile){
                size_t iend = std::min(ilo+ntile,mo_bra.size());
                vecfuncT tmp_mo_bra(mo_bra.begin()+ilo,mo_bra.begin()+iend);
                auto tmp_psif = mul_sparse(world, vket[i], tmp_mo_bra, mul_tol);
                print_size(world, tmp_psif, "tmp_psif before truncation");

                //truncate tmp_mo_bra
                truncate(world, tmp_psif);
                print_size(world, tmp_psif, "tmp_psi_f after truncation");

                //put the results into their final home
                for (size_t i = ilo; i<iend; ++i){
                    psif[i] += tmp_psif[i-ilo];
                }

            }

            //vecfuncT psif = mul_sparse(world, vket[i], mo_bra, mul_tol);
            //auto size=get_size(world,psif);
            //if (world.rank()==0 && printdebug()) print("size of psif after mul_sparse",size);

            //truncate(world, psif);

            //size=get_size(world,psif);
            //if (world.rank()==0 && printdebug()) print("size of psif after truncating",size);

            psif = apply(world, *poisson.get(), psif);
            //size=get_size(world,psif);
            //if (world.rank()==0 && printdebug()) print("size of psif after apply",size);
            truncate(world, psif);
            //size=get_size(world,psif);
            //if (world.rank()==0 && printdebug()) print("size of psif after truncating",size);
            
            // TODO: priority = #coeffs
            for (size_t ilo=0; ilo<mo_ket.size(); ilo+=ntile){
                size_t iend = std::min(ilo+ntile,mo_ket.size());
                vecfuncT tmp_mo_ket(mo_ket.begin()+ilo,mo_ket.begin()+iend);
                vecfuncT tmp_psif(psif.begin()+ilo,psif.begin()+iend);
                // TODO: use matrix_mul_sparse instead, need to implement mul_sparse for
                //       vecfuncT, vecfuncT
                auto tmp_Kf = dot(world, tmp_mo_ket, tmp_psif);
                //auto tmp_Kf = mul_sparse(world, tmp_mo_ket, tmp_psif, mul_tol);

                //print_size(world, tmp_Kf, "tmp_Kf before truncation");
                //tmp_Kf.print_size("tmp_Kf before truncation");
                //truncate
                //truncate(world, tmp_Kf);
                //print_size(world, tmp_Kf, "tmp_Kf after truncation");

                //put the results into their final home
                //for (size_t i = ilo; i<iend; ++i){
                //    Kf[i] += tmp_Kf[i-ilo];
                //}
                Kf[0] += tmp_Kf;
                print_size(world, Kf, "Kf before truncation");
                truncate(world, Kf);
                print_size(world, Kf, "Kf after truncation");
            }
            //Kf[0] +=dot(world, mo_ket, psif);
            //if (world.rank()==0 && printdebug()) Kf[0].print_size("Kf after dot");
            
            //truncate(world, Kf);
            //if (world.rank()==0 && printdebug()) Kf[i].print_size("Kf after truncating");

            return Kf;
        }
    };
};

} /* namespace madness */

#endif /* SRC_APPS_CHEM_EXCHANGEOPERATOR_H_ */
