


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
	typedef Function<T,NDIM> functionT;
	typedef std::vector<functionT> vecfuncT;

    static inline std::atomic<long> apply_timer;
    static inline std::atomic<long> mul2_timer;
    static inline std::atomic<long> mul1_timer; ///< timing

    static void reset_timer() {
        mul1_timer=0l;
        mul2_timer=0l;
        apply_timer=0l;
    }

    static void print_timer(World& world) {
        double t1=double(mul1_timer)*0.001;
        double t2=double(apply_timer)*0.001;
        double t3=double(mul2_timer)*0.001;
        world.gop.sum(t1);
        world.gop.sum(t2);
        world.gop.sum(t3);
        if (world.rank()==0) {
            printf(" cpu time spent in multiply1   %8.2fs\n",t1);
            printf(" cpu time spent in apply       %8.2fs\n",t2);
            printf(" cpu time spent in multiply2   %8.2fs\n",t3);
        }
    }

public:

    class MacroTaskExchangeEfficient : public MacroTaskIntermediate<MacroTaskExchangeEfficient> {

        typedef std::vector<std::shared_ptr<MacroTaskBase> > taskqT;

        static inline std::shared_ptr<vecfuncT> mo_ket, mo_bra;
        static inline std::shared_ptr<real_convolution_3d> poisson;

        std::pair<long,long> row_range, column_range;
        long nocc=0;
        double lo=1.e-4;
        double econv=1.e-6;
        double mul_tol=1.e-7;
        long nrange=1;

        struct accumulate_into_result_op {
            Function<T,NDIM> result;
            accumulate_into_result_op(Function<T,NDIM>& r) {
                result.set_impl(r.get_impl());
            }

            accumulate_into_result_op(Cloud& cloud, World& subworld, const int record) {
                std::shared_ptr<FunctionImpl<T,NDIM> > rimpl;
                cloud.load(subworld,rimpl,record);
                result.set_impl(rimpl);
            }

            void operator()(const Key<NDIM>& key, FunctionNode<T,NDIM>& node) const {
                auto coeffs = const_cast<typename FunctionImpl<T, NDIM>::dcT &>(result.get_impl()->get_coeffs());
                coeffs.send(key, &FunctionNode<T,NDIM>:: template gaxpy_inplace<T,T>, 1.0, node, 1.0);
            }
        };


    public:
        static inline bool async_accumulation=true;
        MacroTaskExchangeEfficient(const std::pair<long,long>& row_range, const std::pair<long,long>& column_range,
                          const long nocc, const double lo, const double econv, const double mul_tol, const long nrange)
                : row_range(row_range)
                , column_range(column_range)
                , nocc(nocc)
                , lo(lo)
                , econv(econv)
                , mul_tol(mul_tol)
                , nrange(nrange){
            this->priority=compute_priority();
        }

        double compute_priority() const {
            long nrow=row_range.second-row_range.first;
            long ncol=column_range.second-column_range.first;
            return double(nrow*ncol);
        }


        void run(World& subworld, Cloud& cloud, taskqT& taskq) {

            if (not poisson) poisson = Exchange<T,NDIM>::set_poisson(subworld,lo,econv);
            // the argument of the exchange operator is the ket vector

            // load bra and ket if not already loaded
            if (not mo_bra.get()) {
                mo_bra.reset(new vecfuncT(nocc));
                for (int i = 0; i < nocc; ++i) (*mo_bra)[i] = cloud.load<functionT>(subworld, i);
            }
            if (not mo_ket.get()) {
                mo_ket.reset(new vecfuncT(nocc));
                for (int i = 0; i < nocc; ++i) (*mo_ket)[i] = cloud.load<functionT>(subworld, i + nocc);
            }

            // make universe-living Kf accessible here in the subworld for result accumulation
            vecfuncT Kf(nocc);
            for (int i=0; i<nocc; ++i) {
                std::shared_ptr<FunctionImpl<T,NDIM> > rimpl;
                cloud.load(subworld,rimpl,i+3*nocc);
                Kf[i].set_impl(rimpl);
            }

            // compute the tile [column_range,row_range], corresponding to bra[nrow], ket[ncolumn]
            vecfuncT bra_batch(mo_bra->begin() + row_range.first, mo_bra->begin() + row_range.second);
            vecfuncT ket_batch(mo_ket->begin() + column_range.first, mo_ket->begin() + column_range.second);

            const double truncate_tol=FunctionDefaults<NDIM>::get_thresh()/nrange*0.01;

            if (row_range == column_range) {
                vecfuncT resultcolumn=compute_symmetric_batch(subworld, bra_batch, ket_batch, truncate_tol);

                if (async_accumulation) {
                    for (int i = column_range.first; i < column_range.second; ++i) {
                        resultcolumn[i - column_range.first].unaryop_node(accumulate_into_result_op(Kf[i]));
                    }
                } else {
                    // store results: columns as columns
                    cloud.store(subworld, resultcolumn, outputrecord(row_range, column_range));
                }

            } else {
                auto [resultcolumn,resultrow]=compute_batch(subworld, bra_batch, ket_batch, truncate_tol);

                if (async_accumulation) {
                    for (int i = column_range.first; i < column_range.second; ++i)
                        resultrow[i - column_range.first].unaryop_node(accumulate_into_result_op(Kf[i]));
                    for (int i = row_range.first; i < row_range.second; ++i)
                        resultcolumn[i - row_range.first].unaryop_node(accumulate_into_result_op(Kf[i]));
                } else {
                    // store results: columns as columns; transpose rows to columns
                    cloud.store(subworld, resultcolumn, outputrecord(row_range, column_range));
                    cloud.store(subworld, resultrow, outputrecord(column_range, row_range));
                }

            }
            subworld.gop.fence();
        }

        /// compute a batch of the exchange matrix, with identical ranges, exploiting the matrix symmetry

        /// \param subworld     the world we're computing in
        /// \param cloud        where to store the results
        /// \param bra_batch    the bra batch of orbitals (including the nuclear correlation factor square)
        /// \param ket_batch    the ket batch of orbitals, also the orbitals to premultiply with
        vecfuncT compute_symmetric_batch(World& subworld, const vecfuncT& bra_batch,
                                         const vecfuncT& ket_batch, const double truncate_tol) const {
            Tensor<double> occ(ket_batch.size());
            occ=1.0;
            double mul_tol=0.0;
            double same=true;
            return Exchange<T,NDIM>::compute_K_tile(subworld,bra_batch,ket_batch,ket_batch,poisson,same,occ,mul_tol,truncate_tol);
        }

        /// compute a batch of the exchange matrix, with non-identical ranges

        /// \param subworld     the world we're computing in
        /// \param cloud        where to store the results
        /// \param bra_batch    the bra batch of orbitals (including the nuclear correlation factor square)
        /// \param ket_batch    the ket batch of orbitals, also the orbitals to premultiply with
        std::pair<vecfuncT, vecfuncT> compute_batch(World& subworld, const vecfuncT& bra_batch,
                                                    const vecfuncT& ket_batch, const double truncate_tol) const {
            // orbital_product is a vector of vectors


            double cpu0 = cpu_time();
            std::vector<vecfuncT> orbital_product=matrix_mul_sparse<T,T,NDIM>(subworld,bra_batch,ket_batch,mul_tol);
            vecfuncT orbital_product_flat=flatten(orbital_product); // convert into a flattened vector
            truncate(subworld, orbital_product_flat);
            double cpu1=cpu_time();
            mul1_timer+=long((cpu1-cpu0)*1000l);

            cpu0=cpu_time();
            vecfuncT Nij = apply(subworld, *poisson.get(), orbital_product_flat);
            truncate(subworld, Nij);
            cpu1=cpu_time();
            apply_timer+=long((cpu1-cpu0)*1000l);

            // accumulate columns:      resultrow(i)=\sum_j j N_ij
            // accumulate rows:      resultcolumn(j)=\sum_i i N_ij
            cpu0=cpu_time();

            // some helper functions
            std::size_t nrow=bra_batch.size();
            std::size_t ncolumn=ket_batch.size();
            auto ij = [&nrow, &ncolumn](const int i, const int j) {return i*ncolumn+j;};

            auto Nslice = [&Nij, &ij, &ncolumn] (const long irow, const Slice s) {
                vecfuncT result;
                MADNESS_CHECK(s.start==0 && s.end==-1 && s.step==1);
                for (std::size_t i=s.start; i<=s.end+ncolumn; ++i) {
                    result.push_back(Nij[ij(irow,i)]);
                }
                return result;
            };
            auto Nslice1 = [&Nij, &ij, &nrow] (const Slice s, const long jcolumn) {
                vecfuncT result;
                MADNESS_CHECK(s.start==0 && s.end==-1 && s.step==1);
                for (int i=s.start; i<=s.end+nrow; ++i) {
                    result.push_back(Nij[ij(i,jcolumn)]);
                }
                return result;
            };

            // corresponds to bra_batch and ket_batch, but without the ncf R^2
            vecfuncT preintegral_row(mo_ket->begin()+row_range.first,mo_ket->begin()+row_range.second);
            vecfuncT preintegral_column(mo_ket->begin()+column_range.first,mo_ket->begin()+column_range.second);

            vecfuncT resultcolumn(nrow);
            for (std::size_t irow=0; irow<nrow; ++irow) {
                resultcolumn[irow]=dot(subworld,preintegral_column,Nslice(irow,_));  // sum over columns result=sum_j ket[j] N[i,j]
            }
            vecfuncT resultrow(ncolumn);
            for (std::size_t icolumn=0; icolumn<ncolumn; ++icolumn) {
                resultrow[icolumn]=dot(subworld,preintegral_row,Nslice1(_,icolumn));  // sum over rows result=sum_i ket[i] N[i,j]
            }
            cpu1=cpu_time();
            mul2_timer+=long((cpu1-cpu0)*1000l);

            truncate(subworld,resultcolumn,truncate_tol,false);
            truncate(subworld,resultrow,truncate_tol,false);
            subworld.gop.fence();

            return std::make_pair(resultcolumn,resultrow);
        }

        void cleanup() {
            if (mo_ket) {
                print("clearing mo_ket");
                mo_ket.reset();
            }
            if (mo_bra) mo_bra.reset();
        }

        /// return the record to store the output based on the batch

        /// \param range        the range of the batch
        /// \param dimension    0: row, 1: column
        /// \return             the record to store the result (positive for column, negative for row)
        static long outputrecord(const std::pair<long,long>& rowrange, const std::pair<long,long>& columnrange) {
            // Cantor pairing
            auto cantor = [] (const long k1, const long k2) {
                return ((k1 + k2) * (k1 + k2 + 1)/2 + k2);
            };
            return cantor(rowrange.first,columnrange.first);
        }

        template <typename Archive>
        void serialize(const Archive& ar) {
            ar & row_range & column_range & nocc & lo & econv & mul_tol & nrange;
        }

        void print_me(std::string s="") const {
            print("K apply task", s, this->stat, "priority",this->get_priority());
        }

    };
public:

    enum Algorithm {small_memory, large_memory, multiworld_efficient};
    Algorithm algorithm_=multiworld_efficient;

    /// default ctor
    Exchange(World& world) : world(world), same_(false) {};

    /// ctor with a conventional calculation
    Exchange(World& world, const SCF* calc, const int ispin);

    /// ctor with a nemo calculation
    Exchange(World& world, const Nemo* nemo, const int ispin);

    /// set the bra and ket orbital spaces, and the occupation

    /// @param[in]	bra		bra space, must be provided as complex conjugate
    /// @param[in]	ket		ket space
    /// @param[in]	occ1	occupation numbers
    void set_parameters(const vecfuncT& bra, const vecfuncT& ket,
            const Tensor<double>& occ1, const double lo1,
            const double econv1=FunctionDefaults<3>::get_thresh()) {
    	mo_bra=copy(world,bra);
    	mo_ket=copy(world,ket);
    	occ=copy(occ1);
    	econv=econv1;
        lo=lo1;
    }

    static auto set_poisson(World& world, const double lo, const double econv) {
        print("setting the poisson operator with lo, econv", lo, econv);
        return std::shared_ptr<real_convolution_3d>(
                CoulombOperatorPtr(world, std::min(lo,1.e-4), econv));
    }

    Function<T,NDIM> operator()(const Function<T,NDIM>& ket) const {
        vecfuncT vket(1,ket);
        vecfuncT vKket=this->operator()(vket);
        return vKket[0];
    }

    /// apply the exchange operator on a vector of functions

    /// note that only one spin is used (either alpha or beta orbitals)
    /// @param[in]  vket       the orbitals |i> that the operator is applied on
    /// @return     a vector of orbitals  K| i>
    vecfuncT operator()(const vecfuncT& vket,const double mul_tol=0.0) const;

    /// compute the matrix element <bra | K | ket>

    /// @param[in]  bra    real_function_3d, the bra state
    /// @param[in]  ket    real_function_3d, the ket state
    T operator()(const Function<T,NDIM>& bra, const Function<T,NDIM>& ket) const {
        return inner(bra,this->operator()(ket));
    }

    /// compute the matrix < vbra | K | vket >

    /// @param[in]  vbra    vector of real_function_3d, the set of bra states
    /// @param[in]  vket    vector of real_function_3d, the set of ket states
    /// @return K_ij
    Tensor<T> operator()(const vecfuncT& vbra, const vecfuncT& vket) const {
        const auto bra_equiv_ket = &vbra == &vket;
        vecfuncT vKket=this->operator()(vket);
        auto result=matrix_inner(world,vbra,vKket,bra_equiv_ket);
        if (world.rank()==0) {
        	print("result matrix in K");
        	print(result);
        }
        return result;
    }

    bool& same() {return same_;}
    bool same() const {return same_;}
    Exchange& same(const bool flag) {
        same_=flag;
        return *this;
    }

    Exchange& set_algorithm(const Algorithm& alg ) {
    	algorithm_=alg;
        return *this;
    }


private:

    /// exchange using macrotasks, i.e. apply K on a function in individual worlds
    vecfuncT K_macrotask_efficient(const vecfuncT& vket, const double mul_tol=0.0) const;

    /// computing the full square of the double sum (over vket and the K orbitals)
    vecfuncT K_small_memory(const vecfuncT& vket, const double mul_tol=0.0) const;

    /// computing the upper triangle of the double sum (over vket and the K orbitals)
    vecfuncT K_large_memory(const vecfuncT& vket, const double mul_tol=0.0) const;

    /// computing the upper triangle of the double sum (over vket and the K orbitals)
    static vecfuncT compute_K_tile(World& world, const vecfuncT& mo_bra, const vecfuncT& mo_ket,
                    const vecfuncT& vket, std::shared_ptr<real_convolution_3d> poisson,
                    const bool same, const Tensor<double>& occ, const double mul_tol=0.0, double trunc_tol=0.0);

    World& world;
    bool same_=false;
    vecfuncT mo_bra, mo_ket;    ///< MOs for bra and ket
    Tensor<double> occ;
    double lo=1.e-4;
    double econv=FunctionDefaults<NDIM>::get_thresh();

};


} /* namespace madness */

#endif /* SRC_APPS_CHEM_EXCHANGEOPERATOR_H_ */
