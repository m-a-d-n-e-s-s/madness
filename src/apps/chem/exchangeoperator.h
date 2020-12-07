


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

public:
	class MacroTaskExchange : public MacroTaskIntermediate<MacroTaskExchange> {

		typedef std::vector<std::shared_ptr<MacroTaskBase> > taskqT;

		static inline std::shared_ptr<vecfuncT> mo_ket, mo_bra;
	public:


		long inputrecord=0;
        long outputrecord=0;
		long nocc=0;
		double lo=1.e-4;
		double econv=1.e-6;
		double mul_tol=1.e-7;

		MacroTaskExchange(const long inputrecord, const long outputrecord,
				const long nocc,
				const double lo, const double econv, const double mul_tol)
			: inputrecord(inputrecord)
			, outputrecord(outputrecord)
			, nocc(nocc)
			, lo(lo)
			, econv(econv)
			, mul_tol(mul_tol) {
		}

		void run(World& subworld, Cloud& cloud, taskqT& taskq) {
	    	auto poisson = std::shared_ptr<real_convolution_3d>(
	    	            CoulombOperatorPtr(subworld, lo, econv));

			// load the K operator argument and the orbitals of K
	    	vecfuncT vf=cloud.load<vecfuncT> (subworld,inputrecord);


//			for (int i=0; i<nocc; ++i) {
//				(*mo_bra)[i]=cloud.load<functionT>(subworld,i);
//				(*mo_ket)[i]=cloud.load<functionT>(subworld,i+nocc);
//			}

	    	if (not mo_bra.get()) {
				mo_bra.reset(new vecfuncT(nocc));
				for (int i=0; i<nocc; ++i) {
					(*mo_bra)[i]=cloud.load<functionT>(subworld,i);
				}
//	    		mo_bra.reset(new vecfuncT(cloud.load<vecfuncT>(subworld,0)));
	    	}
	    	if (not mo_ket.get()) {
				mo_ket.reset(new vecfuncT(nocc));
				for (int i=0; i<nocc; ++i) {
					(*mo_ket)[i]=cloud.load<functionT>(subworld,i+nocc);
				}
//	    		mo_ket.reset(new vecfuncT(cloud.load<vecfuncT>(subworld,1)));
	    	}

	    	double cpu0=cpu_time();
			// psif is a flattened vector (f_i, bra)
			vecfuncT psif;
			const long batchsize=vf.size();
			for (auto f : vf) psif=append(psif,mul_sparse(subworld, f, *mo_bra, mul_tol)); /// was vtol
            truncate(subworld, psif);
			double cpu1=cpu_time();
			double mul1=cpu1-cpu0;

			cpu0=cpu_time();
			psif = apply(subworld, *poisson.get(), psif);
			truncate(subworld, psif);
            cpu1=cpu_time();
            double apply1=cpu1-cpu0;

            cpu0=cpu_time();
			vecfuncT Kf(batchsize);
			auto it=psif.begin();
			for (int i=0; i<batchsize; ++i) {

				vecfuncT psif_slice(it,it+mo_ket->size());
				Kf[i]=dot(subworld,*mo_ket,psif_slice).truncate();
				it+=mo_ket->size();
			}
            cpu1=cpu_time();
            double dot1=cpu1-cpu0;
//	        psif = mul_sparse(subworld, mo_ket[i], psif, mul_tol); /// was vtol
//	        gaxpy(world, 1.0, Kf, occ[i], psif);

			subworld.gop.fence();
			cloud.store(subworld,Kf,outputrecord);
			printf("timings for mul1, apply, dot: %8.2fs %8.2fs %8.2fs\n",mul1,apply1,dot1);
		}

		void cleanup() {
			if (mo_ket) {
				print("clearing mo_ket");
				mo_ket.reset();
			}
			if (mo_bra) mo_bra.reset();
		}

	    template <typename Archive>
	    void serialize(const Archive& ar) {
	    	ar & inputrecord & outputrecord & nocc & lo & econv & mul_tol;
	    }

	    void print_me(std::string s="") const {
	    	print("K apply task", s, inputrecord, outputrecord, this->stat);
	    }

	};

    class MacroTaskExchangeEfficient : public MacroTaskIntermediate<MacroTaskExchangeEfficient> {

        typedef std::vector<std::shared_ptr<MacroTaskBase> > taskqT;

        static inline std::shared_ptr<vecfuncT> mo_ket, mo_bra;
    public:

        std::pair<long,long> row_range, column_range;
        long nocc=0;
        double lo=1.e-4;
        double econv=1.e-6;
        double mul_tol=1.e-7;

        MacroTaskExchangeEfficient(const std::pair<long,long>& row_range, const std::pair<long,long>& column_range,
                          const long nocc,
                          const double lo, const double econv, const double mul_tol)
                : row_range(row_range)
                , column_range(column_range)
                , nocc(nocc)
                , lo(lo)
                , econv(econv)
                , mul_tol(mul_tol) {
        }

        void run(World& subworld, Cloud& cloud, taskqT& taskq) {
            auto poisson = std::shared_ptr<real_convolution_3d>(
                    CoulombOperatorPtr(subworld, lo, econv));

            // the argument of the exchange operator is the ket vector

            // load bra and ket if not already loaded
            if (not mo_bra.get()) {
                mo_bra.reset(new vecfuncT(nocc));
                for (int i=0; i<nocc; ++i) (*mo_bra)[i]=cloud.load<functionT>(subworld,i);
            }
            if (not mo_ket.get()) {
                mo_ket.reset(new vecfuncT(nocc));
                for (int i=0; i<nocc; ++i) (*mo_ket)[i]=cloud.load<functionT>(subworld,i+nocc);
            }

            // compute the tile [column_range,row_range], corresponding to bra[nrow], ket[ncolumn]

            // multiply the functions orbital_product_{ij}(r) = i(r) j(r)
            double cpu0=cpu_time();

            vecfuncT bra_batch(mo_bra->begin()+row_range.first,mo_bra->begin()+row_range.second);
            vecfuncT ket_batch(mo_ket->begin()+column_range.first,mo_ket->begin()+column_range.second);

            // orbital_product is a vector of vectors
            std::vector<vecfuncT> orbital_product=matrix_mul_sparse<T,T,NDIM>(subworld,bra_batch,ket_batch,mul_tol);
            vecfuncT orbital_product_flat=flatten(orbital_product); // convert into a flattened vector
            truncate(subworld, orbital_product_flat);
            double cpu1=cpu_time();
            double mul1=cpu1-cpu0;

            cpu0=cpu_time();
            vecfuncT Nij = apply(subworld, *poisson.get(), orbital_product_flat);
            truncate(subworld, Nij);
            cpu1=cpu_time();
            double apply1=cpu1-cpu0;

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
                for (int i=s.start; i<=s.end+ncolumn; ++i) {
//                    print("Nslice: pushing back element", ij(irow,i), "for row",irow);
                    result.push_back(Nij[ij(irow,i)]);
                }
                return result;
            };
            auto Nslice1 = [&Nij, &ij, &nrow] (const Slice s, const long jcolumn) {
                vecfuncT result;
                MADNESS_CHECK(s.start==0 && s.end==-1 && s.step==1);
                for (int i=s.start; i<=s.end+nrow; ++i) {
//                    print("Nslice: pushing back element", ij(i,jcolumn), "for column",jcolumn);
                    result.push_back(Nij[ij(i,jcolumn)]);
                }
                return result;
            };


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
            double dot1=cpu1-cpu0;

            subworld.gop.fence();

            // store results: columns as columns; transpose rows to columns
            cloud.store(subworld,resultcolumn,outputrecord(row_range,column_range));
            if (row_range!=column_range) {  // store diagonal batches only once
                cloud.store(subworld,resultrow,outputrecord(column_range,row_range));
            }
            print("storing the results for batch",row_range,column_range,
                  "in records",outputrecord(row_range,column_range),
                "and",outputrecord(column_range,row_range));
            printf("timings for mul1, apply, dot: %8.2fs %8.2fs %8.2fs\n",mul1,apply1,dot1);
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
        static long outputrecord(std::pair<long,long>& rowrange, std::pair<long,long>& columnrange) {
            // Cantor pairing
            auto cantor = [] (const long k1, const long k2) {
                return ((k1 + k2) * (k1 + k2 + 1)/2 + k2);
            };
            return cantor(rowrange.first,columnrange.first);
        }

        template <typename Archive>
        void serialize(const Archive& ar) {
            ar & row_range & column_range & nocc & lo & econv & mul_tol;
        }

        void print_me(std::string s="") const {
            print("K apply task", s, this->stat);
        }

    };
public:

    /// default ctor
    Exchange(World& world) : world(world), small_memory_(true), same_(false) {};

    /// ctor with a conventional calculation
    Exchange(World& world, const SCF* calc, const int ispin);

    /// ctor with a nemo calculation
    Exchange(World& world, const Nemo* nemo, const int ispin);

    /// set the bra and ket orbital spaces, and the occupation

    /// @param[in]	bra		bra space, must be provided as complex conjugate
    /// @param[in]	ket		ket space
    /// @param[in]	occ1	occupation numbers
    void set_parameters(const vecfuncT& bra, const vecfuncT& ket,
            const Tensor<double>& occ1, const double lo=1.e-4,
            const double econv=FunctionDefaults<3>::get_thresh()) {
    	mo_bra=copy(world,bra);
    	mo_ket=copy(world,ket);
    	occ=copy(occ1);
    	poisson = std::shared_ptr<real_convolution_3d>(
    	            CoulombOperatorPtr(world, lo, econv));
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

    /// @param[in]  bra    real_funtion_3d, the bra state
    /// @param[in]  ket    real_funtion_3d, the ket state
    T operator()(const Function<T,NDIM>& bra, const Function<T,NDIM>& ket) const {
        return inner(bra,this->operator()(ket));
    }

    /// compute the matrix < vbra | K | vket >

    /// @param[in]  vbra    vector of real_funtion_3d, the set of bra states
    /// @param[in]  vket    vector of real_funtion_3d, the set of ket states
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

    bool& small_memory() {return small_memory_;}
    bool small_memory() const {return small_memory_;}
    Exchange& small_memory(const bool flag) {
        small_memory_=flag;
        return *this;
    }

    bool& same() {return same_;}
    bool same() const {return same_;}
    Exchange& same(const bool flag) {
        same_=flag;
        return *this;
    }

    Exchange& multiworld(const bool flag) {
    	multiworld_=flag;
        return *this;
    }


private:

    /// exchange using macrotasks, i.e. apply K on a function in individual worlds
    vecfuncT K_macrotask(const vecfuncT& vket, const double mul_tol=0.0) const;

    /// exchange using macrotasks, i.e. apply K on a function in individual worlds
    vecfuncT K_macrotask_efficient(const vecfuncT& vket, const double mul_tol=0.0) const;

    /// computing the full square of the double sum (over vket and the K orbitals)
    vecfuncT K_small_memory(const vecfuncT& vket, const double mul_tol=0.0) const;

    /// computing the upper triangle of the double sum (over vket and the K orbitals)
    vecfuncT K_large_memory(const vecfuncT& vket, const double mul_tol=0.0) const;

    World& world;
    bool small_memory_=true;
    bool same_=false;
    vecfuncT mo_bra, mo_ket;    ///< MOs for bra and ket
    Tensor<double> occ;
    std::shared_ptr<real_convolution_3d> poisson;
public:

    bool multiworld_=false;
    bool efficient_=false;
    long ntask_per_subworld=4;

};


} /* namespace madness */

#endif /* SRC_APPS_CHEM_EXCHANGEOPERATOR_H_ */
