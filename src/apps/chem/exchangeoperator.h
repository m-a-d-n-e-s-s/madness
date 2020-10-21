


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

typedef std::vector<real_function_3d> vecfuncT;


template<typename T, std::size_t NDIM>
class Exchange {
	typedef Function<T,NDIM> functionT;
	typedef std::vector<functionT> vecfuncT;

public:
	class MacroTaskExchange : public MacroTaskIntermediate<MacroTaskExchange> {

		typedef std::vector<std::shared_ptr<MacroTaskBase> > taskqT;

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
			vecfuncT mo_bra(nocc), mo_ket(nocc);
			for (int i=0; i<nocc; ++i) {
				mo_bra[i]=cloud.load<functionT>(subworld,i);
				mo_ket[i]=cloud.load<functionT>(subworld,i+nocc);
			}

			// psif is a flattened vector (f_i, bra)
			vecfuncT psif;
			const long batchsize=vf.size();
			for (auto f : vf) psif=append(psif,mul_sparse(subworld, f, mo_bra, mul_tol)); /// was vtol

			truncate(subworld, psif);
			psif = apply(subworld, *poisson.get(), psif);
			truncate(subworld, psif);

			vecfuncT Kf(batchsize);
			auto it=psif.begin();
			for (int i=0; i<batchsize; ++i) {

				vecfuncT psif_slice(it,it+mo_ket.size());
				Kf[i]=dot(subworld,mo_ket,psif_slice).truncate();
			}
//	        psif = mul_sparse(subworld, mo_ket[i], psif, mul_tol); /// was vtol
//	        gaxpy(world, 1.0, Kf, occ[i], psif);

			subworld.gop.fence();
			cloud.store(subworld,Kf,outputrecord);
		}

	    template <typename Archive>
	    void serialize(const Archive& ar) {
	    	ar & inputrecord & outputrecord & nocc & lo & econv & mul_tol;
	    }

	    void print_me(std::string s="") const {
	    	print("K apply task", s, inputrecord, outputrecord, this->stat);
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

    World& world;
    bool small_memory_=true;
    bool same_=false;
    vecfuncT mo_bra, mo_ket;    ///< MOs for bra and ket
    Tensor<double> occ;
    std::shared_ptr<real_convolution_3d> poisson;
public:

    bool multiworld_=false;
    long ntask_per_subworld=4;

};


} /* namespace madness */

#endif /* SRC_APPS_CHEM_EXCHANGEOPERATOR_H_ */
