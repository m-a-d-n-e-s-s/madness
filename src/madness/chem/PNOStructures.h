/*
 * PNOStructures.h
 *
 *  Created on: Oct 24, 2018
 *      Author: kottmanj
 */

#ifndef PAPER_CODE_PNOSTRUCTURES_H_
#define PAPER_CODE_PNOSTRUCTURES_H_

#include <valarray>
#include <string>

#include<madness/chem/PNOParameters.h>
#include<madness/chem/PNOTensors.h>
#include <madness/world/madness_exception.h>
#include <madness.h>

/// this block needs to be known before nonlinsol.h is included (happens when nemo is included)
/// better move this to nonlinsol.h
namespace madness {
/// inner product for std::valarray
double inner(const std::valarray<double>& bra, const std::valarray<double>& ket);
// KAIN allocator for std::valarray
struct valarray_allocator{
	const size_t dim;

	valarray_allocator(const size_t& indim) : dim(indim)  {};

	std::valarray<double> operator()(){
		return std::valarray<double>(dim);
	}
	valarray_allocator operator=(const valarray_allocator &other){
		valarray_allocator tmp(other.dim);
		return tmp;
	}

};
}


#include<madness/chem/nemo.h>
namespace madness {

class ParametrizedExchange {
public:
	ParametrizedExchange(World & world, const Nemo& nemo_, const std::string& type_):
	world(world),
	nemo(nemo_),
	type(type_),
	K(Exchange<double,3>(world, &nemo, 0))
	{

	}

	World& world;
	const Nemo& nemo;
	const std::string type;
	Exchange<double,3> K;

    real_function_3d operator()(const real_function_3d& ket) const {
        vecfuncT vket(1,ket);
        vecfuncT vKket=this->operator()(vket);
        return vKket[0];
    }

	vecfuncT operator ()(const vecfuncT& vket,
			const double& mul_tol = 0.0) const;



};

/// Timer structure
struct MyTimer{
	World& world;
	mutable double wstart_, cstart_;
	mutable double wstop_, cstop_;
	mutable double wtime_, ctime_;
	mutable std::string msg_;

	MyTimer(World& world): world(world), wstart_(-1.0),cstart_(-1.0),wstop_(-1.0),cstop_(-1.0),wtime_(-1.0),ctime_(-1.0){}
	MyTimer start()const{
		world.gop.fence();
		wstart_ = wall_time();
		cstart_ = cpu_time();
		return *this;
	}
	MyTimer stop(){
		wstop_ = wall_time();
		cstop_ = cpu_time();
		wtime_ = wstop_ - wstart_;
		ctime_ = cstop_ - cstart_;
		return *this;
	}
	MyTimer print(const std::string& msg)const{
		msg_=msg;
		if(world.rank()==0){
			printf("timer: %40.40s %8.2fs %8.2fs \n", msg.c_str(), wtime_, ctime_);
		}
		return *this;
	}
	MyTimer print()const{
		if(world.rank()==0){
			printf("timer: %40.40s %8.2fs %8.2fs \n", msg_.c_str(), wtime_, ctime_);
		}
		return *this;
	}
};

/// POD structure for energies
struct PairEnergies{
	PairEnergies(){}
	PairEnergies(const size_t& npairs): eijs(npairs),eijt(npairs),eijs_f12(npairs),eijt_f12(npairs),eij(npairs),energy(0.0),energy_f12(0.0){}

	std::valarray<double> eijs; ///< singlet pair energies (for CIS(D) the GS Part)
	std::valarray<double> eijt; ///< triplet pair energies (for CIS(D) the ES Part)
	std::valarray<double> eijs_f12; ///< singlet f12 pair energies (for CIS(D) the GS Part)
	std::valarray<double> eijt_f12; ///< triplet f12 pair energies (for CIS(D) the ES Part)
	std::valarray<double> eij; ///< total pair energies
	double energy=0.0; ///<total correlation energy (regularized Energy for f12 calculation)
	double energy_f12=0.0; ///< total f12 correlation energy
	double total_energy()const{return energy+energy_f12;}

	PairEnergies operator +=(const PairEnergies& right);

	PairEnergies operator +(const PairEnergies& right) const;

	void update(){
		energy=0.0;
		energy_f12=0.0;
		for(size_t ij=0;ij<eij.size();++ij){
			const double ereg_ij=eijs[ij]+eijt[ij];
			const double ef12_ij=eijs_f12[ij]+eijt_f12[ij];
			energy+=ereg_ij;
			energy_f12+=ef12_ij;
			eij[ij]=ereg_ij+ef12_ij;
		}
	}
};

/// iterates the third index for pair coupling
struct OrbitalIterator{
private:
	size_t i_;
	const size_t start_;
	const size_t stop_;
	bool finished_;
	const size_t freeze_;
public:
	size_t i()const {return i_;}

	OrbitalIterator(const size_t &nocc, const size_t &freeze): i_(0), start_(0),stop_(nocc-freeze),finished_(false),freeze_(freeze) {}

	operator bool() const{ return !finished_;}

	OrbitalIterator& operator ++(){
		if(i_<(stop_-1)) ++i_;
		else finished_=true;
		return *this;
	}

	// gives real index (if some orbtials are frozen)
	std::string name()const {return std::to_string(i_+start_);}

        size_t freeze() const {return freeze_;} // mostly to con the compiler that freeze_ is used
};

/// Iterator over pairs
/// iterates (i,j) from i=start to i<stop and j=i to j<stop
/// iterates like the pno-mp2.cc code
/// For frozen orbitals iteration goes from 0 to nocc-freeze !!! (consistency with pno-mp2.cc code)
/// The name() function gives back the "real" pair number
struct ElectronPairIterator{

	ElectronPairIterator(const size_t& nocc, const size_t& freeze):
		stop_(nocc-freeze),
		freeze_(freeze)
	{MADNESS_ASSERT(start_<stop_);}

	/// check if iteration has to proceed
	operator bool() const {
		return !finished_;
	}

	/// Gives the number of the pair in the valarray of pno-mp2.cc file
	template <typename Int>
	size_t tridx(Int row, Int col) {
		size_t result= (row < col) ? (row + col * (col + 1) / 2)
				: (col + row * (row + 1) / 2);
		return result;//-offset_;
	}

	/// pre-increment operator
	ElectronPairIterator& operator ++();

	/// number of active occupied orbitals
	size_t nocc()const{
		return (stop_-start_);
	}

	/// total number of pairs n * (n + 1) / 2;
	size_t npairs()const{
		const size_t n=stop_-start_;
		return n*(n+1)/2;
	}

	/// Gives back "pairij" (for frozen core i and j are the true indices corresponding to the reference with all orbitals)
	std::string name()const{
		std::string name="pair"+ std::to_string(i_+freeze_) + std::to_string(j_+freeze_);
		return name;
	}

	size_t i()const{ return i_;}
	size_t j()const{ return j_;}
	size_t start()const{ return start_;}
	size_t stop()const{ return stop_;}
	size_t ij()const{return ij_;}
	bool finished()const{return finished_;}
	bool diagonal()const{return (i_==j_);}

private:
	size_t i_=0;	///< current pair index i
	size_t j_=0;	///< current pair index j
	const size_t start_=0; ///< start value for i and j (usually 0)
	const size_t stop_=0;  ///< stop value for i and j (usually number of occ orbitals or occ-nfreeze for frozen core)
	const size_t freeze_=0; ///< number of frozen orbitals (just to print the correct name)
	size_t ij_=0;	///< pair number starting from 00 => 0
	bool finished_=false; ///< true if all pairs where iterated
};

/// POD for CIS excitation
struct CISData{
	CISData(){}
	CISData(const size_t& n, const double& o, const vector_real_function_3d& f) : x(f),number(n),omega(o){}
	vector_real_function_3d x=vector_real_function_3d();
	vector_real_function_3d Kx=vector_real_function_3d();
	vector_real_function_3d Vx=vector_real_function_3d();
	int number=-1;
	double omega=0.0;
	bool initialized()const{return (number>=0 and omega>0.0 and x.size()>0);}
};

/// POD for PNO code
struct PNOPairs{

	struct MemInfo{
		double pno;
		double Kpno;
		double W;
		friend std::ostream& operator <<(std::ostream& os, const MemInfo& mi){
			os << "PNOPairs Memory Information\n";
			os << "Total: " << std::fixed << mi.pno + mi.Kpno + mi.W << " GB \n";
			os << "PNO  : " << std::fixed << mi.pno << " GB \n";
			os << "KPNO : " << std::fixed << mi.Kpno << " GB \n";
			os << "W    : " << std::fixed << mi.W << " GB \n";
			return os;
		}
	};

	typedef std::valarray<vector_real_function_3d> vfT ;
	PNOPairs(const PairType& t, const size_t& n) :  type(t), nocc(n), npairs(n*(n+1)/2) , S_ij_ik(n), S_ij_kj(n) {initialize(n);}
	PNOPairs(const PairType& t, const size_t& n, const CISData& c): type(t), nocc(n), npairs(n*(n+1)/2), cis(c), S_ij_ik(n), S_ij_kj(n)
	{
		initialize(n);
		const size_t n2=cis.x.size();
		MADNESS_ASSERT(nocc==n);
		MADNESS_ASSERT(n==n2);
	}

	void initialize(const std::size_t& nocc);

	const PairType type;							///< type (i.e. MP2_PAIRTYPE, CISPD_PAIRTYPE, later CCSD etc
	const size_t nocc;								///< number of active occupied orbitals
	const size_t npairs;							///< number of Pairs
	CISData cis;									///< CIS excitation structure
	vfT pno_ij;										///< the PNOs for all Pairs
	std::valarray<Tensor<double> > rdm_evals_ij;   	        ///< PNO eigenvalues
	std::valarray<Tensor<double> > t_ij;			///< the amplitudes for all Pairs
	std::valarray<Tensor<double> > F_ij;			///< Fock matrices of the PNOs
	std::valarray<Tensor<double> > W_ij;			///< Fluctuation matrix
	std::valarray<int> maxranks_ij;					///< maxranks for all pairs, negative->unlimited
	std::valarray<bool> frozen_ij;					///< if true then pairs are frozen and not optimized (but still contribute to energy -> not frozen core)
	PairEnergies energies;							///< all Pair Energies
	vfT Kpno_ij; 									///< Exchange Intermediate
	vfT W_ij_i;										///< Fluctuation Potential
	vfT W_ij_j;										///< Fluctuation Potential
        PNOTensors::Tensor_IJ_IK<double> S_ij_ik;					///< PNO overlaps
        PNOTensors::Tensor_IJ_KJ<double> S_ij_kj;					///< PNO overlaps
	mutable MemInfo meminfo;						///< information about the used memory

	/// check if all pairs are empty
	bool empty()const{
		bool empty=true;
		for(const auto& p:pno_ij) empty=(empty and p.empty());
		return empty;
	}

	/// consistency check
	bool is_consistent(std::string& errm) const;
	/// throes exception if consistency check fails fails
	void verify()const{
		std::string errm="";
		if(not is_consistent(errm)) MADNESS_EXCEPTION(errm.c_str(),1);
	}
	// need explicit assignment operator bc of const members
	PNOPairs operator =(const PNOPairs& other);

	// name the pair to store and load on disc
	std::string name(const ElectronPairIterator& it) const;

	// rearrange a valarray to a big vector according to the pair structure of this
	// only return unfrozen pairs
	vector_real_function_3d extract(const vfT& vf) const;
	// reassemble big vector into valarray pair structure of this
	// the inverted operation to extract
	vfT reassemble(const vector_real_function_3d& v){
		vfT result(npairs);
		return reassemble(v,result);
	}
	vfT reassemble(const vector_real_function_3d& v, vfT& result) const;

	void clear_intermediates(const ElectronPairIterator& it);

	MemInfo update_meminfo() const;



};


} /* namespace madness */

#endif /* PAPER_CODE_PNOSTRUCTURES_H_ */
