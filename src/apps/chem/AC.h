/*
 * AC.h
 *
 *  Created on: Dec 13, 2016
 *      Author: msahre
 */

#ifndef SRC_APPS_CHEM_AC_H_
#define SRC_APPS_CHEM_AC_H_

#include <madness/world/MADworld.h>
#include <madness/mra/mra.h>
#include <madness/mra/funcplot.h>
#include <exception>
#include <iterator>
#include <list>
#include <chem/SCF.h>

namespace madness {



/// Needed information about atom to compute asymptotic correction
template<unsigned long int NDIM>
struct atom_information{

	/// Coordinates of nucleus
	Vector<double,NDIM> coord;

	/// Interval boundaries, in which exchange correlation potential (<R1) and 1/r potential (>R2) are used
	double R1;
	double R2;

	/// Nuclear charge
	int charge;

	/// Computes the distance between electron and nucleus
	/// @param[in] elec: position of electron
	/// @param[in] nuc: position of nucleus
	/// @param[out] distance: distance between electron and nucleus
	double get_distance(Vector<double,NDIM> elec, Vector<double,NDIM> nuc) const{
		double distance = 0.0;
		for(unsigned i = 0; i < NDIM; i++){
			distance += (elec[i]-nuc[i])*(elec[i]-nuc[i]);
		}
		distance = sqrt(distance);

		return distance;
	}

    template <typename Archive>
    void serialize(Archive& ar) {
        ar &  coord& R1&  R2& charge;
    }

};

/// Creates vector of atom_information for the given molecule
/// @param[in] molecule: the molecule
/// @param[out] : the vector containing the atom_information for each atom in the molecule
std::vector<atom_information<3> > make_atom_vec(const Molecule& molecule, double R1_,double R2_);

/// Creates vector of atom_information via input file (calls make_atom_vec(const Molecule& molecule))
/// @param[in] input: file with geometry of molecule
/// @param[out] : the vector containing the atom_information for each atom in the molecule
std::vector<atom_information<3> > make_atom_vec(const std::string& input);

/// collection of slater radii needed to calculate R1 and R2; The Journal of Chemical Physics 41, 3199 (1964); doi: 10.1063/1.1725697
/// @param[in] atomic_number: atomic_number of atom
/// @param[out] : slater radius corresponding to atomic_number
double slater_radius(int atomic_number);




/// Contains all the parameters for the asymptotic correction
template<unsigned long int NDIM>
struct ACParameters{
	std::vector< atom_information<NDIM> > atoms_;
	bool use_mult_; /// set true to use multipole approximation
	double e_ion_; /// ionisation energy of corresponding ion without ac
	double eh_; /// energy of the homo without ac
	double R1_; /// lower boundary interpolation area
	double R2_; /// upper boundary for interpolation area
	double dft_coefficient_; /// dft coefficient for hybrid functionals (=1.0 for non hybrid)
	int num_elec_; /// number of electrons
	std::string interpolation_scheme_; /// scheme for interpolation between xc-functional and 1/r-potential


    template <typename Archive>
    void serialize(Archive& ar) {
        ar & atoms_ & e_ion_ & eh_ & R1_ & R2_ & dft_coefficient_ & num_elec_ & interpolation_scheme_;
    }

	ACParameters(){
		use_mult_ = false;
		e_ion_ = 0.0;
		eh_ = 0.0;
		R1_ = 3.0;
		R2_ = 4.0;
		interpolation_scheme_ = "linear";
		num_elec_=-1;
		dft_coefficient_=0.0;
	}

	bool initialize(Molecule molecule, std::string ac_data, double dft_coeff, double tot_charge){
		if(ac_data=="none"){
			return false;
		}else{
		dft_coefficient_ = dft_coeff;

		std::stringstream sac_data(ac_data);
		std::string word;

		while (sac_data >> word) {
			std::transform(word.begin(), word.end(), word.begin(), ::toupper);
			if(word=="MULT"){
				use_mult_ = true;
			} else if(word == "EION") {
				sac_data >> e_ion_;
			} else if(word == "EHOMO"){
				sac_data >> eh_;
			} else if(word == "R1"){
				sac_data >> R1_;
			} else if(word == "R2"){
				sac_data >> R2_;
			} else if(word == "INTERPOL"){
				sac_data >> interpolation_scheme_;
			} else {
				std::cout << "Invalid entry in ac line\n";
				//MADNESS_ASSERT("Invalid ac_data!\n", 1);
			}
		}
		atoms_ = make_atom_vec(molecule, R1_, R2_);
		num_elec_ = molecule.total_nuclear_charge() - tot_charge;
		return true;
		}
	}

	ACParameters(const ACParameters &other)
	: atoms_(other.atoms_), use_mult_(other.use_mult_), e_ion_(other.e_ion_),
	  eh_(other.eh_), R1_(other.R1_), R2_(other.R2_), dft_coefficient_(other.dft_coefficient_),
	  num_elec_(other.num_elec_), interpolation_scheme_(other.interpolation_scheme_) {}

	void print(World& world){
		if(world.rank()==0){
		std::cout << "\nAsymptotic correction parameters\n";
		std::cout << "---------------------------------------\n";
		std::cout << "Atom vector is not empty: " << !atoms_.empty() << "\n";
		std::cout << "Using multipole approximation: " << use_mult_ << "\n";
		std::cout << "Number of electrons = " << num_elec_ << "\n";
		std::cout << "Ionisation energy = " << e_ion_ << "\n";
		std::cout << "E(HOMO) = " << eh_ << "\n";
		std::cout << "R1 = " << R1_ << "\n";
		std::cout << "R2 = " << R2_ << "\n";
		std::cout << "DFT coefficient = " << dft_coefficient_ << "\n";
		std::cout << "Interpolation scheme = " << interpolation_scheme_ << "\n";
		}
	}

	void check(World& world){
		if(atoms_.empty()) { MADNESS_EXCEPTION("Failed to initialize AC object!",1); }
		else if(num_elec_ <0) { MADNESS_EXCEPTION("Negative number of electrons!",1); }
		else if(e_ion_ < 0) { MADNESS_EXCEPTION("Ionisation energy is negative!",1); }
		else if(eh_>0) { MADNESS_EXCEPTION("Energy of homo is positive!",1); }
		else if(R1_ == 0.0 or R2_ == 0.0) std::cout << "\n\nWARNING: R1 or R2 is zero!\n\n";
		else if(interpolation_scheme_ != "linear" and interpolation_scheme_ != "constant") std::cout << "\n\nWARNING: Unknown interpolation scheme, using linear interpolation instead\n\n!";
		else if(world.rank()==0) std::cout << "AC object was initialized succesfully!\n\n";
	}

};

/// Functor for the correction factor, which is multiplied with the exchange correlation potential
template<unsigned long int NDIM>
class int_factor_functor : public FunctionFunctorInterface<double,NDIM>{

public:
	int_factor_functor(){}
	//int_factor_functor(std::vector<atom_information<NDIM> > atoms, std::string interpolation): atoms(atoms), intpol_scheme(interpolation){}

	int_factor_functor(const ACParameters<NDIM>& ac_param) : ac_param_(ac_param) {}

	//int_factor_functor(std::string interpolation) : intpol_scheme(interpolation) {}



	/// Computes correction factor for exchange correlation potential
	/// @param[in] r: position of the electron
	/// @param[out] : correction factor for exchange correlation potential at position r
	double operator ()(const Vector<double, NDIM> &r)const{

		// distance between electron and nuclues
		double d = 0.0;

		bool inside_R1 = false;
		std::vector<atom_information<NDIM> > between;

		// get position of electron relative to nuclei, to compute right correction factor
		for(const auto& x:ac_param_.atoms_){
			d = get_distance(r, x.coord);
			// check if d<R1 for at least one nucleus
			if(d<x.R1) {
				inside_R1 = true;
				break;
			}
			// push_back index of nucleus for which electron is inbetween
			if(x.R1 < d and d < x.R2) {
				between.push_back(x);
			}
		}

		// compute correction factor
		// at least for one nucleus r<R1
		if(inside_R1){
			return 1.0;}
		// all nuclei r>R2
		else if(between.empty()){
			return 0.0;}
		// for N nuclei R1<r<R2
		else {
			return int_factor(r, between);
		}
	}

	/// Computes correction factor for the exchange correlation potential if electron position is between R1 and R2
	/// @param[in] r: position of electron
	/// @param[in] between: vector of the atom_information of all atoms for that electron position is between R1 and R2
	/// @param[out] int_factor: correction factor
	double int_factor(const Vector<double, NDIM> &r, std::vector<atom_information<NDIM> > between) const{
		// Do not forget to add interpolation scheme to ACparameters check()
		if(ac_param_.interpolation_scheme_=="constant"){
			return 0.5;
		}
		else{
			double int_factor = 0.0;
			for(const auto& x:between){
				// linear interpolation scheme between R1 and R2
				int_factor += (- get_distance(r,x.coord) + x.R2)/(x.R2-x.R1);
				//int_factor += 0.5*(erf(-(get_distance(r, x.coord)-5.0))+1.0);
			}
			int_factor = int_factor/(double(between.size()));
			return int_factor;
		}
	}

private:
	/// Computes the distance between electron and nucleus
	/// @param[in] elec: position of electron
	/// @param[in] nuc: position of nucleus
	/// @param[out] distance: distance between electron and nucleus
	double get_distance(Vector<double,NDIM> elec, Vector<double,NDIM> nuc) const{
		double distance = 0.0;
		for(unsigned i = 0; i < NDIM; i++){
			distance += (elec[i]-nuc[i])*(elec[i]-nuc[i]);
		}
		distance = sqrt(distance);

		return distance;
	}

	/// Parameter for the asymtotic correction
	const ACParameters<NDIM> &ac_param_;

};

/// Functor for the 1/r potential to induce the correct asymptotic behaviour of the exchange correlation potential
template<unsigned long int NDIM>
class lr_pot_functor : public FunctionFunctorInterface<double,NDIM>{

public:
	lr_pot_functor(){}
	//lr_pot_functor(std::vector<atom_information<NDIM> > atoms, double dft_coeff_, std::string interpolation): atoms(atoms), dft_coeff(dft_coeff_), intpol_scheme(interpolation){}

	lr_pot_functor(const ACParameters<NDIM> & ac_param): ac_param_(ac_param){}

	/// Computes 1/r potential weighted by a correction factor
	/// @param[in] r: position of the electron (density)
	/// @param[out]: 1/r potential weighted by a correction factor
	double operator ()(const Vector<double, NDIM> &r)const{
		// distance between electron and nuclues
		double d = 0.0;
		bool inside_R1 = false;
		std::vector<atom_information<NDIM> > between;

		// get position of electron relative to nuclei, to compute right correction factor
		for(const auto& x:ac_param_.atoms_){

			d = get_distance(r, x.coord);
			// check if d<R1 for at least one nucleus
			if(d<x.R1) {
				inside_R1 = true;
				break;
			}

			// push_back index of nucleus for which electron is inbetween
			if(x.R1 < d and d < x.R2) {
				between.push_back(x);
			}

		}

		// compute correction factor
		// at least for one nucleus r<R1
		if(inside_R1){return 0.0;} // return 0.0;

		// all nuclei r>R2
		else if(between.empty()){return 1.0*potential(r);} // 1.0*potential(r)+E_ion+E_homo

		// for N nuclei R1<r<R2
		else {return (int_factor(r, between))*potential(r);}
	}

private:
	/// Computes the distance between electron and nucleus
	/// @param[in] elec: position of electron
	/// @param[in] nuc: position of nucleus
	/// @param[out] distance: distance between electron and nucleus
	double get_distance(Vector<double,NDIM> elec, Vector<double,NDIM> nuc) const{
		double distance = 0.0;
		for(unsigned i = 0; i < NDIM; i++){
			distance += (elec[i]-nuc[i])*(elec[i]-nuc[i]);
		}
		distance = sqrt(distance);

		return distance;
	}

	/// Multipole approximation of the 1/r potential for electron coordinate > R1 for all nuclei
	/// @param[in] r: coordinate of the electron
	/// @param[out] : potential at given coordinate r
	double potential(const Vector<double, NDIM> &r) const{
		double pot = 0.0;
		for(const auto& x:ac_param_.atoms_){
			pot += double(x.charge)/get_distance(r, x.coord);
		}
		return (-ac_param_.dft_coefficient_*pot/ac_param_.num_elec_);
	}

	/// Computes correction factor for the 1/r potential if electron position is between R1 and R2
	/// @param[in] r: position of electron
	/// @param[in] between: vector of the atom_information of all atoms for that electron position is between R1 and R2
	/// @param[out] int_factor: correction factor
	double int_factor(const Vector<double, NDIM> &r, std::vector<atom_information<NDIM> > between) const{
		int_factor_functor<NDIM> interpolation(ac_param_);
		double int_factor = interpolation.int_factor(r, between);
		return (1.0 - int_factor);
	}

	/// Parameter for the asymtotic correction
	const ACParameters<NDIM> &ac_param_;
};

/// computes the corrected exchange correlation potential using the multipole approximation

template<unsigned long int NDIM>
struct UnaryOpStructure {

	UnaryOpStructure(const std::shared_ptr<FunctionFunctorInterface<double,NDIM> > f_, const std::shared_ptr<FunctionFunctorInterface<double,NDIM> > f2_) :
		f(f_), f2(f2_),
		cdata(FunctionCommonData<double,NDIM>::get(FunctionDefaults<NDIM>::get_k())){};

	/// computes the corrected potential (weighting with correction factor and adding of 1/r potential)
	/// @param[in] t: uncorrected exchange correlation potential
	/// @param[out] t: corrected exchange correlation potential
    void operator()(const Key<NDIM>& key, Tensor<double>& t) const {
        Tensor<double> intpol(t.ndim(),t.dims());//madness::copy(t);
        Tensor<double> lrpot(t.ndim(),t.dims());
        const Tensor<double>& qp =cdata.quad_x;
        fcube(key,(*f),qp,intpol);
        fcube(key,(*f2),qp,lrpot);

        // multiply xc-functional with interpolation factor
        t.emul(intpol);
        // add 1/r correction to xc-functional
        t= t + lrpot;
    }
    /// shared pointer to object of int_factor_functor
    std::shared_ptr<FunctionFunctorInterface<double,NDIM> > f;
    /// shared pointer to object of lr_pot_functor
    std::shared_ptr<FunctionFunctorInterface<double,NDIM> > f2;
    FunctionCommonData<double,NDIM> cdata;
    template <typename Archive> void serialize(Archive& ar) {}
};

/// computes the corrected exchange correlation potential using the hartree potential

template<unsigned long int NDIM>
struct BinaryOpStructure {
	/// necessary to compile the program without getting a strange error i dont understand
	BinaryOpStructure(): f(NULL), cdata(FunctionCommonData<double,NDIM>::get(FunctionDefaults<NDIM>::get_k())) {}

	BinaryOpStructure(const std::shared_ptr<FunctionFunctorInterface<double,NDIM> > f_) :
		f(f_),
		cdata(FunctionCommonData<double,NDIM>::get(FunctionDefaults<NDIM>::get_k())){};

	BinaryOpStructure(const BinaryOpStructure<NDIM> &other) : f(other.f), cdata(other.cdata) {}

	/// computes the corrected potential (weighting with correction factor and adding of 1/r potential)
	/// @param[in] t: uncorrected exchange correlation potential
	/// @param[out] t: corrected exchange correlation potential
    void operator()(const Key<NDIM>& key, Tensor<double>& result, const Tensor<double>& vxc, const Tensor<double>& vhartree) const {

    	if(f.get()==NULL) MADNESS_EXCEPTION("NULL Pointer in BinaryOpStructure of AC",1);

    	Tensor<double> intpol(result.ndim(),result.dims());//madness::copy(t);
        Tensor<double> lrpot(result.ndim(),result.dims());
        const Tensor<double>& qp =cdata.quad_x;
        fcube(key,(*f),qp,intpol);

        // exchange correlation potential weighted by correction factor
        Tensor<double> tmp1 = madness::copy(vxc).emul(intpol);
        // hartree potential pot weighted by correction factor
        Tensor<double> tmp2 = madness::copy(vhartree).emul(intpol);
        Tensor<double> tmp3 = madness::copy(vhartree);
        // corrected potential
        result = tmp1 - tmp2 + tmp3;
    }

    /// shared pointer to object of int_factor_functor
    std::shared_ptr<FunctionFunctorInterface<double,NDIM> > f;
    FunctionCommonData<double,NDIM> cdata;
    template <typename Archive> void serialize(Archive& ar) {}
};

/// Asymptotic correction for DFT. In the correction the xc-potential is replaced by an 1/r term far away from the nuclei
/// to give the correct asymptotic behavior. Close to the nuclei the standard xc-potential is used. The transition between
/// the different potentials is achieved via a linear interpolation.
/// for more information see: Molecular Physics, Vol. 103, No. 2–3, 20 January–10 February 2005, 413–424;
/// The Journal of Chemical Physics 109, 10180 (1998); doi: 10.1063/1.477711;
/// Molecular Physics, 1999, VOL.97, No. 7, 859-868
template<unsigned long int NDIM>
class AC{
private:

	/// Parameter for the asymtotic correction
	ACParameters<NDIM> ac_param_;
	bool initialized_=false;
	double shift()const{return (-ac_param_.e_ion_-ac_param_.eh_);}

public:
    bool initialized()const{return initialized_;}
    AC() = default;

   AC(const ACParameters<NDIM>& ac_param): ac_param_(ac_param) {}

   AC(World &world, std::shared_ptr<SCF> calc){
	   if(world.rank()==0){
		   initialized_=ac_param_.initialize(calc->molecule, calc->param.ac_data(), 1.0-calc->xc.hf_exchange_coefficient(), calc->param.charge());
	   }
	   world.gop.broadcast_serializable(initialized_,0);
	   world.gop.broadcast_serializable(ac_param_, 0);
//	   ac_param_.print(world);
	   if(calc->param.ac_data()!="none") ac_param_.check(world);
   }

   AC(const AC& other): ac_param_(other.ac_param_){}

   /// /// correction of the exchange correlation potential using the multipole approximation to describe asymptotic behaviour
   /// @param[in] xc_functional: uncorrected (standard) exchange correlation potential
   /// @param[out] : corrected exchange correlation potential
   Function<double,NDIM> apply(Function<double,NDIM> xc_functional)const{
	   MADNESS_ASSERT(initialized());
	   std::cout << "Apply AC Scheme with multipole approximation\n";
	   // fast return
	   if(ac_param_.atoms_.empty()){
		   std::cout << "OR NOT -- EMPTY VECTOR ATOMS!!!\n";
		   return xc_functional;
	   }

	   if(ac_param_.dft_coefficient_ < 0.0) { MADNESS_EXCEPTION("DFT coefficient is negative!",1); }
	   if(ac_param_.dft_coefficient_ == 0.0) { MADNESS_EXCEPTION("DFT coefficient is zero. This is no DFT calculation!\n",1); }

	   // shift of the exchange correlation potential to get the correct asymptotic behaviour
	   xc_functional = xc_functional + shift();
	   // pointer on interpolation factor and long range potential (1/r pot)
	   std::shared_ptr<FunctionFunctorInterface<double, NDIM> > int_ptr(new int_factor_functor<NDIM>(ac_param_));
	   std::shared_ptr<FunctionFunctorInterface<double, NDIM> > lr_pot_ptr(new lr_pot_functor<NDIM>(ac_param_));

		// compute interpolation and long range potential
		UnaryOpStructure<NDIM> U_total(int_ptr,lr_pot_ptr);

		// apply correction on xc functional
		xc_functional.template unaryop<UnaryOpStructure<NDIM> >(U_total);
		//plot_plane( xc_functional, "test_total");
		return xc_functional;
   }

   /// correction of the exchange correlation potential using the hartree potential to describe asymptotic behaviour
   /// @param[in] xc_functional: uncorrected (standard) exchange correlation potential
   /// @param[in] v_hartree: potential to describe asymptotic behaviour
   /// @param[out] : corrected exchange correlation potential
   Function<double,NDIM> apply(Function<double,NDIM> xc_functional, const Function<double,NDIM> v_hartree)const{
	   MADNESS_ASSERT(initialized());
	   if(ac_param_.use_mult_) return apply(xc_functional);
	   //else return apply_potential(xc_functional, v_hartree);
	   std::cout << "Apply AC Scheme with hartree potential\n";
	   // fast return
	   if(ac_param_.atoms_.empty()){
		   std::cout << "OR NOT -- EMPTY VECTOR ATOMS!!!\n";
		   return xc_functional;
	   }
	   // shift of the exchange correlation potential to get the correct asymptotic behaviour
	   xc_functional = xc_functional + shift();

	   // pointer on interpolation factor and long range potential
	   std::shared_ptr<FunctionFunctorInterface<double, NDIM> > int_ptr(new int_factor_functor<NDIM>(ac_param_));

	   // compute interpolation and long range potential
	   BinaryOpStructure<NDIM> U_total(int_ptr);
	   xc_functional = binary_op<double,double,BinaryOpStructure<NDIM>,NDIM >(xc_functional, v_hartree, U_total);

	   return xc_functional;
   }

};

}// end namespace

#endif /* SRC_APPS_CHEM_AC_H_ */
