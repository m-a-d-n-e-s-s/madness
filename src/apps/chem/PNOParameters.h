/*
 * PNOParameters.h
 *
 *  Created on: Sep. 9, 2019
 *      Author: jsk
 */

#ifndef PNOPARAMETERS_H_
#define PNOPARAMETERS_H_

#include <vector>
#include <map>
#include <chem/CalculationParameters.h>

namespace madness {

enum PairType{MP2_PAIRTYPE,CISPD_PAIRTYPE,ALL_PAIRTYPE,NONE_PAIRTYPE,UNKNOWN_PAIRTYPE};
inline std::string type(const PairType& n){ return "PairType";}
std::ostream& operator << (std::ostream& os, const PairType& en);
std::istream& operator >> (std::istream& os, PairType& en);
enum EnergyType{PROJECTED_ENERGYTYPE,HYLLERAAS_ENERGYTYPE,UNKNOWN_ENERGYTYPE};
inline std::string type(const EnergyType& n){ return "EnergyType";}
std::ostream& operator << (std::ostream& os, const EnergyType& en);
std::istream& operator >> (std::istream& os, EnergyType& en);
enum GuessType{PARTIAL_WAVE_GUESSTYPE,FROM_FILE_GUESSTYPE,PREDEFINED_GUESSTYPE,SCF_GUESSTYPE,EXOP_GUESSTYPE,EMPTY_GUESSTYPE,PSI4_GUESSTYPE,UNKNOWN_GUESSTYPE};
inline std::string type(const GuessType& n){ return "GuessType";}
std::ostream& operator << (std::ostream& os, const GuessType& en);
std::istream& operator >> (std::istream& is, GuessType& en);


class ComputeProtocol: public QCCalculationParametersBase {
public:
	using QCCalculationParametersBase::read;
	ComputeProtocol(const QCCalculationParametersBase& param) : QCCalculationParametersBase(param){
		initialize_compute_protocol();
	}

	ComputeProtocol() : QCCalculationParametersBase(){
		initialize_compute_protocol();
	}

	ComputeProtocol(World& world, const std::string& inputfile, const std::string& TAG="computeprotocol") : QCCalculationParametersBase(){
		initialize_compute_protocol();
		read(world,inputfile, TAG);
		set_derived_parameters();
	}

	template<typename T>
	std::vector<T>initialize_vector(const T& one, const T& two)const{
		std::vector<T> result;
		result.push_back(one);
		result.push_back(two);
		return result;
	}

	void initialize_compute_protocol(){
		initialize<std::vector<double> >("thresh",initialize_vector(1.e-3, 1.e-4), "MRA threshold");
		initialize<std::vector<double> >("op_thresh",initialize_vector(1.e-7, 1.e-7), "MRA operator threshold");
		initialize<std::vector<int> >("maxiter_micro",initialize_vector(10, 10), "");
		initialize<std::vector<int> >("maxiter_macro",initialize_vector(10, 1), "");
		initialize<std::vector<double> >("econv_micro",initialize_vector(3.e-4, 1.e-4), "");
		initialize<std::vector<double> >("econv_macro",initialize_vector(3.e-4, 1.e-4), "");
		initialize<std::vector<std::string> >("exop",std::vector<std::string>{"multipole","none"}, "");
	}

	size_t size()const{
		return thresh().size();
	}

	void set_derived_parameters(){

	}

	std::vector<double> thresh()const { return get<std::vector<double> >("thresh");}
	std::vector<double> op_thresh()const { return  get<std::vector<double> >("op_thresh");}
	std::vector<int> maxiter_micro()const { return get<std::vector<int> >("maxiter_micro");}
	std::vector<int> maxiter_macro()const { return get<std::vector<int> >("maxiter_macro");}
	std::vector<double> econv_micro()const { return get<std::vector<double> >("econv_micro");}
	std::vector<double> econv_macro()const { return get<std::vector<double> >("econv_macro");}
	std::vector<std::string> exop()const { return get<std::vector<std::string> >("exop");}

};

class PNOParameters: public QCCalculationParametersBase {
public:

	template<typename T>
	T assign_from_string(const std::string& string)const{
		T result;
		std::stringstream ss(string);
		ss >> result;
		return result;
	}

	PNOParameters() : QCCalculationParametersBase(){
		initialize_pno_parameters();
	}

	PNOParameters(const QCCalculationParametersBase& param) : QCCalculationParametersBase(param){
		initialize_pno_parameters();
	}

	PNOParameters(World& world, const std::string& inputfile, const std::string& TAG="pno") : QCCalculationParametersBase(){
		initialize_pno_parameters();
		QCCalculationParametersBase::read(world,inputfile,TAG);
	}

	void initialize_pno_parameters() {
		initialize<int>("rank_increase", 10 , "maximum rank to increase in every macroiteration");
		initialize<int>("chunk", 100 , "chunk of functions operated on in parrallel when G or K is applied (prevent memory shortage)");
		initialize<bool>("debug",false, "debug mode");
		initialize<std::size_t>("freeze",0, "frozen core approximation");
		initialize<int>("maxrank", 999, "maximal pno rank for all pairs");
		initialize<std::string>("guesstype","exop", "the guesstype: exop (recommended, multiply polynomial excitation operators onto the occupied orbitals), empty (don't compute a guess), scf (use the ao basis from the SCF solver), partial_wave (create s,p,d ..., type functions and place on atoms), predefined (same as partial_wave but for predefined sets) ");
		initialize<std::string>("exop", "multipole", "this string defines which excitation operators can be used, use 'dipole', 'dipole+', 'quadrupole', 'ocopole', or 'multipole', need to set guesstype to 'exop' which is the default ");
		initialize<bool>("exop_trigo", true, "use trigonometric excitation operators ( x --> sin(x) ");
		initialize<std::string>("partial_wave", "", "atom centered partial wave guess of format like '3s2p1d', need to set guesstype to 'partial_wave' ");
		initialize<std::string>("predefined_guess", "", "predefined partial wave guesses of type pvxz with x=d,t,q,5,6, need to set guesstype to 'predefined'");
		initialize<int>("maxiter",10, "maximal number of iterations if no adaptive solver is used (greens function solver)");
		initialize<int>("maxiter_t",100, "maximal number of iterations in the amplitude solver");
		initialize<double>("tpno", 1.e-8, "PNO cutoff threshold");
		initialize<double>("tpno_tight",1.e-10, "PNO cutoff for the first iteration");
		initialize<bool>("canonicalize_pno",true, "canonicalize the pnos before the amplitude solver");
		initialize<double>("thresh", 1.e-4, "MRA threshold");
		initialize<double>("dconv",5.e-4, "convergence of every PNO in the Green's function solver");
		initialize<double>("econv",1.e-4, "convergence of the total energy in the Green's function solver");
		initialize<double>("op_thresh",1.e-6, "MRA operator thresh");
		initialize<std::string>("restart","none", "restart pairs of this type, use 'mp2', 'cispd', 'all' or 'none' ");
		initialize<std::string>("no_compute","none", "do not compute the pairs of this type, use 'mp2', 'cispd', 'all' or 'none' ");
		initialize<std::string>("no_opt","none", "do not optimize the pnos of this type, use 'mp2', 'cispd', 'all' or 'none' ");
		initialize<std::string>("no_guess","none", "guess for this type will be empty, use 'mp2', 'cispd', 'all' or 'none' ");
		initialize<std::string>("adaptive_solver","all", "Use adaptive solver for those pairs, use 'mp2', 'cispd', 'all' or 'none' ");
		initialize<bool>("kain", true, "use KAIN solver in amplitude solver");
		initialize<std::size_t>("kain_subspace", 5 , "subspace size of the KAIN solver");
		initialize<bool>("f12",true, "use explicit correlation");
		initialize<int> ("cispd_number", -1, "CIS(D) excitation numbers, to read in correct functions" );
		initialize<double> ("cispd_energy", 0.0, "CIS energy" );
		initialize<std::string> ("freeze_pairs", "none", "frozen pairs will not be optimized: Expected format 'a b c d' will freeze pairs ab and cd");
		initialize<std::vector<int> >("freeze_pairs_of_orbital",std::vector<int>(), " All pairs which originate from this orbital will not be optimized");
		initialize<std::vector<int> >("active_pairs_of_orbital",std::vector<int>(), " All pairs which originate from this orbital will not be frozen all other pairs will, if this vector is not empty");
		initialize<double>("econv_micro",1.e-4, "Energy convergence for microiterations in adaptive solver");
		initialize<double>("econv_macro",3.e-4, "Energy convergence for macroiterations in adaptive solver");
		initialize<int>("maxiter_micro",10, "Maximum of iterations for every cycle in the adaptive solver");
		initialize<int>("maxiter_macro",10, "Maximum of iterations of cycles in the adaptive solver");
		initialize<bool>("no_opt_in_first_iteration", false, "Do not optimize in the first iteration (then the potentials do not have to be evaluated, use this for large guesses)");
		initialize<std::string>("exchange", "full", "approximate exchange with 'neglect' or xc functional -> same syntax as moldft");
	}

	std::vector<std::pair<int,int> > freeze_pairs()const{
		const std::string str=get<std::string >("freeze_pairs");
		std::vector<std::pair<int,int> > result;
		if (str=="none"){
			return result;
		}
		std::stringstream ss(str);
		int i,j;
		while(ss>>i){
			ss>>j;
			result.push_back(std::make_pair(i,j));
		}
		return result;
	}
	std::string exchange()const {return get<std::string>("exchange");}
	bool exop_trigo()const { return get<bool >("exop_trigo");}
	int rank_increase()const { return get<int >("rank_increase");}
	int chunk()const { return get<int >("chunk");}
	std::vector<std::vector<double> > protocol()const { return get<std::vector<std::vector<double> > >("protocol");}
	bool debug()const { return get<bool >("debug");}
	std::size_t freeze()const { return get<std::size_t >("freeze");}
	int maxrank()const { return get<int >("maxrank");}
	GuessType guesstype()const { return assign_from_string<GuessType>(get<std::string >("guesstype"));}
	std::string exop()const { return get<std::string >("exop");}
	std::map<std::string, std::vector<int> >partial_wave()const {
		// return format atom-name, vector of numbers giving S, P, D, ... functions
		const std::string str=get<std::string >("partial_wave");
		// expected input format: atom-name 3s2p1d atom-name XsYpZd...
		std::stringstream ss(str);
		std::string symbol, pw_string;
		std::map<std::string, std::vector<int> > result;
		while(ss>>symbol){
			ss>>pw_string;
			std::vector<int> numbers(pw_string.size()/2);
			std::vector<char> control = {'s', 'p', 'd', 'f', 'g', 'h', 'i', 'k'};
			for (int i=0; i<pw_string.size()/2; ++i){
				char l = pw_string[2*i];
				char n = pw_string[2*i+1];
				MADNESS_ASSERT(l==control[i]);
				numbers[i]=n;
			}
			result[symbol] = numbers;
		}
		return result;
	}
	std::string predefined_guess()const{ return get<std::string >("predefined_guess");}
	int maxiter()const { return get<int >("maxiter");}
	int maxiter_t()const { return get<int >("maxiter_t");}
	double tpno()const { return get<double >("tpno");}
	double tpno_tight()const { return get<double >("tpno_tight");}
	bool canonicalize_pno()const { return get<bool >("canonicalize_pno");}
	double thresh()const { return get<double >("thresh");}
	double dconv()const { return get<double >("dconv");}
	double econv()const { return get<double >("econv");}
	double op_thresh()const { return get<double >("op_thresh");}
	PairType restart()const {return assign_from_string<PairType>(get<std::string >("restart"));}
	PairType no_compute()const { return assign_from_string<PairType>(get<std::string >("no_compute"));}
	PairType no_opt()const { return assign_from_string<PairType>(get<std::string >("no_opt"));}
	PairType no_guess()const { return assign_from_string<PairType>(get<std::string >("no_guess"));}
	bool kain()const { return get<bool >("kain");}
	std::size_t kain_subspace()const { return get<std::size_t >("kain_subspace");}
	bool f12()const { return get<bool >("f12");}
	std::vector<std::pair<int,double> > cispd()const {
		// workaround new madness data format, can only do one excitation energy each run for now
		const int number = get<int >("cispd_number");
		const double energy = get<double >("cispd_energy");
		if (number < 0){
			return std::vector<std::pair<int,double> >();
		}
		else{
			return std::vector<std::pair<int,double> >(1,std::make_pair(number, energy));
		}
	}
	std::vector<int> freeze_pairs_of_orbital()const { return get<std::vector<int> >("freeze_pairs_of_orbital");}
	std::vector<int> active_pairs_of_orbital()const { return get<std::vector<int> >("active_pairs_of_orbital");}
	PairType adaptive_solver()const { return assign_from_string<PairType>(get<std::string >("adaptive_solver"));}
	double econv_micro()const { return get<double >("econv_micro");}
	double econv_macro()const { return get<double >("econv_macro");}
	int maxiter_micro()const { return get<int >("maxiter_micro");}
	int maxiter_macro()const { return get<int >("maxiter_macro");}
	bool no_opt_in_first_iteration()const { return get<bool >("no_opt_in_first_iteration");}


};

class F12Parameters: public PNOParameters {
public:
	F12Parameters(const PNOParameters& param) : PNOParameters(param){
		initialize_f12_parameters();
	}

	F12Parameters(World& world, const std::string& inputfile, const PNOParameters& param, const std::string& TAG="pno") : PNOParameters(param){
		initialize_f12_parameters();
		QCCalculationParametersBase::read(world,inputfile,TAG);
	}


	void initialize_f12_parameters() {
		initialize<int>("ansatz",2, " (old parameter)");
		initialize<bool>("abs_c",true, " use auxilliary basis on f12Q[Kf12] part of energy (only if energytype is HYLLERAAS_ENERGYTYPE) ");
		initialize<bool>("abs_u",false, " use auxilliary basis on f12QUe part of energy (only if energytype is HYLLERAAS_ENERGYTYPE) ");
		initialize<double>("cabs_thresh",1.e-4, " thresh for auxbasis part in f12 energy ");
		initialize<bool>("external_cabs", false, " use external comp. aux. basis in addition to the pnos as auxbasis");
		initialize<EnergyType>("energytype",PROJECTED_ENERGYTYPE, " the energytype is 'projected' or 'hylleraas' functional projected energies do not need auxilliary bases for the evaluation of the f12 energy");
		initialize<double>("gamma",1.4, "The f12 length scale");
		// broken in madness' new parameter format
		//initialize<std::map<std::string,std::vector<int> > >("auxbas",std::map<std::string,std::vector<int> >() , "create LCAO auxbas, does not work in new parameter format");

	}

	bool f12()const { return get<bool >("f12");}
	int ansatz()const { return get<int >("ansatz");}
	bool abs_c()const { return get<bool >("abs_c");}
	bool abs_u()const { return get<bool >("abs_u");}
	std::size_t freeze()const { return get<std::size_t >("freeze");}
	bool debug()const { return get<bool >("debug");}
	double cabs_thresh()const { return get<double >("cabs_thresh");}
	bool external_cabs()const { return get<bool >("external_cabs");}
	EnergyType energytype()const { return get<EnergyType >("energytype");}
	double gamma()const { return get<double >("gamma");}
	std::map<std::string,std::vector<int> > auxbas()const {
		MADNESS_EXCEPTION("Issues with new parameter structure",1);
	}
	bool test6D()const {return false;}
	bool old_fQc() const {return false;}


};

} /* namespace madness */

#endif /* PNOPARAMETERS_H_ */
