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
#include<madness/chem/CalculationParameters.h>

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

class PNOParameters: public QCCalculationParametersBase {
public:


	std::string get_tag() const override {
		return std::string("pno");
	}

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

	PNOParameters(World& world, const commandlineparser& parser, const std::string& TAG="pno") : QCCalculationParametersBase(){
		initialize_pno_parameters();
		QCCalculationParametersBase::read_input_and_commandline_options(world,parser,TAG);
	}

	PNOParameters(World& world, const commandlineparser& parser, const Molecule& molecule, const std::string& TAG="pno") : QCCalculationParametersBase(){
		initialize_pno_parameters();
		QCCalculationParametersBase::read_input_and_commandline_options(world,parser,TAG);
		set_derived_values(molecule);
	}

	void initialize_pno_parameters() {
		initialize<int>("rank_increase", 15 , "maximum rank to increase in every macroiteration");
		initialize<int>("chunk", 100 , "chunk of functions operated on in parallel when G or K is applied (prevent memory shortage)");
		initialize<bool>("debug",false, "debug mode");
		initialize<std::size_t>("freeze",0, "frozen core approximation");
		initialize<int>("maxrank", 999, "maximal pno rank for all pairs");
		initialize<std::string>("guesstype","exop", "the guesstype: exop (recommended, multiply polynomial excitation operators onto the occupied orbitals), empty (don't compute a guess), scf (use the ao basis from the SCF solver), partial_wave (create s,p,d ..., type functions and place on atoms), predefined (same as partial_wave but for predefined sets) ");
		initialize<std::string>("exop", "multipole", "this string defines which excitation operators can be used, use 'dipole', 'dipole+', 'quadrupole', 'ocopole', or 'multipole', need to set guesstype to 'exop' which is the default ");
		initialize<bool>("exop_trigo", true, "use trigonometric excitation operators ( x --> sin(x) ");
		initialize<std::string>("partial_wave", "", "atom centered partial wave guess of format like 'h-2s1p-o-3s2p1d', need to set guesstype to 'partial_wave' and adaptive_solver to 'none' ");
		initialize<std::string>("predefined", "", "predefined partial wave guesses of type pvxz with x=d,t,q,5,6, need to set guesstype to 'predefined' and adaptive_solver to 'none'");
		initialize<int>("maxiter_t",100, "maximal number of iterations in the amplitude solver");
		initialize<int>("maxiter_micro",10, "Maximum of iterations for every cycle in the adaptive solver");
		initialize<int>("maxiter_macro",10, "Maximum of iterations of cycles in the adaptive solver");
		initialize<double>("tpno", 1.e-8, "PNO cutoff threshold");
		initialize<double>("tpno_tight",1.e-10, "PNO cutoff for the first iteration");
		initialize<bool>("canonicalize_pno",true, "canonicalize the pnos before the amplitude solver");
		initialize<double>("thresh", 1.e-3, "MRA threshold");
		initialize<double>("econv_micro",1.e-3, "Energy convergence for microiterations (Greens function based optimization) in adaptive solver");
		initialize<double>("econv_pairs",1.e-4, "ENergy convergence for individual pairs in adaptive solver. Converged pairs are frozen automatically");
		initialize<double>("econv_macro",1.e-3, "Energy convergence for macroiterations in adaptive solver, no effect if adaptive_solver is deactivated");
		initialize<double>("dconv",1.e-1, "convergence of every PNO in the Green's function solver");
		initialize<double>("op_thresh",1.e-6, "MRA operator thresh");
		initialize<std::string>("restart","none", "restart pairs of this type, use 'mp2', 'cispd', 'all' or 'none' ");
		initialize<std::string>("no_compute","none", "do not compute the pairs of this type, use 'mp2', 'cispd', 'all' or 'none' ");
		initialize<std::string>("no_opt","none", "do not optimize the pnos of this type, use 'mp2', 'cispd', 'all' or 'none' ");
		initialize<std::string>("no_guess","none", "guess for this type will be empty, use 'mp2', 'cispd', 'all' or 'none' ");
		initialize<std::string>("adaptive_solver","all", "Use adaptive solver for those pairs, use 'mp2', 'cispd', 'all' or 'none', works only in combination with guesstype 'exop' ");
		initialize<bool>("kain", true, "use KAIN solver in amplitude solver");
		initialize<std::size_t>("kain_subspace", 5 , "subspace size of the KAIN solver (amplitudes)");
		initialize<bool>("f12",true, "use explicit correlation");
		initialize<int> ("cispd_number", -1, "CIS(D) excitation numbers, to read in correct functions" );
		initialize<double> ("cispd_energy", 0.0, "CIS energy" );
		initialize<std::string> ("freeze_pairs", "none", "frozen pairs will not be optimized: Expected format 'a b c d' will freeze pairs ab and cd");
		initialize<std::vector<int> >("freeze_pairs_of_orbital",std::vector<int>(), " All pairs which originate from this orbital will not be optimized");
		initialize<std::vector<int> >("active_pairs_of_orbital",std::vector<int>(), " All pairs which originate from this orbital will not be frozen all other pairs will, if this vector is not empty");
		initialize<bool>("no_opt_in_first_iteration", false, "Do not optimize in the first iteration (then the potentials do not have to be evaluated, use this for large guesses)");
		initialize<std::string>("exchange", "full", "approximate exchange with 'neglect' or xc functional -> same syntax as moldft");
		initialize<bool>("save_pnos",true, "Save the OBS-PNOs to a file, before and after orthonormalization.");
		initialize<bool>("diagonal", false, "Compute only diagonal PNOs");
	}

	void set_derived_values(const Molecule& molecule) {

		// auto determine freeze parameter for the first two shells
		// deactivate by setting freeze 0 or any other value in the input file
		size_t freeze = 0;
		for(const Atom& atom: molecule.get_atoms()){
			if (atom.atomic_number < 3){
				// no frozen core for H and He
			}else if(atom.atomic_number < 11){
				freeze += 1;
			}else if(atom.atomic_number < 19){
				freeze += 5;
			}
			// beyond this point the parameter is better set manually
		}
		set_derived_value("freeze", freeze);

		set_derived_value("no_guess", get<std::string >("no_compute"));
		set_derived_value("restart", get<std::string>("no_compute"));
		set_derived_value("tpno_tight", 0.01*tpno());

		// set default values for adaptive solver
		if(adaptive_solver()){
			const std::string gt = "exop";
			const std::string ex = "multipole";
			set_derived_value("guesstype", gt);
			set_derived_value("exop", ex);
			set_derived_value("econv_macro", thresh());
			set_derived_value("econv_micro", thresh());
			set_derived_value("econv_pairs", 0.1*thresh());
		}

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
	bool diagonal()const {return get<bool>("diagonal");}
	bool save_pnos()const { return get<bool >("save_pnos");}
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
	std::map<std::string, std::vector<int> >partial_wave(const std::string& key = "partial_wave")const {
		// return format atom-name, vector of numbers giving S, P, D, ... functions
		std::string str=get<std::string >(key);
		std::transform(str.begin(), str.end(), str.begin(), ::tolower);
		// madness parameter format does not allow blancs so we use '-' and transform them here
		std::replace(str.begin(), str.end(), '-', ' ');
		// expected input format: atom-name-3s2p1d atom-name-XsYpZd...
		std::stringstream ss(str);
		std::string symbol, pw_string;
		std::map<std::string, std::vector<int> > result;
		while(ss>>symbol){
			ss>>pw_string;
			std::vector<int> numbers(pw_string.size()/2);
			std::vector<char> control = {'s', 'p', 'd', 'f', 'g', 'h', 'i', 'k'};
			for (size_t i=0; i<pw_string.size()/2; ++i){
				char l = pw_string[2*i+1];
				char n = pw_string[2*i];
				MADNESS_ASSERT(l==control[i]);
				numbers[i]=(n-'0'); // need to convert ascii to int
			}
			result[symbol] = numbers;
		}
		return result;
	}
	std::string predefined_guess()const{ return get<std::string >("predefined");}
	int maxiter()const { return maxiter_micro();}
	int maxiter_t()const { return get<int >("maxiter_t");}
	double tpno()const { return get<double >("tpno");}
	double tpno_tight()const { return get<double >("tpno_tight");}
	bool canonicalize_pno()const { return get<bool >("canonicalize_pno");}
	double thresh()const { return get<double >("thresh");}
	double dconv()const { return get<double >("dconv");}
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
	double econv()const { return econv_micro();}
	double econv_pairs()const { return get<double >("econv_pairs");}
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

	F12Parameters(World& world, const commandlineparser& parser, const PNOParameters& param, const std::string& TAG="pno") : PNOParameters(param){
		initialize_f12_parameters();
		QCCalculationParametersBase::read_input_and_commandline_options(world,parser,TAG);
	}

	std::string get_tag() const override {
		return std::string("f12");
	}


	void initialize_f12_parameters() {
		initialize<bool>("abs_c",true, " use auxilliary basis on f12Q[Kf12] part of energy (only if energytype is HYLLERAAS_ENERGYTYPE). If switched off the part neglected!");
		initialize<bool>("abs_u",false, " use auxilliary basis on f12QUe part of energy (only if energytype is HYLLERAAS_ENERGYTYPE). If switched off the part is computed in full (recommended) ");
		initialize<double>("cabs_thresh",1.e-4, " thresh for auxbasis part in f12 energy ");
		initialize<std::string>("energytype", "projected", " the energytype is 'projected' or 'hylleraas' functional projected energies do not need auxilliary bases for the evaluation of the f12 energy. It's recommended to use projected_energies! For Hylleraas type you need to specify an auxbas from file OR internal");
		initialize<double>("gamma",1.4, "The f12 length scale");
		initialize<std::string>("auxbas", "none", "atom centered partial wave guess of format like 'h-2s1p-o-3s2p1d' ");
		initialize<std::string>("auxbas_file", "none", " use external comp. aux. basis in addition to the pnos as auxbasis. Give the filename as parameter. Give the auxbas in turbomole format. Don't use contractions. If a file is specified the auxbas parameter has no effect");
	}

	bool f12()const { return get<bool >("f12");}
	bool abs_c()const { return get<bool >("abs_c");}
	bool abs_u()const { return get<bool >("abs_u");}
	double cabs_thresh()const { return get<double >("cabs_thresh");}
	std::string auxbas_file()const {
		return get<std::string >("auxbas_file");
	}
	EnergyType energytype()const {
		std::string key = get<std::string>("energytype");
		std::transform(key.begin(), key.end(), key.begin(), ::tolower);
		std::stringstream ss(key);
		EnergyType result = UNKNOWN_ENERGYTYPE;
		ss >> result;
		MADNESS_ASSERT(result != UNKNOWN_ENERGYTYPE);
		return result;
	}
	double gamma()const { return get<double >("gamma");}
	std::map<std::string,std::vector<int> > auxbas()const {
		return partial_wave("auxbas");
	}

	void set_derived_values() {
		if(energytype() == HYLLERAAS_ENERGYTYPE) set_derived_value("abs_c", true);
		if(energytype() == HYLLERAAS_ENERGYTYPE) set_derived_value("abs_u", false);
		set_derived_value("cabs_thresh", thresh());
	}

};

} /* namespace madness */

#endif /* PNOPARAMETERS_H_ */
