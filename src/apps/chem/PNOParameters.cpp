/*
 * PNOParameters.cpp
 *
 *  Created on: Sep. 9, 2019
 *      Author: jsk
 */

#include <PNOParameters.h>

namespace madness {

std::ostream& operator << (std::ostream& os, const PairType& en){
 switch(en){
 case(UNKNOWN_PAIRTYPE):{ os << "unknown"; break;}
 case(MP2_PAIRTYPE):{ os << "mp2"; break;}
 case(CISPD_PAIRTYPE):{ os << "cispd"; break;}
 case(ALL_PAIRTYPE):{ os << "all"; break;}
 case(NONE_PAIRTYPE):{ os << "none"; break;}
 default :{ os << "undefined enum in PairType"; break;}
}
return os;
}

std::istream& operator >> (std::istream& is, PairType& en){
        std::string s;
        is>>s;
        std::transform(s.begin(), s.end(), s.begin(), ::toupper);
        if(s=="UNKNOWN_PAIRTYPE") en=UNKNOWN_PAIRTYPE;
        else if(s=="MP2") en=MP2_PAIRTYPE;
        else if(s=="CISPD") en=CISPD_PAIRTYPE;
        else if(s=="ALL") en=ALL_PAIRTYPE;
        else if(s=="NONE") en=NONE_PAIRTYPE;
        else MADNESS_EXCEPTION(("unknown string to PairType conversion for s="+s).c_str(),1);
        return is;
}
std::ostream& operator << (std::ostream& os, const EnergyType& en){
         switch(en){
         case(UNKNOWN_ENERGYTYPE):{ os << "unknown"; break;}
         case(PROJECTED_ENERGYTYPE):{ os << "projected"; break;}
         case(HYLLERAAS_ENERGYTYPE):{ os << "hylleraas"; break;}
         default :{ os << "undefined enum in EnergyType"; break;}
        }
        return os;
}
std::istream& operator >> (std::istream& is, EnergyType& en){
        std::string s;
        is>>s;
        std::transform(s.begin(), s.end(), s.begin(), ::toupper);
        if(s=="UNKNOWN_ENERGYTYPE") en=UNKNOWN_ENERGYTYPE;
        else if(s=="PROJECTED") en=PROJECTED_ENERGYTYPE;
        else if(s=="HYLLERAAS") en=HYLLERAAS_ENERGYTYPE;
        else MADNESS_EXCEPTION(("unknown string to EnergyType conversion for s="+s).c_str(),1);
        return is;
}
std::ostream& operator << (std::ostream& os, const GuessType& en){
         switch(en){
         case(UNKNOWN_GUESSTYPE):{ os << "unknown"; break;}
         case(PARTIAL_WAVE_GUESSTYPE):{ os << "partial_wave"; break;}
         case(FROM_FILE_GUESSTYPE):{ os << "from_file"; break;}
         case(PREDEFINED_GUESSTYPE):{ os << "predefined"; break;}
         case(SCF_GUESSTYPE):{ os << "scf"; break;}
         case(EXOP_GUESSTYPE):{ os << "exop"; break;}
         case(EMPTY_GUESSTYPE):{ os << "empty"; break;}
         case(PSI4_GUESSTYPE):{ os << "psi4"; break;}
         default :{ os << "undefined enum in GuessType"; break;}
        }
        return os;
}

std::istream& operator >> (std::istream& is, GuessType& en){
	std::string s;
	is>>s;
	std::transform(s.begin(), s.end(), s.begin(), ::toupper);
	if (s=="UNKNOWN") en = UNKNOWN_GUESSTYPE;
	else if (s=="PARTIAL_WAVE") en= PARTIAL_WAVE_GUESSTYPE;
	else if (s=="FROM_FILE") en= FROM_FILE_GUESSTYPE;
	else if (s=="PREDEFINED") en= PREDEFINED_GUESSTYPE;
	else if (s=="SCF") en= SCF_GUESSTYPE;
	else if (s=="PARTIAL_WAVE") en= PARTIAL_WAVE_GUESSTYPE;
	else if (s=="EXOP") en= EXOP_GUESSTYPE;
	else if (s=="EMPTY") en= EMPTY_GUESSTYPE;
	else if (s=="PSI4") en= PSI4_GUESSTYPE;
	else MADNESS_EXCEPTION(("unknown string to GuessType conversion for s="+s).c_str(),1);
	return is;

}

} /* namespace madness */
