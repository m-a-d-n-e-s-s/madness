#include <iostream>
#include <misc/madexcept.h>

/// \file madexcept.cc 

namespace madness {
	/// Print a MadnessException to the stream (for human consumption)
	std::ostream& operator <<(std::ostream& out, const MadnessException& e) {
	        out << "MadnessException : '";
	        if (e.msg) out << "msg=" << e.msg << " : "; 
	        if (e.assertion) out << "assertion=" << e.assertion << " : ";
	        out << "value=" << e.value << " : ";
	        if (e.line) out << "line=" << e.line << " : ";
	        if (e.function) out << "function=" << e.function << " : ";
	        if (e.filename) out << "filename='" << e.filename << "'";
	        out << std::endl;
	        
	        return out;
	}
}
