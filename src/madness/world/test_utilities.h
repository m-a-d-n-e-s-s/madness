/*
 * test_utilities.h
 *
 *  Created on: 15 May 2019
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_TEST_UTILITIES_H_
#define SRC_APPS_CHEM_TEST_UTILITIES_H_

#include<chem/CalculationParameters.h>

namespace madness {

/// small class for pretty printing of test output
struct test_output {
	test_output(std::string line) {
		int ncharacter=line.size();
		std::cout << line ;
		if (line.size()<70) std::cout << std::string(70-ncharacter, ' ' );
		logger << std::scientific << std::setprecision(8) ;
	}

	std::stringstream logger;

	void print_and_clear_log() {
		std::cout << logger.str() << std::endl;
		logger.clear();
	}

	int end(bool success) {
		if (success) std::cout << "\033[32m"   << "passed " << "\033[0m" << std::endl;
		if (not success) {
			std::cout << "\033[31m"   << "failed " << "\033[0m" << std::endl;
			print_and_clear_log();
		}
//		MADNESS_ASSERT(success);
		return (success) ? 0 : 1;
	}
};


    /// will write a test input and remove it from disk upon destruction
    struct write_test_input {

        double eprec=1.e-3; // was 1e-4 ... trying to make test faster

        std::string filename_;
        write_test_input() : filename_("testinput") {}

        write_test_input(const CalculationParameters& param, const std::string& mol="lih") : filename_("test_MO_input") {
            std::ofstream of(filename_);
            write_to_test_input("dft",&param,of);
            write_molecule_to_test_input(mol,of);
            of.close();
        }

        ~write_test_input() {
            std::remove(filename_.c_str());
        }

        std::string filename() const {return filename_;}

        static std::ostream& write_to_test_input(const std::string groupname, const QCCalculationParametersBase* param, std::ostream& of) {
            of << groupname << endl;
            of << param->print_to_string(true);
            of << "end\n";
            return of;
        }

        static std::ostream& write_molecule_to_test_input(std::string mol, std::ostream& of) {
            if (mol=="lih") {
                of << "geometry\n";
                of << "Li 0.0    0.0 0.0\n";
                of << "H  1.4375 0.0 0.0\n";
                of << "end\n";
            } else if (mol=="hf") {
                //double eprec=1.e-5; // trying to make test faster
                of << "geometry\n";
                of << "F  0.1    0.0 0.2\n";
                of << "H  1.4375 0.0 0.0\n";
                of << "end\n";
            }
            return of;
        }

    };


}

#endif /* SRC_APPS_CHEM_TEST_UTILITIES_H_ */
