/*
 * test_utilities.h
 *
 *  Created on: 15 May 2019
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_TEST_UTILITIES_H_
#define SRC_APPS_CHEM_TEST_UTILITIES_H_


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

}

#endif /* SRC_APPS_CHEM_TEST_UTILITIES_H_ */
