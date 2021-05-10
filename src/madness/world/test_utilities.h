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


/// write an input file to disk and remove upon destruction

/**
 * usage: write data group mp3 with some parameters
  	std::string inputlines=R"input(mp3
			econv 1.e-4
			dconv 1.e-4
 			# econv 1.e-3
			maxiter 12# asd
			ncf (slater,1.2)
			localize no
			LocAl CanON
			end)input";
	inputfile ifile("input1",inputlines);
 */
struct test_inputfile {
	std::string fname;
	bool keepfile=false;
	test_inputfile(const std::string filename, const std::string lines) {
		fname=filename;
		std::ofstream myfile;
		myfile.open (fname);
		myfile << lines << std::endl;
		myfile.close();
	}

	~test_inputfile() {
		if (not keepfile) remove(fname.c_str());
	}
};

}

#endif /* SRC_APPS_CHEM_TEST_UTILITIES_H_ */
