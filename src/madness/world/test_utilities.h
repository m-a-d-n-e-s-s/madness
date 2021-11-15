/*
 * test_utilities.h
 *
 *  Created on: 15 May 2019
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_TEST_UTILITIES_H_
#define SRC_APPS_CHEM_TEST_UTILITIES_H_


namespace madness {
extern void set_cout_to_ostream(std::ostream& os);
extern void set_cout_to_terminal(std::ostream& os);


/// small class for pretty printing of test output
struct test_output {
	test_output(std::string line, const bool verbose=false) : verbose(verbose) {
		int ncharacter=line.size();
		std::cout << line ;
		if (line.size()<70) std::cout << std::string(70-ncharacter, ' ' );
		logger << std::scientific << std::setprecision(8) ;
        set_cout_to_logger();
	}

    ~test_output() {
        set_cout_to_terminal();
    }


	void print_and_clear_log() {
        set_cout_to_terminal();
		std::cout << logger.str() << std::endl;
		logger.clear();
	}

	int end(bool success) {
        set_cout_to_terminal();
		if (success) std::cout << "\033[32m"   << "passed " << "\033[0m" << std::endl;
		if (not success) {
			std::cout << "\033[31m"   << "failed " << "\033[0m" << std::endl;
			print_and_clear_log();
		}
//		MADNESS_ASSERT(success);
		return (success) ? 0 : 1;
	}

    void set_cout_to_logger() {
        if (verbose) return;
        stream_buffer_cout = std::cout.rdbuf();
        std::streambuf* stream_buffer_file = logger.rdbuf();
        cout.rdbuf(stream_buffer_file);
    }

    void set_cout_to_terminal() {
        if (verbose) return;
        cout.rdbuf(stream_buffer_cout);
    }

    std::streambuf* stream_buffer_cout;
    std::stringstream logger;
    bool verbose=false;

};


}

#endif /* SRC_APPS_CHEM_TEST_UTILITIES_H_ */
