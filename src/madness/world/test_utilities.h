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
    /// @param[in]  use as if (world.rank()==0) to avoid printing on all ranks
	test_output(std::string line, bool print=true) {
        if (print) std::cout << ltrim_to_length(line,70);
		logger << std::scientific << std::setprecision(8) ;
        time_begin=cpu_time();
        time_last_checkpoint=time_begin;
        set_cout_to_logger();
	}

    static std::string ltrim_to_length(std::string line, long length=70) {
        int ncharacter=line.size();
        if (line.size()<size_t(length)) line+= std::string(length-ncharacter, ' ' );
        return line;
    }

    ~test_output() {
        set_cout_to_terminal(false);
    }

	void print_and_clear_log() {
        set_cout_to_terminal();
		std::cout << logger.str() << std::endl;
		logger.clear();
	}

    void checkpoint(double error, double tol,
                    std::string message, double time=-1.0) {
        bool use_logger=cout_set_to_logger;
        set_cout_to_terminal(false);
        bool success=error<tol;
        final_success = success and final_success;
        if (not have_checkpoints) print("");    // first checkpoint
        have_checkpoints=true;
        std::cout << "  " << ltrim_to_length(message,66);
        double time1=cpu_time()-time_last_checkpoint;
        time_last_checkpoint=cpu_time();
        print_success_fail(std::cout,success,time1,error);
        if (not success) {
            print_and_clear_log();
        }
        if (use_logger) set_cout_to_logger();
    }

    void checkpoint(double value, double reference, double tol,
                    std::string message, double time=-1.0) {
        bool use_logger=cout_set_to_logger;
        set_cout_to_terminal(false);
        double error=fabs(value-reference);
        bool success=error<tol;
        final_success = success and final_success;
        if (not have_checkpoints) print("");    // first checkpoint
        have_checkpoints=true;
        std::cout << "  " << ltrim_to_length(message,66);
        double time1=cpu_time()-time_last_checkpoint;
        time_last_checkpoint=cpu_time();
        print_success_fail(std::cout,success,time1,error);
        if (not success) {
            print_and_clear_log();
        }
        if (use_logger) set_cout_to_logger();
    }

    void checkpoint(bool success, std::string message, double time=-1.0) {
        bool use_logger=cout_set_to_logger;
        set_cout_to_terminal(false);
        final_success = success and final_success;
        if (not have_checkpoints) print("");    // first checkpoint
        have_checkpoints=true;
        std::cout << "  " << ltrim_to_length(message,66);
        double time1=cpu_time()-time_last_checkpoint;
        time_last_checkpoint=cpu_time();
        print_success_fail(std::cout,success,time1,-1.0);
        if (not success) {
            print_and_clear_log();
        }
        if (use_logger) set_cout_to_logger();
    }

    void print_success_fail(std::ostream& os, bool success, double time, double error) {

        if (success) os << "\033[32m"   << "passed " << "\033[0m";
        else os << "\033[31m"   << "failed " << "\033[0m";
        if (time>0) {
            std::stringstream ss;
            ss<< " in " << std::fixed << std::setprecision(1) << time << "s";
            os << ss.str();
        }
        if (error>=0.0) os << " error " << error;
        os << std::endl;
    }

	int end(bool success=true) {
        set_cout_to_terminal(false);
        if (have_checkpoints) std::cout << ltrim_to_length("--> final result -->",70);
        success = success and final_success;
        double time_end=cpu_time();
        print_success_fail(std::cout,success,time_end-time_begin,-1.0);
        if (not success) print_and_clear_log();
		return (success) ? 0 : 1;
	}

    void set_cout_to_logger() {
        if (cout_set_to_logger) return;
        cout_set_to_logger=true;
        stream_buffer_cout = std::cout.rdbuf();
        std::streambuf* stream_buffer_file = logger.rdbuf();
        std::cout.rdbuf(stream_buffer_file);
    }

    /// newline for use by user, not for internal use (e.g. checkpoint())
    void set_cout_to_terminal(bool newline=true) {
        if (not cout_set_to_logger) return;
        if (cout_set_to_logger) {
            std::cout.rdbuf(stream_buffer_cout);
        }
        cout_set_to_logger=false;
        if (newline) std::cout << std::endl;
    }

    std::stringstream logger;
    bool get_final_success() const {return final_success;}
private:

    bool final_success=true;
    bool cout_set_to_logger=false;          // do not change this directly!
    bool have_checkpoints=false;
    std::streambuf* stream_buffer_cout;
    double time_begin=0.0;
    double time_last_checkpoint=0.0;
};


}

#endif /* SRC_APPS_CHEM_TEST_UTILITIES_H_ */
