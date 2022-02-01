//
// Created by Florian Bischoff on 5/27/21.
// Modified by Adrian Hurtado on 2/01/22
//

#ifndef MADNESS_WRITE_TEST_INPUT_H
#define MADNESS_WRITE_TEST_INPUT_H

#include <chem/CalculationParameters.h>

#include "response_parameters.h"

namespace madness {

    /// will write a test input and remove it from disk upon destruction
    struct write_test_input {

        double eprec = 1.e-4;// was 1e-4 ... trying to make test faster

        std::string filename_;
        bool keepfile = true;

        write_test_input() : filename_("moldft.in") {}

        explicit write_test_input(const CalculationParameters& param, const std::string& filename,
                                  const std::string& mol = "lih")
            : filename_(filename) {
            std::ofstream of(filename_);
            write_to_test_input("dft", &param, of);
            write_molecule_to_test_input(mol, of);
            of.close();
        }

        ~write_test_input() {
            if (not keepfile) std::remove(filename_.c_str());
        }

        std::string filename() const { return filename_; }

        static std::ostream& write_to_test_input(const std::string& groupname,
                                                 const QCCalculationParametersBase* param,
                                                 std::ostream& of) {
            of << groupname << endl;
            of << param->print_to_string(true);
            of << "end\n";
            return of;
        }

        static std::ostream& write_molecule_to_test_input(std::string mol, std::ostream& of) {
            if (mol == "lih") {
                of << "geometry\n";
                of << "Li 0.0    0.0 0.0\n";
                of << "H  1.4375 0.0 0.0\n";
                of << "end\n";
            } else if (mol == "hf") {
                //double eprec=1.e-5; // trying to make test faster
                of << "geometry\n";
                of << "F  0.1    0.0 0.2\n";
                of << "H  1.4375 0.0 0.0\n";
                of << "end\n";
            }
            return of;
        }
    };
    /// will write a test input and remove it from disk upon destruction
    struct write_response_input {

        double eprec = 1.e-4;// was 1e-4 ... trying to make test faster

        std::string filename_;
        bool keepfile = true;

        write_response_input() : filename_("response_test") {}

        explicit write_response_input(const ResponseParameters& param, const std::string& filename)
            : filename_(filename) {
            std::ofstream of(filename_);
            write_to_test_input("response", &param, of);
            of.close();
        }

        ~write_response_input() {
            if (not keepfile) std::remove(filename_.c_str());
        }

        std::string filename() const { return filename_; }

        static std::ostream& write_to_test_input(const std::string& groupname,
                                                 const ResponseParameters* param,
                                                 std::ostream& of) {
            of << groupname << endl;
            of << param->print_to_string(true);
            of << "end\n";
            return of;
        }
    };


}// namespace madness


#endif//MADNESS_WRITE_TEST_INPUT_H
