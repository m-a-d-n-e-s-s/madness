//
// Created by Florian Bischoff on 5/27/21.
// Modified by Adrian Hurtado on 2/01/22
//

#ifndef MADNESS_WRITE_RESPONSE_INPUT_H
#define MADNESS_WRITE_RESPONSE_INPUT_H

#include <madness/chem/CalculationParameters.h>

#include "response_parameters.h"

namespace molresponse {

    /// will write a test input and remove it from disk upon destruction
    struct write_test_input {

        double eprec = 1.e-4;// was 1e-4 ... trying to make test faster

        std::string filename_;
        std::string molecule_path;
        bool keepfile = true;

        write_test_input() : filename_("moldft.in") {}

        explicit write_test_input(const CalculationParameters &param,  std::string filename,  std::string mol_path)
            : filename_(std::move(filename)), molecule_path(std::move(mol_path)) {
            std::ofstream of(filename_);
            write_to_test_input("dft", &param, of);
            write_molecule_to_test_input(molecule_path, of);
            of.close();
        }

        ~write_test_input() {
            if (not keepfile) std::remove(filename_.c_str());
        }

        std::string filename() const { return filename_; }

        static std::ostream &write_to_test_input(const std::string &groupname, const QCCalculationParametersBase *param,
                                                 std::ostream &of) {
            of << groupname << endl;
            of << param->print_to_string(true);
            of << "end\n";
            return of;
        }

        static std::ostream &write_molecule_to_test_input(std::string mol_path, std::ostream &of) {


            std::cout << mol_path << "\n";
            std::ifstream mol_file(mol_path);
            std::string line;
            std::string no_orient_line;
            int i = 0;
            while (getline(mol_file, line)) {
                if (i == 1) {
                    no_orient_line = "no_orient true";
                    std::cout << no_orient_line << "\n";
                    of << no_orient_line << "\n";
                }
                std::cout << line << "\n";
                of << line << "\n";
                i++;
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

        explicit write_response_input(const ResponseParameters &param, const std::string &filename)
            : filename_(filename) {
            std::ofstream of(filename_);
            write_to_test_input("response", &param, of);
            of.close();
        }

        ~write_response_input() {
            if (not keepfile) std::remove(filename_.c_str());
        }

        std::string filename() const { return filename_; }

        static std::ostream &write_to_test_input(const std::string &groupname, const ResponseParameters *param,
                                                 std::ostream &of) {
            of << groupname << endl;
            of << param->print_to_string(true);
            of << "end\n";
            return of;
        }
    };


}// namespace molresponse


#endif//MADNESS_WRITE_RESPONSE_INPUT_H
