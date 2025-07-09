//
// Created by Florian Bischoff on 08.07.25.
//

#ifndef RESULTS_H
#define RESULTS_H

#include <madness/external/nlohmann_json/json.hpp>
#include <madness/tensor/tensor_json.hpp>
#include<string>

//* base and derived classes for holding results of a calculation
namespace madness {

    class ResultsBase {

    public:
        ResultsBase() = default;
        virtual ~ResultsBase() = default;

        /// serialize the results to a JSON object
        virtual nlohmann::json to_json() const = 0;
    };

    /// holds metadata of the calculation

    /// create right before the calculation starts, stop() must be called after the calculation is finished
    class MetaDataResults : public ResultsBase {
    public:
        MetaDataResults(World& world) {
            time_begin=wall_time();
            mpi_size= world.size();
        }

        double time_begin=0.0;
        double time_end=0.0;
        std::string finished_at="";
        std::string git_hash="";
        long mpi_size=-1;

        void stop() {
            time_end=wall_time();
            finished_at=time_tag();
        }

        nlohmann::json to_json() const override {
            nlohmann::json j;
            // compute timing on-the-fly unless they have been set
            if (time_end==0.0) {
                j["elapsed_time"] = wall_time() - time_begin;
                j["finished_at"] = time_tag();
            } else {
                j["elapsed_time"] = time_end - time_begin;
                j["finished_at"] = finished_at;
            }
            j["git_hash"] = git_hash;
            j["host"] = std::string(HOST_SYSTEM);
            j["nthreads"] = ThreadPool::size();
            j["mpi_size"] = mpi_size;
            return j;
        }

    private:
        /// borrowed from Adrian's MolDFTLib
        std::string time_tag() const {
            auto print_time = std::chrono::system_clock::now();
            auto in_time_t = std::chrono::system_clock::to_time_t(print_time);
            std::stringstream ss;
            ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X");
            return ss.str();
        }
    };

    /// holds convergence results of the calculation
    class ConvergenceResults : public ResultsBase {
        double converged_for_thresh = 1.e10;
        double converged_for_dconv = 1.e10;

    public:
        ConvergenceResults() = default;

        /// construct from JSON
        ConvergenceResults(const nlohmann::json& j) {
            converged_for_thresh = j.value("converged_for_thresh", 1.e10);
            converged_for_dconv = j.value("converged_for_dconv", 1.e10);
        }

        /// assignment operator from JSON
        ConvergenceResults& operator=(const nlohmann::json& j) {
            converged_for_thresh = j.value("converged_for_thresh", 1.e10);
            converged_for_dconv = j.value("converged_for_dconv", 1.e10);
            return *this;
        }

        ConvergenceResults& set_converged_thresh(double thresh) {
            converged_for_thresh = thresh;
            return *this;
        }

        ConvergenceResults& set_converged_dconv(double dconv) {
            converged_for_dconv = dconv;
            return *this;
        }

        nlohmann::json to_json() const override {
            nlohmann::json j;
            j["converged_for_thresh"] = converged_for_thresh;
            j["converged_for_dconv"] = converged_for_dconv;
            return j;
        }
    };

    class PropertyResults : public ResultsBase {
    public:
        PropertyResults() = default;

        /// construct from JSON
        PropertyResults(const nlohmann::json& j) {
            energy= j.value("energy", 0.0);
            if (j.count("dipole")==1) dipole = tensor_from_json<double>(j["dipole"]);
            if (j.count("gradient")==1) gradient= tensor_from_json<double>(j["gradient"]);
        }

        double energy = 0.0;
        Tensor<double> dipole;
        Tensor<double> gradient;

        nlohmann::json to_json() const override {
            nlohmann::json j;
            j["energy"] = energy;
            j["dipole"] = tensor_to_json(dipole);
            j["gradient"] = tensor_to_json(gradient);
            return j;
        }

    };

    class SCFResults: public ResultsBase {
    public:
        Tensor<double> aeps;
        Tensor<double> beps;
        Tensor<double> afock;
        Tensor<double> bfock;

        SCFResults() = default;

        /// construct from JSON
        SCFResults(const nlohmann::json& j) {
            if (j.count("scf_eigenvalues_a")>0)
                aeps = tensor_from_json<double>(j["scf_eigenvalues_a"]);
            if (j.count("scf_eigenvalues_b") > 0)
                beps = tensor_from_json<double>(j["scf_eigenvalues_b"]);
            if (j.count("scf_fock_b") > 0)
                bfock = tensor_from_json<double>(j["scf_fock_b"]);
            if (j.count("scf_fock_b") > 0)
                bfock = tensor_from_json<double>(j["scf_fock_b"]);
        }

        nlohmann::json to_json() const override {
            nlohmann::json j;
            j["scf_eigenvalues_a"] = tensor_to_json(aeps);
            j["scf_fock_a"] = tensor_to_json(afock);
            if (beps.size() > 0) j["scf_eigenvalues_b"] = tensor_to_json(beps);
            if (bfock.size() > 0) j["scf_fock_b"] = tensor_to_json(bfock);
            return j;
        }
    };

    class CISResults: public ResultsBase {
    public:
        struct excitation_info {
            std::string irrep; // irreducible representation
            double omega; // excitation energy in Hartree
            double current_error; // error in the excitation energy
            double oscillator_strength_length; // oscillator strength
            double oscillator_strength_velocity; // oscillator strength
        };
        std::vector<excitation_info> excitations;
        long nfreeze=-1;

        CISResults() = default;

        /// construct from JSON
        CISResults(const nlohmann::json& j) {
            if (j.count("cis_excitations") > 0) {
                for (const auto& ex : j["cis_excitations"]) {
                    excitation_info ei;
                    ei.irrep = ex.value("irrep", "");
                    ei.omega = ex.value("omega", 0.0);
                    ei.current_error = ex.value("current_error", 0.0);
                    ei.oscillator_strength_length = ex.value("oscillator_strength_length", 0.0);
                    ei.oscillator_strength_velocity = ex.value("oscillator_strength_velocity", 0.0);
                    excitations.push_back(ei);
                }
            }
            nfreeze = j.value("nfreeze", -1);
        }

        nlohmann::json to_json() const override {
            nlohmann::json j;
            for (const auto& ex : excitations) {
                nlohmann::json ex_json;
                ex_json["irrep"] = ex.irrep;
                ex_json["omega"] = ex.omega;
                ex_json["current_error"] = ex.current_error;
                ex_json["oscillator_strength_length"] = ex.oscillator_strength_length;
                ex_json["oscillator_strength_velocity"] = ex.oscillator_strength_velocity;
                j["cis_excitations"].push_back(ex_json);
            }
            j["nfreeze"] = nfreeze;
            return j;
        }

    };

    class CC2Results: public CISResults {
    public:
        CC2Results() = default;

        /// construct from JSON
        CC2Results(const nlohmann::json& j) : CISResults(j) {
            mp2_correlation_energy_ = j.value("mp2_correlation_energy", 0.0);
            cc2_correlation_energy_ = j.value("cc2_correlation_energy", 0.0);
        }


        nlohmann::json to_json() const override {
            nlohmann::json j;
            j = CISResults::to_json();
            j["mp2_correlation_energy"] = mp2_correlation_energy_;
            j["cc2_correlation_energy"] = cc2_correlation_energy_;
            return j;
        }

        double mp2_correlation_energy_;
        double cc2_correlation_energy_;

    };

    class ZnemoResults: public SCFResults {
    public:

        double B=0.0; // B value for the Znemo calculation

        ZnemoResults() = default;
        /// construct from JSON
        ZnemoResults(const nlohmann::json& j) : SCFResults(j) {
            B = j.value("B", 0.0);
        }

        nlohmann::json to_json() const override {
            nlohmann::json j;
            j = SCFResults::to_json();
            j["B"] = B;
            return j;
        }

    };

    class OEPResults: public SCFResults {
    public:
        std::string model="oaep"; // model used for the OEP calculation"
        double drho=0.0; // delta rho =difference to reference (=HF?) density
        double dvir14=0.0; // diagnostic parameter
        double dvir17=0.0; // diagnostic parameter
        double x_local=0.0; // local exchange energy
        double x_hf=0.0; // HF exchange energy

        OEPResults() = default;

        /// construct from JSON
        OEPResults(const nlohmann::json& j) : SCFResults(j) {
            model = j.value("model", "oaep");
            drho = j.value("drho", 0.0);
            dvir14 = j.value("dvir14", 0.0);
            dvir17 = j.value("dvir17", 0.0);
            x_local = j.value("x_local", 0.0);
            x_hf = j.value("x_hf", 0.0);
        }

        nlohmann::json to_json() const override {
            nlohmann::json j;
            j = SCFResults::to_json();
            j["model"] = model;
            j["drho"] = drho;
            j["dvir14"] = dvir14;
            j["dvir17"] = dvir17;
            j["x_local"] = x_local;
            j["x_hf"] = x_hf;
            return j;
        }

    };

    class ResponseResults: public ResultsBase {
        virtual nlohmann::json to_json() const {
            MADNESS_EXCEPTION("to_json not implemented for ConvergenceResults", 1);
        }

    };


}
#endif //RESULTS_H
