#include <madness.h>
#include <madness/chem/SCF.h>
#include <madness/chem/SCFOperators.h>
#include <madness/chem/nemo.h>
#include <madness/chem/write_test_input.h>
#include <madness/misc/info.h>
#include <madness/world/test_utilities.h>
#include <ParameterManager.hpp>

class MemoryParameters : public QCCalculationParametersBase {
public:
    static constexpr char const* tag = "memory";

    MemoryParameters(World& world, const commandlineparser& parser) : MemoryParameters() {
        QCCalculationParametersBase::read_input_and_commandline_options(world, parser, "memory");
        set_derived_values();
    }

    MemoryParameters(const MemoryParameters& mp) = default;

    MemoryParameters() {
        initialize<std::string>("memory_algorithm", "default",
                                "preset for memory algorithm to use", {
                                    "default", "node_replicated_target", "small_memory", "large_memory"
                                });
        initialize<std::string>("cloud_storage", "Function",
                                "Store function or pointer to target function in exchange operator", {
                                    "Function", "FunctionPointer", "FunctionViaPointer"
                                });
        initialize<std::string>("cloud_distribution", "RankReplicated", "distribution of cloud container",
                                {"Distributed", "NodeReplicated", "RankReplicated"});
        initialize<std::string>("target_distribution", "NodeReplicated", "distribution of target functions",
                                {"Distributed", "NodeReplicated", "RankReplicated"});
    }

    std::string memory_algorithm() const { return get<std::string>("memory_algorithm"); }
    MacroTaskInfo::StoragePolicy cloud_storage() const { return MacroTaskInfo::policy_to_string(get<std::string>("cloud_storage"));}
    DistributionType cloud_distribution() const { return from_string(get<std::string>("cloud_distribution")); }
    DistributionType target_distribution() const { return from_string(get<std::string>("target_distribution")); }

    std::string get_tag() const override { return tag; }

    void set_derived_values() {
        if (memory_algorithm()=="default") {
            set_derived_value("cloud_storage",std::string("Function"));
            set_derived_value("cloud_distribution",std::string("NodeReplicated"));
            set_derived_value("target_distribution",std::string("Distributed"));
        } else if (memory_algorithm()=="node_replicated_target") {
            set_derived_value("cloud_storage",std::string("FunctionPointer"));
            set_derived_value("cloud_distribution",std::string("RankReplicated"));
            set_derived_value("target_distribution",std::string("NodeReplicated"));
        } else if (memory_algorithm()=="small_memory") {
            set_derived_value("cloud_storage",std::string("FunctionViaPointer"));
            set_derived_value("cloud_distribution",std::string("RankReplicated"));
            set_derived_value("target_distribution",std::string("Distributed"));
        } else if (memory_algorithm()=="large_memory") {
            set_derived_value("cloud_storage",std::string("Function"));
            set_derived_value("cloud_distribution",std::string("RankReplicated"));
            set_derived_value("target_distribution",std::string("Distributed"));
        }
    }
};

using namespace madness;

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    {
        startup(world, argc, argv, true);
        commandlineparser parser(argc, argv);

        srand(time(nullptr));
        if (world.rank() == 0) print(info::print_revision_information());

        if (not parser.key_exists("geometry")) parser.set_keyval("geometry", "water2");
        MemoryParameters memparam(world, parser);

        Params pm(world, parser);
        pm.get<CalculationParameters>().set_derived_value<int>("maxiter", 1);
        // memparam.print("memory","end");

        auto exchange_alg=Exchange<double,3>::string2algorithm(pm.get<CalculationParameters>().hfexalg());
        MacroTaskInfo info;
        info.storage_policy=memparam.cloud_storage();
        info.ptr_target_distribution_policy=memparam.target_distribution();
        info.cloud_distribution_policy=memparam.cloud_distribution();

        double cpu0 = cpu_time();

        SCF calc(world, pm.get<CalculationParameters>(), pm.get<Molecule>());
        if (world.rank() == 0) {
            calc.param.print("dft", "end");
            memparam.print("memory","end");
            calc.molecule.print();
        }

        // prepare orbitals
        calc.set_protocol<3>(world, 1e-4);
        MolecularEnergy me(world, calc);
        me.value(calc.molecule.get_all_coords());

        Exchange<double, 3> K = Exchange<double, 3>(world, &calc, 0);

        if (world.size() > 1) {
            LoadBalanceDeux<3> lb(world);
            for (unsigned int i = 0; i < calc.amo.size(); ++i) {
                lb.add_tree(calc.amo[i], lbcost<double, 3>(1.0, 8.0), false);
            }
            world.gop.fence();
            FunctionDefaults<3>::redistribute(world, lb.load_balance( calc.param.loadbalparts()));
            world.gop.fence();
        }

        double cpu1 = cpu_time();
        if (world.rank() == 0) printf("\ntimings for preparation   %8.2fs\n", cpu1 - cpu0);

        // compute reference number
        cpu0 = cpu1;
        K.set_algorithm(Exchange<double, 3>::large_memory);
        K.set_symmetric(true);
        K.set_printlevel(20);
        //        const vecfuncT reference = K(calc.amo);
        //        cpu1 = cpu_time();
        //        double norm = norm2(world, reference);
        //        if (world.rank() == 0)
        //            printf("timings exchange operator no multiworld largemem   %8.2fs, norm %.15e\n", cpu1 - cpu0, norm);

        // gather all results in a vector
        nlohmann::json all_results = nlohmann::json::array();


        //        for (auto exchange_alg : {
        //                 Exchange<double, 3>::multiworld_efficient_row,
        //                 Exchange<double, 3>::multiworld_efficient,
        //                 Exchange<double, 3>::large_memory,
        //                 Exchange<double, 3>::small_memory
        //             }) {
        //            for (auto info : MacroTaskInfo::get_all_presets()) {
        {
            {
                K.set_macro_task_info(info);
                if (world.rank() == 0) {
                    print("MacroTaskInfo:");
                    print(info);
                    print("Algorithm:", exchange_alg);
                }


                cpu0 = cpu_time();
                K.set_algorithm(exchange_alg);
                vecfuncT tmp = K(calc.amo);
                cpu1 = cpu_time();
                nlohmann::json result = K.statistics;
                result["time"] = cpu1 - cpu0;
                double err = 0.0;
                // err= norm2(world, reference - tmp);
                result["error"] = err;
                all_results.push_back(result);

                if (world.rank() == 0) {
                    printf("timings exchange operator                          %8.2fs, error %.2e\n", cpu1 - cpu0, err);
                    print(all_results.dump(4));
                }
            }
        }

        // print out all_results into a file
        if (world.rank() == 0) {
            std::string prefix = pm.get<CalculationParameters>().prefix();
            prefix+=".hfexalg_"+calc.param.hfexalg()
                    +".cloud_"+to_string(memparam.cloud_storage())
                    +".cloud_"+to_string(memparam.cloud_distribution())
                    +".target_"+to_string(memparam.target_distribution());
            std::ofstream file(prefix + ".benchmark.json");
            file << std::setw(4) << all_results << std::endl;
            file.close();
        }
    }

    print("\n all done -- thumbs up :-)\n");
    madness::finalize();
    return 0;
}
