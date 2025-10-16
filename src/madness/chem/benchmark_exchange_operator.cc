#include <madness.h>
#include <madness/chem/SCF.h>
#include <madness/chem/SCFOperators.h>
#include <madness/chem/nemo.h>
#include <madness/chem/write_test_input.h>
#include <madness/misc/info.h>
#include <madness/world/test_utilities.h>
#include <ParameterManager.hpp>


using namespace madness;

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    {
        startup(world, argc, argv, true);
        commandlineparser parser(argc, argv);

        srand(time(nullptr));
        if (world.rank() == 0) print(info::print_revision_information());

        if (not parser.key_exists("geometry")) parser.set_keyval("geometry", "water2");

        Exchange<double,3>::Algorithm alg;
        if (not parser.key_exists("algorithm")) {
            print("no algorithm specified, using multiworld_efficient_row");
            parser.set_keyval("algorithm", "multiworld_efficient_row");
        } else {
            if (parser.value("algorithm")=="small") alg=Exchange<double,3>::small_memory;
            else if (parser.value("algorithm")=="large") alg=Exchange<double,3>::large_memory;
            else if (parser.value("algorithm")=="multiworld_efficient") alg=Exchange<double,3>::multiworld_efficient;
            else if (parser.value("algorithm")=="multiworld_efficient_row") alg=Exchange<double,3>::multiworld_efficient_row;
            else {
                MADNESS_EXCEPTION("no valid algorithm specified, use small, large, multiworld_efficient, multiworld_efficient_row",1);
            }
        }
        MacroTaskInfo info=MacroTaskInfo::preset("default");
        if (not parser.key_exists("cloud_preset")) {
            print("no cloud_preset specified, using default");
            parser.set_keyval("cloud_preset", "default");
        } else {
            if (parser.value("cloud_preset")=="default") info=MacroTaskInfo::preset("default");
            else if (parser.value("cloud_preset")=="small_memory") info=MacroTaskInfo::preset("small_memory");
            else if (parser.value("cloud_preset")=="large_memory") info=MacroTaskInfo::preset("large_memory");
            else if (parser.value("cloud_preset")=="node_replicated_target") info=MacroTaskInfo::preset("node_replicated_target");
            else {
                MADNESS_EXCEPTION("no valid cloud preset specified, use default, small_memory, large_memory",1);
            }
        }

        Params pm(world, parser);
        pm.get<CalculationParameters>().set_derived_value<int>("maxiter",1);

        double cpu0 = cpu_time();

        SCF calc(world, pm.get<CalculationParameters>(), pm.get<Molecule>());
        if (world.rank() == 0) {
            calc.param.print("", "");
            calc.molecule.print();
        }

        calc.set_protocol<3>(world, 1.e-4);
        MolecularEnergy me(world, calc);
        me.value(calc.molecule.get_all_coords());
        Exchange<double, 3> K = Exchange<double, 3>(world, &calc, 0);

        if (world.size() > 1) {
            LoadBalanceDeux<3> lb(world);
            for (unsigned int i = 0; i < calc.amo.size(); ++i) {
                lb.add_tree(calc.amo[i], lbcost<double, 3>(1.0, 8.0), false);
            }
            world.gop.fence();
            FunctionDefaults<3>::redistribute(world, lb.load_balance(
                                                  calc.param.loadbalparts()));
            // 6.0 needs retuning after param.vnucextra
            world.gop.fence();
        }

        vecfuncT tmp, tmp1;

        double cpu1 = cpu_time();
        if (world.rank() == 0) printf("\ntimings for preparation   %8.2fs\n", cpu1 - cpu0);

        auto amo=calc.amo;
        print_size(world,amo,"amo");

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


//        for (auto alg : {
//                 Exchange<double, 3>::multiworld_efficient_row,
//                 Exchange<double, 3>::multiworld_efficient,
//                 Exchange<double, 3>::large_memory,
//                 Exchange<double, 3>::small_memory
//             }) {
//            for (auto info : MacroTaskInfo::get_all_presets()) {
        {
            {
                K.set_macro_task_info(info);
                if (world.rank()==0) {
                    print("MacroTaskInfo:");
                    print(info);
                    print("Algorithm:", alg);
                }


                cpu0 = cpu_time();
                K.set_algorithm(alg);
                tmp = K(calc.amo);
                cpu1 = cpu_time();
                nlohmann::json result=K.statistics;
                result["time"] = cpu1 - cpu0;
                double err=0.0;
                // err= norm2(world, reference - tmp);
                result["error"] = err;
                all_results.push_back(result);

                if (world.rank() == 0) {
                    printf("timings exchange operator                          %8.2fs, error %.2e\n", cpu1 - cpu0, err);
                    print(K.statistics.dump(4));
                }
            }
        }

        // print out all_results into a file
        if (world.rank()==0) {
            std::string prefix=pm.get<CalculationParameters>().prefix();
            std::ofstream file(prefix+".benchmark_exchange_operator_results.json");
            file << std::setw(4) << all_results << std::endl;
            file.close();
        }
    }

    print("\n all done -- thumbs up :-)\n");
    madness::finalize();
    return 0;
}
