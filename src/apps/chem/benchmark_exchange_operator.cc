
//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness.h>
#include <chem/SCFOperators.h>
#include <chem/SCF.h>
#include <chem/nemo.h>
#include <madness/world/test_utilities.h>
#include <chem/write_test_input.h>
#include <madness/misc/gitinfo.h>

using namespace madness;

int main(int argc, char** argv) {
    commandlineparser parser(argc,argv);
    madness::initialize(argc, argv);
    {
        madness::World world(SafeMPI::COMM_WORLD);
        {

            world.gop.fence();
            startup(world, argc, argv, true);
            srand(time(nullptr));
            if (world.rank() == 0) print(info::print_revision_information());

            commandlineparser parser(argc, argv);
            if (not parser.key_exists("structure")) parser.set_keyval("structure","water2");
            CalculationParameters param;
            param.set_user_defined_value("econv",1.e-6);
            std::vector<double> protocol;   // empty protocol
            param.set_user_defined_value("protocol",protocol);
            param.set_user_defined_value("molecular_structure",parser.value("structure"));

            write_test_input inputfile(param);
            inputfile.keepfile=true;
            parser.set_keyval("input",inputfile.filename());
            world.gop.fence();

            double cpu0 = cpu_time();

            SCF calc(world, parser);
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
                        calc.param.loadbalparts())); // 6.0 needs retuning after param.vnucextra
                world.gop.fence();
            }

            vecfuncT tmp, tmp1;
            double err;

            double cpu1 = cpu_time();
            if (world.rank() == 0) printf("\ntimings for preparation   %8.2fs\n", cpu1 - cpu0);

            /*
             * serialize pointer to funcimpl
             * make a universe function in subworld by function::set_impl
             * const_cast<dcT&>.send(key,&FunctionNode<T,NDIM>::gaxpy, gaxpy_args)
             *
            */

//    // split world into subworlds
//    std::shared_ptr<World> subworld_ptr=MacroTaskQ::create_worlds(world,world.size());
//    World& subworld=*subworld_ptr;
//
//    print("world.size()   ",world.size(), world.id());
//    print("subworld.size()",subworld.size(), subworld.id());
//
//    // input function: a gaussian
//    FunctionDefaults<3>::set_default_pmap(subworld);
//    real_function_3d source=real_factory_3d(subworld).functor([](const coord_3d& r){return exp(-3.0*r.normf()*r.normf());});
//    print("source function lives in world",source.get_impl()->get_world().id());
//
//    // output function: an empty function
//    FunctionDefaults<3>::set_default_pmap(world);
//    real_function_3d result=real_factory_3d(world).functor([](const coord_3d& r){return exp(-3.0*r.normf()*r.normf());});
//
//    Cloud cloud(world);
//    cloud.store(world,result.get_impl(),1);
//
//    struct transfer_op {
//        Function<double,3> result;
//        World& subworld;
//        transfer_op(Function<double,3>& r, World& subworld) : subworld(subworld) {
//            result.set_impl(r.get_impl());
//        }
//        transfer_op(Cloud& cloud, World& subworld) : subworld(subworld) {
//            std::shared_ptr<FunctionImpl<double,3> > rimpl;
//            cloud.load(subworld,rimpl,1);
//            result.set_impl(rimpl);
//        }
//
//        void operator()(const Key<3>& key, FunctionNode<double,3>& node) const {
//            FunctionImpl<double,3>::dcT&  coeffs=const_cast<FunctionImpl<double,3>::dcT&>(result.get_impl()->get_coeffs());
//            coeffs.send(key, &FunctionNode<double,3>:: template gaxpy_inplace<double,double>, 1.0, node, 1.0);
//
//        }
//    };
//    source.unaryop_node(transfer_op(result,subworld));
//    print("result function lives in world",result.get_impl()->get_world().id());
//    double n2=result.norm2();
//    double n1=source.norm2();
//    world.gop.fence();
//    print("norm(source), norm(result)", n1, n2);
//
//    MADNESS_EXCEPTION("\n\nending now\n\n",1);


            cpu0 = cpu1;
            K.set_algorithm(Exchange<double, 3>::large_memory);
            K.set_symmetric(true);
            K.set_printlevel(10);
            const vecfuncT reference = K(calc.amo);
            cpu1 = cpu_time();
            double norm = norm2(world, reference);
            if (world.rank() == 0)
                printf("timings exchange operator no multiworld largemem   %8.2fs, norm %.15e\n", cpu1 - cpu0, norm);

            cpu0 = cpu1;
            K.set_algorithm(Exchange<double, 3>::multiworld_efficient);
            tmp = K(calc.amo);
            cpu1 = cpu_time();
            err = norm2(world, reference - tmp);
            if (world.rank() == 0)
                printf("timings exchange operator efficient                %8.2fs, error %.2e\n", cpu1 - cpu0, err);

            cpu0 = cpu1;
            K.set_algorithm(Exchange<double, 3>::small_memory);
            tmp = K(calc.amo);
            cpu1 = cpu_time();
            err = norm2(world, reference - tmp);
            if (world.rank() == 0)
                printf("timings exchange operator no multiworld smallmem   %8.2fs, error %.2e\n", cpu1 - cpu0, err);
        }
        world.gop.fence();
        world.gop.fence();
        world.gop.fence();
    }

    madness::finalize();
    return 0;
}
