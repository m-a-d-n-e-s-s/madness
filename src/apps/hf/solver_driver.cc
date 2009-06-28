#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include "electronicstructureapp.h"
#include "solver.h"

using namespace madness;

int main(int argc, char** argv)
{
    initialize(argc, argv);

    World world(MPI::COMM_WORLD);

    try {
        // Load info for MADNESS numerical routines
        startup(world,argc,argv);
        FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap(world)));


        // Process 0 reads input information and broadcasts
        ElectronicStructureApp app(world, "input");

        // Warm and fuzzy for the user
        if (world.rank() == 0) {
            print("\n\n");
            print(" MADNESS Hartree-Fock and Density Functional Theory Program");
            print(" ----------------------------------------------------------\n");
            print("\n");
            app.entity().print();
            print("\n");
            //app.params().print(world);
        }

        app.make_nuclear_potential();
        app.initial_guess();
        ElectronicStructureParams params = app.params();
        Function<double,3> vnucrhon = app.vnucrhon();
        vecfuncT orbs = app.orbitals();
        std::vector<double> eigs;
        std::vector< Function< std::complex<double>,3> > phis;
        std::vector<double> tmpe = app.eigs();
        print(tmpe);
        int neps = eigs.size();
        for (int i = 0; i < neps; i++)
        {
          phis.push_back(orbs[i]);
          eigs.push_back(tmpe[i]);
        }

        Solver<double,3> dftcalc(world, vnucrhon, app.orbitals(), app.eigs(), 
                                 app.kpoints(), app.occs(), app.params(),
                                 app.entity());
        dftcalc.solve();
        world.gop.fence();

    } catch (const MPI::Exception& e) {
        //        print(e);
        error("caught an MPI exception");
    } catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    } catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    } catch (const char* s) {
        print(s);
        error("caught a string exception");
    } catch (char* s) {
        print(s);
        error("caught a string exception");
    } catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    } catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    } catch (...) {
        error("caught unhandled exception");
    }

    finalize();

    return 0;
}
