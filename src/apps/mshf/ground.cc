#include "input.h"
#include "densities.h"
#include "loadbalance.h"
#include "states.h"
#include "potential.h"
#include "energies.h"
#include "output.h"
#include "orthonorm.h"
#include "iterate.h"
#include "ground.h"

double ttt, sss;
double details = 1.0;

#define START_TIMER world.gop.fence(); ttt=wall_time(); sss=cpu_time()
#define END_TIMER(msg) ttt=wall_time()-ttt; sss=cpu_time()-sss; if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt)


// Directory for checkpointing files
//------------------
char record[] = "./";
char bck_record[] = "./";


// Calculates nuclear ground state via Hartree-Fock iterations
void ground_state(World& world, const double A,
                                const double Z,
                                const double length,
                                const int initial,
                                const int boundary,
                                const int jellium,
                                const int spinorbit,
                                const int meff,
                                const int screening,
                                const double screenl,
                                const int avg_pot,
                                const int avg_lap,
                                const int avg_wav,
                                const int lap_comp,
                                const double prec,
                                const int vtk_output,
                                const int txt_output,
                                const int use_infile,
                                double knumber,
                                double thresh,
                                double chi,
                                double tol,
                                double brad,
                                int IO_nodes,
                                const real_tensorT b,
                                const real_tensorT bp,
                                const double alpha,
                                const double k_fn,
                                const double k_fp,
                                const int project,
                                const int timing)
{
    double time_old     = 0.0;          // Checkpointed wall_time
    const double L      = length/2.0;      // Half of box size

    double chi_in       = chi;
    double brad_in      = brad;
    double tol_in       = tol;

    if (world.rank() == 0) {
        print("Mass number:                     ", A);
        print("Charge number:                   ", Z);
        print("Total box size:                  ", length);
        print("Additional precision factor:     ", prec);
        print("Density laplacian via method:    ", lap_comp);
        print(" ");

        if (boundary == 1) print("- Periodic boundary conditions");
        else if (boundary == 0) print("- Free boundary conditions");
        else {
            print ("Error: Boundary condition must be either 1 (periodic) or 0 (free).");
            assert(0);
        }

        if (jellium == 1) print("- Jellium approximation");
        else if (jellium == 0) print("- No jellium approximation");
        else {
            print ("Error: Jellium switch must be 1 (with jellium) or 0 (no jellium).");
            assert(0);
        }

        if (spinorbit == 1) print("- Spin-orbit potential included");
        else if (spinorbit == 0) print("- No spin-orbit potential");
        else {
            print ("Error: Spin-orbit swithc must be 1 (with spinorbit) or 0 (no spinorbit).");
            assert(0);
        }
        if (meff == 1) print("- Effective mass potential included");
        else if (meff == 0) print("- No effective mass potential");
        else {
            print ("Error: Effective mass switch must be 1 (with eff. mass) or 0 (no eff. mass).");
            assert(0);
        }

        if (screening == 1) {
            print("- Coulomb screening");
            print("  Screening length:                ", screenl);
        }
        else if (screening == 0) print("- No Coulomb screening");
        else {
            print ("Error: Coulomb screening switch must be 1 (with screening) or 0 (no screening).");
            assert(0);
        }

        if (avg_pot == 1) print("- Potential mixing");
        else if (avg_pot == 0) {}
        else {
            print ("Error: Potential mixing switch must be 1 (with mixing) or 0 (no mixing).");
            assert(0);
        }

        if (avg_lap == 1) print("- Laplacian mixing");
        else if (avg_lap == 0) {}
        else {
            print ("Error: Laplacian mixing switch must be 1 (with mixing) or 0 (no mixing).");
            assert(0);
        }

        if (avg_wav == 1) print("- State mixing");
        else if (avg_wav == 0) {}
        else {
            print ("Error: State mixing switch must be 1 (mixing) or 0 (no mixing).");
            assert(0);
        }
    }

    if (world.rank() == 0) {std::cout << " " << std::endl;}
    if (world.rank() == 0) {std::cout << "I/O with" << " " << IO_nodes << " " << "nodes" << std::endl;}
    if (world.rank() == 0) {std::cout << " " << std::endl;}

    // Check, if there is a checkpointed file with updated values for above parameters
    FILE *tfile;
    char trecord_name1[100];
    sprintf(trecord_name1, "%s/checkpoint_t.00000",record);
    tfile=fopen(trecord_name1,"r");
    if (tfile!=NULL) {
        if (world.rank() == 0) {print("Checkpointed files are present - read simulation information");}
        char trecord_name2[100];
        sprintf(trecord_name2, "%s/checkpoint_t",record);
        archive::ParallelInputArchive tin(world, trecord_name2, IO_nodes);
        tin & thresh;
        tin & knumber;
        tin & brad_in;
        tin & chi_in;
        tin & tol_in;
        tin & time_old;
        tin.close();
        fclose(tfile);
    }

    // Use values from checkpointed file (ignore_chk == 0) or input file (ignore_chk == 1)
    if (use_infile == 0){
        brad    = brad_in;
        chi     = chi_in;
        tol     = tol_in;
    }
    else if (use_infile == 1) { }
    else {if (world.rank() == 0) {print("Error; check value for 'use_infile'. Must be either 0 or 1");}}

    // FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap< Key<3> >(world)));
    std::cout.precision(6);
    FunctionDefaults<3>::set_k(knumber);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_initial_level(5);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_cubic_cell(-L, L);

    if (boundary == 1) {
        FunctionDefaults<3>::set_bc(BC_PERIODIC);
        FunctionDefaults<3>::set_truncate_mode(1);
    }

    FunctionDefaults<3>::set_apply_randomize(false);
    FunctionDefaults<3>::set_project_randomize(false);
    FunctionDefaults<3>::set_truncate_on_project(true);

    std::ofstream outfile;

    comp_vecfuncT psi_n(A-Z), psi_p(Z);                                         // Neutron and proton states
    comp_vecfuncT psi_pu = zero_functions<double_complex,3>(world,Z);           // Proton spin-up states
    comp_vecfuncT psi_pd = zero_functions<double_complex,3>(world,Z);           // Proton spin-down states
    comp_vecfuncT psi_nu = zero_functions<double_complex,3>(world,A-Z);         // Neutron spin-up states
    comp_vecfuncT psi_nd = zero_functions<double_complex,3>(world,A-Z);         // Neutron spin-down states

    real_tensorT energy_n(A-Z);                                                 // Energy of states
    real_tensorT energy_p(Z);                                                   // Energy of states

    real_tensorT nx(A-Z), ny(A-Z), nz(A-Z);
    real_tensorT px(Z),   py(Z),   pz(Z);

    real_functionT U_n, U_p;                                            // Neutron and proton potentials
    real_functionT U_n_old, U_p_old;                                    // Old neutron and proton potentials
    real_functionT rho_p, rho_n, rho;                                   // Number densities of neutrons and protons
    real_functionT lap_n_old, lap_p_old;                                // Old lapalcians of number densities 

    double delta_psi     = -10.0;                                       // Maximum change is single particle states
    double delta_psiold  =  20.0;                                       // Maximum change in previous timestep 
    double BE       =   0.0;                                            // Initialization of binding energy
    bool file_check = false;                                            // Marker for checkpointed files


    // Initialization from functions or text file (only done when there is no checkpointed file)
    //--------------------

    // Check if a checkpoint file exists. Only use one file to check for existance
    file_check=false;
    if (world.rank() == 0) {
        char grecord_name1[100];
        sprintf(grecord_name1, "%s/checkpoint_g.00000",record);
        FILE *file;
        file=fopen(grecord_name1,"r");
        if (file==NULL) {file_check = false;}
        else {file_check = true; fclose(file);}
    }
    world.gop.broadcast(file_check);

    // If there is no checkpoint file, initialize wavefunctions. Otherwise read from checkpoint file
    if (file_check == false) {
        if (world.rank() == 0) { print(""); print(""); print("Initialize single particle states");}

        if (initial == 1) {
            if (world.rank() == 0) {print("Gaussians");}
            make_MD(world, psi_n, psi_p, A, Z, L, prec);
            if (world.rank() == 0) {print("Done");}
        }
        else if (initial == 2) {
            if (world.rank() == 0) {print("Harmonic Oscillator");}
            make_HO(world, psi_p, A);
            make_HO(world, psi_n, A);
            if (world.rank() == 0) {print("Done");}
        }
        else if (initial == 3) {
            if (world.rank() == 0) {print("Plane Waves");}
            make_Fermi(world, psi_p, A, Z, L, prec);
            loadbalance_v1(world, psi_p);
            world.gop.fence();
            make_Fermi(world, psi_n, A, Z, L, prec);
            loadbalance_v1(world, psi_n);
            world.gop.fence();
            if (world.rank() == 0) {print("Done");}
        }
        else {
            if (world.rank() == 0) {print("Error: Initalization parameter 'initial' must be 1-3");}
            assert(0);
        }

        truncate2(world, psi_n, psi_p, prec);
        world.gop.fence();
        normalize_2v(world, psi_n, psi_p);
        world.gop.fence();

        spin_split(world, psi_p, psi_pu, psi_pd); // Split into spin-up and down
        spin_split(world, psi_n, psi_nu, psi_nd); // Split into spin-up and down

        world.gop.fence();
        psi_p.clear(); psi_n.clear();
   }

  // ------------------------------
  // End of Initialization



    if (world.rank() == 0) {
        print("");print("");
        print("Main iteration starts");
        print("--------------");
    }


    int iter      = -1;
    int pindex    =  1;
    int terminate =  0;
    while (terminate == 0) {
        if ((sqrt(delta_psi * delta_psi) < 1.e-8) && (knumber >= 9)) {
            if (world.rank() == 0) {print(delta_psi, knumber);}
            terminate = 1;
            world.gop.fence();
        }
        else {
            // Read-in of checkpointed files
            if (file_check == true) {
                if (timing == 1) {START_TIMER;}
                if (world.rank() == 0) {print(" "); print("Reading checkpoint");}

                char grecord_name2[100];
                sprintf(grecord_name2, "%s/checkpoint_g",record);
                archive::ParallelInputArchive gin(world, grecord_name2, IO_nodes);
                gin & iter;
                gin & pindex;
                gin & delta_psi;
                gin & delta_psiold;
                gin & energy_n;
                gin & energy_p;
                gin & U_n;
                gin & U_p;
                gin & lap_n_old;
                gin & lap_p_old;
                gin.close();
                world.gop.fence();

                if(world.rank() == 0) {print("Reading neutrons");}
                char nrecord_name[100];
                sprintf(nrecord_name, "%s/checkpoint_n", record);
                archive::ParallelInputArchive nind(world, nrecord_name, IO_nodes);
                for (unsigned int i=0; i< psi_nu.size(); i++) {nind & psi_nu[i]; nind & psi_nd[i];}
                nind.close();
                world.gop.fence();

                if (world.rank() == 0) {print("Reading protons");}
                char precord_name[100];
                sprintf(precord_name, "%s/checkpoint_p", record);
                archive::ParallelInputArchive pind(world, precord_name, IO_nodes);
                for (unsigned int i=0; i< psi_pu.size(); i++) {pind & psi_pu[i]; pind & psi_pd[i];}
                pind.close();
                world.gop.fence();
                file_check = false;
                if (world.rank() == 0) {print("Reading done"); print(" ");}
                if (timing == 1) {END_TIMER("Reading Input");}
            }

            iter++;
            delta_psi = -1.0;
            double time = wall_time() + time_old;
            if(world.rank() == 0) {
                print("Iteration Nr.:    ", iter);
                print("Wavelet number:   ", knumber);
                print("Trunc. threshold: ", thresh);
                print("Mixing parameter: ", chi);
                if(brad > 0.0) {print("Smoothing radius: ", brad);}
                else {print("No smoothing");}
                print("Wall-time:        ", time);
                print(" ");
            }

            if( BE == 0.0) {
                rho_p = ndensity(world, psi_pu, psi_pd);
                rho_n = ndensity(world, psi_nu, psi_nd);
            }

            truncate2(world, psi_pu, psi_pd, prec);
            truncate2(world, psi_nu, psi_nd, prec);
            world.gop.fence();    

            // Make local potential
            if (world.rank() == 0) {print("Make potential");}
            if (iter > 0) {
                U_p_old = copy(U_p);
                U_n_old = copy(U_n);
            }

            if (timing == 1) {START_TIMER;}
            Potential(world, psi_pu, psi_pd, psi_nu, psi_nd, energy_p, energy_n, U_p, U_n, rho_p,
                      rho_n, rho, iter, BE, delta_psi, tol, brad, U_p_old, U_n_old, lap_p_old, lap_n_old,
                      L, jellium, spinorbit, screening, screenl, avg_pot, avg_lap,
                      lap_comp, prec, vtk_output, txt_output, b, bp, alpha, k_fn, k_fp, timing);
            if (timing == 1) {END_TIMER("Potential");}

            if (world.rank() == 0) {
                print(" ");
                print("Binding energy: ", BE); print(" ");
                print("Calculate new wave functions");
            }

            // Iterate neutrons
            if (world.rank() == 0) {print("Neutrons: ");}
            if (timing == 1) {START_TIMER;}
            iterate(world, U_n, psi_nu, psi_nd, energy_n, delta_psi, rho, rho_n, iter,
                    chi, tol, spinorbit, meff, avg_wav, prec, b, bp, k_fn);
            if (timing == 1) {END_TIMER("Iterate neutrons");}
            world.gop.fence();

            if (timing == 1) {START_TIMER;}
            truncate2(world, psi_nu, psi_nd, prec);
            world.gop.fence();
            loadbalance_v2(world, psi_nu, psi_nd);
            world.gop.fence();
            if (timing == 1) {END_TIMER("loadbalance_v2");}

            // Iterate protons 
            if (world.rank() == 0) {print("Protons: ");}
            if (timing == 1) {START_TIMER;}
            iterate(world, U_p, psi_pu, psi_pd, energy_p, delta_psi, rho, rho_p, iter,
                    chi, tol, spinorbit, meff, avg_wav, prec, b, bp, k_fn);
            if (timing == 1) {END_TIMER("Iterate protons");}
            world.gop.fence();

            if (timing == 1) {START_TIMER;}
            truncate2(world, psi_pu, psi_pd, prec);
            world.gop.fence();
            loadbalance_v2(world, psi_pu, psi_pd);
            world.gop.fence();
            if (timing == 1) {END_TIMER("loadbalance_v2");}

            U_p_old.clear(); U_n_old.clear();
            world.gop.fence();

            // Print out maximum wavefunction change
            if (world.rank() == 0) {
                print(" "); print(" ");
                print("max. single particle state change: ",delta_psi);
            }

            // Output file to keep track of delta_psi
            if (timing == 1) {START_TIMER;}
            if (iter%1 == 0 || iter == 0) {
                if (world.rank() == 0) {
                    outfile.open("data_log.txt", std::ios::app);
                    outfile << iter <<" "<< delta_psi <<" "<< BE <<" "<< thresh <<" "<< time << std::endl;
                    outfile.close();
                }
            }
            if (timing == 1) {END_TIMER("Output data_log.txt");}
            world.gop.fence();  // FENCE gif

            // Project to higher wavelet number
            if ((delta_psi <= 2.0 * A * thresh * prec && knumber < 10) && project == 1) {
                if (timing == 1) {START_TIMER;}
                if (world.rank() == 0) {print("Project");}

                thresh *= 1.e-1;
                knumber += 1;
                if (knumber == 8) {brad = -1.0;}

                FunctionDefaults<3>::set_k(knumber);
                FunctionDefaults<3>::set_thresh(thresh);
                reconstruct(world, psi_pu);
                reconstruct(world, psi_pd);
                reconstruct(world, psi_nu);
                reconstruct(world, psi_nd);
                U_n.reconstruct();
                U_p.reconstruct();
                rho_n.reconstruct();
                rho_p.reconstruct();
                lap_n_old.reconstruct();
                lap_p_old.reconstruct();

                for (unsigned int i = 0; i < psi_pu.size(); i++) {
                    psi_pu[i] = madness::project(psi_pu[i], knumber, thresh);
                    psi_pd[i] = madness::project(psi_pd[i], knumber, thresh);
                }
                truncate2(world, psi_pu, psi_pd, prec);
                world.gop.fence();
                normalize_ud(world, psi_pu, psi_pd);
                world.gop.fence();

                for (unsigned int i = 0; i < psi_nu.size(); i++) {
                    psi_nu[i] = madness::project(psi_nu[i], knumber, thresh);
                    psi_nd[i] = madness::project(psi_nd[i], knumber, thresh);
                }
                truncate2(world, psi_nu, psi_nd, prec);
                world.gop.fence();
                normalize_ud(world, psi_nu, psi_nd);
                world.gop.fence();

                U_n = madness::project(U_n, knumber, thresh);
                U_p = madness::project(U_p, knumber, thresh);

                rho_n = madness::project(rho_n, knumber, thresh);
                rho_p = madness::project(rho_p, knumber, thresh);

                lap_n_old = madness::project(lap_n_old, knumber, thresh);
                lap_p_old = madness::project(lap_p_old, knumber, thresh);

                pindex++;
                world.gop.fence();

                if (timing == 1) {END_TIMER("Reconstruct");}
                rho = rho_p + rho_n;
            }

            if (timing == 1) {START_TIMER;}
            if (iter%1 == 0) {
                if (iter%2 == 0 && iter != 0) {
                    if (world.rank() == 0) { print(" "); print("Checkpoint"); print(" "); }
                    char trecord_name[100];
                    sprintf(trecord_name, "%s/checkpoint_t", bck_record);
                    archive::ParallelOutputArchive tout(world, trecord_name, IO_nodes);
                    tout & thresh;
                    tout & knumber;
                    tout & brad;
                    tout & chi;
                    tout & tol;
                    tout & time;
                    tout.close();
                    world.gop.fence();

                    char grecord_name[100];
                    sprintf(grecord_name, "%s/checkpoint_g", bck_record);
                    archive::ParallelOutputArchive gout(world, grecord_name, IO_nodes);
                    gout & iter;
                    gout & pindex;
                    gout & delta_psi;
                    gout & delta_psiold;
                    gout & energy_n;
                    gout & energy_p;
                    gout & U_n;
                    gout & U_p;
                    gout & lap_n_old;
                    gout & lap_p_old;
                    gout.close();
                    world.gop.fence();

                    char nrecord_name[100];
                    sprintf(nrecord_name, "%s/checkpoint_n", bck_record);
                    archive::ParallelOutputArchive noutc(world, nrecord_name, IO_nodes);
                    for (unsigned int i = 0; i < psi_nu.size(); i++) {noutc & psi_nu[i]; noutc & psi_nd[i];}
                    noutc.close();
                    world.gop.fence();

                    char precord_name[100];
                    sprintf(precord_name, "%s/checkpoint_p", bck_record);
                    archive::ParallelOutputArchive poutc(world, precord_name, IO_nodes);
                    for (unsigned int i = 0; i < psi_pu.size(); i++) {poutc & psi_pu[i]; poutc & psi_pd[i];}
                    poutc.close();
                    world.gop.fence();
                }
                else {
                    if (world.rank() == 0) {print(" "); print("Checkpoint"); print(" ");}
                    char trecord_name[100];
                    sprintf(trecord_name, "%s/checkpoint_t", record);
                    archive::ParallelOutputArchive tout(world, trecord_name, IO_nodes);
                    tout & thresh;
                    tout & knumber;
                    tout & brad;
                    tout & chi;
                    tout & tol;
                    tout & time;
                    tout.close();
                    world.gop.fence();

                    char grecord_name[100];
                    sprintf(grecord_name, "%s/checkpoint_g", record);
                    archive::ParallelOutputArchive gout(world, grecord_name, IO_nodes);
                    gout & iter;
                    gout & pindex;
                    gout & delta_psi;
                    gout & delta_psiold;
                    gout & energy_n;
                    gout & energy_p;
                    gout & U_n;
                    gout & U_p;
                    gout & lap_n_old;
                    gout & lap_p_old;
                    gout.close();
                    world.gop.fence();

                    char nrecord_name[100];
                    sprintf(nrecord_name, "%s/checkpoint_n", record);
                    archive::ParallelOutputArchive noutc(world, nrecord_name, IO_nodes);
                    for (unsigned int i = 0; i < psi_nu.size(); i++) { noutc & psi_nu[i]; noutc & psi_nd[i];}
                    noutc.close();
                    world.gop.fence();

                    char precord_name[100];
                    sprintf(precord_name, "%s/checkpoint_p", record);
                    archive::ParallelOutputArchive poutc(world, precord_name, IO_nodes);
                    for (unsigned int i = 0; i < psi_pu.size(); i++) { poutc & psi_pu[i]; poutc & psi_pd[i]; }
                    poutc.close();
                    world.gop.fence();
                }
            }
            if (timing == 1) {END_TIMER("Checkpointing");}

            delta_psiold = delta_psi;

            if (world.rank() == 0) {
                print("____________________________________________________");
                print("");
                print("");
            }
            world.gop.fence();
        }
    }
    world.gop.fence();
}







