//
// Created by Florian Bischoff on 11/1/21.
//



#include<madness/mra/mra.h>
#include<madness/mra/vmra.h>
#include<madness/chem/localizer.h>
#include<madness/chem/molecularbasis.h>
#include<madness/chem/MolecularOrbitals.h>
#include<madness/chem/potentialmanager.h>
#include<madness/chem/correlationfactor.h>
#include<madness/world/test_utilities.h>
#include<madness/chem/SCFOperators.h>
#include<madness/chem/SCF.h>
#include<madness/world/timing_utilities.h>
#include<nemo.h>
#include<write_test_input.h>

#include<vector>

using namespace madness;

/*
 * test different localization methods including core-valence separation
 */

double fock_trace(const Tensor<double>& fock) {
    MADNESS_CHECK(fock.ndim()==2);
    MADNESS_CHECK(fock.dim(0)==fock.dim(1));
    double result=0.0;
    for (int i=0; i<fock.dim(0); ++i) result+=fock(i,i);
    return result;
}

/// distance between v1 and v2
double dist(const Vector<double,3> v1, const Vector<double,3> v2) {
    return (v1-v2).normf();
}

/// angle in degrees between v1--v2--v3 (v2 being the center)
double bend(const Vector<double,3> v1, const Vector<double,3> v2, const Vector<double,3> v3) {
    Vector<double,3> w1=v1-v2;
    Vector<double,3> w2=v3-v2;
    double denom=inner(w1,w2);
    double num=w1.normf()*w2.normf();
    return acos(denom/num)/constants::pi*180.0;
}

bool is_close(const Vector<double,3>& v1, const Vector<double,3>& v2, const double thresh=1.e-3) {
    return ((v1-v2).normf()<thresh);
}
bool is_aligned(const Vector<double,3>& v1, const Vector<double,3>& v2, const Vector<double,3>& v3,
                const double thresh=1.e-3) {
    const Vector<double,3> a1=(v1-v2)*(1.0/(v1-v2).normf());
    const Vector<double,3> a2=(v1-v3)*(1.0/(v1-v3).normf());
    return (fabs(fabs(inner(a1,a2))-1.0)<thresh);
}

/// check if v3 is orthogonal on v1-v2 midbond
bool is_centered_normal(const Vector<double,3>& v1, const Vector<double,3>& v2, const Vector<double,3>& v3,
                   const double thresh=1.e-3) {
    const Vector<double,3> midbond=v2+0.5*(v1-v2);
    const Vector<double,3> normal=v1-v2;
    double inplane=fabs(inner(v3-midbond,normal));
    return (inplane<thresh);
}


/// test localized orbitals: should be pointing towards the edges of a tetrahedron
bool test_ne_boys(World& world, const Nemo& nemo) {
    test_output tout("testing ne boys localization");

    MolecularOrbitals<double,3> mos;
    mos.update_mos_and_eps(nemo.get_calc()->amo,nemo.get_calc()->aeps);
    mos.set_all_orbitals_occupied();
    mos.recompute_localize_sets();
    mos.pretty_print("initial orbitals after a few SCF iterations");

    Localizer localizer(world,nemo.get_calc()->aobasis,nemo.molecule(),nemo.get_calc()->ao);

    localizer.set_method("boys");
    localizer.set_enforce_core_valence_separation(true);
    localizer.print_info();
    Tensor<double> fock1 =nemo.compute_fock_matrix(mos.get_mos(), mos.get_occ());

    auto lmo=localizer.localize(mos, fock1, true);
    Tensor<double> fock2 =nemo.compute_fock_matrix(lmo.get_mos(), lmo.get_occ());
    lmo.pretty_print("cv-localized MOs");
    print(fock2);
    std::vector<Vector<double,3>> center=lmo.compute_center(nemo.get_ncf_ptr()->square());
    for (size_t i=0; i<lmo.get_mos().size(); ++i) {
        print("center of mo",i,center[i][0],center[i][1],center[i][2]);
    }
    print("distances ne tetrahedron",dist(center[1],center[2]),
          dist(center[1],center[3]),
          dist(center[1],center[4]),
          dist(center[2],center[3]),
          dist(center[2],center[4]),
          dist(center[3],center[4]));
    print("angles ne tetrahedron",
          bend(center[1],center[0],center[2]),
          bend(center[1],center[0],center[3]),
          bend(center[1],center[0],center[4]),
          bend(center[2],center[0],center[3]),
          bend(center[2],center[0],center[4]),
          bend(center[3],center[0],center[4]));

    bool success=true;
    // distances depend on the ao basis..
//    double refdistance=localizer.get_aobasis().get_name()=="STO-3G" ? 0.7158 :  0.51877 ;
//    success=success and (fabs(dist(center[1],center[2])-refdistance)<0.001);
//    success=success and (fabs(dist(center[1],center[3])-refdistance)<0.001);
//    success=success and (fabs(dist(center[1],center[4])-refdistance)<0.001);
//    success=success and (fabs(dist(center[2],center[3])-refdistance)<0.001);
//    success=success and (fabs(dist(center[2],center[4])-refdistance)<0.001);
//    success=success and (fabs(dist(center[3],center[4])-refdistance)<0.001);

    success=success and (fabs(bend(center[1],center[0],center[2])-109.4712)<0.001);
    success=success and (fabs(bend(center[1],center[0],center[3])-109.4712)<0.001);
    success=success and (fabs(bend(center[1],center[0],center[4])-109.4712)<0.001);
    success=success and (fabs(bend(center[2],center[0],center[3])-109.4712)<0.001);
    success=success and (fabs(bend(center[2],center[0],center[4])-109.4712)<0.001);
    success=success and (fabs(bend(center[3],center[0],center[4])-109.4712)<0.001);

    tout.end(success);
    return success;

}

/// test localized orbitals

/// tests invariance of the Fock matrix trace, core-valence separation and location of orbitals on bonds
bool test_ethylene(World& world, const Nemo& nemo, const std::string geometry="ethylene") {
    test_output tout("testing localization on: "+geometry);
//    tout.set_cout_to_terminal();
    bool success=true;

    // number of orbitals per bond
    std::map<std::pair<int,int>,int> bonds;
    std::map<std::pair<int,int>,int> bananabonds;
    if (geometry=="ethylene") {
        bonds[{0,2}]=0; // C1--H1
        bonds[{0,3}]=0; // C1--H2
        bonds[{1,4}]=0; // C2--H3
        bonds[{1,5}]=0; // C2--H4
        bonds[{0,1}]=0; // C1--C2
        bananabonds[{0,1}]=0; // C1--C2
    } else if (geometry=="methane") {
        bonds[{0, 1}] = 0; // C1--H1
        bonds[{0, 2}] = 0; // C1--H2
        bonds[{0, 3}] = 0; // C1--H3
        bonds[{0, 4}] = 0; // C1--H4
    }

    Molecule mol=nemo.molecule();
    MolecularOrbitals<double,3> mos;
    mos.update_mos_and_eps(nemo.get_calc()->amo,nemo.get_calc()->aeps);
    mos.set_all_orbitals_occupied();
    mos.recompute_localize_sets();
    mos.pretty_print("initial orbitals after a few SCF iterations");

    Tensor<double> fock =nemo.compute_fock_matrix(mos.get_mos(), mos.get_occ());

//    Tensor<double> fock, overlap;
    double trace_canonical=fock_trace(fock);
    print("initial fock matrix");
    print(fock);

    Localizer localizer(world,nemo.get_calc()->aobasis,nemo.molecule(),nemo.get_calc()->ao);

//    for (std::string method : {"boys","pm","new"}) {
    for (std::string method : {"boys","new"}) {
        for (bool enforce_cv : {true, false}) {

            for (auto& bond : bonds) bond.second=0;
            for (auto& bond : bananabonds) bond.second=0;
            localizer.set_method(method);
            localizer.set_enforce_core_valence_separation(enforce_cv);
            localizer.set_metric(nemo.get_ncf_ptr()->function());
            localizer.print_info();

            auto lmo=localizer.localize(mos, fock, true);
//            if (enforce_cv) lmo.print_cubefiles("mo_ethylene"+method,mol.cubefile_header());
            nemo.get_calc()->amo=lmo.get_mos();

            Tensor<double> fock2 = nemo.compute_fock_matrix(lmo.get_mos(), lmo.get_occ());
            print("final fock matrix for",method);
            print(fock2);
            if (enforce_cv) {
                bool success1=Localizer::check_core_valence_separation(fock2, lmo.get_localize_sets());
                tout.checkpoint(success1,"core-valence separation for "+method);
                success=success and success1;
            }
            double trace_local=fock_trace(fock2);
            print("canonical trace",trace_canonical);
            print("local trace    ",trace_local);
            print("difference, rel. error, thresh ",fabs(trace_local-trace_canonical),
                  fabs(trace_local-trace_canonical)/trace_canonical,FunctionDefaults<3>::get_thresh());
            bool success1=(fabs(trace_local-trace_canonical)/trace_local<FunctionDefaults<3>::get_thresh());
            tout.checkpoint(success1,"trace of the fock matrix for "+method);
            success=success and success1;


            std::vector<Vector<double,3>> center=lmo.compute_center(nemo.get_ncf_ptr()->square());
            for (auto c : center) {
                for (auto& bond : bonds) {
                    int i=bond.first.first;
                    int j=bond.first.second;
                    Vector<double,3> iatom=mol.get_atom(i).get_coords();
                    Vector<double,3> jatom=mol.get_atom(j).get_coords();
                    if (is_close(c,iatom) or is_close(c,jatom)) continue; // ignore core orbitals
                    if (is_aligned(iatom,jatom,c)) bond.second++;
                }
                for (auto& bond : bananabonds) {
                    int i=bond.first.first;
                    int j=bond.first.second;
                    Vector<double,3> iatom=mol.get_atom(i).get_coords();
                    Vector<double,3> jatom=mol.get_atom(j).get_coords();
                    if (is_centered_normal(iatom,jatom,c) and
                            (not is_aligned(iatom,jatom,c))) bond.second++;
                }
            }
            mol.print();
            for (auto bond : bonds) {
                print("found ",bond.second,"bonds between atoms",bond.first.first,bond.first.second);
            }
            for (auto bond : bananabonds) {
                print("found ",bond.second,"banana bonds between atoms",bond.first.first,bond.first.second);
            }
            bool success2=true;
            if (geometry=="ethylene") {
                success2=success2 and bonds[{0,2}]==1; // C1--H1
                success2=success2 and bonds[{0,3}]==1; // C1--H2
                success2=success2 and bonds[{1,4}]==1; // C2--H3
                success2=success2 and bonds[{1,5}]==1; // C2--H4
                if (method=="boys") success2=success2 and bananabonds[{0,1}]==2; // C1-C2
                if (method=="new" or method=="pm") success2=success2 and bonds[{0,1}]==2;
            } else if (geometry=="methane") {
                success2=success2 and bonds[{0,1}]==1; // C--H1
                success2=success2 and bonds[{0,2}]==1; // C--H2
                success2=success2 and bonds[{0,3}]==1; // C--H3
                success2=success2 and bonds[{0,4}]==1; // C--H4
            }
            tout.checkpoint(success2,"center of the local orbitals for "+method);
            success=success and success2;

        }
    }
    tout.end(success);
    return success;

}

int main(int argc, char **argv) {
    int result = 0;
    {
        World& world = madness::initialize(argc, argv);
        startup(world, argc, argv);
        print("entering test_localizer");
        timer t(world);
        commandlineparser parser(argc,argv);
        parser.print_map();
        const int k = 8;
        const double thresh = 1.e-5;
//        const double L = 20.0;
        const double L = 50.0;
        FunctionDefaults<3>::set_cubic_cell(-L, L);
        FunctionDefaults<3>::set_thresh(thresh);
        FunctionDefaults<3>::set_k(k);

        std::string geometry="ethylene";
        Nemo::NemoCalculationParameters nemo_param;
        CalculationParameters param;
        param.set_user_defined_value<std::vector<double>>("protocol",{1.e-5});
        param.set_user_defined_value("k",8);
        param.set_user_defined_value("econv",1.e-4);
        param.set_user_defined_value("maxiter",0);
        param.set_user_defined_value("xc",std::string("lda"));
        param.set_user_defined_value("localize",std::string("canon"));
        param.set_user_defined_value("print_level",2);
//        param.set_user_defined_value("ncf",std::pair<std::string,double>("none",0.0));
        write_test_input test_input(param);
        parser.set_keyval("input",test_input.filename());
        parser.set_keyval("geometry","source_type=library");
        parser.set_keyval("geometry","source_name="+geometry);
        Molecule molecule(world, parser);
        Nemo nemo(world,param,nemo_param,molecule);
        param.print("dft","end");

        test_output tout1(geometry+" prep calculation");
        nemo.value();
        tout1.end(true);

        bool success=test_ethylene(world,nemo,geometry);
        if (not success) result++;

//        parser.set_keyval("structure","ne");
        parser.set_keyval("geometry","source_type = library");
        parser.set_keyval("geometry","source_name = ne");
        Nemo nemo1(world,parser);
        test_output tout("ne prep calculation");
        nemo1.value();
        tout.end(true);
        success=test_ne_boys(world,nemo1);
        if (not success) result++;
        print("result", result);
        t.end("finished localizer test");
    }
    madness::finalize();
    return result;

}
