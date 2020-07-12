/*
 * pno_integrals.hpp
 *
 *  Created on: Feb. 28, 2020
 *      Author: jsk


 CURRENTLY TO PLAY AROUND
 A LOT IS HARDCODED FOR H2
 SO BE CAREFUL WITH DIFFERENT INPUTS


 */


//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <iomanip>
#include <madness/mra/vmra.h>
#include <chem/SCF.h>
#include <chem/nemo.h>
#include <chem/PNO.h>

using namespace madness;

// DEFINE PARAMETER TAGS FOR THE INPUT FILE
const std::string TAG_PNO = "pno";
const std::string TAG_F12 = "f12";
const std::string TAG_CP = "computeprotocol";

// this needs to be added to include
#include "NumCpp.hpp"

//std::vector<real_function_3d> gram_schmidt(World world, std::vector<real_function_3d> &in_basis){
  //auto out_basis = in_basis;
  //out_basis[0] = in_basis[0];
  //std::cout << "norm of 0: " << std::sqrt(out_basis[0].inner(out_basis[0])) << std::endl;
  //for (int i=1; i<in_basis.size(); ++i) {
      //out_basis[i] = in_basis[i];
      //for (int j=0; j<i; ++j) {
          //Q = madness::QProjector<double, 3> (world, out_basis[j]);
          //out_basis[i] = Q(out_basis[i]);
      //}
  //std::cout << "norm of " << i << ": " << std::sqrt(out_basis[i].inner(out_basis[i])) << std::endl;
  //}

  //return out_basis;
//}


int main(int argc, char** argv) {
    madness::initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    if (world.rank() == 0) printf("starting at time %.1f\n", wall_time());
    const double time_start = wall_time();
    std::cout.precision(6);

    startup(world,argc,argv,true);
    print_meminfo(world.rank(), "startup");

    // Get the name of the input file (if given)
    const std::string input = (argc > 1) ? argv[1] : "input";

    // get the orthogonalization method (cholesky, canonicalize, gs -> Gram-Schmidt, none) and basis_size
    const std::string orthogonalization = (argc > 2) ? argv[2] : "cholesky";
    const int basis_size = (argc > 3) ? std::atoi(argv[3]) : 10;

    const bool only_diag = (argc > 4) ? bool(std::atoi(argv[4])) : false;

    if(world.rank()==0){
        std::cout << "\n\n";
        std::cout << "-------------------------------------------------------------------------------------\n";
        std::cout << "SOLVING MRA-PNO-F12 as described in \n";
        std::cout << "J.S. Kottmann, F.A. Bischoff, E.F. Valeev\n";
        std::cout << "Direct determination of optimal pair-natural orbitals in a real-space representation:\n";
        std::cout << "the second-order Møller-Plesset energy\n";
        std::cout << "Journal of Chemical Physics ... 2020\n";
        std::cout << "-------------------------------------------------------------------------------------\n";
        std::cout << "\n\n";

        std::cout << "This script will run PNO-MP2 and print out tensors in binary\n";
        std::cout << "Call as: pno_integrals inputfile orthogonalization basis_size";
        std::cout << "input is " << input << "\n";
        std::cout << "orthogonalization is " << orthogonalization << "\n";
        std::cout << "basis size is " << basis_size << "\n";
        std::cout << "only diag is " << only_diag << "\n";
    }




    // Compute the SCF Reference
    const double time_scf_start = wall_time();
    std::shared_ptr<SCF> calc(new SCF(world, input));
    Nemo nemo(world, calc, input);
    nemo.get_calc()->param.print();
    const double scf_energy = nemo.value();
    if (world.rank() == 0) print("nemo energy: ", scf_energy);
    if (world.rank() == 0) printf(" at time %.1f\n", wall_time());
    const double time_scf_end = wall_time();

    // Compute MRA-PNO-MP2-F12
    const double time_pno_start = wall_time();
    PNOParameters parameters(world,input,nemo.get_calc()->molecule,TAG_PNO);
    F12Parameters paramf12(world, input, parameters, TAG_F12);
    PNO pno(world, nemo, parameters, paramf12);
    std::vector<PNOPairs> all_pairs;
    pno.solve(all_pairs);
    const double time_pno_end = wall_time();


    if(world.rank()==0){
        std::cout << std::setfill(' ');
        std::cout << "\n\n\n";
        std::cout << "--------------------------------------------------\n";
        std::cout << "MRA-PNO-MP2-F12 ended \n";
        std::cout << "--------------------------------------------------\n";
        std::cout << std::setw(25) << "time scf" << " = " << time_scf_end - time_scf_start << "\n";
        std::cout << std::setw(25) << "energy scf" << " = " << scf_energy << "\n";
        std::cout << "--------------------------------------------------\n";
    }
    double mp2_energy = 0.0;
    std::cout<< std::setw(25) << "time pno" << " = " << time_pno_end - time_pno_start << "\n";
    for(const auto& pairs: all_pairs){
        if(pairs.type == MP2_PAIRTYPE){
            mp2_energy = pairs.energies.total_energy();
        }
        std::pair<size_t, size_t> ranks= pno.get_average_rank(pairs.pno_ij);
        if(world.rank()==0){
            std::string name;
            std::stringstream ss;
            ss << pairs.type;
            ss >> name;
            std::cout<< std::setw(25) << "energy "+name << " = " << pairs.energies.total_energy() << "\n";
            std::cout<< std::setw(25) << "average pno rank " + name << " = " << ranks.first << "\n";
            std::cout<< std::setw(25) << "max pno rank " + name << " = " << ranks.second << "\n";
        }
    }
    if(world.rank()==0 and mp2_energy != 0.0){
        std::cout << "--------------------------------------------------\n";
            std::cout<< std::setw(25) << "energy(total)" << " = " << scf_energy + mp2_energy << "\n";
            std::cout << "--------------------------------------------------\n";
            std::cout << "\n\n\n";
    }

    // compute orthogonalized mp2 basis and print out hamiltonian tensors
    std::cout << std::setprecision(8);
    std::cout << std::fixed;
    std::cout << std::showpos;

    const bool canonicalize = orthogonalization == "canonicalize";
    const bool orthogonalize = orthogonalization != "none";
    const double h_thresh = 1.e-4;
    const double thresh = 1.e-4;//parameters.thresh();
    const auto amo = nemo.get_calc()->amo;

    if(world.rank()==0) std::cout << "Tightening thresholds to 1.e-6 for post-processing\n";
    FunctionDefaults<3>::set_thresh(1.e-6);

    for(auto& pairs: all_pairs){
        const auto& pno_ij = pairs.pno_ij;
        const auto& rdm_evals = pairs.rdm_evals_ij;
        const bool is_gs = pairs.type == MP2_PAIRTYPE;
        std::string name = "gs";

        vecfuncT reference = amo;
        if (not is_gs){
            const auto& x = pairs.cis.x;
            reference.insert(reference.end(), x.begin(), x.end());
            name = "ex" + std::to_string(pairs.cis.number);
        }

        std::vector<real_function_3d> all_basis_functions;// = nemo.get_calc()->amo;

        std::vector<double> occ;
        // collect PNOs from all pairs and sort by occupation number
        for(ElectronPairIterator it=pno.pit();it;++it){
            if (only_diag and not it.diagonal()){
                std::cout << "skipping pair " << it.name() << "\n";
                continue;
            }else{
                const auto& pair = pno_ij[it.ij()];
                all_basis_functions.insert(all_basis_functions.end(), pair.begin(), pair.end());
                for (auto ii=0; ii<rdm_evals[it.ij()].size();++ii){
                    occ.push_back(rdm_evals[it.ij()][ii]);
                }
            }
        }
        for (auto i=0; i<rdm_evals.size(); ++i){
            for (auto ii=0; ii<rdm_evals[i].size();++ii){
                occ.push_back(rdm_evals[i][ii]);
            }
        }
        std::vector<std::pair<double, real_function_3d> > zipped;
        for (auto i=0; i< all_basis_functions.size(); ++i){
            zipped.push_back(std::make_pair(occ[i], all_basis_functions[i]));
        }

        std::sort(zipped.begin(), zipped.end(), [](const auto& i, const auto& j) { return i.first > j.first; });

        std::vector<double> unzipped_first;
        std::vector<real_function_3d> unzipped_second;
        for (auto i=0; i<basis_size;++i){
            unzipped_first.push_back(zipped[i].first);
            unzipped_second.push_back(zipped[i].second);
        }
        occ = unzipped_first;
        all_basis_functions = unzipped_second;

        if(world.rank()==0){
            std::cout << "all used occupation numbers:\n" << occ << "\n";
        }

        // reference projector (not fullfilled for CIS)
        madness::QProjector<double, 3> Q(world, reference);
        all_basis_functions = Q(all_basis_functions);

        // compute overlap for cholesky decomposition
        const auto S = madness::matrix_inner(world, all_basis_functions, all_basis_functions, true);
        if(world.rank()==0) std::cout << "Overlap Matrix of all PNOs:\n";
        for (int i=0;i<all_basis_functions.size();++i){
            for (int j=0;j<all_basis_functions.size();++j){
                std::cout << S(i,j) << " ";
            }
            if(world.rank()==0) std::cout << "\n";
        }
        auto gop = std::shared_ptr < real_convolution_3d > (madness::CoulombOperatorPtr(world, 1.e-6, parameters.op_thresh()));

        
        auto basis = all_basis_functions;

        if (orthogonalize){

            //basis = madness::orthonormalize_rrcd(all_basis_functions, 1.e-5);
            //Use standard cd, since pivoting swaps PNOs around
            if (orthogonalization == "cholesky") {    
                basis = madness::orthonormalize_cd(all_basis_functions);
                if(world.rank()==0) std::cout << "Basis size after global Cholesky: " << basis.size() << "\n";
            }

            //do  gram-schmidt
            else if (orthogonalization == "gs") {
                std::cout << "orthonormalize...\n";
                // u_0 = v_0
                double norm_i = std::sqrt(basis[0].inner(basis[0]));
                basis[0].scale(1.0/norm_i);
                // Gram-Schmidt iterations
                for (int i=1; i<basis.size(); ++i) {
                    //basis[i] = in_basis[i];
                    for (int j=0; j<i; ++j) {
                        std::vector<real_function_3d> basis_j(1);
                        basis_j[0] = basis[j];
                        Q = madness::QProjector<double, 3> (world, basis_j);
                        basis[i] = Q(basis[i]);
                    }
                norm_i = std::sqrt(basis[i].inner(basis[i]));
                basis[i].scale(1.0/norm_i);
                }
            }

            if(world.rank()==0) std::cout << "Basis size after Gram-Schmidt: " << basis.size() << "\n";

        }

        if(world.rank()==0) std::cout << "Adding Reference orbitals\n";
        const auto amo = nemo.get_calc()->amo;
        basis.insert(basis.begin(), reference.begin(), reference.end());



       // |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

        auto fop = std::shared_ptr < real_convolution_3d > (SlaterF12OperatorPtr(world, paramf12.gamma(), 1.e-6, parameters.op_thresh()));
	// auto slaterop = SlaterOperator(world, paramf12.gamma(), 1.e-6, parameters.op_thresh());
        fop->is_slaterf12 = false; // make sure it is the correct formulation of the exponent
    

       /* Here (or 3 lines further up), read in {cabs} basis if required
         * -> read in {cabs}
         * -> orthogonalize {cabs}
         * -> {all_basis_functions} = {all_PNOs} + {cabs}
         * -> compute integrals based on {all_basis_functions}
         */
        // > define Q here (reference proj now!)
        // > different routine for orthonormalization? -> GS?
        // > check whether we only project out PNO's used in PNO-basis!
        bool corr_cabs = true; // todo: read from params file
        auto pno_plus_ref = basis;
        if (corr_cabs) {
            std::vector<real_function_3d> cabs;
            cabs = pno.f12.read_cabs_from_file(paramf12.auxbas_file()); // sadly, F12Potential f12 is private member of pno...
            if (not cabs.empty()){
                MyTimer time2 = MyTimer(world).start();
                // Project out reference 
                //cabs = Q(cabs); 
                // Project out {pno} + ref
                madness::QProjector<double, 3> Qpno(world, pno_plus_ref);
                cabs = Qpno(cabs);
                // Orthonormalize {cabs}
                cabs = orthonormalize_cd(cabs);
                //for (ElectronPairIterator it = pno.pit(); it; ++it) {
                     //right now this will make the same guess for all pairs
                    //const vector_real_function_3d tmp=guess_virtuals(param.abs);
                    //madness::QProjector<double, 3> Qpno(world, pno[it.ij()]);
                    //cabs = Qpno(cabs);
                    //const vector_real_function_3d tmp = Qpno(cabs);
                    //abs_ij[it.ij()] = tmp;
                //}
               time2.stop().print("Make pair specific ABS from PNOS and " + std::to_string(cabs.size()) + " functions");
            }
            else if (cabs.empty()) {
                std::cout << "Complaining..." << std::endl; // todo: raise exception or so
            }
        
            // Merge {cabs} + {pno}
            // necessary?
            if(world.rank()==0) std::cout << "Adding {cabs} to {pno+ref}.\n";
            if(world.rank()==0) std::cout << "Size before: " << basis.size() << ".\n";
            basis.insert(basis.end(), cabs.begin(), cabs.end());
            if(world.rank()==0) std::cout << "Size after: " << basis.size() << ".\n";

        
        }
          
          
        
        // ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| 




        madness::Tensor<double> g(basis.size(), basis.size(), basis.size(), basis.size()); // using mulliken notation since thats more efficient to compute here: Tensor is (pq|g|rs) = <pr|g|qs>
        madness::Tensor<double> f(basis.size(), basis.size(), basis.size(), basis.size());
        // madness::Tensor<double> slater(basis.size(), basis.size(), basis.size(), basis.size());


        if(canonicalize){
            if(world.rank()==0) std::cout << "canonicalizing!\n";
            auto F = madness::Fock(world, &nemo);
            const auto Fmat = F(basis, basis);
            Tensor<double> U, evals;
            syev(Fmat, U, evals);
            basis = madness::transform(world, basis, U);
        }

        std::vector<vecfuncT> PQ;
        for (const auto& x : basis){
            PQ.push_back(madness::truncate(x*basis,thresh));
        }
        std::vector<vecfuncT> GPQ;
        std::vector<vecfuncT> FPQ;
        //std::vector<vecfuncT> slaterPQ;
        for (const auto& x : basis){
            GPQ.push_back(madness::truncate(madness::apply(world, *gop, madness::truncate(x*basis,thresh)), thresh));
            FPQ.push_back(madness::truncate(madness::apply(world, *fop, madness::truncate(x*basis,thresh)), thresh));
            //slaterPQ.push_back(madness::truncate(madness::apply(world, slaterop, madness::truncate(x*basis,thresh)), thresh));
        }

        auto J = madness::Coulomb(world, &nemo);
        auto K = madness::Exchange<double, 3>(world, &nemo, 0);
        auto Jmat = J(basis, basis);
        auto Kmat = K(basis, basis);

        // delete me ---------------------------------------------
        for (auto p=0; p<basis.size(); p++){
            for (auto q=0; q<basis.size(); q++){
		    std::cout << Jmat(p,q) << "\t";
	    }
	    std::cout << std::endl;
	 }
        for (auto p=0; p<basis.size(); p++){
            for (auto q=0; q<basis.size(); q++){
		    std::cout << Kmat(p,q) << "\t";
	    }
	    std::cout << std::endl;
	 }
        // delete me ---------------------------------------------

        int non_zero=0, non_zero_f=0;
        if(world.rank() ==0 ) std::cout << "Compute G Tensor:\n";
        for (auto p=0; p<basis.size(); p++){
            for (auto q=0; q<basis.size(); q++){
                if (PQ[p][q].norm2() < h_thresh) continue;
                for (auto r=0; r<basis.size(); r++){
                    for (auto s=0; s<basis.size(); s++){
                        if (GPQ[r][s].norm2() < h_thresh) continue;
                        else{
                            g(p,q,r,s) = PQ[p][q].inner(GPQ[r][s]);
                            if(canonicalize and p==q){
                                g(p,q,r,s) += Kmat(r,s) - Jmat(r,s);
                            }else if(canonicalize and r==s){
                                g(p,q,r,s) += Kmat(p,q) - Jmat(p,q);
                            }
                            if(std::fabs(g(p,q,r,s)) > h_thresh ){
                                if(world.rank()==0 and basis.size() < 3) std::cout << " g " << p  << " "<<  q  << " "<<  r  << " "<<  s << " = " << g(p,q,r,s) << "\n";
                                ++non_zero;
                            }
                        }
                        if (FPQ[r][s].norm2() < h_thresh) continue;
                        else{
                            f(p,q,r,s) = PQ[p][q].inner(FPQ[r][s]);
                            //slater(p,q,r,s) = PQ[p][q].inner(slaterPQ[r][s]);
                            if(std::fabs(f(p,q,r,s)) > h_thresh ){
                                if(world.rank()==0 and basis.size() < 3) std::cout << " f " << p  << " "<<  q  << " "<<  r  << " "<<  s << " = " << g(p,q,r,s) << "\n";
                                ++non_zero_f;
                            }
                        }
                    }
                }
            }
        }
        int non_zero_h = 0;

        Tensor<double> h;
        if(canonicalize){
            auto F = madness::Fock(world, &nemo);
            h = F(basis, basis);
        }else{
            auto T = madness::Kinetic<double, 3>(world);
            auto V = madness::Nuclear(world, &nemo);
            h = T(basis,basis) + V(basis,basis);
        }
        for (int i=0;i<basis.size();++i){
            for (int j=0;j<basis.size();++j){
                if(std::fabs(h(i,j)) > h_thresh){
                    if(world.rank()==0 and basis.size() < 3) std::cout << " h " << i << " "<< j << " "<< h(i,j) << "\n";
                    ++non_zero_h;
                }
            }
        }

        if(world.rank()==0) {
            std::cout << "non zero elements:\n g : " << non_zero;
            std::cout << "\n f :" << non_zero_f <<"\n h : " << non_zero_h << "\n";
        }

        h = h.flat();
        nc::NdArray<double> hh(h.ptr(), h.size(), 1);
        hh.tofile(name+"_hcore.bin", "");

        g = g.flat();
        nc::NdArray<double> gg(g.ptr(), g.size(), 1);
        gg.tofile(name+"_gtensor.bin", "");

        f = f.flat();
        nc::NdArray<double> ff(f.ptr(), f.size(), 1);
        ff.tofile(name+"_f12tensor.bin", "");
        //slater = slater.flat();
        //nc::NdArray<double> ss(slater.ptr(), slater.size(), 1);
        //ss.tofile(name+"_f12_tensor.bin", "");

        if (not orthogonalize){
            auto S = madness::matrix_inner(world, basis, basis, true);
            S = S.flat();
            nc::NdArray<double> gg(S.ptr(), S.size(), 1);
            gg.tofile(name+"_overlap.bin", "");
        }

        auto Fop =  madness::Fock(world, &nemo);
        auto F = Fop(basis, basis);
        std::cout << "F\n" << F << "\n";

        const auto Stest = madness::matrix_inner(world, basis, basis, true);
        std::cout << "Overlap over whole basis\n" << Stest << "\n";

    }




    world.gop.fence();
    if (world.rank() == 0) printf("finished at time %.1f\n", wall_time());


    print_stats(world);
    finalize();

    return 0;

}
