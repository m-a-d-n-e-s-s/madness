//
// Created by Florian Bischoff on 9/28/23.
//

#include<madness.h>
#include<test_utilities.h>


using namespace madness;


int test_conversion(World& world) {
    real_function_2d f=real_factory_2d(world).functor([](const coord_2d& r) {return exp(-r.normf());});
    real_function_2d f1=real_factory_2d(world).functor([](const coord_2d& r) {return exp(-inner(r,r));});
    real_function_2d f2=real_factory_2d(world).functor([](const coord_2d& r) {return inner(r,r)*exp(-2.0*r.normf());});

    real_function_2d f3=copy(f2);
    f3.get_impl()->gaxpy_inplace_reconstructed(1.0,*f1.get_impl(),1.0,false);
    real_function_2d f3ref=real_factory_2d(world).functor([](const coord_2d& r) {return inner(r,r)*exp(-2.0*r.normf()) + exp(-inner(r,r));});
    double f3norm=f3ref.norm2();
    MADNESS_CHECK_THROW(f3.get_impl()->verify_tree_state_local(),"tree state get_tree_state(is invalid");
    MADNESS_CHECK_THROW(get_tree_state(f3)==redundant_after_merge,"tree state is invalid");

    f.print_size("f");
    f.reconstruct();
    double fnorm=f.norm2();
    double f1norm=f1.norm2();
    std::vector<real_function_2d> vf={f1,f2,f1};
    std::vector<double> vfnorm=norm2s(world,vf);
    real_function_2d ref;
    double norm=fnorm;

    std::vector<TreeState> states={reconstructed,				///< s coeffs at the leaves only
                                   compressed, 				///< d coeffs in internal nodes, s and d coeffs at the root
                                   nonstandard, 				///< s and d coeffs in internal nodes
                                   nonstandard_with_leaves, 	///< like nonstandard, with s coeffs at the leaves
                                   redundant};//,					///< s coeffs everywhere
    std::vector<TreeState> additional_states={redundant_after_merge};

    long k=FunctionDefaults<2>::get_k();



    // check on node coeffs size and return the norm
    auto check_nodes_have_coeffs = [](const real_function_2d& arg, const long k, const bool do_leaf) {
        bool correct_k=true;
        double norm=0.0;
        for (const auto& datum : arg.get_impl()->get_coeffs()) {
            const auto& node=datum.second;
            if (do_leaf == node.is_leaf()) {
                if (k>0) correct_k = correct_k and (node.has_coeff() and node.coeff().dim(0)==k);
                if (node.has_coeff()) norm+=std::pow(node.coeff().normf(),2);
            }
        }
        return std::make_pair(correct_k, sqrt(norm));
    };

    auto check_is_reconstructed = [&](const real_function_2d& arg) {
        auto [correct_k_leaf, norm_leaf]=check_nodes_have_coeffs(arg,k,true);
        auto [correct_k_interior, norm_interior]=check_nodes_have_coeffs(arg,0,false);
        bool correct_norm=(std::abs(norm_leaf-norm)<1.e-10) and (std::abs(norm_interior)<1.e-10);
        return correct_k_interior and correct_k_leaf and correct_norm and (arg.tree_size()==ref.tree_size());
    };

    auto check_is_compressed = [&](const real_function_2d& arg) {
        auto [correct_k_leaf, norm_leaf]=check_nodes_have_coeffs(arg,0,true);
        auto [correct_k_interior, norm_interior]=check_nodes_have_coeffs(arg,2*k,false);
        bool correct_norm=(std::abs(norm_leaf)<1.e-10) and (std::abs(norm_interior-norm)<1.e-10);
        return correct_k_interior and correct_k_leaf and correct_norm and (arg.tree_size()==ref.tree_size());
    };

    auto check_is_nonstandard = [&](const real_function_2d& arg) {
        auto [correct_k_leaf, norm_leaf]=check_nodes_have_coeffs(arg,0,true);
        auto [correct_k_interior, norm_interior]=check_nodes_have_coeffs(arg,2*k,false);
        bool correct_norm=norm_leaf<1.e-12;
        return correct_k_interior and correct_k_leaf and correct_norm and (arg.tree_size()==ref.tree_size());
    };

    auto check_is_nonstandard_with_leaves = [&](const real_function_2d& arg) {
        auto [correct_k_leaf, norm_leaf]=check_nodes_have_coeffs(arg,k,true);
        auto [correct_k_interior, norm_interior]=check_nodes_have_coeffs(arg,2*k,false);
        bool correct_norm=(std::abs(norm_leaf-norm)<1.e-10);
        return correct_k_interior and correct_k_leaf and correct_norm and (arg.tree_size()==ref.tree_size());
    };

    auto check_is_redundant = [&](const real_function_2d& arg) {
        auto [correct_k_leaf, norm_leaf]=check_nodes_have_coeffs(arg,k,true);
        auto [correct_k_interior, norm_interior]=check_nodes_have_coeffs(arg,k,false);
        bool correct_norm=(std::abs(norm_leaf-norm)<1.e-10);
        return correct_k_interior and correct_k_leaf and correct_norm and (arg.tree_size()==ref.tree_size());
    };

    auto check_tree_state = [&](const real_function_2d& arg, const TreeState state) {
        MADNESS_CHECK_THROW(arg.get_impl()->verify_tree_state_local(),"tree state is invalid");
        if (state==reconstructed) return check_is_reconstructed(arg);
        if (state==compressed) return check_is_compressed(arg);
        if (state==nonstandard) return check_is_nonstandard(arg);
        if (state==nonstandard_with_leaves) return check_is_nonstandard_with_leaves(arg);
        if (state==redundant) return check_is_redundant(arg);
        print("unknown state");
        return false;
    };

    auto vector_check_tree_state = [&](const std::vector<real_function_2d>& arg, const TreeState state) {
        for (int i=0; i<arg.size(); ++i) {
            ref=vf[i];
            norm=vfnorm[i];
            auto a=arg[i];
            bool ok=check_tree_state(a,state);
            if (not ok) return false;
        }
        return true;
    };

    print("f is reconstructed      ", check_is_reconstructed(f));
    print("f is compressed         ", check_is_compressed(f));
    print("f is nonstandard        ", check_is_nonstandard(f));
    print("f is nonstandard_leaves ", check_is_nonstandard_with_leaves(f));
    print("f is redundant          ", check_is_redundant(f));

    test_output t("testing tree state conversion");
    t.set_cout_to_terminal();

    // convert from all redundant_after_merge state to reconstructed
    for (const auto& ff : {f3}) {
        ref=f3ref;
        norm=f3norm;
        auto fcopy=copy(ff);
        std::stringstream ss_initial;
        ss_initial << get_tree_state(fcopy);
        fcopy=change_tree_state(fcopy,reconstructed);
        world.gop.fence();
        bool success=check_tree_state(fcopy,reconstructed);
        t.checkpoint(success,"conversion to reconstructed from "+ss_initial.str());
    }
    norm=fnorm;



    ref=f;
    // convert to all states from reconstructed
    for (auto state : states) {
        auto fcopy=copy(f);
        fcopy.change_tree_state(state);
        bool success=check_tree_state(fcopy,state);
        std::stringstream ss_initial;
        ss_initial << state;
        t.checkpoint(success,"conversion from reconstructed to "+ss_initial.str());
        print("f is",state, check_tree_state(fcopy,state));
    }

    for (auto initial_state : states) {
        auto fcopy=copy(f);
        fcopy.change_tree_state(initial_state);
        bool success=check_tree_state(fcopy,initial_state);
        std::stringstream ss_initial;
        ss_initial << initial_state;
        t.checkpoint(success,"initial conversion from reconstructed to "+ss_initial.str());
        for (auto final_state : states) {
            auto ffcopy=copy(fcopy);
            MADNESS_CHECK(ffcopy.get_impl()->get_tree_state()==initial_state);
            ffcopy.change_tree_state(final_state);
            success=check_tree_state(ffcopy,final_state);
            std::stringstream ss_final;
            ss_final << final_state;
            t.checkpoint(success,"conversion from "+ss_initial.str()+" to "+ss_final.str());
        }
    }


    print_header2("testing vector tree state conversion");
    // repeat for vectors of functions



    // convert to all states from reconstructed
    for (auto state : states) {
        auto fcopy=copy(world,vf);
        fcopy=change_tree_state(fcopy,state,false);
        world.gop.fence();
        bool success=vector_check_tree_state(fcopy,state);
        std::stringstream ss_initial;
        ss_initial << state;
        t.checkpoint(success,"conversion from reconstructed to "+ss_initial.str());
        print("f is",state, vector_check_tree_state(fcopy,state));
    }

    for (auto initial_state : states) {
        auto fcopy=copy(world,vf);
        change_tree_state(fcopy,initial_state);
        bool success=vector_check_tree_state(fcopy,initial_state);
        std::stringstream ss_initial;
        ss_initial << initial_state;
        t.checkpoint(success,"initial conversion from reconstructed to "+ss_initial.str());
        for (auto final_state : states) {
            auto ffcopy=copy(world,fcopy);
//            MADNESS_CHECK(ffcopy.get_impl()->get_tree_state()==initial_state);
            change_tree_state(ffcopy,final_state,false);
            world.gop.fence();
            success=vector_check_tree_state(ffcopy,final_state);
            std::stringstream ss_final;
            ss_final << final_state;
            t.checkpoint(success,"conversion from "+ss_initial.str()+" to "+ss_final.str());
        }
    }



    return t.end();


}


int main(int argc, char **argv) {
    madness::World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);
    FunctionDefaults<2>::set_thresh(1.e-6);
    FunctionDefaults<2>::set_k(6);
    FunctionDefaults<2>::set_cubic_cell(-20,20);

    int success = 0;

    success+=test_conversion(world);

    madness::finalize();
    return 0;
}
