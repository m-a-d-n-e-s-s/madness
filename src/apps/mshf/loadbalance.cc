
//////////////////////////////////////////////////////////////////////
////////////////////// DATA LOAD BALANCING ///////////////////////////
////////////////////// see examples/dataloadbal.cc for details ///////
//////////////////////////////////////////////////////////////////////

#include "loadbalance.h"

// Load balance: three real functions
void loadbalance_xxp(World& world, real_functionT& p_den,
                                   real_functionT& n_den,
                                   real_functionT& den3)
{   
    LoadBalanceDeux<3> lb(world);
    lb.add_tree(p_den, LBCost<double,3>(1.0,1.0), true);
    world.gop.fence(); 
    lb.add_tree(n_den, LBCost<double,3>(1.0,1.0), true);
    world.gop.fence(); 
    lb.add_tree(den3,  LBCost<double,3>(1.0,1.0), true);
    world.gop.fence();
    FunctionDefaults<3>::redistribute(world, lb.load_balance(2.0, false));
}


// Load-balance: one vector of complex functions
void loadbalance_v1(World& world, comp_vecfuncT& u1)
{   
    LoadBalanceDeux<3> lb(world);
    for (unsigned int i=0; i<u1.size(); i++) {
        lb.add_tree(u1[i], LBCost<double_complex,3>(1.0,1.0), true);
    }
    world.gop.fence();
    FunctionDefaults<3>::redistribute(world, lb.load_balance(2.0, false));
}


// Load-balance: two vectors of complex functions
void loadbalance_v2(World& world,
                    comp_vecfuncT& u1,
                    comp_vecfuncT& u2)
{   
    LoadBalanceDeux<3> lb(world);
    for (unsigned int i=0; i<u1.size(); i++) {
        lb.add_tree(u1[i], LBCost<double_complex,3>(1.0,1.0), true);
    }
    world.gop.fence();
    
    for (unsigned int i=0; i<u2.size(); i++) {
        lb.add_tree(u2[i], LBCost<double_complex,3>(1.0,1.0), true);
    }
    world.gop.fence();
    FunctionDefaults<3>::redistribute(world, lb.load_balance(2.0, false));
}


// Load-balance: three vectors of complex functions
void loadbalance_v3(World& world,
                    comp_vecfuncT& u1,
                    comp_vecfuncT& u2,
                    comp_vecfuncT& u3)
{   
    LoadBalanceDeux<3> lb(world);
    for (unsigned int i=0; i<u1.size(); i++) {
        lb.add_tree(u1[i], LBCost<double_complex,3>(1.0,1.0), true);
    }
    world.gop.fence();
    for (unsigned int i=0; i<u2.size(); i++) {
        lb.add_tree(u2[i], LBCost<double_complex,3>(1.0,1.0), true);
    }
    world.gop.fence();
    for (unsigned int i=0; i<u3.size(); i++) {
        lb.add_tree(u3[i], LBCost<double_complex,3>(1.0,1.0), true);
    }
    world.gop.fence();
    FunctionDefaults<3>::redistribute(world, lb.load_balance(2.0, false));
}


// Load-balance: four vectors of complex functions 
void loadbalance_v4(World& world,
                    comp_vecfuncT& u1,
                    comp_vecfuncT& u2,
                    comp_vecfuncT& u3,
                    comp_vecfuncT& u4)
{   
    LoadBalanceDeux<3> lb(world);
    for (unsigned int i=0; i<u1.size(); i++) {
        lb.add_tree(u1[i], LBCost<double_complex,3>(1.0,1.0), true);
    }
    world.gop.fence();
    for (unsigned int i=0; i<u2.size(); i++) {
        lb.add_tree(u2[i], LBCost<double_complex,3>(1.0,1.0), true);
    }
    world.gop.fence();
    for (unsigned int i=0; i<u3.size(); i++) {
        lb.add_tree(u3[i], LBCost<double_complex,3>(1.0,1.0), true);
    }
    world.gop.fence();
    for (unsigned int i=0; i<u4.size(); i++) {
        lb.add_tree(u4[i], LBCost<double_complex,3>(1.0,1.0), true);
    }
    world.gop.fence();
    FunctionDefaults<3>::redistribute(world, lb.load_balance(2.0, false));
}


// Load-balance: six vectors of complex functions
void loadbalance_v6(World& world,
                    comp_vecfuncT& u1,
                    comp_vecfuncT& u2,
                    comp_vecfuncT& u3,
                    comp_vecfuncT& v1,
                    comp_vecfuncT& v2,
                    comp_vecfuncT& v3)
{   
    LoadBalanceDeux<3> lb(world);
    for (unsigned int i=0; i<u1.size(); i++) {
        lb.add_tree(u1[i], LBCost<double_complex,3>(1.0,1.0), true);
    }
    world.gop.fence();
    for (unsigned int i=0; i<u2.size(); i++) {
        lb.add_tree(u2[i], LBCost<double_complex,3>(1.0,1.0), true);
    }
    world.gop.fence();
    for (unsigned int i=0; i<u3.size(); i++) {
        lb.add_tree(u3[i], LBCost<double_complex,3>(1.0,1.0), true);
    }
    world.gop.fence();
    for (unsigned int i=0; i<v1.size(); i++) {
        lb.add_tree(v1[i], LBCost<double_complex,3>(1.0,1.0), true);
    }
    world.gop.fence();
    for (unsigned int i=0; i<v2.size(); i++) {
        lb.add_tree(v2[i], LBCost<double_complex,3>(1.0,1.0), true);
    }
    world.gop.fence();
    for (unsigned int i=0; i<v3.size(); i++) {
        lb.add_tree(v3[i], LBCost<double_complex,3>(1.0,1.0), true);
    }
    world.gop.fence();
    FunctionDefaults<3>::redistribute(world, lb.load_balance(2.0, false));
}

