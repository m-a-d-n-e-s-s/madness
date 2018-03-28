#ifndef LOADBALANCE_H
#define LOADBALANCE_H

#include "input.h"

// Load Balance
template <typename T, int NDIM>
struct LBCost
{
    double leaf_value;
    double parent_value;
    LBCost(double leaf_value=1.0, double parent_value=0.0)
        : leaf_value(leaf_value), parent_value(parent_value){}
    double operator()(const Key<NDIM>& key, const FunctionNode<T,NDIM>& node) const
    {
        if (key.level() <= 1){
        return 100.0*(leaf_value+parent_value);
        }
        else if (node.is_leaf()) { return leaf_value; }
        else { return parent_value; }
    }
};

// Load balance: three real functions
extern void loadbalance_xxp(World& world, real_functionT& p_den,
                                   real_functionT& n_den,
                                   real_functionT& den3);

// Load-balance: one vector of complex functions
extern void loadbalance_v1(World& world, comp_vecfuncT& u1);

// Load-balance: two vectors of complex functions
extern void loadbalance_v2(World& world,
                    comp_vecfuncT& u1,
                    comp_vecfuncT& u2);

// Load-balance: three vectors of complex functions
extern void loadbalance_v3(World& world,
                    comp_vecfuncT& u1,
                    comp_vecfuncT& u2,
                    comp_vecfuncT& u3);

// Load-balance: four vectors of complex functions
extern void loadbalance_v4(World& world,
                    comp_vecfuncT& u1,
                    comp_vecfuncT& u2,
                    comp_vecfuncT& u3,
                    comp_vecfuncT& u4);

// Load-balance: six vectors of complex functions
extern void loadbalance_v6(World& world,
                    comp_vecfuncT& u1,
                    comp_vecfuncT& u2,
                    comp_vecfuncT& u3,
                    comp_vecfuncT& v1,
                    comp_vecfuncT& v2,
                    comp_vecfuncT& v3);

#endif

