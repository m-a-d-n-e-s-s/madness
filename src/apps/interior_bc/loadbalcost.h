// load balancing structure lifted from dataloadbal.cc
struct LBCost {
    double leaf_value;
    double parent_value;
    LBCost(double leaf_value = 1.0, double parent_value = 1.0)
        : leaf_value(leaf_value), parent_value(parent_value) {}

    double operator() (const Key<3> &key, const FunctionNode<double, 3> &node)
        const {

        if(key.level() <= 1) {
            return 100.0*(leaf_value+parent_value);
        }
        else if(node.is_leaf()) {
            return leaf_value;
        }
        else {
            return parent_value;
        }
    }
};
