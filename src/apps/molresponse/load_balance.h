#ifndef SRC_APPS_molresponse_LOAD_BALANCE_H
#define SRC_APPS_molresponse_LOAD_BALANCE_H
#include "TDDFT.h"
// Needed for rebalancing
template <typename T, int NDIM>
struct lbcost {
  double leaf_value;
  double parent_value;
  explicit lbcost(double leaf_value = 1.0, double parent_value = 0.0)
      : leaf_value(leaf_value), parent_value(parent_value) {}
  double operator()(const Key<NDIM>& key, const FunctionNode<T, NDIM>& node) const {
    if (key.level() < 1) {
      return 100.0 * (leaf_value + parent_value);
    } else if (node.is_leaf()) {
      return leaf_value;
    } else {
      return parent_value;
    }
  }
};

#endif
