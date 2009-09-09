#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/qmprop.h>
#include <mra/operator.h>
#include <constants.h>
#include <tensor/vmath.h>
#include <mra/lbdeux.h>

using namespace madness;

template <int NDIM>
class MyCostFn {
private:
  static SharedPtr<ApplyTime<NDIM> > apply_time;
  static double wt_factor;
public:
//    MyCostFn(SharedPtr<ApplyTime<NDIM> > app_ptr = SharedPtr<ApplyTime<NDIM> > ()) :
//      apply_time(app_ptr)
//      , wt_factor(100.0) {
//      if (apply_time) {
//        // find smallest non-zero entry in apply_time
//        // set wt_factor to order of magnitude of that entry
//      }
//    }

  static void set_defaults() {
    apply_time = SharedPtr<ApplyTime<NDIM> > ();
    wt_factor = 100.0;
  }

  template <typename T>
  static Cost cost_fun(const Key<NDIM> & key, const FunctionNode<T,NDIM>& node) {
    return actual_cost_fun(key);
  }
  static Cost actual_cost_fun(const Key<NDIM> & key) {
    if (!(apply_time)) return 0;
    double f = apply_time->get(key);
    return (Cost) (f*wt_factor);
  }
//   static long cost_fun(const Key<NDIM> & key, const FunctionNode<double_complex,NDIM>& node) {
//     return cost_fun(key);
//   }
//   static long cost_fun(const Key<NDIM> & key, const FunctionNode<double,NDIM>& node) {
//     return cost_fun(key);
//   }

  static void set_apply_time(SharedPtr<ApplyTime<NDIM> > apt) {
    apply_time = apt;
  }

  static void set_wt_factor() {
    // find smallest non-zero entry in apply_time
    // set wt_factor to order of magnitude of that entry
  }
  static void set_wt_factor(double w) {
    wt_factor = w;
  }

};
