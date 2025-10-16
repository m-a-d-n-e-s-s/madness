/**
 \file test_vectormacrotask.cc
 \brief Tests the \c MacroTaskQ and MacroTask classes
 \ingroup mra

 For more information also see macrotask.h

 The user defines a macrotask (an example is found in test_vectormacrotask.cc), the tasks are
 lightweight and carry only bookkeeping information, actual input and output are stored in a
 cloud (see cloud.h)

 The user-defined MacroTask is derived from MacroTaskOperationBase and must implement the call operator
 method. A heterogeneous task queue is possible.

*/

#include <madness/mra/mra.h>
#include <madness/world/cloud.h>
#include <madness/mra/macrotaskq.h>
#include <madness/world/world.h>
#include <madness/world/timing_utilities.h>
#include <madness/mra/macrotaskpartitioner.h>

#include "test_utilities.h"


using namespace madness;
using namespace archive;

struct slater {
    double a=1.0;

    slater(double aa) : a(aa) {}

    template<std::size_t NDIM>
    double operator()(const Vector<double,NDIM> &r) const {
        double x = inner(r,r);
        return exp(-a *sqrt(x));
    }
};

struct gaussian {
    double a;

    gaussian() : a() {};

    gaussian(double aa) : a(aa) {}

    template<std::size_t NDIM>
    double operator()(const Vector<double,NDIM> &r) const {
        double r2 = inner(r,r);
        return exp(-a * r2);//*abs(sin(abs(2.0*x))) *cos(y);
    }
};



template<typename T, std::size_t NDIM>
class MicroTask : public MacroTaskOperationBase {
public:
    // you need to define the exact argument(s) of operator() as tuple
    typedef std::tuple<const Function<T,NDIM> &, const double &,
            const std::vector<Function<T,NDIM>> &> argtupleT;

    // you need to define the result type
    // resultT must implement gaxpy(alpha, result, beta, contribution)
    // with resultT result, contribution;
    using resultT = std::vector<Function<T,NDIM>>;

    // you need to define an empty constructor for the result
    resultT allocator(World &world, const argtupleT &argtuple) const {
        std::size_t n = std::get<2>(argtuple).size();
        resultT result = zero_functions_auto_tree_state<T,NDIM>(world, n);
        return result;
    }

    resultT operator()(const Function<T,NDIM> &f1, const double &arg2,
                       const std::vector<Function<T,NDIM>> &f2) const {
        World& world=f1.world();
        if (world.rank()==0) print("in MicroTask");
        SeparatedConvolution<double,NDIM> op=BSHOperator<NDIM>(world,0,1.e-4,1.e-5);
        return arg2 * f1 * apply(world,op,f2);
    }
};


class MicroTask1 : public MacroTaskOperationBase{
public:
    // you need to define the result type
    // resultT must implement gaxpy(alpha, result, beta, contribution)
    // with resultT result, contribution;
    typedef real_function_3d resultT;

    // you need to define the exact argument(s) of operator() as tuple
    typedef std::tuple<const real_function_3d &, const double &,
            const std::vector<real_function_3d> &> argtupleT;

    resultT allocator(World &world, const argtupleT &argtuple) const {
        resultT result = real_factory_3d(world).compressed();
        return result;
    }

    resultT operator()(const real_function_3d &f1, const double &arg2,
                       const std::vector<real_function_3d> &f2) const {
        auto result = arg2 * f1 * inner(f2, f2);
        return result;
    }
};

class MicroTask2 : public MacroTaskOperationBase{
public:
    // you need to define the result type
    // resultT must implement gaxpy(alpha, result, beta, contribution)
    // with resultT result, contribution;
    typedef std::vector<real_function_3d> resultT;

    // you need to define the exact argument(s) of operator() as tuple
    typedef std::tuple<const std::vector<real_function_3d> &, const double &,
            const std::vector<real_function_3d> &> argtupleT;

    resultT allocator(World &world, const argtupleT &argtuple) const {
        std::size_t n = std::get<2>(argtuple).size();
        resultT result = zero_functions_compressed<double, 3>(world, n);
        return result;
    }

    resultT operator()(const std::vector<real_function_3d>& f1, const double &arg2,
                       const std::vector<real_function_3d>& f2) const {
        // World &world = f1[0].world();
        // // won't work because of nested loop over f1
        // if (batch.input[0]==batch.input[1]) {
        //     return arg2 * f1 * inner(f1,f2);
        // }

        // // won't work because result batches are the same as f1 batches
        // return arg2*f2*inner(f1,f1);

        // will work
        return arg2*f1*inner(f2,f2);

    }
};

class VectorOfScalarTask : public MacroTaskOperationBase{
public:
    // you need to define the result type
    // resultT must implement gaxpy(alpha, result, beta, contribution)
    // with resultT result, contribution;
    typedef std::vector<ScalarResult<double>> resultT;

    // you need to define the exact argument(s) of operator() as tuple
    typedef std::tuple<const std::vector<real_function_3d> &, const double &,
            const std::vector<real_function_3d> &> argtupleT;

    resultT allocator(World &world, const argtupleT &argtuple) const {
        std::size_t n = std::get<0>(argtuple).size();
        return scalar_result_vector<double>(world,n);
    }

    resultT operator()(const std::vector<real_function_3d>& f1, const double &arg2,
                       const std::vector<real_function_3d>& f2) const {
        World &world = f1[0].world();
        auto result=scalar_result_vector<double>(world,f1.size());
        for (int i=0; i<f1.size(); ++i) result[i]=double(i+batch.input[0].begin);

        return result;
    }
};

class ScalarTask : public MacroTaskOperationBase{
public:
    // you need to define the result type
    // resultT must implement gaxpy(alpha, result, beta, contribution)
    // with resultT result, contribution;
    // typedef std::shared_ptr<ScalarResultImpl<double>> resultT;
    typedef ScalarResult<double> resultT;

    // you need to define the exact argument(s) of operator() as tuple
    typedef std::tuple<const std::vector<real_function_3d> &> argtupleT;

    resultT allocator(World &world, const argtupleT &argtuple) const {
        return resultT(world);
    }

    resultT operator()(const std::vector<real_function_3d>& f1) const {
        World &world = f1[0].world();
        resultT result(world);
        result=double(f1.size());
        return result;
    }
};

class TupleTask : public MacroTaskOperationBase{
public:
    // you need to define the result type
    // resultT must implement gaxpy(alpha, result, beta, contribution)
    // with resultT result, contribution;
    typedef std::tuple<vector_real_function_3d,vector_real_function_3d> resultT;

    // you need to define the exact argument(s) of operator() as tuple
    typedef std::tuple<const std::vector<real_function_3d> &> argtupleT;

    resultT allocator(World &world, const argtupleT &argtuple) const {
        std::size_t n=std::get<0>(argtuple).size();
        auto v1=zero_functions_compressed<double,3>(world,n);
        auto v2=zero_functions_compressed<double,3>(world,n);
        return std::make_tuple(v1,v2);
    }

    resultT operator()(const std::vector<real_function_3d>& f1) const {
        World &world = f1[0].world();
        auto result1=copy(world,f1);
        auto result2=f1+f1;
        auto result=std::make_tuple(result1,result2);
        return result;
    }
};

/// this task won't do anything, is mostly to check if the combinations compile
template<typename elementT, typename elementR>
class MixedTupleTask : public MacroTaskOperationBase{
public:
    // you need to define the result type
    // resultT must implement gaxpy(alpha, result, beta, contribution)
    // with resultT result, contribution;
    typedef std::tuple<elementT,elementR> resultT;

    // you need to define the exact argument(s) of operator() as tuple
    typedef std::tuple<const std::vector<real_function_3d> &> argtupleT;

    // some allocators
    template<typename T>
    typename std::enable_if<std::is_same<T, Function<double,3>>::value, T>::type
    allocator(World& world, std::size_t n) const {
        return Function<double,3>(world);
    }

    template<typename T>
    typename std::enable_if<std::is_same<T, std::vector<real_function_3d>>::value, T>::type
    allocator(World& world, std::size_t n) const {
        return zero_functions_compressed<double,3>(world,n);
    }

    template<typename T>
    typename std::enable_if<std::is_same<T, std::vector<ScalarResult<double>>>::value, T>::type
    allocator(World& world, std::size_t n) const {
        return scalar_result_vector<double>(world,n);
    }

    template<typename T>
    typename std::enable_if<std::is_same<T, ScalarResult<double>>::value, T>::type
    allocator(World& world, std::size_t n) const {
        return ScalarResult<double>(world);
    }

    resultT allocator(World &world, const argtupleT &argtuple) const {
        std::size_t n=std::get<0>(argtuple).size();
        auto v1=allocator<elementT>(world,n);
        auto v2=allocator<elementR>(world,n);
        return std::make_tuple(v1,v2);
    }

    resultT operator()(const std::vector<real_function_3d>& f1) const {
        World &world = f1[0].world();
        std::size_t n=f1.size();
        auto v1=allocator<elementT>(world,n);
        auto v2=allocator<elementR>(world,n);
        return std::make_tuple(v1,v2);
    }
};



template<typename T, std::size_t NDIM>
int check_vector(World& universe, const std::vector<Function<T,NDIM>> &ref,
                    const std::vector<Function<T,NDIM>> &test,
                    const std::string msg) {
    double norm_ref = norm2(universe, ref);
    double norm_test = norm2(universe, test);
    double error = norm2(universe, ref - test);
    bool success = error/norm_ref < 1.e-10;
    if (universe.rank() == 0) {
        print("norm ref, test, diff", norm_ref, norm_test, error);
        if (success) print("test", msg, " \033[32m", "passed ", "\033[0m");
        else print("test", msg, " \033[31m", "failed \033[0m ");
    }
    return (success) ? 0 : 1;
};


int check(World& universe, const double &ref, const double &test, const std::string msg) {
    double error = (ref - test);
    bool success = error/ref < 1.e-10;
    if (universe.rank() == 0) {
                print("norm ref, test, diff", ref, test, error);
        if (success) print("test", msg, " \033[32m", "passed ", "\033[0m");
        else print("test", msg, " \033[31m", "failed \033[0m ");
    }
    return (success) ? 0 : 1;
}

int check(World& universe, const real_function_3d &ref, const real_function_3d &test, const std::string msg) {
    double norm_ref = ref.norm2();
    double norm_test = test.norm2();
    double error = (ref - test).norm2();
    bool success = error/norm_ref < 1.e-10;
    if (universe.rank() == 0) {
                print("norm ref, test, diff", norm_ref, norm_test, error);
        if (success) print("test", msg, " \033[32m", "passed ", "\033[0m");
        else print("test", msg, " \033[31m", "failed \033[0m ");
    }
    return (success) ? 0 : 1;
};

template<typename T, std::size_t NDIM>
int test_immediate(World& universe, const std::vector<Function<T,NDIM>>& v3,
                   const std::vector<Function<T,NDIM>>& ref) {
    test_output t1("testing immediate execution",universe.rank()==0);
    // t1.set_cout_to_terminal();
    for (std::string preset : MacroTaskInfo::get_all_preset_names()) {
        MicroTask<T,NDIM> t;
        MacroTask task_immediate(universe, t, MacroTaskQFactory(universe).preset(preset));
        task_immediate.set_debug(true);
        std::vector<Function<T,NDIM>> v = task_immediate(v3[0], 2.0, v3);
        double error=norm2(universe, v-ref);
        std::stringstream ss;
        ss << "immediate execution -- storage policy " << preset;
        t1.checkpoint(error,1.e-10,ss.str());
    }
    return t1.end();
}

int test_deferred(World& universe, const std::vector<real_function_3d>& v3,
                   const std::vector<real_function_3d>& ref) {
    test_output t1("testing deferred execution",universe.rank()==0);
    t1.set_cout_to_terminal();
    for (const auto& policy_name : MacroTaskInfo::get_all_preset_names()) {
        auto policy=MacroTaskInfo::preset(policy_name);
        print(policy);
        auto taskq = std::make_shared<MacroTaskQ>(MacroTaskQFactory(universe).set_policy(policy));
        taskq->set_printlevel(20);
        MicroTask<double,3> t;
        MacroTask task(universe, t, taskq);
        task.set_debug(true);
        std::vector<real_function_3d> f2a = task(v3[0], 2.0, v3);
        taskq->print_taskq();
        // auto [global_size,global_memsize,min_memsize,max_memsize,max_record_size]=taskq->cloud.get_size(universe);
        nlohmann::json stats=taskq->cloud.gather_memory_statistics(universe);
        double global_memsize=stats["memory_global"].template get<double>();


        if (policy.storage_policy==MacroTaskInfo::StoreFunction) {
            t1.checkpoint(global_memsize>1.e3,"cloud global memory size > 1 kB");
        } else if (policy.storage_policy==MacroTaskInfo::StoreFunctionViaPointer) {
            t1.checkpoint(global_memsize<1.e3,"cloud global memory size < 1 kB");
        }
        taskq->run_all();
        taskq->cloud.print_timings(taskq->cloud.gather_timings(universe));
        taskq->cloud.clear_timings();
        double error=norm2(universe, f2a-ref);
        std::stringstream ss;
        ss << "deferred execution -- storage policy " << policy_name;
        t1.checkpoint(error,1.e-9,ss.str());
    }
    return t1.end();
}

int test_twice(World& universe, const std::vector<real_function_3d>& v3,
                  const std::vector<real_function_3d>& ref) {
    if (universe.rank() == 0) print("\nstarting Microtask twice (check caching)\n");
    MacroTaskInfo largemem=MacroTaskInfo::preset("large_memory");
    auto taskq = std::shared_ptr<MacroTaskQ>(new MacroTaskQ(MacroTaskQFactory(universe).set_policy(largemem)));
    taskq->set_printlevel(3);
    MicroTask<double,3> t;
    MacroTask task(universe, t, taskq);
    std::vector<real_function_3d> f2a1 = task(v3[0], 2.0, v3);
    std::vector<real_function_3d> f2a2 = task(v3[0], 2.0, v3);
    taskq->print_taskq();
    taskq->run_all();
    taskq->cloud.print_timings(taskq->cloud.gather_timings(universe));
    int success=0;
    success += check_vector(universe, ref, f2a1, "taske twice a");
    success += check_vector(universe, ref, f2a2, "taske twice b");
    return success;
}

int test_task1(World& universe, const std::vector<real_function_3d>& v3) {
    if (universe.rank()==0) print("\nstarting Microtask1\n");
    MicroTask1 t1;
    real_function_3d ref_t1 = t1(v3[0], 2.0, v3);
    MacroTask task1(universe, t1);
    task1.set_debug(true);
    real_function_3d ref_t2 = task1(v3[0], 2.0, v3);
    int success = check(universe,ref_t1, ref_t2, "task1 immediate");
    return success;
}

/// each task accumulates into the same result
int test_scalar_task(World& universe, const std::vector<real_function_3d>& v3) {
    if (universe.rank()==0) print("\nstarting ScalarTask\n");
    ScalarTask t1;
    ScalarResult<double> result = t1(v3);
    double ref_t1=result.get();
    print("result reference",ref_t1);

    MacroTask task1(universe, t1);
    auto result2= task1(v3);
    double result_t1=result2.get();
    print("result macro",result_t1);

    int success = check(universe,ref_t1, result_t1, "ScalarTask");
    return success;
}


int test_vector_of_scalar_task(World& universe, const std::vector<real_function_3d>& v3) {
    if (universe.rank()==0) print("\nstarting VectorOfScalarTask\n");
    VectorOfScalarTask t1;
    std::vector<ScalarResult<double>> result = t1(v3, 2.0, v3);

    MacroTask task1(universe, t1);
    auto result2= task1(v3, 2.0, v3);

    double error=1.e-15;
    for (int i=0; i<result.size(); ++i) error+=fabs(result[i].get()-result2[i].get());
    int success = check(universe,1.e-15, error, "vector of scalar task");
    return success;
}

/// each task accumulates into the same result
int test_tuple_of_vectors(World& universe, const std::vector<real_function_3d>& v3) {
    if (universe.rank()==0) print("\nstarting TupleOfVectorsTask\n");
    TupleTask t1;
    auto [ref1,ref2]=t1(v3);
    double n1=norm2(universe,ref1);
    double n2=norm2(universe,ref2);
    print("result reference",n1,n2);


    MacroTask task1(universe, t1);
    auto [res1,res2] = task1(v3);
    double m1=norm2(universe,res1);
    double m2=norm2(universe,res2);
    print("result macro",m1,m2);

    int success = check(universe,n1,m1, "task1 immediate");
    return success;
}


// this will only test compilation
int test_mixed_tuple(World& universe, const std::vector<real_function_3d>& v3) {
    if (universe.rank()==0) print("\nstarting MixedTupleTask\n");

    typedef real_function_3d type1;
    typedef std::vector<real_function_3d> type2;
    typedef ScalarResult<double> type3;
    typedef std::vector<ScalarResult<double>> type4;
    MixedTupleTask<type1,type1> t11;
    MixedTupleTask<type1,type2> t12;
    MixedTupleTask<type1,type3> t13;
    MixedTupleTask<type1,type4> t14;
    MixedTupleTask<type2,type1> t21;
    MixedTupleTask<type2,type2> t22;
    MixedTupleTask<type2,type3> t23;
    MixedTupleTask<type2,type4> t24;
    MixedTupleTask<type3,type1> t31;
    MixedTupleTask<type3,type2> t32;
    MixedTupleTask<type3,type3> t33;
    MixedTupleTask<type3,type4> t34;
    MixedTupleTask<type4,type1> t41;
    MixedTupleTask<type4,type2> t42;
    MixedTupleTask<type4,type3> t43;
    MixedTupleTask<type4,type4> t44;
    auto result11=t11(v3);
    auto result12=t12(v3);
    auto result13=t13(v3);
    auto result14=t14(v3);
    auto result21=t21(v3);
    auto result22=t22(v3);
    auto result23=t23(v3);
    auto result24=t24(v3);
    auto result31=t31(v3);
    auto result32=t32(v3);
    auto result33=t33(v3);
    auto result34=t34(v3);
    auto result41=t41(v3);
    auto result42=t42(v3);
    auto result43=t43(v3);
    auto result44=t44(v3);
    return 0;
}

int test_2d_partitioning(World& universe, const std::vector<real_function_3d>& v3) {
    if (universe.rank() == 0) print("\nstarting 2d partitioning");
    auto policy=MacroTaskInfo::preset("large_memory");
    auto taskq = std::shared_ptr<MacroTaskQ>(new MacroTaskQ(MacroTaskQFactory(universe).set_policy(policy)));
    taskq->set_printlevel(3);
    MicroTask2 t;
    auto ref=t(v3,2.0,v3);
    t.partitioner->set_dimension(2);
    MacroTask task(universe, t, taskq);
    std::vector<real_function_3d> f2a = task(v3, 2.0, v3);
    taskq->print_taskq();
    taskq->run_all();
    taskq->cloud.print_timings(taskq->cloud.gather_timings(universe));
    taskq->cloud.clear_timings();
    int success=check_vector(universe,ref,f2a,"test 2d partitioning");
    return success;
}

int main(int argc, char **argv) {
    madness::World &universe = madness::initialize(argc, argv);
    startup(universe, argc, argv, true);
    FunctionDefaults<3>::set_thresh(1.e-5);
    FunctionDefaults<3>::set_k(9);
    FunctionDefaults<3>::set_cubic_cell(-20,20);

    int success = 0;

    universe.gop.fence();
    int nworld = universe.size();
    if (universe.rank() == 0) print("creating nworld", nworld, universe.id());

    {
        // execution in a taskq, result will be complete only after the taskq is finished
        real_function_3d f1 = real_factory_3d(universe).functor(slater(1.0));
        real_function_3d i2 = real_factory_3d(universe).functor(slater(2.0));
        real_function_3d i3 = real_factory_3d(universe).functor(slater(2.0));
        std::vector<real_function_3d> v3;
        for (int i=0; i<20; ++i) v3.push_back(real_factory_3d(universe).functor(slater(sqrt(double(i)))));

        // test a low-rank vector
        FunctionDefaults<2>::set_tensor_type(TT_2D);
        std::vector<real_function_2d> v2;
        for (int i=0; i<20; ++i) v2.push_back(real_factory_2d(universe).functor(slater(sqrt(double(i)))));

        timer timer1(universe);
        MicroTask<double,3> t;
        std::vector<real_function_3d> ref3 = t(v3[0], 2.0, v3);
        timer1.tag("direct execution");

        MicroTask<double,2> t2;
        std::vector<real_function_2d> ref2 = t2(v2[0], 2.0, v2);
        timer1.tag("direct execution low-rank");

        // success+=test_immediate<double,3>(universe,v3,ref3);
        // timer1.tag("immediate taskq execution");

        success+=test_immediate<double,2>(universe,v2,ref2);
        timer1.tag("immediate taskq execution with low-rank tensors");

        success+=test_scalar_task(universe,v3);
        timer1.tag("scalar task execution");

        success+=test_vector_of_scalar_task(universe,v3);
        timer1.tag("vector of scalar task execution");

        success+=test_tuple_of_vectors(universe,v3);
        timer1.tag("vector of tuples task execution");

        success+=test_mixed_tuple(universe,v3);
        timer1.tag("mixed tuple task execution");

        success+=test_deferred(universe,v3,ref3);
        timer1.tag("deferred taskq execution");

        success+=test_twice(universe,v3,ref3);
        timer1.tag("executing a task twice");

        success+=test_task1(universe,v3);
        timer1.tag("task1 immediate execution");

        success+=test_2d_partitioning(universe,v3);
        timer1.tag("2D partitioning");

        if (universe.rank() == 0) {
            if (success==0) print("\n --> all tests \033[32m", "passed ", "\033[0m\n");
            else print("\n --> all tests \033[31m", "failed \033[0m \n");
        }
    }

    madness::finalize();
    return success;
}

template<> volatile std::list<detail::PendingMsg> WorldObject<MacroTaskQ>::pending = std::list<detail::PendingMsg>();
template<> Spinlock WorldObject<MacroTaskQ>::pending_mutex(0);

//template <> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, MacroTask, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
//template <> Spinlock WorldObject<WorldContainerImpl<long, MacroTask, madness::Hash<long> > >::pending_mutex(0);

template<> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
template<> Spinlock WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending_mutex(
        0);
