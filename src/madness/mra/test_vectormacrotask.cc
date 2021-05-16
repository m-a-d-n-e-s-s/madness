/**
 \file test_vectormacrotask.h
 \brief Tests the \c MacroTaskQ and MacroTask classes
 \ingroup mra

 For more information also see macrotask.h

 The user defines a macrotask (an example is found in test_vectormacrotask.cc), the tasks are
 lightweight and carry only bookkeeping information, actual input and output are stored in a
 cloud (see cloud.h)

 The user-defined MacroTask is derived from MacroTaskIntermediate and must implement the run()
 method. A heterogeneous task queue is possible.

 When running the taskq with Functions make sure the process map is reset before and after the
 task execution:

 	    FunctionDefaults<NDIM>::set_default_pmap(taskq.get_subworld());
	    taskq.run_all(task_list);
		FunctionDefaults<NDIM>::set_default_pmap(universe);

 Also make sure the subworld objects are properly destroyed after the execution

*/

/*
 * for_each i
 *      result[i] = op(arg[i],other_args)
 *
 * for_each i
 *      result[j] = op(arg[i],other_args)
 *
 */



#include <madness/mra/mra.h>
#include <madness/world/cloud.h>
#include <madness/mra/macrotaskq.h>
#include <madness/world/world.h>
#include <madness/world/timing_utilities.h>


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


static std::atomic<int> atomic_idx;


template<typename ... Ts>
constexpr auto decay_types(std::tuple<Ts...> const &)
-> std::tuple<std::remove_cv_t<std::remove_reference_t<Ts>>...>;

template<typename T>
using decay_tuple = decltype(decay_types(std::declval<T>()));

template<typename Q>
struct is_vector : std::false_type {
};
template<typename Q>
struct is_vector<std::vector<Q>> : std::true_type {
};

template<typename>
struct is_madness_function_vector : std::false_type {
};

template<typename T, std::size_t NDIM>
struct is_madness_function_vector<std::vector<typename madness::Function<T, NDIM>>> : std::true_type {
};


template<size_t I = 0, typename tupleT>
constexpr std::size_t get_index_of_first_vector_argument() {

    typedef decay_tuple<tupleT> argtupleT;   // removes const, &, etc

    if constexpr(I >= std::tuple_size_v<tupleT>) {
        // Last case, if nothing is left to iterate, then exit the function
        MADNESS_EXCEPTION("there is no madness function vector argument in the list, cannot partition the tasks",1);
        return I;
    } else {
        // Print the tuple and go to next element
        using typeT = typename std::tuple_element<I,argtupleT>::type;       // use decay types for determining a vector
        if constexpr (is_madness_function_vector<typeT>::value) {
            return I;
        } else {
            // Going for next element.
            return get_index_of_first_vector_argument<I + 1,tupleT>();
        }
    }
}


class Batch {
    friend class OrbitalPartitioner;
public:
    std::size_t begin, end; ///< first and first past last index [begin,end)

    Batch(const std::size_t& begin, const std::size_t& end) : begin(begin), end(end) {}

    /// select the relevant vector elements from the argument tuple
    template<typename tupleT>
    tupleT operator()(const tupleT &arg) const {

        constexpr std::size_t index=get_index_of_first_vector_argument<0,tupleT>();
        auto v=std::get<index>(arg);
        decltype(v) v_batch;
        std::copy(v.begin()+begin,v.begin()+end,std::back_inserter(v_batch));
        tupleT batched_arg=arg;
        std::get<index>(batched_arg)=v_batch;
        return batched_arg;
    }

    friend std::ostream& operator<<(std::ostream& os, const Batch& batch) {
        std::stringstream ss;
        ss << "[" << std::setw(3) << batch.begin << ", " << batch.end << ")";
        os << ss.str();
        return os;
    }

};

class OrbitalPartitioner {
    friend class Batch;
public:
    typedef std::list<Batch> partitionT;
    std::size_t min_batch_size;
    std::size_t max_batch_size=10;
    std::size_t nsubworld;
    std::string policy="guided";

    OrbitalPartitioner(long nsubworld, long min_batch_size=5, std::string policy="guided")
            : min_batch_size(min_batch_size), nsubworld(nsubworld), policy(policy) {}

    template<typename tupleT>
    partitionT partition_tasks(const tupleT &argtuple) const {

        auto v=get_first_vector_argument(argtuple);
        partitionT result;

        if (policy == "guided") {
            std::size_t begin=0;
            std::size_t end=0;
            while (end<v.size()) {
                end+=std::min(max_batch_size,std::max(min_batch_size,((v.size()-end)/nsubworld)));
                end=std::min(end,v.size());
                result.push_back(Batch(begin,end));
                begin=end;
            }
        } else {
            std::string msg="unknown partitioning policy: "+policy;
            MADNESS_EXCEPTION(msg.c_str(),1);
        }
        return result;
    }


    /// find the first madness function vector argument in tuple and return it
    template<size_t I = 0, typename... Ts>
    static constexpr auto get_first_vector_argument(std::tuple<Ts...> argtuple) {

        typedef decay_tuple<typename std::tuple<Ts...>> argtupleT;   // removes const, &, etc

        if constexpr(I == sizeof...(Ts)) {
            // Last case, if nothing is left to iterate, then exit the function
            MADNESS_EXCEPTION("there is no madness function vector argument in the list, cannot partition the tasks",1);
            return double(0.);
        } else {
            // Print the tuple and go to next element
            using typeT = typename std::tuple_element<I,argtupleT>::type;       // use decay types for determining a vector
            if constexpr (is_madness_function_vector<typeT>::value) {
                return std::get<I>(argtuple);
            } else {
                // Going for next element.
                return get_first_vector_argument<I + 1>(argtuple);
            }
        }
    }
};


template<typename taskT>
class MacroTask_2G {
    using partitionT = OrbitalPartitioner::partitionT;
    typedef typename taskT::resultT resultT;
    typedef typename taskT::argtupleT argtupleT;
    typedef Cloud::recordlistT recordlistT;
    taskT task;

public:

    MacroTask_2G(World &world, taskT &task, std::shared_ptr<MacroTaskQ> taskq_ptr = 0)
            : world(world), task(task), taskq_ptr(taskq_ptr) {}

    template<typename ... Ts>
    resultT operator()(const Ts &... args) {

        const bool immediate_execution = (not taskq_ptr);
        if (not taskq_ptr) taskq_ptr.reset(new MacroTaskQ(world, world.size()));

        auto argtuple = std::tie(args...);
        static_assert(std::is_same<decltype(argtuple), argtupleT>::value, "type or number of arguments incorrect");

        // partition the argument vector into batches
        OrbitalPartitioner op(taskq_ptr->get_nsubworld());
        partitionT partition = op.partition_tasks(argtuple);

        // store input and output: output being a pointer to a universe function (vector)
        recordlistT inputrecords = taskq_ptr->cloud.store(world, argtuple);
        auto[outputrecords, result] =prepare_output(taskq_ptr->cloud, argtuple);

        // create tasks and add them to the taskq
        MacroTaskBase::taskqT vtask;
        for (const Batch &batch : partition) {
            vtask.push_back(
                    std::shared_ptr<MacroTaskBase>(new MacroTaskInternal(task, batch, inputrecords, outputrecords)));
        }
        taskq_ptr->add_tasks(vtask);

        if (immediate_execution) taskq_ptr->run_all(vtask);

        return result;
    }

private:

    World &world;
    std::shared_ptr<MacroTaskQ> taskq_ptr;

    /// prepare the output of the macrotask: WorldObjects must be created in the universe
    std::pair<recordlistT, resultT> prepare_output(Cloud &cloud, const argtupleT &argtuple) {
        static_assert(is_madness_function<resultT>::value || is_madness_function_vector<resultT>::value);
        resultT result = task.allocator(world, argtuple);
        recordlistT outputrecords;
        if constexpr (is_madness_function<resultT>::value) {
            outputrecords += cloud.store(world, result.get_impl().get()); // store pointer to FunctionImpl
        } else if constexpr (is_vector<resultT>::value) {
            outputrecords += cloud.store(world, get_impl(result));
        } else {
            MADNESS_EXCEPTION("\n\n  unknown result type in prepare_input ", 1);
        }
        return std::make_pair(outputrecords, result);
    }

    class MacroTaskInternal : public MacroTaskIntermediate<MacroTask_2G> {

        typedef decay_tuple<typename taskT::argtupleT> argtupleT;   // removes const, &, etc
        typedef typename taskT::resultT resultT;
        Batch batch;
        recordlistT inputrecords;
        recordlistT outputrecords;
    public:
        taskT task;

        MacroTaskInternal(const taskT &task, const Batch &batch,
                          const recordlistT &inputrecords, const recordlistT &outputrecords)
                : task(task), batch(batch), inputrecords(inputrecords), outputrecords(outputrecords) {
            static_assert(is_madness_function<resultT>::value || is_madness_function_vector<resultT>::value);
        }


        virtual void print_me(std::string s="") const {
            print("this is task",typeid(task).name(),"with batch", batch,"priority",this->get_priority());
        }

        virtual void print_me_as_table(std::string s="") const {
            std::stringstream ss;
            std::string name=typeid(task).name();
            std::size_t namesize=name.size();
            name += std::string(20-namesize,' ');

            ss  << name
                << std::setw(10) << batch
                << std::setw(5) << this->get_priority() << "        "
                <<this->stat ;
            print(ss.str());
        }

        void run(World &subworld, Cloud &cloud, MacroTaskBase::taskqT &taskq) {

            const argtupleT argtuple = cloud.load<argtupleT>(subworld, inputrecords);
            const argtupleT batched_argtuple = batch(argtuple);

            resultT result_tmp = std::apply(task, batched_argtuple);

            resultT result = get_output(subworld, cloud, argtuple);       // lives in the universe
            if constexpr (is_madness_function<resultT>::value) {
                result_tmp.compress();
                result += result_tmp;
            } else if constexpr(is_madness_function_vector<resultT>::value) {
                compress(subworld, result_tmp);
                resultT tmp1=task.allocator(subworld,argtuple);
                std::copy(result_tmp.begin(),result_tmp.end(),tmp1.begin()+batch.begin);
                result += tmp1;
            }

        };

        resultT get_output(World &subworld, Cloud &cloud, const argtupleT &argtuple) {
            resultT result;
            if constexpr (is_madness_function<resultT>::value) {
                typedef std::shared_ptr<typename resultT::implT> impl_ptrT;
                result.set_impl(cloud.load<impl_ptrT>(subworld, outputrecords));
            } else if constexpr (is_madness_function_vector<resultT>::value) {
                typedef std::shared_ptr<typename resultT::value_type::implT> impl_ptrT;
                std::vector<impl_ptrT> rimpl = cloud.load<std::vector<impl_ptrT>>(subworld, outputrecords);
                result.resize(rimpl.size());
                set_impl(result, rimpl);
            } else {
                MADNESS_EXCEPTION("unknown result type in get_output", 1);
            }
            return result;
        }

    };

};

class MicroTask {
public:
    // you need to define the exact argument(s) of operator() as tuple
    typedef std::tuple<const real_function_3d &, const double &,
            const std::vector<real_function_3d> &> argtupleT;

    using resultT = std::vector<real_function_3d>;

    // you need to define an empty constructor for the result
    // resultT must implement operator+=(const resultT&)
    resultT allocator(World &world, const argtupleT &argtuple) const {
        std::size_t n = std::get<2>(argtuple).size();
        resultT result = zero_functions_compressed<double, 3>(world, n);
        return result;
    }

    std::pair<long, hashT> partition(const argtupleT &argtuple) {
        return std::make_pair<long, hashT>(0, 0);
    }

    resultT operator()(const real_function_3d &f1, const double &arg2,
                       const std::vector<real_function_3d> &f2) const {
        World& world=f1.world();
        real_convolution_3d op=CoulombOperator(world,1.e-4,1.e-5);
        return arg2 * f1 * apply(world,op,f2);
    }

    real_function_3d density;
};


class MicroTask1 {
public:
    // you need to define the result type
    // resultT must implement operator+=(const resultT&)
    typedef real_function_3d resultT;

    // you need to define the exact argument(s) of operator() as tuple
    typedef std::tuple<const real_function_3d &, const double &,
            const std::vector<real_function_3d> &> argtupleT;

    resultT allocator(World &world, const argtupleT &argtuple) const {
        resultT result = real_factory_3d(world).compressed();
        return result;
    }

    std::pair<long, hashT> partition(const argtupleT &argtuple) {
        return std::make_pair<long, hashT>(0, 0);
    }

    resultT operator()(const real_function_3d &f1, const double &arg2,
                       const std::vector<real_function_3d> &f2) const {
        World &world = f1.world();
        Derivative<double, 3> D(world, 1);
        auto result = arg2 * f1 * inner(apply(world,D,f2), f2);
        return result;
    }
};

int test_orbital_partitioner(World &world) {
    std::tuple<double, int, std::vector<real_function_4d>, std::vector<real_function_3d>> tuple;
    auto v=OrbitalPartitioner::get_first_vector_argument(tuple);
    constexpr std::size_t index=get_index_of_first_vector_argument<0,decltype(tuple)>();
    print("test_orbital_partitioner, index ",index);

    return 1;
}


int main(int argc, char **argv) {
    madness::World &universe = madness::initialize(argc, argv);
    startup(universe, argc, argv);
    FunctionDefaults<3>::set_thresh(1.e-5);
    FunctionDefaults<3>::set_k(9);
    FunctionDefaults<3>::set_cubic_cell(-20,20);


    universe.gop.fence();
    int nworld = universe.size();
    if (universe.rank() == 0) print("creating nworld", nworld, universe.id());

    // TODO: compute batches
    // TODO: serialize member variables of tasks
    // TODO: pretty-print cloud content/input/output records
    {

        auto check_vector = [&](const std::vector<real_function_3d> &ref, const std::vector<real_function_3d> &test,
                                const std::string msg) {
            double norm_ref = norm2(universe, ref);
            double norm_test = norm2(universe, test);
            double error = norm2(universe, ref - test);
            bool success = error < 1.e-10;
            if (universe.rank() == 0) {
//                print("norm ref, test, diff", norm_ref, norm_test, error);
                if (success) print("test", msg, " \033[32m", "passed ", "\033[0m");
                else print("test", msg, " \033[31m", "failed \033[0m ");
            }
            return success;
        };
        auto check = [&](const real_function_3d &ref, const real_function_3d &test, const std::string msg) {
            double norm_ref = ref.norm2();
            double norm_test = test.norm2();
            double error = (ref - test).norm2();
            bool success = error < 1.e-10;
            if (universe.rank() == 0) {
//                print("norm ref, test, diff", norm_ref, norm_test, error);
                if (success) print("test", msg, " \033[32m", "passed ", "\033[0m");
                else print("test", msg, " \033[31m", "failed \033[0m ");
            }
            return success;
        };

        // execution in a taskq, result will be complete only after the taskq is finished
        real_function_3d f1 = real_factory_3d(universe).functor(slater(1.0));
        real_function_3d i2 = real_factory_3d(universe).functor(slater(2.0));
        real_function_3d i3 = real_factory_3d(universe).functor(slater(2.0));
        std::vector<real_function_3d> v2 = {2.0 * f1, i2};
        std::vector<real_function_3d> v3;
        for (int i=0; i<20; ++i) v3.push_back(real_factory_3d(universe).functor(slater(sqrt(double(i)))));


        timer timer1(universe);
        MicroTask t;
        std::vector<real_function_3d> ref_t = t(f1, 2.0, v3);
        timer1.tag("direct exection");

        MacroTask_2G task_immediate(universe, t);
        std::vector<real_function_3d> f2b = task_immediate(f1, 2.0, v3);
        timer1.tag("immediate taskq exection");

        auto taskq = std::shared_ptr<MacroTaskQ>(new MacroTaskQ(universe, nworld));
        taskq->set_printlevel(3);
        MacroTask_2G task(universe, t, taskq);
        std::vector<real_function_3d> f2a = task(f1, 2.0, v3);
        taskq->print_taskq();
        taskq->run_all();
        taskq->cloud.print_timings(universe);
        taskq->cloud.clear_timings();
        timer1.tag("deferred taskq execution");

        if (universe.rank()==0) print("\nstarting Microtask twice (check caching)\n");
        std::vector<real_function_3d> f2a1 = task(f1, 2.0, v3);
        std::vector<real_function_3d> f2a2 = task(f1, 2.0, v3);
        taskq->print_taskq();
        taskq->run_all();
        taskq->cloud.print_timings(universe);
        timer1.tag("executing a task twice");

        if (universe.rank()==0) print("\nstarting Microtask1\n");


        timer t2(universe);
        MicroTask1 t1;
        real_function_3d ref_t1 = t1(f1, 2.0, v3);
        t2.tag("immediate execution");

        auto taskq2 = std::shared_ptr<MacroTaskQ>(new MacroTaskQ(universe, nworld));
        taskq2->set_printlevel(3);
        MacroTask_2G task1(universe, t1, taskq2);
        real_function_3d r1 = task1(f1, 2.0, v3);
        taskq2->run_all();
        t2.tag("deferred execution");

        bool success = true;
        success = success && check(ref_t1, r1, "task1");
        success = success && check_vector(ref_t, f2a, "task again");
        success = success && check_vector(ref_t, f2b, "task immediate");

        if (universe.rank() == 0) {
            if (success) print("\n --> all tests \033[32m", "passed ", "\033[0m\n");
            else print("\n --> all tests \033[31m", "failed \033[0m \n");
        }
    }
    test_orbital_partitioner(universe);

    madness::finalize();
    return 0;
}

template<> volatile std::list<detail::PendingMsg> WorldObject<MacroTaskQ>::pending = std::list<detail::PendingMsg>();
template<> Spinlock WorldObject<MacroTaskQ>::pending_mutex(0);

//template <> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, MacroTask, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
//template <> Spinlock WorldObject<WorldContainerImpl<long, MacroTask, madness::Hash<long> > >::pending_mutex(0);

template<> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
template<> Spinlock WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending_mutex(
        0);
