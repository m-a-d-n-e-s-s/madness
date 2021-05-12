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




#include <madness/mra/mra.h>
#include <madness/world/cloud.h>
#include <madness/mra/macrotaskq.h>
#include <madness/world/world.h>


using namespace madness;
using namespace archive;

struct gaussian {
	double a;
	gaussian() : a() {};
	gaussian(double aa) : a(aa) {}
	double operator()(const coord_4d& r) const {
		double x=r[0], y=r[1], z=r[2], aa=r[3];
		return exp(-a*(x*x + y*y + z*z * aa*aa));//*abs(sin(abs(2.0*x))) *cos(y);
	}
    double operator()(const coord_3d& r) const {
        double x=r[0], y=r[1], z=r[2];
        return exp(-a*(x*x + y*y + z*z ));//*abs(sin(abs(2.0*x))) *cos(y);
    }
};


static std::atomic<int> atomic_idx;

/// a macro task for top-level operations requiring a world
//class MacroTask : public MacroTaskIntermediate<MacroTask> {
//public:
//	int inputrecord=0;
//	int outputrecord=1000;
//	int idx=0;
//	double exponent=1.0;
//
//	MacroTask(double e) : idx(atomic_idx++), exponent(e) {}
//
//    template <typename Archive>
//    void serialize(const Archive& ar) {
//    	ar & inputrecord & outputrecord & idx;
//    }
//
//    /// implement this
//	void run(World& world, Cloud& cloud, taskqT& taskq) {
//		{
//			// read the input data
//			if (world.rank()==0) print("task",idx,"is reading from record",inputrecord);
//			Function<double,4> f=real_factory_4d(world);
//			cloud.load(world,f,inputrecord);
//
//			// do the work
//			Function<double,4> g=real_factory_4d(world).functor(gaussian(exponent));
//			Function<double,4> f2=square(f)+g;
//			Derivative<double,4> D(world,1);
//			Function<double,4> df2=(D(f2)).truncate();
//			Function<double,4> df3=(D(df2)).truncate();
//			double trace=df2.trace();
//			world.gop.fence();
//
//			// store the result data
//			if (world.rank()==0) print("task",idx,"is writing to record", outputrecord);
//			cloud.store(world,f2,outputrecord);
//		}
//		world.gop.fence();
//	}
//
//
//    void print_me(std::string s="") const {
//    	print("task",s, idx,this,this->stat);
//    }
//
//};

template <typename ... Ts>
constexpr auto decay_types (std::tuple<Ts...> const &)
-> std::tuple<std::remove_cv_t<std::remove_reference_t<Ts>>...>;

template <typename T>
using decay_tuple = decltype(decay_types(std::declval<T>()));

template <typename Q> struct is_vector: std::false_type {};
template <typename Q> struct is_vector<std::vector<Q>>: std::true_type {};

template <typename>
struct is_madness_function_vector: std::false_type {};

template <typename T, std::size_t NDIM>
struct is_madness_function_vector<std::vector<typename madness::Function<T, NDIM>>>: std::true_type {};


template<typename taskT>
class MacroTask_2G {
    typedef std::pair<long,hashT> batchT;
    typedef std::list<batchT> partitionT;
    typedef typename taskT::resultT resultT;
    typedef typename taskT::argtupleT argtupleT;
    typedef Cloud::recordlistT recordlistT;
    taskT task;

public:

    MacroTask_2G(World& world, taskT& task, std::shared_ptr<MacroTaskQ> taskq_ptr=0)
            : world(world), task(task),taskq_ptr(taskq_ptr) {}

    template<typename ... Ts>
    resultT operator()(const Ts& ... args) {

        const bool immediate_execution = (not taskq_ptr);
        if (not taskq_ptr) taskq_ptr.reset(new MacroTaskQ(world,world.size()));

        auto argtuple=std::tie(args...);
        static_assert(std::is_same<decltype(argtuple),argtupleT>::value,"type or number of arguments incorrect");

        // TODO: do the partitioning
        partitionT partition;//=OrbitalPartitioner::partition_tasks(arg1);
        partition.push_back(batchT(0,1));

        recordlistT inputrecords=taskq_ptr->cloud.store(world,argtuple);
        auto [outputrecords,result] =prepare_output(taskq_ptr->cloud,argtuple);

        MacroTaskBase::taskqT vtask;
        for (const batchT& batch : partition) {
            vtask.push_back(std::shared_ptr<MacroTaskBase>(new MacroTaskInternal(task,batch,inputrecords,outputrecords)));
        }
        taskq_ptr->add_tasks(vtask);

        if (immediate_execution) taskq_ptr->run_all(vtask);

        return result;
    }

private:

    World& world;
    std::shared_ptr<MacroTaskQ> taskq_ptr;

    /// prepare the output of the macrotask: WorldObjects must be created in the universe
    std::pair<recordlistT, resultT> prepare_output(Cloud& cloud, const argtupleT& argtuple) {
        static_assert(is_madness_function<resultT>::value || is_madness_function_vector<resultT>::value);
        resultT result=task.allocator(world,argtuple);
        recordlistT outputrecords;
        if constexpr (is_madness_function<resultT>::value) {
            outputrecords+=cloud.store(world, result.get_impl().get()); // store pointer to FunctionImpl
        } else if constexpr (is_vector<resultT>::value) {
            outputrecords+=cloud.store(world,result.size());
            for (auto& r : result) outputrecords+= cloud.store(world, r.get_impl().get());
        } else {
            MADNESS_EXCEPTION("\n\n  unknown result type in prepare_input ",1);
        }
        return std::make_pair(outputrecords, result);
    }

    class MacroTaskInternal : public MacroTaskIntermediate<MacroTask_2G> {

        typedef decay_tuple<typename taskT::argtupleT> argtupleT;   // removes const, &, etc
        typedef typename taskT::resultT resultT;
        batchT batch;
        recordlistT inputrecords;
        recordlistT outputrecords;
    public:
        taskT task;
        MacroTaskInternal(const taskT& task, const batchT& batch,
                          const recordlistT& inputrecords, const recordlistT& outputrecords)
                          : task(task), batch(batch), inputrecords(inputrecords), outputrecords(outputrecords) {}

        void run(World& subworld, Cloud& cloud, MacroTaskBase::taskqT& taskq) {

            argtupleT argtuple=cloud.load<argtupleT>(subworld,inputrecords);
            resultT result_tmp;

            constexpr std::size_t narg=(std::tuple_size<argtupleT>::value);
            if constexpr (narg==1) result_tmp=task(std::get<0>(argtuple));
            if constexpr (narg==2) result_tmp=task(std::get<0>(argtuple),std::get<1>(argtuple));
            if constexpr (narg==3) result_tmp=task(std::get<0>(argtuple),std::get<1>(argtuple),
                                               std::get<2>(argtuple));
            if constexpr (narg==4) result_tmp=task(std::get<0>(argtuple),std::get<1>(argtuple),
                                               std::get<2>(argtuple),std::get<3>(argtuple));
            if constexpr (narg==5) result_tmp=task(std::get<0>(argtuple),std::get<1>(argtuple),
                                               std::get<2>(argtuple),std::get<3>(argtuple),
                                               std::get<4>(argtuple));
            if constexpr (narg==6) result_tmp=task(std::get<0>(argtuple),std::get<1>(argtuple),
                                               std::get<2>(argtuple),std::get<3>(argtuple),
                                               std::get<4>(argtuple),std::get<5>(argtuple));
            static_assert(is_madness_function<resultT>::value || is_madness_function_vector<resultT>::value);
            if constexpr (is_madness_function<resultT>::value) result_tmp.compress();
            else if constexpr(is_madness_function_vector<resultT>::value) compress(subworld,result_tmp);

            resultT result=get_output(subworld, cloud, argtuple);       // lives in the universe
            result+=result_tmp;
        };

        resultT get_output(World& subworld, Cloud& cloud, const argtupleT& argtuple) {
            resultT result;
            if constexpr (is_madness_function<resultT>::value) {
                typedef std::shared_ptr<typename resultT::implT> impl_ptrT;
                result.set_impl(cloud.load<impl_ptrT>(subworld,outputrecords));
            } else if constexpr (is_madness_function_vector<resultT>::value){
                typedef std::shared_ptr<typename resultT::value_type::implT> impl_ptrT;
                std::size_t n=cloud.load<std::size_t>(subworld,outputrecords);
                result.resize(n);
                std::vector<impl_ptrT> rimpl(n);
                for (std::size_t i=0; i<n; ++i) {
                    result[i].set_impl(cloud.load<impl_ptrT>(subworld,outputrecords));
                }
            } else {
                MADNESS_EXCEPTION("unknown result type in get_output",1);
            }
            return result;
        }

    };

};

class MicroTask {
public:
    // you need to define the exact argument(s) of operator() as tuple
    typedef std::tuple<const real_function_3d&, const double&,
            const std::vector<real_function_3d>&> argtupleT;

    using resultT = std::vector<real_function_3d>;

    // you need to define an empty constructor for the result
    // resultT must implement operator+=(const resultT&)
    resultT allocator(World& world, const argtupleT& argtuple) const {
        std::size_t n=std::get<2>(argtuple).size();
        resultT result=zero_functions_compressed<double,3>(world,n);
        return result;
    }

    std::pair<long,hashT> partition(const argtupleT& argtuple) {
        return std::make_pair<long,hashT>(0,0);
    }

    resultT operator()(const real_function_3d& f1, const double& arg2,
            const std::vector<real_function_3d>& f2) const {
        return arg2*f1*f2;
    }

    real_function_3d density;
};


class MicroTask1 {
public:
    // you need to define the result type
    // resultT must implement operator+=(const resultT&)
    typedef real_function_3d resultT;

    // you need to define the exact argument(s) of operator() as tuple
    typedef std::tuple<const real_function_3d&, const double&,
            const std::vector<real_function_3d>&> argtupleT;

    resultT allocator(World& world, const argtupleT& argtuple) const {
       resultT result=real_factory_3d(world).compressed();
       return result;
    }
    std::pair<long,hashT> partition(const argtupleT& argtuple) {
        return std::make_pair<long,hashT>(0,0);
    }

    resultT operator()(const real_function_3d& f1, const double& arg2,
                       const std::vector<real_function_3d>& f2) const {
        print("Hello world from MicroTask 1");
        World& world=f1.world();
        Derivative<double,3> D(world,1);
        return arg2*f1*inner(D(f1),f1);
    }
};


/// similar to MacroTask, for testing heterogeneous task queues
//class MacroTask1 : public MacroTaskIntermediate<MacroTask1> {
//
//public:
//	int idx=0;
//	int inputrecord=0;
//	int outputrecord=1000;
//
//	MacroTask1(): idx(atomic_idx++) {}
//
//	void run(World& world, Cloud& cloud, taskqT& taskq) {
//		Function<double,4> f=cloud.load<Function<double,4> > (world,inputrecord);
//		Function<double,4> g=real_factory_4d(world).functor(gaussian(2.0));
//		Function<double,4> f2=square(f)+g;
//		world.gop.fence();
//		cloud.store(world,g,outputrecord);
//	}
//
//    template <typename Archive>
//    void serialize(const Archive& ar) {
//    	ar & inputrecord & outputrecord & idx;
//    }
//
//    void print_me(std::string s="") const {
//    	print("macro task 1",s, idx,this,this->stat);
//    }
//
//};



int main(int argc, char** argv) {
    madness::World& universe = madness::initialize(argc, argv);
    startup(universe,argc,argv);
    FunctionDefaults<4>::set_thresh(1.e-9);
    FunctionDefaults<4>::set_k(7);

    universe.gop.fence();
    int nworld=universe.size();
    if (universe.rank()==0) print("creating nworld",nworld);

    // TODO: compute batches
    // TODO: cache input functions
    // TODO: serialize member variables of tasks
    // TODO: pretty-print cloud content/input/output records
    {

        auto check_vector = [&](const std::vector<real_function_3d>& ref, const std::vector<real_function_3d>& test, const std::string msg)
        {
            double norm_ref=norm2(universe,ref);
            double norm_test=norm2(universe,test);
            double error=norm2(universe,ref-test);
            bool success=error<1.e-10;
            if (universe.rank()==0) {
//                print("norm ref, test, diff", norm_ref, norm_test, error);
                if (success) print("test",msg," \033[32m", "passed ", "\033[0m");
                else print("test",msg," \033[31m", "failed \033[0m ");
            }
            return success;
        };
        auto check = [&](const real_function_3d& ref, const real_function_3d& test, const std::string msg)
                {
            double norm_ref=ref.norm2();
            double norm_test=test.norm2();
            double error=(ref-test).norm2();
            bool success=error<1.e-10;
            if (universe.rank()==0) {
//                print("norm ref, test, diff", norm_ref, norm_test, error);
                if (success) print("test",msg," \033[32m", "passed ", "\033[0m");
                else print("test",msg," \033[31m", "failed \033[0m ");
            }
            return success;
        };

        // execution in a taskq, result will be complete only after the taskq is finished
        real_function_3d f1=real_factory_3d(universe).functor(gaussian(1.0));
        real_function_3d i2=real_factory_3d(universe).functor(gaussian(2.0));
        real_function_3d i3=real_factory_3d(universe).functor(gaussian(2.0));
        std::vector<real_function_3d> v2={2.0*f1,i2};
        std::vector<real_function_3d> v3={2.0*f1,i3,i2};

        MicroTask t;
        MicroTask1 t1;
        std::vector<real_function_3d> ref_t=t(f1,2.0,v2);
        real_function_3d ref_t1=t1(f1,2.0,v3);

        auto taskq=std::shared_ptr<MacroTaskQ> (new MacroTaskQ(universe,nworld));
//        taskq->set_printlevel(10);
        MacroTask_2G task(universe,t,taskq);
        MacroTask_2G task_immediate(universe,t);
        MacroTask_2G task1(universe,t1,taskq);

        std::vector<real_function_3d> f2=task(f1,2.0, v2);
        std::vector<real_function_3d> f2a=task(f1,2.0, v2);
        std::vector<real_function_3d> f2b=task_immediate(f1,2.0, v2);

        real_function_3d f2_1=task1(f1,2.0, v3);

        taskq->run_all();

        bool success=true;
        success = success && check_vector(ref_t,f2,"task1");
        success = success && check_vector(ref_t,f2a,"task1 again");
        success = success && check_vector(ref_t,f2b,"task1 immediate");
        success = success && check(ref_t1,f2_1,"task2");

        if (universe.rank()==0) {
            if (success) print("\n --> all tests \033[32m", "passed ", "\033[0m\n");
            else print("\n --> all tests \033[31m", "failed \033[0m \n");
        }
    }

    madness::finalize();
    return 0;
}

template <> volatile std::list<detail::PendingMsg> WorldObject<MacroTaskQ>::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<MacroTaskQ>::pending_mutex(0);

//template <> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, MacroTask, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
//template <> Spinlock WorldObject<WorldContainerImpl<long, MacroTask, madness::Hash<long> > >::pending_mutex(0);

template <> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending_mutex(0);
