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

#include <madness/mra/macrotaskq.h>
#include <madness/world/cloud.h>
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


template <typename>
struct is_madness_function_vector: std::false_type {};

template <typename T, std::size_t NDIM>
struct is_madness_function_vector<std::vector<typename std::is_same<T,typename madness::Function<T, NDIM>>>>: std::true_type {};

template<typename taskT>
class MacroTask_2G {
    typedef std::pair<long,hashT> batchT;
    typedef std::list<batchT> partitionT;
    typedef typename taskT::resultT resultT;
    typedef typename taskT::argtupleT argtupleT;
    typedef Cloud::recordlistT recordlistT;
    taskT task;

public:
    MacroTask_2G(World& world, taskT& task) : world(world), task(task), taskq_ptr() {}
    MacroTask_2G(World& world, taskT& task, std::shared_ptr<MacroTaskQ> taskq_ptr) : world(world), task(task), taskq_ptr(taskq_ptr) {}

    template<typename ... Ts>
    resultT operator()(const Ts& ... args) {

        // a MacroTask cannot be called twice, because the output will be overwritten
        // feel free to add this functionality if you need
        if (outputrecords.size()!=0) {
            print("\n\nU can call a MacroTask only once because otherwise the output ",
                  "will be overwritten.\nFeel free to add this functionality\n");
            MADNESS_EXCEPTION("missing feature",1);
        }

        auto allargs=std::tie(args...);
        static_assert(std::is_same<decltype(allargs),argtupleT>::value,"type or number of arguments incorrect");

        // TODO: do the partitioning
        partitionT partition;//=OrbitalPartitioner::partition_tasks(arg1);
        partition.push_back(batchT(0,1));

        inputrecords=taskq_ptr->cloud.store(world,allargs);
        resultT result=prepare_output(taskq_ptr->cloud);

        MacroTaskBase::taskqT vtask;
        for (const batchT& batch : partition) {
            vtask.push_back(std::shared_ptr<MacroTaskBase>(new MacroTaskInternal(task,batch,inputrecords,outputrecords)));
        }

        if (not deferred_execution())  {
            MacroTaskQ taskq(world, world.size());
            taskq.add_tasks(vtask);
            taskq.run_all(vtask);
        } else {
            taskq_ptr->add_tasks(vtask);
        }

        return result;
    }

private:

    World& world;
    std::shared_ptr<MacroTaskQ> taskq_ptr;
    recordlistT inputrecords;
    recordlistT outputrecords;


    bool deferred_execution() const {return (taskq_ptr) ? true : false;}

    // TODO: return bookkeeping information, serialize
    /// prepare the output of the macrotask: WorldObjects must be created in the universe
    resultT prepare_output(Cloud& cloud) {
        // TODO: generalize this
        static_assert(is_madness_function<resultT>::value);
        if constexpr (is_madness_function<resultT>::value) {
            resultT result(world);
            result.compress();
            outputrecords+=cloud.store(world, result.get_impl().get()); // store pointer to FunctionImpl
            return result;
        } else {
            resultT result;
            cloud.store(world, result, 0); // store result
            return result;
        }
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

//            argtupleT  argtuple(subworld);
//            argtupleT argtuple=get_input<argtupleT>(subworld,cloud,inputrecords);
//            cloud.set_debug(true);
            argtupleT argtuple=cloud.load<argtupleT>(subworld,inputrecords);
            resultT result=get_output(subworld, cloud);
            print("result.id in run()",result.get_impl()->id());


            constexpr std::size_t narg=(std::tuple_size<argtupleT>::value);
            if constexpr (narg==1) result+=task(std::get<0>(argtuple));
            if constexpr (narg==2) result+=task(std::get<0>(argtuple),std::get<1>(argtuple));
            if constexpr (narg==3) result+=task(std::get<0>(argtuple),std::get<1>(argtuple),
                                                std::get<2>(argtuple));
            if constexpr (narg==4) result+=task(std::get<0>(argtuple),std::get<1>(argtuple),
                                                std::get<2>(argtuple),std::get<3>(argtuple));
            if constexpr (narg==5) result+=task(std::get<0>(argtuple),std::get<1>(argtuple),
                                                std::get<2>(argtuple),std::get<3>(argtuple),
                                                std::get<4>(argtuple));
            if constexpr (narg==6) result+=task(std::get<0>(argtuple),std::get<1>(argtuple),
                                                std::get<2>(argtuple),std::get<3>(argtuple),
                                                std::get<4>(argtuple),std::get<5>(argtuple));
            MADNESS_CHECK(is_madness_function<resultT>::value);
        };

        resultT get_output(World& subworld, Cloud& cloud) {
            resultT result;
            if constexpr (is_madness_function<resultT>::value) {
                typedef std::shared_ptr<typename resultT::implT> impl_ptrT;
                impl_ptrT rimpl;
                rimpl=cloud.load<impl_ptrT>(subworld,outputrecords);
                result.set_impl(rimpl);
            } else {
                MADNESS_CHECK(is_madness_function<resultT>::value);
            }

            return result;
        }

    };

};

class MicroTask {
public:
    // you need to define the result type
    // resultT must implement operator+=(const resultT&)
    typedef real_function_3d resultT;

    // you need to define the exact argument(s) of operator() as tuple
    typedef std::tuple<const real_function_3d&, const double&,
            const std::vector<real_function_3d>&> argtupleT;

    std::pair<long,hashT> partition(const argtupleT& argtuple) {
        return std::make_pair<long,hashT>(0,0);
    }

    resultT operator()(const real_function_3d& f1, const double& arg2,
            const std::vector<real_function_3d>& f2) const {
        double norm=f1.norm2();
        print("Hello world from MicroTask",norm);
        World& world=f1.world();
        resultT result= arg2*f1*dot(world,f2,f2);
        double normresult=result.norm2();
        print("normresult in task",normresult);
        return result;
    }

    real_function_3d density;

    void store_member_variables_to_cloud() const {

    }

};


class MicroTask1 {
public:
    // you need to define the result type
    // resultT must implement operator+=(const resultT&)
    typedef real_function_3d resultT;

    // you need to define the exact argument(s) of operator() as tuple
    typedef std::tuple<const real_function_3d&, const double&,
            const std::vector<real_function_3d>&> argtupleT;

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
//    madness::World& universe = madness::initialize(argc,argv);
    madness::World& universe = madness::initialize(argc, argv);
    startup(universe,argc,argv);
    FunctionDefaults<4>::set_thresh(1.e-9);
    FunctionDefaults<4>::set_k(7);

    std::cout << "Hello from " << universe.rank() << std::endl;
    universe.gop.fence();
//    int nworld=std::min(int(universe.size()),int(3));
    int nworld=universe.size();
//    int nworld=2;
    if (universe.rank()==0) print("creating nworld",nworld);

//    if (0) {
//    	MacroTaskQ taskq(universe,nworld);
//
//		long ntask=15;
//
//		Function<double,4> f=real_factory_4d(universe);
//		f=f+1.0;
//		taskq.cloud.store(universe,f,0);
//		universe.gop.fence();
//
//		MacroTaskBase::taskqT vtask;
//		for (int i=0; i<ntask; ++i) {
//			MacroTask task(double(ntask-i)/2.0);
//			task.inputrecord=0;
//			task.outputrecord=i+1;
//			vtask.push_back(std::shared_ptr<MacroTaskBase>(new MacroTask(task)));
//		}
//		universe.gop.fence();
//
//		// try heterogeneous task list
//	    vtask.push_back(std::shared_ptr<MacroTaskBase>(new MacroTask1()));
//
//	    FunctionDefaults<4>::set_default_pmap(taskq.get_subworld());
//	    taskq.run_all(vtask);
//		FunctionDefaults<4>::set_default_pmap(universe);
//
//		taskq.cloud.print_timings(universe);
//	//    MacroTask<real_function_4d, dataT> task;
//	//    task.idx=3;
//	//    std::vector<Function<double,4> > result=taskq.map(task,vdata);
//
//	//    print_size(universe,result,"result after map");
//    }
	universe.gop.fence();

    {
        // execution in a taskq, result will be complete only after the taskq is finished
        real_function_3d f1=real_factory_3d(universe).functor(gaussian(1.0));
        real_function_3d i2=real_factory_3d(universe).functor(gaussian(2.0));
        std::vector<real_function_3d> v2={2.0*f1,i2};

        MicroTask t;
        MicroTask1 t1;
        real_function_3d f3=t(f1,2.0,v2);
        real_function_3d f3_1=t1(f1,2.0,v2);

        auto taskq=std::shared_ptr<MacroTaskQ> (new MacroTaskQ(universe,nworld));
        taskq->set_printlevel(10);
        MacroTask_2G task(universe,t,taskq);
        MacroTask_2G task1(universe,t1,taskq);

        real_function_3d f2=task(f1,2.0, v2);
        real_function_3d f2_1=task1(f1,2.0, v2);
        double norm2a=(f2).norm2();

        print("before running all tasks");
	    FunctionDefaults<3>::set_default_pmap(taskq->get_subworld());
	    taskq->set_printlevel(10);
        universe.gop.fence();

        taskq->run_all();
		FunctionDefaults<3>::set_default_pmap(universe);
        universe.gop.fence();
        double norm2=(f2).norm2();
        double norm3=(f3).norm2();
        double error=(f3-f2).norm2();
        double norm_f2_1=(f2_1).norm2();
        double norm_f3_1=(f3_1).norm2();
        double error1=(f3_1-f2_1).norm2();
        print("norm2/2a/3, error",norm2, norm2a, norm3, error);
        print("norm2_1/3_1, error", norm_f2_1,norm_f3_1, error1);
        if (error<1.e-10) print("test_vectormacrotask \033[32m"  ,"passed ", "\033[0m");
        else print("test_vectormacrotask \033[31m", "failed \033[0m ");
        if (error1<1.e-10) print("test_vectormacrotask \033[32m"  ,"passed ", "\033[0m");
        else print("test_vectormacrotask \033[31m", "failed \033[0m ");
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
