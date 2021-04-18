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
class MacroTask : public MacroTaskIntermediate<MacroTask> {
public:
	int inputrecord=0;
	int outputrecord=1000;
	int idx=0;
	double exponent=1.0;

	MacroTask(double e) : idx(atomic_idx++), exponent(e) {}

    template <typename Archive>
    void serialize(const Archive& ar) {
    	ar & inputrecord & outputrecord & idx;
    }

    /// implement this
	void run(World& world, Cloud& cloud, taskqT& taskq) {
		{
			// read the input data
			if (world.rank()==0) print("task",idx,"is reading from record",inputrecord);
			Function<double,4> f=real_factory_4d(world);
			cloud.load(world,f,inputrecord);

			// do the work
			Function<double,4> g=real_factory_4d(world).functor(gaussian(exponent));
			Function<double,4> f2=square(f)+g;
			Derivative<double,4> D(world,1);
			Function<double,4> df2=(D(f2)).truncate();
			Function<double,4> df3=(D(df2)).truncate();
			double trace=df2.trace();
			world.gop.fence();

			// store the result data
			if (world.rank()==0) print("task",idx,"is writing to record", outputrecord);
			cloud.store(world,f2,outputrecord);
		}
		world.gop.fence();
	}


    void print_me(std::string s="") const {
    	print("task",s, idx,this,this->stat);
    }

};


template<typename taskT>
class MacroTask_2G {

    typedef std::pair<long,hashT> batchT;
    typedef std::list<batchT> partitionT;
    typedef typename taskT::resultT resultT;
    taskT task;
public:
    class MacroTaskInternal : public MacroTaskIntermediate<MacroTask> {

        typedef typename taskT::argT1 argT1;
        typedef typename taskT::argT2 argT2;
        typedef typename taskT::resultT resultT;
        batchT batch;
    public:
        taskT task;
        MacroTaskInternal(const taskT& task, const batchT& batch) : task(task), batch(batch) {}

        void run(World& subworld, Cloud& cloud, MacroTaskBase::taskqT& taskq) {

            argT1 arg1=get_input<argT1>(subworld,cloud);
            argT2 arg2=argT2(2);//=get_input<argT2>(subworld,cloud);
            resultT result=get_output(subworld, cloud);
            result+=task(arg1,arg2);
        };

        resultT get_output(World& subworld, Cloud& cloud) {
            resultT result;
            std::shared_ptr<FunctionImpl<double,3> > rimpl;
            cloud.load(subworld,rimpl,0);
            result.set_impl(rimpl);
            return result;
        }

        template<typename argT>
        argT get_input(World& subworld, Cloud& cloud) {
            hashT inputrecord=batch.second;
            argT arg;
            cloud.load(subworld,arg,inputrecord);
            return arg;
        }
    };

public:
    MacroTask_2G(World& world, taskT& task) : world(world), task(task), taskq_ptr() {}
    MacroTask_2G(World& world, taskT& task, std::shared_ptr<MacroTaskQ> taskq_ptr) : world(world), task(task), taskq_ptr(taskq_ptr) {}

    template<typename argT1, typename argT2>
    resultT operator()(const argT1& arg1, const argT2& arg2) {
        partitionT partition;//=OrbitalPartitioner::partition_tasks(arg1);
        partition.push_back(batchT(0,1));

        prepare_input(partition,arg1,arg2);
        resultT result=prepare_output(taskq_ptr->cloud);

        MacroTaskBase::taskqT vtask;
        for (const batchT& batch : partition) {
            vtask.push_back(std::shared_ptr<MacroTaskBase>(new MacroTaskInternal(task,batch)));
        }

        if (not deferred_execution())  {
            MacroTaskQ taskq(world, world.size());
            taskq.add_tasks(vtask);
            taskq.run_all(vtask);
            world.gop.fence();
        } else {
            taskq_ptr->add_tasks(vtask);
        }
        return result;
    }


private:
    bool deferred_execution() const {return (taskq_ptr) ? true : false;}

    template<typename argT1, typename argT2>
    void prepare_input(const partitionT& partition, const argT1& arg1, const argT2& arg2) {
        hashT inputrecord=partition.front().second;
//        auto storage=std::tie(arg1,arg2);
//        taskq_ptr->cloud.store(world,storage,inputrecord);
        if (std::is_constructible<argT1,World&>::value) taskq_ptr->cloud.store(world,arg1,inputrecord);
        if (std::is_constructible<argT2,World&>::value) taskq_ptr->cloud.store(world,arg2,inputrecord);
    }

    resultT prepare_output(Cloud& cloud) {
        resultT result(world);
        cloud.store(world, result.get_impl(), 0); // store pointer to FunctionImpl
        return result;
    }


    World& world;
    std::shared_ptr<MacroTaskQ> taskq_ptr;

};

class MicroTask {
public:
    typedef real_function_3d resultT;
    typedef real_function_3d argT1;
    typedef double argT2;
    resultT operator()(const real_function_3d& arg1, const double arg2) const {
        print("Hello world from MicroTask", arg2);
        myusleep(2.e6);
        return arg2*arg1;
    }
};


/// similar to MacroTask, for testing heterogeneous task queues
class MacroTask1 : public MacroTaskIntermediate<MacroTask1> {

public:
	int idx=0;
	int inputrecord=0;
	int outputrecord=1000;

	MacroTask1(): idx(atomic_idx++) {}

	void run(World& world, Cloud& cloud, taskqT& taskq) {
		Function<double,4> f=cloud.load<Function<double,4> > (world,inputrecord);
		Function<double,4> g=real_factory_4d(world).functor(gaussian(2.0));
		Function<double,4> f2=square(f)+g;
		world.gop.fence();
		cloud.store(world,g,outputrecord);
	}

    template <typename Archive>
    void serialize(const Archive& ar) {
    	ar & inputrecord & outputrecord & idx;
    }

    void print_me(std::string s="") const {
    	print("macro task 1",s, idx,this,this->stat);
    }

};




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
    if (universe.rank()==0) print("creating nworld",nworld);

    if (0) {
    	MacroTaskQ taskq(universe,nworld);

		long ntask=15;

		Function<double,4> f=real_factory_4d(universe);
		f=f+1.0;
		taskq.cloud.store(universe,f,0);
		universe.gop.fence();

		MacroTaskBase::taskqT vtask;
		for (int i=0; i<ntask; ++i) {
			MacroTask task(double(ntask-i)/2.0);
			task.inputrecord=0;
			task.outputrecord=i+1;
			vtask.push_back(std::shared_ptr<MacroTaskBase>(new MacroTask(task)));
		}
		universe.gop.fence();

		// try heterogeneous task list
	    vtask.push_back(std::shared_ptr<MacroTaskBase>(new MacroTask1()));

	    FunctionDefaults<4>::set_default_pmap(taskq.get_subworld());
	    taskq.run_all(vtask);
		FunctionDefaults<4>::set_default_pmap(universe);

		taskq.cloud.print_timings(universe);
	//    MacroTask<real_function_4d, dataT> task;
	//    task.idx=3;
	//    std::vector<Function<double,4> > result=taskq.map(task,vdata);

	//    print_size(universe,result,"result after map");
    }
	universe.gop.fence();

    {
        // execution in a taskq, result will be complete only after the taskq is finished
        real_function_3d f1=real_factory_3d(universe).functor(gaussian(1.0));

        MicroTask t;
        real_function_3d f3=t(f1,2.0);

        auto taskq=std::shared_ptr<MacroTaskQ> (new MacroTaskQ(universe,nworld));
        taskq->set_printlevel(10);
        MacroTask_2G<MicroTask> task(universe,t,taskq);
        real_function_3d f2=task(f1,2.0);
        double norm2a=(f2).norm2();

        print("before running all tasks");
        taskq->run_all();
        universe.gop.fence();
        double norm2=(f2).norm2();
        double norm3=(f3).norm2();
        double error=(f3-f2).norm2();
        print("norm2/2a/3, error",norm2, norm2a, norm3, error);
    }

    madness::finalize();
    return 0;
}

template <> volatile std::list<detail::PendingMsg> WorldObject<MacroTaskQ>::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<MacroTaskQ>::pending_mutex(0);

template <> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, MacroTask, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<WorldContainerImpl<long, MacroTask, madness::Hash<long> > >::pending_mutex(0);

template <> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending_mutex(0);
