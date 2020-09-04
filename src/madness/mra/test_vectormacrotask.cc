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

    {
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

    madness::finalize();
    return 0;
}

template <> volatile std::list<detail::PendingMsg> WorldObject<MacroTaskQ>::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<MacroTaskQ>::pending_mutex(0);

template <> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, MacroTask, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<WorldContainerImpl<long, MacroTask, madness::Hash<long> > >::pending_mutex(0);

template <> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending_mutex(0);
