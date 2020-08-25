#ifndef MAD_PARALLEL_DC_ARCHIVE_H_INCLUDED
#define MAD_PARALLEL_DC_ARCHIVE_H_INCLUDED

#include <madness/world/MADworld.h>
#include <madness/world/worlddc.h>
#include <madness/world/vector_archive.h>

namespace madness {



    namespace archive {
        
        class ContainerRecordOutputArchive : public BaseOutputArchive {
            using keyT = long;
            using containerT = WorldContainer<keyT,std::vector<unsigned char>>;
            World& subworld;
            keyT key;
            containerT& dc; // lifetime???
            std::vector<unsigned char> v;
            VectorOutputArchive ar;
            
        public:

            ContainerRecordOutputArchive(World& subworld, containerT& dc, const keyT& key)
                : subworld(subworld)
                , key(key)
                , dc(dc)
                , v()
                , ar(v)
            {}
            
            ~ContainerRecordOutputArchive()
            {
                close();
            }
            
            template <class T>
            inline
            typename std::enable_if< madness::is_trivially_serializable<T>::value, void >::type
            store(const T* t, long n) const {
                MADNESS_CHECK(subworld.rank() == 0);
                ar.store(t, n);
            }
            
            void open() {}
            
            void flush() {}
            
            void close() {
                if (subworld.rank() == 0) dc.replace(key,v);
            }
        };
        
        class ContainerRecordInputArchive : public BaseInputArchive {
            using keyT = long;
            using containerT = WorldContainer<keyT,std::vector<unsigned char>>;
            ProcessID rank;
            std::vector<unsigned char> v;
            VectorInputArchive ar;
            
        public:
            ContainerRecordInputArchive(World& subworld, containerT& dc, const keyT& key)
                : rank(subworld.rank())
                , v()
                , ar(v)
            {
                if (rank==0) {
                    containerT::const_iterator it = dc.find(key).get();
                    if (it != dc.end()) {
                        v = it->second;
                    }
                    else {
                    	std::cout << "key " << key << " in world " << subworld.id()
                    			<< "dc.world " << dc.get_world().id() << std::endl;
                        MADNESS_EXCEPTION("record not found", key);
                    }
                }
            }
            
            ~ContainerRecordInputArchive()
            {}
            
            template <class T>
            inline
            typename std::enable_if< madness::is_trivially_serializable<T>::value, void >::type
            load(T* t, long n) const {
                MADNESS_CHECK(rank == 0);
                ar.load(t,n);
            }
            
            void open() {}
            
            void flush() {}
            
            void close() {}
        };
        
        void xxxtest(World& world) {
            //std::vector<unsigned char> v;
            WorldContainer<long,std::vector<unsigned char>> dc(world);
            if (world.rank() == 0) {
                //VectorOutputArchive ar(v);
                ContainerRecordOutputArchive ar(world, dc, 99);
                int a=1, b=7;
                ar & a & b;
            }
            if (world.rank() == 0) {
                //VectorInputArchive ar(v);
                ContainerRecordInputArchive ar(world, dc, 99);
                int a, b;
                ar & a & b;
                std::cout << "I read " << a << " " << b << std::endl;
            }
            dc.get_world().gop.fence();
        }
    }


    class Cloud  {

    	madness::WorldContainer<long,std::vector<unsigned char> > container;
    	bool debug=false;

    public:
    	/// @param[in]	w	the universe world
    	Cloud(madness::World& w) : container(w) {
    		MADNESS_ASSERT(w.id()==0);
    	}

    	void set_debug(bool value) {
    		debug=value;
    	}

    	template<typename T>
    	void store(madness::World& world, const T& source, const int record) {
    		if (debug) std::cout << "storing " << typeid(T).name() << " to record " << record << std::endl;

    		// scope is important because of destruction ordering of world objects and fence
    		{
    			madness::archive::ContainerRecordOutputArchive ar(world,container,record);
    			madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(world, ar);
    			par & source;
    		}
    		world.gop.fence();
    		if (debug) std::cout << "done with storing; container size " << container.size() << std::endl;
    	}

    	// call this load if your argument has a constructor taking world as argument (e.g. WorldContainer)
    	template<typename T>
    	typename std::enable_if<std::is_constructible<T,World&>::value, T>::type
    	load(madness::World& world, const int record) {
            if (debug) std::cout << "loading " << typeid(T).name() << " to world " << world.id() << " from record " << record << std::endl;
    		T target(world);
            madness::archive::ContainerRecordInputArchive ar(world,container,record);
            madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(world, ar);
            par & target;
            return target;
    	}

    	// call this load for simple types
    	template<typename T>
    	typename std::enable_if<!std::is_constructible<T,World&>::value, T>::type
    	load(madness::World& world, const int record) {
            if (debug) std::cout << "loading " << typeid(T).name() << " to world " << world.id() << " from record " << record << std::endl;
    		T target;
            madness::archive::ContainerRecordInputArchive ar(world,container,record);
            madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(world, ar);
            par & target;
            return target;
    	}

    	// call this load to pass in an already constructed argument
    	template<typename T>
    	void load(madness::World& world, T& target, const int record) {
            if (debug) std::cout << "loading " << typeid(T).name() << " to world " << world.id() << " from record " << record << std::endl;
            madness::archive::ContainerRecordInputArchive ar(world,container,record);
            madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(world, ar);
            par & target;
    	}

    };


	struct test_cloud_pod {
		int i;
		double d;
		madness::WorldContainer<int,double> c;
		test_cloud_pod(madness::World& world) : i(), d(), c(world) {}
		template<typename Archive>
		void serialize(Archive& ar) {
			ar & i & d &c;
		}
	};

    void test_cloud(madness::World& universe) {
    	typedef madness::WorldContainer<int,double> dcT;

    	madness::Cloud cloud(universe);
    	cloud.set_debug(true);

    	std::cout << "dcT::is_constructible " << std::is_constructible<dcT,World&>::value << std::endl;
    	std::cout << "int::is_constructible " << std::is_constructible<int,World>::value << std::endl;

        // Make subworlds
        int color = universe.rank() % 2;
        int nsubworld=2;
        SafeMPI::Intracomm comm = universe.mpi.comm().Split(color, universe.rank() / 2);

        madness::World subworld(comm);
        universe.gop.fence();
        subworld.gop.fence();

        // store simple types
        {
        	cloud.store(universe,int(1),0);
        	int i=cloud.load<int>(subworld,0);
        	int j=0;
        	cloud.load(subworld,j,0);
        	std::cout << "j= " << j << std::endl;
        }
        {
        	cloud.store(universe,double(1),0);
        	int i=cloud.load<double>(subworld,0);
        }

        // WorldObjects
        {
        	dcT w(universe);
        	cloud.store(universe,w,0);
        	dcT w1=cloud.load<dcT>(subworld,0);
        	dcT w2(subworld);
        	cloud.load(subworld,w2,0);
        }

//        {
//        	test_cloud_pod p(universe), p2(subworld);
//    		p.c=dcT(universe);
//    		cloud.store(universe,p,0);
//
//    		test_cloud_pod p1=cloud.load<test_cloud_pod>(subworld,0);
//    		cloud.load(subworld,p2,0);
//
//    		cloud.store(universe,std::make_tuple(color,p.c),1);
//    		int i3; dcT c3(subworld);
//    		auto p3=std::tie(i3,c3);
//    		cloud.load(subworld,p3,1);
//        }
        {


    		// prepare universe functions
    		dcT a(universe);
    		for (int i=0; i<10; ++i) a.replace(i,0);
    		std::cout << "a worldid " << a.get_world().id() << std::endl;
    		for (int i=0; i<10; ++i) std::cout << "a " << a.find(i).get()->first << " " << a.find(i).get()->second  << std::endl;
    		std::cout << "a.size (local to universe rank" << universe.rank() <<") " << a.size() << std::endl;
    		cloud.store(universe,a,1);

    		// do subworld work
    		dcT sub_a=cloud.load<dcT>(subworld,1);	// sub_a is now replicated
    		for (int i=0; i<10; ++i) sub_a.replace(i,double(i)*subworld.id());
    		for (int i=0; i<10; ++i) std::cout << "sub_a " << sub_a.find(i).get()->first << " " << sub_a.find(i).get()->second  << std::endl;
    		subworld.gop.fence();
    		cloud.store(subworld,sub_a,color);		// sub_a is stored on different records


    		// collect data into universe
    		dcT b0=cloud.load<dcT>(universe,0);		// b0 is from subworld 1, b1 is from subworld 1000
    		std::cout << "b0.size (local to universe rank" << universe.rank() <<") " << b0.size() << std::endl;
    		dcT b1=cloud.load<dcT>(universe,1);
    		std::cout << "b1.size (local to universe rank" << universe.rank() <<") " << b1.size() << std::endl;

    		universe.gop.fence();
    		for (int i=0; i<10; ++i) std::cout << "b0 " << b0.find(i).get()->first << " " << b0.find(i).get()->second  << std::endl;
    		universe.gop.fence();
    		for (int i=0; i<10; ++i) std::cout << "b1 " << b1.find(i).get()->first << " " << b1.find(i).get()->second  << std::endl;

    		universe.gop.fence();

    	}
        universe.gop.fence();
    	subworld.gop.fence();

    }
}

#endif
