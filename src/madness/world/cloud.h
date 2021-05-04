
/**
 \file cloud.h
 \brief Declares the \c Cloud class for storing data and transfering them between worlds
 \ingroup world

*/


#ifndef SRC_MADNESS_WORLD_CLOUD_H_
#define SRC_MADNESS_WORLD_CLOUD_H_


#include <madness/world/parallel_dc_archive.h>
#include<iomanip>

namespace madness {



/// cloud class

/// store and load data from the cloud into arbitrary worlds
///
/// Note that storing and load is always done only by rank=0 of the respective world and
/// the data is communicated between the processes of that world.
/// If data is to be collected into the universe from different subworlds (after forking worlds)
/// make sure to first load the subworlds data into a universe container, and do the reduction/
/// further handling afterwards.
///
/// Note that there must be a fence after the destruction of subworld containers, as in:
///
///  create subworlds
///  {
///      dcT(subworld)
///      do work
///  }
///  subworld.gop.fence();
class Cloud  {

	madness::WorldContainer<long,std::vector<unsigned char> > container;
	bool debug=false;
	bool dofence=true;

public:
	mutable std::atomic<long> reading_time;	// in ms
	mutable std::atomic<long> writing_time;	// in ms
    using keyT=madness::archive::ContainerRecordOutputArchive::keyT;


    template<typename T>
    struct recordlist {
        std::list<T> list;
        recordlist() : list() {};
        recordlist(const T& key) : list{key} {};
        recordlist(const recordlist& other) : list(other.list) {};

        recordlist& operator+=(const recordlist& list2) {
            for (auto& l2 : list2.list) list.push_back(l2);
            return *this;
        }

        recordlist& operator+=(const T& key) {
            list.push_back(key);
            return *this;
        }

        T pop_front_and_return() {
            T key=list.front();
            list.pop_front();
            return key;
        }

        std::size_t size() const {
            return list.size();
        }

        friend std::ostream& operator<<(std::ostream& os, const recordlist& arg) {
            using namespace madness::operators;
            os << arg.list ;
            return os;
        }

    };
    using recordlistT = recordlist<keyT>;
    using lookuplistT = std::map<keyT,std::string>;

    lookuplistT lookuplist;
    template<typename T>
    void add_to_lookuplist(const keyT& key, const T& arg) {
        lookuplist[key]=typeid(arg).name();
    }
    void print_lookuplist() const {
        print("lookuplist of the cloud");
        print("record/key    typeid");
        for (const auto l : lookuplist) printf("%12zu %s\n",l.first,l.second.c_str());
    }

private:
    template <typename T> struct is_vector: std::false_type {};
    template <typename T> struct is_vector<std::vector<T>>: std::true_type {};

    template <typename> struct is_tuple: std::false_type {};
    template <typename ...T> struct is_tuple<std::tuple<T...>>: std::true_type {};

    template <typename> struct is_madness_function_vector: std::false_type {};
    template <typename T, std::size_t NDIM> struct is_madness_function_vector<std::vector<Function<T,NDIM>>> : std::true_type {};

    template<typename T> struct is_world_object : std::false_type {};
    template<typename T> struct is_world_object<madness::WorldObject<T>> : std::true_type {};

    template<typename T> struct is_funcimpl_ptr : std::false_type {};
    template<typename T, std::size_t NDIM> struct is_funcimpl_ptr<madness::FunctionImpl<T,NDIM>*> : std::true_type {};

//    template<typename T> using is_world_constructible=std::is_constructible<T, World&>;
//    template<typename T> using is_vector_of_world_constructible=is_vector<std::vector<std::is_constructible<T, World&>>>;



    struct storetimer {
		double cpu0;
		std::atomic<long>& wtime;
		storetimer(std::atomic<long>& writetime) : cpu0(cpu_time()), wtime(writetime) {}
		~storetimer() {
			wtime+=long((cpu_time()-cpu0)*1000l);
		}
	};
	struct loadtimer {
		double cpu0;
		std::atomic<long>& rtime;
		loadtimer(std::atomic<long>& readtime) : cpu0(cpu_time()), rtime(readtime)  {}
		~loadtimer() {
			rtime+=long((cpu_time()-cpu0)*1000l);
		}
	};


public:
	/// @param[in]	universe	the universe world
	Cloud(madness::World& universe) : container(universe), reading_time(0l), writing_time(0l){}

	void set_debug(bool value) {
		debug=value;
	}

	void set_fence(bool value) {
		dofence=value;
	}

	void print_timings(World& universe) const {
		double rtime=double(reading_time);
		double wtime=double(writing_time);
		universe.gop.sum(rtime);
		universe.gop.sum(wtime);
		if (universe.rank()==0) {
		    double precision = std::cout.precision();
            std::cout << std::fixed << std::setprecision(1);
			print("cloud storing cpu time", wtime*0.001);
			print("cloud reading cpu time", rtime*0.001,std::defaultfloat);
			std::cout << std::setprecision(precision) << std::scientific;
		}
	}


    /// compute the record of the argument

    /// distinguish 3 cases for the argument: Function, vector<Function>, everything else
    /// TODO: fix hash value for vector of functions,
    template<typename T>
    static keyT compute_record(const T& arg) {
        if constexpr (is_madness_function_vector<T>::value) {
            print("hurray");
        }
        if constexpr (is_madness_function<T>::value) {
//            print("hashing madness::Function");
            return hash_value(arg.get_impl()->id());
        }
        else if constexpr (is_world_object<T>::value) {
//            print("hashing WorldObject");
            return hash_value(arg.id());
        }
        else if constexpr (is_funcimpl_ptr<T>::value) {
//            print("hashing pointer to FunctionImpl");
            return hash_value(arg->id());
        }
        else if constexpr (is_vector<T>::value) {
            if constexpr (is_madness_function<typename T::value_type>::value) {
//                print("hashing vector of madness::Function");
                return hash_value(arg.front().get_impl()->id());
            }
        }
        else {
//            print("hashing simple ", typeid(T).name());
            return keyT(hash_value(arg));
        }
    }


    template<typename T>
	recordlistT store(madness::World& world, const std::vector<T>& source) {
		storetimer t(writing_time);
        keyT record= compute_record(source);
		if (debug)
			std::cout << "storing vector of " << typeid(T).name() << " of size " << source.size() << " to record " << record << std::endl;

		// scope is important because of destruction ordering of world objects and fence
		{
			madness::archive::ContainerRecordOutputArchive ar(world,container,record);
			madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(world, ar);
			par & source.size();
            add_to_lookuplist(record,source);
			for (const auto& s : source) par & s;
		}
		if (dofence) world.gop.fence();
		if (debug) std::cout << "done with vector storing; container size " << container.size() << std::endl;
        return recordlistT{record};
	}

    template<typename T, std::size_t NDIM>
    recordlistT store(madness::World& world, const Function<T,NDIM>& source) {
        storetimer t(writing_time);
        keyT record=compute_record(source);
        if (debug)
            std::cout << "storing object of " << typeid(T).name() << " of size " << source.size() << " to record " << record << std::endl;

        // scope is important because of destruction ordering of world objects and fence
        {
            madness::archive::ContainerRecordOutputArchive ar(world,container,record);
            madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(world, ar);
            par & source;
            add_to_lookuplist(record,source);
        }
        if (dofence) world.gop.fence();
        if (debug) std::cout << "done with function storing; container size " << container.size() << std::endl;
        return recordlistT{record};
    }

    template<typename T, std::size_t NDIM>
    recordlistT store(madness::World& world, const std::vector<Function<T,NDIM>>& source) {
        storetimer t(writing_time);
        keyT record= compute_record(source.size());
        recordlistT l(record);
        if (debug)
            std::cout << "storing vector of " << typeid(T).name() << " of size " << source.size() << " to record " << record << std::endl;

        // scope is important because of destruction ordering of world objects and fence
        {
            madness::archive::ContainerRecordOutputArchive ar(world,container,record);
            madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(world, ar);
            par & source.size();
            add_to_lookuplist(record,source.size());
        }
        for (auto s : source) l+=store(world,s);
        if (dofence) world.gop.fence();
        if (debug) std::cout << "done with vector storing; container size " << container.size() << std::endl;
        return l;
    }

    template<typename T>
    recordlistT store(madness::World& world, const T& source) {
        storetimer t(writing_time);
        auto record= compute_record(source);
        if (debug)
            std::cout << "storing object of " << typeid(T).name() << " to record " << record << std::endl;

        // scope is important because of destruction ordering of world objects and fence
        {
            madness::archive::ContainerRecordOutputArchive ar(world,container,record);
            madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(world, ar);
            par & source;
            add_to_lookuplist(record,source);
        }
        if (dofence) world.gop.fence();
        if (debug) std::cout << "done with vector storing; container size " << container.size() << std::endl;
        return recordlistT{record};
    }

    /// store a tuple in multiple records
    template<typename... Ts>
    recordlistT store(World& world, const std::tuple<Ts...>& input) {
        recordlistT v;
        auto storeaway=[&](const auto& arg) {
            v+=this->store(world,arg);
        };
        auto l=[&](Ts const &... arg) {
            ((storeaway(arg)),...);
        };
        std::apply(l,input);
        return v;
    }

    // call this to store a vector of world-constructible objects
    template<typename vecT>
    typename std::enable_if<is_vector<vecT>::value && std::is_constructible<typename vecT::value_type,World&>::value, recordlistT>::type
    store(madness::World& world, const vecT& target) const {
        loadtimer t(reading_time);
        keyT record=compute_record(target);
        recordlistT recordlist(record);
        typedef typename vecT::value_type T;
        if (debug) std::cout << "loading vector of type " << typeid(T).name() << " to world " << world.id() << " from record " << record << std::endl;
        madness::archive::ContainerRecordInputArchive ar(world,container,record);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(world, ar);
        std::size_t sz;
        par & sz;
        add_to_lookuplist(record,target);
        for (auto& i : target) recordlist+=store(world,i);
        if (dofence) world.gop.fence();
    }

	// call this load a vector of world-constructible objects
	template<typename vecT>
	typename std::enable_if<is_vector<vecT>::value && std::is_constructible<typename vecT::value_type,World&>::value, vecT>::type
	load(madness::World& world, recordlistT& recordlist) const {
		typedef typename vecT::value_type T;
        keyT record=recordlist.pop_front_and_return();
        if (debug) std::cout << "loading vector of type " << typeid(T).name() << " to world " << world.id() << " from record " << record << std::endl;
        madness::archive::ContainerRecordInputArchive ar(world,container,record);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(world, ar);
        std::size_t sz;
        par & sz;
        vecT target(sz);
        for (std::size_t i=0; i<sz; ++i) target[i]=load<T>(world,recordlist);
        if (dofence) world.gop.fence();
        return target;
	}

    // call this load a vector of simple objects
    template<typename vecT>
    typename std::enable_if<is_vector<vecT>::value && !std::is_constructible<typename vecT::value_type,World&>::value, vecT>::type
    load(madness::World& world, recordlistT& recordlist) const {
        typedef typename vecT::value_type T;
        keyT record=recordlist.pop_front_and_return();
        if (debug) std::cout << "loading vector of type " << typeid(T).name() << " to world " << world.id() << " from record " << record << std::endl;
        madness::archive::ContainerRecordInputArchive ar(world,container,record);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(world, ar);
        std::size_t sz;
        par & sz;
        vecT target(sz);
        for (std::size_t i=0; i<sz; ++i) par & target[i];
        if (dofence) world.gop.fence();
        return target;
    }

    // call this load a simple object
    template<typename T>
    typename std::enable_if<!is_vector<T>::value && !is_tuple<T>::value && !std::is_constructible<T, World&>::value, T>::type
    load(madness::World& world, recordlistT& recordlist) const {
//        ShowType<T> dummy;
        loadtimer t(reading_time);
        keyT record=recordlist.pop_front_and_return();
        if (debug) std::cout << "loading type " << typeid(T).name() << " to world " << world.id() << " from record " << record << std::endl;
        madness::archive::ContainerRecordInputArchive ar(world,container,record);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(world, ar);
        T target;
        par & target;
        if (dofence) world.gop.fence();
        return target;
    }

    // call this load a world-constructible object
    template<typename T>
    typename std::enable_if<!is_vector<T>::value && !is_tuple<T>::value && std::is_constructible<T, World&>::value, T>::type
    load(madness::World& world, recordlistT& recordlist) const {
        loadtimer t(reading_time);
        keyT record=recordlist.pop_front_and_return();
        if (debug) std::cout << "loading type " << typeid(T).name() << " to world " << world.id() << " from record " << record << std::endl;
        madness::archive::ContainerRecordInputArchive ar(world,container,record);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(world, ar);
        T target(world);
        par & target;
        if (dofence) world.gop.fence();
        return target;
    }

    template<typename T>
    typename std::enable_if<is_tuple<T>::value, T>::type
    load(madness::World& world, recordlistT& recordlist) const {
        if (debug) std::cout << "loading tuple of type " << typeid(T).name() << " to world " << world.id() << std::endl;
        T target;
        std::apply([&](auto &&... args) { ((args=load<typename std::remove_reference<decltype(args)>::type>(world,recordlist)), ...); }, target);
        return target;
    }

};

} /* namespace madness */

#endif /* SRC_MADNESS_WORLD_CLOUD_H_ */
