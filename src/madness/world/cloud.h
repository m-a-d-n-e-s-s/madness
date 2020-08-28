
/**
 \file cloud.h
 \brief Declares the \c Cloud class for storing data and transfering them between worlds
 \ingroup world

*/


#ifndef SRC_MADNESS_WORLD_CLOUD_H_
#define SRC_MADNESS_WORLD_CLOUD_H_

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
	/// @param[in]	w	the universe world
	Cloud(madness::World& w) : container(w) {
		MADNESS_ASSERT(w.id()==0);
	}

	void set_debug(bool value) {
		debug=value;
	}

	template<typename T>
	void store(madness::World& world, const T& source, const int record) {
		if (debug)
			std::cout << "storing " << typeid(T).name() << " to record " << record << std::endl;

		// scope is important because of destruction ordering of world objects and fence
		{
			madness::archive::ContainerRecordOutputArchive ar(world,container,record);
			madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(world, ar);
			par & source;
		}
		if (dofence) world.gop.fence();
		if (debug) std::cout << "done with storing; container size " << container.size() << std::endl;
	}

	// call this load if your argument has a constructor taking world as argument (e.g. WorldContainer)
	template<typename T>
	typename std::enable_if<std::is_constructible<T,World&>::value, T>::type
	load(madness::World& world, const int record) {
		T target(world);
		load(world, target, record);
        return target;
	}


	// call this load for simple types
	template<typename T>
	typename std::enable_if<!std::is_constructible<T,World&>::value, T>::type
	load(madness::World& world, const int record) {
		T target;
		load(world, target, record);
        return target;
	}

	// call this load to pass in an already constructed argument
	template<typename T>
	void load(madness::World& world, T& target, const int record) {
        if (debug) std::cout << "loading " << typeid(T).name() << " to world " << world.id() << " from record " << record << std::endl;
        madness::archive::ContainerRecordInputArchive ar(world,container,record);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(world, ar);
        par & target;
        if (dofence) world.gop.fence();
	}

};

} /* namespace madness */

#endif /* SRC_MADNESS_WORLD_CLOUD_H_ */
