
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

	/// vector helper function
	template <typename T, typename _ = void>
	struct is_vector {
	    static const bool value = false;
	};
	template <typename T>
	struct is_vector< T, typename std::enable_if<
		  std::is_same<T,std::vector< typename T::value_type,
			   typename T::allocator_type > >::value>::type>
	{
	    static const bool value = true;
	};

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


	template<typename T>
	void store(madness::World& world, const std::vector<T>& source, const int record) {
		if (debug)
			std::cout << "storing vector of " << typeid(T).name() << " of size " << source.size() << " to record " << record << std::endl;

		// scope is important because of destruction ordering of world objects and fence
		{
			madness::archive::ContainerRecordOutputArchive ar(world,container,record);
			madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(world, ar);
			par & source.size();
			for (const auto& s : source) par & s;
		}
		if (dofence) world.gop.fence();
		if (debug) std::cout << "done with vector storing; container size " << container.size() << std::endl;
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

	// call this load to pass in an already constructed argument (not a vector)
	template<typename T>
	typename std::enable_if<!is_vector<T>::value, void>::type
	load(madness::World& world, T& target, const int record) {
        if (debug) std::cout << "loading " << typeid(T).name() << " to world " << world.id() << " from record " << record << std::endl;
        madness::archive::ContainerRecordInputArchive ar(world,container,record);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(world, ar);
        par & target;
        if (dofence) world.gop.fence();
	}

	// call this load to pass in an uninitialized vector, default-constructible
	template<typename vecT>
	typename std::enable_if<is_vector<vecT>::value && !std::is_constructible<typename vecT::value_type,World&>::value, void>::type
	load(madness::World& world, vecT& target, const int record) {
		typedef typename vecT::value_type T;
        if (debug) std::cout << "loading vector of type " << typeid(T).name() << " to world " << world.id() << " from record " << record << std::endl;
        madness::archive::ContainerRecordInputArchive ar(world,container,record);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(world, ar);
        std::size_t sz;
        par & sz;
        target.resize(sz);
        for (std::size_t i=0; i<sz; ++i) {
        	target[i]=T();
            par & target[i];
        }
        if (dofence) world.gop.fence();
	}

	// call this load to pass in an uninitialized vector, world-constructible
	template<typename vecT>
	typename std::enable_if<is_vector<vecT>::value && std::is_constructible<typename vecT::value_type,World&>::value, void>::type
	load(madness::World& world, vecT& target, const int record) {
		typedef typename vecT::value_type T;
        if (debug) std::cout << "loading vector of type " << typeid(T).name() << " to world " << world.id() << " from record " << record << std::endl;
        madness::archive::ContainerRecordInputArchive ar(world,container,record);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(world, ar);
        std::size_t sz;
        par & sz;
        target.resize(sz);
        for (std::size_t i=0; i<sz; ++i) {
        	target[i]=T(world);
            par & target[i];
        }
        if (dofence) world.gop.fence();
	}

};

} /* namespace madness */

#endif /* SRC_MADNESS_WORLD_CLOUD_H_ */
