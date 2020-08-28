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


}

#endif
