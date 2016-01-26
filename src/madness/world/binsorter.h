
#ifndef MADNESS_WORLD_BINSORTER_H__INCLUDED
#define MADNESS_WORLD_BINSORTER_H__INCLUDED

#include <madness/world/MADworld.h>
#include <madness/world/world_object.h>
#include <madness/world/madness_exception.h>
#include <vector>

namespace madness {

    /// A parallel bin sort across MPI processes
    template <typename T, typename inserterT>
    class BinSorter : public WorldObject< BinSorter<T,inserterT> > {
        World* pworld;
        inserterT inserter;
        std::size_t bufsize;
        std::vector<T>* bins;
        
        void flush(int owner) {
            if (bins[owner].size())
                this->send(owner, &BinSorter<T,inserterT>::sorter,  bins[owner]);
            bins[owner].clear();
        }
        
        void sorter(const std::vector<T>& v) {
            for(typename std::vector<T>::const_iterator it = v.begin(); it != v.end(); ++it) {
                inserter(*it);
            }
        }
        
    public:
        /// Constructs the sorter object 

        /// @param[in] world The world object that must persist during the existence of this object
        /// @param[in] inserter User provides this routine to process an item of data on remote end
        /// @param[in] bufsize Size of bin (in units of T) ... default value is as large as possible
        BinSorter(World& world, inserterT inserter, int bufsize=0)
            : WorldObject<BinSorter>(world)
            , pworld(&world)
            , inserter(inserter)
            , bufsize(bufsize)
            , bins(new std::vector<T>[world.size()])
        {
            // bufsize ... max from AM buffer size is about 512K/sizeof(T)
            // bufsize ... max from total buffer use is about 1GB/sizeof(T)/P
            if (bufsize <= 0) this->bufsize = std::min((RMI::max_msg_len()-1024)/sizeof(T),(1<<30)/(world.size()*sizeof(T)));

            // for (int i=0; i<world.size(); i++) {
            //     bins[i].reserve(bufsize); // Not a good idea on large process counts unless truly all to all?
            // }

            //print("binsorter bufsize is", this->bufsize, this->bufsize*sizeof(T));
            WorldObject< BinSorter<T,inserterT> >::process_pending();
        }
        
        virtual ~BinSorter() {
            for (int i=0; i<pworld->size(); i++) {
                MADNESS_ASSERT(bins[i].size() == 0);
            }
            delete [] bins;
        }
        
        /// Invoke to complete the sort, flush all buffers, and ensure communication/processing is complete
        void finish() {
            for (int i=0; i<pworld->size(); i++) {
                flush(i);
            }
            pworld->gop.fence();
        }
        
        /// Application calls this to add a value to the bin for process p
        void insert(ProcessID p, const T& value) {
            bins[p].push_back(value);

            // More intelligent buffer management would look at total use and flush
            // largest buffers rather than using a fixed buffersize per process
            if (bins[p].size() >= bufsize) flush(p);
        }
    };
}

#endif // MADNESS_WORLD_BINSORTER_H__INCLUDED
