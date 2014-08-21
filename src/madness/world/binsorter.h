#include <madness/world/world.h>
#include <madness/world/worldobj.h>
#include <madness/world/worldexc.h>
#include <vector>

namespace madness {

    /// A parallel bin sort across MPI processes
    template <typename T, typename inserterT>
    class BinSorter : public WorldObject< BinSorter<T,inserterT> > {
        World* pworld;
        inserterT inserter;
        const std::size_t bufsize;
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
        /// @param[in] bufize Size of bin (in units of T)
        BinSorter(World& world, inserterT inserter, int bufsize=1024)
            : WorldObject<BinSorter>(world)
            , pworld(&world)
            , inserter(inserter)
            , bufsize(bufsize)
            , bins(new std::vector<T>[world.size()])
        {
            for (int i=0; i<world.size(); i++) {
                bins[i].reserve(bufsize); // Maybe not a good idea on large process counts?
            }
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
        
        /// Application calls this to add a value to the bin for process owner
        void insert(ProcessID p, const T& value) {
            bins[p].push_back(value);
            if (bins[p].size() >= bufsize) flush(p);
        }
    };
}
      
    

