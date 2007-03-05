#include <map>
#include <iostream>
#include <misc/print.h>
#include <misc/communicator.h>
#include <mra/mra.h>
#include <misc/shared_ptr.h>

namespace madness {    
    namespace distmap {
        /// Key<n> wraps an array of ulongs for use as key into std::map 
        template <int n>
        class Key {
            typedef unsigned long ulong;
            ulong key[n];
            public:
                Key() {};
                Key(const ulong k[n]) {
                    for (int i=0; i<n; i++) key[i] = k[i];
                };
                Key(const Key<n>& other) {
                    for (int i=0; i<n; i++) key[i] = other.key[i];
                };
                Key<n>& operator=(const Key<n>& other) {
                    for (int i=0; i<n; i++) key[i] = other.key[i];
                    return *this;
                };
                bool operator<(const Key<n>& other) const {
                    for (int i=0; i<n; i++) {
                        if (key[i] < other.key[i]) return true;
                        else if (key[i] > other.key[i]) return false;
                    }
                    return false; // They are equal
                };
        };    
        
        class DistMapData {
        private:
            std::vector<bool> a;
            ProcessID rank;
        public:
            DistMapData() : a(madness::FunctionNode::size), rank(-1) {};
            DistMapData(ProcessID rank, std::vector<bool> a) : a(a), rank(rank) {};
            DistMapData(const DistMapData& d) : a(d.a), rank(d.rank) {};
            inline ProcessID get_rank() {return rank;};
            inline bool isactive(int ind) {return a[ind];};
        };
        
        typedef SharedPtr<DistMapData> DistMapDataPtr;
        typedef Key<4> key4T;
        typedef std::map< Key<4>, DistMapDataPtr > map4T;
        typedef std::pair< Key<4>, DistMapDataPtr > pair4T;
        typedef unsigned long ulong;
        
        map4T dir4;
        void clear4() {dir4.clear();};
        
        void set4(const ulong k[4], ProcessID rank, std::vector<bool>& a) {
            dir4.insert(pair4T(key4T(k),DistMapDataPtr(new DistMapData(rank,a))));
        }
        
        /// Find (n,l) of closest parent with coefficients and owning MPI process
        
        /// This is the C++ equivalent of the old Python routine sock_it_to_me.
        ///
        /// k contains [n,lx,ly,lz,ind] and the client wants coeffcients
        /// from either this node or a parent.  If such a node exists, k is set
        /// to the located values [n',lx',ly',lz',ind] and the MPI rank of the
        /// the owning process is returned.  Otherwise, -1 is returned to
        /// indicate that no coefficients were found which implies that the
        /// tree is locally more deeply refined.
        ProcessID get4_with_coeff(ulong k[5]) {
	    int n=k[0];
            ulong lx=k[1], ly=k[2], lz=k[3];
            int ind = int(k[4]);
            while (n>=0) {
                map4T::iterator test = dir4.find(key4T(k));
                if (test == dir4.end() || !test->second->isactive(ind)) {
                    n--;
                    lx = lx >> 1;
                    ly = ly >> 1;
                    lz = lz >> 1;
                }
                else {
                    k[0]=n; k[1]=lx; k[2]=ly; k[3]=lz;
                    return test->second->get_rank();
                }
            }
            return -1;
        }

        /// Find MPI process that logically owns (n,l)
        
        /// k contains [n,lx,ly,lz] and the client wants to locate the 
        /// MPI process that logically owns it, even if it does not 
        /// yet exist.
        ProcessID get4(ulong k[5]) {
	    int  n=k[0];
            ulong lx=k[1], ly=k[2], lz=k[3];
            while (n>=0) {
                map4T::iterator test = dir4.find(key4T(k));
                if (test == dir4.end()) {
                    n--;
                    lx = lx >> 1;
                    ly = ly >> 1;
                    lz = lz >> 1;
                }
                else {
                    return test->second->get_rank();
                }
            }
            return -1;
        }
    }
}

using namespace madness;
using namespace std;

int main(int argc, char** argv) {
    Communicator comm = startup(argc,argv);
    int dir4_set = comm.am_register();

    return 0;
}
