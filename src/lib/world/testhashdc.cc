#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <world/world.h>

using namespace madness;
using namespace std;

// User defined key class that we don't want to modify
struct Key {
    int k;
    Key() : k(-1) {}
    Key(int k) : k(k) {}

    bool operator==(const Key& b) const {
        return k==b.k;
    }
};

ostream& operator<<(ostream&s, const Key& key) {
    s << "Key(" << key.k << ")";
    return s;
}

// Make the key serialiable using non-intrusive mechanism
namespace madness {
    namespace archive {
        template <class Archive>
        struct ArchiveSerializeImpl<Archive,Key> {
            static inline void serialize(const Archive& ar, Key& obj) {
                ar & obj.k;
            }
        };
    }
}

// Make the key hashable using non-intrusive mechanism
unsigned int hash(const Key& key) {
    return key.k;
}

int main(int argc, char** argv) {
    initialize(argc,argv);
    World world(MPI::COMM_WORLD);

    WorldContainer<Key,double> fred(world);

    fred.replace(Key(99),99.0);

    cout << fred.find(Key(99)) << endl;

    finalize();
    return 0;
}
