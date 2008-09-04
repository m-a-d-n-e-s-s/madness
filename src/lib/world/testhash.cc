#include <iostream>
#include <world/worldhash.h>

using namespace std;
using namespace HASH_MAP_NAMESPACE;

int main() {
   hash_map<int,int> h;

   h.insert(pair<int,int>(1,2));
   cout << h.find(1)->second << endl;

   return 0;
}
