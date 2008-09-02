#include <iostream>
#include <ext/hash_map>

using namespace std;
using namespace __gnu_cxx;

int main() {
   hash_map<int,int> h;

   h.insert(pair<int,int>(1,2));
   cout << h.find(1)->second << endl;

   return 0;
}
