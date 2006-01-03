#include <iostream>
#include <vector>
using namespace std;

template <typename T>
class Ftypes {
public:
  enum {id=-1};
};

template<> class Ftypes<int> {public: enum {id=1};};
template<> class Ftypes<double> {public: enum {id=2};};

class F {
private:
  vector<double> d;
  vector<int> i;
public:
  template <class T> int get() {
    return Ftypes<T>::id;
  };
};

int main() {
  vector<int> ii;
  cout << ii.size() << endl;
  F f;
  cout << f.get<double>() << endl;
  cout << f.get<int>() << endl;
  return 0;
}


