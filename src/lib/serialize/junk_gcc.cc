#include <iostream>

using namespace std;

template <class T> 
struct typeinfo {
  static const int cookie = 255;
};

template <> struct typeinfo<int> {
  static const int cookie = 0;
};

int main() {
  cout << typeinfo<int>::cookie << endl; // Specialized
  cout << typeinfo<double>::cookie << endl; // Default
  return 0;
}

