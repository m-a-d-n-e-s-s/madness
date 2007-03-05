#include <iostream>
#include <vector>
using namespace std;

class A {
    int val;
public:
    A() : val(0) {
        cout << "Default constructor" << endl;
    };

    A(int val) : val(val) {
        cout << "Value constructor " << val << endl;
    };

    A(const A& a) : val(a.val) {
        cout << "Copy constructor " << val << endl;
    };

    A& operator=(const A& a) {
        cout << "Assignment " << val << " " << a.val << endl;
        val = a.val;
        return *this;
    };

    void set(int a) {val = a;};

    ~A() {
        cout << "Destructor" << endl;
    }
};

int main() {
    cout << "Making vector(3)" << endl;
    vector<A> v(3);
    cout << "Finished making vector" << endl;

    cout << "Assigning vector values" << endl;
    for (int i=0; i<(int)v.size(); i++) v[i].set(i+1);
    cout << "Finished assigning vector values" << endl;

    cout << "Making empty vector" << endl;
    vector<A> u;
    cout << "Finished making empty vector" << endl;

    cout << "Assigning to empty vector from vector(3)" << endl;
    u = v;
    cout << "Finished assigning to empty vector" << endl;

    cout << "Reassigning vector values" << endl;
    for (int i=0; i<(int)v.size(); i++) v[i].set(i+5);
    cout << "Finished reassigning vector values" << endl;

    cout << "Assigning to existing vector from vector(3)" << endl;
    u = v;
    cout << "Finished assigning to existing vector" << endl;

    cout << "Vector copy constructor from vector(3)" << endl;
    vector<A> p(v);
    cout << "Finished vector copy constructor from vector(3)" << endl;

    return 0;
}




        

