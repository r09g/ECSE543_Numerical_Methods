#include <iostream>
#include <memory>
#include "basic.h"
#include "Matrix.h"

using namespace std;

int main(){
    unique_ptr<int> ptr = make_unique<int>(1);
    *ptr = 1;
    cout << *ptr << endl;

    return 0;
}
