#include <iostream>
#include "basic.h"

using namespace std;

int main(){
    int a[4] = {1,2,3,4};
    double b[4] = {1.1111,2,3,4};

    int sum = basic::sum(a);


    cout << sum << endl;

    return 0;
}
