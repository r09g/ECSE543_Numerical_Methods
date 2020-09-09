/*
    Author: Raymond Yang
    Date: 2020/09/09
*/

#include <iostream>
#include "Matrix.h"
#include "basic.h"

using namespace std;

int main(){

    int a[] = {1,2,3,2,3,4,3,4,5};
    int b[] = {-1,-1};
    Matrix<int> m1(3,3,a);
    Matrix<int> m2(1,2,b);
    Matrix<int> m3 = m1*m2;

    m3.show();
    return 0;
}
