/*
    Author: Raymond Yang
    Date: 2020/09/09
*/

#include <iostream>
#include "Matrix.h"
#include "basic.h"

using namespace std;

int main(){

    int a[] = {1,2,3,2};
    int b[] = {0,0,0,0};
    Matrix<int> m1(2,2,a);
    Matrix<int> m2(2,2,b);
    Matrix<int> m3 = m1.div(m2);

    m3.show();
    return 0;
}
