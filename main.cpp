/*
    Author: Raymond Yang
    Date: 2020/09/09
*/

#include <iostream>
#include <fstream>
#include <ctime>
#include "Matrix.h"
#include "Basic.h"

using namespace std;

int main(){
    clock_t start = clock();
    
    Matrix<double> mat1 = Matrix<double>::read_mat("test.csv");

    Matrix<double> mat2(mat1);
    mat2.transpose();

    Matrix<double> mat3 = mat1 * mat2;

    cout << mat3.mean() << endl;

    double duration = (clock() - start) / (double)(CLOCKS_PER_SEC);
    cout << "Executed in " << duration << "s" << endl;
    return 0;
}
