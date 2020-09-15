/******************************************************************************/
/* Name: main.cpp                                                             */
/* Date: 2020/09/10                                                           */
/* Author: Raymond Yang                                                       */
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <ctime>
#include <cassert>
#include <stdlib.h>
#include <time.h>
#include "Matrix.h"
#include "Basic.h"

using namespace std;

void test(){
    int num = 1000000;
    int* arr1 = new int[num];
    for(int i = 0; i < num; i++){
        arr1[i] = rand() % 100 - 50;
    }

    Matrix<int> mat1 = Matrix<int>(1000,1000,arr1);
    mat1.write_mat("test.csv");

    Matrix<double> mat2 = (Matrix<double>::read_mat("test.csv"));
    Matrix<int> mat3 = mat2.to_int();
    
    delete [] arr1;
    return;
}

int main(){

    srand(time(NULL));

    clock_t start = clock();

    test();

    double duration = (clock() - start) / (double)(CLOCKS_PER_SEC);
    cout << "Executed in " << duration << "s" << endl;

    return 0;
}
