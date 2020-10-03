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
#include "Matrix.h"
#include "Basic.h"
#include "matrix_solver.h"
#include "LRN.h"

#include "A1.cpp"

using namespace std;

int FLAG = 0;

int main(){
    cout << endl;
    srand(time(NULL));
    clock_t tic = clock();
    cout << "ECSE 543 Numerical Methods Assignment 1" << endl;

    try{
        A1 a1 = A1();
    }catch(const char* msg){
        cout << msg << endl;
    }

    cout << "\nERROR FLAG: " << FLAG << endl;

    clock_t toc = clock();
    double duration = (toc - tic) / (double)(CLOCKS_PER_SEC);
    cout << "Executed in " << duration << "s" << endl;
    return 0;
}
