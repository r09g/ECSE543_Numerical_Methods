/******************************************************************************/
/* Name: main.cpp                                                             */
/* Date: 2020/09/10                                                           */
/* Author: Raymond Yang                                                       */
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <windows.h>
#include <stdlib.h>
#include <math.h>
#include "Matrix.h"
#include "Basic.h"
#include "matrix_solver.h"
#include "LRN.h"
#include "FDM.h"
#include "A1.cpp"

using namespace std;

int FLAG = 0;  // to be used at task assignment level or higher

int main(){
    cout << endl;
    srand(time(NULL));
    LARGE_INTEGER freq, tic, toc;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&tic);

    cout << ">>> ECSE 543 Numerical Methods Assignment 1 <<<" << endl;
    
    // solve assignment questions here
    try{
        A1 a1 = A1();
    }catch(const char* msg){
        FLAG -= 1;
        cout << msg << endl;
    }

    cout << "\nERROR FLAG: " << FLAG << endl;  // 0 = no error

    QueryPerformanceCounter(&toc);
    double duration = 1.0*(toc.QuadPart - tic.QuadPart)/freq.QuadPart;
    cout << "Executed in " << duration << "s" << endl << endl;
    return EXIT_SUCCESS;
}
