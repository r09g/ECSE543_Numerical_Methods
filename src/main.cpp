/******************************************************************************/
/* Name: main.cpp                                                             */
/* Date: 2020/09/10                                                           */
/* Author: Raymond Yang                                                       */
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <chrono>
#include <stdlib.h>
#include <math.h>
#include "Matrix.h"
#include "Basic.h"
#include "Matrix_Solver.h"
#include "LRN.h"
#include "FDM.h"
#include "Poly_CF.h"
#include "A1.cpp"
#include "A2.cpp"
#include "A3.cpp"
#include "NLEq.h"

using namespace std;

int FLAG = 0;  // to be used at task assignment level or higher

int main(){
    cout << endl;
    srand(time(NULL));
    auto start = chrono::high_resolution_clock::now(); 
    
    // solve assignment questions here
    try{
	A1 a1 = A1();
	A2 a2 = A2();
        A3 a3 = A3();

    }catch(const char* msg){
        FLAG -= 1;
        cout << msg << endl;
    }

    cout << "\nERROR FLAG: " << FLAG << endl;  // 0 = no error

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration
            <double, std::milli>(end - start).count()/1e3;
    cout << "Executed in " << duration << "s" << endl << endl;
    return EXIT_SUCCESS;
}
