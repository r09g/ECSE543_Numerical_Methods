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
#include "matrix_solver.h"
#include "LRN.h"

using namespace std;

int FLAG = 0;

void Q1(){

    Matrix<double> A = Matrix<double>::SSPD_mat(10);
    double x1[] = {1,2,3,4,5,6,7,8,9,10};
    Matrix<double> x(10,1,x1);
    Matrix<double> b = A*x; 

    Matrix_Solver::choleski_solve(&A, &b);

    for(int i = 0; i < 10; i++){
        if(std::abs(b.get(i) - x1[i]) > 1e-10){
            cout << b.get(i) << " _ " << x1[i] << endl;
        }
    }

    b.show();

    return;
}

int main(){
    cout << endl;
    srand(time(NULL));
    clock_t start = clock();

    try{
        // LRN<double> circuit("./data/A1Q1d_test_circuit_5.csv");
        double circ[] = {
            1,0,0,20,-10,
            1,2,0,10,0,
            1,3,0,10,0,
            2,3,0,30,0,
            2,0,0,30,0,
            3,0,0,30,0
        };
        LRN<double> circuit(4,6,circ);
        circuit.choleski_solve();
        circuit.get_v().show();
    }catch(const char* msg){
        cout << msg << endl;
    }

    cout << "\nERROR FLAG: " << FLAG << endl;

    double duration = (clock() - start) / (double)(CLOCKS_PER_SEC);
    cout << "Executed in " << duration << "s" << endl;
    return 0;
}
