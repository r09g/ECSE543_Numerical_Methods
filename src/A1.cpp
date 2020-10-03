/******************************************************************************/
/* Name: A1.cpp                                                               */
/* Date: 2020/10/03                                                           */
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

using namespace std;

extern int FLAG;

class A1{
    public:
    A1();

    // Questions
    void Q1();
    void Q2();
    void Q3();

};

// -----------------------------------------------------------------------------
// Constructor & Destructor
// -----------------------------------------------------------------------------

/*  Executes all questions  */
A1::A1(){
    this->Q1();
    this->Q2();
    this->Q3();
    return;
}

// -----------------------------------------------------------------------------
// Member Function
// -----------------------------------------------------------------------------

/*  Question 1  */
void A1::Q1(){

    cout << "\n--------------------" << endl;
    cout << "Solving A1 Q1..." << endl ;
    cout << "--------------------" << endl;

    // ----- Part c -----
    Matrix<int> results(1,9);
    int num_correct = 0;
    for(int size = 2; size <= 10; size++){ 
        Matrix<double> A = Matrix<double>::SSPD_mat(10);
        double x1[] = {1,2,3,4,5,6,7,8,9,10};
        Matrix<double> x(10,1,x1);
        Matrix<double> b = A*x; 

        Matrix_Solver::cholesky_solve(&A, &b);
        bool skip = false;
        for(int i = 0; i < size; i++){
            if(std::abs(b.get(i) - x1[i]) > 1e-10){
                skip = true;
                break;
            }
        }

        if(!skip){
            results.set(0, size - 2, 1);
            num_correct++;
        }
    }

    if(num_correct == 9){
        cout << "Cholesky & SSPD Test Passed." << endl;
    } else {
        cout << "Cholesky & SSPD Test Failed: " << (9 - num_correct) 
            << " errors." << endl;
        results.show();
    }

    // ----- Part d -----
    results = Matrix<int>(1,5);
    num_correct = 0;
    for(int i = 0; i < 5; i++){
        std::string filepath = std::string("./data/A1Q1d_test_circuit_")
            + std::to_string(i + 1) + std::string(".csv");
        LRN<double> cct(filepath);
        cct.cholesky_solve();
        Matrix<double> node_v = cct.get_v();

        switch(i){
            case 0:
                if(std::abs(node_v.get(0) - 5) < 1e-10){
                    num_correct++;
                    results.set(0, 0, 1);
                }
                break;
            case 1:
                if(std::abs(node_v.get(0) - 50) < 1e-10){
                    num_correct++;
                    results.set(0, 1, 1);
                }
                break;
            case 2:
                if(std::abs(node_v.get(0) - 55) < 1e-10){
                    num_correct++;
                    results.set(0, 2, 1);
                }
                break;
            case 3:
                if(std::abs(node_v.get(0) - 20) < 1e-10
                    && std::abs(node_v.get(1) - 35) < 1e-10){
                    num_correct++;
                    results.set(0, 3, 1);
                }
                break;
            default:
                if(std::abs(node_v.get(0) - 5) < 1e-10
                    && std::abs(node_v.get(1) - 3.75) < 1e-10
                    && std::abs(node_v.get(2) - 3.75) < 1e-10){
                    num_correct++;
                    results.set(0, 4, 1);
                }
        }

    }

    if(num_correct == 5){
        cout << "Circuit Solver Test Passed." << endl;
    } else {
        cout << "Circuit Solver Test Failed: " << (5 - num_correct) 
            << " errors." << endl;
        results.show();
    }

    cout << "\nA1 Q1 Solved." << endl;

    return;
}

/*  Question 2  */
void A1::Q2(){

    cout << "\n--------------------" << endl;
    cout << "Solving A1 Q2..." << endl;
    cout << "--------------------" << endl;

    // ----- Part a -----
    Matrix<double> req1(1,14);
    for(int N = 2; N <= 15; N++){
        LRN<double> cct;
        cct.nxn_res_mesh(N);
        cct.cholesky_solve();
        Matrix<double> node_V = cct.get_v();
        Matrix<double> branch_I = cct.get_i();
        req1.set(0, N - 2, -(node_V.get((N+1)*(N+1)-2))/(branch_I.get(0)));
    }

    cout << "Non-banded:" << endl;
    cout << "a) The equivalent resistances are:" << endl;
    req1.show();

    // ----- Part b -----
    Matrix<double> perf1(1,14);
    for(int N = 2; N <= 15; N++){
        LRN<double> cct;
        cct.nxn_res_mesh(N);

        clock_t start = clock();
        cct.cholesky_solve();
        clock_t end = clock();

        perf1.set(0, N - 2, (end - start)/(double)(CLOCKS_PER_SEC));
    }

    cout << "b) The performance is:" << endl;
    perf1.show();

    // ----- Part c -----
    Matrix<double> req2(1,14);
    Matrix<double> perf2(1,14);
    for(int N = 2; N <= 15; N++){
        LRN<double> cct;
        cct.nxn_res_mesh(N);

        clock_t start = clock();
        cct.cholesky_solve_banded(N+2);  // for this problem only: b = N+2
        clock_t end = clock();

        Matrix<double> node_V = cct.get_v();
        Matrix<double> branch_I = cct.get_i();
        req2.set(0, N - 2, -(node_V.get((N+1)*(N+1)-2))/(branch_I.get(0)));
        perf2.set(0, N - 2, (end - start)/(double)(CLOCKS_PER_SEC));
    }

    cout << "Banded:" << endl;
    cout << "c) The performance is:" << endl;
    perf2.show();

    // ----- Part c -----
    Matrix<double> req3(1,14);
    for(int N = 2; N <= 15; N++){
        LRN<double> cct;
        cct.nxn_res_mesh(N);
        cct.cholesky_solve_banded(N+2);  // for this problem only: b = N+2
        Matrix<double> node_V = cct.get_v();
        Matrix<double> branch_I = cct.get_i();
        req3.set(0, N - 2, -(node_V.get((N+1)*(N+1)-2))/(branch_I.get(0)));
    }
    
    cout << "d) The equivalent resistances are:" << endl;
    req3.show();

    cout << "\nA1 Q2 Solved." << endl;
    return;
}

/*  Question 3  */
void A1::Q3(){

    return;
}













