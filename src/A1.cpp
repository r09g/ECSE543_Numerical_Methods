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
    A1(int question=-1);

    // Questions
    void Q1();
    void Q2();
    void Q3();

};

// -----------------------------------------------------------------------------
// Constructor & Destructor
// -----------------------------------------------------------------------------

/*  Executes all questions  */
A1::A1(int question){
    switch(question){
        case 1:
            this->Q1();
            break;
        case 2:
            this->Q2();
            break;
        case 3:
            this->Q3();
            break;
        default:
            this->Q1();
            this->Q2();
            this->Q3();
            break;
    }
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
    {
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
    }

    // ----- Part d -----
    {
        Matrix<int> results = Matrix<int>(1,5);
        int num_correct = 0;
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
    {
        Matrix<double> req1(1,14);
        for(int N = 2; N <= 15; N++){
            LRN<double> cct;
            cct.nxn_res_mesh(N);
            cct.cholesky_solve();
            Matrix<double> node_V = cct.get_v();
            Matrix<double> branch_I = cct.get_i();
            req1.set(0, N - 2, -(node_V.get((N+1)*(N+1)-2))/(branch_I.get(0)));
        }

        cout << "---> Non-banded:" << endl;
        cout << "a)\nThe equivalent resistances are:" << endl;
        req1.show();
    }

    // ----- Part b -----
    {
        Matrix<double> perf1(1,14);
        for(int N = 2; N <= 15; N++){
            LRN<double> cct;
            cct.nxn_res_mesh(N);

            clock_t start = clock();
            cct.cholesky_solve();
            clock_t end = clock();

            perf1.set(0, N - 2, (end - start)/(double)(CLOCKS_PER_SEC));
        }

        cout << "b)\nThe performance is:" << endl;
        perf1.show();
    }

    // ----- Part c -----
    {
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

        cout << "---> Banded:" << endl;
        cout << "c)\nThe performance is:" << endl;
        perf2.show();
    }

    // ----- Part c -----
    {
        Matrix<double> req3(1,14);
        for(int N = 2; N <= 15; N++){
            LRN<double> cct;
            cct.nxn_res_mesh(N);
            cct.cholesky_solve_banded(N+2);  // for this problem only: b = N+2
            Matrix<double> node_V = cct.get_v();
            Matrix<double> branch_I = cct.get_i();
            req3.set(0, N - 2, -(node_V.get((N+1)*(N+1)-2))/(branch_I.get(0)));
        }
        
        cout << "d)\nThe equivalent resistances are:" << endl;
        req3.show();
    }

    cout << "\nA1 Q2 Solved." << endl;
    return;
}

/*  Question 3  */
void A1::Q3(){
    
    cout << "\n--------------------" << endl;
    cout << "Solving A1 Q3..." << endl;
    cout << "--------------------" << endl;

    // ----- Part a -----
    {
        // set up problem
        double width = 0.1;
        double h = 0.000625; 
        double omega = 1.5;
        double tol = 1e-5;

        FDM<> fdm(width/h + 1, width/h + 1, h);
        // set boundaries
        // lower Dirichlet bound
        fdm.set(0, 0, 0, width/h, 0, false);
        // left Dirichlet bound
        fdm.set(0, width/h, 0, 0, 0, false);
        // top right Dirichlet bound
        fdm.set(0.08/h, width/h, 0.06/h, width/h, 110, false);
        
        // run solver
        fdm.SOR(omega, tol);
        // fdm.get_phi().write_mat("solution.csv");
        
        cout << "---> SOR:" << endl;
        cout << "a)\nFor h = " << h << " and omega = " << omega << ":" << endl;
        cout << "Nodes: " << (fdm.get_n_row()*fdm.get_n_col()) << endl;
        cout << "Iterations: " << fdm.num_itr << std::endl;
        cout << "Time taken: " << fdm.duration << "s" << std::endl;
    }

    // ----- Part b -----
    {
        // set up problem
        double width = 0.1;
        double h = 0.02; 
        double tol = 1e-5;
        Matrix<> omega_vals(1, 10);
        Matrix<> iterations(1, 10);
        Matrix<> node_phi(1, 10);

        for(int i = 0; i < 10; i++){
            double omega = 1.05 + i*0.1;
            FDM<> fdm(width/h + 1, width/h + 1, h);
            // set boundaries
            fdm.set(0, 0, 0, width/h, 0, false);
            fdm.set(0, width/h, 0, 0, 0, false);
            fdm.set(0.08/h, width/h, 0.06/h, width/h, 110, false);
            
            try{
                fdm.SOR(omega, tol);
            }catch(const char* msg){
                FLAG -= 1;
                cout << endl << msg  << endl << endl;
            }
            
            omega_vals.set(0, i, omega);
            iterations.set(0, i, fdm.num_itr);
            node_phi.set(0, i, fdm.get_phi(0.06/h, 0.04/h));
        }
        
        cout << "b)\nFor h = " << h << endl;
        cout << "The omega values are:" << endl;
        omega_vals.show();
        cout << "The number of iterations are:" << endl;
        iterations.show();
        cout << "The node potential at (0.06, 0.04) is:" << endl;
        node_phi.show();
    }

    // ----- Part c -----
    {
        // set up problem
        double width = 0.1;
        double omega = 1.35; 
        double tol = 1e-5;
        Matrix<> h_vals(1,5);
        Matrix<> iterations = Matrix<>(1,5);
        Matrix<> node_phi = Matrix<>(1,5);

        for(int i = 0; i < 5; i++){
            double h = (double)(0.02)/pow(2,i);
            FDM<> fdm(width/h + 1, width/h + 1, h);
            // set boundaries
            fdm.set(0, 0, 0, width/h, 0, false);
            fdm.set(0, width/h, 0, 0, 0, false);
            fdm.set(0.08/h, width/h, 0.06/h, width/h, 110, false);
            
            try{
                fdm.SOR(omega, tol);
            }catch(const char* msg){
                FLAG -= 1;
                cout << endl << msg  << endl << endl;
            }
            
            h_vals.set(0, i, 1/h);
            iterations.set(0, i, fdm.num_itr);
            node_phi.set(0, i, fdm.get_phi(0.06/h, 0.04/h));
        }
        
        cout << "c)\nFor omega = " << omega << endl;
        cout << "The 1/h values are:" << endl;
        h_vals.show();
        cout << "The number of iterations are:" << endl;
        iterations.show();
        cout << "The node potential at (0.06, 0.04) is:" << endl;
        node_phi.show();
    }

    // ----- Part d -----
    {
        // set up problem
        double width = 0.1;
        double tol = 1e-5;
        Matrix<> h_vals = Matrix<>(1,5);
        Matrix<> iterations = Matrix<>(1,5);
        Matrix<> node_phi = Matrix<>(1,5);

        for(int i = 0; i < 5; i++){
            double h = (double)(0.02)/pow(2,i);
            FDM<> fdm(width/h + 1, width/h + 1, h);
            // set boundaries
            fdm.set(0, 0, 0, width/h, 0, false);
            fdm.set(0, width/h, 0, 0, 0, false);
            fdm.set(0.08/h, width/h, 0.06/h, width/h, 110, false);
            
            try{
                fdm.jacobi(tol);
            }catch(const char* msg){
                FLAG -= 1;
                cout << endl << msg << endl << endl;
            }
            
            h_vals.set(0, i, 1/h);
            iterations.set(0, i, fdm.num_itr);
            node_phi.set(0, i, fdm.get_phi(0.06/h, 0.04/h));
        }
        
        cout << "---> Jacobi:" << endl;
        cout << "d)\nThe 1/h values are:" << endl;
        h_vals.show();
        cout << "The number of iterations are:" << endl;
        iterations.show();
        cout << "The node potential at (0.06, 0.04) is:" << endl;
        node_phi.show();
    }

    // ----- Part e -----
    {
        // set up problem
        double width = 0.1;
        double omega = 1.35;
        double tol = 1e-5;

        double v_values[] = 
            {0,0.03,0.0525,0.0575,0.059,0.06,0.061,0.0625,0.0675,0.08,0.1};
        double h_values[] = 
            {0,0.02,0.0325,0.0375,0.039,0.04,0.041,0.0425,0.055,0.08,0.1};
        Matrix<> h_lines(1,11,h_values);
        Matrix<> v_lines(1,11,v_values);

        FDM<> fdm(h_lines, v_lines);
        fdm.set(0, 0, 0, 10, 0, false);
        fdm.set(0, 10, 0, 0, 0, false);
        fdm.set(9, 10, 5, 10, 110, false);
        try{    
            fdm.SOR(omega, tol);
        } catch(const char* msg){
            FLAG -= 1;
            cout << endl << msg << endl << endl;
        }

        cout << "---> SOR:" << endl;
        cout << "e)\nFor omega = " << omega << ":" << endl;
        cout << "Nodes: " << (fdm.get_n_row()*fdm.get_n_col()) << endl;
        cout << "Iterations: " << fdm.num_itr << std::endl;
        cout << "Time taken: " << fdm.duration << "s" << std::endl;
        cout << "The node potential at (0.06, 0.04) is:" 
            << fdm.get_phi(5,5) << endl;

        cout << "\nA1 Q3 Solved." << endl;
    }

    cout << endl << ">>> End of ECSE 543 Assignment 1 <<<" << endl;
    return;
}












