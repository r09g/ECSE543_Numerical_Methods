/******************************************************************************/
/* Name: A3.cpp                                                               */
/* Date: 2020/12/09                                                           */
/* Author: Raymond Yang                                                       */
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "Matrix.h"
#include "NLEq.h"

extern int FLAG;

using namespace std;

class A3{
    private:
        // questions
        void Q1();
        void Q2();

    public:
        A3(int question=-1);

};

// -----------------------------------------------------------------------------
// Constructor & Destructor
// -----------------------------------------------------------------------------

/*  Executes all questions  */
A3::A3(int question){
    std::cout << ">>> ECSE 543 Numerical Methods Assignment 3 <<<" << std::endl;
    switch(question){
        case 1:
            this->Q1();
            break;
        case 2:
            this->Q2();
            break;
        default:
            this->Q1();
            this->Q2();
            break;
    }
    return;
}

// -----------------------------------------------------------------------------
// Member Function
// -----------------------------------------------------------------------------

/*  Question 1  */
void A3::Q1(){

    cout << "\n--------------------" << endl;
    cout << "Solving A3 Q1..." << endl;
    cout << "--------------------" << endl;

    // ----- Part a -----
    {
        Matrix<> b6 = Matrix<>::read_mat("./data/A3/b6.csv");
        Matrix<> h6 = Matrix<>::read_mat("./data/A3/h6.csv");
        
        Poly_CF<> fit(h6, b6);
        Matrix<> coeffs = fit.Lagrange_WD();

        cout << "a) The polynomial coefficients are" << endl;
        coeffs.show();
    }

    // ----- Part b -----
    {
        Matrix<> b6 = Matrix<>::read_mat("./data/A3/b6s.csv");
        Matrix<> h6 = Matrix<>::read_mat("./data/A3/h6s.csv");

        Poly_CF<> fit(h6, b6);
        Matrix<> coeffs = fit.Lagrange_WD();

        cout << "b) The polynomial coefficients are" << endl;
        coeffs.show();
    }

    // ----- Part c -----
    {
        Matrix<> b6 = Matrix<>::read_mat("./data/A3/b6s.csv");
        Matrix<> h6 = Matrix<>::read_mat("./data/A3/h6s.csv");
        Matrix<> x = Matrix<>::read_mat("./data/A3/h_1c.csv");

        Poly_CF<> fit(h6, b6);
        Matrix<> y = fit.Hermite_3(x);
        y.write_mat("./data/A3/b_1c.csv");
        cout << "c) The results are stored in ./data/A3/b_1c.csv" << endl;
        
    }

    // ----- Part e -----
    {
        Matrix<> b = Matrix<>::read_mat("./data/A3/b.csv");
        Matrix<> h = Matrix<>::read_mat("./data/A3/h.csv");

        double A = 1e-4;
        double lg = 0.5e-2;
        double lc = 30e-2;
        double N = 800;
        double I = 10;
        double u0 = 1.25663706e-6;

        auto H = [&](double psi) -> double{
            double B = psi/A;
            int len = b.get_n_row()*b.get_n_col();
            if(B < b.get(0)){
                B = b.get(0);
            } else if(B > b.get(len - 1)) {
                B = b.get(len - 1);
            } else {}

            for(int i = 0; i < len - 1; i++){
                if(B >= b.get(i) && B <= b.get(i+1)){
                    double L1 = (B - b.get(i+1))/(b.get(i) - b.get(i+1));
                    double L2 = (B - b.get(i))/(b.get(i+1) - b.get(i));
                    return h.get(i)*L1 + h.get(i+1)*L2;
                }
            }

            return -1.0;  // impossible
        };

        auto dH = [&](double psi) -> double{
            double B = psi/A;
            int len = b.get_n_row()*b.get_n_col();
            if(B < b.get(0)){
                B = b.get(0);
            } else if(B > b.get(len - 1)) {
                B = b.get(len - 1);
            } else {}

            for(int i = 0; i < len - 1; i++){
                if(B >= b.get(i) && B <= b.get(i+1)){
                    return (h.get(i+1) - h.get(i))/(b.get(i+1) - b.get(i));
                }
            }

            return -1.0;  // impossible
        };   

        auto f = [&](double psi) -> double{
            return lg/u0/A*psi + lc*H(psi) - N*I;
        };

        auto df = [&](double psi) -> double{
            return lg/u0/A + lc/A*dH(psi);
        };

        Matrix<> results = NLEq::Newton_Raphson(f, df);
        cout << "e) The Newton-Raphson results are:" << endl;
        cout << "Psi = " << (results.get(0)) << endl;
        cout << "Iterations: " << (results.get(1)) << endl;
        cout << "Final Error: " << (results.get(2)) << endl;
    }

    // ----- Part f -----
    {
        Matrix<> b = Matrix<>::read_mat("./data/A3/b.csv");
        Matrix<> h = Matrix<>::read_mat("./data/A3/h.csv");

        double A = 1e-4;
        double lg = 0.5e-2;
        double lc = 30e-2;
        double N = 800;
        double I = 10;
        double u0 = 1.25663706e-6;

        auto H = [&](double psi) -> double {
            double B = psi/A;
            int len = b.get_n_row()*b.get_n_col();
            if(B < b.get(0)){
                B = b.get(0);
            } else if(B > b.get(len - 1)) {
                B = b.get(len - 1);
            } else {}

            for(int i = 0; i < len - 1; i++){
                if(B >= b.get(i) && B <= b.get(i+1)){
                    double L1 = (B - b.get(i+1))/(b.get(i) - b.get(i+1));
                    double L2 = (B - b.get(i))/(b.get(i+1) - b.get(i));
                    return h.get(i)*L1 + h.get(i+1)*L2;
                }
            }

            return -1.0;  // impossible
        };

        auto f = [&](double psi) -> double {
            return 1e-8*(lg/u0/A*psi + lc*H(psi) - N*I);
        };

        Matrix<> results = NLEq::successive_substitution(f);
        cout << "f) The Successive Substitution results are:" << endl;
        cout << "Psi = " << (results.get(0)) << endl;
        cout << "Iterations: " << (results.get(1)) << endl;
        cout << "Final Error: " << (results.get(2)) << endl;
    }

    cout << "\nA3 Q1 Solved." << endl;
    return;
}

/*  Question 2  */
void A3::Q2(){

    cout << "\n--------------------" << endl;
    cout << "Solving A3 Q2..." << endl;
    cout << "--------------------" << endl;

    // ----- Part b -----
    {
        double E = 200e-3;
        double R = 512;
        double IsA = 0.8e-6;
        double IsB = 1.1e-6;
        double kT_q = 25e-3;

        auto f = [&](Matrix<> v) -> Matrix<> {
            Matrix<> fout(2, 1);
            fout.set(0, 0, (E - v.get(0))/R - IsA*(exp((v.get(0) - v.get(1))/
                (kT_q)) - 1));
            fout.set(1, 0, IsA*(exp((v.get(0) - v.get(1))/(kT_q)) - 1) - 
                IsB*(exp(v.get(1)/(kT_q)) - 1));
            return fout;
        };

        auto df = [&](Matrix<> v) -> Matrix<> {
            Matrix<> dfout(2, 2);
            dfout.set(0, 0, (-1/R) - IsA/kT_q*exp((v.get(0) - v.get(1))/kT_q));
            dfout.set(0, 1, IsA/kT_q*exp((v.get(0) - v.get(1))/kT_q));
            dfout.set(1, 0, IsA/kT_q*exp((v.get(0) - v.get(1))/kT_q));
            dfout.set(1, 1, -IsA/kT_q*exp((v.get(0) - v.get(1))/kT_q) - 
                IsB/kT_q*exp(v.get(1)/kT_q));
            return dfout;
        };

        Matrix<> results = NLEq::Newton_Raphson(f, df, 2);
        
        cout << "b) The Newton-Raphson solutions are " << 
            "[node V, Iterations, Final Error]" << endl;
        results.show();
    }

    cout << "\nA3 Q2 Solved." << endl;
    cout << endl << ">>> End of ECSE 543 Assignment 3 <<<" << endl;
    return;
}











