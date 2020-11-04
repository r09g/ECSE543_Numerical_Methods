/******************************************************************************/
/* Name: A2.cpp                                                               */
/* Date: 2020/11/02                                                           */
/* Author: Raymond Yang                                                       */
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <chrono>
#include <stdlib.h>
#include "Matrix.h"
#include "FDM.h"

#define E0 8.85418782e-12  // vacuum permittivity
extern int FLAG;

using namespace std;

class A2{
    public:
    A2(int question=-1);

    // Questions
    void Q2();
    void Q3();

};

// -----------------------------------------------------------------------------
// Constructor & Destructor
// -----------------------------------------------------------------------------

/*  Executes all questions  */
A2::A2(int question){
    switch(question){
        case 2:
            this->Q2();
            break;
        case 3:
            this->Q3();
            break;
        default:
            this->Q2();
            this->Q3();
            break;
    }
    return;
}

// -----------------------------------------------------------------------------
// Member Function
// -----------------------------------------------------------------------------

/*  Question 2  */
void A2::Q2(){

    cout << "\n--------------------" << endl;
    cout << "Solving A2 Q2..." << endl;
    cout << "--------------------" << endl;

    // ----- Part c -----
    {
        // import solution
        Matrix<> s2d_sln = Matrix<>::read_mat("./data/A2/output.csv");
        Matrix<> S_con = Matrix<>::read_mat("./data/A2/S_con.csv");

        double h = 0.02;
        double width = 0.1;
        double V = 110;
        FDM<> fdm(width/h + 1, width/h + 1, h);
        fdm.set(0, 0, 0, width/h, 0, false);
        fdm.set(0, width/h, 0, 0, 0, false);
        fdm.set(0.08/h, width/h, 0.06/h, width/h, 110, false);

        for(int row = 0; row < s2d_sln.get_n_row(); row++){
            fdm.set_phi(s2d_sln.get(row, 1)/h, s2d_sln.get(row, 0)/h, 
                s2d_sln.get(row, 2));
        }

        double E = 4*E0*fdm.calc_W(S_con);
        double C = 2*E/(110*110);
        cout << "c) The total capacitance per unit length is: " << C << " F/m" 
            << endl; 
    }

    cout << "\nA2 Q2 Solved." << endl;
    return;
}

/*  Question 3  */
void A2::Q3(){
    
    cout << "\n--------------------" << endl;
    cout << "Solving A2 Q3..." << endl;
    cout << "--------------------" << endl;

    // ----- Part a -----
    {
        double width = 0.1;
        double h = 0.02; 

        FDM<> fdm(width/h + 1, width/h + 1, h);
        // set boundaries
        // lower Dirichlet bound
        fdm.set(0, 0, 0, width/h, 0, false);
        // left Dirichlet bound
        fdm.set(0, width/h, 0, 0, 0, false);
        // top right Dirichlet bound
        fdm.set(0.08/h, width/h, 0.06/h, width/h, 110, false);

        fdm.CG();
        // fdm.get_phi().show();
    }


    cout << "\nA2 Q3 Solved." << endl;
    cout << endl << ">>> End of ECSE 543 Assignment 2 <<<" << endl;
    return;
}













