/******************************************************************************/
/* Name: NLEq.h                                                               */
/* Date: 2020/12/09                                                           */
/* Author: Raymond Yang                                                       */
/* Description: This namespace is a group of nonlinear equation solvers       */
/******************************************************************************/

#ifndef __NLEq__
#define __NLEq__

#include <functional>
#include "Matrix.h"
#include "Matrix_Solver.h"


namespace NLEq{
    Matrix<> Newton_Raphson(std::function<double(double)> f, 
        std::function<double(double)> df, double tol = -1);

    Matrix<> successive_substitution(std::function<double(double)> f, 
        double tol = -1);

    Matrix<> Newton_Raphson(std::function<Matrix<>(Matrix<>)> f, 
        std::function<Matrix<>(Matrix<>)> df, int nvar, double tol = -1);
};

// -----------------------------------------------------------------------------

/*
    Newton-Raphson method for solving nonlinear equations

    Input:
    f: the objective function (double) -> (double)
    df: the derivative function (double) -> (double)
    tol: tolerance for stopping condition

    Returns: column vector with [solution, iterations, final error]
*/
Matrix<> NLEq::Newton_Raphson(std::function<double(double)> f, 
    std::function<double(double)> df, double tol){

    if(tol <= 0){
        tol = 1e-6;
    }

    double psi = 0;
    double itr = 0;
    double fval = f(psi);
    double dfval = df(psi);
    double f_init = fval;

    while(std::abs(fval/f_init) > tol){
        psi = psi - (fval/dfval);
        fval = f(psi);
        dfval = df(psi);
        itr++;
    }
    
    Matrix<> out(3, 1);
    out.set(0, psi);
    out.set(1, itr);
    out.set(2, fval);

    return out;
}

/*
    Successive Substitution method for solving nonlinear equations

    Input:
    f: the objective function (double) -> (double)
    tol: tolerance for stopping condition

    Returns: column vector with [solution, iterations, final error]
*/
Matrix<> NLEq::successive_substitution(std::function<double(double)> f, 
    double tol){

    if(tol <= 0){
        tol = 1e-6;
    }

    double psi = 0;
    double itr = 0;
    double fval = f(psi);
    double f_init = fval;

    while(std::abs(fval/f_init) > tol){
        psi = psi - (fval);
        fval = f(psi);
        itr++;
    }
    
    Matrix<> out(3, 1);
    out.set(0, psi);
    out.set(1, itr);
    out.set(2, fval);

    return out;
}

/*
    Newton-Raphson method for solving systems of nonlinear equations

    Input:
    f: the objective function (Matrix) -> (Matrix)
    df: the derivative function (Matrix) -> (Matrix)
    nvar: number of unknown variables
    tol: tolerance for stopping condition

    Returns: Matrix with columns [solution, iterations, final error]
*/
Matrix<> NLEq::Newton_Raphson(std::function<Matrix<>(Matrix<>)> f, 
    std::function<Matrix<>(Matrix<>)> df, int nvar, double tol){

    if(tol <= 0){
        tol = 1e-6;
    }



    Matrix<> v(nvar, 1);
    v.set(0, 0.1);
    double itr = 0;
    bool sw = false;

    Matrix<> fval = f(v);
    Matrix<> dfval = df(v);
    Matrix<> dv(nvar, 1);
    Matrix<> f_init = fval;

    try{
        Matrix<> temp(nvar, 1);
        Matrix_Solver::cholesky(dfval, &temp);
    } catch(const char* msg) {
        sw = true;
    }

    while(!(((fval.div(f_init)).abs()) < tol)){
        if(sw){
            Matrix_Solver::cholesky_solve(-1.0*dfval, fval, &dv);
        } else {
            Matrix_Solver::cholesky_solve(dfval, -1.0*fval, &dv);
        }
        v = v + dv;
        fval = f(v);
        dfval = df(v);
        itr++;
    }
    
    Matrix<> out(nvar, 3);
    out.set(0, 0, v);
    out.set(0, 1, itr);
    out.set(0, 2, fval);

    return out;
}


#endif