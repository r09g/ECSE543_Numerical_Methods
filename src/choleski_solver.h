/******************************************************************************/
/* Name: choleski_solver.h                                                    */
/* Description: functions for assignment 1 of ECSE 543 Numerical Methods      */
/* Date: 2020/09/10                                                           */
/* Author: Raymond Yang                                                       */
/******************************************************************************/

#ifndef __CHOLESKI_SOLVER__
#define __CHOLESKI_SOLVER__

#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Matrix.h"

extern int FLAG;

/*  
    Choleski Decomposition

    Formula: A = L*L'; solves for L
    A: square, symmetric, and positive definite matrix
    *L: decomposed lower triangular matrix
*/
template <typename T> 
void choleski(Matrix<T> A, Matrix<T>* L){
    int row = A.get_n_row();
    int col = A.get_n_col();
    for(int j = 0; j < row; j++){
        if(A.get(j, j) <= 0){
            FLAG -= 1;
            throw "choleski Error: A is not P.D.";
        }
        L->set(j, j, sqrt(A.get(j, j)));
        for(int i = j + 1; i < row; i++){
            L->set(i, j, A.get(i, j)/L->get(j, j));
            for(int k = j + 1; k <= i; k++){
                A.set(i, k, A.get(i, k) - L->get(i, j)*L->get(k, j));
            }
        }
    }
    return;
}

/*  In-place computation version of choleski function  */
template <typename T> 
void choleski(Matrix<T>* A){
    int row = A->get_n_row();
    int col = A->get_n_col();
    for(int j = 0; j < row; j++){
        if(A->get(j, j) <= 0){
            FLAG -= 1;
            throw "choleski Error: A is not P.D.";
        }
        A->set(j, j, sqrt(A->get(j, j)));
        for(int i = j + 1; i < row; i++){
            A->set(i, j, A->get(i, j)/A->get(j, j));
            for(int k = j + 1; k <= i; k++){
                A->set(i, k, A->get(i, k) - A->get(i, j)*A->get(k, j));
            }
        }
    }
    // set upper right (excluding diagonal) to 0
    for(int j = 1; j < row; j++){
        for(int i = 0; i < j; i++){
            A->set(i, j, 0);
        }
    }
    return;
}

/*  
    Forward Elimination

    Formula: L*y = b; solves for y
    L: lower triangular matrix
    b: output vector
    *y: unknown vector
*/
template <typename T>
void forward_elimination(Matrix<T>& L, Matrix<T> b, Matrix<T>* y){
    int length = b.get_n_row();
    for(int j = 0; j < length; j++){
        y->set(j, b.get(j)/L.get(j, j));
        for(int i = j + 1; i < length; i++){
            b.set(i, b.get(i) - L.get(i, j)*y->get(j));
        }
    }
    return;
}

/*  In-place computation version of forward_elimination  */
template <typename T>
void forward_elimination(Matrix<T>& L, Matrix<T>* b){
    int length = b->get_n_row();
    for(int j = 0; j < length; j++){
        b->set(j, b->get(j)/L.get(j, j));
        for(int i = j + 1; i < length; i++){
            b->set(i, b->get(i) - L.get(i, j)*b->get(j));
        }
    }
    return;
}

/*  
    Elimination 

    Cholesky decomposition and forward elimination combined

    Formula: A*x = b -> L*y = b; solves for L and y
    A: square, symmetric, P.D. 
    b: output vector
    *L: lower triangular matrix
    *y: unknown vector
*/
template <typename T>
void elimination(Matrix<T> A, Matrix<T> b, Matrix<T>* L, Matrix<T>* y){

    int row = A.get_n_row();
    int col = A.get_n_col();
    for(int j = 0; j < row; j++){
        if(A.get(j, j) <= 0){
            FLAG -= 1;
            throw "choleski Error: A is not P.D.";
        }
        L->set(j, j, sqrt(A.get(j, j)));
        y->set(j, y->get(j)/L->get(j, j));
        for(int i = j + 1; i < row; i++){
            L->set(i, j, A.get(i, j)/L->get(j, j));
            b.set(i, b.get(i) - L->get(i, j) * y->get(j));
            for(int k = j + 1; k <= i; k++){
                A.set(i, k, A.get(i, k) - L->get(i, j)*L->get(k, j));
            }
        }
    }
    return;
}

/*  In-place computation version of elimination  */
template <typename T>
void elimination(Matrix<T>* A, Matrix<T>* b){
    int row = A->get_n_row();
    int col = A->get_n_col();
    for(int j = 0; j < row; j++){
        if(A->get(j, j) <= 0){
            FLAG -= 1;
            throw "choleski Error: A is not P.D.";
        }
        A->set(j, j, sqrt(A->get(j, j)));
        b->set(j, b->get(j)/A->get(j, j));
        for(int i = j + 1; i < row; i++){
            A->set(i, j, A->get(i, j)/A->get(j, j));
            b->set(i, b->get(i) - A->get(i, j) * b->get(j));
            for(int k = j + 1; k <= i; k++){
                A->set(i, k, A->get(i, k) - A->get(i, j)*A->get(k, j));
            }
        }
    }
        // set upper right (excluding diagonal) to 0
    for(int j = 1; j < row; j++){
        for(int i = 0; i < j; i++){
            A->set(i, j, 0);
        }
    }
    return;
}

/*  
    Back Substitution 

    Formula: Ux = y; solves for x
    U: Upper triangular matrix
    y: output vector
    *x: unknown vector
*/
template <typename T>
void back_substitution(Matrix<T>& U, Matrix<T> y, Matrix<T>* x){
    int length = y.get_n_row();
    for(int j = length - 1; j >= 0; j--){
        x->set(j, y.get(j)/U.get(j, j));
        for(int i = j - 1; i >= 0; i--){
            y.set(i, y.get(i) - U.get(i, j)*x->get(j));
        }
    }
    return;
}

/*  In-place computation version of back_substitution  */
template <typename T>
void back_substitution(Matrix<T>& U, Matrix<T>* y){
    int length = y->get_n_row();
    for(int j = length - 1; j >= 0; j--){
        y->set(j, y->get(j)/U.get(j, j));
        for(int i = j - 1; i >= 0; i--){
            y->set(i, y->get(i) - U.get(i, j)*y->get(j));
        }
    }
    return;
}

/*  
    Solve Ax = b using Choleski Decomposition 

    A: square, symmetric, P.D.
    b: output vector
    *x: unknown vector
*/
template <typename T>
void solve_choleski(Matrix<T> A, Matrix<T> b, Matrix<T>* x){
    elimination(&A, &b);
    A.transpose();
    back_substitution(A, &b);
    *x = b;
    return;
}

/*  In-place computation version of solve_choleski  */
template <typename T>
void solve_choleski(Matrix<T>* A, Matrix<T>* b){
    elimination(A, b);
    A->transpose();
    back_substitution(*A, b);
    return;
}

/*  
    Create square symmetric P.D. matrix 
    
    Implemented in Matrix.h

    n: matrix dimension (n, n)
    Returns: matrix
*/

#endif


