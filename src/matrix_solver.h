/******************************************************************************/
/* Name: Matrix_Solver.h                                                      */
/* Description: functions for assignment 1 of ECSE 543 Numerical Methods      */
/* Date: 2020/09/10                                                           */
/* Author: Raymond Yang                                                       */
/******************************************************************************/

#ifndef __MATRIX_SOLVER__
#define __MATRIX_SOLVER__

#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Matrix.h"

namespace Matrix_Solver{

    template <typename T> void cholesky(Matrix<T> A, Matrix<T>* L);
    template <typename T> void cholesky(Matrix<T>* A);
    template <typename T> 
        void forward_elimination(Matrix<T>& L, Matrix<T> b, Matrix<T>* y);
    template <typename T> void forward_elimination(Matrix<T>& L, Matrix<T>* b);
    template <typename T>
        void elimination(Matrix<T> A, Matrix<T> b, Matrix<T>* L, Matrix<T>* y);
    template <typename T> void elimination(Matrix<T>* A, Matrix<T>* b);
    template <typename T> 
        void back_substitution(Matrix<T>& U, Matrix<T> y, Matrix<T>* x);
    template <typename T> void back_substitution(Matrix<T>& U, Matrix<T>* y);
    template <typename T>  
        void cholesky_solve(Matrix<T> A, Matrix<T> b, Matrix<T>* x);
    template <typename T> void cholesky_solve(Matrix<T>* A, Matrix<T>* b);
    template <typename T> int find_HBW(Matrix<T>* A);
    template <typename T> void cholesky_banded(Matrix<T>* A, int HBW=-1);
    template <typename T> 
        void elimination_banded(Matrix<T>* A, Matrix<T>* b, int HBW=-1);
    template <typename T>
        void cholesky_solve_banded(Matrix<T>* A, Matrix<T>* b, int HBW=-1);

    template <typename T> Matrix<T> CG_solve(Matrix<T> A, Matrix<T> b, 
        int itr = -1, Matrix<T> IC = Matrix<T>(0));

}

/*  
    Cholesky Decomposition

    Formula: A = L*L'; solves for L
    A: square, symmetric, and positive definite matrix
    *L: decomposed lower triangular matrix
*/
template <typename T> 
void Matrix_Solver::cholesky(Matrix<T> A, Matrix<T>* L){
    int row = A.get_n_row();
    int col = A.get_n_col();
    for(int j = 0; j < row; j++){
        if(A.get(j, j) <= 0){
            throw "cholesky Error: A is not P.D.";
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

/*  In-place computation version of cholesky function  */
template <typename T> 
void Matrix_Solver::cholesky(Matrix<T>* A){
    int row = A->get_n_row();
    int col = A->get_n_col();
    for(int j = 0; j < row; j++){
        if(A->get(j, j) <= 0){
            throw "cholesky Error: A is not P.D.";
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
void Matrix_Solver::forward_elimination(Matrix<T>& L, Matrix<T> b, Matrix<T>* y){
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
void Matrix_Solver::forward_elimination(Matrix<T>& L, Matrix<T>* b){
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
void Matrix_Solver::elimination(Matrix<T> A, Matrix<T> b, Matrix<T>* L, Matrix<T>* y){

    int row = A.get_n_row();
    int col = A.get_n_col();
    for(int j = 0; j < row; j++){
        if(A.get(j, j) <= 0){
            throw "cholesky Error: A is not P.D.";
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
void Matrix_Solver::elimination(Matrix<T>* A, Matrix<T>* b){
    int row = A->get_n_row();
    int col = A->get_n_col();
    for(int j = 0; j < row; j++){
        if(A->get(j, j) <= 0){
            throw "cholesky Error: A is not P.D.";
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
void Matrix_Solver::back_substitution(Matrix<T>& U, Matrix<T> y, Matrix<T>* x){
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
void Matrix_Solver::back_substitution(Matrix<T>& U, Matrix<T>* y){
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
    Solve Ax = b using Cholesky Decomposition 

    A: square, symmetric, P.D.
    b: output vector
    *x: unknown vector
*/
template <typename T>
void Matrix_Solver::cholesky_solve(Matrix<T> A, Matrix<T> b, Matrix<T>* x){
    elimination(&A, &b);
    A.transpose();
    back_substitution(A, &b);
    *x = b;
    return;
}

/*  In-place computation version of cholesky_solve  */
template <typename T>
void Matrix_Solver::cholesky_solve(Matrix<T>* A, Matrix<T>* b){
    elimination(A, b);
    A->transpose();
    back_substitution(*A, b);
    return;
}

/*  Find the half-bandwidth of A  */
template <typename T>
int Matrix_Solver::find_HBW(Matrix<T>* A){
    int num_row = A->get_n_row();
    int HBW = 0;
    for(int i = 0; i < num_row; i++){
        int j = num_row - 1;
        while(j >= i){
            if(A->get(i, j) == 0){
                j--;
            } else {
                break;
            }
        }
        HBW = std::max(HBW, j - i + 1);
    }
    return HBW;
}

/*    
    Banded Cholesky Decomposition (In-Place)

    Formula: A = L*L'; solves for L
    *A: square, symmetric, and positive definite matrix
    HBW: half bandwidth. Will automatically determine if not provided.
*/
template <typename T> 
void Matrix_Solver::cholesky_banded(Matrix<T>* A, int HBW){
    if(HBW == -1){
        HBW = find_HBW(A);
    }
    int row = A->get_n_row();
    int col = A->get_n_col();
    for(int j = 0; j < row; j++){
        if(A->get(j, j) <= 0){
            throw "cholesky Error: A is not P.D.";
        }
        A->set(j, j, sqrt(A->get(j, j)));
        for(int i = j + 1; i < std::min(row, j + HBW); i++){
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
    Elimination (In-Place)

    Banded Cholesky decomposition and forward elimination combined

    Formula: A*x = b -> L*y = b; solves for L and y
    *A: square, symmetric, P.D. 
    *b: output vector
    HBW: half bandwidth. Will automatically determine if not provided.
*/
template <typename T>
void Matrix_Solver::elimination_banded(Matrix<T>* A, Matrix<T>* b, int HBW){
    if(HBW == -1){
        HBW = find_HBW(A);
    }
    int row = A->get_n_row();
    int col = A->get_n_col();
    for(int j = 0; j < row; j++){
        if(A->get(j, j) <= 0){
            throw "cholesky Error: A is not P.D.";
        }
        A->set(j, j, sqrt(A->get(j, j)));
        b->set(j, b->get(j)/A->get(j, j));
        for(int i = j + 1; i < std::min(row, j + HBW); i++){
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
    Solve Ax = b using Banded Cholesky Decomposition (In-Place)

    *A: square, symmetric, P.D.
    *b: output vector
    HBW: half bandwidth. Will automatically determine if not provided.
*/
template <typename T>
void Matrix_Solver::cholesky_solve_banded(Matrix<T>* A, Matrix<T>* b, int HBW){
    elimination_banded(A, b, HBW);
    A->transpose();
    back_substitution(*A, b);
    return;
}

/*  
    Solve Ax = b using Conjugate Gradient (Unpreconditioned)

    A: square, symmetric, P.D.
    b: output vector
    tol: tolerance for stopping condition
    
    Returns: solution vector
*/
template <typename T> 
Matrix<T> Matrix_Solver::CG_solve(Matrix<T> A, Matrix<T> b, int itr, 
    Matrix<T> IC){
    Matrix<T> x(b.get_n_row(), 1);
    Matrix<T> r(b.get_n_row(), 1);
    Matrix<T> p(b.get_n_row(), 1);
    double alpha;
    double beta;
    // initial condition
    if(IC.get_n_row() != b.get_n_row() || IC.get_n_col() != 1){
        x = Matrix<T>(b.get_n_row(), 1);
    } else {
        x = IC;
    }
    itr = (itr < 0) ? b.get_n_row() : itr;
    Matrix<> r_evo(b.get_n_row(), itr + 1);

    // initial guess
    r = b - A*x;
    p = r;
    // iteration
    int count = 0;
    while(count <= itr){
        r_evo.set(0, count, r);
        alpha = (transpose(p)*r).get(0)/(transpose(p)*A*p).get(0);
        x = x + alpha*p;
        r = b - A*x;
        beta = -1*(transpose(p)*A*r).get(0)/(transpose(p)*A*p).get(0);
        p = r + beta*p;
        count++;
    }
    
    r_evo.write_mat("./data/A2/r_results.csv");
    return x;
}



#endif
