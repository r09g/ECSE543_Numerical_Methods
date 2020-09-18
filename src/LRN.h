/******************************************************************************/
/* Name: LRN.h                                                              */
/* Description: Linear Resistive Network                                      */
/* Date: 2020/09/10                                                           */
/* Author: Raymond Yang                                                       */
/******************************************************************************/

#ifndef __LRN__
#define __LRN__

#include "Matrix.h"
#include "matrix_solver.h" 

template <typename T>
class LRN {
    private:
        Matrix<T> A_;
        Matrix<T> J_;
        Matrix<T> R_;
        Matrix<T> E_;
        Matrix<T> v_;
        Matrix<T> i_;

    public:
        LRN(const std::string& filepath);
        ~LRN();

        // Getters
        Matrix<T> get_A() const;
        Matrix<T> get_J() const;
        Matrix<T> get_R() const;
        Matrix<T> get_E() const;
        Matrix<T> get_v() const;
        Matrix<T> get_i() const;
        
        // Circuit Solvers
        void choleski_solve();

};

// -----------------------------------------------------------------------------
// Constructor & Destructor
// -----------------------------------------------------------------------------

/*  
    Create a LRN from .csv description  
    
    .csv format:
    #nodes, #branches
    N+, N-, J, R, E
    ...

    Convention:
               <--- i_k
        |------------0     
    |-------|        
    ^       R        +
    ^       |
    I       -       v_k
    ^       V
    ^       +        -
    |-------|        
        |------------0

    KVL: (I + i_k)*R - V = v_k   

    GND node number must be 0
    I flowing out is +; I flowing in is -

    A, J, R, E: Pass by reference these matrices to receive output
    filepath: string filepath
    &A: reduced incidence matrix
    &J: I source vector
    &R: R vector
    &E: B source vector    
*/
template <typename T>
LRN<T>::LRN(const std::string& filepath){
    std::ifstream file;
    file.open(filepath);

    int n_nodes, n_branches;
    std::string line, output;
    getline(file, line);
    std::stringstream ss(line);
    getline(ss, output, ',');
    n_nodes = stoi(output);
    getline(ss, output, ',');
    n_branches = stoi(output);
    if(n_nodes <  2 || n_branches < 2){
        throw "LRN::LRN ERROR: Too few nodes/branches"; 
    }  
    
    this->A_ = Matrix<T>(n_nodes, n_branches);
    this->J_ = Matrix<T>(n_branches, 1);
    this->R_ = Matrix<T>(n_branches, 1);
    this->E_ = Matrix<T>(n_branches, 1);
    
    for(int i = 0; i < n_branches; i++){
        getline(file,line);
        std::stringstream ss(line);
        getline(ss, output, ',');
        this->A_.set(stoi(output), i, 1);
        getline(ss, output, ',');
        this->A_.set(stoi(output), i, -1);
        getline(ss, output, ',');
        this->J_.set(i, stod(output));
        getline(ss, output, ',');
        this->R_.set(i, stod(output));
        getline(ss, output, ',');
        this->E_.set(i, stod(output));
    }
    this->A_ = this->A_.get(1, n_nodes - 1, 0, n_branches - 1);

    file.close();
    return;
}

/*  Default destructor used  */
template <typename T>
LRN<T>::~LRN() {}

// -----------------------------------------------------------------------------
// Member Functions
// -----------------------------------------------------------------------------

/*  Solve LRN using Choleski Decomposition  */
template <typename T>
void LRN<T>::choleski_solve(){
    Matrix<T> Y(this->A_.get_n_col(), this->A_.get_n_col());
    for(int i = 0; i < this->A_.get_n_col(); i++){
        Y.set(i, i, 1/this->R_.get(i));
    }

    Matrix<T> A = this->A_*Y*transpose(this->A_);
    Matrix<T> b = this->A_*(this->J_ - Y*this->E_);
    solve_choleski(&A, &b);
    
    this->v_ = b;
    this->i_ = Y*(transpose(this->A_)*this->v_) + Y*this->E_ - this->J_;
    return;
}

/*  Get the reduced incidence matrix  */
template <typename T>
Matrix<T> LRN<T>::get_A() const{
    return this->A_;
}

/*  Get the current source vector  */
template <typename T>
Matrix<T> LRN<T>::get_J() const{
    return this->J_;
}

/*  Get the branch resistance vector  */
template <typename T>
Matrix<T> LRN<T>::get_R() const{
    return this->R_;
}

/*  Get the voltage source vector  */
template <typename T>
Matrix<T> LRN<T>::get_E() const{
    return this->E_;
}

/*  Get the node voltages vector  */
template <typename T>
Matrix<T> LRN<T>::get_v() const{
    return this->v_;
}

/*  Get the branch current vector  */
template <typename T>
Matrix<T> LRN<T>::get_i() const{
    return this->i_;
}


#endif








