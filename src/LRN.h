/******************************************************************************/
/* Name: LRN.h                                                                */
/* Description: Linear Resistive Network                                      */
/* Date: 2020/09/17                                                           */
/* Author: Raymond Yang                                                       */
/******************************************************************************/

#ifndef __LRN__
#define __LRN__

#include <math.h>
#include "Matrix.h"
#include "matrix_solver.h" 

extern int FLAG;

template <typename T>
class LRN {
    private:
        Matrix<T> A_;   // reduced incidence matrix
        Matrix<T> J_;   // current source vector
        Matrix<T> R_;   // resistance vector
        Matrix<T> E_;   // voltage source vector
        Matrix<T> v_;   // node voltages
        Matrix<T> i_;   // branch currents
        Matrix<T> Y_;   // diagonal conductance matrix

        // helper function
        void compute_Y();

    public:
        LRN();
        LRN(const std::string& filepath);
        LRN(int n_nodes, int n_branches, T* data);
        ~LRN();

        // Getters
        Matrix<T> get_A() const;
        Matrix<T> get_J() const;
        Matrix<T> get_R() const;
        Matrix<T> get_E() const;
        Matrix<T> get_v() const;
        Matrix<T> get_i() const;
        Matrix<T> get_Y() const;
        
        // Generate NxN Resistor Mesh
        void nxn_res_mesh(int N, bool write_csv=false, double r_ohms=10e3);

        // Circuit Solvers
        void cholesky_solve();
        void cholesky_solve_banded(int HBW=-1);

};

// -----------------------------------------------------------------------------
// Constructor & Destructor
// -----------------------------------------------------------------------------

/*  Default constructor  */
template <typename T>
LRN<T>::LRN(){}

/*  
    Create a LRN from .csv description  
    
    A, J, R, E: Pass by reference these matrices to receive output
    filepath: string filepath
    &A: reduced incidence matrix
    &J: I source vector
    &R: R vector
    &E: B source vector

    .csv format:
    Ground node is numbered 0.
    -----------------
    #nodes, #branches
    N+, N-, J, R, E
    ...
    -----------------

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
        FLAG -= 1;
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
    this->compute_Y();

    file.close();
    return;
}

/*  Specify circuit in array  */
template <typename T>
LRN<T>::LRN(int n_nodes, int n_branches, T* data){
    this->A_ = Matrix<T>(n_nodes, n_branches);
    this->J_ = Matrix<T>(n_branches, 1);
    this->R_ = Matrix<T>(n_branches, 1);
    this->E_ = Matrix<T>(n_branches, 1);

    for(int i = 0; i < n_branches; i++){
        this->A_.set(data[5*i + 0], i, 1);
        this->A_.set(data[5*i + 1], i, -1);
        this->J_.set(i, data[5*i + 2]);
        this->R_.set(i, data[5*i + 3]);
        this->E_.set(i, data[5*i + 4]);
    }

    this->A_ = this->A_.get(1, n_nodes - 1, 0, n_branches - 1);
    this->compute_Y();
    return;
}

/*  Default destructor used  */
template <typename T>
LRN<T>::~LRN() {}

// -----------------------------------------------------------------------------
// Member Functions
// -----------------------------------------------------------------------------

/*  Compute the diagonal conductance matrix  */
template <typename T>
void LRN<T>::compute_Y(){
    this->Y_ = Matrix<T>(this->A_.get_n_col(), this->A_.get_n_col());
    for(int i = 0; i < this->A_.get_n_col(); i++){
        this->Y_.set(i, i, 1/this->R_.get(i));
    }
    return;
}

/*  
    Generate NxN resistor mesh.

    Resistors have uniform resistance r_ohms = 10k by default.
    Test source: 1A current source with Norton equivlent resistance = r_ohms.
*/
template <typename T>
void LRN<T>::nxn_res_mesh(int N, bool write_csv, double r_ohms){
    if(N <= 0){
        throw "LRN::nxn_res_mesh ERROR: Invalid mesh size";
    }
    int n_nodes = (N + 1)*(N + 1);
    int n_branches = 2*N*(N + 1) + 1;
    
    Matrix<T> cct_file(n_branches + 1, 5);
    cct_file.set(0, 0, n_nodes);
    cct_file.set(0, 1, n_branches);

    // iterate through branches
    for(int k = 0; k < n_branches; k++){
        if(k == 0){  // test source branch
            cct_file.set(k + 1, 0, n_nodes - 1);
            cct_file.set(k + 1, 1, 0);
            cct_file.set(k + 1, 2, 1);
            cct_file.set(k + 1, 3, r_ohms);
            cct_file.set(k + 1, 4, 0);
        } else if(k <= N*(N + 1)) {  // horizontal branches
            cct_file.set(k + 1, 0, k + ceil((double)(k)/N) - 1);
            cct_file.set(k + 1, 1, k + ceil((double)(k)/N) - 2);
            cct_file.set(k + 1, 2, 0);
            cct_file.set(k + 1, 3, r_ohms);
            cct_file.set(k + 1, 4, 0);
        } else {  // vertical branches
            cct_file.set(k + 1, 0, k - N*(N + 1) - 1);
            cct_file.set(k + 1, 1, k - N*N);
            cct_file.set(k + 1, 2, 0);
            cct_file.set(k + 1, 3, r_ohms);
            cct_file.set(k + 1, 4, 0);
        }
    }

    if(write_csv == true){
        std::string filepath = std::string("./") + std::to_string(N) 
            + std::string("x") + std::to_string(N) 
            + std::string("_R_Mesh_CCT.csv");
        cct_file.write_mat(filepath);
    }
    
    cct_file = cct_file.get(1, cct_file.get_n_row() - 1, 
        0, cct_file.get_n_col() - 1);
    LRN<T> mesh_cct(n_nodes, n_branches, cct_file.get());
    *this = mesh_cct;
    return;
}

/*  Solve LRN using Cholesky Decomposition  */
template <typename T>
void LRN<T>::cholesky_solve(){
    Matrix<T> A = this->A_*this->Y_*transpose(this->A_);
    Matrix<T> b = this->A_*(this->J_ - this->Y_*this->E_);
    Matrix_Solver::cholesky_solve(&A, &b);
    
    this->v_ = b;
    this->i_ = this->Y_*(transpose(this->A_)*this->v_) 
        + this->Y_*this->E_ - this->J_;
    return;
}

/*  Solve LRN using Cholesky Decomposition  */
template <typename T>
void LRN<T>::cholesky_solve_banded(int HBW){
    Matrix<T> A = this->A_*this->Y_*transpose(this->A_);
    Matrix<T> b = this->A_*(this->J_ - this->Y_*this->E_);
    Matrix_Solver::cholesky_solve_banded(&A, &b, HBW);
    
    this->v_ = b;
    this->i_ = this->Y_*(transpose(this->A_)*this->v_) 
        + this->Y_*this->E_ - this->J_;
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

/*  Get the diagonal conductance matrix  */
template <typename T>
Matrix<T> LRN<T>::get_Y() const{
    return this->Y_;
}

#endif








