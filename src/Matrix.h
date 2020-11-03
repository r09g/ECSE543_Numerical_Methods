/******************************************************************************/
/* Name: Matrix.h                                                             */
/* Date: 2020/09/10                                                           */
/* Author: Raymond Yang                                                       */
/******************************************************************************/

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <new>
#include <cstdlib>
#include <type_traits>
#include "Basic.h"

template <class T = double>
class Matrix {
    private:
        // parameters
        int n_row;
        int n_col;
        T* data = NULL;

        // memory allocation
        T* allocate(int n_row, int n_col);
        void deallocate();

        // deep copy
        Matrix<T> deep_copy() const;

    public:
        // constructors
        Matrix();
        Matrix(int n);
        Matrix(int n_row, int n_col);
        Matrix(int n_row, int n_col, T val);
        Matrix(int n_row, int n_col, const T* data);
        Matrix(const Matrix<T>& mat);

        // show entire matrix
        void show() const;
        
        // row, col to index conversion
        int get_index(int row, int col) const;

        // create special matrices
        static Matrix<T> zero_mat(int n);
        static Matrix<T> zero_mat(int n_row, int n_col);
        static Matrix<T> identity_mat(int n);
        static Matrix<T> rand_mat(int n);
        static Matrix<T> rand_mat(int n_row, int n_col);
        static Matrix<T> SSPD_mat(int n);

        // getters and setters
        T* get() const;
        void set(const T* data, unsigned int length);
        void set(T*&& data);
        T get(int row, int col) const;
        T& get_ref(int row, int col);
        void set(int row, int col, T value);
        T get(int index);
        void set(int index, T value);
        Matrix<T> get(int row1, int row2, int col1, int col2) const;
        void set(int row1, int row2, int col1, int col2, T value);
        void set(int row, int col, const Matrix<T> mat);
        Matrix<T> get_row(int row) const;
        Matrix<T> get_col(int col) const;
        int get_n_row() const;
        int get_n_col() const;

        // copy (duplicate) matrix
        Matrix<T>& operator= (const Matrix<T>& mat);
        Matrix<T>& operator= (const T* data);
        static T* deep_copy(const T* data, unsigned int length);

        // matrix operations
        void transpose();
        Matrix<T> mul(const Matrix<T>& mat);
        Matrix<T> div(const Matrix<T>& mat);
        Matrix<T> operator+ (const Matrix<T>& mat);
        Matrix<T> operator- (const Matrix<T>& mat);
        Matrix<T> operator* (const Matrix<T>& mat);
        Matrix<T> mul(const T value);
        Matrix<T> div(const T value);
        T sum();
        T mean();

        // check matrix characteristics
        bool is_square();
        bool is_symmetric();
        bool operator== (const Matrix& mat);

        // matrix data type conversions
        Matrix<int> to_int();
        Matrix<double> to_double();
        Matrix<float> to_float();

        // I/O to Excel
        static Matrix<T> read_mat(const std::string& filepath);
        void write_mat(const std::string& filepath);

        // destructor
        ~Matrix();
};

// -----------------------------------------------------------------------------
// Dummy Namespace
// -----------------------------------------------------------------------------

/*  Used for calling non-member functions over member functions  */
namespace Matrix_Dummy{

    /*  creates a matrix containing the transpose of mat  */
    template <class T>
    Matrix<T> transpose(Matrix<T> mat){
        mat.transpose();
        return mat;
    }

}

// -----------------------------------------------------------------------------
// Static Member-Functions
// -----------------------------------------------------------------------------

/*  creates a deep copy of the array  */
template <class T>
T* Matrix<T>::deep_copy(const T* data, unsigned int length){
    T* new_data = new T[length]();
    for(int i = 0; i < length; i++){
        new_data[i] = data[i];
    }
    return new_data;
}

/*  reads and returns the matrix from a csv  */
template <class T>
Matrix<T> Matrix<T>::read_mat(const std::string& filepath){
    
    std::ifstream file;
    file.open(filepath);
    
    int n_row, n_col;
    std::string line, output;
    getline(file, line);
    std::stringstream ss(line);
    getline(ss, output, ',');
    n_row = stoi(output);
    getline(ss, output, ',');
    n_col = stoi(output);  
    
    Matrix<T> new_matrix(n_row, n_col);
    int row = 0;
    int col = 0;
    while(getline(file,line)){
        std::stringstream ss(line);
        col = 0;
        while(getline(ss, output, ',')){
            new_matrix.set(row, col, stod(output));
            col++;
            if(col >= n_col){
                break;
            }
        }
        row++;
        if(row >= n_row){
            break;
        }
    }
    file.close();
    return new_matrix;
}

/*  creates a square zero matrix of size n  */
template <class T>
Matrix<T> Matrix<T>::zero_mat(int n){
    return Matrix<T>(n);
}

/*  creates a square zero matrix of shape (n_row, n_col)  */
template <class T>
Matrix<T> Matrix<T>::zero_mat(int n_row, int n_col){
    return Matrix<T>(n_row, n_col);
}

/*  creates an identity matrix of size n  */
template <class T>
Matrix<T> Matrix<T>::identity_mat(int n){
    Matrix<T> mat(n);
    for(int i = 0; i < n; i++){
        mat.set(i,i,1);
    }
    return mat;
}

/*  
    creates a random square matrix of shape (n, n)
    values range in [-n^2/2,n^2/2] 
*/
template <class T>
Matrix<T> Matrix<T>::rand_mat(int n){
    return rand_mat(n, n);
}

/*  
    creates a random matrix of shape (n_row, n_col)
    values range in [-n^2/2,n^2/2]
*/
template <class T>
Matrix<T> Matrix<T>::rand_mat(int n_row, int n_col){
    Matrix<T> mat(n_row, n_col);
    for(int i = 0; i < n_row*n_col; i++){
        mat.set(i, rand() % (n_row*n_col + 1) - n_row*n_col/2);
    }
    return mat;
}

/*  
    creates an square, symmetric, and positive definite matrix of shape (n, n)
*/
template <class T>
Matrix<T> Matrix<T>::SSPD_mat(int n){
    Matrix<T> mat(n, n);
    for(int row = 0; row < n; row++){
        for(int col = 0; col <= row; col++){
            mat.set(row, col, rand() % (n*n + 1) - n*n/2);
        }
    }
    for(int diag = 0; diag < n; diag++){
        mat.set(diag, diag, std::abs(mat.get(diag, diag)) + 1);
    }  
    mat = mat * Matrix_Dummy::transpose(mat);
    return mat;
}

// -----------------------------------------------------------------------------
// Constructors  & Destructor
// -----------------------------------------------------------------------------

/*  creates a square matrix of zeros with size n  */
template <class T>
Matrix<T>::Matrix(int n): n_row(n), n_col(n) {
    this->data = allocate(n, n);
    return;
}

/*  creates an empty matrix  */
template <class T>
Matrix<T>::Matrix(): n_row(0), n_col(0) {
    this->data = allocate(n_row, n_col);
    return;
}

/*  creates a (n_row, n_col) matrix with elements initialized to 0  */
template <class T>
Matrix<T>::Matrix(int n_row, int n_col): n_row(n_row), n_col(n_col) {
    this->data = allocate(n_row, n_col);
    return;
}

/*  creates a (n_row, n_col) matrix with elements initialized to val  */
template <class T>
Matrix<T>::Matrix(int n_row, int n_col, T val): n_row(n_row), n_col(n_col) {
    this->data = allocate(n_row, n_col);
    for(int i = 0; i < n_row*n_col; i++){
        this->data[i] = val;
    }
    return;
}

/*  creates a (n_row, n_col) matrix initialized to "data"
    values copied over  */
template <class T>
Matrix<T>::Matrix(int n_row, int n_col, const T* data){
    this->n_row = n_row;
    this->n_col = n_col;
    this->data = allocate(this->n_row, this->n_col);
    for(int i = 0; i < n_row*n_col; i++){
        this->data[i] = data[i];
    }
    return;
}

/*  copy constructor; returns deep copy of "mat"  */
template <class T>
Matrix<T>::Matrix(const Matrix<T>& mat): n_row(mat.n_row), n_col(mat.n_col){
    this->data = deep_copy(mat.data, this->n_row * this->n_col); 
    return;
}

/*  destructor; calls deallocate()  */
template <class T>
Matrix<T>::~Matrix(){
    this->deallocate();
    // std::cout << "deleted" << std::endl;
    return;
}

// -----------------------------------------------------------------------------
// Member-Functions
// -----------------------------------------------------------------------------

/*  allocates memory for data array with n_row*n_col elements
    returns a pointer  */
template <class T>
T* Matrix<T>::allocate(int n_row, int n_col){
    T* new_data = NULL;
    if(this->n_row != 0 && this->n_col != 0){
        try{
            new_data = new T[this->n_row * this->n_col]();
        } catch(std::bad_alloc& e){
            std::cout << e.what() << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    // std::cout << "created" << std::endl;
    return new_data;
}

/*  deallocates memory assigned to data array  */
template <class T>
void Matrix<T>::deallocate(){
    delete [] this->data;
    this->data = NULL;
    // std::cout << "deallocated" << std::endl;
    return;
}

/*  prints a matrix in terminal  */
template <class T>
void Matrix<T>::show() const {
    if(this->n_row == 0 && this->n_col == 0){
        std::cout << "0 [ ]" << std::endl;
    }
    for(int i = 0; i < this->n_row; i++){
        std::cout << i << " [ ";
        for(int j = 0; j < this->n_col; j++){
            std::cout << this->get(i,j) << ",";
        }
        std::cout << " ]" << std::endl;
    }   
    return;
}

/*  get a pointer to a deep copy of the data array  */
template <class T>
T* Matrix<T>::get() const {
    return deep_copy(this->data, this->n_col*this->n_row);
}

/*  set a new data array to the Matrix by deep copy  */
template <class T>
void Matrix<T>::set(const T* data, unsigned int length){
    if(this->data != NULL){
        this->deallocate();
    }
    this->data = deep_copy(data, length);
    return;
}

/*  set a new data array to the Matrix by move-semantics  */
template <class T>
void Matrix<T>::set(T*&& data){
    if(this->data != NULL){
        this->deallocate();
    }
    this->data = data;
    data = NULL;
    return;
}

/*  get the (row, col) element from the Matrix  */
template <class T>
T Matrix<T>::get(int row, int col) const {
    return this->data[this->get_index(row,col)];
}

/*  get the reference to (row, col) element from the Matrix  */
template <class T>
T& Matrix<T>::get_ref(int row, int col){
    return this->data[this->get_index(row,col)];
}

/*  set the (row, col) element to "value"  */
template <class T>
void Matrix<T>::set(int row, int col, T value){
    this->data[this->get_index(row,col)] = value;
    return;
}

/*  get the i-th element in the data array  */
template <class T>
T Matrix<T>::get(int index){
    return this->data[index];
}

/*  set the i-th element in the data array  */
template <class T>
void Matrix<T>::set(int index, T value){
    this->data[index] = value;
    return;
}

/*  get a copy of part of the matrix (row1:row2, col1:col2) */
template <class T>
Matrix<T> Matrix<T>::get(int row1, int row2, int col1, int col2) const {
    Matrix<T> new_mat(row2 - row1 + 1, col2 - col1 + 1);
    for(int i = 0; i < new_mat.n_row; i++){
        for(int j = 0; j < new_mat.n_col; j++){
            new_mat.data[new_mat.get_index(i,j)] 
                = this->get(row1 + i, col1 + j);
        }
    }
    return new_mat;
}

/*  set the same value to part of the Matrix (row1:row2, col1:col2)  */
template <class T>
void Matrix<T>::set(int row1, int row2, int col1, int col2, T value){
    for(int row = row1; row <= row2; row++){
        for(int col = col1; col <= col2; col++){
            this->set(row, col, value);
        }
    }
    return;
}

/*  set part of the Matrix to "mat", starting at (row, col)  */
template <class T>
void Matrix<T>::set(int row, int col, const Matrix<T> mat){
    for(int i = 0; i < mat.n_row; i++){
        for(int j = 0; j < mat.n_col; j++){
            this->data[this->get_index(row + i, col + j)]
                = mat.data[i*mat.n_col + j];
        }
    }
    return;
}

/*  get the row at (row,:)  */
template <class T>
Matrix<T> Matrix<T>::get_row(int row) const {
    Matrix<T> new_mat(1, this->n.col);
    for(int i = 0; i < this->n.col; i++){
        new_mat.data[i] = this->get(row, i);
    }
    return new_mat;
}

/*  get the column at (:, col)  */
template <class T>
Matrix<T> Matrix<T>::get_col(int col) const {
    Matrix<T> new_mat(this->n.col, 1);
    for(int i = 0; i < this->n.row; i++){
        new_mat.data[i] = this->get(i, col);
    }
    return new_mat;
}

/*  convert (row, col) to index for data array  */
template <class T>
int Matrix<T>::get_index(int row, int col) const{
    return row*this->n_col + col;
}

/*  return the total number of rows  */
template <class T>
int Matrix<T>::get_n_row() const{
    return this->n_row;
}

/*  return the total number of columns  */
template <class T>
int Matrix<T>::get_n_col() const{
    return this->n_col;
}

/*  deallocate current data array in Matrix and assign data array with values
    from "mat"; returns self  */
template <class T>
Matrix<T>& Matrix<T>::operator= (const Matrix<T>& mat){
    this->n_col = mat.n_col;
    this->n_row = mat.n_row;
    if(this->data != NULL){
        this->deallocate();
    }
    this->data = allocate(this->n_row, this->n_col);
    for(int i = 0; i < mat.n_row * mat.n_col; i++){
        this->data[i] = mat.data[i];
    }
    return *this;
}

/*  overwrite current data array with "data"; returns self  */
template <class T>
Matrix<T>& Matrix<T>::operator= (const T* data){
    for(int i = 0; i < this->n_row * this->n_col; i++){
        this->data[i] = data[i];
    }
    return *this;
}

/*  converts Matrix to its transpose  */
template <class T>
void Matrix<T>::transpose(){
    T* new_data = new T[this->n_col * this->n_row]();
    for(int i = 0; i < this->n_row; i++){
        for(int j = 0; j < this->n_col; j++){
            new_data[j*this->n_row + i] = this->get(i,j);
        }
    }
    this->deallocate();
    this->data = new_data;
    std::swap(this->n_row, this->n_col);
    return;
}

/*  creates a deep copy of the current Matrix  */
template <class T> 
Matrix<T> Matrix<T>::deep_copy() const {
    Matrix<T> new_mat(this->n_row, this->n_col);
    new_mat.data = deep_copy(this->data, this->n_row * this->n_col);
    return new_mat;
}

/*  element-wise multiplication  */
template <class T>
Matrix<T> Matrix<T>::mul(const Matrix<T>& mat){
    Matrix<T> new_mat(this->n_row, this->n_col);
    for(int i = 0; i < this->n_row * this->n_col; i++){
        new_mat.data[i] = this->data[i] * mat.data[i];
    }
    return new_mat;
}

/*  element-wise division  */
template <class T>
Matrix<T> Matrix<T>::div(const Matrix<T>& mat){
    Matrix<T> new_mat(this->n_row, this->n_col);
    for(int i = 0; i < this->n_row * this->n_col; i++){
        if(mat.data[i] == 0){
            throw "Matrix::div ERROR: Divide by zero";
        }
        new_mat.data[i] = this->data[i] / mat.data[i];
    }
    return new_mat;
}

/*  element-wise addition  */
template <class T>
Matrix<T> Matrix<T>::operator+ (const Matrix<T>& mat){
    Matrix<T> new_mat(this->n_row, this->n_col);
    for(int i = 0; i < this->n_row * this->n_col; i++){
        new_mat.data[i] = this->data[i] + mat.data[i];
    }
    return new_mat;
}

/*  element-wise subtraction  */
template <class T>
Matrix<T> Matrix<T>::operator- (const Matrix<T>& mat){
    Matrix<T> new_mat(this->n_row, this->n_col);
    for(int i = 0; i < this->n_row * this->n_col; i++){
        new_mat.data[i] = this->data[i] - mat.data[i];
    }
    return new_mat;
}

/*  matrix multiplication  */
template <class T>
Matrix<T> Matrix<T>::operator* (const Matrix<T>& mat){
    if(this->n_col != mat.n_row){
        throw "Matrix::operator* ERROR: Dimension error";
    }
    Matrix<T> output(this->n_row, mat.n_col);
    for(int i = 0; i < this->n_row; i++){
        for(int j = 0; j < mat.n_col; j++){
            Matrix<T> mat1 = this->get(i,i,0,this->n_col-1);
            Matrix<T> mat2 = mat.get(0,mat.n_row-1,j,j);
            output.data[output.get_index(i,j)] 
                = Basic::dot(mat1.data, mat2.data, this->n_col);
        }
    }
    return output;
}

/*  sum of all elements in matrix  */
template <class T>
T Matrix<T>::sum(){
    return Basic::sum(this->data,this->n_col * this->n_row);
}

/*  mean of all elements in matrix  */
template <class T>
T Matrix<T>::mean(){
    return Basic::mean(this->data,this->n_col * this->n_row);
}

/*  check if matrix is square; true -> square  */
template <class T>
bool Matrix<T>::is_square(){
    return this->n_col == this->n_row;
}

/*  check if matrix is symmetric; true -> symmetric  */
template <class T>
bool Matrix<T>::is_symmetric(){
    if(this->is_square() == false){
        return false;
    }
    for(int i = 0; i < this->n_row; i++){
        for(int j = 0; j < i; j++){
            if(this->get(i,j) == this->get(j,i)){
                continue;
            } else {
                return false;
            }
        }
    }
    return true;
}

/*  check if two matrices are identical  */
template <class T>
bool Matrix<T>::operator== (const Matrix& mat){
    if(this->n_col != mat.n_col || this->n_row != mat.n_row){
        return false;
    }
    for(int i = 0; i < this->n_row * this->n_col; i++){
        if(std::abs(this->data[i] - mat.data[i]) < 1e-10){
            continue;
        } else {
            return false;
        }
    }
    return true;
}

/*  write current matrix to csv file  */
template <class T>
void Matrix<T>::write_mat(const std::string& filepath){
    std::stringstream ss;
    std::ofstream file;
    file.open(filepath);
    ss << this->n_row << "," << this->n_col << std::endl;
    for(int i = 0; i < this->n_col; i++){
        for(int j = 0; j < this->n_row; j++){
            ss << this->get(i,j) << ",";
        } 
        ss << std::endl;
    }
    file << ss.str();
    file.close();
    return;
}

/*  convert all elements to int type  */
template <class T>
Matrix<int> Matrix<T>::to_int(){
    Matrix<int> new_mat(this->n_row, this->n_col);
    new_mat.set(Basic::to_int(this->data, this->n_row * this->n_col));
    return new_mat;
}

/*  convert all elements to double type  */
template <class T>
Matrix<double> Matrix<T>::to_double(){
    Matrix<double> new_mat(this->n_row, this->n_col);
    new_mat.set(Basic::to_double(this->data, this->n_row * this->n_col));
    return new_mat;
}

/*  convert all elements to float type  */
template <class T>
Matrix<float> Matrix<T>::to_float(){
    Matrix<float> new_mat(this->n_row, this->n_col);
    new_mat.set(Basic::to_float(this->data, this->n_row * this->n_col));
    return new_mat;
}

// -----------------------------------------------------------------------------
// Non-Member Functions
// -----------------------------------------------------------------------------

/*  overload +: matrix + value  */
template <class T>
Matrix<T> operator+ (Matrix<T> mat, T value){
    for(int i = 0; i < mat.get_n_col() * mat.get_n_row(); i++){
        mat.set(i, mat.get(i) + value); 
    }
    return mat;
}

/*  overload +: value + matrix  */
template <class T>
Matrix<T> operator+ (T value, Matrix<T> mat){
    for(int i = 0; i < mat.get_n_col() * mat.get_n_row(); i++){
        mat.set(i, mat.get(i) + value); 
    }
    return mat;
}

/*  overload -: matrix - value  */
template <class T>
Matrix<T> operator- (Matrix<T> mat, T value){
    for(int i = 0; i < mat.get_n_col() * mat.get_n_row(); i++){
        mat.set(i, mat.get(i) - value); 
    }
    return mat;
}

/*  overload -: value - matrix  */
template <class T>
Matrix<T> operator- (T value, Matrix<T> mat){
    for(int i = 0; i < mat.get_n_col() * mat.get_n_row(); i++){
        mat.set(i, mat.get(i) - value); 
    }
    return mat;
}

/*  overload *: matrix * value  */
template <class T>
Matrix<T> operator* (Matrix<T> mat, T value){
    for(int i = 0; i < mat.get_n_col() * mat.get_n_row(); i++){
        mat.set(i, mat.get(i) * value); 
    }
    return mat;
}

/*  overload *: value * matrix  */
template <class T>
Matrix<T> operator* (T value, Matrix<T> mat){
    for(int i = 0; i < mat.get_n_col() * mat.get_n_row(); i++){
        mat.set(i, mat.get(i) * value); 
    }
    return mat;
}

/*  overload /: value / matrix  */
template <class T>
Matrix<T> operator/ (T value, Matrix<T> mat){
    for(int i = 0; i < mat.get_n_col() * mat.get_n_row(); i++){
        if(value == 0){
            throw "Matrix::operator/ ERROR: Division by zero";
        }
        mat.set(i, mat.get(i) / value); 
    }
    return mat;
}

/*  overload /: value / matrix  */
template <class T>
Matrix<T> operator/ (Matrix<T> mat, T value){
    for(int i = 0; i < mat.get_n_col() * mat.get_n_row(); i++){
        if(value == 0){
            throw "Matrix::operator/ ERROR: Division by zero";
        }
        mat.set(i, mat.get(i) / value); 
    }
    return mat;
}

/*  creates a matrix containing the transpose of mat  */
template <class T>
Matrix<T> transpose(Matrix<T> mat){
    mat.transpose();
    return mat;
}

#endif