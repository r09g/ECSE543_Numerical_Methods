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
#include <cassert>
#include <new>
#include "Basic.h"

template <typename T>
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
        Matrix(int n_row, int n_col);
        Matrix(int n_row, int n_col, const T* data);
        Matrix(const Matrix<T>& mat);

        // show entire matrix
        void show() const;
        
        // row, col to index conversion
        int get_index(int row, int col) const;

        // getters and setters
        T* get() const;
        void set(const T* data, unsigned int length);
        void set(T*&& data);
        T get(int row, int col) const;
        void set(int row, int col, T value);
        T get(int index);
        void set(int index, T value);
        Matrix<T> get(int row1, int row2, int col1, int col2) const;
        void set(int row, int col, const Matrix<T>& mat);
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
        static Matrix<double>& read_mat(const std::string& filepath);
        void write_mat(const std::string& filepath);

        // destructor
        ~Matrix();
};

// -----------------------------------------------------------------------------
// Static Member-Functions
// -----------------------------------------------------------------------------
template <typename T>
T* Matrix<T>::deep_copy(const T* data, unsigned int length){
    T* new_data = new T[length]();
    for(int i = 0; i < length; i++){
        new_data[i] = data[i];
    }
    return new_data;
}

template <typename T>
Matrix<double>& Matrix<T>::read_mat(const std::string& filepath){
    
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
    
    Matrix<double>* new_matrix = new Matrix<double>(n_row, n_col);
    int row = 0;
    int col = 0;
    while(getline(file,line)){
        std::stringstream ss(line);
        col = 0;
        while(getline(ss, output, ',')){
            new_matrix->set(row, col, stod(output));
            col++;
        }
        row++;
    }
    file.close();
    return *new_matrix;
}

// -----------------------------------------------------------------------------
// Constructors  & Destructor
// -----------------------------------------------------------------------------
template <typename T>
Matrix<T>::Matrix(): n_row(0), n_col(0) {
    this->data = allocate(n_row, n_col);
    return;
}

template <typename T>
Matrix<T>::Matrix(int n_row, int n_col): n_row(n_row), n_col(n_col) {
    this->data = allocate(n_row, n_col);
    return;
}

template <typename T>
Matrix<T>::Matrix(int n_row, int n_col, const T* data){
    this->n_row = n_row;
    this->n_col = n_col;
    this->data = allocate(this->n_row, this->n_col);
    for(int i = 0; i < n_row*n_col; i++){
        this->data[i] = data[i];
    }
    return;
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T>& mat): n_row(mat.n_row), n_col(mat.n_col){
    this->data = deep_copy(mat.data, this->n_row * this->n_col); 
    return;
}

template <typename T>
Matrix<T>::~Matrix(){
    this->deallocate();
    // std::cout << "deleted" << std::endl;
    return;
}

// -----------------------------------------------------------------------------
// Member-Functions
// -----------------------------------------------------------------------------
template <typename T>
T* Matrix<T>::allocate(int n_row, int n_col){
    T* new_data = NULL;
    try{
        new_data = new T[this->n_row * this->n_col]();
    } catch(std::bad_alloc& e){
        std::cout << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    // std::cout << "created" << std::endl;
    return new_data;
}

template <typename T>
void Matrix<T>::deallocate(){
    delete [] this->data;
    this->data = NULL;
    // std::cout << "deallocated" << std::endl;
    return;
}

template <typename T>
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

template <typename T>
T* Matrix<T>::get() const {
    return deep_copy(this->data, this->n_col*this->n_row);
}

template <typename T>
void Matrix<T>::set(const T* data, unsigned int length){
    if(this->data != NULL){
        this->deallocate();
    }
    this->data = deep_copy(data, length);
    return;
}

template <typename T>
void Matrix<T>::set(T*&& data){
    if(this->data != NULL){
        this->deallocate();
    }
    this->data = data;
    data = NULL;
    return;
}

template <typename T>
T Matrix<T>::get(int row, int col) const {
    return this->data[this->get_index(row,col)];
}

template <typename T>
void Matrix<T>::set(int row, int col, T value){
    this->data[this->get_index(row,col)] = value;
    return;
}

template <typename T>
T Matrix<T>::get(int index){
    return this->data[index];
}

template <typename T>
void Matrix<T>::set(int index, T value){
    this->data[index] = value;
    return;
}

template <typename T>
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

template <typename T>
void Matrix<T>::set(int row, int col, const Matrix<T>& mat){
    for(int i = 0; i < mat.n_row; i++){
        for(int j = 0; j < mat.n_col; j++){
            this->data[this->get_index(row + i, col + j)]
                = mat.data[i*mat.n_col + j];
        }
    }
    return;
}

template <typename T>
Matrix<T> Matrix<T>::get_row(int row) const {
    Matrix<T> new_mat(1, this->n.col);
    for(int i = 0; i < this->n.col; i++){
        new_mat.data[i] = this->get(row, i);
    }
    return new_mat;
}

template <typename T>
Matrix<T> Matrix<T>::get_col(int col) const {
    Matrix<T> new_mat(this->n.col, 1);
    for(int i = 0; i < this->n.row; i++){
        new_mat.data[i] = this->get(i, col);
    }
    return new_mat;
}

template <typename T>
int Matrix<T>::get_index(int row, int col) const{
    return row*this->n_col + col;
}

template <typename T>
int Matrix<T>::get_n_row() const{
    return this->n_row;
}

template <typename T>
int Matrix<T>::get_n_col() const{
    return this->n_col;
}

template <typename T>
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

template <typename T>
Matrix<T>& Matrix<T>::operator= (const T* data){
    for(int i = 0; i < this->n_row * this->n_col; i++){
        this->data[i] = data[i];
    }
    return *this;
}

template <typename T>
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

template <typename T> 
Matrix<T> Matrix<T>::deep_copy() const {
    Matrix<T> new_mat(this->n_row, this->n_col);
    new_mat.data = deep_copy(this->data, this->n_row * this->n_col);
    return new_mat;
}

template <typename T>
Matrix<T> Matrix<T>::mul(const Matrix<T>& mat){
    Matrix<T> new_mat(this->n_row, this->n_col);
    for(int i = 0; i < this->n_row * this->n_col; i++){
        new_mat.data[i] = this->data[i] * mat.data[i];
    }
    return new_mat;
}

template <typename T>
Matrix<T> Matrix<T>::div(const Matrix<T>& mat){
    Matrix<T> new_mat(this->n_row, this->n_col);
    for(int i = 0; i < this->n_row * this->n_col; i++){
        assert(mat.data[i] != 0);
        new_mat.data[i] = this->data[i] / mat.data[i];
    }
    return new_mat;
}

template <typename T>
Matrix<T> Matrix<T>::operator+ (const Matrix<T>& mat){
    Matrix<T> new_mat(this->n_row, this->n_col);
    for(int i = 0; i < this->n_row * this->n_col; i++){
        new_mat.data[i] = this->data[i] + mat.data[i];
    }
    return new_mat;
}

template <typename T>
Matrix<T> Matrix<T>::operator- (const Matrix<T>& mat){
    Matrix<T>* new_mat(this->n_row, this->n_col);
    for(int i = 0; i < this->n_row * this->n_col; i++){
        new_mat.data[i] = this->data[i] - mat.data[i];
    }
    return new_mat;
}

template <typename T>
Matrix<T> Matrix<T>::operator* (const Matrix<T>& mat){
    assert(this->n_col == mat.n_row);
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

template <typename T>
T Matrix<T>::sum(){
    return Basic::sum(this->data,this->n_col * this->n_row);
}

template <typename T>
T Matrix<T>::mean(){
    return Basic::mean(this->data,this->n_col * this->n_row);
}

template <typename T>
bool Matrix<T>::is_square(){
    return this->n_col == this->n_row;
}

template <typename T>
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

template <typename T>
bool Matrix<T>::operator== (const Matrix& mat){
    if(this->n_col != mat.n_col || this->n_row != mat.n_row){
        return false;
    }
    for(int i = 0; i < this->n_row * this->n_col; i++){
        if(this->data[i] == mat.data[i]){
            continue;
        } else {
            return false;
        }
    }
    return true;
}

template <typename T>
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

template <typename T>
Matrix<int> Matrix<T>::to_int(){
    Matrix<int> new_mat(this->n_row, this->n_col);
    new_mat.set(Basic::to_int(this->data, this->n_row * this->n_col));
    return new_mat;
}

template <typename T>
Matrix<double> Matrix<T>::to_double(){
    Matrix<double> new_mat(this->n_row, this->n_col);
    new_mat.set(Basic::to_double(this->data, this->n_row * this->n_col));
    return new_mat;
}

template <typename T>
Matrix<float> Matrix<T>::to_float(){
    Matrix<float> new_mat(this->n_row, this->n_col);
    new_mat.set(Basic::to_float(this->data, this->n_row * this->n_col));
    return new_mat;
}

// -----------------------------------------------------------------------------
// Non-Member Functions
// -----------------------------------------------------------------------------
template <typename T>
Matrix<T> operator+ (Matrix<T> mat, T value){
    for(int i = 0; i < mat.get_n_col() * mat.get_n_row(); i++){
        mat.set(i, mat.get(i) + value); 
    }
    return mat;
}

template <typename T>
Matrix<T> operator+ (T value, Matrix<T> mat){
    for(int i = 0; i < mat.get_n_col() * mat.get_n_row(); i++){
        mat.set(i, mat.get(i) + value); 
    }
    return mat;
}

template <typename T>
Matrix<T> operator- (Matrix<T> mat, T value){
    for(int i = 0; i < mat.get_n_col() * mat.get_n_row(); i++){
        mat.set(i, mat.get(i) - value); 
    }
    return mat;
}

template <typename T>
Matrix<T> operator- (T value, Matrix<T> mat){
    for(int i = 0; i < mat.get_n_col() * mat.get_n_row(); i++){
        mat.set(i, mat.get(i) - value); 
    }
    return mat;
}

template <typename T>
Matrix<T> operator* (Matrix<T> mat, T value){
    for(int i = 0; i < mat.get_n_col() * mat.get_n_row(); i++){
        mat.set(i, mat.get(i) * value); 
    }
    return mat;
}

template <typename T>
Matrix<T> operator* (T value, Matrix<T> mat){
    for(int i = 0; i < mat.get_n_col() * mat.get_n_row(); i++){
        mat.set(i, mat.get(i) * value); 
    }
    return mat;
}

template <typename T>
Matrix<T> operator/ (T value, Matrix<T> mat){
    for(int i = 0; i < mat.get_n_col() * mat.get_n_row(); i++){
        assert(value != 0);
        mat.set(i, mat.get(i) / value); 
    }
    return mat;
}

template <typename T>
Matrix<T> operator/ (Matrix<T> mat, T value){
    for(int i = 0; i < mat.get_n_col() * mat.get_n_row(); i++){
        assert(value != 0);
        mat.set(i, mat.get(i) / value); 
    }
    return mat;
}

#endif