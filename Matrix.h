/*
    Author: Raymond Yang
    Date: 2020/09/09
*/

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include "Basic.h"

template <typename T>
class Matrix {
    private:
        // parameters
        int n_row;
        int n_col;
        T* data = NULL;

        // memory allocation
        T* allocate_memory(int n_row, int n_col);

        // deep copy
        Matrix<T>* deep_copy() const;

    public:
        // constructors
        Matrix(int n_row, int n_col);
        Matrix(int n_row, int n_col, const T* data);
        Matrix(const Matrix<T>& matrix);

        // show entire matrix
        void show() const;
        
        // row, col to index conversion
        int get_index(int row, int col) const;

        // getters and setters
        T* get() const;
        void set(const T* data, unsigned int length);
        T get(int row, int col) const;
        void set(int row, int col, T value);
        Matrix<T>* get(int row1, int row2, int col1, int col2) const;
        void set(int row, int col, const Matrix<T>& matrix);
        int get_n_row() const;
        int get_n_col() const;

        // copy (duplicate) matrix
        Matrix<T>& operator= (const Matrix<T>& matrix);
        static T* deep_copy(const T* data, unsigned int length);

        // matrix operations
        void transpose();
        Matrix<T>& mul(const Matrix<T>& matrix);
        Matrix<T>& div(const Matrix<T>& matrix);
        Matrix<T>& operator+ (const Matrix<T>& matrix);
        Matrix<T>& operator- (const Matrix<T>& matrix);
        Matrix<T>& operator* (const Matrix<T>& matrix);
        Matrix<T>& mul(const T value);
        Matrix<T>& div(const T value);
        Matrix<T>& operator+ (const T value);
        Matrix<T>& operator- (const T value);
        Matrix<T>& operator* (const T value);
        T sum();
        T mean();

        // check matrix characteristics
        bool is_square();
        bool is_symmetric();
        bool operator== (const Matrix& matrix);

        // matrix data type conversions
        Matrix<T>& to_int();
        Matrix<T>& to_double();
        Matrix<T>& to_float();

        // I/O to Excel
        static Matrix<double>& read_mat(const std::string& filepath);
        void write_mat(const std::string& filepath);

        // destructor
        ~Matrix();
};

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

template <typename T>
Matrix<T>::Matrix(int n_row, int n_col): n_row(n_row), n_col(n_col) {
    this->data = allocate_memory(n_row, n_col);
    return;
}

template <typename T>
Matrix<T>::Matrix(int n_row, int n_col, const T* data){
    this->n_row = n_row;
    this->n_col = n_col;
    this->data = allocate_memory(this->n_row, this->n_col);
    for(int i = 0; i < n_row*n_col; i++){
        this->data[i] = data[i];
    }
    return;
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T>& matrix){
    this->n_col = matrix.n_col;
    this->n_row = matrix.n_row;
    this->data = deep_copy(matrix.data, this->n_row * this->n_col); 
    return;
}

template <typename T>
Matrix<T>::~Matrix(){
    delete [] this->data;
    return;
}

template <typename T>
T* Matrix<T>::allocate_memory(int n_row, int n_col){
    T* new_data = new T[this->n_row * this->n_col]();
    return new_data;
}

template <typename T>
void Matrix<T>::show() const {
    for(int i = 0; i < this->n_row; i++){
        std::cout << "[ ";
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
        delete [] this->data;
    }
    this->data = deep_copy(data, length);
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
Matrix<T>* Matrix<T>::get(int row1, int row2, int col1, int col2) const {
    Matrix<T>* new_matrix = new Matrix<T>(row2 - row1 + 1, col2 - col1 + 1);
    for(int i = 0; i < new_matrix->n_row; i++){
        for(int j = 0; j < new_matrix->n_col; j++){
            new_matrix->data[new_matrix->get_index(i,j)] 
                = this->get(row1 + i, col1 + j);
        }
    }
    return new_matrix;
}

template <typename T>
void Matrix<T>::set(int row, int col, const Matrix<T>& matrix){
    for(int i = 0; i < matrix.n_row; i++){
        for(int j = 0; j < matrix.n_col; j++){
            this->data[this->get_index(row + i, col + j)]
                = matrix.data[i*matrix.n_col + j];
        }
    }
    return;
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
Matrix<T>& Matrix<T>::operator= (const Matrix<T>& matrix){
    this->n_col = matrix.n_col;
    this->n_row = matrix.n_row;
    if(this->data != NULL){
        delete [] this->data;
    }
    this->data = allocate_memory(this->n_row, this->n_col);
    for(int i = 0; i < matrix.n_row * matrix.n_col; i++){
        this->data[i] = matrix.data[i];
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
    delete [] this->data;
    this->data = new_data;
    std::swap(this->n_row, this->n_col);
    return;
}

template <typename T> 
Matrix<T>* Matrix<T>::deep_copy() const {
    Matrix<T>* new_matrix = new Matrix<T>(this->n_row, this->n_col);
    new_matrix->data = deep_copy(this->data, this->n_row * this->n_col);
    return new_matrix;
}

template <typename T>
Matrix<T>& Matrix<T>::mul(const Matrix<T>& matrix){
    Matrix<T>* new_matrix = new Matrix<T>(this->n_row, this->n_col);
    for(int i = 0; i < this->n_row * this->n_col; i++){
        new_matrix->data[i] = this->data[i] * matrix.data[i];
    }
    return *new_matrix;
}

template <typename T>
Matrix<T>& Matrix<T>::div(const Matrix<T>& matrix){
    Matrix<T>* new_matrix = new Matrix<T>(this->n_row, this->n_col);
    for(int i = 0; i < this->n_row * this->n_col; i++){
        new_matrix->data[i] = this->data[i] / matrix.data[i];
    }
    return *new_matrix;
}

template <typename T>
Matrix<T>& Matrix<T>::operator+ (const Matrix<T>& matrix){
    Matrix<T>* new_matrix = new Matrix<T>(this->n_row, this->n_col);
    for(int i = 0; i < this->n_row * this->n_col; i++){
        new_matrix->data[i] = this->data[i] + matrix.data[i];
    }
    return *new_matrix;
}

template <typename T>
Matrix<T>& Matrix<T>::operator- (const Matrix<T>& matrix){
    Matrix<T>* new_matrix = new Matrix<T>(this->n_row, this->n_col);
    for(int i = 0; i < this->n_row * this->n_col; i++){
        new_matrix->data[i] = this->data[i] - matrix.data[i];
    }
    return *new_matrix;
}

template <typename T>
Matrix<T>& Matrix<T>::mul(const T value){
    Matrix<T>* new_matrix = new Matrix<T>(this->n_row, this->n_col);
    for(int i = 0; i < this->n_row * this->n_col; i++){
        new_matrix->data[i] = this->data[i] * value;
    }
    return *new_matrix;
}

template <typename T>
Matrix<T>& Matrix<T>::div(const T value){
    Matrix<T>* new_matrix = new Matrix<T>(this->n_row, this->n_col);
    for(int i = 0; i < this->n_row * this->n_col; i++){
        new_matrix->data[i] = this->data[i] / value;
    }
    return *new_matrix;
}

template <typename T>
Matrix<T>& Matrix<T>::operator+ (const T value){
    Matrix<T>* new_matrix = new Matrix<T>(this->n_row, this->n_col);
    for(int i = 0; i < this->n_row * this->n_col; i++){
        new_matrix->data[i] = this->data[i] + value;
    }
    return *new_matrix;
}

template <typename T>
Matrix<T>& Matrix<T>::operator- (const T value){
    Matrix<T>* new_matrix = new Matrix<T>(this->n_row, this->n_col);
    for(int i = 0; i < this->n_row * this->n_col; i++){
        new_matrix->data[i] = this->data[i] - value;
    }
    return *new_matrix;
}

template <typename T>
Matrix<T>& Matrix<T>::operator* (const T value){
    Matrix<T>* new_matrix = new Matrix<T>(this->n_row, this->n_col);
    for(int i = 0; i < this->n_row * this->n_col; i++){
        new_matrix->data[i] = this->data[i] * value;
    }
    return *new_matrix;
}

template <typename T>
Matrix<T>& Matrix<T>::operator* (const Matrix<T>& matrix){
    if(this->n_col != matrix.n_row){
        std::cout << "Dimension error" << std::endl;
        exit(EXIT_FAILURE);    
    }

    Matrix<T>* output = new Matrix<T>(this->n_row, matrix.n_col);
    for(int i = 0; i < this->n_row; i++){
        for(int j = 0; j < matrix.n_col; j++){
            Matrix<T>* mat1 = this->get(i,i,0,this->n_col-1);
            Matrix<T>* mat2 = matrix.get(0,matrix.n_row-1,j,j);
            output->data[output->get_index(i,j)] 
                = basic::dot(mat1->data, mat2->data, this->n_col);
            delete mat1;
            delete mat2;
        }
    }
    return *output;
}

template <typename T>
T Matrix<T>::sum(){
    return basic::sum(this->data,this->n_col * this->n_row);
}

template <typename T>
T Matrix<T>::mean(){
    return basic::mean(this->data,this->n_col * this->n_row);
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
bool Matrix<T>::operator== (const Matrix& matrix){
    if(this->n_col != matrix.n_col || this->n_row != matrix.n_row){
        return false;
    }
    for(int i = 0; i < this->n_row * this->n_col; i++){
        if(this->data[i] == matrix.data[i]){
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
Matrix<T>& Matrix<T>::to_int(){
    Matrix<int>* new_matrix = new Matrix<int>(this->n_row, this->n_col);
    new_matrix->data = basic::to_int(this->data, this->n_row * this->n_col);
    return *new_matrix;
}

template <typename T>
Matrix<T>& Matrix<T>::to_double(){
    Matrix<double>* new_matrix = new Matrix<double>(this->n_row, this->n_col);
    new_matrix->data = basic::to_double(this->data, this->n_row * this->n_col);
    return *new_matrix;
}

template <typename T>
Matrix<T>& Matrix<T>::to_float(){
    Matrix<float>* new_matrix = new Matrix<float>(this->n_row, this->n_col);
    new_matrix->data = basic::to_float(this->data, this->n_row * this->n_col);
    return *new_matrix;
}



#endif