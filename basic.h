/*
    Author: Raymond Yang
    Date: 2020/09/09
*/

#ifndef __BASIC__
#define __BASIC__

#include <iostream>

namespace basic{
    
    template <typename T, size_t N> T sum(const T (&data)[N]);
    template <typename T> T sum(const T* data, unsigned int length);
    template <typename T, size_t N> double mean(const T (&data)[N]);
    template <typename T, size_t N> T product(const T (&data)[N]);
    
    // dot product (overloaded)
    template <typename T, typename C, size_t N>
    double dot(const T (&data1)[N], const C (&data2)[N]);
    template <typename T, typename C>
    double dot(const T* data1, const C* data2, unsigned int length);
};

// -----------------------------------------------------------------------------

template <typename T, size_t N>
T basic::sum(const T (&data)[N]){
    T* result = new T();
    for(int i = 0; i < N; i++){
        *result += data[i];
    }
    return *result;
}

template <typename T> 
T basic::sum(const T* data, unsigned int length){
    T* result = new T();
    for(int i = 0; i < length; i++){
        *result += data[i];
    }
    return *result;
}

template <typename T, size_t N>
double& basic::mean(const T (&data)[N]){
    double* result = new double();
    *result = (double)(basic::sum(data)) / N;
    return *result;
}

template <typename T, size_t N>
T basic::product(const T (&data)[N]){
    T* result = new T();
    for(int i = 0; i < N; i++){
        *result *= data[i];
    }
    return *result;
}

template <typename T, typename C, size_t N>
double basic::dot(const T (&data1)[N], const C (&data2)[N]){
    double* result = new double(0.0);
    for(int i = 0; i < N; i++){
        *result += data1[i] * data2[i];
    }
    return *result;
}

template <typename T, typename C>
double basic::dot(const T* data1, const C* data2, unsigned int length){
    double* result = new double(0.0);
    for(int i = 0; i < length; i++){
        *result += data1[i] * data2[i];
    }
    return *result;
}


#endif