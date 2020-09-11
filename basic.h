/******************************************************************************/
/* Name: Basic.h                                                              */
/* Date: 2020/09/10                                                           */
/* Author: Raymond Yang                                                       */
/******************************************************************************/

#ifndef __BASIC__
#define __BASIC__

#include <iostream>
#include <cmath>

namespace Basic{

    template <typename T, size_t N> T sum(const T (&data)[N]);
    template <typename T> T sum(const T* data, unsigned int length);
    
    template <typename T, size_t N> double mean(const T (&data)[N]);
    template <typename T> double mean(const T* data, unsigned int length);
    
    template <typename T, size_t N> T product(const T (&data)[N]);
    template <typename T> T product(const T* data, unsigned int length);
    
    template <typename T, typename C, size_t N>
    double dot(const T (&data1)[N], const C (&data2)[N]);
    
    template <typename T, typename C>
    double dot(const T* data1, const C* data2, unsigned int length);
    
    template <typename T, size_t N> double* to_double(const T (&data)[N]);
    template <typename T> double* to_double(const T* data, unsigned int length);
    
    template <typename T, size_t N> int* to_int(const T (&data)[N]);
    template <typename T> int* to_int(const T* data, unsigned int length);
    
    template <typename T, size_t N> float* to_float(const T (&data)[N]);
    template <typename T> float* to_float(const T* data, unsigned int length);
    
};

// -----------------------------------------------------------------------------

template <typename T, size_t N>
T Basic::sum(const T (&data)[N]){
    T result = (T)(0.0);
    for(int i = 0; i < N; i++){
        result += data[i];
    }
    return result;
}

template <typename T> 
T Basic::sum(const T* data, unsigned int length){
    T result = (T)(0.0);
    for(int i = 0; i < length; i++){
        result += data[i];
    }
    return result;
}

template <typename T, size_t N>
double Basic::mean(const T (&data)[N]){
    double result = 0.0;
    result = (double)(Basic::sum(data)) / N;
    return result;
}

template <typename T>
double Basic::mean(const T* data, unsigned int length){
    double result = 0.0;
    result = (double)(Basic::sum(data, length)) / length;
    return result;
}

template <typename T, size_t N>
T Basic::product(const T (&data)[N]){
    T result = (T)(0.0);
    for(int i = 0; i < N; i++){
        result *= data[i];
    }
    return result;
}

template <typename T>
T Basic::product(const T* data, unsigned int length){
    T result = (T)(0.0);
    for(int i = 0; i < length; i++){
        result *= data[i];
    }
    return result;
}

template <typename T, typename C, size_t N>
double Basic::dot(const T (&data1)[N], const C (&data2)[N]){
    double result = 0.0;
    for(int i = 0; i < N; i++){
        result += data1[i] * data2[i];
    }
    return result;
}

template <typename T, typename C>
double Basic::dot(const T* data1, const C* data2, unsigned int length){
    double result = 0.0;
    for(int i = 0; i < length; i++){
        result += data1[i] * data2[i];
    }
    return result;
}

template <typename T, size_t N> 
double* Basic::to_double(const T (&data)[N]){
    double* output = new double[N]();
    for(int i = 0; i < N; i++){
        output[i] = (double)(data[i]);
    }
    return output;
}

template <typename T> 
double* Basic::to_double(const T* data, unsigned int length){
    double* output = new double[length]();
    for(int i = 0; i < length; i++){
        output[i] = (double)(data[i]);
    }
    return output;
}

template <typename T, size_t N> 
int* Basic::to_int(const T (&data)[N]){
    int* output = new int[N]();
    for(int i = 0; i < N; i++){
        output[i] = std::round(data[i]);
    }
    return output;
}

template <typename T> 
int* Basic::to_int(const T* data, unsigned int length){
    int* output = new int[length]();
    for(int i = 0; i < length; i++){
        output[i] = std::round(data[i]);
    }
    return output;
}

template <typename T, size_t N> 
float* Basic::to_float(const T (&data)[N]){
    float* output = new float[N]();
    for(int i = 0; i < N; i++){
        output[i] = (float)(data[i]);
    }
    return output;
}

template <typename T> 
float* Basic::to_float(const T* data, unsigned int length){
    float* output = new float[length]();
    for(int i = 0; i < length; i++){
        output[i] = (float)(data[i]);
    }
    return output;
}


#endif