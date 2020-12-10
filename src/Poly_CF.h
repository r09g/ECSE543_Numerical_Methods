/******************************************************************************/
/* Name: Poly_CF.h                                                            */
/* Description: Polynomial Curve Fitting                                      */
/* Date: 2020/12/09                                                           */
/* Author: Raymond Yang                                                       */
/******************************************************************************/

#ifndef __Poly_CF__
#define __Poly_CF__

#include "Matrix.h"
#include "Matrix_Solver.h"

template <typename T = double>
class Poly_CF {
    private:
        Matrix<T> xdata;
        Matrix<T> ydata;
        int ndata;

        // helper functions
        Matrix<T> diff(int type = -1);

    public:
        Poly_CF(Matrix<T> xdata, Matrix<T> ydata);
        ~Poly_CF();

        // getters and setters
        Matrix<T> get_x();
        Matrix<T> get_y();
        int get_n();
        
        // Fitting Methods
        Matrix<T> Lagrange_WD();
        Matrix<T> Hermite_3(Matrix<T> x_pts);

        // 
        
};

// -----------------------------------------------------------------------------
// Constructor & Destructor
// -----------------------------------------------------------------------------

/*  Constructor. Input data must be 1-D.  */
template <typename T>
Poly_CF<T>::Poly_CF(Matrix<T> xdata, Matrix<T> ydata): xdata(xdata), 
    ydata(ydata){

    this->ndata = (xdata.get_n_row())*(xdata.get_n_col());

    if(this->ndata < 2){
        throw "Poly_CF::Poly_CF ERROR: Less than 2 data points";
    }

    if(xdata.get_n_col() != 1 && xdata.get_n_row() != 1 || 
        ydata.get_n_col() != 1 && ydata.get_n_row() != 1){
        throw "Poly_CF::Poly_CF ERROR: Data is not 1-D"; 
    }
    
    if((xdata.get_n_col())*(xdata.get_n_row()) !=
        (ydata.get_n_col())*(ydata.get_n_row())){
        throw "Poly_CF::Poly_CF ERROR: x and y have different sizes"; 
    }
    return;
}

/*  Destructor  */
template <typename T>
Poly_CF<T>::~Poly_CF(){}

// -----------------------------------------------------------------------------
// Member Functions
// -----------------------------------------------------------------------------

/*  
    Compute the discrete derivative of the data  
    
    Input: type of derivative.
        0: (i+1) and (i)
        1: (i+1) and (i-1)
        2: inverse y-weighted (i+1) and (i-1)
        default: type 0

    Returns: column vector of derivatives at each data point
*/
template <typename T>
Matrix<T> Poly_CF<T>::diff(int type){
    Matrix<T> df(this->ndata, 1);
    if(type == -1 || type == 0){  // (i+1) and (i)
        for(int i = 0; i < this->ndata - 1; i++){
            df.set(i, (this->ydata.get(i+1) - this->ydata.get(i))/
                (this->xdata.get(i+1) - this->xdata.get(i)));
        }
        df.set(this->ndata - 1, df.get(this->ndata - 2));
    } else if(type == 1) {  // (i+1) and (i-1)
        if(this->ndata < 3){
            throw "Poly_CF::diff ERROR: Less than 3 data points";
        }
        for(int i = 0; i < this->ndata; i++){
            int prev = i - 1;
            if(prev < 0){
                prev = 0;
            }
            int post = i + 1;
            if(post > this->ndata - 1){
                post = this->ndata - 1;
            }
            df.set(i, (this->ydata.get(post) - this->ydata.get(prev))/
                (this->xdata.get(post) - this->xdata.get(prev)));
        }
    } else {  //  inverse y-weighted (i+1) and (i-1)
        df.set(0, (this->ydata.get(1) - this->ydata.get(0))/
            (this->xdata.get(1) - this->xdata.get(0)));
        int end = this->ndata - 1;
        df.set(end, (this->ydata.get(end) - this->ydata.get(end-1))/
            (this->xdata.get(end) - this->xdata.get(end-1)));

        for(int i = 1; i < this->ndata - 1; i++){
            double w1 = (this->ydata.get(i+1) - this->ydata.get(i))/
                (this->ydata.get(i+1) - this->ydata.get(i-1));
            double w2 = (this->ydata.get(i) - this->ydata.get(i-1))/
                (this->ydata.get(i+1) - this->ydata.get(i-1));

            df.set(i, w1*(this->ydata.get(i) - this->ydata.get(i-1))/
                (this->xdata.get(i) - this->xdata.get(i-1)) + w2*
                (this->ydata.get(i+1) - this->ydata.get(i))/
                (this->xdata.get(i+1) - this->xdata.get(i)));
        }
    }
    return df;
}

/*  Get x data  */
template <typename T>
Matrix<T> Poly_CF<T>::get_x(){
    return this->xdata;
}

/*  Get y data  */
template <typename T>
Matrix<T> Poly_CF<T>::get_y(){
    return this->ydata;
}

/*  Get total number of data points given  */
template <typename T>
int Poly_CF<T>::get_n(){
    return this->ndata;
}

/*  
    Lagrange Polynomial Whole Domain Fitting

    Returns: Column vector of coefficients representing polynomial from lowest
             order to highest order x^0 + ... + x^(n-1).
*/
template <typename T>
Matrix<T> Poly_CF<T>::Lagrange_WD(){
    Matrix<T> G(this->ndata);
    Matrix<T> b(this->ndata, 1);
    for(int i = 0; i < this->ndata; i++){
        for(int j = 0; j < this->ndata; j++){
            G.set(i, j, ((this->xdata)^(i+j)).sum());
        }
        b.set(i, (((this->xdata)^(i)).mul(this->ydata)).sum());
    }

    Matrix_Solver::cholesky_solve(&G, &b);
    return b;
}

/*  
    Cubic Hermite Polynomial Sub-domain Fitting

    Input: vector of given data points + interpolation points

    Returns: Column vector of y data points for curve.
*/
template <typename T>
Matrix<T> Poly_CF<T>::Hermite_3(Matrix<T> x_pts){
    
    auto u1 = [](T x, T x1, T x2){
        return (1 - 2*(x - x1)/(x1 - x2))*pow((x - x2)/(x1 - x2),2);
    };  
    auto u2 = [](T x, T x1, T x2){
        return (1 - 2*(x - x2)/(x2 - x1))*pow((x - x1)/(x2 - x1),2);
    };  
    auto v1 = [](T x, T x1, T x2){
        return (x - x1)*pow((x - x2)/(x1 - x2),2);
    };  
    auto v2 = [](T x, T x1, T x2){
        return (x - x2)*pow((x - x1)/(x2 - x1),2);;
    }; 
    
    Matrix<T> df = diff(1);

    Matrix<T> y(x_pts.get_n_col() * x_pts.get_n_row(), 1);
    int i = 0;
    for(int sd = 0; sd < this->ndata - 1; sd++){
        int i_next = x_pts.find(this->xdata.get(sd+1));
        double x1 = this->xdata.get(sd);
        double x2 = this->xdata.get(sd+1);
        for(i; i <= i_next; i++){
            y.set(i, this->ydata.get(sd)*u1(x_pts.get(i), x1, x2) +
                this->ydata.get(sd+1)*u2(x_pts.get(i), x1, x2) + 
                df.get(sd)*v1(x_pts.get(i), x1, x2) + 
                df.get(sd+1)*v2(x_pts.get(i), x1, x2)); 
        }

    } 

    return y;
}


#endif








