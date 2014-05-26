#ifndef MISC_CC
#define MISC_CC

#define RAND_A1 40014
#define RAND_M1 2147483563
#define RAND_Q1 53668
#define RAND_R1 12211

#define RAND_A2 40692
#define RAND_M2 2147483399
#define RAND_Q2 52744
#define RAND_R2 3791

#define RAND_SCALE1 (1.0 / RAND_M1)

#include <cmath>
#include <stdlib.h>
#include <stdio.h>

const double pi = 3.1415927;
const double au = 1.496e13;
const double pc = 3.086e18;
const double R_sun = 6.955e10;
const double M_sun = 1.98892e33;
const double L_sun = 3.84e33;
const double c_l = 2.99792458e10;
const double h = 6.6260755e-27;
const double k = 1.380658e-16;
const double sigma = 5.67051e-5;

/* Get a random number between 0 and 1. */

static int seed1 = 1, seed2 = 1;

double random_number() {
    //return ((double) rand() / (RAND_MAX));
    int k, result;

    k = seed1 / RAND_Q1;
    seed1 = RAND_A1 * (seed1 - k * RAND_Q1) - k * RAND_R1;
    if (seed1 < 0) seed1 += RAND_M1;

    k = seed2 / RAND_Q2;
    seed2 = RAND_A2 * (seed2 - k * RAND_Q2) - k * RAND_R2;
    if (seed2 < 0) seed2 += RAND_M2;

    result = seed1 - seed2;
    if (result < 1) result += RAND_M1 - 1;

    return result * (double) RAND_SCALE1;
};

/* Calculate the blackbody function for a given frequency and temperature. */

double planck_function(double nu, double T) {
    return 2.0*h*nu*nu*nu/(c_l*c_l)*1.0/(exp(h*nu/(k*T))-1.0);
};

/* Integrate an array y over array x using the trapezoidal method. */

double integrate(double *y, double *x, int nx) {
    double sum = 0.0;

    for (int i=0; i<nx-1; i++)
        sum += 0.5*(y[i+1]+y[i])*(x[i+1]-x[i]);

    return sum;
};

/* Test whether two values are equal within a given tolerance. */

bool equal(double x, double y, double tol) {
    if (fabs(x-y) < fabs(y)*tol)
        return true;
    else
        return false;
};

/* Find the cell in an array in which the given value is located using a 
   tree. */

int find_in_arr(double val, double *arr, int n) {
    int lmin = 0;
    int lmax = n-1;
    bool not_found = true;
    int l;

    while (not_found) {
        int ltest = (lmax-lmin)/2+lmin;

        if ((val >= arr[ltest]) && (val <= arr[ltest+1])) {
            l = ltest;
            not_found = false;
        }
        else {
            if (val < arr[ltest])
                lmax = ltest;
            else
                lmin = ltest;
        }
    }

    return l;
};

/* Find the cell in an array in which the given value is located. This 
 * function overloads the previous one by allowing you to search a smaller 
 * portion of the array. */

int find_in_arr(double val, double *arr, int lmin, int lmax) {
    int l;

    for (int i=lmin; i <= lmax; i++) {
        if ((val >= arr[i]) && (val <= arr[i+1]))
            l = i;
    }

    return l;
};

/* Find the cell in an array in which the given value is located, but in this
 * case the array is cyclic. */

int find_in_periodic_arr(double val, double *arr, int n, int lmin, int lmax) {
    int l;

    for (int i=lmin; i <= lmax; i++) {
        int index = (i+n)%(n);
        if ((val >= arr[index]) && (val <= arr[index+1])) {
            l = index;
        }
    }

    return l;
}

void swap (double* x, double* y) {
        double temp;
        temp = *x;
        *x = *y;
        *y = temp;
}

void bubbleSort(double arr [], int size) {
    int last = size - 2;
    int isChanged = 1, k;

    while ((last >= 0) && isChanged) {
        isChanged = 0;
        for (k = 0; k <= last; k++)
            if (arr[k] > arr[k+1]) {
                swap (&arr[k], &arr[k+1]);
                isChanged = 1;
            }
        last--;
    }
};

/* Create an empty 2-dimensional array. */

double **create2DArr(int nx, int ny) {
    double **arr = new double*[nx];
    for (int i=0; i<nx; i++)
        arr[i] = new double[ny];

    return arr;
};

/* Create an empty 3-dimensional array. */

double ***create3DArr(int nx, int ny, int nz) {
    double ***arr = new double**[nx];
    for (int i=0; i<nx; i++) {
        arr[i] = new double*[ny];
        for (int j=0; j<ny; j++)
            arr[i][j] = new double[nz];
    }

    return arr;
};

/* Create a 3-dimensional array filled with a particular value. */

double ***create3DArrValue(int nx, int ny, int nz, int value) {
    double ***arr = new double**[nx];
    for (int i=0; i<nx; i++) {
        arr[i] = new double*[ny];
        for (int j=0; j<ny; j++) {
            arr[i][j] = new double[nz];
            for (int k=0; k<nz; k++)
                arr[i][j][k] = value;
        }
    }

    return arr;
};

/* Set the value of a 3-dimensional array to a constant value. */

void set3DArrValue(double ***arr, double value, int nx, int ny, int nz) {
    for (int i=0; i<nx; i++)
        for (int j=0; j<ny; j++)
            for (int k=0; k<nz; k++)
                arr[i][j][k] = value;
}

/* Set one 3-dimensional array equal to another 3-dimensional array, element
 * by element. */

void equate3DArrs(double ***arr1, double ***arr2, int nx, int ny, int nz) {
    for (int i=0; i<nx; i++)
        for (int j=0; j<ny; j++)
            for (int k=0; k<nz; k++)
                arr1[i][j][k] = arr2[i][j][k];
}

double delta(double x1, double x2) {
    double d1 = x1/x2;
    double d2 = x2/x1;

    double delt = d1;
    if (d2 > d1) delt = d2;

    return delt;
}

double quantile(double ***R, double p, int nx, int ny, int nz) {
    double *Rline = new double[nx*ny*nz];

    for (int i=0; i<nx; i++)
        for (int j=0; j<ny; j++)
            for (int k=0; k<nz; k++)
                Rline[i*ny*nz+j*nz+k] = R[i][j][k];

    bubbleSort(Rline, nx*ny*nz);

    double quant = Rline[int(p*nx*ny*nz)];

    delete Rline;

    return quant;
}

bool converged(double ***newArr, double ***oldArr, double ***reallyoldArr,
        int n1, int n2, int n3) {
    double Qthresh = 2.0;
    double Delthresh = 1.1;
    double p = 0.99;

    double ***R = create3DArr(n1, n2, n3);
    double ***Rold = create3DArr(n1, n2, n3);

    for (int i=0; i<n1; i++) {
        for (int j=0; j<n2; j++) {
            for (int k=0; k<n3; k++) {
                R[i][j][k] = delta(oldArr[i][j][k],newArr[i][j][k]);
                Rold[i][j][k] = delta(reallyoldArr[i][j][k],newArr[i][j][k]);
            }
        }
    }

    double Q = quantile(R,p,n1,n2,n3);
    double Qold = quantile(Rold,p,n1,n2,n3);
    printf("%f   %f\n", Q, Qold);

    double Del = delta(Qold,Q);
    printf("%f\n", Del);

    bool conv = ((Q < Qthresh) && (Del < Delthresh));

    delete R;
    delete Rold;

    return conv;
}

#endif
