#include "misc.h"

/* Get a random number between 0 and 1. */

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
    double value = 2.0*h*nu*nu*nu/(c_l*c_l)*1.0/(exp(h*nu/(k*T))-1.0);

    if (std::isnan(value))
        return 0;
    else
        return value;
};

double planck_function_derivative(double nu, double T) {
    double value = (-2.0*h*nu*nu*nu*nu)/(c_l*c_l*k*T*T) / 
        (exp(h*nu/(k*T))-1.0) / (1. - exp(-h*nu/(k*T)));

    if (isnan(value))
        return 0;
    else
        return value;
};

/* Integrate an array y over array x using the trapezoidal method. */

double integrate(double *y, double *x, int nx) {
    double sum = 0.0;

    for (int i=0; i<nx-1; i++)
        sum += 0.5*(y[i+1]+y[i])*(x[i+1]-x[i]);

    return sum;
};

/* Cumulatively integrate. */

double *cumulative_integrate(double *y, double *x, int nx) {
    double sum = 0;
    double* cum_sum = new double[nx];

    cum_sum[0] = sum;
    for (int i = 1; i < nx; i++) {
        sum += 0.5*(y[i]+y[i-1])*(x[i]-x[i-1]);
        cum_sum[i] = sum;
    }

    for (int i = 0; i < nx; i++)
        cum_sum[i] /= sum;

    return cum_sum;
}

/* Take the derivative of an array. */

double *derivative(double *y, double *x, int nx) {
    double *result = new double[nx-1];

    for (int i = 0; i < nx-1; i++)
        result[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]);

    return result;
}

double **derivative2D_ax0(double **y, double *x, int nx, int ny) {
    double **result = create2DArr(nx-1, ny);

    for (int i = 0; i < nx-1; i++)
        for (int j = 0; j < ny; j++)
            result[i][j] = (y[i+1][j] - y[i][j]) / (x[i+1] - x[i]);

    return result;
}

/* Test whether two values are equal within a given tolerance. */

bool equal(double x, double y, double tol) {
    if (fabs(x-y) < fabs(y)*tol)
        return true;
    else
        return false;
};

bool equal_zero(double x, double tol) {
    if (fabs(x) < tol)
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
    int l = -1;

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

void delete2DArr(double **arr, int nx, int ny) {
    for (int i = 0; i < nx; i++)
        delete[] arr[i];
    delete[] arr;
}

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

/* Delete a 3-dimensional array. */

void delete3DArr(double ***arr, int nx, int ny, int nz) {
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++)
            delete[] arr[i][j];
        delete[] arr[i];
    }
    delete[] arr;
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

/* Create an empty 4-dimensional array. */

std::vector<double***> create4DArr(int nx, int ny, int nz, int nq) {
    std::vector<double***> arr;
    for (int i=0; i<nx; i++)
        arr.push_back(create3DArr(ny, nz, nq));

    return arr;
};

void delete4DArr(std::vector<double***> arr, int nx, int ny, int nz, int nq) {
    for (int i = 0; i < nx; i++)
        delete3DArr(arr[i], ny, nz, nq);
    arr.clear();
}

void delete4DArr(double ****arr, int nx, int ny, int nz, int nq) {
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            for (int k=0; k<nz; k++)
                delete[] arr[i][j][k];
            delete[] arr[i][j];
        }
        delete[] arr[i];
    }
    delete[] arr;
};

/* Create a 4-dimensional array filled with a particular value. */

double ****create4DArrValue(int nx, int ny, int nz, int nq, double value) {
    double ****arr = new double***[nx];
    for (int i=0; i<nx; i++) {
        arr[i] = new double**[ny];
        for (int j=0; j<ny; j++) {
            arr[i][j] = new double*[nz];
            for (int k=0; k<nz; k++) {
                arr[i][j][k] = new double[nq];
                for (int l=0; l<nq; l++)
                    arr[i][j][k][l] = value;
            }
        }
    }

    return arr;
};

/* Set the value of a 4-dimensional array to a constant value. */

void set4DArrValue(std::vector<double***> arr, double value, int nx, int ny, 
        int nz, int nq) {
    for (int i=0; i<nx; i++)
        for (int j=0; j<ny; j++)
            for (int k=0; k<nz; k++)
                for (int l=0; l<nq; l++)
                    arr[i][j][k][l] = value;
}

void set4DArrValue(double**** arr, double value, int nx, int ny, int nz, 
        int nq) {
    for (int i=0; i<nx; i++)
        for (int j=0; j<ny; j++)
            for (int k=0; k<nz; k++)
                for (int l=0; l<nq; l++)
                    arr[i][j][k][l] = value;
}

/* Set one 4-dimensional array equal to another 4-dimensional array, element
 * by element. */

void equate4DArrs(std::vector<double***> arr1, std::vector<double***> arr2, 
        int nx, int ny, int nz, int nq) {
    for (int i=0; i<nx; i++)
        for (int j=0; j<ny; j++)
            for (int k=0; k<nz; k++)
                for (int l=0; l<nq; l++)
                    arr1[i][j][k][l] = arr2[i][j][k][l];
}

double delta(double x1, double x2) {
    double d1 = x1/x2;
    double d2 = x2/x1;

    double delt = d1;
    if (d2 > d1) delt = d2;

    return delt;
}

double quantile(std::vector<double***> R, double p, int nx, int ny, int nz, 
        int nq) {
    double *Rline = new double[nx*ny*nz*nq];

    for (int i=0; i<nx; i++)
        for (int j=0; j<ny; j++)
            for (int k=0; k<nz; k++)
                for (int l=0; l<nq; l++)
                    Rline[i*ny*nz*nq+j*nz*nq+k*nq+l] = R[i][j][k][l];

    bubbleSort(Rline, nx*ny*nz*nq);

    double quant = Rline[int(p*nx*ny*nz*nq)];

    delete[] Rline;

    return quant;
}

bool converged(std::vector<double***> newArr, std::vector<double***> oldArr, 
        std::vector<double***> reallyoldArr, int n1, int n2, int n3, int n4) {
    double Qthresh = 2.0;
    double Delthresh = 1.1;
    double p = 0.99;

    std::vector<double***> R = create4DArr(n1, n2, n3, n4);
    std::vector<double***> Rold = create4DArr(n1, n2, n3, n4);

    for (int i=0; i<n1; i++) {
        for (int j=0; j<n2; j++) {
            for (int k=0; k<n3; k++) {
                for (int l=0; l<n4; l++) {
                    R[i][j][k][l] = delta(oldArr[i][j][k][l], 
                            newArr[i][j][k][l]);
                    Rold[i][j][k][l] = delta(reallyoldArr[i][j][k][l],
                            newArr[i][j][k][l]);
                }
            }
        }
    }

    double Q = quantile(R,p,n1,n2,n3,n4);
    double Qold = quantile(Rold,p,n1,n2,n3,n4);
    printf("%f   %f\n", Q, Qold);

    double Del = delta(Qold,Q);
    printf("%f\n", Del);

    bool conv = ((Q < Qthresh) && (Del < Delthresh));

    delete4DArr(R, n1, n2, n3, n4);
    delete4DArr(Rold, n1, n2, n3, n4);

    return conv;
}
