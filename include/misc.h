#ifndef MISC_H
#define MISC_H

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
#include <vector>

const double pi = 3.14159265;
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

double random_number();

/* Calculate the blackbody function for a given frequency and temperature. */

double planck_function(double nu, double T);

double planck_function_derivative(double nu, double T);

/* Integrate an array y over array x using the trapezoidal method. */

double integrate(double *y, double *x, int nx);

double* cumulative_integrate(double *y, double *x, int nx);

double* derivative(double *y, double *x, int nx);

double **derivative2D_ax0(double **y, double *x, int nx, int ny);

/* Define what amounts to a tiny value. */

const double EPSILON = 1.0e-6;

/* Test whether two values are equal within a given tolerance. */

bool equal(double x, double y, double tol);

bool equal_zero(double x, double tol);

/* Find the cell in an array in which the given value is located using a 
   tree. */

int find_in_arr(double val, double *arr, int n);

/* Find the cell in an array in which the given value is located. This 
 * function overloads the previous one by allowing you to search a smaller 
 * portion of the array. */

int find_in_arr(double val, double *arr, int lmin, int lmax);

/* Find the cell in an array in which the given value is located, but in this
 * case the array is cyclic. */

int find_in_periodic_arr(double val, double *arr, int n, int lmin, int lmax);

void swap (double* x, double* y);

void bubbleSort(double arr [], int size);

/* Create an empty 2-dimensional array. */

double **create2DArr(int nx, int ny);

void delete2DArr(double **arr, int nx, int ny);

/* Create an empty 3-dimensional array. */

double ***create3DArr(int nx, int ny, int nz);

/* Create a 3-dimensional array filled with a particular value. */

double ***create3DArrValue(int nx, int ny, int nz, int value);

void delete3DArr(double ***arr, int nx, int ny, int nz);

/* Set the value of a 3-dimensional array to a constant value. */

void set3DArrValue(double ***arr, double value, int nx, int ny, int nz);

/* Set one 3-dimensional array equal to another 3-dimensional array, element
 * by element. */

void equate3DArrs(double ***arr1, double ***arr2, int nx, int ny, int nz);

/* Create an empty 4-dimensional array. */

std::vector<double***> create4DArr(int nx, int ny, int nz, int nq);

void delete4DArr(std::vector<double***> arr, int nx, int ny, int nz, int nq);

/* Create a 3-dimensional array filled with a particular value. */

//double ****create4DArrValue(int nx, int ny, int nz, int nq, double value);

/* Set the value of a 4-dimensional array to a constant value. */

void set4DArrValue(std::vector<double***> arr, double value, int nx, int ny, 
        int nz, int nq);

void set4DArrValue(double**** arr, double value, int nx, int ny, 
        int nz, int nq);

/* Set one 4-dimensional array equal to another 4-dimensional array, element
 * by element. */

void equate4DArrs(std::vector<double***> arr1, std::vector<double***> arr2, 
        int nx, int ny, int nz, int nq);

double delta(double x1, double x2);

double quantile(std::vector<double***> R, double p, int nx, int ny, int nz, 
        int nq);

bool converged(std::vector<double***> newArr, std::vector<double***> oldArr, 
        std::vector<double***> reallyoldArr, int n1, int n2, int n3, int n4);

#endif
