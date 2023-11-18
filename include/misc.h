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

#ifdef _OPENMP
#include <omp.h>
#endif

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
const double k_B = 1.380658e-16;
const double sigma = 5.67051e-5;
const double arcsec = 4.8481e-6;
const double Jy = 1.0e-23;
const double h_p = 6.6261e-27;
const double m_p = 1.6726e-24;

/* Get a random number between 0 and 1. */

static int seed1 = 1, seed2 = 1;
#pragma omp threadprivate(seed1, seed2)

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

std::vector<double*> create2DVecArr(int nx, int ny);
void delete2DVecArr(std::vector<double*> arr, int nx, int ny);

void equate2DVecArrs(std::vector<double*> arr1, std::vector<double*> arr2, 
        int nx, int ny);

/* Copy the values in one 2D array into another array. */

void equate2DArrs(double *arr1, double *arr2, int nx);

/* Functions to assess the convergence of a Lucy iteration. */

double delta(double x1, double x2);

double quantile(std::vector<double***> R, double p, int nx, int ny, int nz, 
        int nq);

bool converged(std::vector<double***> newArr, std::vector<double***> oldArr, 
        std::vector<double***> reallyoldArr, int n1, int n2, int n3, int n4);

#endif
