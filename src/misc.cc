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

double random_number(Kokkos::Random_XorShift64_Pool<> *random_pool) {
    auto generator = random_pool->get_state();
    double x = generator.drand(0., 1.);
    random_pool->free_state(generator);
    return x;
};

Kokkos::View<double*> view_from_array(py::array_t<double> arr) {
    auto arr_buf = arr.request();
    int n = arr_buf.shape[0];

    Kokkos::View<double*> result("result", n);
    auto arr_access = arr.unchecked<1>();
    for (size_t i = 0; i < n; i++) result(i) = arr_access[i];

    return result;
}

Kokkos::View<int*> view_from_array(py::array_t<int> arr) {
    auto arr_buf = arr.request();
    int n = arr_buf.shape[0];

    Kokkos::View<int*> result("result", n);
    auto arr_access = arr.unchecked<1>();
    for (size_t i = 0; i < n; i++) result(i) = arr_access[i];

    return result;
}

template<typename Ta, typename Tv>
py::array_t<Ta> array_from_view(Kokkos::View<Tv> v, int ndim, std::vector<size_t> extents) {
    std::vector<size_t> strides;
    for (int i=0; i < ndim; i++) {
        strides.push_back(sizeof(Ta));
        for (int j = i+1; j<ndim; j++)
            strides[i] *= extents[j];
    }
    py::array_t<Ta> arr = py::array_t<Ta>(py::buffer_info(
            v.data(),                               /* Pointer to buffer */
            sizeof(Ta),                          /* Size of one scalar */
            py::format_descriptor<Ta>::format(), /* Python struct-style format descriptor */
            ndim,                                      /* Number of dimensions */
            extents,                                /* Buffer dimensions */
            strides
        )
    );

    return arr;
};

py::array_t<double> array_from_view(Kokkos::View<double*> v, int ndim, std::vector<size_t> extents) {
    std::vector<size_t> strides;
    for (int i=0; i < ndim; i++) {
        strides.push_back(sizeof(double));
        for (int j = i+1; j<ndim; j++)
            strides[i] *= extents[j];
    }
    py::array_t<double> arr = py::array_t<double>(py::buffer_info(
            v.data(),                               /* Pointer to buffer */
            sizeof(double),                          /* Size of one scalar */
            py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
            ndim,                                      /* Number of dimensions */
            extents,                                /* Buffer dimensions */
            strides
        )
    );

    return arr;
};

py::array_t<double> array_from_view(Kokkos::View<double**> v, int ndim, std::vector<size_t> extents) {
    std::vector<size_t> strides;
    for (int i=0; i < ndim; i++) {
        strides.push_back(sizeof(double));
        for (int j = i+1; j<ndim; j++)
            strides[i] *= extents[j];
    }
    py::array_t<double> arr = py::array_t<double>(py::buffer_info(
            v.data(),                               /* Pointer to buffer */
            sizeof(double),                          /* Size of one scalar */
            py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
            ndim,                                      /* Number of dimensions */
            extents,                 /* Buffer dimensions */
            strides
        )
    );

    return arr;
};

py::array_t<int> array_from_view(Kokkos::View<int*> v, int ndim, std::vector<size_t> extents) {
    std::vector<size_t> strides;
    for (int i=0; i < ndim; i++) {
        strides.push_back(sizeof(int));
        for (int j = i+1; j<ndim; j++)
            strides[i] *= extents[j];
    }
    py::array_t<int> arr = py::array_t<int>(py::buffer_info(
            v.data(),                               /* Pointer to buffer */
            sizeof(int),                          /* Size of one scalar */
            py::format_descriptor<int>::format(), /* Python struct-style format descriptor */
            ndim,                                      /* Number of dimensions */
            extents,                                /* Buffer dimensions */
            strides
        )
    );

    return arr;
};

/* Calculate the blackbody function for a given frequency and temperature. */

double planck_function(double nu, double T) {
    double value = 2.0*h*nu*nu*nu/(c_l*c_l)*1.0/(exp(h*nu/(k_B*T))-1.0);

    if (std::isnan(value))
        return 0;
    else
        return value;
};

double planck_function_derivative(double nu, double T) {
    double value = (-2.0*h*nu*nu*nu*nu)/(c_l*c_l*k_B*T*T) / 
        (exp(h*nu/(k_B*T))-1.0) / (1. - exp(-h*nu/(k_B*T)));

    if (std::isnan(value))
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

double integrate(Kokkos::View<double*> y, Kokkos::View<double*> x, int nx) {
    double sum = 0.0;

    for (int i=0; i<nx-1; i++)
        sum += 0.5*(y(i+1)+y(i))*(x(i+1)-x(i));

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

Kokkos::View<double*> cumulative_integrate(Kokkos::View<double*> y, Kokkos::View<double*> x, int nx) {
    double sum = 0;
    Kokkos::View<double*> cum_sum("cum_sum", nx);

    cum_sum(0) = sum;
    for (int i = 1; i < nx; i++) {
        sum += 0.5*(y(i)+y(i-1))*(x(i)-x(i-1));
        cum_sum(i) = sum;
    }

    for (int i = 0; i < nx; i++)
        cum_sum(i) /= sum;

    return cum_sum;
}

/* Take the derivative of an array. */

double *derivative(double *y, double *x, int nx) {
    double *result = new double[nx-1];

    for (int i = 0; i < nx-1; i++)
        result[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]);

    return result;
}

Kokkos::View<double*> derivative(Kokkos::View<double*> y, Kokkos::View<double*> x, int nx) {
    Kokkos::View<double*> result("result", nx-1);

    for (int i = 0; i < nx-1; i++)
        result(i) = (y(i+1) - y(i)) / (x(i+1) - x(i));

    return result;
}

double **derivative2D_ax0(double **y, double *x, int nx, int ny) {
    double **result = create2DArr(nx-1, ny);

    for (int i = 0; i < nx-1; i++)
        for (int j = 0; j < ny; j++)
            result[i][j] = (y[i+1][j] - y[i][j]) / (x[i+1] - x[i]);

    return result;
}

Kokkos::View<double**> derivative2D_ax0(Kokkos::View<double**> y, Kokkos::View<double*> x, int nx, int ny) {
    Kokkos::View<double**> result("result", nx-1, ny);

    for (int i = 0; i < nx-1; i++)
        for (int j = 0; j < ny; j++)
            result(i,j) = (y(i+1,j) - y(i,j)) / (x(i+1) - x(i));

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

int find_in_arr(double val, Kokkos::View<double*> arr, int n) {
    int lmin = 0;
    int lmax = n-1;
    bool not_found = true;
    int l;

    while (not_found) {
        int ltest = (lmax-lmin)/2+lmin;

        if ((val >= arr(ltest)) && (val <= arr(ltest+1))) {
            l = ltest;
            not_found = false;
        }
        else {
            if (val < arr(ltest))
                lmax = ltest;
            else
                lmin = ltest;
        }
    }

    return l;
};

int find_in_arr(double val, double* arr, int n) {
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

int find_in_arr(double val, Kokkos::View<double*> arr, int lmin, int lmax) {
    int l;

    for (int i=lmin; i <= lmax; i++) {
        if ((val >= arr(i)) && (val <= arr(i+1)))
            l = i;
    }

    return l;
};

int find_in_arr(double val, double* arr, int lmin, int lmax) {
    int l;

    for (int i=lmin; i <= lmax; i++) {
        if ((val >= arr[i]) && (val <= arr[i+1]))
            l = i;
    }

    return l;
};

/* Find the cell in an array in which the given value is located, but in this
 * case the array is cyclic. */

int find_in_periodic_arr(double val, Kokkos::View<double*> arr, int n, int lmin, int lmax) {
    int l = -1;

    for (int i=lmin; i <= lmax; i++) {
        int index = (i+n)%(n);
        if ((val >= arr(index)) && (val <= arr(index+1))) {
            l = index;
        }
    }

    return l;
}

int find_in_periodic_arr(double val, double* arr, int n, int lmin, int lmax) {
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

std::vector<double*> create2DVecArr(int nx, int ny) {
    std::vector<double*> arr;
    for (int i=0; i<nx; i++)
        arr.push_back(new double[ny]);

    return arr;
};

void delete2DVecArr(std::vector<double*> arr, int nx, int ny) {
    for (int i = 0; i < nx; i++)
        delete[] arr[i];
    arr.clear();
}

void equate2DVecArrs(std::vector<double*> arr1, std::vector<double*> arr2, 
        int nx, int ny) {
    for (int i=0; i<nx; i++)
        for (int j=0; j<ny; j++)
            arr1[i][j] = arr2[i][j];
}

double delta(double x1, double x2) {
    double d1 = x1/x2;
    double d2 = x2/x1;

    double delt = d1;
    if (d2 > d1) delt = d2;

    return delt;
}

double quantile(double* R, double p, int nx, int ny, int nz, int nq) {
    double *Rline = new double[nx*ny*nz*nq];

    for (int i=0; i<nx*ny*nz*nq; i++)
        Rline[i] = R[i];

    bubbleSort(Rline, nx*ny*nz*nq);

    double quant = Rline[int(p*nx*ny*nz*nq)];

    delete[] Rline;

    return quant;
}

bool converged(Kokkos::View<double****> newArr, Kokkos::View<double****> oldArr, 
        Kokkos::View<double****> reallyoldArr, int n1, int n2, int n3, int n4) {
    double Qthresh = 2.0;
    double Delthresh = 1.1;
    double p = 0.99;

    double* R = new double[n1*n2*n3*n4];
    double* Rold = new double[n1*n2*n3*n4];

    for (int i=0; i<n1; i++) {
        for (int j=0; j<n2; j++) {
            for (int k=0; k<n3; k++) {
                for (int l=0; l<n4; l++) {
                    R[i*n2*n3*n4 + j*n3*n4 + k*n4 + j] = delta(oldArr(i,j,k,l), newArr(i,j,k,l));
                    Rold[i*n2*n3*n4 + j*n3*n4 + k*n4 + j] = delta(reallyoldArr(i,j,k,l), newArr(i,j,k,l));
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

    delete[] R; delete[] Rold;

    return conv;
}
