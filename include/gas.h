#ifndef GAS_H
#define GAS_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <stdlib.h>
#include <cmath>
#include "misc.h"

namespace py = pybind11;

struct Gas {
    double mu;

    int nlevels;
    int *levels;
    double *energies;
    double *weights;
    int *J;

    py::array_t<int> _levels;
    py::array_t<double> _energies;
    py::array_t<double> _weights;
    py::array_t<int> _J;

    int ntransitions;
    int *transitions;
    int *up;
    int *low;
    double *A;
    double *nu;
    double *Eu;

    py::array_t<int> _transitions;
    py::array_t<int> _up;
    py::array_t<int> _low;
    py::array_t<double> _A;
    py::array_t<double> _nu;
    py::array_t<double> _Eu;

    int ntemp;
    double *temp;
    double *Z;
    double *dZdT;

    Gas(double _mu, py::array_t<int> __levels, py::array_t<double> __energies, 
            py::array_t<double> __weights, py::array_t<int> __J, 
            py::array_t<int> __transitions, py::array_t<int> __up, 
            py::array_t<int> __low, py::array_t<double> __A, 
            py::array_t<double> __nu, py::array_t<double> __Eu);

    ~Gas();

    double partition_function(double T);
};

#endif
