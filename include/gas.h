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
    Kokkos::View<int*> levels{"levels", 0};
    Kokkos::View<double*> energies{"energies", 0};
    Kokkos::View<double*> weights{"weights", 0};
    Kokkos::View<int*> J{"J", 0};

    py::array_t<int> _levels;
    py::array_t<double> _energies;
    py::array_t<double> _weights;
    py::array_t<int> _J;

    int ntransitions;
    Kokkos::View<int*> transitions{"transitions", 0};
    Kokkos::View<int*> up{"up", 0};
    Kokkos::View<int*> low{"low", 0};
    Kokkos::View<double*> A{"A", 0};
    Kokkos::View<double*> nu{"nu", 0};
    Kokkos::View<double*> Eu{"Eu", 0};

    py::array_t<int> _transitions;
    py::array_t<int> _up;
    py::array_t<int> _low;
    py::array_t<double> _A;
    py::array_t<double> _nu;
    py::array_t<double> _Eu;

    int ntemp;
    Kokkos::View<double*> temp{"temp", 0};
    Kokkos::View<double*> Z{"Z", 0};
    Kokkos::View<double*> dZdT{"dZdT", 0};

    Gas(double _mu, py::array_t<int> __levels, py::array_t<double> __energies, 
            py::array_t<double> __weights, py::array_t<int> __J, 
            py::array_t<int> __transitions, py::array_t<int> __up, 
            py::array_t<int> __low, py::array_t<double> __A, 
            py::array_t<double> __nu, py::array_t<double> __Eu);

    ~Gas();

    double partition_function(double T);
};

#endif
