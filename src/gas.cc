#include "gas.h"

/* Functions to set up the dust. */

Gas::Gas(double _mu, py::array_t<int> __levels, py::array_t<double> __energies, 
        py::array_t<double> __weights, py::array_t<int> __J, 
        py::array_t<int> __transitions, py::array_t<int> __up, 
        py::array_t<int> __low, py::array_t<double> __A, 
        py::array_t<double> __nu, py::array_t<double> __Eu) {

    mu = _mu;

    /*_levels = __levels; _energies = __energies; _weights = __weights; _J = __J;
    _transitions = __transitions; _up = __up; _low = __low; _A = __A;
    _nu = __nu; _Eu = __Eu;

    // Load the array buffers to get the proper setup info.*/

    auto _levels_buf = __levels.request(); 
    /*auto _energies_buf = __energies.request();
    auto _weights_buf = __weights.request();
    auto _J_buf = __J.request();*/

    auto _transitions_buf = __transitions.request(); 
    /*auto _up_buf = __up.request();
    auto _low_buf = __low.request();
    auto _A_buf = __A.request();
    auto _nu_buf = __nu.request();
    auto _Eu_buf = __Eu.request();

    if (_levels_buf.ndim != 1 || _energies_buf.ndim != 1 || 
            _weights_buf.ndim != 1 || _J_buf.ndim != 1 || 
            _transitions_buf.ndim != 1 || _up_buf.ndim != 1 || 
            _low_buf.ndim != 1 || _A_buf.ndim != 1 || _nu_buf.ndim != 1 
            || _Eu_buf.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one");

    // Now get the correct format.*/

    nlevels = _levels_buf.shape[0];
    /*levels = (int *) _levels_buf.ptr;
    energies = (double *) _energies_buf.ptr;
    weights = (double *) _weights_buf.ptr;
    J = (int *) _J_buf.ptr;*/

    ntransitions = _transitions_buf.shape[0];
    /*transitions = (int *) _transitions_buf.ptr;
    up = (int *) _up_buf.ptr;
    low = (int *) _low_buf.ptr;
    A = (double *) _A_buf.ptr;
    nu = (double *) _nu_buf.ptr;
    Eu = (double *) _Eu_buf.ptr;*/

    Kokkos::resize(levels, nlevels);
    Kokkos::deep_copy(levels, view_from_array(__levels));
    Kokkos::resize(energies, nlevels);
    Kokkos::deep_copy(energies, view_from_array(__energies));
    Kokkos::resize(weights, nlevels);
    Kokkos::deep_copy(weights, view_from_array(__weights));
    Kokkos::resize(J, nlevels);
    Kokkos::deep_copy(J, view_from_array(__J));

    _levels = array_from_view(levels, 1, {(size_t) nlevels});
    _energies = array_from_view(energies, 1, {(size_t) nlevels});
    _weights = array_from_view(weights, 1, {(size_t) nlevels});
    _J = array_from_view(J, 1, {(size_t) nlevels});

    Kokkos::resize(transitions, ntransitions);
    Kokkos::deep_copy(transitions, view_from_array(__transitions));
    Kokkos::resize(up, ntransitions);
    Kokkos::deep_copy(up, view_from_array(__up));
    Kokkos::resize(low, ntransitions);
    Kokkos::deep_copy(low, view_from_array(__low));
    Kokkos::resize(A, ntransitions);
    Kokkos::deep_copy(A, view_from_array(__A));
    Kokkos::resize(nu, ntransitions);
    Kokkos::deep_copy(nu, view_from_array(__nu));
    Kokkos::resize(Eu, ntransitions);
    Kokkos::deep_copy(Eu, view_from_array(__Eu));

    _transitions = array_from_view(transitions, 1, {(size_t) ntransitions});
    _up = array_from_view(up, 1, {(size_t) ntransitions});
    _low = array_from_view(low, 1, {(size_t) ntransitions});
    _A = array_from_view(A, 1, {(size_t) ntransitions});
    _nu = array_from_view(nu, 1, {(size_t) ntransitions});
    _Eu = array_from_view(Eu, 1, {(size_t) ntransitions});

    // Calculate the partition function.

    ntemp = 1000;
    Kokkos::resize(temp, ntemp);
    Kokkos::resize(Z, ntemp);

    for (int i = 0; i < ntemp; i++) {
        temp(i) = pow(10.,-1.+i*6./(ntemp-1));
        Z(i) = 0;

        for (int j = 0; j < nlevels; j++)
            Z(i) += weights(j)*exp(-h_p*c_l*energies(j) / (k_B * temp(i)));
    }

    Kokkos::resize(dZdT, ntemp-1);
    Kokkos::deep_copy(dZdT, derivative(Z, temp, ntemp));
}

Gas::~Gas() {
}

double Gas::partition_function(double T) {
    int n = find_in_arr(T,temp,ntemp);

    double partition_function = dZdT(n)*(T-temp(n))+Z(n);

    return partition_function;
}
