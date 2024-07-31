#ifndef STAR_H
#define STAR_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <Kokkos_Random.hpp>

#include <stdlib.h>
#include "vector.h"
#include "dust.h"
#include "photon.h"
#include "misc.h"
#include "params.h"
#include "source.h"

namespace py = pybind11;

struct Star : public Source {
    double mass;
    double radius;
    double temperature;
    double luminosity;

    Kokkos::Random_XorShift64_Pool<> *random_pool;

    Star(double x, double y, double z, double _mass, double _radius, 
            double _temperature);

    ~Star();

    void set_blackbody_spectrum(py::array_t<double> nu);

    Kokkos::View<double*> random_nu_CPD{"random_nu_CPD", 0};

    Photon *emit(int nphot);
    Photon *emit(double _nu, double _dnu, int nphot);
    Ray *emit_ray(Kokkos::View<double*> _nu, int _nnu, double _pixelsize, \
            Vector<double, 3> _n, int nphot);
    Ray *emit_ray(Kokkos::View<double*> _nu, int nnu, Vector<double, 3> _n, int nphot);
    double intercept_distance(Photon *P);
    double random_nu();
    double flux(double freq);

    void reemit(Photon *P);
};

#endif
