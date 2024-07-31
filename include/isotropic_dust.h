#ifndef ISOTROPIC_DUST_H
#define ISOTROPIC_DUST_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <stdlib.h>
#include <cmath>
#include "misc.h"
#include "photon.h"
#include "dust.h"

namespace py = pybind11;

struct IsotropicDust : public Dust {
    IsotropicDust(py::array_t<double> lam, py::array_t<double> kabs, 
            py::array_t<double> ksca) : Dust(lam, kabs, ksca) {};

    /*IsotropicDust(int _nlam, double *_nu, double *_lam, double *_kabs, 
            double *_ksca, double *_kext, double *_albedo) : 
        Dust(_nlam, _nu, _lam, _kabs, _ksca, _kext, _albedo) {};*/
        
    void scatter(Photon *P);
};

#endif
