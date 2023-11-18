#ifndef DUST_H
#define DUST_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <stdlib.h>
#include <cmath>
#include "misc.h"
#include "photon.h"

namespace py = pybind11;

struct Dust {
    int nlam;
    double *nu;
    double *lam;
    double *kabs;
    double *ksca;
    double *kext;
    double *albedo;
    double *dkextdnu;
    double *dalbedodnu;

    py::array_t<double> _nu;
    py::array_t<double> _lam;
    py::array_t<double> _kabs;
    py::array_t<double> _ksca;
    py::array_t<double> _kext;
    py::array_t<double> _albedo;

    int ntemp;
    double *temp;
    double *planck_opacity;
    double *dplanck_opacity_dT;
    double *rosseland_extinction;
    double *drosseland_extinction_dT;
    double **random_nu_CPD;
    double **random_nu_CPD_bw;
    double **drandom_nu_CPD_dT;
    double **drandom_nu_CPD_bw_dT;

    Dust(py::array_t<double> lam, py::array_t<double> kabs, 
            py::array_t<double> ksca);

    Dust(int _nlam, double *_nu, double *_lam, double *_kabs, double *_ksca, 
            double *_kext, double *_albedo);

    ~Dust();
        
    void set_lookup_tables();

    virtual void scatter(Photon *P);
    void absorb(Photon *P, double T, bool bw);
    void absorb_mrw(Photon *P, double T, bool bw);

    double random_nu(double T, bool bw);
    double opacity(double freq);
    double albdo(double freq);
    double planck_mean_opacity(double T);
    double rosseland_mean_extinction(double T);
};

#endif
