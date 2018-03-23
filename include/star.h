#ifndef STAR_H
#define STAR_H

#include <stdlib.h>
#include "vector.h"
#include "dust.h"
#include "photon.h"
#include "misc.h"
#include "params.h"
#include "source.h"

struct Star : public Source {
    double mass;
    double radius;
    double temperature;
    double luminosity;

    Star(double x, double y, double z, double _mass, double _radius, 
            double _temperature);

    void set_blackbody_spectrum(int _nnu, double *_nu, double *_Bnu, 
            double _luminosity, double *_random_nu_CPD);

    double *random_nu_CPD;

    Photon *emit(int nphot);
    Photon *emit(double _nu, double _dnu, int nphot);
    Ray *emit_ray(double _nu, double _dnu, double _pixelsize, \
            Vector<double, 3> _n, int nphot);
    double intercept_distance(Photon *P);
    double random_nu();
    double flux(double freq);
};

#endif
