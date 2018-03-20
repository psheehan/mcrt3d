#ifndef SOURCE_H
#define SOURCE_H

#include <stdlib.h>
#include "vector.h"
#include "dust.h"
#include "photon.h"
#include "misc.h"
#include "params.h"

struct Source {
    Vector<double, 3> r;
    double *nu;
    double *Bnu;
    int nnu;

    virtual Photon *emit(int nphot);
    virtual Photon *emit(double _nu, double _dnu, int nphot);
    virtual double intercept_distance(Photon *P);
};

#endif
