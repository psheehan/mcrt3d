#ifndef SOURCES_CC
#define SOURCES_CC

#include <stdlib.h>
#include "vector.cc"
#include "dust.cc"
#include "photon.cc"
#include "misc.cc"
#include "params.cc"

struct Source {
    Vector<double, 3> r;
    double *nu;
    double *Bnu;
    int nnu;

    virtual Photon *emit(int nphot);
    virtual double intercept_distance(Photon *P);
};

/* Emit a photon from the source. */

Photon *Source::emit(int nphot) {
    return new Photon();
}

double Source::intercept_distance(Photon *P) {
    return 0;
}

#endif
