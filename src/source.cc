#include "source.h"

Source::~Source() {
}

/* Emit a photon from the source. */

Photon *Source::emit(int nphot) {
    return new Photon();
}

Photon *Source::emit(double _nu, double _dnu, int nphot) {
    return new Photon();
}

Ray *Source::emit_ray(Kokkos::View<double*> _nu, int _nnu, double _pixelsize, \
        Vector<double, 3> _n, int nphot) {
    return new Ray();
}

Ray *Source::emit_ray(Kokkos::View<double*> _nu, int _nnu, Vector<double, 3> _n, 
        int nphot) {
    return new Ray();
}

double Source::intercept_distance(Photon *P) {
    return 0;
}

double Source::flux(double freq) {
    return 0;
}

double Source::random_nu() {
    return 0;
}

void Source::reemit(Photon *P) {
}
