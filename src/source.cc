#include "source.h"

/* Emit a photon from the source. */

Photon *Source::emit(int nphot) {
    return new Photon();
}

Photon *Source::emit(double _nu, double _dnu, int nphot) {
    return new Photon();
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
