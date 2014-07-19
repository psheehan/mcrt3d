#ifndef PHOTON_CC
#define PHOTON_CC

#include "vector.cc"

struct Photon {
    double energy;
    double nu;
    Vector<double, 3> r;
    Vector<double, 3> n;
    Vector<double, 3> invn;
    double *current_kext, *current_albedo;
    Vector<int, 3> l;

    double rad, phi, theta;

    /* Move the photon a distance s along its direction vector. */

    void move(double s) {
        r += s*n;
    }

    /* Clean up the photon to remove any pointers. */

    void clean() {
        delete current_kext;
        delete current_albedo;
    }
};

struct Ray : public Photon {
    double tau;
    double intensity;
    double pixel_size;
    bool pixel_too_large;
};

#endif
