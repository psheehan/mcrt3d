#ifndef PHOTON_H
#define PHOTON_H

#include "vector.h"

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

    void move(double s);

    /* Clean up the photon to remove any pointers. */

    void clean();
};

struct Ray : public Photon {
    double tau;
    double intensity;
    double pixel_size;
    bool pixel_too_large;
};

#endif
