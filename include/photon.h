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

    double event_count;
    double same_cell_count;

    double rad, phi, theta;

    /* Clean up the photon to remove any pointers. */

    ~Photon();

    /* Move the photon a distance s along its direction vector. */

    void move(double s);
};

struct Ray : public Photon {
    double tau;
    double intensity;
    double pixel_size;
    bool pixel_too_large;
};

#endif
