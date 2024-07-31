#ifndef PHOTON_H
#define PHOTON_H

#include "vector.h"
#include "misc.h"

struct Photon {
    double energy;
    double nu;
    Vector<double, 3> r;
    Vector<double, 3> n;
    Vector<double, 3> invn;
    Vector<double, 3> nframe;
    Kokkos::View<double*> current_kext{"current_kext", 0}, current_albedo{"current_albedo", 0};
    Vector<int, 3> l;
    int cell_index;

    double event_count;
    double same_cell_count;

    double rad, phi, theta;
    double sin_theta, cos_theta;
    double sin_phi, cos_phi;

    int ithread;

    /* Clean up the photon to remove any pointers. */

    ~Photon();

    /* Move the photon a distance s along its direction vector. */

    void move(double s);
};

struct Ray : public Photon {
    Kokkos::View<double**> current_kext{"current_kext", 0, 0}, current_albedo{"current_albedo", 0, 0};

    int nnu;
    int ndust;

    Kokkos::View<double*> nu{"nu", 0};
    Kokkos::View<double*> tau{"tau", 0};
    Kokkos::View<double*> intensity{"intensity", 0};
    double pixel_size;
    bool pixel_too_large;

    ~Ray();
};

#endif
