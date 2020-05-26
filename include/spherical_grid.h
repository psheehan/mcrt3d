#ifndef SPHERICAL_GRID_H
#define SPHERICAL_GRID_H

#include "grid.h"
#include "vector.h"
#include "photon.h"

struct SphericalGrid : public Grid {
    SphericalGrid(int _n1, int _n2, int _n3, int _nw1, int _nw2, int _nw3, 
            double *_w1, double *_w2, double *_w3, double *_volume) : 
        Grid(_n1, _n2, _n3, _nw1, _nw2, _nw3, _w1, _w2, _w3, _volume) {};

    double next_wall_distance(Photon *P);
    double outer_wall_distance(Photon *P);
    double minimum_wall_distance(Photon *P);
    Vector<int, 3> photon_loc(Photon *P);
    Vector<double, 3> random_location_in_cell(int ix, int iy, int iz);
    bool in_grid(Photon *P);
};

#endif
