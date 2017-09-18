#ifndef CARTESIAN_GRID_H
#define CARTESIAN_GRID_H

#include "grid.h"
#include "vector.cc"
#include "photon.h"

struct CartesianGrid : public Grid {
    double next_wall_distance(Photon *P);
    double outer_wall_distance(Photon *P);
    double minimum_wall_distance(Photon *P);
    Vector<int, 3> photon_loc(Photon *P);
    bool in_grid(Photon *P);
};

#endif
