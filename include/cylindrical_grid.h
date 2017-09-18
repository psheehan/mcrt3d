#ifndef CYLINDRICAL_GRID_H
#define CYLINDRICAL_GRID_H

#include "grid.h"
#include "vector.h"
#include "photon.h"

struct CylindricalGrid : public Grid {
    double next_wall_distance(Photon *P);
    double outer_wall_distance(Photon *P);
    double minimum_wall_distance(Photon *P);
    Vector<int, 3> photon_loc(Photon *P);
    bool in_grid(Photon *P);
};

#endif
