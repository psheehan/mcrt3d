#ifndef SPHERICAL_GRID_H
#define SPHERICAL_GRID_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "grid.h"
#include "vector.h"
#include "photon.h"

struct SphericalGrid : public Grid {
    SphericalGrid(py::array_t<double> r, py::array_t<double> theta,
            py::array_t<double> phi);

    /*SphericalGrid(int _n1, int _n2, int _n3, int _nw1, int _nw2, int _nw3, 
            double *_w1, double *_w2, double *_w3, double *_volume,
            bool _mirror_symmetry) : 
        Grid(_n1, _n2, _n3, _nw1, _nw2, _nw3, _w1, _w2, _w3, _volume, 
                _mirror_symmetry) {};*/

    ~SphericalGrid();

    py::array_t<double> r;
    py::array_t<double> theta;
    py::array_t<double> phi;

    Kokkos::View<double*> sin_w2{"sin_w2", 0};
    Kokkos::View<double*> cos_w2{"cos_w2", 0};
    Kokkos::View<double*> neg_mu{"neg_mu", 0};
    Kokkos::View<double*> tan_w2{"tan_w2", 0};
    Kokkos::View<double*> sin_w3{"sin_w3", 0};
    Kokkos::View<double*> cos_w3{"cos_w3", 0};
    Kokkos::View<double*> sin_tol_w2{"sin_tol_w2", 0};
    Kokkos::View<double*> cos_tol_w2{"cos_tol_w2", 0};

    double next_wall_distance(Photon *P);
    double outer_wall_distance(Photon *P);
    double minimum_wall_distance(Photon *P);
    double smallest_wall_size(Photon *P);
    double smallest_wall_size(Ray *R);
    double grid_size();
    Vector<int, 3> photon_loc(Photon *P);
    void photon_loc_mrw(Photon *P);
    Vector<double, 3> random_location_in_cell(int ix, int iy, int iz);
    bool in_grid(Photon *P);
};

#endif
