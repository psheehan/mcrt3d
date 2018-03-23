#ifndef GRID_H
#define GRID_H

#include <cmath>
#include <vector>
#include "pymangle.h"
#include "vector.h"
#include "dust.h"
#include "isotropic_dust.h"
#include "source.h"
#include "photon.h"
#include "misc.h"
#include "params.h"

struct Grid {
    int n1;
    int n2;
    int n3;
    int nw1;
    int nw2;
    int nw3;
    double *w1;
    double *w2;
    double *w3;

    std::vector<double***> dens;
    std::vector<double***> energy;
    std::vector<double***> temp;
    std::vector<double***> mass;
    std::vector<double****> scatt;
    double ***volume;

    int nspecies;
    std::vector<Dust*> dust;

    int nsources;
    std::vector<Source*> sources;
    double total_lum;

    Params *Q;

    int ny;
    double *y;
    double *f;
    double *dydf;

    Grid(int _n1, int _n2, int _n3, int _nw1, int _nw2, int _nw3, 
            double *_w1, double *_w2, double *_w3, double *_volume);

    void add_density(double *_dens, double *_temp, double *_mass, 
            Dust *D);
    void add_source(Source *S);
    void set_mrw_tables(double *y, double *f, double *dydf, int ny);
    void initialize_scattering_array();
    void deallocate_scattering_array();

    Photon *emit(int iphot);

    virtual double next_wall_distance(Photon *P);
    virtual double outer_wall_distance(Photon *P);
    virtual double minimum_wall_distance(Photon *P);

    void propagate_photon_full(Photon *P);
    void propagate_photon(Photon *P, double tau, bool absorb);
    void propagate_photon_mrw(Photon *P);
    void propagate_ray(Ray *R);
    void propagate_ray_from_source(Ray *R);

    void absorb(Photon *P, int idust);
    void scatter(Photon *P, int idust);

    virtual Vector<int, 3> photon_loc(Photon *P);
    virtual bool in_grid(Photon *P);

    void update_grid(Vector<int, 3> l);
    void update_grid();

    double cell_lum(Vector<int, 3> l);
};

#endif
