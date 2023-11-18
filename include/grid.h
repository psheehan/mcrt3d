#ifndef GRID_H
#define GRID_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <signal.h>

#include <cmath>
#include <vector>
#include <algorithm>
#include "pymangle.h"
#include "vector.h"
#include "dust.h"
#include "isotropic_dust.h"
#include "source.h"
#include "star.h"
#include "photon.h"
#include "misc.h"
#include "params.h"
#include "gas.h"

namespace py = pybind11;

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

    bool mirror_symmetry;

    py::array_t<double> _w1;
    py::array_t<double> _w2;
    py::array_t<double> _w3;
    py::array_t<double> _volume;

    py::list _dens;
    py::list _temp;
    py::list _scatt;

    std::vector<double*> dens;
    std::vector<std::vector<double*>> energy;
    std::vector<double*> temp;
    std::vector<double*> mass;
    std::vector<double*> rosseland_mean_extinction;
    std::vector<double*> planck_mean_opacity;
    std::vector<double*> luminosity;
    std::vector<double*> scatt;
    double *volume;
    double *uses_mrw;

    std::vector<double*> number_dens;
    std::vector<double*> gas_temp;
    std::vector<double*> velocity;
    std::vector<double*> microturbulence;

    py::list _number_dens;
    py::list _gas_temp;
    py::list _microturbulence;
    py::list _velocity;

    #ifdef _OPENMP
    std::vector<omp_lock_t*> lock;
    #endif

    int nspecies;
    std::vector<Dust*> dust;
    py::list _dust;

    int ngases;
    std::vector<Gas *> gas;
    py::list _gas;
    std::vector<int> include_lines;
    std::vector<double*> level_populations;
    std::vector<double*> alpha_line;
    std::vector<double*> inv_gamma_thermal;


    int nsources;
    std::vector<Source*> sources;
    py::list _sources;
    double total_lum;

    Params *Q;

    Grid(py::array_t<double> w1, py::array_t<double> w2, 
            py::array_t<double> w3);

    Grid(int _n1, int _n2, int _n3, int _nw1, int _nw2, int _nw3, 
            double *_w1, double *_w2, double *_w3, double *_volume,
            bool _mirror_symmetry);

    ~Grid();

    //void add_density(double *_dens, double *_temp, double *_mass, 
    //        Dust *D);
    void add_density(py::array_t<double>, Dust *d);
    void add_number_density(py::array_t<double>, py::array_t<double>, 
            py::array_t<double>, Gas *g);

    //void add_source(Source *S);
    void add_star(double x, double y, double z, double _mass, double _radius, 
            double _temperature);

    void add_scattering_array(py::array_t<double> _scatt, int nthreads);
    void initialize_scattering_array(int nthreads);
    void deallocate_scattering_array(int start);
    void collapse_scattering_array();

    void add_energy_arrays(int nthreads);
    void deallocate_energy_arrays();

    void initialize_luminosity_array();
    void initialize_luminosity_array(double nu);
    void deallocate_luminosity_array();

    Photon *emit(int iphot);
    Photon *emit(double _nu, double _dnu, int photons_per_source);

    virtual Vector<double, 3> random_location_in_cell(int ix, int iy, int iz);

    virtual double next_wall_distance(Photon *P);
    virtual double outer_wall_distance(Photon *P);
    virtual double minimum_wall_distance(Photon *P);
    virtual double smallest_wall_size(Photon *P);
    virtual double smallest_wall_size(Ray *R);
    virtual double grid_size();

    void propagate_photon_full(Photon *P);
    void propagate_photon(Photon *P, double tau, bool absorb);
    void propagate_photon_scattering(Photon *P);
    void propagate_photon_mrw(Photon *P);
    void propagate_ray(Ray *R);
    void propagate_ray_from_source(Ray *R);

    void absorb(Photon *P, int idust);
    void absorb_mrw(Photon *P, int idust);
    void scatter(Photon *P, int idust);

    void random_dir_mrw(Photon *P);

    virtual Vector<int, 3> photon_loc(Photon *P);
    virtual void photon_loc_mrw(Photon *P);
    virtual bool in_grid(Photon *P);
    virtual bool on_and_parallel_to_wall(Photon *P);

    void update_grid(Vector<int, 3> l, int cell_index);
    void update_grid();

    double cell_lum(int idust, int cell_index);
    double cell_lum(int idust, int cell_index, double nu);

    Vector<double, 3> vector_velocity(int igas, Photon *P);
    double maximum_velocity(int igas);
    double maximum_gas_temperature(int igas);
    double maximum_microturbulence(int igas);
    double line_profile(int igas, int iline, int itrans, int cell_index, 
            double nu);
    void calculate_level_populations(int igas, int iline);
    void set_tgas_eq_tdust();
    void select_lines(py::array_t<double> lam);
    void deselect_lines();
};

#endif
