#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdio.h>
#include <cmath>

#include "params.h"
#include "dust.h"
#include "isotropic_dust.h"
#include "grid.h"
#include "cartesian_grid.h"
#include "cylindrical_grid.h"
#include "spherical_grid.h"
#include "source.h"
#include "star.h"
#include "camera.h"
#include "misc.h"
#include "photon.h"

namespace py = pybind11;

struct MCRT {
    Grid *G;
    Params *Q;
    Camera *C;

    py::list images;
    py::list spectra;

    //MCRT(Grid *G, Params *Q);
    MCRT();
    ~MCRT();

    void set_cartesian_grid(py::array_t<double> x, py::array_t<double> y,
            py::array_t<double> z);
    void set_cylindrical_grid(py::array_t<double> r, py::array_t<double> phi,
            py::array_t<double> z);
    void set_spherical_grid(py::array_t<double> r, py::array_t<double> theta,
            py::array_t<double> phi);

    void thermal_mc(int nphot, bool bw, bool use_mrw, double mrw_gamma,
            bool verbose, int nthreads);
    void scattering_mc(py::array_t<double> scatt, int nphot, bool verbose, 
            bool save, int nthreads);

    void mc_iteration(int nthreads);

    void run_image(py::array_t<double> lam, int nx, int ny, double pixel_size, 
            int nphot, double incl, double pa, double dpc, int nthreads);
    void run_unstructured_image(py::array_t<double> lam, int nx, int ny, 
            double pixel_size, int nphot, double incl, double pa, double dpc, 
            int nthreads);
    void run_spectrum(py::array_t<double> lam, int nphot, double incl, 
            double pa, double dpc, int nthreads);
};
