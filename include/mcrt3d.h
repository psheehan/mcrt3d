#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <stdio.h>
#include <cmath>
#include "params.h"
#include "isotropic_dust.h"
#include "grid.h"
#include "cartesian_grid.h"
#include "cylindrical_grid.h"
#include "spherical_grid.h"
#include "star.h"
#include "pymangle.h"
//#include "timer.c"
#include "camera.h"

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

    void thermal_mc(int nphot, bool bw, bool use_mrw, int mrw_gamma,
            bool verbose);
    void scattering_mc(int nphot, bool verbose);

    void mc_iteration();

    void run_image(int nx, int ny, double pixel_size, py::array_t<double> lam, 
            int nphot, double incl, double pa, double dpc);
    void run_spectrum(py::array_t<double> lam, int nphot, double incl, 
            double pa, double dpc);
};
