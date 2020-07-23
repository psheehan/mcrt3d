#ifndef CAMERA_H
#define CAMERA_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "vector.h"
#include "photon.h"
#include "grid.h"
#include "misc.h"
#include "params.h"
#include "misc.h"

namespace py = pybind11;

struct Image {
    double r;
    double incl;
    double pa;

    py::array_t<double> _x;
    py::array_t<double> _y;
    py::array_t<double> _intensity;
    py::array_t<double> _lam;
    py::array_t<double> _nu;

    double *x;
    double *y;
    Vector<double, 3> i;
    Vector<double, 3> ex;
    Vector<double, 3> ey;
    Vector<double, 3> ez;
    double *lam;
    double *nu;
    double ***intensity;
    double pixel_size;
    int nx;
    int ny;
    int nnu;

    Image(int nx, int ny, double pixel_size, py::array_t<double> lam, 
            double incl, double pa, double dpc);
};

struct Spectrum {
    double r;
    double incl;
    double pa;

    py::array_t<double> _intensity;
    py::array_t<double> _lam;
    py::array_t<double> _nu;

    Vector<double, 3> i;
    Vector<double, 3> ex;
    Vector<double, 3> ey;
    Vector<double, 3> ez;
    double *lam;
    double *nu;
    double *intensity;
    double pixel_size;
    int nnu;

    Spectrum(py::array_t<double> lam, double incl, double pa,
        double dpc);
};

struct Camera {
    Image* image;
    Grid* G;
    Params *Q;

    Camera(Grid *_G, Params *_Q);

    Image* make_image(int nx, int ny, double pixel_size, 
            py::array_t<double> lam, double incl, double pa, double dpc);
    Spectrum* make_spectrum(py::array_t<double> lam, double incl, 
            double pa, double dpc);

    Ray* emit_ray(double x, double y, double pixel_size, double nu);
    double raytrace_pixel(double x, double y, double pixel_size, double nu, 
            int count);
    double raytrace(double x, double y, double pixel_size, double nu);

    void raytrace_sources(Image *I);
};

#endif
