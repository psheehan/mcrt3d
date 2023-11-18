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
    py::array_t<double> _x;
    py::array_t<double> _y;
    py::array_t<double> _intensity;
    py::array_t<double> _lam;
    py::array_t<double> _nu;

    double *x;
    double *y;
    double *lam;
    double *nu;
    double *intensity;

    double pixel_size;
    int nx;
    int ny;
    int nnu;

    Image(int nx, int ny, double pixel_size, py::array_t<double> lam);

    ~Image();
};

struct UnstructuredImage {
    py::array_t<double> _x;
    py::array_t<double> _y;
    py::array_t<double> _intensity;
    py::array_t<double> _lam;
    py::array_t<double> _nu;

    std::vector<double> x;
    std::vector<double> y;
    double* lam;
    double* nu;
    std::vector<std::vector<double>> intensity;

    double pixel_size;
    int nx;
    int ny;
    int nnu;

    UnstructuredImage(int nx, int ny, double pixel_size, 
            py::array_t<double> lam);
    UnstructuredImage(int nr, int nphi, double rmin, double rmax, 
            py::array_t<double> lam);
};

struct Spectrum {
    py::array_t<double> _intensity;
    py::array_t<double> _lam;
    py::array_t<double> _nu;

    double *lam;
    double *nu;
    double *intensity;
    double pixel_size;
    int nnu;

    Spectrum(py::array_t<double> lam);
};

struct Camera {
    Grid* G;
    Params *Q;

    double r;
    double incl;
    double pa;

    Vector<double, 3> i;
    Vector<double, 3> ex;
    Vector<double, 3> ey;
    Vector<double, 3> ez;

    Camera(Grid *_G, Params *_Q);

    Image* make_image(int nx, int ny, double pixel_size, 
            py::array_t<double> lam, double incl, double pa, double dpc, 
            int nthreads);

    UnstructuredImage* make_unstructured_image(int nx, int ny, 
            double pixel_size, py::array_t<double> lam, double incl, 
            double pa, double dpc, int nthreads);

    UnstructuredImage* make_circular_image(int nr, int nphi, 
            py::array_t<double> lam, double incl, double pa, double dpc, 
            int nthreads);

    Spectrum* make_spectrum(py::array_t<double> lam, double incl, 
            double pa, double dpc, int nthreads);

    void set_orientation(double incl, double pa, double dpc);

    Ray* emit_ray(double x, double y, double pixel_size, double *nu, int nnu);
    double* raytrace_pixel(double x, double y, double pixel_size, double *nu, 
            int nnu, int count);
    void raytrace_pixel(UnstructuredImage *image, int ix, double pixel_size); 
    double* raytrace(double x, double y, double pixel_size, double *nu, 
            int nnu, bool unstructured, bool *pixel_too_large);

    void raytrace_sources(Image *I, int nthreads);
    void raytrace_sources(UnstructuredImage *I);
};

#endif
