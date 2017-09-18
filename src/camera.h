#ifndef CAMERA_H
#define CAMERA_H

#include "vector.cc"
#include "photon.h"
#include "grid.h"
#include "misc.h"
#include "params.h"

struct Image {
    double *x;
    double *y;
    Vector<double, 3> i;
    Vector<double, 3> ex;
    Vector<double, 3> ey;
    Vector<double, 3> ez;
    double *nu;
    double ***intensity;
    double pixel_size;
    int nx;
    int ny;
    int nnu;

    Image(double r, double incl, double pa, double *_x, double *_y, 
            double *_intensity, int _nx, int _ny, double *_nu, 
            double _pixel_size, int _nnu);
};

struct Camera {
    Image* image;
    Grid* G;
    Params *Q;

    Camera(Grid *_G, Params *_Q);

    void make_image(Image *I);
    Ray* emit_ray(double x, double y, double pixel_size, double nu);
    double raytrace_pixel(double x, double y, double pixel_size, double nu, 
            int count);
    double raytrace(double x, double y, double pixel_size, double nu);
};

#endif
