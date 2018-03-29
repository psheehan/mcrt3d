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

struct MCRT {
    Grid *G;
    Params *Q;
    Camera *C;

    MCRT(Grid *G, Params *Q);

    void thermal_mc();
    void scattering_mc();
    void mc_iteration();
    void run_image(Image *I);
    void run_spectrum(Spectrum *S);
};
