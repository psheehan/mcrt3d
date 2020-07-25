#include "camera.h"

Image::Image(int _nx, int _ny, double _pixel_size, py::array_t<double> __lam,
            double _incl, double _pa, double _dpc) {
    // Start by setting up the appropriate Python arrays.

    nx = _nx; ny = _ny;
    _x = py::array_t<double>(nx);
    _y = py::array_t<double>(ny);
    _lam = __lam;

    auto _x_buf = _x.request(); auto _y_buf = _y.request(); 
    auto _lam_buf = _lam.request();

    if (_lam_buf.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one");

    // Now get the correct format.
    
    x = (double *) _x_buf.ptr; 
    y = (double *) _y_buf.ptr; 
    lam = (double *) _lam_buf.ptr;

    nnu = _lam_buf.shape[0];

    // Set up the x and y values properly.

    pixel_size = _pixel_size * arcsec * _dpc*pc;

    for (int i = 0; i < nx; i++)
        x[i] = (i - nx/2)*pixel_size;
    for (int i = 0; i < ny; i++)
        y[i] = (i - ny/2)*pixel_size;

    // Set up the frequency array.

    _nu = py::array_t<double>(nnu);
    auto _nu_buf = _nu.request();

    nu = (double *) _nu_buf.ptr;

    for (int i = 0; i < nnu; i++)
        nu[i] = c_l / (lam[i]*1.0e-4);

    // Set up the volume of each cell.

    _intensity = py::array_t<double>(nx*ny*nnu);
    _intensity.resize({nx, ny, nnu});

    auto _intensity_buf = _intensity.request();
    intensity = pymangle(nx, ny, nnu, (double *) _intensity_buf.ptr);

    set3DArrValue(intensity, 0., nx, ny, nnu);
    
    // Set viewing angle parameters.

    r = _dpc*pc;
    incl = _incl * pi/180.;
    pa = _pa * pi/180.;

    double phi = -pi/2 - pa;

    i[0] = r*sin(incl)*cos(phi);
    i[1] = r*sin(incl)*sin(phi);
    i[2] = r*cos(incl);

    ex[0] = -sin(phi);
    ex[1] = cos(phi);
    ex[2] = 0.0;

    ey[0] = -cos(incl)*cos(phi);
    ey[1] = -cos(incl)*sin(phi);
    ey[2] = sin(incl);

    ez[0] = -sin(incl)*cos(phi);
    ez[1] = -sin(incl)*sin(phi);
    ez[2] = -cos(incl);
}

Image::~Image() {
    freepymangle(intensity);
}

Spectrum::Spectrum(py::array_t<double> __lam, double _incl, double _pa, 
        double _dpc) {
    // Start by setting up the appropriate Python arrays.

    _lam = __lam;
    auto _lam_buf = _lam.request();

    if (_lam_buf.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one");

    // Now get the correct format.
    
    lam = (double *) _lam_buf.ptr;

    nnu = _lam_buf.shape[0];

    // Set up the frequency array.

    _nu = py::array_t<double>(nnu);
    auto _nu_buf = _nu.request();

    nu = (double *) _nu_buf.ptr;

    for (int i = 0; i < nnu; i++)
        nu[i] = c_l / (lam[i]*1.0e-4);

    // Set up the volume of each cell.

    _intensity = py::array_t<double>(nnu);

    auto _intensity_buf = _intensity.request();
    intensity = (double *) _intensity_buf.ptr;

    for (int i = 0; i < nnu; i++)
        intensity[i] = 0;
    
    // Set viewing angle parameters.

    r = _dpc*pc;
    incl = _incl * pi/180.;
    pa = _pa * pi/180.;

    double phi = -pi/2 - pa;

    i[0] = r*sin(incl)*cos(phi);
    i[1] = r*sin(incl)*sin(phi);
    i[2] = r*cos(incl);

    ex[0] = -sin(phi);
    ex[1] = cos(phi);
    ex[2] = 0.0;

    ey[0] = -cos(incl)*cos(phi);
    ey[1] = -cos(incl)*sin(phi);
    ey[2] = sin(incl);

    ez[0] = -sin(incl)*cos(phi);
    ez[1] = -sin(incl)*sin(phi);
}

Camera::Camera(Grid *_G, Params *_Q) {
    G = _G;
    Q = _Q;
}

Image *Camera::make_image(int nx, int ny, double pixel_size, 
        py::array_t<double> lam, double incl, double pa, double dpc) {

    // Set up the image.

    Image *I = new Image(nx, ny, pixel_size, lam, incl, pa, dpc);
    image = I;

    // Now go through and raytrace.

    for (int i=0; i<image->nnu; i++)
    {
        Q->inu = i;
        for (int j=0; j<image->nx; j++)
            for (int k=0; k<image->ny; k++) {
                if (Q->verbose) printf("%d   %d\n", j, k);
                image->intensity[j][k][i] = raytrace_pixel(image->x[j], 
                        image->y[k], image->pixel_size, image->nu[i], 0) * 
                        image->pixel_size * image->pixel_size / 
                        (image->r * image->r)/ Jy;
            }
    }

    // Also raytrace the sources.

    raytrace_sources(image);

    // And return.

    return I;
}

Spectrum *Camera::make_spectrum(py::array_t<double> lam, double incl,
        double pa, double dpc) {
    // Set up parameters for the image.

    int nx = 100;
    int ny = 100;

    double pixel_size = G->grid_size()*1.1/(dpc*pc)/nx/arcsec;

    // Set up and create an image.

    Image *image = make_image(nx, ny, pixel_size, lam, incl, pa, dpc);

    // Sum the image intensity.
    Spectrum *S = new Spectrum(lam, incl, pa, dpc);

    for (int i=0; i<image->nnu; i++)
        for (int j=0; j<image->nx; j++)
            for (int k=0; k<image->ny; k++) {
                S->intensity[i] += image->intensity[j][k][i];
            }

    // Delete the parts of the image we no longer need.
    delete image;

    // And return the spectrum.

    return S;
}

Ray *Camera::emit_ray(double x, double y, double pixel_size, double nu) {
    Ray *R = new Ray();

    R->tau = 0.0;
    R->intensity = 0.0;
    R->pixel_size = pixel_size;
    R->pixel_too_large = false;

    R->nu = nu;

    double *current_kext = new double[G->nspecies];
    double *current_albedo = new double[G->nspecies];

    for (int j=0; j<G->nspecies; j++) {
        current_kext[j] = G->dust[j]->opacity(R->nu);
        current_albedo[j] = G->dust[j]->albdo(R->nu);
    }

    R->current_kext = current_kext;
    R->current_albedo = current_albedo;

    R->r = image->i + x*image->ex + y*image->ey;
    R->n = image->ez;

    if (equal_zero(R->n[0],EPSILON)) R->n[0] = 0.;
    if (equal_zero(R->n[1],EPSILON)) R->n[1] = 0.;
    if (equal_zero(R->n[2],EPSILON)) R->n[2] = 0.;

    R->invn[0] = 1.0/R->n[0];
    R->invn[1] = 1.0/R->n[1];
    R->invn[2] = 1.0/R->n[2];

    R->l[0] = -1;
    R->l[1] = -1;
    R->l[2] = -1;

    //R->l = G->photon_loc(R, false);
    
    return R;
}

double Camera::raytrace_pixel(double x, double y, double pixel_size, 
        double nu, int count) {
    //printf("%d\n", count);
    double intensity = raytrace(x, y, pixel_size, nu);

    count++;

    if ((intensity < 0)) { // && (count < 1)) {
        double intensity1 = raytrace_pixel(x-pixel_size/4, y-pixel_size/4, 
                pixel_size/2, nu, count);
        double intensity2 = raytrace_pixel(x-pixel_size/4, y+pixel_size/4, 
                pixel_size/2, nu, count);
        double intensity3 = raytrace_pixel(x+pixel_size/4, y-pixel_size/4, 
                pixel_size/2, nu, count);
        double intensity4 = raytrace_pixel(x+pixel_size/4, y+pixel_size/4, 
                pixel_size/2, nu, count);

        return (intensity1+intensity2+intensity3+intensity4)/4;
    }
    else
        return intensity;
}

double Camera::raytrace(double x, double y, double pixel_size, double nu) {
    /* Emit a ray from the given location. */
    Ray *R = emit_ray(x, y, pixel_size, nu);

    /* Move the ray onto the grid boundary */
    double s = G->outer_wall_distance(R);
    if (Q->verbose) printf("%7.4f   %7.4f   %7.4f\n", R->r[0]/au, 
            R->r[1]/au, R->r[2]/au);
    if (Q->verbose) printf("%7.4f\n", s/au);

    if (s != HUGE_VAL) {
        R->move(s);

        if (Q->verbose) printf("%7.4f   %7.4f   %7.4f\n", R->r[0]/au, 
                R->r[1]/au, R->r[2]/au);
        R->l = G->photon_loc(R);

        /* Check whether this photon happens to fall on a wall and is traveling
         * along that wall. */

        if (G->on_and_parallel_to_wall(R)) {
            return -1.0;
        }

        /* Move the ray through the grid, calculating the intensity as 
         * you go. */
        if (Q->verbose) printf("\n");
        if (G->in_grid(R))
            G->propagate_ray(R);
        if (Q->verbose) printf("\n");

        /* Check whether the run was successful or if we need to sub-pixel 
         * to get a good intensity measurement. */
        double intensity = R->intensity;
        if (R->pixel_too_large)
            intensity = -1.0;

        /* Clean up the ray. */
        delete R;

        return intensity;
    }
    else {
        delete R; // Make sure the Ray instance is cleaned up.

        return 0.0;
    }
}

void Camera::raytrace_sources(Image *I) {
    Image *image = I;

    for (int isource=0; isource < G->nsources; isource++) {
        for (int inu=0; inu < image->nnu; inu++) {
            double dnu = image->nu[inu+1] - image->nu[inu];

            for (int iphot=0; iphot < 1000; iphot++) {
                // Emit the ray.
                Ray *R = G->sources[isource]->emit_ray(image->nu[inu], dnu, 
                        image->pixel_size, image->ez, 1000);

                R->l = G->photon_loc(R);

                // Get the appropriate dust opacities.

                double *current_kext = new double[G->nspecies];
                double *current_albedo = new double[G->nspecies];

                for (int j=0; j<G->nspecies; j++) {
                    current_kext[j] = G->dust[j]->opacity(R->nu);
                    current_albedo[j] = G->dust[j]->albdo(R->nu);
                }

                R->current_kext = current_kext;
                R->current_albedo = current_albedo;

                // Now propagate the ray.
                G->propagate_ray_from_source(R);

                // Now bin the photon into the right cell.
                double ximage = R->r * image->ey;
                double yimage = R->r * image->ex;

                int ix = int(image->nx * (ximage + image->x[image->nx-1]) / 
                        (2*image->x[image->nx-1]) + 0.5);
                int iy = int(image->ny * (yimage + image->y[image->ny-1]) / 
                        (2*image->y[image->ny-1]) + 0.5);

                // Finally, add the energy into the appropriate cell.
                image->intensity[ix][iy][inu] += R->intensity;

                // And clean up the Ray.
                delete R;
            }
        }
    }
}
