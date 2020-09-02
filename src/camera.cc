#include "camera.h"

Image::Image(int _nx, int _ny, double _pixel_size, py::array_t<double> __lam) {
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

    pixel_size = _pixel_size;

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
}

Image::~Image() {
    freepymangle(intensity);
}

UnstructuredImage::UnstructuredImage(int _nx, int _ny, double _pixel_size, 
        py::array_t<double> __lam) {
    // Start by setting up the appropriate Python arrays.

    nx = _nx; ny = _ny;
    _lam = __lam;

    auto _lam_buf = _lam.request();

    if (_lam_buf.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one");

    // Now get the correct format.
    
    lam = (double *) _lam_buf.ptr;

    nnu = _lam_buf.shape[0];

    // Set up the x and y values properly.

    pixel_size = _pixel_size;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            x.push_back((i - nx/2)*pixel_size + (random_number()-0.5)*
                    pixel_size/4.);
            y.push_back((j - ny/2)*pixel_size + (random_number()-0.5)*
                    pixel_size/4.);
            intensity.push_back(std::vector<double>());
        }
    }

    // Set up the frequency array.

    _nu = py::array_t<double>(nnu);
    auto _nu_buf = _nu.request();

    nu = (double *) _nu_buf.ptr;

    for (int i = 0; i < nnu; i++)
        nu[i] = c_l / (lam[i]*1.0e-4);
}

Spectrum::Spectrum(py::array_t<double> __lam) {
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
}

Camera::Camera(Grid *_G, Params *_Q) {
    G = _G;
    Q = _Q;
}

void Camera::set_orientation(double _incl, double _pa, double _dpc) {
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

Image *Camera::make_image(int nx, int ny, double pixel_size, 
        py::array_t<double> lam, double incl, double pa, double dpc) {

    // Set the camera orientation.
    set_orientation(incl, pa, dpc);

    // Set up the image.

    Image *image = new Image(nx, ny, pixel_size*arcsec*dpc*pc, lam);

    // Now go through and raytrace.

    for (int j=0; j<image->nx; j++)
        for (int k=0; k<image->ny; k++) {
            if (Q->verbose) printf("%d   %d\n", j, k);

            double *intensity = raytrace_pixel(image->x[j], image->y[k], 
                    image->pixel_size, image->nu, image->nnu, 0);

            for (int i = 0; i < image->nnu; i++)
                image->intensity[j][k][i] = intensity[i] * 
                    image->pixel_size * image->pixel_size / (r * r)/ Jy;

            delete[] intensity;
        }

    // Also raytrace the sources.

    raytrace_sources(image);

    // And return.

    return image;
}

UnstructuredImage *Camera::make_unstructured_image(int nx, int ny, 
        double pixel_size, py::array_t<double> lam, double incl, double pa, 
        double dpc) {

    // Set the camera orientation.
    set_orientation(incl, pa, dpc);

    // Set up the image.

    UnstructuredImage *image = new UnstructuredImage(nx, ny, pixel_size*arcsec*
            dpc*pc, lam);

    // Now go through and raytrace.

    for (int j=0; j < nx*ny; j++)
        raytrace_pixel(image, j, image->pixel_size);

    // Also raytrace the sources.

    //raytrace_sources(image);

    // Now that thats done, create an intensity Python array to populate.
    //
    image->_x = py::array_t<double>(image->x.size());
    image->_y = py::array_t<double>(image->y.size());
    image->_intensity = py::array_t<double>(image->x.size()*image->nnu);
    image->_intensity.resize({(int)image->x.size(), image->nnu});

    auto _x_buf = image->_x.request();
    auto _y_buf = image->_y.request();
    auto _intensity_buf = image->_intensity.request();

    double *_x_arr = (double *) _x_buf.ptr;
    double *_y_arr = (double *) _y_buf.ptr;
    double *_intensity_arr = (double *) _intensity_buf.ptr;

    for (int i = 0; i < image->x.size(); i++) {
        _x_arr[i] = image->x[i];
        _y_arr[i] = image->y[i];

        for (int j = 0; j < image->nnu; j++)
            _intensity_arr[i*image->nnu + j] = image->intensity[i][j];
    }

    // And return.

    return image;
}

Spectrum *Camera::make_spectrum(py::array_t<double> lam, double incl,
        double pa, double dpc) {

    // Set up parameters for the image.

    int nx = 100;
    int ny = 100;

    double pixel_size = G->grid_size()*1.1/(dpc*pc)/arcsec/nx;

    // Set up and create an image.

    Image *image = make_image(nx, ny, pixel_size, lam, incl, pa, dpc);

    // Sum the image intensity.
    Spectrum *S = new Spectrum(lam);

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

Ray *Camera::emit_ray(double x, double y, double pixel_size, double *nu, 
        int nnu) {
    Ray *R = new Ray();

    R->nnu = nnu;
    R->tau = new double[nnu];
    R->intensity = new double[nnu];
    for (int i = 0; i < nnu; i++) {
        R->tau[i] = 0; R->intensity[i] = 0;
    }

    R->pixel_size = pixel_size;
    R->pixel_too_large = false;

    R->nu = nu;

    R->ndust = G->nspecies;
    R->current_kext = create2DArr(G->nspecies, nnu);
    R->current_albedo = create2DArr(G->nspecies, nnu);

    for (int j=0; j<G->nspecies; j++) {
        for (int k = 0; k < nnu; k++) {
            R->current_kext[j][k] = G->dust[j]->opacity(R->nu[k]);
            R->current_albedo[j][k] = G->dust[j]->albdo(R->nu[k]);
        }
    }

    R->r = i + x*ex + y*ey;
    R->n = ez;

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

double* Camera::raytrace_pixel(double x, double y, double pixel_size, 
        double *nu, int nnu, int count) {
    //printf("%d\n", count);
    bool subpixel = false;

    double *intensity = raytrace(x, y, pixel_size, nu, nnu, false, &subpixel);

    count++;

    if (subpixel) { // && (count < 1)) {
        double *intensity1 = raytrace_pixel(x-pixel_size/4, y-pixel_size/4, 
                pixel_size/2, nu, nnu, count);
        double *intensity2 = raytrace_pixel(x-pixel_size/4, y+pixel_size/4, 
                pixel_size/2, nu, nnu, count);
        double *intensity3 = raytrace_pixel(x+pixel_size/4, y-pixel_size/4, 
                pixel_size/2, nu, nnu, count);
        double *intensity4 = raytrace_pixel(x+pixel_size/4, y+pixel_size/4, 
                pixel_size/2, nu, nnu, count);

        for (int i = 0; i < nnu; i++) {
            intensity[i] = (intensity1[i]+intensity2[i]+intensity3[i]+
                    intensity4[i])/4;
        }

        delete[] intensity1; delete[] intensity2; delete[] intensity3;
        delete[] intensity4;
    }

    return intensity;
}

void Camera::raytrace_pixel(UnstructuredImage *image, int ix, 
        double pixel_size) {

    bool subpixel = false;

    // Raytrace all frequencies.

    double *intensity = raytrace(image->x[ix], image->y[ix], pixel_size, 
            image->nu, image->nnu, true, &subpixel);

    for (int i=0; i<image->nnu; i++)
        image->intensity[ix].push_back(intensity[i]);

    delete[] intensity;
    
    // Split the cell into four and raytrace again.
    if (subpixel) {
        int nxy = image->x.size();

        image->x.push_back(image->x[ix] - pixel_size/4 + (random_number()-0.5)*
                pixel_size/8);
        image->x.push_back(image->x[ix] - pixel_size/4 + (random_number()-0.5)*
                pixel_size/8);
        image->x.push_back(image->x[ix] + pixel_size/4 + (random_number()-0.5)*
                pixel_size/8);
        image->x.push_back(image->x[ix] + pixel_size/4 + (random_number()-0.5)*
                pixel_size/8);

        image->y.push_back(image->y[ix] - pixel_size/4 + (random_number()-0.5)*
                pixel_size/8);
        image->y.push_back(image->y[ix] + pixel_size/4 + (random_number()-0.5)* 
                pixel_size/8);
        image->y.push_back(image->y[ix] - pixel_size/4 + (random_number()-0.5)* 
                pixel_size/8);
        image->y.push_back(image->y[ix] + pixel_size/4 + (random_number()-0.5)* 
                pixel_size/8);

        image->intensity.push_back(std::vector<double>());
        image->intensity.push_back(std::vector<double>());
        image->intensity.push_back(std::vector<double>());
        image->intensity.push_back(std::vector<double>());

        raytrace_pixel(image, nxy+0, pixel_size/2);
        raytrace_pixel(image, nxy+1, pixel_size/2);
        raytrace_pixel(image, nxy+2, pixel_size/2);
        raytrace_pixel(image, nxy+3, pixel_size/2);
    }
}

double* Camera::raytrace(double x, double y, double pixel_size, double *nu, 
        int nnu, bool unstructured, bool *pixel_too_large) {
    /* Emit a ray from the given location. */
    Ray *R = emit_ray(x, y, pixel_size, nu, nnu);

    /* Create an intensity array for the result to go into. */
    double *intensity = new double[nnu];

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

        if (G->on_and_parallel_to_wall(R) and not unstructured) {
            for (int i = 0; i < nnu; i++)
                intensity[i] = -1.0;

            *pixel_too_large = true;

            return intensity;
        }

        /* Move the ray through the grid, calculating the intensity as 
         * you go. */
        if (Q->verbose) printf("\n");
        if (G->in_grid(R))
            G->propagate_ray(R);
        if (Q->verbose) printf("\n");

        /* Check whether the run was successful or if we need to sub-pixel 
         * to get a good intensity measurement. */
        for (int i = 0; i < nnu; i++)
            intensity[i] = R->intensity[i];

        *pixel_too_large = R->pixel_too_large;

        /* Clean up the ray. */
        delete R;

        return intensity;
    }
    else {
        delete R; // Make sure the Ray instance is cleaned up.

        for (int i = 0; i < nnu; i++)
            intensity[i] = 0.0;

        return intensity;
    }
}

void Camera::raytrace_sources(Image *image) {

    for (int isource=0; isource < G->nsources; isource++) {
        for (int iphot=0; iphot < 1000; iphot++) {
            // Emit the ray.
            Ray *R = G->sources[isource]->emit_ray(image->nu, image->nnu, 
                    image->pixel_size, ez, 1000);

            R->l = G->photon_loc(R);

            // Get the appropriate dust opacities.

            R->ndust = G->nspecies;
            R->current_kext = create2DArr(G->nspecies, image->nnu);
            R->current_albedo = create2DArr(G->nspecies, image->nnu);

            for (int j=0; j<G->nspecies; j++) {
                for (int k = 0; k < image->nnu; k++) {
                    R->current_kext[j][k] = G->dust[j]->opacity(R->nu[k]);
                    R->current_albedo[j][k] = G->dust[j]->albdo(R->nu[k]);
                }
            }

            // Now propagate the ray.
            G->propagate_ray_from_source(R);

            // Now bin the photon into the right cell.
            double ximage = R->r * ey;
            double yimage = R->r * ex;

            int ix = int(image->nx * (ximage + image->x[image->nx-1]) / 
                    (2*image->x[image->nx-1]) + 0.5);
            int iy = int(image->ny * (yimage + image->y[image->ny-1]) / 
                    (2*image->y[image->ny-1]) + 0.5);

            // Finally, add the energy into the appropriate cell.
            for (int inu=0; inu < image->nnu; inu++) {
                image->intensity[ix][iy][inu] += R->intensity[inu] *
                    image->pixel_size * image->pixel_size / (r * r)/ Jy;
            }

            // And clean up the Ray.
            delete R;
        }
    }
}

void Camera::raytrace_sources(UnstructuredImage *image) {

    for (int isource=0; isource < G->nsources; isource++) {
        int nxy = image->x.size();

        for (int iphot=0; iphot < 1; iphot++) {
            // Emit the ray.
            Ray *R = G->sources[isource]->emit_ray(image->nu, image->nnu, 
                    ez, 1);

            R->l = G->photon_loc(R);

            // Get the appropriate dust opacities.

            R->ndust = G->nspecies;
            R->current_kext = create2DArr(G->nspecies, image->nnu);
            R->current_albedo = create2DArr(G->nspecies, image->nnu);

            for (int j=0; j<G->nspecies; j++) {
                for (int k = 0; k < image->nnu; k++) {
                    R->current_kext[j][k] = G->dust[j]->opacity(R->nu[k]);
                    R->current_albedo[j][k] = G->dust[j]->albdo(R->nu[k]);
                }
            }

            // Now propagate the ray.
            G->propagate_ray_from_source(R);

            // Now bin the photon into the right cell.
            double ximage = R->r * ey;
            double yimage = R->r * ex;

            image->x.push_back(ximage);
            image->y.push_back(yimage);
            image->intensity.push_back(std::vector<double>());

            // Finally, add the energy into the appropriate cell.
            for (int inu=0; inu < image->nnu; inu++) {
                image->intensity[nxy+iphot].push_back(R->intensity[inu]);
            }

            // And clean up the Ray.
            delete R;
        }
    }
}
