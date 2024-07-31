#include "camera.h"

Image::Image(int _nx, int _ny, double _pixel_size, py::array_t<double> __lam) {
    // Start by setting up the appropriate Python arrays.

    nx = _nx; ny = _ny;
    Kokkos::resize(x, nx);
    Kokkos::resize(y, ny);

    // Set up the x and y values properly.

    pixel_size = _pixel_size;
    
    for (int i = 0; i < nx; i++)
        x(i) = (i - nx/2)*pixel_size;
    for (int i = 0; i < ny; i++)
        y(i) = (i - ny/2)*pixel_size;

    // Now get the correct format.
    
    _x = array_from_view(x, 1, {(size_t) nx});
    _y = array_from_view(y, 1, {(size_t) ny});
    
    // Set up the wavelength array.

    auto _lam_buf = __lam.request();
    nnu = _lam_buf.shape[0];

    Kokkos::resize(lam, nnu);
    Kokkos::deep_copy(lam, view_from_array(__lam));
    _lam = __lam;

    // Set up the frequency array.

    Kokkos::resize(nu, nnu);
    for (int i = 0; i < nnu; i++)
        nu(i) = c_l / (lam(i)*1.0e-4);
    _nu = array_from_view(nu, 1, {(size_t) nnu});

    // Set up the volume of each cell.

    Kokkos::resize(intensity, nx*ny*nnu);
    for (int i = 0; i < nx*ny*nnu; i++)
        intensity(i) = 0.;

    _intensity = array_from_view(intensity, 3, {(size_t) nx, (size_t) ny, (size_t) nnu});
}

Image::~Image() {
}

UnstructuredImage::UnstructuredImage(int _nx, int _ny, double _pixel_size, 
        py::array_t<double> __lam) {
    // Start by setting up the appropriate Python arrays.

    nx = _nx; ny = _ny;
    _lam = __lam;

    auto _lam_buf = _lam.request();

    if (_lam_buf.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one");

    nnu = _lam_buf.shape[0];

    // Now get the correct format.

    Kokkos::resize(lam, nnu);
    lam = view_from_array(_lam);
    _lam = __lam;

    // Set up the x and y values properly.

    pixel_size = _pixel_size;

    Kokkos::resize(x, nx*ny);
    Kokkos::resize(y, nx*ny);
    Kokkos::resize(intensity, nx*ny, nnu);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            x(i*ny+j) = (i - nx/2)*pixel_size + (random_number()-0.5)*
                    pixel_size/10000;
            y(i*ny+j) = (j - ny/2)*pixel_size + (random_number()-0.5)*
                    pixel_size/10000;
            
            for (int k = 0; k < nnu; k++)
                intensity(i*ny+j,k) = 0.;
        }
    }

    // Set up the frequency array.

    Kokkos::resize(nu, nnu);
    for (int i = 0; i < nnu; i++)
        nu(i) = c_l / (lam(i)*1.0e-4);
    _nu = array_from_view(nu, 1, {(size_t) nnu});
}

UnstructuredImage::UnstructuredImage(int _nr, int _nphi, double rmin, 
        double rmax, py::array_t<double> __lam) {
    // Start by setting up the appropriate Python arrays.

    nx = _nr; ny = _nphi;
    _lam = __lam;

    auto _lam_buf = _lam.request();

    if (_lam_buf.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one");

    nnu = _lam_buf.shape[0];

    // Now get the correct format.

    Kokkos::resize(lam, nnu);
    Kokkos::deep_copy(lam, view_from_array(_lam));
    _lam = __lam;

    // Set up the x and y values properly.

    pixel_size = 0;

    double logrmin = log10(rmin);
    double logrmax = log10(rmax);
    double dlogr = (logrmax - logrmin) / (nx - 1);

    Kokkos::resize(x, nx*ny+1);
    Kokkos::resize(y, nx*ny+1);
    Kokkos::resize(intensity, nx*ny+1, nnu);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double r = pow(10, logrmin + dlogr*i);
            double phi = 2*pi * (j + 0.5) / ny;
            x(i*ny+j) = r*cos(phi);
            y(i*ny+j) = r*sin(phi);
            
            for (int k = 0; k < nnu; k++)
                intensity(i*ny+j,k) = 0.;
        }
    }
    x(x.extent(0)-1) = 0.;
    y(x.extent(0)-1) = 0.;

    // Set up the frequency array.

    Kokkos::resize(nu, nnu);
    for (int i = 0; i < nnu; i++)
        nu(i) = c_l / (lam(i)*1.0e-4);
    _nu = array_from_view(nu, 1, {(size_t) nnu});
}

Spectrum::Spectrum(py::array_t<double> __lam) {
    // Start by setting up the appropriate Python arrays.

    _lam = __lam;
    auto _lam_buf = _lam.request();

    if (_lam_buf.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one");

    nnu = _lam_buf.shape[0];

    // Now get the correct format.

    Kokkos::resize(lam, nnu);
    Kokkos::deep_copy(lam, view_from_array(_lam));
    _lam = __lam;

    // Set up the frequency array.

    Kokkos::resize(nu, nnu);
    for (int i = 0; i < nnu; i++)
        nu(i) = c_l / (lam(i)*1.0e-4);
    _nu = array_from_view(nu, 1, {(size_t) nnu});

    // Set up the volume of each cell.

    Kokkos::resize(intensity, nnu);
    for (int i = 0; i < nnu; i++)
        intensity(i) = 0;
    _intensity = array_from_view(intensity, 1, {(size_t) nnu});
}

Camera::Camera(Grid *_G, Params *_Q) {
    G = _G;
    Q = _Q;
    random_pool = new Kokkos::Random_XorShift64_Pool<>(/*seed=*/12345);
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
        py::array_t<double> lam, double incl, double pa, double dpc, 
        int nthreads) {

    // Set the camera orientation.
    set_orientation(incl, pa, dpc);

    // Set up the image.

    Image *image = new Image(nx, ny, pixel_size*arcsec*dpc*pc, lam);

    // Now go through and raytrace.

    Kokkos::parallel_for("RaytraceRectangularImage", 
            Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0},{image->nx,image->ny}),
            [=] (int64_t j, int64_t k) {
    //for (int j=0; j<image->nx; j++)
    //    for (int k=0; k<image->ny; k++) {
        if (Q->verbose) printf("%d   %d\n", j, k);

        Kokkos::View<double*> intensity(raytrace_pixel(image->x(j), image->y(k), 
                image->pixel_size, image->nu, image->nnu, 0));

        for (int i = 0; i < image->nnu; i++)
            image->intensity(k*image->ny*image->nnu + j*image->nnu + i) = 
                intensity(i) * image->pixel_size * image->pixel_size / 
                (r * r)/ Jy;
    });

    // Also raytrace the sources.

    raytrace_sources(image, nthreads);

    // And return.

    image->_intensity = array_from_view(image->intensity, 3, {(size_t) image->nx, (size_t) image->ny, (size_t) image->nnu});

    return image;
}

UnstructuredImage *Camera::make_unstructured_image(int nx, int ny, 
        double pixel_size, py::array_t<double> lam, double incl, double pa, 
        double dpc, int nthreads) {

    // Set the camera orientation.
    set_orientation(incl, pa, dpc);

    // Set up the image.

    UnstructuredImage *image = new UnstructuredImage(nx, ny, pixel_size*arcsec*
            dpc*pc, lam);

    // Now go through and raytrace.

    seed1 = int(time(NULL));
    seed2 = int(time(NULL));

    Kokkos::parallel_for(nx*ny, [=] (int64_t j) {
        raytrace_pixel(image, j, image->pixel_size);
    });

    // Also raytrace the sources.

    //raytrace_sources(image);

    // Now that thats done, create an intensity Python array to populate.
    
    image->_x = array_from_view(image->x, 1, {image->x.extent(0)});
    image->_y = array_from_view(image->y, 1, {image->y.extent(0)});
    image->_intensity = array_from_view(image->intensity, 2, {image->intensity.extent(0), image->intensity.extent(1)});

    // And return.

    return image;
}

UnstructuredImage *Camera::make_circular_image(int nr, int nphi, 
        py::array_t<double> lam, double incl, double pa, double dpc, 
        int nthreads) {

    // Set the camera orientation.
    set_orientation(incl, pa, dpc);

    // Set up the image.

    UnstructuredImage *image = new UnstructuredImage(nr, nphi, 
            G->w1[1]*0.05, 0.5*G->grid_size()*1.05, lam);

    // Now go through and raytrace.

    seed1 = int(time(NULL));
    seed2 = int(time(NULL));

    Kokkos::parallel_for(image->x.extent(0), [=] (int64_t j) {
        raytrace_pixel(image, j, image->pixel_size);
    });

    // Also raytrace the sources.

    //raytrace_sources(image);

    // Now that thats done, create an intensity Python array to populate.

    image->_x = array_from_view(image->x, 1, {image->x.extent(0)});
    image->_y = array_from_view(image->y, 1, {image->y.extent(0)});
    image->_intensity = array_from_view(image->intensity, 2, {image->intensity.extent(0), image->intensity.extent(1)});

    // And return.

    return image;
}

Spectrum *Camera::make_spectrum(py::array_t<double> lam, double incl,
        double pa, double dpc, int nthreads) {

    // Set up parameters for the image.

    int nx = 100;
    int ny = 100;

    double pixel_size = G->grid_size()*1.1/(dpc*pc)/arcsec/nx;

    // Set up and create an image.

    Image *image = make_image(nx, ny, pixel_size, lam, incl, pa, dpc, 
            nthreads);

    // Sum the image intensity.
    Spectrum *S = new Spectrum(lam);

    for (int i=0; i<image->nnu; i++) {
        for (int j=0; j<image->nx; j++) {
            for (int k=0; k<image->ny; k++) {
                S->intensity(i) += image->intensity(j*image->nnu*image->ny 
                    + k*image->nnu + i);
            }
        }
    }

    S->_intensity = array_from_view(S->intensity, 1, {(size_t) S->nnu});

    // Delete the parts of the image we no longer need.
    delete image;

    // And return the spectrum.

    return S;
}

Ray *Camera::emit_ray(double x, double y, double pixel_size, Kokkos::View<double*> nu, 
        int nnu) {
    Ray *R = new Ray();

    R->nnu = nnu;
    Kokkos::resize(R->nu, nnu);
    Kokkos::resize(R->tau, nnu);
    Kokkos::resize(R->intensity, nnu);
    for (int i = 0; i < nnu; i++) {
        R->nu(i) = nu(i);
        R->tau(i) = 0; R->intensity(i) = 0;
    }

    R->pixel_size = pixel_size;
    R->pixel_too_large = false;

    R->ndust = G->nspecies;
    Kokkos::resize(R->current_kext, G->nspecies, nnu);
    Kokkos::resize(R->current_albedo, G->nspecies, nnu);

    for (int j=0; j<G->nspecies; j++) {
        for (int k = 0; k < nnu; k++) {
            R->current_kext(j,k) = G->dust[j]->opacity(R->nu(k));
            R->current_albedo(j,k) = G->dust[j]->albdo(R->nu(k));
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

Kokkos::View<double*> Camera::raytrace_pixel(double x, double y, double pixel_size, 
        Kokkos::View<double*> nu, int nnu, int count) {
    //printf("%d\n", count);
    bool subpixel = false;

    Kokkos::View<double*> intensity(raytrace(x, y, pixel_size, nu, nnu, false, &subpixel));

    count++;

    if (subpixel) { // && (count < 1)) {
        Kokkos::View<double*> intensity1(raytrace_pixel(x-pixel_size/4, y-pixel_size/4, 
                pixel_size/2, nu, nnu, count));
        Kokkos::View<double*> intensity2(raytrace_pixel(x-pixel_size/4, y+pixel_size/4, 
                pixel_size/2, nu, nnu, count));
        Kokkos::View<double*> intensity3(raytrace_pixel(x+pixel_size/4, y-pixel_size/4, 
                pixel_size/2, nu, nnu, count));
        Kokkos::View<double*> intensity4(raytrace_pixel(x+pixel_size/4, y+pixel_size/4, 
                pixel_size/2, nu, nnu, count));

        for (int i = 0; i < nnu; i++) {
            intensity(i) = (intensity1(i)+intensity2(i)+intensity3(i)+
                    intensity4(i))/4;
        }
    }

    return intensity;
}

void Camera::raytrace_pixel(UnstructuredImage *image, int ix, 
        double pixel_size) {
    bool subpixel = false;

    // Raytrace all frequencies.

    Kokkos::View<double*> intensity(raytrace(image->x(ix), image->y(ix), pixel_size, 
            image->nu, image->nnu, true, &subpixel));

    for (int i=0; i<image->nnu; i++)
        image->intensity(ix,i) = intensity(i);
    
    // Split the cell into four and raytrace again.
    if (subpixel) {
        int nxy;
        //#pragma omp critical
        //{
        nxy = image->x.extent(0);

        Kokkos::resize(image->x, image->x.extent(0)+4);
        Kokkos::resize(image->y, image->y.extent(0)+4);
        Kokkos::resize(image->intensity, image->intensity.extent(0)+4, image->intensity.extent(1));

        for (int i=0; i < 4; i++) {
            image->x(nxy+i) = image->x(ix) + pow(-1,i)*pixel_size/4 + (random_number(random_pool)-0.5)*
                pixel_size/10000;
            image->y(nxy+i) = image->y(ix) + pow(-1, (int) i < 2)*pixel_size/4 + (random_number(random_pool)-0.5)*
                pixel_size/10000;
        }
        //}

        raytrace_pixel(image, nxy+0, pixel_size/2);
        raytrace_pixel(image, nxy+1, pixel_size/2);
        raytrace_pixel(image, nxy+2, pixel_size/2);
        raytrace_pixel(image, nxy+3, pixel_size/2);
    }
}

Kokkos::View<double*> Camera::raytrace(double x, double y, double pixel_size, Kokkos::View<double*> nu, 
        int nnu, bool unstructured, bool *pixel_too_large) {
    /* Emit a ray from the given location. */
    Ray *R = emit_ray(x, y, pixel_size, nu, nnu);

    /* Create an intensity array for the result to go into. */
    Kokkos::View<double*> intensity("intensity", nnu);

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
                intensity(i) = -1.0;

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
            intensity(i) = R->intensity(i);

        *pixel_too_large = R->pixel_too_large;

        /* Clean up the ray. */
        delete R;

        return intensity;
    }
    else {
        delete R; // Make sure the Ray instance is cleaned up.

        for (int i = 0; i < nnu; i++)
            intensity(i) = 0.0;

        return intensity;
    }
}

void Camera::raytrace_sources(Image *image, int nthreads) {

    #pragma omp parallel for num_threads(nthreads) schedule(guided) collapse(2)
    for (int isource=0; isource < G->nsources; isource++) {
        for (int iphot=0; iphot < 1000; iphot++) {
            // Emit the ray.
            Ray *R = G->sources[isource]->emit_ray(image->nu, image->nnu, 
                    image->pixel_size, ez, 1000);

            R->l = G->photon_loc(R);

            // Get the appropriate dust opacities.

            R->ndust = G->nspecies;
            Kokkos::resize(R->current_kext, G->nspecies, image->nnu);
            Kokkos::resize(R->current_albedo, G->nspecies, image->nnu);

            for (int j=0; j<G->nspecies; j++) {
                for (int k = 0; k < image->nnu; k++) {
                    R->current_kext(j,k) = G->dust[j]->opacity(R->nu(k));
                    R->current_albedo(j,k) = G->dust[j]->albdo(R->nu(k));
                }
            }

            // Now propagate the ray.
            G->propagate_ray_from_source(R);

            // Now bin the photon into the right cell.
            double ximage = R->r * ey;
            double yimage = R->r * ex;

            int ix = int(image->nx * (ximage + image->x(image->nx-1)) / 
                    (2*image->x(image->nx-1)) + 0.5);
            int iy = int(image->ny * (yimage + image->y(image->ny-1)) / 
                    (2*image->y(image->ny-1)) + 0.5);

            // Finally, add the energy into the appropriate cell.
            #pragma omp critical
            {
            for (int inu=0; inu < image->nnu; inu++) {
                image->intensity(ix*image->ny*image->nnu + iy*image->nnu + inu)
                    += R->intensity(inu) * image->pixel_size * 
                    image->pixel_size / (r * r)/ Jy;
            }
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
            Kokkos::resize(R->current_kext, G->nspecies, image->nnu);
            Kokkos::resize(R->current_albedo, G->nspecies, image->nnu);

            for (int j=0; j<G->nspecies; j++) {
                for (int k = 0; k < image->nnu; k++) {
                    R->current_kext(j,k) = G->dust[j]->opacity(R->nu(k));
                    R->current_albedo(j,k) = G->dust[j]->albdo(R->nu(k));
                }
            }

            // Now propagate the ray.
            G->propagate_ray_from_source(R);

            // Now bin the photon into the right cell.
            double ximage = R->r * ey;
            double yimage = R->r * ex;

            Kokkos::resize(image->x, image->x.extent(0)+1);
            Kokkos::resize(image->y, image->y.extent(0)+1);
            Kokkos::resize(image->intensity, image->intensity.extent(0)+1, image->intensity.extent(1));

            // Finally, add the energy into the appropriate cell.
            for (int inu=0; inu < image->nnu; inu++) {
                image->intensity(nxy+iphot,inu) = R->intensity(inu);
            }

            // And clean up the Ray.
            delete R;
        }
    }
}
