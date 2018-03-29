#include "camera.h"

Image::Image(double _r, double _incl, double _pa, double *_x, double *_y, 
            double *_intensity, int _nx, int _ny, double *_nu, 
            double _pixel_size, int _nnu) {

    r = _r;
    incl = _incl;
    pa = _pa;

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

    x = _x;
    y = _y;
    intensity = pymangle(_nx, _ny, _nnu, _intensity);
    nx = _nx;
    ny = _ny;
    nu = _nu;
    nnu = _nnu;

    pixel_size = _pixel_size;
}

Image::Image(double _r, double _incl, double _pa, double *_x, double *_y, 
            double ***_intensity, int _nx, int _ny, double *_nu, 
            double _pixel_size, int _nnu) {

    r = _r;
    incl = _incl;
    pa = _pa;

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

    x = _x;
    y = _y;
    intensity = _intensity;
    nx = _nx;
    ny = _ny;
    nu = _nu;
    nnu = _nnu;

    pixel_size = _pixel_size;
}

Spectrum::Spectrum(double _r, double _incl, double _pa, double *_intensity, 
        double *_nu, double _pixel_size, int _nnu) {

    r = _r;
    incl = _incl;
    pa = _pa;

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

    intensity = _intensity;
    nu = _nu;
    nnu = _nnu;

    pixel_size = _pixel_size;
}

Camera::Camera(Grid *_G, Params *_Q) {
    G = _G;
    Q = _Q;
}

void Camera::make_image(Image *I) {
    image = I;

    for (int i=0; i<image->nnu; i++)
    {
        Q->inu = i;
        for (int j=0; j<image->nx; j++)
            for (int k=0; k<image->ny; k++) {
                if (Q->verbose) printf("%d   %d\n", j, k);
                image->intensity[j][k][i] = raytrace_pixel(image->x[j], 
                        image->y[k], image->pixel_size, image->nu[i], 0);
            }
    }

    raytrace_sources(image);
}

void Camera::make_spectrum(Spectrum *S) {
    // Set up parameters for the image.
    int nx = 100;
    int ny = 100;

    double *x = new double[nx];
    double *y = new double[ny];

    for (int i=0; i<nx; i++)
        x[i] = (i - nx/2)*S->pixel_size;
    for (int i=0; i<ny; i++)
        y[i] = (i - ny/2)*S->pixel_size;
    
    double ***intensity = create3DArrValue(nx, ny, S->nnu, 0.);

    // Set up and create an image.
    Image *image = new Image(S->r, S->incl, S->pa, x, y, intensity,  
            nx, ny, S->nu, S->pixel_size, S->nnu);

    // Run the image.
    make_image(image);

    // Sum the image intensity.
    for (int i=0; i<image->nnu; i++)
        for (int j=0; j<image->nx; j++)
            for (int k=0; k<image->ny; k++) {
                S->intensity[i] += image->intensity[j][k][i];
            }

    // Delete the parts of the image we no longer need.
    delete[] x;
    delete[] y;
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
    else
        return 0.0;
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
                        (2*image->x[image->nx-1]));
                int iy = int(image->ny * (yimage + image->y[image->ny-1]) / 
                        (2*image->y[image->ny-1]));

                // Finally, add the energy into the appropriate cell.
                image->intensity[ix][iy][inu] += R->intensity;
            }
        }
    }
}
