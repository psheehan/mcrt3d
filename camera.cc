#include "vector.cc"
#include "photon.cc"
#include "grid.cc"
#include "misc.cc"

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
            double ***_intensity, int _nx, int _ny, double *_nu, 
            double _pixel_size, int _nnu) {

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
};

struct Camera {
    Image* image;
    Grid* G;

    void make_image();
    Ray* emit_ray(double x, double y, double pixel_size, double nu);
    double raytrace_pixel(double x, double y, double pixel_size, double nu);
    double raytrace(double x, double y, double pixel_size, double nu);
};

void Camera::make_image() {
    for (int i=0; i<image->nnu; i++)
        for (int j=0; j<image->nx; j++)
            for (int k=0; k<image->ny; k++)
                image->intensity[j][k][i] = raytrace_pixel(image->x[j], 
                        image->y[k], image->pixel_size, image->nu[i]);
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
        current_kext[j] = G->dust_species[j].opacity(R->nu);
        current_albedo[j] = G->dust_species[j].albdo(R->nu);
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

    R->l = G->photon_loc(R, false);
    
    return R;
}

double Camera::raytrace_pixel(double x, double y, double pixel_size, 
        double nu) {
    double intensity = raytrace(x, y, pixel_size, nu);

    if (intensity < 0) {
        double intensity1 = raytrace_pixel(x-pixel_size/4, y-pixel_size/4, 
                pixel_size/2, nu);
        double intensity2 = raytrace_pixel(x-pixel_size/4, y+pixel_size/4, 
                pixel_size/2, nu);
        double intensity3 = raytrace_pixel(x+pixel_size/4, y-pixel_size/4, 
                pixel_size/2, nu);
        double intensity4 = raytrace_pixel(x+pixel_size/4, y+pixel_size/4, 
                pixel_size/2, nu);

        return (intensity1+intensity2+intensity3+intensity4)/4;
    }
    else
        return intensity;
}

double Camera::raytrace(double x, double y, double pixel_size, double nu) {
    Ray *R = emit_ray(x, y, pixel_size, nu);

    double s;
    do {
        s = G->next_wall_distance(R, false);

        if (s == HUGE_VAL)
            break;
        else {
            R->move(s);

            R->l = G->photon_loc(R, false);
        }
    } while ((!equal(s,HUGE_VAL,1.0e-6)) && (!G->in_grid(R)));

    if (G->in_grid(R))
        G->propagate_ray(R, false);

    double intensity = R->intensity;
    if (R->pixel_too_large)
        intensity = -1.0;

    delete R;

    return intensity;
}
