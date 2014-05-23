#ifndef GRID_CC
#define GRID_CC

#include <cmath>
#include "vector.cc"
#include "dust.cc"
#include "sources.cc"
#include "photon.cc"

struct Grid {
    int n1;
    int n2;
    int n3;
    int nw1;
    int nw2;
    int nw3;
    double *w1;
    double *w2;
    double *w3;
    double *dw1;
    double *dw2;
    double *dw3;
    double ***dens;
    double ***temp;
    double ***mass;
    Dust *dust_species;
    int ***dust;
    Source *sources;
    int nspecies;
    double total_lum;

    Photon *emit();
    virtual double next_wall_distance(Photon *P);
    void propagate_photon(Photon *P, double tau, double ***pcount, bool lucy,
            bool verbose);
    void propagate_ray(Ray *R, bool verbose);
    void absorb(Photon *P, bool lucy);
    void isoscatt(Photon *P);
    virtual Vector<int, 3> photon_loc(Photon *P, bool verbose);
    virtual bool in_grid(Photon *P);
    virtual bool in_grid_raytracing(Ray *R);
    void update_grid(int nphot, Vector<int, 3> l, double ***pcount);
    void update_grid(int nphot, double ***pcount);
    double cell_lum(Vector<int, 3> l);
};

/* Emit a photon from the grid. */

Photon *Grid::emit() {
    Photon *P = sources[0].emit(nspecies, dust_species);
    P->l[0] = -1;
    P->l[1] = -1;
    P->l[2] = -1;

    P->l = photon_loc(P, false);

    return P;
}

/* Linker function to the dust absorb function. */

void Grid::absorb(Photon *P, bool lucy) {
    dust_species[dust[P->l[0]][P->l[1]][P->l[2]]].absorb(P,
            temp[P->l[0]][P->l[1]][P->l[2]], !lucy, dust_species, nspecies);
}

/* Linker function to the dust isoscatt function. */

void Grid::isoscatt(Photon *P) {
    dust_species[dust[P->l[0]][P->l[1]][P->l[2]]].isoscatt(P);
}

/* Propagate a photon through the grid a distance equivalent to tau. */

void Grid::propagate_photon(Photon *P, double tau, double ***pcount, 
        bool lucy, bool verbose) {

    int i = 0;
    while ((tau > 1.0e-6) && (in_grid(P))) {
        double s1 = next_wall_distance(P);

        double s2 = tau/(P->current_kext[dust[P->l[0]][P->l[1]][P->l[2]]]*
                dens[P->l[0]][P->l[1]][P->l[2]]);

        double s = s1;
        if (s2 < s1) s = s2;

        if (verbose) {
            printf("%2i  %7.4f  %i  %7.4f  %7.4f  %7.4f\n", i, tau, P->l[0],
                    P->r[0]/au, s1*P->n[0]/au, s2*P->n[0]/au);
            printf("%14i  %7.4f  %7.4f  %7.4f\n", P->l[1], P->r[1]/au, 
                    s1*P->n[1]/au, s2*P->n[1]/au);
            printf("%14i  %7.4f  %7.4f  %7.4f\n", P->l[2], P->r[2]/au, 
                    s1*P->n[2]/au, s2*P->n[2]/au);
        }

        if (lucy) {
            pcount[P->l[0]][P->l[1]][P->l[2]] += s*
                P->current_kext[dust[P->l[0]][P->l[1]][P->l[2]]]*
                dens[P->l[0]][P->l[1]][P->l[2]];
        }

        tau -= s*P->current_kext[dust[P->l[0]][P->l[1]][P->l[2]]]*
                dens[P->l[0]][P->l[1]][P->l[2]];

        P->move(s);

        if (s1 < s2) P->l = photon_loc(P, verbose);
        i++;
    }
}

/* Propagate a ray through the grid for raytracing. */

void Grid::propagate_ray(Ray *R, bool verbose) {

    int i=0;
    do {
        if (R->l[0] < 0) R->l[0] = 0;
        if (R->l[1] < 0) R->l[1] = 0;
        if (R->l[2] < 0) R->l[2] = 0;
        if (R->l[0] >= nw1-1) R->l[0] = nw1-2;
        if (R->l[1] >= nw2-1) R->l[1] = nw2-2;
        if (R->l[2] >= nw3-1) R->l[2] = nw3-2;

        if (dw1[R->l[0]] < R->pixel_size) {
            R->pixel_too_large = true;
            break;
        }

        double s = next_wall_distance(R);

        double tau_cell = s*R->current_kext[dust[R->l[0]][R->l[1]][R->l[2]]]*
                dens[R->l[0]][R->l[1]][R->l[2]];

        double intensity_cell = (1.0-exp(-tau_cell))*planck_function(R->nu,
                temp[R->l[0]][R->l[1]][R->l[2]]);

        if (verbose) {
            printf("%2i  %7.4f  %i  %7.4f  %7.4f\n", i, tau_cell, 
                    R->l[0], R->r[0]/au, s*R->n[0]/au);
            printf("%11.1e  %i  %7.4f  %7.4f\n", R->intensity, R->l[1], 
                    R->r[1]/au, s*R->n[1]/au);
            printf("%11.4f  %i  %7.4f  %7.4f\n", R->tau, R->l[2], R->r[2]/au, 
                    s*R->n[2]/au);
        }

        R->intensity += intensity_cell*exp(-R->tau);
        R->tau += tau_cell;

        R->move(s);

        R->l = photon_loc(R, verbose);

        i++;
    } while (in_grid(R));
}

/* Calculate the distance between the photon and the nearest wall. */

double Grid::next_wall_distance(Photon *P) {
    return 0.0;
}

/* Determine which cell the photon is in. */

Vector<int, 3> Grid::photon_loc(Photon *P, bool verbose) {
    return Vector<int, 3>();
}

/* Check whether a photon is in the boundaries of the grid. */

bool Grid::in_grid(Photon *P) {
    return true;
}

/* Determine whether a ray is in the boundaries of the grid to start the 
 * raytracing. */

bool Grid::in_grid_raytracing(Ray *R) {
    return true;
}

/* Update the temperature in a cell given the number of photons that have 
 * been absorbed in the cell. */

void Grid::update_grid(int nphot, Vector<int, 3> l, double ***pcount) {
    bool not_converged = true;

    while (not_converged) {
        double T_old = temp[l[0]][l[1]][l[2]];

        temp[l[0]][l[1]][l[2]]=pow(pcount[l[0]][l[1]][l[2]]*
            sources[0].luminosity/(4*sigma*nphot*
            dust_species[dust[l[0]][l[1]][l[2]]].
            planck_mean_opacity(temp[l[0]][l[1]][l[2]])*
            mass[l[0]][l[1]][l[2]]),0.25);

        if ((fabs(T_old-temp[l[0]][l[1]][l[2]])/T_old < 1.0e-2))
            not_converged = false;
    }
}

void Grid::update_grid(int nphot, double ***pcount) {
    for (int i=0; i<nw1-1; i++)
        for (int j=0; j<nw2-1; j++)
            for (int k=0; k<nw3-1; k++)
                update_grid(nphot,Vector<int, 3>(i,j,k),pcount);
}

/* Calculate the luminosity of the cell indicated by l. */

double Grid::cell_lum(Vector<int, 3> l) {
    return 4*mass[l[0]][l[1]][l[2]]*dust_species[dust[l[0]][l[1]][l[2]]].
        planck_mean_opacity(temp[l[0]][l[1]][l[2]])*sigma*
        pow(temp[l[0]][l[1]][l[2]],4);
}

#endif
