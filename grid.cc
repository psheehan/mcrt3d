#ifndef GRID_CC
#define GRID_CC

#include <cmath>
#include "vector.cc"
#include "dust.cc"
#include "sources.cc"
#include "photon.cc"
#include "misc.cc"

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
    double ****dens;
    double ****energy;
    double ****temp;
    double ****mass;
    double ***volume;
    int nspecies;
    Dust *dust_species;
    int nsources;
    Source *sources;
    double total_lum;

    Photon *emit(int nphot, int iphot);
    virtual double next_wall_distance(Photon *P, bool verbose);
    virtual double outer_wall_distance(Photon *P, bool verbose);
    void propagate_photon_full(Photon *P, int nphot, 
            bool bw, bool verbose);
    void propagate_photon(Photon *P, double tau, bool absorb,
            bool bw, bool verbose, int nphot);
    void propagate_ray(Ray *R, bool verbose);
    void absorb(Photon *P, int idust, bool bw);
    void isoscatt(Photon *P, int idust);
    virtual Vector<int, 3> photon_loc(Photon *P, bool verbose);
    virtual bool in_grid(Photon *P);
    void update_grid(int nphot, Vector<int, 3> l);
    void update_grid(int nphot);
    double cell_lum(Vector<int, 3> l);
};

/* Emit a photon from the grid. */

Photon *Grid::emit(int nphot, int iphot) {
    /* Cycle through the various stars, having them emit photons one after 
     * another. This way each source will get the same nuber of photons 
     * +/- 1. */
    int isource = 0;
    int photons_per_source = nphot;
    if (nsources > 1) {
        int isource = fmod(iphot, nsources);
        int photons_per_source = int(nphot/nsources);
    }

    Photon *P = sources[isource].emit(photons_per_source, nspecies, 
            dust_species);

    /* Check the photon's location in the grid. */
    P->l = photon_loc(P, false);

    return P;
}

/* Linker function to the dust absorb function. */

void Grid::absorb(Photon *P, int idust, bool bw) {
    dust_species[idust].absorb(P, temp[idust][P->l[0]][P->l[1]][P->l[2]], bw, 
            dust_species, nspecies);

    // Check the photon's location again because there's a small chance that 
    // the photon was absorbed on a wall, and if it was we may need to update
    // which cell it is in if the direction has changed.
    P->l = photon_loc(P, false);
}

/* Linker function to the dust isoscatt function. */

void Grid::isoscatt(Photon *P, int idust) {
    dust_species[idust].isoscatt(P);

    // Check the photon's location again because there's a small chance that 
    // the photon was absorbed on a wall, and if it was we may need to update
    // which cell it is in if the direction has changed.
    P->l = photon_loc(P, false);
}

/* Propagate a photon through the grid until it escapes. */

void Grid::propagate_photon_full(Photon *P, int nphot, 
        bool bw, bool verbose) {
    while (in_grid(P)) {
        // Determin the optical depth that the photon can travel until it's
        // next interaction.
        double tau = -log(1-random_number());

        // Figure out what that next action is, absorption or scattering. This
        // is figured out early for the sake of the continuous absorption
        // method.
        double albedo;
        int idust;
        if (nspecies == 1) {
            albedo = P->current_albedo[0];
            idust = 0;
        }
        else {
            double ksca_tot = 0;
            double *ksca_cum = new double[nspecies];
            double kext_tot = 0;
            for (int i=0; i<nspecies; i++) {
                ksca_tot += P->current_albedo[i]*P->current_kext[i]*
                    dens[i][P->l[0]][P->l[1]][P->l[2]];
                ksca_cum[i] = ksca_tot;
                kext_tot += P->current_kext[i]*
                    dens[i][P->l[0]][P->l[1]][P->l[2]];
            }
            albedo = ksca_tot / kext_tot;

            double rand = random_number();
            for (int i=0; i<nspecies; i++) {
                if (rand < ksca_cum[i] / kext_tot) {
                    idust = i;
                    break;
                }
            }
            delete ksca_cum;
        }

        bool absorb_photon = random_number() > albedo;

        // Move the photon to the point of it's next interaction.
        propagate_photon(P, tau, absorb_photon, bw, verbose, nphot);

        // If the photon is still in the grid when it reaches it's 
        // destination...
        if (in_grid(P)) {
            // If the next interaction is absorption...
            if (absorb_photon) {
                absorb(P, idust, bw);
                // If we've asked for verbose output, print some info.
                if (verbose) {
                    printf("Absorbing photon at %i  %i  %i\n", P->l[0],
                            P->l[1], P->l[2]);
                    printf("Absorbed in a cell with temperature: %f\n",
                            temp[idust][P->l[0]][P->l[1]][P->l[2]]);
                    printf("Re-emitted with direction: %f  %f  %f\n",
                            P->n[0], P->n[1], P->n[2]);
                    printf("Re-emitted with frequency: %e\n", P->nu);
                }
            }
            // Otherwise, scatter the photon.
            else {
                isoscatt(P, idust);
                // If we've asked for verbose output, print some info.
                if (verbose) {
                    printf("Scattering photon at cell  %i  %i  %i\n",
                            P->l[0], P->l[1], P->l[2]);
                    printf("Scattered with direction: %f  %f  %f\n",
                            P->n[0], P->n[1], P->n[2]);
                }
            }
        }
    }
}

/* Propagate a photon through the grid a distance equivalent to tau. */

void Grid::propagate_photon(Photon *P, double tau,
        bool absorb, bool bw, bool verbose, int nphot) {

    bool absorbed_by_source = false;
    int i = 0;
    while ((tau > 0) && (in_grid(P))) {
        // Calculate the distance to the next wall.
        double s1 = next_wall_distance(P, verbose);

        // Calculate how far the photon can go with the current tau.
        double alpha = 0;
        for (int idust = 0; idust<nspecies; idust++)
            alpha += P->current_kext[idust]*
                dens[idust][P->l[0]][P->l[1]][P->l[2]];

        double s2 = tau/alpha;

        // Determine whether to move to the next wall or to the end of tau.
        double s = s1;
        if (s2 < s) s = s2;

        // Calculate how far the photon can go before running into a source.
        for (int isource=0; isource<nsources; isource++) {
            double s3 = sources[isource].intercept_distance(P);

            if (s3 < s) {
                s = s3;
                absorbed_by_source = true;
            }
        }


        // Continuously absorb the photon's energy, if the end result of the
        // current trajectory is absorption.
        if (absorb) {
            for (int idust=0; idust<nspecies; idust++)
                energy[idust][P->l[0]][P->l[1]][P->l[2]] += P->energy*
                    s*P->current_kext[idust]*
                    dens[idust][P->l[0]][P->l[1]][P->l[2]];
            // If we're doing a Bjorkman & Wood simulation, update the cell to
            // find its new temperature.
            if (bw) {
                update_grid(nphot,P->l);
            }
        }

        // Remvove the tau we've used up with this stepl
        tau -= s*alpha;

        // Move the photon to it's new position.
        P->move(s);

        // If the photon moved to the next cell, update it's location.
        if (s1 < s2) P->l = photon_loc(P, verbose);
        i++;

        // If we've asked for verbose, print some information out.
        if (verbose) {
            printf("%2i  %7.4f  %i  %7.4f  %7.4f  %7.4f\n", i, tau, P->l[0],
                    P->r[0]/au, s1*P->n[0]/au, s2*P->n[0]/au);
            printf("%14i  %7.4f  %7.4f  %7.4f\n", P->l[1], P->r[1]/au, 
                    s1*P->n[1]/au, s2*P->n[1]/au);
            printf("%14i  %7.4f  %7.4f  %7.4f\n", P->l[2], P->r[2]/au, 
                    s1*P->n[2]/au, s2*P->n[2]/au);
        }

        // If the distance to the star is the shortest distance, kill the 
        // photon.
        if (absorbed_by_source) {
            P->l[0] = nw1;
            P->l[1] = nw2;
            P->l[2] = nw3;
        }

        // Kill the photon if it bounces around too many times...
        if (i > 1000) {
            tau = -1.0;
            printf("!!!!!!! ERROR - Killing photon because it seems to be stuck.\n");
        }
    }
}

/* Propagate a ray through the grid for raytracing. */

void Grid::propagate_ray(Ray *R, bool verbose) {

    int i=0;
    do {
        if (volume[R->l[0]][R->l[1]][R->l[2]] < 
                pi*R->pixel_size*R->pixel_size*R->pixel_size/6.) {
            R->pixel_too_large = true;
            break;
        }

        double s = next_wall_distance(R, verbose);

        double tau_cell = 0;
        double intensity_cell = 0;
        for (int idust=0; idust<nspecies; idust++) {
            double tau = s*R->current_kext[idust]*
                dens[idust][R->l[0]][R->l[1]][R->l[2]];

            tau_cell += tau;
            intensity_cell += (1.0-exp(-tau))*planck_function(R->nu,
                    temp[idust][R->l[0]][R->l[1]][R->l[2]]);
        }

        if (verbose) {
            printf("%2i  %7.5f  %i  %7.4f  %7.4f\n", i, tau_cell, 
                    R->l[0], R->r[0]/au, s*R->n[0]/au);
            printf("%11.1e  %i  %7.4f  %7.4f\n", R->intensity, R->l[1], 
                    R->r[1]/au, s*R->n[1]/au);
            printf("%11.5f  %i  %7.4f  %7.4f\n", R->tau, R->l[2], R->r[2]/au, 
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

double Grid::next_wall_distance(Photon *P, bool verbose) {
    return 0.0;
}

/* Calculate the distance between the photon and the outermost wall. */

double Grid::outer_wall_distance(Photon *P, bool verbose) {
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

/* Update the temperature in a cell given the number of photons that have 
 * been absorbed in the cell. */

void Grid::update_grid(int nphot, Vector<int, 3> l) {
    bool not_converged = true;

    for (int idust=0; idust<nspecies; idust++) {
        while (not_converged) {
            double T_old = temp[idust][l[0]][l[1]][l[2]];

            temp[idust][l[0]][l[1]][l[2]]=pow(energy[idust][l[0]][l[1]][l[2]]/
                (4*sigma*dust_species[idust].\
                planck_mean_opacity(temp[idust][l[0]][l[1]][l[2]])*
                mass[idust][l[0]][l[1]][l[2]]),0.25);

            // Make sure that there is a minimum temperature that the grid can
            // get to.
            if (temp[idust][l[0]][l[1]][l[2]] < 0.1) 
                temp[idust][l[0]][l[1]][l[2]] = 0.1;

            if ((fabs(T_old-temp[idust][l[0]][l[1]][l[2]])/T_old < 1.0e-2))
                not_converged = false;
        }
    }
}

void Grid::update_grid(int nphot) {
    for (int i=0; i<nw1-1; i++)
        for (int j=0; j<nw2-1; j++)
            for (int k=0; k<nw3-1; k++)
                update_grid(nphot, Vector<int, 3>(i,j,k));
}

/* Calculate the luminosity of the cell indicated by l. */

double Grid::cell_lum(Vector<int, 3> l) {
    return 4*mass[0][l[0]][l[1]][l[2]]*dust_species[0].
        planck_mean_opacity(temp[0][l[0]][l[1]][l[2]])*sigma*
        pow(temp[0][l[0]][l[1]][l[2]],4);
}

#endif
