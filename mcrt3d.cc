#include <stdio.h>
#include <cmath>
#include "grid.cc"
#include "cartesian_grid.cc"
#include "cylindrical_grid.cc"
#include "spherical_grid.cc"
#include "pymangle.cc"
#include "timer.c"
#include "camera.cc"

struct MCRT {
    Grid *G;

    void thermal_mc(int nphot, bool bw);
    void thermal_mc_bw(int nphot);
    void thermal_mc_lucy(int nphot);
    void lucy_iteration(int nphot, double ***pcount);
    void propagate_photon_bw(Photon *P, double ***pcount, int nphot);
    void propagate_photon_lucy(Photon *P, double ***pcount, int nphot);
};

/* Run a Monte Carlo simulation to calculate the temperature throughout the 
 * grid. */

void MCRT::thermal_mc(int nphot, bool bw) {
    if (bw)
        thermal_mc_bw(nphot);
    else
        thermal_mc_lucy(nphot);
}

void MCRT::thermal_mc_bw(int nphot) {
    double ***pcount = create3DArrValue(G->n1, G->n2, G->n3, 0);
    bool verbose = false;

    for (int i=0; i<nphot; i++) {
        if (fmod(i+1,nphot/10) == 0) printf("%i\n",i+1);

        Photon *P = G->emit();

        if (verbose) {
            printf("Emitting photon # %i\n", i);
            printf("Emitted with direction: %f  %f  %f\n", P->n[0], P->n[1], 
                    P->n[2]);
            printf("Emitted from a cell with temperature: %f\n", 
                    G->temp[P->l[0]][P->l[1]][P->l[2]]);
            printf("Emitted with frequency: %e\n", P->nu);
        }

        propagate_photon_bw(P, pcount, nphot);

        P->clean();
        delete P;
        if (verbose) printf("Photon has escaped the grid.\n\n");
    }
}

void MCRT::thermal_mc_lucy(int nphot) {
    double ***pcount = create3DArrValue(G->n1, G->n2, G->n3, 0);
    double ***told = create3DArr(G->n1,G->n2,G->n3);
    double ***treallyold = create3DArr(G->n1,G->n2,G->n3);

    int maxniter = 20;

    equate3DArrs(told, G->temp, G->n1, G->n2, G->n3);

    int i = 1;
    while (i <= maxniter) {
        printf("Starting iteration # %i \n\n", i);

        equate3DArrs(treallyold, told, G->n1, G->n2, G->n3);
        equate3DArrs(told, G->temp, G->n1, G->n2, G->n3);

        lucy_iteration(nphot, pcount);

        G->update_grid(nphot, pcount);
        set3DArrValue(pcount, 0.0, G->n1, G->n2, G->n3);

        if (i > 2)
            if (converged(G->temp, told, treallyold, G->n1, G->n2, G->n3))
                i = maxniter;

        i++;
        printf("\n");
    }

    delete told;
    delete treallyold;
    delete pcount;
}

void MCRT::lucy_iteration(int nphot, double ***pcount) {
    for (int i=0; i<nphot; i++) {
        if (fmod(i+1,nphot/10) == 0) printf("%i\n",i+1);

        Photon *P = G->emit();

        propagate_photon_lucy(P, pcount, nphot);

        P->clean();
        delete P;
    }
}


/* The Bjorkman & Wood method for calculating the temperature throughout the
 * grid. */

void MCRT::propagate_photon_bw(Photon *P, double ***pcount, int nphot) {
    bool verbose = false;

    while (G->in_grid(P)) {
        double tau = -log(1-random_number());

        bool absorb = random_number() > P->current_albedo[G->dust[P->l[0]]
                [P->l[1]][P->l[2]]];

        G->propagate_photon(P, tau, pcount, absorb, true, verbose, nphot);

        if (G->in_grid(P)) {
            if (absorb) {
                //pcount[P->l[0]][P->l[1]][P->l[2]]++;
                //G->update_grid(nphot,P->l,pcount);
                G->absorb(P, false);
                if (verbose) {
                    printf("Absorbing photon at %i  %i  %i\n", P->l[0], 
                            P->l[1], P->l[2]);
                    printf("Absorbed in a cell with temperature: %f\n", 
                            G->temp[P->l[0]][P->l[1]][P->l[2]]);
                    printf("Re-emitted with direction: %f  %f  %f\n", 
                            P->n[0], P->n[1], P->n[2]);
                    printf("Re-emitted with frequency: %e\n", P->nu);
                }
            }
            else {
                G->isoscatt(P);
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

/* The Lucy method for calculating temperature throughout the grid. */

void MCRT::propagate_photon_lucy(Photon *P, double ***pcount, int nphot) {
    while (G->in_grid(P)) {
        double tau = -log(1-random_number());

        bool absorb = random_number() > P->current_albedo[G->dust[P->l[0]]
                                [P->l[1]][P->l[2]]];

        G->propagate_photon(P, tau, pcount, absorb, false, false, nphot);

        if (G->in_grid(P)) {
            if (absorb) {
                G->absorb(P, true);
            }
            else {
                G->isoscatt(P);
            }
        }
    }
}

/* Functions to link to Python. */

extern "C" {
    MCRT* new_mcrt() {
        return new MCRT();
    }

    void set_Grid(MCRT *me, Grid *G) {
        me->G = G;
    }

    /* Functions to set up the grid. */

    CartesianGrid* new_CartesianGrid() {
        return new CartesianGrid();
    }

    CylindricalGrid* new_CylindricalGrid() {
        return new CylindricalGrid();
    }

    SphericalGrid* new_SphericalGrid() {
        return new SphericalGrid();
    }

    void set_walls(Grid *G, int n1, int n2, int n3, int nw1, int nw2, int nw3, 
            double *w1, double *w2, double *w3, double *dw1, double *dw2, 
            double *dw3) {

        G->n1 = n1;
        G->n2 = n2;
        G->n3 = n3;
        G->nw1 = nw1;
        G->nw2 = nw2;
        G->nw3 = nw3;
        G->w1 = w1;
        G->w2 = w2;
        G->w3 = w3;
        G->dw1 = dw1;
        G->dw2 = dw2;
        G->dw3 = dw3;
    }

    void set_physical_properties(Grid *G, double *_dens, double *_temp,
            double *_mass, int *_dust) {

        double ***dens = pymangle(G->n1, G->n2, G->n3, _dens);
        double ***temp = pymangle(G->n1, G->n2, G->n3, _temp);
        double ***mass = pymangle(G->n1, G->n2, G->n3, _mass);
        int ***dust = pymangle(G->n1, G->n2, G->n3, _dust);

        G->dens = dens;
        G->temp = temp;
        G->mass = mass;
        G->dust = dust;
    }

    void create_dust_species_array(Grid *G, int nspecies) {
        G->nspecies = nspecies;

        G->dust_species = new Dust[nspecies];
    }

    void set_dust_species(Grid *G, Dust *D, int index) {
        G->dust_species[index] = *D;
    }

    void create_sources_array(Grid *G, int nsources) {
        G->sources = new Source[nsources];
    }

    void set_sources(Grid *G, Source *S, int index) {
        G->sources[index] = *S;
    }

    /* Functions to set up the dust. */

    Dust* new_Dust() {
        return new Dust();
    }

    void set_optical_properties(Dust *D, int nlam, double *nu, double *lam, 
            double *kabs, double *ksca, double *kext, double *albedo) {
        
        D->nlam = nlam;
        D->nu = nu;
        D->lam = lam;
        D->kabs = kabs;
        D->ksca = ksca;
        D->kext = kext;
        D->albedo = albedo;
    }

    void set_lookup_tables(Dust *D, int ntemp, double *temp, 
            double *planck_opacity, double *int_dBnu_knu, 
            double *dplanck_opacity_dT, double *dint_dBnu_knu_dT,
            double *dkextdnu, double *dalbedodnu, double *_Bnu, double *_dBnu,
            double *_dBnudT, double *_ddBnudT, double *_int_Bnu_knu_nu,
            double *_int_dBnu_knu_nu) {

        D->ntemp = ntemp;
        D->temp = temp;
        D->planck_opacity = planck_opacity;
        D->int_dBnu_knu = int_dBnu_knu;
        D->dplanck_opacity_dT = dplanck_opacity_dT;
        D->dint_dBnu_knu_dT = dint_dBnu_knu_dT;
        D->dkextdnu = dkextdnu;
        D->dalbedodnu = dalbedodnu;

        double **Bnu = pymangle(D->ntemp, D->nlam, _Bnu);
        double **dBnu = pymangle(D->ntemp, D->nlam, _dBnu);
        double **dBnudT = pymangle(D->ntemp, D->nlam, _dBnudT);
        double **ddBnudT = pymangle(D->ntemp, D->nlam, _ddBnudT);
        double **int_Bnu_knu_nu = pymangle(D->ntemp, D->nlam, _int_Bnu_knu_nu);
        double **int_dBnu_knu_nu = pymangle(D->ntemp, D->nlam, 
                _int_dBnu_knu_nu);

        D->Bnu = Bnu;
        D->dBnu = dBnu;
        D->dBnudT = dBnudT;
        D->ddBnudT = ddBnudT;
        D->int_Bnu_knu_nu = int_Bnu_knu_nu;
        D->int_dBnu_knu_nu = int_dBnu_knu_nu;
    }

    /* Functions to set up the sources. */

    Source* new_Source() {
        return new Source();
    }

    void set_parameters(Source *S, double x, double y, double z, double mass, 
            double radius, double temperature) {

        S->r[0] = x;
        S->r[1] = y;
        S->r[2] = z;
        S->mass = mass;
        S->radius = radius;
        S->temperature = temperature;
    }

    void set_blackbody_spectrum(Source *S, int nnu, double *nu, double *Bnu, 
            double luminosity) {

        S->nnu = nnu;
        S->nu = nu;
        S->Bnu = Bnu;
        S->luminosity = luminosity;
    }

    /* Function to test that everything is being set up correctly. */

    void test(Grid *G) {
        printf("%i\n", G->nw1);
    }

    /* Function to run a thermal monte carlo simulation. */

    void run_thermal_mc(MCRT *me, int nphot, bool bw) {
        me->thermal_mc(nphot, bw);
    }

    /* Function to create an image. */

    Image *new_Image(double r, double incl, double pa, double *x, double *y, 
            double *_intensity, int nx, int ny, double *nu, double pixel_size, 
            int nnu) {

        double ***intensity = pymangle(nx, ny, nnu, _intensity);

        return new Image(r, incl, pa, x, y, intensity, nx, ny, nu, 
                pixel_size, nnu);
    }

    /* Functions to operate the camera. */

    Camera *new_Camera() {
        return new Camera();
    }

    void set_camera_grid(Camera *C, Grid *G) {
        C->G = G;
    }

    void make_image(Camera *C, Image *I) {
        C->image = I;

        C->make_image();
    }
}
