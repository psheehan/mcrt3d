#include <stdio.h>
#include <cmath>
#include <omp.h>
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
    void thermal_mc_omp(int nphot, bool bw, int nthreads);
    void bw_iteration(int nphot, double ***pcount);
    void lucy_iteration(int nphot, double ***pcount);

    bool converged(double ***temp, double ***tempold, double ***tempreallyold);
    double delta(double x1, double x2);
    double quantile(double ***R, double p);
};

/* Run a Monte Carlo simulation to calculate the temperature throughout the 
 * grid. */

void MCRT::thermal_mc(int nphot, bool bw) {
    double ***pcount = new double**[G->nw1-1];
    for (int i=0; i<G->nw1-1; i++) {
        pcount[i] = new double*[G->nw2-1];
        for (int j=0; j<G->nw2-1; j++) {
            pcount[i][j] = new double[G->nw3-1];
            for (int k=0; k<G->nw3-1; k++)
                pcount[i][j][k] = 0;
        }
    }

    if (bw) 
        bw_iteration(nphot,pcount);
    else {
        double ***told = create3DArr(G->nw1-1,G->nw2-1,G->nw3-1);
        double ***treallyold = create3DArr(G->nw1-1,G->nw2-1,G->nw3-1);

        int maxniter = 20;

        for (int i=0; i<G->nw1-1; i++)
            for (int j=0; j<G->nw2-1; j++)
                for (int k=0; k<G->nw3-1; k++)
                    told[i][j][k] = G->temp[i][j][k];

        int i = 1;
        while (i <= maxniter) {
            printf("Starting iteration # %i \n\n", i);

            for (int l=0; l<G->nw1-1; l++) {
                for (int j=0; j<G->nw2-1; j++) {
                    for (int k=0; k<G->nw3-1; k++) {
                        treallyold[l][j][k] = told[l][j][k];
                        told[l][j][k] = G->temp[l][j][k];
                    }
                }
            }

            lucy_iteration(nphot,pcount);

            for (int j=0; j<G->nw1-1; j++) {
                for (int k=0; k<G->nw2-1; k++) {
                    for (int l=0; l<G->nw3-1; l++) {
                        G->update_grid(nphot,Vector<int, 3>(j,k,l),pcount);
                        pcount[j][k][l] = 0.0;
                    }
                }
            }

            if (i > 2)
                if (converged(G->temp,told,treallyold))
                    i = maxniter;

            i++;
            printf("\n");
        }

        delete told;
        delete treallyold;
    }
    delete pcount;
}

/* Run a Monte Carlo simulation to calculate the temperature throughout the 
 * grid. */

void MCRT::thermal_mc_omp(int nphot, bool bw, int nthreads) {
    double ****pcount = new double***[nthreads];
    for (int l=0; l<nthreads; l++) {
        pcount[l] = new double**[G->nw1-1];
        for (int i=0; i<G->nw1-1; i++) {
            pcount[l][i] = new double*[G->nw2-1];
            for (int j=0; j<G->nw2-1; j++) {
                pcount[l][i][j] = new double[G->nw3-1];
                for (int k=0; k<G->nw3-1; k++)
                    pcount[l][i][j][k] = 0;
            }
        }
    }

    if (bw) 
        printf("Parallel Bjorkman & Wood Simulation not yet supported.\n");
    else {
        double ***told = create3DArr(G->nw1-1,G->nw2-1,G->nw3-1);
        double ***treallyold = create3DArr(G->nw1-1,G->nw2-1,G->nw3-1);

        int maxniter = 20;

        for (int i=0; i<G->nw1-1; i++)
            for (int j=0; j<G->nw2-1; j++)
                for (int k=0; k<G->nw3-1; k++)
                    told[i][j][k] = G->temp[i][j][k];

        int i = 1;
        while (i <= maxniter) {
            printf("Starting iteration # %i \n\n", i);

            for (int l=0; l<G->nw1-1; l++) {
                for (int j=0; j<G->nw2-1; j++) {
                    for (int k=0; k<G->nw3-1; k++) {
                        treallyold[l][j][k] = told[l][j][k];
                        told[l][j][k] = G->temp[l][j][k];
                    }
                }
            }

            omp_set_num_threads(nthreads);
#pragma omp parallel
{
            lucy_iteration(nphot/nthreads,pcount[omp_get_thread_num()]);
}

            double ***pcount_tot = create3DArr(G->nw1-1,G->nw2-1,G->nw3-1);

            for (int j=0; j<G->nw1-1; j++) {
                for (int k=0; k<G->nw2-1; k++) {
                    for (int l=0; l<G->nw3-1; l++) {
                        pcount_tot[j][k][l] = 0.0;
                        for (int m=0; m<nthreads; m++) {
                            pcount_tot[j][k][l] += pcount[m][j][k][l];
                            pcount[m][j][k][l] = 0.0;

                        }
                        G->update_grid(nphot,Vector<int, 3>(j,k,l),pcount_tot);
                    }
                }
            }

            delete pcount_tot;

            if (i > 2)
                if (converged(G->temp,told,treallyold))
                    i = maxniter;

            i++;
            printf("\n");
        }

        delete told;
        delete treallyold;
    }
    delete pcount;
}

/* The Bjorkman & Wood method for calculating the temperature throughout the
 * grid. */

void MCRT::bw_iteration(int nphot, double ***pcount) {

    bool verbose = false;

    TCREATE(moo); TCLEAR(moo);
    for (int i=0; i<nphot; i++) {
        if (fmod(i+1,nphot/10) == 0) printf("%i\n",i+1);

        /*if (i == 9)
            verbose = true;
        else
            verbose = false;*/
        TSTART(moo);
        Photon *P = G->emit();
        TSTOP(moo);

        if (verbose) {
            printf("Emitting photon # %i\n", i);
            printf("Emitted with direction: %f  %f  %f\n", P->n[0], P->n[1], 
                    P->n[2]);
            printf("Emitted from a cell with temperature: %f\n", 
                    G->temp[P->l[0]][P->l[1]][P->l[2]]);
            printf("Emitted with frequency: %e\n", P->nu);
        }

        while (G->in_grid(P)) {
            double tau = -log(1-random_number());

            G->propagate_photon(P, tau, pcount, false, verbose);

            if (G->in_grid(P)) {
                if (random_number() < P->current_albedo[G->dust[P->l[0]]
                        [P->l[1]][P->l[2]]]) {
                    G->isoscatt(P);
                    if (verbose) {
                        printf("Scattering photon at cell  %i  %i  %i\n", 
                                P->l[0], P->l[1], P->l[2]);
                        printf("Scattered with direction: %f  %f  %f\n", 
                                P->n[0], P->n[1], P->n[2]);
                    }
                }
                else {
                    pcount[P->l[0]][P->l[1]][P->l[2]]++;
                    G->update_grid(nphot,P->l,pcount);
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
            }
        }
        P->clean();
        delete P;
        if (verbose) printf("Photon has escaped the grid.\n\n");
    }
    printf("%f\n", TGIVE(moo));
}

/* The Lucy method for calculating temperature throughout the grid. */

void MCRT::lucy_iteration(int nphot, double ***pcount) {

    for (int i=0; i<nphot; i++) {
        if (fmod(i+1,nphot/10) == 0) printf("%i\n",i+1);

        Photon *P = G->emit();

        while (G->in_grid(P)) {
            double tau = -log(1-random_number());

            G->propagate_photon(P, tau, pcount, true, false);

            if (G->in_grid(P)) {
                if (random_number() < P->current_albedo[G->dust[P->l[0]]
                        [P->l[1]][P->l[2]]]) {
                    G->isoscatt(P);
                }
                else {
                    G->absorb(P, true);
                }
            }
        }
        P->clean();
        delete P;
    }
}

bool MCRT::converged(double ***temp, double ***tempold, 
        double ***tempreallyold) {

    double Qthresh = 2.0;
    double Delthresh = 1.1;
    double p = 0.99;

    double ***R = new double**[G->nw1-1];
    double ***Rold = new double**[G->nw1-1];
    for (int i=0; i<G->nw1-1; i++) {
        R[i] = new double*[G->nw2-1];
        Rold[i] = new double*[G->nw2-1];
        for (int j=0; j<G->nw2-1; j++) {
            R[i][j] = new double[G->nw3-1];
            Rold[i][j] = new double[G->nw3-1];
            for (int k=0; k<G->nw3-1; k++) {
                R[i][j][k] = delta(tempold[i][j][k],temp[i][j][k]);
                Rold[i][j][k] = delta(tempreallyold[i][j][k],temp[i][j][k]);
            }
        }
    }

    double Q = quantile(R,p);
    double Qold = quantile(Rold,p);
    printf("%f   %f\n", Q, Qold);

    double Del = delta(Qold,Q);
    printf("%f\n", Del);

    bool conv = ((Q < Qthresh) && (Del < Delthresh));

    delete R;
    delete Rold;

    return conv;
}

double MCRT::delta(double x1, double x2) {
    double d1 = x1/x2;
    double d2 = x2/x1;

    double delt = d1;
    if (d2 > d1) delt = d2;

    return delt;
}

double MCRT::quantile(double ***R, double p) {
    double *Rline = new double[(G->nw1-1)*(G->nw2-1)*(G->nw3-1)];

    for (int i=0; i<G->nw1-1; i++)
        for (int j=0; j<G->nw2-1; j++)
            for (int k=0; k<G->nw3-1; k++)
                Rline[i*(G->nw2-1)*(G->nw3-1)+j*(G->nw3-1)+k] = R[i][j][k];

    bubbleSort(Rline, (G->nw1-1)*(G->nw2-1)*(G->nw3-1));

    double quant = Rline[int(p*(G->nw1-1)*(G->nw2-1)*(G->nw3-1))];

    delete Rline;

    return quant;
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

    void set_walls(Grid *G, int nw1, int nw2, int nw3, double *w1, double *w2, 
            double *w3, double *dw1, double *dw2, double *dw3) {

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

        double ***dens = pymangle(G->nw1-1, G->nw2-1, G->nw3-1, _dens);
        double ***temp = pymangle(G->nw1-1, G->nw2-1, G->nw3-1, _temp);
        double ***mass = pymangle(G->nw1-1, G->nw2-1, G->nw3-1, _mass);
        int ***dust = pymangle(G->nw1-1, G->nw2-1, G->nw3-1, _dust);

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

    void run_thermal_mc_omp(MCRT *me, int nphot, bool bw, int nthreads) {
        me->thermal_mc_omp(nphot, bw, nthreads);
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
