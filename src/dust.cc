#ifndef DUST_CC
#define DUST_CC

#include <stdlib.h>
#include <cmath>
#include "misc.cc"
#include "photon.cc"

struct Dust {
    int nlam;
    double *nu;
    double *lam;
    double *kabs;
    double *ksca;
    double *kext;
    double *albedo;
    double *dkextdnu;
    double *dalbedodnu;

    int ntemp;
    double *temp;
    double *planck_opacity;
    double *dplanck_opacity_dT;
    double *rosseland_extinction;
    double *drosseland_extinction_dT;
    double **random_nu_CPD;
    double **random_nu_CPD_bw;
    double **drandom_nu_CPD_dT;
    double **drandom_nu_CPD_bw_dT;

    void scatter(Photon *P);
    void absorb(Photon *P, double T, bool bw);

    double random_nu(double T, bool bw);
    double opacity(double freq);
    double albdo(double freq);
    double planck_mean_opacity(double T);
    double rosseland_mean_extinction(double T);

    int find_freq_bin(double freq);
};

/* Scatter a photon isotropically off of dust. */

void Dust::scatter(Photon *P) {
    double cost = -1+2*random_number();
    double sint = sqrt(1-pow(cost,2));
    double phi = 2*pi*random_number();

    P->n[0] = sint*cos(phi);
    P->n[1] = sint*sin(phi);
    P->n[2] = cost;
    P->invn[0] = 1.0/P->n[0];
    P->invn[1] = 1.0/P->n[1];
    P->invn[2] = 1.0/P->n[2];
}

/* Absorb and then re-emit a photon from dust. */

void Dust::absorb(Photon *P, double T, bool bw) {
    double cost = -1+2*random_number();
    double sint = sqrt(1-pow(cost,2));
    double phi = 2*pi*random_number();

    P->n[0] = sint*cos(phi);
    P->n[1] = sint*sin(phi);
    P->n[2] = cost;
    P->invn[0] = 1.0/P->n[0];
    P->invn[1] = 1.0/P->n[1];
    P->invn[2] = 1.0/P->n[2];

    P->nu = random_nu(T,bw);
}

/* Calculate a random frequency for a photon. */

double Dust::random_nu(double T, bool bw) {
    double freq, CPD;

    int i = find_in_arr(T,temp,ntemp);

    double ksi = random_number();

    for (int j=0; j < nlam; j++) {
        if (bw)
            CPD = drandom_nu_CPD_bw_dT[i][j] * (T - temp[i]) + 
                random_nu_CPD_bw[i][j];
        else
            CPD = drandom_nu_CPD_dT[i][j] * (T - temp[i]) + 
                random_nu_CPD[i][j];

        if (CPD > ksi) {
            freq = random_number() * (nu[j+1] - nu[j]) + nu[j];
            break;
        }
    }

    return freq;
}

/* Calculate the opacity of a dust grain at a specific frequency. */

double Dust::opacity(double freq) {
    int l = find_freq_bin(freq);

    double opacity = dkextdnu[l]*(freq-nu[l])+kext[l];

    return opacity;
};

/* Calculate the albedo of a dust grain at a specific frequency. */

double Dust::albdo(double freq) {
    int l = find_freq_bin(freq);

    double albdo = dalbedodnu[l]*(freq-nu[l])+albedo[l];

    return albdo;
}

/* Calculate the Planck Mean Opacity for a dust grain at a given temperature. */

double Dust::planck_mean_opacity(double T) {
    int n = find_in_arr(T,temp,ntemp);

    double planck_mean_opacity = dplanck_opacity_dT[n]*(T-temp[n])+
        planck_opacity[n];

    return planck_mean_opacity;
}

/* Calculate the Rosseland Mean Extinction for a dust grain at a given 
 * temperature. */

double Dust::rosseland_mean_extinction(double T) {
    int n = find_in_arr(T,temp,ntemp);

    double rosseland_mean_extinction = drosseland_extinction_dT[n]*(T-temp[n])+
        rosseland_extinction[n];

    return rosseland_mean_extinction;
}

/* Determine which frequency bin a given frequency is in. */

int Dust::find_freq_bin(double freq) {
    int lmin = 0;
    int lmax = nlam-1;
    int l = 0;
    bool not_found = true;

    while (not_found) {
        int ltest = (lmax - lmin)/2 + lmin;

        if ((freq > nu[ltest+1]) && (freq <= nu[ltest])) {
            l = ltest;
            not_found = false;
        }
        else {
            if (freq < nu[ltest])
                lmin = ltest;
            else
                lmax = ltest;
        }
    }

    return l;
};

#endif
