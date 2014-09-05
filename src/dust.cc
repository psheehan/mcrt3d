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

    int ntemp;
    double *temp;
    double *planck_opacity;
    double *dplanck_opacity_dT;
    double *rosseland_extinction;
    double *drosseland_extinction_dT;

    double *int_dBnu_knu;
    double *dint_dBnu_knu_dT;
    double **Bnu;
    double **dBnu;
    double **dBnudT;
    double **ddBnudT;

    double *dkextdnu;
    double *dalbedodnu;

    void scatter(Photon *P);
    void absorb(Photon *P, double T, bool bw);

    double random_nu(double T, bool bw);
    double opacity(double freq);
    double albdo(double freq);
    double planck_mean_opacity(double T);
    double rosseland_mean_extinction(double T);
    double intdBnuknu(double T);

    double *Bnu_arr(double T);
    double *dBnu_arr(double T);

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
    double freq;

    double *F;
    double norm;
    if (bw) {
        F = dBnu_arr(T);
        norm = 1.0/intdBnuknu(T);
    }
    else {
        F = Bnu_arr(T);
        norm = -pi/(planck_mean_opacity(T)*sigma*pow(T,4));
    }

    double ksi = random_number();

    double tot = 0.0;
    for (int i=0; i<nlam-2; i++) {
        tot += 0.5*(nu[nlam-i-1]-nu[nlam-i-2])*(kext[nlam-i-1]*F[nlam-i-1]+
                kext[nlam-i-2]*F[nlam-i-2]);

        double Prob = tot*norm;

        if (Prob > ksi) {
            freq = random_number()*(nu[nlam-i-1]-nu[nlam-i-2])+
                nu[nlam-i-2];
            break;
        }
    }

    delete F;

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

/* Calculate the integral of d(B_nu)/dT * k_nu over all frequencies for a given
   temperature. */

double Dust::intdBnuknu(double T) {
    int n = find_in_arr(T,temp,ntemp);

    double int_dBnuknu = dint_dBnu_knu_dT[n]*(T-temp[n])+int_dBnu_knu[n];

    return int_dBnuknu;
}

/* Calculate the blackbody function as a function of frequency for a given
   temperature. */

double *Dust::Bnu_arr(double T) {
    double *B_nu = new double[nlam];

    int i = find_in_arr(T,temp,ntemp);

    for (int j=0; j < nlam; j++) 
        B_nu[j] = dBnudT[i][j]*(T-temp[i])+Bnu[i][j];

    return B_nu;
}

/* Calculate the derivative of the blackbody function as a function of 
 * frequency for a given temperature. */

double *Dust::dBnu_arr(double T) {
    double *dB_nu = new double[nlam];

    int i = find_in_arr(T,temp,ntemp);

    for (int j=0; j < nlam; j++) 
        dB_nu[j] = ddBnudT[i][j]*(T-temp[i])+dBnu[i][j];

    return dB_nu;
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
