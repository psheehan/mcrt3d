#ifndef DUST_H
#define DUST_H

#include <stdlib.h>
#include <cmath>
#include "pymangle.h"
#include "misc.h"
#include "photon.h"

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

    Dust(int _nlam, double *_nu, double *_lam, double *_kabs, double *_ksca, 
            double *_kext, double *_albedo);
        
    void set_lookup_tables(int _ntemp, double *_temp, 
            double *_planck_opacity, double *_rosseland_extinction, 
            double *_dplanck_opacity_dT, double *_drosseland_extinction_dT,
            double *_dkextdnu, double *dalbedodnu, double *_random_nu_CPD, 
            double *_random_nu_CPD_bw, double *_drandom_nu_CPD_dT, 
            double *_drandom_nu_CPD_bw_dT);

    virtual void scatter(Photon *P);
    void absorb(Photon *P, double T, bool bw);

    double random_nu(double T, bool bw);
    double opacity(double freq);
    double albdo(double freq);
    double planck_mean_opacity(double T);
    double rosseland_mean_extinction(double T);
};

#endif
