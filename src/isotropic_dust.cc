#include "isotropic_dust.h"

/*IsotropicDust::IsotropicDust(int _nlam, double *_nu, double *_lam, \
        double *_kabs, double *_ksca, double *_kext, double *_albedo) {
    Dust(_nlam, _nu, _lam, _kabs, _ksca, _kext, _albedo);
}*/
        
/* Scatter a photon isotropically off of dust. */

void IsotropicDust::scatter(Photon *P) {
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
