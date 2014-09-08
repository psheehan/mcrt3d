#ifndef ISOTROPIC_DUST_CC
#define ISOTROPIC_DUST_CC

#include <stdlib.h>
#include <cmath>
#include "misc.cc"
#include "photon.cc"
#include "dust.cc"

struct IsotropicDust : public Dust {
    void scatter(Photon *P);
};

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
