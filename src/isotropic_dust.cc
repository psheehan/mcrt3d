#include "isotropic_dust.h"

void IsotropicDust::scatter(Photon *P) {
    double cost = -1+2*random_number(random_pool);
    double sint = sqrt(1-pow(cost,2));
    double phi = 2*pi*random_number(random_pool);

    P->n[0] = sint*cos(phi);
    P->n[1] = sint*sin(phi);
    P->n[2] = cost;
    P->invn[0] = 1.0/P->n[0];
    P->invn[1] = 1.0/P->n[1];
    P->invn[2] = 1.0/P->n[2];
}
