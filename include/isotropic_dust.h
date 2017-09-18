#ifndef ISOTROPIC_DUST_H
#define ISOTROPIC_DUST_H

#include <stdlib.h>
#include <cmath>
#include "misc.h"
#include "photon.h"
#include "dust.h"

struct IsotropicDust : public Dust {
    IsotropicDust(int _nlam, double *_nu, double *_lam, double *_kabs, 
            double *_ksca, double *_kext, double *_albedo) : 
        Dust(_nlam, _nu, _lam, _kabs, _ksca, _kext, _albedo) {};
        
    void scatter(Photon *P);
};

#endif
