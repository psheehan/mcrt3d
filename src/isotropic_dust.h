#ifndef ISOTROPIC_DUST_H
#define ISOTROPIC_DUST_H

#include <stdlib.h>
#include <cmath>
#include "misc.h"
#include "photon.h"
#include "dust.h"

struct IsotropicDust : public Dust {
    void scatter(Photon *P);
};

#endif
