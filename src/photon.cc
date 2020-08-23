#include "photon.h"

/* Clean up the photon to remove any pointers. */

Photon::~Photon() {
    delete[] current_kext; delete[] current_albedo;
}

/* Move the photon a distance s along its direction vector. */

void Photon::move(double s) {
    r += s*n;
}

/* Clean up the ray to properly remove pointers. */

Ray::~Ray() {
    for (int i = 0; i < ndust; i++) {
        delete[] current_kext[i]; delete[] current_albedo[i];
    }
    delete[] current_kext; delete[] current_albedo;
    delete[] intensity; delete[] tau;
}

