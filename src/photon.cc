#include "photon.h"

/* Move the photon a distance s along its direction vector. */

void Photon::move(double s) {
    r += s*n;
}

/* Clean up the photon to remove any pointers. */

void Photon::clean() {
    delete current_kext;
    delete current_albedo;
}
