#include "photon.h"

/* Clean up the photon to remove any pointers. */

Photon::~Photon() {
    delete[] current_kext; delete[] current_albedo;
}

/* Move the photon a distance s along its direction vector. */

void Photon::move(double s) {
    r += s*n;
}
