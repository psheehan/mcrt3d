#include "photon.h"

/* Clean up the photon to remove any pointers. */

Photon::~Photon() {
}

/* Move the photon a distance s along its direction vector. */

void Photon::move(double s) {
    r += s*n;
}

/* Clean up the ray to properly remove pointers. */

Ray::~Ray() {
}

