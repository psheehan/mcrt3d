#include "cylindrical_grid.h"

CylindricalGrid::CylindricalGrid(py::array_t<double> _r, 
            py::array_t<double> _phi, py::array_t<double> _z) : 
        Grid(_r, _phi, _z) {

    // Set up the x, y, and z arrays.

    r = py::array_t<double>(n1);
    phi = py::array_t<double>(n2);
    z = py::array_t<double>(n3);

    // Load the array buffers to get the proper setup info.

    auto r_buf = r.request(); auto phi_buf = phi.request(); 
    auto z_buf = z.request();

    // Now get the correct values.

    double *__r = (double *) r_buf.ptr;
    for (int i = 0; i < n1; i++) __r[i] = 0.5 * (w1[i+1] + w1[i]);

    double *__phi = (double *) phi_buf.ptr;
    for (int i = 0; i < n2; i++) __phi[i] = 0.5 * (w2[i+1] + w2[i]);

    double *__z = (double *) z_buf.ptr;
    for (int i = 0; i < n3; i++) __z[i] = 0.5 * (w3[i+1] + w3[i]);
    
    // Set up the volume of each cell.

    auto _volume_buf = _volume.request();
    double *__volume = (double *) _volume_buf.ptr;

    for (int i = 0; i < n1; i++)
        for (int j = 0; j < n2; j++)
            for (int k = 0; k < n3; k++)
                __volume[i*n2*n3 + j*n3 + k] = (w1[i+1] - w1[i]) * 
                    (w2[j+1] - w2[j]) * (w3[k+1] - w3[k]);

    // Finally, set up the mirror symmetry.

     mirror_symmetry = false;
}

/* Calculate the distance between the photon and the nearest wall. */

double CylindricalGrid::next_wall_distance(Photon *P) {

    //double r = sqrt(P->r[0]*P->r[0]+P->r[1]*P->r[1]);
    double r = P->rad;

    // Calculate the distance to the intersection with the next radial wall.
    
    double a = P->n[0]*P->n[0]+P->n[1]*P->n[1];
    double b = P->r[0]*P->n[0]+P->r[1]*P->n[1];

    double s = HUGE_VAL;
    for (int i=P->l[0]; i <= P->l[0]+1; i++) {
        if (r == w1[i]) {
            double sr1 = (-b + fabs(b))/a;
            if ((sr1 < s) && (sr1 > 0) && (not equal_zero(sr1/
                    (P->rad*(w2[P->l[1]+1]-w2[P->l[1]])),EPSILON))) s = sr1;
            double sr2 = (-b - fabs(b))/a;
            if ((sr2 < s) && (sr2 > 0) && (not equal_zero(sr2/
                    (P->rad*(w2[P->l[1]+1]-w2[P->l[1]])),EPSILON))) s = sr2;
        }
        else {
            double c = r*r - w1[i]*w1[i];
            double d = b*b - a*c;

            if (d >= 0) {
                double sr1 = (-b + sqrt(d))/a;
                if ((sr1 < s) && (sr1 > 0)) s = sr1;
                double sr2 = (-b - sqrt(d))/a;
                if ((sr2 < s) && (sr2 > 0)) s = sr2;
            }
        }
    }

    // Calculate the distance to the intersection with the next phi wall.
    
    if (nw2 != 2) {
        //double phi = fmod(atan2(P->r[1],P->r[0])+2*pi,2*pi);
        double phi = P->phi;

        for (int i=P->l[1]; i <= P->l[1]+1; i++) {
            if (phi != w2[i]) {
                double c = P->r[0]*sin(w2[i])-P->r[1]*cos(w2[i]);
                double d = P->n[0]*sin(w2[i])-P->n[1]*cos(w2[i]);

                double sp = -c/d;

                if ((sp < s) && (sp > 0)) s = sp;
            }
        }
    }

    // Calculate the distance to intersection with the nearest z wall.
    
    double sz1 = (w3[P->l[2]] - P->r[2])*P->invn[2];
    if ((sz1 < s) && (sz1 > 0)) s = sz1;
    double sz2 = (w3[P->l[2]+1] - P->r[2])*P->invn[2];
    if ((sz2 < s) && (sz2 > 0)) s = sz2;
    
    return s;
}

/* Calculate the distance between the photon and the outermost wall. */

double CylindricalGrid::outer_wall_distance(Photon *P) {
    double s = 0;

    double r = sqrt(P->r[0]*P->r[0]+P->r[1]*P->r[1]);

    // Calculate the distance to the intersection with the next radial wall.
    
    if (r >= w1[nw1-1]) {
        double sr = HUGE_VAL;

        double a = P->n[0]*P->n[0]+P->n[1]*P->n[1];
        double b = P->r[0]*P->n[0]+P->r[1]*P->n[1];
        double c = r*r - w1[nw1-1]*w1[nw1-1];
        double d = b*b - a*c;

        if (d >= 0) {
            double sr1 = (-b + sqrt(d))/a;
            if ((sr1 < sr) && (sr1 > 0)) sr = sr1;
            double sr2 = (-b - sqrt(d))/a;
            if ((sr2 < sr) && (sr2 > 0)) sr = sr2;

            if (sr > s) s = sr;
        }
    }

    // Calculate the distance to intersection with the nearest z wall.
    
    if (P->n[2] != 0) {
        if (mirror_symmetry) {
            if (P->r[2] <= -w3[nw3-1]) {
                double sz = (-w3[nw3-1] - P->r[2])*P->invn[2];
                if (sz > s) s = sz;
            }
            else if (P->r[2] >= w3[nw3-1]) {
                double sz = (w3[nw3-1] - P->r[2])*P->invn[2];
                if (sz > s) s = sz;
            }
        }
        else {
            if (P->r[2] <= w3[0]) {
                double sz = (w3[0] - P->r[2])*P->invn[2];
                if (sz > s) s = sz;
            }
            else if (P->r[2] >= w3[nw3-1]) {
                double sz = (w3[nw3-1] - P->r[2])*P->invn[2];
                if (sz > s) s = sz;
            }
        }
    }

    Vector<double, 3> newr = P->r + s*P->n;
    double newtwodr = sqrt(newr[0]*newr[0] + newr[1]*newr[1]);

    if (Q->verbose) printf("%20.17f   %7.4f   %7.4f\n", newr[0]/au, newr[1]/au, 
            newr[2]/au);
    if (Q->verbose) printf("%20.17f\n", newtwodr/au);

    if (equal(newtwodr, w1[nw1-1], EPSILON)) newtwodr = w1[nw1-1];
    if (mirror_symmetry) {
        if (equal(newr[2],w3[0],EPSILON)) newr[2] = w3[0];
        else if (equal(newr[2],w3[nw3-1],EPSILON)) newr[2] = w3[nw3-1];
    }
    else {
        if (equal(newr[2],-w3[nw3-1],EPSILON)) newr[2] = -w3[nw3-1];
        else if (equal(newr[2],w3[nw3-1],EPSILON)) newr[2] = w3[nw3-1];
    }

    if (Q->verbose) printf("%20.17f   %7.4f   %7.4f\n", newr[0]/au, newr[1]/au, 
            newr[2]/au);
    if (Q->verbose) printf("%20.17f\n", newtwodr/au);

    if (mirror_symmetry) {
        if ((newr[2] < -w3[nw3-1]) || (newr[2] > w3[nw3-1]) || 
                (newtwodr > w1[nw1-1])) {
            s = HUGE_VAL;
        }
    }
    else {
        if ((newr[2] < w3[0]) || (newr[2] > w3[nw3-1]) || 
                (newtwodr > w1[nw1-1])) {
            s = HUGE_VAL;
        }
    }

    return s;
}

/* Calculate the smallest absolute distance to the nearest wall. */

double CylindricalGrid::minimum_wall_distance(Photon *P) {
    double r = P->rad;

    // Calculate the distance to the nearest radial wall.
    
    double s = HUGE_VAL;
    for (int i=P->l[0]; i <= P->l[0]+1; i++) {
        double sr = fabs(r - w1[i]);
        if (sr < s) s = sr;
    }

    // Calculate the distance to the nearest phi wall.
    
    if (nw2 != 2) {
        double phi = P->phi;

        for (int i=P->l[1]; i <= P->l[1]+1; i++) {
            Vector<double, 3> phi_hat = Vector<double, 3>(-sin(w2[i]), 
                    cos(w2[i]), 0);

            double sp = fabs(phi_hat * P->r);
            if (sp < s) s = sp;
        }
    }

    // Calculate the distance to the nearest z wall.
    
    double sz1 = fabs(w3[P->l[2]] - P->r[2]);
    if (sz1 < s) s = sz1;
    double sz2 = fabs(w3[P->l[2]+1] - P->r[2]);
    if (sz2 < s) s = sz2;
    
    return s;
}

/* Calculate the size of the grid. */

double CylindricalGrid::grid_size() {
    double rw1_max = w1[nw1-1];
    double rw3_max = std::max(abs(w3[0]), abs(w3[nw3-1]));

    return 2*sqrt(rw1_max*rw1_max + rw3_max*rw3_max);
}

/* Determine which cell the photon is in. */

Vector<int, 3> CylindricalGrid::photon_loc(Photon *P) {
    Vector<int, 3> l;

    /* Find the location in the radial grid. */
    
    double pi = 3.14159265;
    //double r = sqrt(P->r[0]*P->r[0]+P->r[1]*P->r[1]);
    //double phi = fmod(atan2(P->r[1],P->r[0])+2*pi,2*pi);
    P->rad = sqrt(P->r[0]*P->r[0]+P->r[1]*P->r[1]);
    /* If P->rad = 0 we need to be careful about how we calculate phi. */
    if (P->rad == 0) {
        P->phi = fmod(P->phi+pi,2*pi);
        P->l[1] = -1;
    }
    else
        P->phi = fmod(atan2(P->r[1],P->r[0])+2*pi,2*pi);
    double r = P->rad;
    double phi = P->phi;

    double gnx = cos(phi);
    double gny = sin(phi);
    if (equal_zero(gnx, EPSILON)) gnx = 0.;
    if (equal_zero(gny, EPSILON)) gny = 0.;
    
    // Check if we are using mirror symmetry and we're in the southern
    // hemisphere. If so, we need to flip.
    
    if (mirror_symmetry) {
        if (P->r[2] < 0) {
            P->r[2] *= -1;
            P->n[2] *= -1;
        }

        if (equal_zero(P->r[2], EPSILON) and P->n[2] < 0)
            P->n[2] *= -1;
    }

    /* Determine which cell the photon is currently in. */

    if (r >= w1[nw1-1])
        l[0] = n1-1;
    else if (r <= w1[0])
        l[0] = 0;
    else {
        if (P->l[0] == -1)
            l[0] = find_in_arr(r,w1,nw1);
        else {
            int lower = P->l[0]-1;
            if (lower < 0) lower = 0;
            int upper = P->l[0]+1;
            if (upper > n1-1) upper = n1;
            
            l[0] = find_in_arr(r,w1,lower,upper);
        }
    }

    /* Because of floating point errors it may be the case that the photon 
     * should be on the wall exactly, but is not exactly on the wall. We
     * need to put the photon exactly on the wall. */

    if (equal(r,w1[l[0]],EPSILON))
        r = w1[l[0]];
    else if (equal(r,w1[l[0]+1],EPSILON))
        r = w1[l[0]+1];

    /* Finally, update which cell the photon is in based on the direction it
     * is going. */

    double nr = P->n[0]*gnx + P->n[1]*gny;
    if (equal_zero(nr, EPSILON)) nr = 0.0;
    
    if ((r == w1[l[0]]) && (nr < 0))
        l[0] -= 1;
    else if ((r == w1[l[0]+1]) && (nr >= 0))
        l[0] += 1;

    // Find the location in the phi grid.
    
    if (nw2 == 2)
        l[1] = 0;
    else {
        if (P->l[1] == -1)
            l[1] = find_in_arr(phi,w2,nw2);
        else
            l[1] = find_in_periodic_arr(phi,w2,n2,P->l[1]-1,P->l[1]+1);

        /* Check whether the photon is supposed to be exactly on the cell
         * wall. Floating point errors may keep it from being exactly on the
         * wall, and we need to fix that. */

        if (equal(phi,w2[l[1]],EPSILON))
            phi = w2[l[1]];
        else if (equal(phi,w2[l[1]+1],EPSILON))
            phi = w2[l[1]+1];

        /* Update which cell the photon is in depending on the 
         * direction it is going. */

        double gnx = -sin(phi);
        double gny = cos(phi);
        if (equal_zero(gnx, EPSILON)) gnx = 0.;
        if (equal_zero(gny, EPSILON)) gny = 0.;
        
        if ((phi == w2[l[1]]) && (P->n[0]*gnx+P->n[1]*gny <= 0))
            l[1] -= 1;
        else if ((phi == w2[l[1]+1]) && (P->n[0]*gnx+P->n[1]*gny >= 0))
            l[1] += 1;
        l[1] = (l[1]+n2)%(n2);

        /* Finally, if you are at phi = 0, but going towards negative phi, 
         * you should set phi = 2*pi. */

        if ((phi == 0) && (l[1] == n2-1))
            phi = w2[l[1]+1];
    }

    /* Since we may have updated r and phi to be exactly on the grid cell
     * walls, change the photon position slightly to reflect this. */

    P->r[0] = r * cos(phi);
    P->r[1] = r * sin(phi);
    P->rad = r;
    P->phi = phi;

    // Find the location in the z grid.
    
    if (P->r[2] >= w3[nw3-1])
        l[2] = n3-1;
    else if (P->r[2] <= w3[0])
        l[2] = 0;
    else {
        if (P->l[2] == -1)
            l[2] = find_in_arr(P->r[2],w3,nw3);
        else {
            int lower = P->l[2]-1;
            if (lower < 0) lower = 0;
            int upper = P->l[2]+1;
            if (upper > nw3-1) upper = nw3-1;
            
            l[2] = find_in_arr(P->r[2],w3,lower,upper);
        }
    }
    
    if ((P->r[2] == w3[l[2]]) && (P->n[2] < 0))
        l[2] -= 1;
    else if ((P->r[2] == w3[l[2]+1]) && (P->n[2] > 0))
        l[2] += 1;
    
    return l;
}

/* Check whether a photon is on a wall and going parallel to it. */

bool CylindricalGrid::on_and_parallel_to_wall(Photon *P) {
    if (P->r[2] == w3[P->l[2]] and P->n[2] == 0)
        return true;
    if (P->r[2] == w3[P->l[2]+1] and P->n[2] == 0)
        return true;

    return false;
}

/* Randomly generate a photon location within a cell. */
 
Vector<double, 3> CylindricalGrid::random_location_in_cell(int ix, int iy, 
        int iz) {
    double r = w1[ix] + random_number() * (w1[ix+1] - w1[ix]);
    double phi = w2[iy] + random_number() * (w2[iy+1] - w1[iy]);
    double z = w3[iz] + random_number() * (w3[iz+1] - w1[iz]);

    double x = r * cos(phi);
    double y = r * sin(phi);

    return Vector<double, 3>(x, y, z);
}

/* Check whether a photon is in the boundaries of the grid. */

bool CylindricalGrid::in_grid(Photon *P) {

    if ((P->l[0] >= n1) || (P->l[0] < 0) || (P->l[2] >= n3) || (P->l[2] < 0))
        return false;
    else
        return true;
}
