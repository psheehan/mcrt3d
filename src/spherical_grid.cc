#include "spherical_grid.h"

SphericalGrid::SphericalGrid(py::array_t<double> _r, 
            py::array_t<double> _theta, py::array_t<double> _phi) : 
        Grid(_r, _theta, _phi) {

    // Set up the x, y, and z arrays.

    r = py::array_t<double>(n1);
    theta = py::array_t<double>(n2);
    phi = py::array_t<double>(n3);

    // Load the array buffers to get the proper setup info.

    auto r_buf = r.request(); auto theta_buf = theta.request(); 
    auto phi_buf = phi.request();

    // Now get the correct values.

    double *__r = (double *) r_buf.ptr;
    for (int i = 0; i < n1; i++) __r[i] = 0.5 * (w1[i+1] + w1[i]);

    double *__theta = (double *) theta_buf.ptr;
    for (int i = 0; i < n2; i++) __theta[i] = 0.5 * (w2[i+1] + w2[i]);

    double *__phi = (double *) phi_buf.ptr;
    for (int i = 0; i < n3; i++) __phi[i] = 0.5 * (w3[i+1] + w3[i]);
    
    // Check for mirror symmetry.

    int volume_scale = 1;
    if (equal_zero(cos(w2[nw2-1]), EPSILON))
    {
        mirror_symmetry = true;
        volume_scale = 2;
    }
    else
        mirror_symmetry = false;

    // Set up the volume of each cell.

    auto _volume_buf = _volume.request();
    double *__volume = (double *) _volume_buf.ptr;

    for (int i = 0; i < n1; i++)
        for (int j = 0; j < n2; j++)
            for (int k = 0; k < n3; k++)
                __volume[i*n2*n3 + j*n3 + k] = (w1[i+1]*w1[i+1]*w1[i+1] - 
                        w1[i]*w1[i]*w1[i]) * (cos(w2[j]) - cos(w2[j+1])) * 
                        (w3[k+1] - w3[k]) / 3 * volume_scale;
}

/* Calculate the distance between the photon and the nearest wall. */

double SphericalGrid::next_wall_distance(Photon *P) {

    //double r = P->r.norm();
    double r = P->rad;

    // Calculate the distance to the intersection with the next radial wall.
    
    double b = P->r*P->n;

    double s = HUGE_VAL;
    for (int i=P->l[0]; i <= P->l[0]+1; i++) {
        if (r == w1[i]) {
            double sr1 = -b + fabs(b);
            if ((sr1 < s) && (sr1 > 0) && (not equal_zero(sr1/
                    (P->rad*(w2[P->l[1]+1]-w2[P->l[1]])),EPSILON))) s = sr1;
            double sr2 = -b - fabs(b);
            if ((sr2 < s) && (sr2 > 0) && (not equal_zero(sr2/
                    (P->rad*(w2[P->l[1]+1]-w2[P->l[1]])),EPSILON))) s = sr2;
        }
        else {
            double c = r*r - w1[i]*w1[i];
            double d = b*b - c;

            if (d >= 0) {
                double sr1 = -b + sqrt(d);
                if ((sr1 < s) && (sr1 > 0)) s = sr1;
                double sr2 = -b - sqrt(d);
                if ((sr2 < s) && (sr2 > 0)) s = sr2;
            }
        }
    }

    // Calculate the distance to the intersection with the next theta wall.
    
    if (nw2 != 2) {
        double theta = P->theta;
        
        for (int i=P->l[1]; i <= P->l[1]+1; i++) {
            if (equal_zero(cos(w2[i]),1.0e-10)) {
                double st1 = -P->r[2]*P->invn[2];
                if (equal_zero(st1/(P->rad*(w2[P->l[1]+1]-w2[P->l[1]])),
                        EPSILON)) st1 = 0;
                if ((st1 < s) && (st1 > 0)) s = st1;
            }
            else {
                double a = P->n[0]*P->n[0]+P->n[1]*P->n[1]-P->n[2]*P->n[2]*
                    pow(tan(w2[i]),2);
                double b = 2*(P->r[0]*P->n[0]+P->r[1]*P->n[1]-P->r[2]*P->n[2]*
                    pow(tan(w2[i]),2));

                //if (theta == w2[i]) {
                if (equal(sin(theta),sin(w2[i]),1.0e-10)) {
                    double st1 = (-b + fabs(b))/(2*a);
                    if ((st1 < s) && (st1 > 0)) s = st1;
                    double st2 = (-b - fabs(b))/(2*a);
                    if ((st2 < s) && (st2 > 0)) s = st2;
                }
                else {
                    double c = P->r[0]*P->r[0]+P->r[1]*P->r[1]-P->r[2]*P->r[2]*
                        pow(tan(w2[i]),2);
                    double d = b*b-4*a*c;

                    if (d >= 0) {
                        double st1 = (-b + sqrt(d))/(2*a);
                        if ((st1 < s) && (st1 > 0)) s = st1;
                        double st2 = (-b - sqrt(d))/(2*a);
                        if ((st2 < s) && (st2 > 0)) s = st2;
                    }
                }
            }
        }
    }

    // Calculate the distance to intersection with the nearest phi wall.
    
    if (nw3 != 2) {
        double phi = P->phi;
        
        for (int i=P->l[2]; i <= P->l[2]+1; i++) {
            if (phi != w3[i]) {
                double c = P->r[0]*sin(w3[i])-P->r[1]*cos(w3[i]);
                double d = P->n[0]*sin(w3[i])-P->n[1]*cos(w3[i]);

                double sp = -c/d;

                if ((sp < s) && (sp > 0)) s = sp;
            }
        }
    }
    
    return s;
}

/* Calculate the distance between the photon and the outermost wall. */

double SphericalGrid::outer_wall_distance(Photon *P) {

    double r = P->r.norm();

    double s = HUGE_VAL;

    double b = P->r*P->n;
    double c = r*r - w1[nw1-1]*w1[nw1-1];
    double d = b*b - c;

    if (d >= 0) {
        double sr1 = -b + sqrt(d);
        if ((sr1 < s) && (sr1 > 0)) s = sr1;
        double sr2 = -b - sqrt(d);
        if ((sr2 < s) && (sr2 > 0)) s = sr2;
    }

    return s;
}

/* Calculate the smallest absolute distance to the nearest wall. */

double SphericalGrid::minimum_wall_distance(Photon *P) {
    double r = P->rad;

    // Calculate the distance to the nearest radial wall.
    
    double s = HUGE_VAL;
    for (int i=P->l[0]; i <= P->l[0]+1; i++) {
        double sr = fabs(r - w1[i]);
        if (sr < s) s = sr;
    }

    // Calculate the distance to the nearest theta wall.
    
    if (nw2 != 2) {
        //double theta = P->theta;
        double rho = sqrt(P->r[0]*P->r[0] + P->r[1]*P->r[1]);
        
        for (int i=P->l[1]; i <= P->l[1]+1; i++) {
            double st = fabs(P->r[3]*sin(w2[i]) - rho*cos(w3[i]));
            if (st < s) s = st;
        }
    }

    // Calculate the distance to the nearest phi wall.
    
    if (nw3 != 2) {
        //double phi = P->phi;
        
        for (int i=P->l[2]; i <= P->l[2]+1; i++) {
            Vector<double, 3> phi_hat = Vector<double, 3>(-sin(w3[i]), 
                    cos(w3[i]), 0);

            double sp = fabs(phi_hat * P->r);
            if (sp < s) s = sp;
        }
    }
    
    return s;
}

/* Calculate the smallest distance across the cell. */

double SphericalGrid::smallest_wall_size(Photon *P) {

    double s = fabs(w1[P->l[0]+1] - w1[P->l[0]]);

    if (nw2 != 2) {
        double r = w1[P->l[0]];
        if (w1[P->l[0]] == 0)
            r = w1[P->l[0]+1]*0.5;

        double st = fabs(r*(w2[P->l[1]+1] - w2[P->l[1]]));
        if (st < s) s = st;
    }
    
    if (nw3 != 2) {
        double r = w1[P->l[0]];
        if (w1[P->l[0]] == 0)
            r = w1[P->l[0]+1]*0.5;

        double sint = fmin(sin(w2[P->l[1]]), sin(w2[P->l[1]+1]));
        if (equal_zero(sint, EPSILON))
            sint = sin(0.5*(w2[P->l[1]] + w2[P->l[1]+1]));

        double sp = fabs(r * sint * (w3[P->l[2]+1] - w3[P->l[2]]));
        if (sp < s) s = sp;
    }
    
    return s;
}

double SphericalGrid::smallest_wall_size(Ray *R) {

    double s = fabs(w1[R->l[0]+1] - w1[R->l[0]]);

    if (nw2 != 2) {
        double r = w1[R->l[0]];
        if (w1[R->l[0]] == 0)
            r = w1[R->l[0]+1]*0.5;
        double st = fabs(r*fmax(0.05, w2[R->l[1]+1] - w2[R->l[1]]));
        if (st < s) s = st;
    }

    if (nw3 != 2) {
        double r = w1[R->l[0]];
        if (w1[R->l[0]] == 0)
            r = w1[R->l[0]+1]*0.5;
        double sint = fmin(sin(w2[R->l[1]]), sin(w2[R->l[1]+1]));
        if (equal_zero(sint, EPSILON))
            sint = sin(0.5*(w2[R->l[1]] + w2[R->l[1]+1]));
        double sp = fabs(r * sint * fmax(0.05, w3[R->l[2]+1] - w3[R->l[2]]));
        if (sp < s) s = sp;
    }

    return s;
}

/* Calculate the size of the grid. */

double SphericalGrid::grid_size() {
    return 2*w1[nw1-1];
}

/* Determine which cell the photon is in. */

Vector<int, 3> SphericalGrid::photon_loc(Photon *P) {
    Vector<int, 3> l;

    double pi = 3.14159265;
    P->rad = P->r.norm();
    /* If P->rad = 0 we need to be careful about how we calculate theta and 
     * phi. */
    if (P->rad == 0) {
        P->theta = pi - P->theta;
        P->l[1] = -1;
        P->phi = fmod(P->phi + pi, 2*pi);
        P->l[2] = -1;
    }
    else {
        P->theta = acos(P->r[2]/P->rad);
        P->phi = fmod(atan2(P->r[1],P->r[0])+2*pi,2*pi);
    }
    double r = P->rad;
    double theta = P->theta;
    double phi = P->phi;

    // Check if we are using mirror symmetry and we're in the southern
    // hemisphere. If so, we need to flip.
    
    if (mirror_symmetry) {
        if (theta > pi/2) {
            theta = pi - theta;
            P->n[2] *= -1;
        }

        if (equal_zero(cos(theta), EPSILON) and P->n[2] < 0)
            P->n[2] *= -1;
    }

    // Find the location in the radial grid.
    
    double gnx = sin(theta)*cos(phi);
    double gny = sin(theta)*sin(phi);
    double gnz = cos(theta);
    if (equal_zero(gnx, EPSILON)) gnx = 0.;
    if (equal_zero(gny, EPSILON)) gny = 0.;
    if (equal_zero(gnz, EPSILON)) gnz = 0.;
    
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
    if ((r == w1[l[0]]) && (P->n[0]*gnx+P->n[1]*gny+P->n[2]*gnz < 0))
        l[0] -= 1;
    else if ((r == w1[l[0]+1]) && (P->n[0]*gnx+P->n[1]*gny+P->n[2]*gnz >= 0))
        l[0] += 1;

    // Find the location in the theta grid.
    
    if (nw2 == 2)
        l[1] = 0;
    else {
        if (theta >= w2[nw2-1])
            l[1] = n2-1;
        else if (theta <= w2[0])
            l[1] = 0;
        else {
            if (P->l[1] == -1)
                l[1] = find_in_arr(theta,w2,nw2);
            else {
                int lower = P->l[1]-1;
                if (lower < 0) lower = 0;
                int upper = P->l[1]+1;
                if (upper > n2-1) upper = n2-1;
                
                l[1] = find_in_arr(theta,w2,lower,upper);
            }
            if (l[1] == n2) l[1] = n2-1;
        }

        /* Because of floating point errors it may be the case that the photon 
         * should be on the wall exactly, but is not exactly on the wall. We
         * need to put the photon exactly on the wall. */

        if (equal(theta,w2[l[1]],EPSILON))
            theta = w2[l[1]];
        else if (equal(theta,w2[l[1]+1],EPSILON))
            theta = w2[l[1]+1];

        /* Update which cell the photon is in based on the direction it
         * is going. */

        double gnx = cos(theta)*cos(phi);
        double gny = cos(theta)*sin(phi);
        double gnz = -sin(theta);
        if (equal_zero(gnx, EPSILON)) gnx = 0.;
        if (equal_zero(gny, EPSILON)) gny = 0.;
        if (equal_zero(gnz, EPSILON)) gnz = 0.;
        
        if ((theta == w2[l[1]]) && (P->n[0]*gnx+P->n[1]*gny+P->n[2]*gnz < 0))
            l[1] -= 1;
        else if ((theta == w2[l[1]+1]) && (P->n[0]*gnx+P->n[1]*gny+P->n[2]*gnz >= 0))
            l[1] += 1;

        /* Finally, if you somehow end up with l[1] = -1 or l[1] = n2, change
         * those to the correct values because you can't escape the grid in the
         * theta direction. */

        if (l[1] == -1) l[1] = 0;
        if (l[1] == n2) l[1] = n2-1;
    }

    // Find the location in the phi grid.
    
    if (nw3 == 2)
        l[2] = 0;
    else {
        if (P->l[2] == -1)
            l[2] = find_in_arr(phi,w3,nw3);
        else
            l[2] = find_in_periodic_arr(phi,w3,n3,P->l[2]-1,P->l[2]+1);

        if (l[2] == -1) find_in_arr(phi,w3,nw3);

        /* Check whether the photon is supposed to be exactly on the cell
         * wall. Floating point errors may keep it from being exactly on the
         * wall, and we need to fix that. */

        if (equal(phi,w3[l[2]],EPSILON))
            phi = w3[l[2]];
        else if (equal(phi,w3[l[2]+1],EPSILON))
            phi = w3[l[2]+1];

        /* Update which cell the photon is in depending on the 
         * direction it is going. */

        double gnx = -sin(phi);
        double gny = cos(phi);
        double gnz = 0.0;
        if (equal_zero(gnx, EPSILON)) gnx = 0.;
        if (equal_zero(gny, EPSILON)) gny = 0.;
        if (equal_zero(gnz, EPSILON)) gnz = 0.;
        
        if ((phi == w3[l[2]]) && (P->n[0]*gnx+P->n[1]*gny <= 0))
            l[2] -= 1;
        else if ((phi == w3[l[2]+1]) && (P->n[0]*gnx+P->n[1]*gny >= 0))
            l[2] += 1;
        l[2] = (l[2]+n3)%(n3);

        /* Finally, if you are at phi = 0, but going towards negative phi, 
         * you should set phi = 2*pi. */

        if ((phi == 0) && (l[2] == n3-1))
            phi = w3[l[2]+1];
    }

    /* Since we may have updated r, theta and phi to be exactly on the grid 
     * cell walls, change the photon position slightly to reflect this. */

    P->r[0] = r * sin(theta) * cos(phi);
    P->r[1] = r * sin(theta) * sin(phi);
    P->r[2] = r * cos(theta);
    P->rad = r;
    P->theta = theta;
    P->phi = phi;
    
    return l;
}

/* Randomly generate a photon location within a cell. */
 
Vector<double, 3> SphericalGrid::random_location_in_cell(int ix, int iy, 
        int iz) {
    double r = w1[ix] + random_number() * (w1[ix+1] - w1[ix]);
    double theta = w2[iy] + random_number() * (w2[iy+1] - w2[iy]);
    double phi = w3[iz] + random_number() * (w3[iz+1] - w3[iz]);

    double x = r * sin(theta) * cos(phi);
    double y = r * sin(theta) * sin(phi);
    double z = r * cos(theta);

    return Vector<double, 3>(x, y, z);
}

/* Check whether a photon is in the boundaries of the grid. */

bool SphericalGrid::in_grid(Photon *P) {

    /*double r = P->r.norm();

    if ((r >= w1[nw1-1]) || (equal(r,w1[nw1-1],EPSILON)))
        return false;
    else if ((r <= w1[0]) || (equal(r,w1[0],EPSILON)))
        return false; */
    if ((P->l[0] >= n1) || (P->l[0] < 0))
        return false;
    else
        return true;
}
