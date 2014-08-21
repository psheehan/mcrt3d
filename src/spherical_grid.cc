#ifndef SPHERICAL_GRID_CC
#define SPHERICAL_GRID_CC

#include "grid.cc"
#include "vector.cc"
#include "photon.cc"

struct SphericalGrid : public Grid {
    double next_wall_distance(Photon *P);
    double outer_wall_distance(Photon *P);
    double minimum_wall_distance(Photon *P);
    Vector<int, 3> photon_loc(Photon *P);
    bool in_grid(Photon *P);
};

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
            if ((sr1 < s) && (sr1 > 0)) s = sr1;
            double sr2 = -b - fabs(b);
            if ((sr2 < s) && (sr2 > 0)) s = sr2;
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
        if (r == w1[i]) {
            double sr = fabs(r - w1[i]);
            if (sr < s) s = sr;
        }
    }

    // Calculate the distance to the nearest theta wall.
    
    if (nw2 != 2) {
        double theta = P->theta;
        
        for (int i=P->l[1]; i <= P->l[1]+1; i++) {
        }
    }

    // Calculate the distance to the nearest phi wall.
    
    if (nw3 != 2) {
        double phi = P->phi;
        
        for (int i=P->l[2]; i <= P->l[2]+1; i++) {
            if (phi != w3[i]) {
                Vector<double, 3> phi_hat = Vector<double, 3>(-sin(w3[i]), 
                        cos(w3[i]), 0);

                double sp = fabs(phi_hat * P->r);
                if (sp < s) s = sp;
            }
        }
    }
    
    return s;
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

    // Find the location in the radial grid.
    
    double gnx = sin(theta)*cos(phi);
    double gny = sin(theta)*sin(phi);
    double gnz = cos(theta);
    
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

    if (equal(r,w1[l[0]],1.0e-6))
        r = w1[l[0]];
    else if (equal(r,w1[l[0]+1],1.0e-6))
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

        if (equal(theta,w2[l[1]],1.0e-6))
            theta = w2[l[1]];
        else if (equal(theta,w2[l[1]+1],1.0e-6))
            theta = w2[l[1]+1];

        /* Update which cell the photon is in based on the direction it
         * is going. */

        double gnx = cos(theta)*cos(phi);
        double gny = cos(theta)*sin(phi);
        double gnz = -sin(theta);
        
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

        /* Check whether the photon is supposed to be exactly on the cell
         * wall. Floating point errors may keep it from being exactly on the
         * wall, and we need to fix that. */

        if (equal(phi,w3[l[2]],1.0e-6))
            phi = w3[l[2]];
        else if (equal(phi,w3[l[2]+1],1.0e-6))
            phi = w3[l[2]+1];

        /* Update which cell the photon is in depending on the 
         * direction it is going. */

        double gnx = -sin(phi);
        double gny = cos(phi);
        double gnz = 0.0;
        
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

/* Check whether a photon is in the boundaries of the grid. */

bool SphericalGrid::in_grid(Photon *P) {

    /*double r = P->r.norm();

    if ((r >= w1[nw1-1]) || (equal(r,w1[nw1-1],1.0e-6)))
        return false;
    else if ((r <= w1[0]) || (equal(r,w1[0],1.0e-6)))
        return false; */
    if ((P->l[0] >= n1) || (P->l[0] < 0))
        return false;
    else
        return true;
}

#endif
