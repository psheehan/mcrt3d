#ifndef SPHERICAL_GRID_CC
#define SPHERICAL_GRID_CC

#include "grid.cc"
#include "vector.cc"
#include "photon.cc"

struct SphericalGrid : public Grid {
    double next_wall_distance(Photon *P);
    Vector<int, 3> photon_loc(Photon *P, bool verbose);
    bool in_grid(Photon *P);
};

/* Calculate the distance between the photon and the nearest wall. */

double SphericalGrid::next_wall_distance(Photon *P) {

    double r = P->r.norm();

    // Calculate the distance to the intersection with the next radial wall.
    
    double b = P->r*P->n;

    int l = P->l[0]-1;
    if (l < 0) l = 0;
    int u = P->l[0]+2;
    if (u >= nw1) u = nw1-1;
    
    double s = HUGE_VAL;
    for (int i=l; i <= u; i++) {
        double c = r*r - w1[i]*w1[i];
        double d = b*b - c;

        if (d >= 0) {
            double sr1 = -b + sqrt(d);
            double sr2 = -b - sqrt(d);

            if ((sr1 < s) && (sr1 > 1.0e-6*dw1[i])) s = sr1;
            if ((sr2 < s) && (sr2 > 1.0e-6*dw1[i])) s = sr2;
        }
    }

    // Calculate the distance to the intersection with the next theta wall.
    
    if (nw2 != 2) {
        l = P->l[1]-1;
        if (l < 0) l = 0;
        u = P->l[1]+2;
        if (u >= nw2) u = nw2-1;
        
        for (int i=l; i <= u; i++) {
            double a = P->n[0]*P->n[0]+P->n[1]*P->n[1]-P->n[2]*P->n[2]*
                pow(tan(w2[i]),2);
            double b = 2*(P->r[0]*P->n[0]+P->r[1]*P->n[1]-P->r[2]*P->n[2]*
                pow(tan(w2[i]),2));
            double c = P->r[0]*P->r[0]+P->r[1]*P->r[1]-P->r[2]*P->r[2]*
                pow(tan(w2[i]),2);
            double d = b*b-4*a*c;

            if (d >= 0) {
                double st1 = (-b + sqrt(d))/(2*a);
                double st2 = (-b - sqrt(d))/(2*a);

                if ((st1 < s) && (st1 > 1.0e-6*r*dw2[i])) s = st1;
                if ((st2 < s) && (st2 > 1.0e-6*r*dw2[i])) s = st2;
            }
        }
    }

    // Calculate the distance to intersection with the nearest phi wall.
    
    if (nw3 != 2) {
        l = P->l[2]-1;
        u = P->l[2]+2;

        double twodr = sqrt(P->r[0]*P->r[0]+P->r[1]*P->r[1]);
        
        for (int i=l; i <= u; i++) {
            int index = (i+(nw3-1))%(nw3-1);
            double c = P->r[0]*sin(w3[index])-P->r[1]*cos(w3[index]);
            double d = P->n[0]*sin(w3[index])-P->n[1]*cos(w3[index]);

            double sp = -c/d;

            if ((sp < s) && (sp > 1.0e-6*twodr*dw3[index])) s = sp;
        }
    }
    
    return s;
}

/* Determine which cell the photon is in. */

Vector<int, 3> SphericalGrid::photon_loc(Photon *P, bool verbose) {
    Vector<int, 3> l;

    // Find the location in the radial grid.
    
    double pi = 3.1415926;
    double r = P->r.norm();
    double theta = acos(P->r[2]/r);
    double phi = fmod(atan2(P->r[1],P->r[0])+2*pi,2*pi);

    double gnx = sin(theta)*cos(phi);
    double gny = sin(theta)*sin(phi);
    double gnz = cos(theta);
    
    if (r >= w1[nw1-1])
        l[0] = nw1-2;
    else if (r <= w1[0])
        l[0] = 0;
    else {
        if (P->l[0] == -1)
            l[0] = find_in_arr(r,w1,nw1);
        else {
            int lower = P->l[0]-1;
            if (lower < 0) lower = 0;
            int upper = P->l[0]+2;
            if (upper > nw1-1) upper = nw1-1;
            
            l[0] = find_in_arr(r,w1,lower,upper);
        }
    }
    
    if ((equal(r,w1[l[0]],1.0e-6)) && 
            (P->n[0]*gnx+P->n[1]*gny+P->n[2]*gnz <= 0))
        l[0] -= 1;
    else if ((equal(r,w1[l[0]+1],1.0e-6)) && 
            (P->n[0]*gnx+P->n[1]*gny+P->n[2]*gnz >= 0))
        l[0] += 1;

    // Find the location in the theta grid.
    
    if (nw2 == 2)
        l[1] = 0;
    else {
        if (theta >= w2[nw2-1])
            l[1] = nw2-2;
        else if (theta <= w2[0])
            l[1] = 0;
        else {
            if (P->l[1] == -1)
                l[1] = find_in_arr(theta,w2,nw2);
            else {
                int lower = P->l[1]-1;
                if (lower < 0) lower = 0;
                int upper = P->l[1]+2;
                if (upper > nw2-1) upper = nw2-1;
                
                l[1] = find_in_arr(theta,w2,lower,upper);
            }
            if (l[1] == nw2-1) l[1]--;
        }

        double gnx = cos(theta)*cos(phi);
        double gny = cos(theta)*sin(phi);
        double gnz = -sin(theta);
        
        if ((equal(theta,w2[l[1]],1.0e-3)) && 
                (P->n[0]*gnx+P->n[1]*gny+P->n[2]*gnz < 0) && (l[1] != 0))
            l[1] -= 1;
        else if ((equal(theta,w2[l[1]+1],1.0e-3)) && 
                (P->n[0]*gnx+P->n[1]*gny+P->n[2]*gnz >= 0) && (l[1] != nw2-2))
            l[1] += 1;
    }

    // Find the location in the phi grid.
    
    if (nw3 == 2)
        l[2] = 0;
    else {
        if (phi >= w3[nw3-1])
            l[2] = 0;
        else if (phi <= w3[0])
            l[2] = nw3-2;
        else {
            if (P->l[2] == -1)
                l[2] = find_in_arr(phi,w3,nw3);
            else {
                int lower = P->l[2]-1;
                int upper = P->l[2]+2;
                
                l[2] = find_in_periodic_arr(phi,w3,nw3,lower,upper);
            }
            if (l[2] == nw3-1) l[2] = 0;
        }

        double gnx = -sin(phi);
        double gny = cos(phi);
        double gnz = 0.0;
        
        if ((equal(phi,w3[l[2]],1.0e-3)) && 
                (P->n[0]*gnx+P->n[1]*gny+P->n[2]*gnz <= 0) && (l[2] != 0))
            l[2] -= 1;
        else if ((equal(phi,w3[l[2]+1],1.0e-3)) && 
                (P->n[0]*gnx+P->n[1]*gny+P->n[2]*gnz >= 0) && (l[2] != nw3-2))
            l[2] += 1;
    }
    
    return l;
}

/* Check whether a photon is in the boundaries of the grid. */

bool SphericalGrid::in_grid(Photon *P) {

    double r = P->r.norm();

    if ((r >= w1[nw1-1]) || (equal(r,w1[nw1-1],1.0e-6)))
        return false;
    else if ((r <= w1[0]) || (equal(r,w1[0],1.0e-6)))
        return false;
    else
        return true;
}

#endif
