#ifndef CYLINDRICAL_GRID_CC
#define CYLINDRICAL_GRID_CC

#include "grid.cc"
#include "vector.cc"
#include "photon.cc"

struct CylindricalGrid : public Grid {
    double next_wall_distance(Photon *P);
    Vector<int, 3> photon_loc(Photon *P, bool verbose);
    bool in_grid(Photon *P);
};

/* Calculate the distance between the photon and the nearest wall. */

double CylindricalGrid::next_wall_distance(Photon *P) {

    double r = sqrt(P->r[0]*P->r[0]+P->r[1]*P->r[1]);

    // Calculate the distance to the intersection with the next radial wall.
    
    double a = P->n[0]*P->n[0]+P->n[1]*P->n[1];
    double b = P->r[0]*P->n[0]+P->r[1]*P->n[1];

    int l = P->l[0]-1;
    if (l < 0) l = 0;
    int u = P->l[0]+2;
    if (u >= nw1) u = nw1-1;
    
    double s = HUGE_VAL;
    for (int i=l; i <= u; i++) {
        double c = r*r - w1[i]*w1[i];
        double d = b*b - a*c;

        if (d >= 0) {
            double sr1 = (-b + sqrt(d))/a;
            double sr2 = (-b - sqrt(d))/a;

            if ((sr1 < s) && (sr1 > 1.0e-6*dw1[i])) s = sr1;
            if ((sr2 < s) && (sr2 > 1.0e-6*dw1[i])) s = sr2;
        }
    }

    // Calculate the distance to the intersection with the next phi wall.
    
    if (nw2 != 2) {
        l = P->l[1]-1;
        u = P->l[1]+2;
        
        for (int i=l; i <= u; i++) {
            int index = (i+(nw2-1))%(nw2-1);
            double c = P->r[0]*sin(w2[index])-P->r[1]*cos(w2[index]);
            double d = P->n[0]*sin(w2[index])-P->n[1]*cos(w2[index]);

            double sp = -c/d;

            if ((sp < s) && (sp > 1.0e-6*r*dw2[index])) s = sp;
        }
    }

    // Calculate the distance to intersection with the nearest z wall.
    
    l = P->l[2]-1;
    if (l < 0) l = 0;
    u = P->l[2]+2;
    if (u >= nw3) u = nw3-1;
    
    for (int i=l; i <= u; i++) {
        double sz = (w3[i] - P->r[2])*P->invn[2];
        if ((sz < s) && (sz > 1.0e-6*dw3[i])) s = sz;
    }
    
    return s;
}

/* Determine which cell the photon is in. */

Vector<int, 3> CylindricalGrid::photon_loc(Photon *P, bool verbose) {
    Vector<int, 3> l;

    // Find the location in the radial grid.
    
    double pi = 3.14159265;
    double r = sqrt(P->r[0]*P->r[0]+P->r[1]*P->r[1]);
    double phi = fmod(atan2(P->r[1],P->r[0])+2*pi,2*pi);

    double gnx = cos(phi);
    double gny = sin(phi);
    
    if (r >= w1[nw1-1])
        l[0] = nw1-1;
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
    
    if ((equal(r,w1[l[0]],1.0e-6)) && (P->n[0]*gnx+P->n[1]*gny <= 0))
        l[0] -= 1;
    else if ((equal(r,w1[l[0]+1],1.0e-6)) && (P->n[0]*gnx+P->n[1]*gny >= 0))
        l[0] += 1;

    // Find the location in the phi grid.
    
    if (nw2 == 2)
        l[1] = 0;
    else {
        if (phi >= w2[nw2-1])
            l[1] = 0;
        else if (phi <= w2[0])
            l[1] = nw2-2;
        else {
            if (P->l[1] == -1)
                l[1] = find_in_arr(phi,w2,nw2);
            else {
                int lower = P->l[1]-1;
                int upper = P->l[1]+2;
                
                l[1] = find_in_periodic_arr(phi,w2,nw3,lower,upper);
            }
            if (l[1] == nw2-1) l[1] = 0;
        }

        double gnx = -sin(phi);
        double gny = cos(phi);
        
        if ((equal(phi,w2[l[1]],1.0e-3)) && (P->n[0]*gnx+P->n[1]*gny <= 0))
            l[1] -= 1;
        else if ((equal(phi,w2[l[1]+1],1.0e-3)) && 
                (P->n[0]*gnx+P->n[1]*gny >= 0))
            l[1] += 1;
        l[1] = (l[1]+nw2-1)%(nw2-1);
    }

    // Find the location in the z grid.
    
    if (P->r[2] >= w3[nw3-1])
        l[2] = nw3-1;
    else if (P->r[2] <= w3[0])
        l[2] = 0;
    else {
        if (P->l[2] == -1)
            l[2] = find_in_arr(P->r[2],w3,nw3);
        else {
            int lower = P->l[2]-1;
            if (lower < 0) lower = 0;
            int upper = P->l[2]+2;
            if (upper > nw3-1) upper = nw3-1;
            
            l[2] = find_in_arr(P->r[2],w3,lower,upper);
        }
    }
    
    if ((equal(P->r[2],w3[l[2]],1.0e-6)) && (P->n[2] <= 0))
        l[2] -= 1;
    else if ((equal(P->r[2],w3[l[2]+1],1.0e-6)) && (P->n[2] >= 0))
        l[2] += 1;
    
    return l;
}

/* Check whether a photon is in the boundaries of the grid. */

bool CylindricalGrid::in_grid(Photon *P) {

    double r = sqrt(P->r[0]*P->r[0]+P->r[1]*P->r[1]);

    if ((r >= w1[nw1-1]) || (P->r[2] >= w3[nw3-1]) || 
            (equal(r,w1[nw1-1],1.0e-6)) || (equal(P->r[2],w3[nw3-1],1.0e-6)))
        return false;
    else if ((r <= w1[0]) || (equal(r,w1[0],1.0e-6)))
        return false;
    else if ((P->r[2] <= w3[0]) || (equal(P->r[2],w3[0],1.0e-6)))
        return false;
    else
        return true;
}

#endif
