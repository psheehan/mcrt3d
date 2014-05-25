#ifndef CARTESIAN_GRID_CC
#define CARTESIAN_GRID_CC

#include "grid.cc"
#include "vector.cc"
#include "photon.cc"

struct CartesianGrid : public Grid {
    double next_wall_distance(Photon *P);
    Vector<int, 3> photon_loc(Photon *P, bool verbose);
    bool in_grid(Photon *P);
    bool in_grid_raytracing(Ray *R);
};

/* Calculate the distance between the photon and the nearest wall. */

double CartesianGrid::next_wall_distance(Photon *P) {
    double s = HUGE_VAL;

    for (int i=P->l[0]; i <= P->l[0]+1; i++) {
        double sx = (w1[i] - P->r[0])*P->invn[0];
        if ((sx < s) && (sx > 0)) s = sx;
    }
    
    for (int i=P->l[1]; i <= P->l[1]+1; i++) {
        double sy = (w2[i] - P->r[1])*P->invn[1];
        if ((sy < s) && (sy > 0)) s = sy;
    }
    
    for (int i=P->l[2]; i <= P->l[2]+1; i++) {
        double sz = (w3[i] - P->r[2])*P->invn[2];
        if ((sz < s) && (sz > 0)) s = sz;
    }
    
    return s;
}

/* Determine which cell the photon is in. */

Vector<int, 3> CartesianGrid::photon_loc(Photon *P, bool verbose) {
    Vector<int, 3> l;
    
    if (P->r[0] >= w1[nw1-1])
        l[0] = n1-1;
    else if (P->r[0] <= w1[0])
        l[0] = 0;
    else {
        if (P->l[0] == -1)
            l[0] = find_in_arr(P->r[0],w1,nw1);
        else {
            int lower = P->l[0]-1;
            if (lower < 0) lower = 0;
            int upper = P->l[0]+1;
            if (upper > n1-1) upper = n1-1;
            
            l[0] = find_in_arr(P->r[0],w1,lower,upper);
        }
    }
    
    if ((P->r[0] == w1[l[0]]) && (P->n[0] <= 0))
        l[0] -= 1;
    else if ((P->r[0] == w1[l[0]+1]) && (P->n[0] >= 0))
        l[0] += 1;
    
    if (P->r[1] >= w2[nw2-1])
        l[1] = n2-1;
    else if (P->r[1] <= w2[0])
        l[1] = 0;
    else {
        if (P->l[1] == -1)
            l[1] = find_in_arr(P->r[1],w2,nw2);
        else {
            int lower = P->l[1]-1;
            if (lower < 0) lower = 0;
            int upper = P->l[1]+1;
            if (upper > n2-1) upper = n2-1;
            
            l[1] = find_in_arr(P->r[1],w2,lower,upper);
        }
    }
    
    if ((P->r[1] == w2[l[1]]) && (P->n[1] <= 0))
        l[1] -= 1;
    else if ((P->r[1] == w2[l[1]+1]) && (P->n[1] >= 0))
        l[1] += 1;
    
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
    
    if ((P->r[2] == w3[l[2]]) && (P->n[2] <= 0))
        l[2] -= 1;
    else if ((P->r[2] == w3[l[2]+1]) && (P->n[2] >= 0))
        l[2] += 1;
    
    return l;
}

/* Check whether a photon is in the boundaries of the grid. */

bool CartesianGrid::in_grid(Photon *P) {
    /*if ((P->r[0] >= w1[nw1-1]) || (P->r[1] >= w2[nw2-1]) ||
        (P->r[2] >= w3[nw2-1]) || (P->r[0] <= w1[0]) ||
        (P->r[1] <= w2[0]) || (P->r[2] <= w3[0]))
        return false;*/
    if ((P->l[0] >= n1) || (P->l[0] < 0) || (P->l[1] >= n2) || (P->l[1] < 0) ||
            (P->l[2] >= n3) || (P->l[2] < 0))
        return false;
    else
        return true;
}

bool CartesianGrid::in_grid_raytracing(Ray *R) {
    if (((R->r[0] >= w1[0]) || (equal(R->r[0],w1[0],1.0e-6))) &&
        ((R->r[1] >= w2[0]) || (equal(R->r[1],w2[0],1.0e-6))) &&
        ((R->r[2] >= w3[0]) || (equal(R->r[2],w3[0],1.0e-6))) &&
        ((R->r[0] <= w1[nw1-1]) || (equal(R->r[0],w1[nw1-1],1.0e-6))) &&
        ((R->r[1] <= w2[nw2-1]) || (equal(R->r[1],w2[nw2-1],1.0e-6))) &&
        ((R->r[2] <= w3[nw3-1]) || (equal(R->r[2],w3[nw3-1],1.0e-6))))
        return true;
    else
        return false;
}

#endif
