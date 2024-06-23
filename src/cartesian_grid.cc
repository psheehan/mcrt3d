#include "cartesian_grid.h"

CartesianGrid::CartesianGrid(py::array_t<double> _x, py::array_t<double> _y,
            py::array_t<double> _z) : Grid(_x, _y, _z) {

    // Set up the x, y, and z arrays.

    x = py::array_t<double>(n1);
    y = py::array_t<double>(n2);
    z = py::array_t<double>(n3);

    // Load the array buffers to get the proper setup info.

    auto x_buf = x.request(); auto y_buf = y.request(); 
    auto z_buf = z.request();

    // Now get the correct values.

    double *__x = (double *) x_buf.ptr;
    for (int i = 0; i < n1; i++) __x[i] = 0.5 * (w1[i+1] + w1[i]);

    double *__y = (double *) y_buf.ptr;
    for (int i = 0; i < n2; i++) __y[i] = 0.5 * (w2[i+1] + w2[i]);

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

double CartesianGrid::next_wall_distance(Photon *P) {
    double s = HUGE_VAL;

    double sx1 = (w1[P->l[0]] - P->r[0])*P->invn[0];
    if ((sx1 < s) && (sx1 > 0)) s = sx1;
    double sx2 = (w1[P->l[0]+1] - P->r[0])*P->invn[0];
    if ((sx2 < s) && (sx2 > 0)) s = sx2;
    
    double sy1 = (w2[P->l[1]] - P->r[1])*P->invn[1];
    if ((sy1 < s) && (sy1 > 0)) s = sy1;
    double sy2 = (w2[P->l[1]+1] - P->r[1])*P->invn[1];
    if ((sy2 < s) && (sy2 > 0)) s = sy2;
    
    double sz1 = (w3[P->l[2]] - P->r[2])*P->invn[2];
    if ((sz1 < s) && (sz1 > 0)) s = sz1;
    double sz2 = (w3[P->l[2]+1] - P->r[2])*P->invn[2];
    if ((sz2 < s) && (sz2 > 0)) s = sz2;
    
    return s;
}

/* Calculate the distance between the photon and the outermost wall. */

double CartesianGrid::outer_wall_distance(Photon *P) {
    double s = 0;

    if (P->n[0] != 0) {
        if (P->r[0] <= w1[0]) {
            double sx = (w1[0] - P->r[0])*P->invn[0];
            if (sx > s) s = sx;
        }
        else if (P->r[0] >= w1[nw1-1]) {
            double sx = (w1[nw1-1] - P->r[0])*P->invn[0];
            if (sx > s) s = sx;
        }
    }
    if (Q->verbose) printf("%7.4f\n", s/au);
    
    if (P->n[1] != 0) {
        if (P->r[1] <= w2[0]) {
            double sy = (w2[0] - P->r[1])*P->invn[1];
            if (sy > s) s = sy;
        }
        else if (P->r[1] >= w2[nw2-1]) {
            double sy = (w2[nw2-1] - P->r[1])*P->invn[1];
            if (sy > s) s = sy;
        }
    }
    if (Q->verbose) printf("%7.4f\n", s/au);
    
    if (P->n[2] != 0) {
        if (P->r[2] <= w3[0]) {
            double sz = (w3[0] - P->r[2])*P->invn[2];
            if (sz > s) s = sz;
        }
        else if (P->r[2] >= w3[nw3-1]) {
            double sz = (w3[nw3-1] - P->r[2])*P->invn[2];
            if (sz > s) s = sz;
        }
    }
    if (Q->verbose) printf("%7.4f\n", s/au);

    Vector<double, 3> newr = P->r + s*P->n;

    if (Q->verbose) printf("%20.17f   %7.4f   %7.4f\n", newr[0]/au, newr[1]/au, 
            newr[2]/au);

    if (equal(newr[0],w1[0],EPSILON)) newr[0] = w1[0];
    else if (equal(newr[0],w1[nw1-1],EPSILON)) newr[0] = w1[nw1-1];
    if (equal(newr[1],w2[0],EPSILON)) newr[1] = w2[0];
    else if (equal(newr[1],w2[nw2-1],EPSILON)) newr[1] = w2[nw2-1];
    if (equal(newr[2],w3[0],EPSILON)) newr[2] = w3[0];
    else if (equal(newr[2],w3[nw3-1],EPSILON)) newr[2] = w3[nw3-1];

    if (Q->verbose) printf("%20.17f   %7.4f   %7.4f\n", newr[0]/au, newr[1]/au, 
            newr[2]/au);
    if ((newr[0] < w1[0]) || (newr[0] > w1[nw1-1]) || (newr[1] < w2[0]) ||
            (newr[1] > w2[nw2-1]) || (newr[2] < w3[0]) || (newr[2] > w3[nw3-1]))
        s = HUGE_VAL;
    if (Q->verbose) printf("%7.4f\n", s/au);
    
    return s;
}

/* Calculate the smallest absolute distance to the nearest wall. */

double CartesianGrid::minimum_wall_distance(Photon *P) {
    double s = HUGE_VAL;

    double sx1 = fabs(w1[P->l[0]] - P->r[0]);
    if (sx1 < s) s = sx1;
    double sx2 = fabs(w1[P->l[0]+1] - P->r[0]);
    if (sx2 < s) s = sx2;
    
    double sy1 = fabs(w2[P->l[1]] - P->r[1]);
    if (sy1 < s) s = sy1;
    double sy2 = fabs(w2[P->l[1]+1] - P->r[1]);
    if (sy2 < s) s = sy2;
    
    double sz1 = fabs(w3[P->l[2]] - P->r[2]);
    if (sz1 < s) s = sz1;
    double sz2 = fabs(w3[P->l[2]+1] - P->r[2]);
    if (sz2 < s) s = sz2;
    
    return s;
}

/* Calculate the smallest distance across the cell. */

double CartesianGrid::smallest_wall_size(Photon *P) {

    double s = fabs(w1[P->l[0]+1] - w1[P->l[0]]);
    
    double sy = fabs(w2[P->l[0]+1] - w2[P->l[0]]);
    if (sy < s) s = sy;
    
    double sz = fabs(w3[P->l[0]+1] - w3[P->l[0]]);
    if (sz < s) s = sz;
    
    return s;
}

double CartesianGrid::smallest_wall_size(Ray *R) {

    // Use the cell volume as an estimator of the average size of a cell.

    double cell_volume = volume[R->l[0]*n2*n3 + R->l[1]*n3 + R->l[2]];

    double s = pow(cell_volume, 1./3);

    return s;
}

/* Calculate the size of the grid. */

double CartesianGrid::grid_size() {
    double rw1_max = std::max(abs(w1[0]), abs(w1[nw1-1]));
    double rw2_max = std::max(abs(w2[0]), abs(w2[nw2-1]));
    double rw3_max = std::max(abs(w3[0]), abs(w3[nw3-1]));

    return 2*sqrt(rw1_max*rw1_max + rw2_max*rw2_max + rw3_max*rw3_max);
}

/* Determine which cell the photon is in. */

Vector<int, 3> CartesianGrid::photon_loc(Photon *P) {
    Vector<int, 3> l;
    
    // Determine which cell the photon is currently in.

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

    /* Because of floating point errors it may be the case that the photon 
     * should be on the wall exactly, but is not exactly on the wall. We
     * need to put the photon exactly on the wall. */

    if (equal(P->r[0], w1[l[0]], EPSILON))
        P->r[0] = w1[l[0]];
    else if (equal(P->r[0], w1[l[0]+1], EPSILON))
        P->r[0] = w1[l[0]+1];

    /* Finally, update which cell the photon is in based on the direction it
     * is going. */
    
    if ((P->r[0] == w1[l[0]]) && (P->n[0] < 0))
        l[0] -= 1;
    else if ((P->r[0] == w1[l[0]+1]) && (P->n[0] > 0))
        l[0] += 1;
    
    // Determine which cell the photon is currently in.

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
    
    /* Because of floating point errors it may be the case that the photon 
     * should be on the wall exactly, but is not exactly on the wall. We
     * need to put the photon exactly on the wall. */

    if (equal(P->r[1], w2[l[1]], EPSILON))
        P->r[1] = w2[l[1]];
    else if (equal(P->r[1], w2[l[1]+1], EPSILON))
        P->r[1] = w2[l[1]+1];

    /* Finally, update which cell the photon is in based on the direction it
     * is going. */
    
    if ((P->r[1] == w2[l[1]]) && (P->n[1] < 0))
        l[1] -= 1;
    else if ((P->r[1] == w2[l[1]+1]) && (P->n[1] > 0))
        l[1] += 1;
    
    // Determine which cell the photon is currently in.

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
            if (upper > n3-1) upper = n3-1;
            
            l[2] = find_in_arr(P->r[2],w3,lower,upper);
        }
    }
    
    /* Because of floating point errors it may be the case that the photon 
     * should be on the wall exactly, but is not exactly on the wall. We
     * need to put the photon exactly on the wall. */

    if (equal(P->r[2], w3[l[2]], EPSILON))
        P->r[2] = w3[l[2]];
    else if (equal(P->r[2], w3[l[2]+1], EPSILON))
        P->r[2] = w3[l[2]+1];

    /* Finally, update which cell the photon is in based on the direction it
     * is going. */
    
    if ((P->r[2] == w3[l[2]]) && (P->n[2] < 0))
        l[2] -= 1;
    else if ((P->r[2] == w3[l[2]+1]) && (P->n[2] > 0))
        l[2] += 1;
    
    /* Also calculate n in the coordinate system frame. */

    P->nframe = P->n;

    /* And the cell index. */

    P->cell_index = l[0]*n2*n3 + l[1]*n3 + l[2];

    return l;
}

/* Update extra position parameters like rad and theta during MRW. */

void CartesianGrid::photon_loc_mrw(Photon *P) {
    /* Nothing needed here because we don't use rad, theta, or phi in
     * cartesian coordinates. */
}

/* Check whether a photon is on a wall and going parallel to it. */

bool CartesianGrid::on_and_parallel_to_wall(Photon *P) {
    if (P->r[0] == w1[P->l[0]] and P->n[0] == 0)
        return true;
    if (P->r[0] == w1[P->l[0]+1] and P->n[0] == 0)
        return true;

    if (P->r[1] == w2[P->l[1]] and P->n[1] == 0)
        return true;
    if (P->r[1] == w2[P->l[1]+1] and P->n[1] == 0)
        return true;

    if (P->r[2] == w3[P->l[2]] and P->n[2] == 0)
        return true;
    if (P->r[2] == w3[P->l[2]+1] and P->n[2] == 0)
        return true;

    return false;
}

/* Randomly generate a photon location within a cell. */
 
Vector<double, 3> CartesianGrid::random_location_in_cell(int ix, int iy, 
        int iz) {
    double x = w1[ix] + random_number(random_pool) * (w1[ix+1] - w1[ix]);
    double y = w2[iy] + random_number(random_pool) * (w2[iy+1] - w2[iy]);
    double z = w3[iz] + random_number(random_pool) * (w3[iz+1] - w3[iz]);

    return Vector<double, 3>(x, y, z);
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

