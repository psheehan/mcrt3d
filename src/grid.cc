#include "grid.h"

// Initialize the grid from numpy arrays directly.

Grid::Grid(py::array_t<double> __w1, py::array_t<double> __w2,
            py::array_t<double> __w3) {

    _w1 = __w1; _w2 = __w2; _w3 = __w3;

    // Load the array buffers to get the proper setup info.

    auto _w1_buf = _w1.request(); auto _w2_buf = _w2.request(); 
    auto _w3_buf = _w3.request();

    if (_w1_buf.ndim != 1 || _w2_buf.ndim != 1 || _w3_buf.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one");

    // Now get the correct format.
    
    nw1 = _w1_buf.shape[0], nw2 = _w2_buf.shape[0], nw3 = _w3_buf.shape[0];
    n1 = _w1_buf.shape[0]-1, n2 = _w2_buf.shape[0]-1, n3 = _w3_buf.shape[0]-1;
    w1 = (double *) _w1_buf.ptr; 
    w2 = (double *) _w2_buf.ptr; 
    w3 = (double *) _w3_buf.ptr;

    // Set up the volume of each cell.

    _volume = py::array_t<double>(n1*n2*n3);
    _volume.resize({n1, n2, n3});

    auto _volume_buf = _volume.request();
    double *__volume = (double *) _volume_buf.ptr;

    volume = pymangle(n1, n2, n3, __volume);

    // Make sure the number of species and sources are set correctly.

    nspecies = 0; nsources = 0;

    // Initialize the uses_mrw array.

    uses_mrw = create3DArrValue(n1, n2, n3, -1);
}

/* Functions to set up the grid. */

Grid::Grid(int _n1, int _n2, int _n3, int _nw1, int _nw2, int _nw3, 
        double *_w1, double *_w2, double *_w3, double *_volume, 
        bool _mirror_symmetry) {

    n1 = _n1;
    n2 = _n2;
    n3 = _n3;
    nw1 = _nw1;
    nw2 = _nw2;
    nw3 = _nw3;
    w1 = _w1;
    w2 = _w2;
    w3 = _w3;

    mirror_symmetry = _mirror_symmetry;

    volume = pymangle(n1, n2, n3, _volume);
}

Grid::~Grid() {
    // Free the physical parameters;

    freepymangle(volume);

    for (int i = 0; i < nspecies; i++) {
        freepymangle(dens[i]);
        freepymangle(temp[i]);
        delete3DArr(mass[i], n1, n2, n3);
        delete3DArr(energy[i], n1, n2, n3);
        delete3DArr(rosseland_mean_extinction[i], n1, n2, n3);
        delete3DArr(planck_mean_opacity[i], n1, n2, n3);
    }
    dens.clear(); temp.clear(); mass.clear(); energy.clear(); dust.clear();

    // Deallocate the uses_mrw array.

    delete3DArr(uses_mrw, n1, n2, n3);

    // Make sure the scattering array is deallocated.

    if (scatt.size() > 0)
        deallocate_scattering_array();

    // Clear the sources.
    
    for (int i = 0; i < nsources; i++) {
        delete sources[i];
    }
    sources.clear();
}

//void Grid::add_density(double *_dens, double *_temp, double *_mass, 
//        Dust *D) {
void Grid::add_density(py::array_t<double> ___dens, Dust *D) {

    // Deal with the Python code.

    auto _dens_buf = ___dens.request();

    _dens.append(___dens);

    // Create temperature array.

    py::array_t<double> ___temp = py::array_t<double>(n1*n2*n3);
    ___temp.resize({n1, n2, n3});

    _temp.append(___temp);

    auto _temp_buf = ___temp.request();

    // Now send to C++ useful things.

    double ***__dens = pymangle(n1, n2, n3, (double *) _dens_buf.ptr);
    double ***__temp = pymangle(n1, n2, n3, (double *) _temp_buf.ptr);
    double ***__mass = create3DArrValue(n1, n2, n3, 0);
    double ***__energy = create3DArrValue(n1, n2, n3, 0);
    double ***__rosseland_mean_extinction = create3DArrValue(n1, n2, n3, 0);
    double ***__planck_mean_opacity = create3DArrValue(n1, n2, n3, 0);

    dens.push_back(__dens);
    temp.push_back(__temp);
    mass.push_back(__mass);
    energy.push_back(__energy);
    rosseland_mean_extinction.push_back(__rosseland_mean_extinction);
    planck_mean_opacity.push_back(__planck_mean_opacity);

    // Initialize their values.

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            for (int k = 0; k < n3; k++) {
                __temp[i][j][k] = 0.1;
                __mass[i][j][k] = __dens[i][j][k] * volume[i][j][k];
                __rosseland_mean_extinction[i][j][k] = 
                    D->rosseland_mean_extinction(__temp[i][j][k]);
                __planck_mean_opacity[i][j][k] = D->planck_mean_opacity(
                        __temp[i][j][k]);
            }
        }
    }

    // Add the dust to the list of dust classes.

    dust.push_back(D);
    _dust.append(D);
    nspecies++;
}

//void Grid::add_source(Source *S) {
void Grid::add_star(double x, double y, double z, double _mass, double _radius, 
            double _temperature) {
    Star *S = new Star(x, y, z, _mass, _radius, _temperature);

    //TODO: calculate the lookup tables somehow!
    
    sources.push_back(S);
    _sources.append(S);
    nsources++;
}

void Grid::set_mrw_tables(double *_y, double *_f, double *_dydf, int _ny) {
    y = _y;
    f = _f;
    dydf = _dydf;
    ny = _ny;
}

/*void Grid::add_scattering_array(double *_scatt, int nnu) {
    double ****__scatt = pymangle(n1, n2, n3, nnu, _scatt);

    scatt.push_back(__scatt);
}*/

void Grid::initialize_scattering_array() {
    // Create a 4D scattering array for each of the dust components.

    for (int idust = 0; idust<nspecies; idust++) {
        // Create temperature array in Numpy.

        py::array_t<double> ___scatt = py::array_t<double>(n1*n2*n3*Q->nnu);
        ___scatt.resize({n1, n2, n3, Q->nnu});

        _scatt.append(___scatt);

        // Get the buffer and create an array useful for C++.

        auto _scatt_buf = ___scatt.request();

        scatt.push_back(pymangle(n1, n2, n3, Q->nnu, 
                    (double *) _scatt_buf.ptr));

        set4DArrValue(scatt[idust], 0., n1, n2, n3, Q->nnu);
    }
}

void Grid::deallocate_scattering_array() {
    for (int idust = 0; idust < (int) scatt.size(); idust++)
        freepymangle(scatt[idust]);
    scatt.clear();
}

void Grid::initialize_luminosity_array() {
    total_lum = 0.;

    for (int idust = 0; idust<nspecies; idust++) {
        luminosity.push_back(create3DArrValue(n1, n2, n3, 0.));

        for (int ix = 0; ix<n1; ix++) {
            for (int iy = 0; iy<n2; iy++) {
                for (int iz = 0; iz<n3; iz++) {
                    luminosity[idust][ix][iy][iz] = cell_lum(idust, ix, iy, iz);

                    total_lum += luminosity[idust][ix][iy][iz];
                }
            }
        }
    }
}

void Grid::initialize_luminosity_array(double nu) {
    total_lum = 0.;

    for (int idust = 0; idust<nspecies; idust++) {
        luminosity.push_back(create3DArrValue(n1, n2, n3, 0.));

        for (int ix = 0; ix<n1; ix++) {
            for (int iy = 0; iy<n2; iy++) {
                for (int iz = 0; iz<n3; iz++) {
                    luminosity[idust][ix][iy][iz] = cell_lum(idust, ix, iy, 
                            iz, nu);

                    total_lum += luminosity[idust][ix][iy][iz];
                }
            }
        }
    }

    for (int idust = 0; idust<nspecies; idust++) {
        for (int ix = 0; ix<n1; ix++) {
            for (int iy = 0; iy<n2; iy++) {
                for (int iz = 0; iz<n3; iz++) {
                }
            }
        }
    }
}

void Grid::deallocate_luminosity_array() {
    for (int idust = 0; idust<nspecies; idust++) {
        delete3DArr(luminosity[idust], n1, n2, n3);
    }
    luminosity.clear();
}

/* Emit a photon from the grid. */

Photon *Grid::emit(int iphot) {
    /* Cycle through the various stars, having them emit photons one after 
     * another. This way each source will get the same nuber of photons 
     * +/- 1. */
    int isource = 0;
    int photons_per_source = Q->nphot;
    if (Q->scattering) {
        isource = fmod(iphot, nsources+1);
        photons_per_source = int(Q->nphot/(nsources+1));
    }
    else if (nsources > 1) {
        isource = fmod(iphot, nsources);
        photons_per_source = int(Q->nphot/nsources);
    }

    Photon *P;
    if (Q->scattering) {
        if (isource == nsources)
            P = emit(Q->scattering_nu[Q->inu], Q->dnu, photons_per_source);
        else
            P = sources[isource]->emit(Q->scattering_nu[Q->inu], Q->dnu, 
                    photons_per_source);
    }
    else
        P = sources[isource]->emit(photons_per_source);

    /* Calculate kext and albedo at the photon's current frequency for all
     * dust species. */

    P->current_kext = new double[nspecies];
    P->current_albedo = new double[nspecies];
    for (int i=0; i<nspecies; i++) {
        P->current_kext[i] = dust[i]->opacity(P->nu);
        P->current_albedo[i] = dust[i]->albdo(P->nu);
    }

    /* Check the photon's location in the grid. */
    P->l = photon_loc(P);

    return P;
}

Photon *Grid::emit(double _nu, double _dnu, int nphot) {
    Photon *P = new Photon();

    int idust, ix, iy, iz;

    // Get the cell that randomly emits.

    double rand = random_number();
    double cum_lum = 0.;

    for (idust = 0; idust<nspecies; idust++) {
        for (ix = 0; ix<n1; ix++) {
            for (iy = 0; iy<n2; iy++) {
                for (iz = 0; iz<n3; iz++) {
                    cum_lum += luminosity[idust][ix][iy][iz] / total_lum;

                    if (cum_lum >= rand)
                        break;
                }
                if (cum_lum >= rand)
                    break;
            }
            if (cum_lum >= rand)
                break;
        }
        if (cum_lum >= rand)
            break;
    }

    // Now set up the photon.

    P->r = random_location_in_cell(ix, iy, iz);

    double theta = pi*random_number();
    double phi = 2*pi*random_number();

    P->n[0] = sin(theta)*cos(phi);
    P->n[1] = sin(theta)*sin(phi);
    P->n[2] = cos(theta);

    P->invn[0] = 1.0/P->n[0];
    P->invn[1] = 1.0/P->n[1];
    P->invn[2] = 1.0/P->n[2];
    P->l[0] = -1;
    P->l[1] = -1;
    P->l[2] = -1;

    P->nu = _nu;

    P->energy = total_lum / nphot;

    return P;
}

/* Linker function to the dust absorb function. */

void Grid::absorb(Photon *P, int idust) {
    dust[idust]->absorb(P, temp[idust][P->l[0]][P->l[1]][P->l[2]], Q->bw);

    // Update the photon's arrays of kext and albedo since P->nu has changed
    // upon absorption.
    for (int i=0; i<nspecies; i++) {
        P->current_kext[i] = dust[i]->opacity(P->nu);
        P->current_albedo[i] = dust[i]->albdo(P->nu);
    }

    // Check the photon's location again because there's a small chance that 
    // the photon was absorbed on a wall, and if it was we may need to update
    // which cell it is in if the direction has changed.
    P->l = photon_loc(P);
}

void Grid::absorb_mrw(Photon *P, int idust) {
    dust[idust]->absorb_mrw(P, temp[idust][P->l[0]][P->l[1]][P->l[2]], Q->bw);

    // Update the photon's arrays of kext and albedo since P->nu has changed
    // upon absorption.
    for (int i=0; i<nspecies; i++) {
        P->current_kext[i] = dust[i]->opacity(P->nu);
        P->current_albedo[i] = dust[i]->albdo(P->nu);
    }

    // Check the photon's location again because there's a small chance that 
    // the photon was absorbed on a wall, and if it was we may need to update
    // which cell it is in if the direction has changed.
    P->l = photon_loc(P);
}

/* Linker function to the dust scatter function. */

void Grid::scatter(Photon *P, int idust) {
    dust[idust]->scatter(P);

    // Check the photon's location again because there's a small chance that 
    // the photon was absorbed on a wall, and if it was we may need to update
    // which cell it is in if the direction has changed.
    P->l = photon_loc(P);
}

/* Propagate a photon through the grid until it escapes. */

void Grid::propagate_photon_full(Photon *P) {
    P->same_cell_count = 0;

    while (in_grid(P)) {
        // Determin the optical depth that the photon can travel until it's
        // next interaction.
        double tau = -log(1-random_number());

        // Figure out what that next action is, absorption or scattering. This
        // is figured out early for the sake of the continuous absorption
        // method.
        double albedo;
        int idust;
        if (nspecies == 1) {
            albedo = P->current_albedo[0];
            idust = 0;
        }
        else {
            double ksca_tot = 0;
            double *ksca_cum = new double[nspecies];
            double kext_tot = 0;
            for (int i=0; i<nspecies; i++) {
                ksca_tot += P->current_albedo[i]*P->current_kext[i]*
                    dens[i][P->l[0]][P->l[1]][P->l[2]];
                ksca_cum[i] = ksca_tot;
                kext_tot += P->current_kext[i]*
                    dens[i][P->l[0]][P->l[1]][P->l[2]];
            }
            albedo = ksca_tot / kext_tot;

            double rand = random_number();
            for (int i=0; i<nspecies; i++) {
                if (rand < ksca_cum[i] / kext_tot) {
                    idust = i;
                    break;
                }
            }
            delete[] ksca_cum;
        }

        bool absorb_photon = random_number() > albedo;

        // Move the photon to the point of it's next interaction.
        if (Q-> scattering)
            propagate_photon(P, tau, false);
        else
            propagate_photon(P, tau, absorb_photon);

        // If the photon is still in the grid when it reaches it's 
        // destination...
        if (in_grid(P)) {
            P->event_count += 1;
            // If the next interaction is absorption...
            if (absorb_photon) {
                absorb(P, idust);
                // If we've asked for verbose output, print some info.
                if (Q->verbose) {
                    printf("Absorbing photon at %i  %i  %i\n", P->l[0],
                            P->l[1], P->l[2]);
                    printf("Absorbed in a cell with temperature: %f\n",
                            temp[idust][P->l[0]][P->l[1]][P->l[2]]);
                    printf("Re-emitted with direction: %f  %f  %f\n",
                            P->n[0], P->n[1], P->n[2]);
                    printf("Re-emitted with frequency: %e\n", P->nu);
                }
            }
            // Otherwise, scatter the photon.
            else {
                // Now scatter the photon.
                scatter(P, idust);
                // If we've asked for verbose output, print some info.
                if (Q->verbose) {
                    printf("Scattering photon at cell  %i  %i  %i\n",
                            P->l[0], P->l[1], P->l[2]);
                    printf("Scattered with direction: %f  %f  %f\n",
                            P->n[0], P->n[1], P->n[2]);
                }
            }

            // If we've enabled MRW, try running an MRW step.
            // Need to do the in_grid check again because the absorption
            // or scattering may have been on the outermost wall, putting
            // it outside of the grid.
            if (Q->use_mrw && in_grid(P) && ((P->same_cell_count > 400) || 
                        (uses_mrw[P->l[0]][P->l[1]][P->l[2]] > 0))) {
                double alpha = 0.0;
                for (int idust = 0; idust<nspecies; idust++)
                    alpha += rosseland_mean_extinction[idust][P->l[0]][P->l[1]]
                            [P->l[2]] * dens[idust][P->l[0]][P->l[1]][P->l[2]];

                double dmin = smallest_wall_size(P);

                if (alpha * dmin > 20) {
                    uses_mrw[P->l[0]][P->l[1]][P->l[2]] = 1.0;

                    // Propagate the photon.
                    propagate_photon_mrw(P);

                    // Update the temperature after absorbing all that energy.
                    if (Q->bw) update_grid(P->l);

                    // Absorb the photon
                    int idust = 0;
                    if (nspecies == 1) {
                        idust = 0;
                    }
                    else {
                        double ksca_tot = 0;
                        double *ksca_cum = new double[nspecies];
                        double kext_tot = 0;
                        for (int i=0; i<nspecies; i++) {
                            ksca_tot += P->current_albedo[i]*P->current_kext[i]*
                                dens[i][P->l[0]][P->l[1]][P->l[2]];
                            ksca_cum[i] = ksca_tot;
                            kext_tot += P->current_kext[i]*
                                dens[i][P->l[0]][P->l[1]][P->l[2]];
                        }

                        double rand = random_number();
                        for (int i=0; i<nspecies; i++) {
                            if (rand < ksca_cum[i] / kext_tot) {
                                idust = i;
                                break;
                            }
                        }
                        delete[] ksca_cum;
                    }

                    absorb_mrw(P, idust);
                }
            }
        }
    }
}

/* Propagate a photon through the grid a distance equivalent to tau. */

void Grid::propagate_photon(Photon *P, double tau, bool absorb) {

    bool absorbed_by_source = false;
    int i = 0;
    while ((tau > EPSILON) && (in_grid(P))) {
        // Calculate the distance to the next wall.
        double s1 = next_wall_distance(P);

        // Calculate how far the photon can go with the current tau.
        double alpha = 0;
        for (int idust = 0; idust<nspecies; idust++)
            alpha += P->current_kext[idust]*
                dens[idust][P->l[0]][P->l[1]][P->l[2]];

        double s2 = tau/alpha;

        // Determine whether to move to the next wall or to the end of tau.
        double s = s1;
        if (s2 < s) {
            s = s2;
            P->same_cell_count++;
        } else
            P->same_cell_count = 0;

        // Calculate how far the photon can go before running into a source.
        int isource_intercept;
        for (int isource=0; isource<nsources; isource++) {
            double s3 = sources[isource]->intercept_distance(P);

            if (s3 < s) {
                s = s3;
                absorbed_by_source = true;
                isource_intercept = isource;
            }
        }

        // Continuously absorb the photon's energy, if the end result of the
        // current trajectory is absorption.
        if (absorb) {
            for (int idust=0; idust<nspecies; idust++)
                energy[idust][P->l[0]][P->l[1]][P->l[2]] += P->energy*
                    s*P->current_kext[idust]*
                    dens[idust][P->l[0]][P->l[1]][P->l[2]];
            // If we're doing a Bjorkman & Wood simulation, update the cell to
            // find its new temperature.
            if (Q->bw) {
                update_grid(P->l);
            }
        }

        // Remvove the tau we've used up with this stepl
        tau -= s*alpha;

        // Move the photon to it's new position.
        P->move(s);

        // If the photon moved to the next cell, update it's location.
        if (s1 < s2) P->l = photon_loc(P);
        i++;

        // If we've asked for verbose, print some information out.
        if (Q->verbose) {
            printf("%2i  %7.4f  %i  %7.4f  %7.4f  %7.4f\n", i, tau, P->l[0],
                    P->r[0]/au, s1*P->n[0]/au, s2*P->n[0]/au);
            printf("%14i  %7.4f  %7.4f  %7.4f\n", P->l[1], P->r[1]/au, 
                    s1*P->n[1]/au, s2*P->n[1]/au);
            printf("%14i  %7.4f  %7.4f  %7.4f\n", P->l[2], P->r[2]/au, 
                    s1*P->n[2]/au, s2*P->n[2]/au);
        }

        // If the distance to the star is the shortest distance, kill the 
        // photon.
        if (absorbed_by_source) {
            sources[isource_intercept]->reemit(P);
            P->l = photon_loc(P);
        }

        // Kill the photon if it bounces around too many times...
        if (i > 1000) {
            tau = -1.0;
            printf("!!!!!!! ERROR - Killing photon because it seems to be stuck.\n");
        }
    }
}

/* Propagate a photon through the grid for a scattering simulation. */

void Grid::propagate_photon_scattering(Photon *P) {
    double total_tau_abs = 0.;

    //while (in_grid(P) && P->energy > EPSILON) {
    while (in_grid(P) && total_tau_abs < 30.) {
        // Determin the optical depth that the photon can travel until it's
        // next interaction.
        double tau = -log(1-random_number());

        // Move the photon to the point of it's next interaction.
        while ((tau > EPSILON) && (in_grid(P))) {
            // Calculate the distance to the next wall.
            double s1 = next_wall_distance(P);

            // Calculate how far the photon can go with the current tau.
            double alpha_abs = 0;
            double alpha_scat = 0;
            for (int idust = 0; idust<nspecies; idust++) {
                alpha_abs += P->current_kext[idust]* 
                    (1.-P->current_albedo[idust])*
                    dens[idust][P->l[0]][P->l[1]][P->l[2]];
                alpha_scat += P->current_kext[idust]*P->current_albedo[idust]*
                    dens[idust][P->l[0]][P->l[1]][P->l[2]];
            }

            double s2 = tau/alpha_scat;

            // Determine whether to move to the next wall or to the end of tau.
            double s = s1;
            if (s2 < s) s = s2;

            // Add some of the energy to the scattering array.
            for (int idust=0; idust<nspecies; idust++)
                scatt[idust][P->l[0]][P->l[1]][P->l[2]][Q->inu] +=
                P->energy * s * P->current_albedo[idust] /
                (4*pi * mass[idust][P->l[0]][P->l[1]][P->l[2]]/
                dens[idust][P->l[0]][P->l[1]][P->l[2]]);

            // Absorb some of the photon's energy.
            P->energy *= exp(-s*alpha_abs);

            // Remvove the tau we've used up with this stepl
            tau -= s*alpha_scat;

            // Add to the total absorbed tau.
            total_tau_abs += s*alpha_abs;

            // Move the photon to it's new position.
            P->move(s);

            // If the photon moved to the next cell, update it's location.
            if (s1 < s2) P->l = photon_loc(P);
        }

        // If the photon is still in the grid when it reaches it's 
        // destination...
        if (in_grid(P)) {
            // Now scatter the photon.
            scatter(P, 0);
        }
    }
}

/* Get a random direction for the photon. */

void Grid::random_dir_mrw(Photon *P) {
    // Adjust the photon's direction if we are doing MRW.
    double cost = -1+2*random_number();
    double sint = sqrt(1-pow(cost,2));
    double phi = 2*pi*random_number();

    P->n[0] = sint*cos(phi);
    P->n[1] = sint*sin(phi);
    P->n[2] = cost;
    P->invn[0] = 1.0/P->n[0];
    P->invn[1] = 1.0/P->n[1];
    P->invn[2] = 1.0/P->n[2];
}

/* Propagate a photon using the MRW method for high optical depths. */

void Grid::propagate_photon_mrw(Photon *P) {
    // Calculate the Rosseland mean opacity.
    double alpha = 0.0;
    for (int idust = 0; idust<nspecies; idust++)
        alpha += rosseland_mean_extinction[idust][P->l[0]][P->l[1]][P->l[2]] * 
                dens[idust][P->l[0]][P->l[1]][P->l[2]];

    // Calculate the energy threshold at which to stop and re-calculate the
    // Rosseland mean opacity.
    double energy_threshold = 0.;
    for (int idust = 0; idust<nspecies; idust++) {
        if (temp[idust][P->l[0]][P->l[1]][P->l[2]] > 3.) {
            energy_threshold += energy[idust][P->l[0]][P->l[1]][P->l[2]];
        } else {
            energy_threshold = HUGE_VAL;
            break;
        }
    }
    energy_threshold *= (1. + 0.3);

    // Start walking the photon.
    while (true) {
        double dmin = minimum_wall_distance(P);

        if (alpha * dmin > Q->mrw_gamma) {
            if (Q->verbose) {
                printf("Doing MRW step...\n");
            }

            // Calculate the radius of the sphere that we will move the photon 
            // to somewhere on.
            double R_0 = 0.99 * dmin;
            if (Q->verbose) {
                printf("%f %f\n", dmin/au, R_0/au);
            }

            // Calculate the actual distance that the photon traveled through 
            // diffusion.
            double s = 0.5 * R_0 * R_0 * alpha; //Per RADMC-3D, this is the 
                                                //average distance.

            if (Q->verbose) {
                printf("%f\n", s/au);
            }
            // Add the energy absorbed into the grid.
            for (int idust=0; idust<nspecies; idust++)
                energy[idust][P->l[0]][P->l[1]][P->l[2]] += P->energy*
                        s*planck_mean_opacity[idust][P->l[0]][P->l[1]][P->l[2]]*
                        dens[idust][P->l[0]][P->l[1]][P->l[2]];

            // Move the photon to the edge of the sphere.
            P->move(R_0);

            // Pick a random direction to be traveling in next.
            random_dir_mrw(P);

        } else {
            double tau = -log(1-random_number());

            // Calculate the distance to the next wall.
            double s1 = next_wall_distance(P);

            // Calculate how far the photon can go with the current tau.
            double s2 = tau/alpha;

            // Determine whether to move to the next wall or to the end of tau.
            double s = s1;
            if (s2 < s) s = s2;

            // Continuously absorb the photon's energy, if the end result of the
            // current trajectory is absorption.
            for (int idust=0; idust<nspecies; idust++)
                energy[idust][P->l[0]][P->l[1]][P->l[2]] += P->energy*
                        s*planck_mean_opacity[idust][P->l[0]][P->l[1]][P->l[2]]*
                        dens[idust][P->l[0]][P->l[1]][P->l[2]];

            // Move the photon to it's new position.
            P->move(s);

            // If the photon moved to the next cell, update it's location.
            if (s1 < s2) {
                P->same_cell_count = 0;
                break;
            } else
                // "Absorb" to get a new direction.
                random_dir_mrw(P);
        }

        double energy_tot = 0.;
        for (int idust=0; idust<nspecies; idust++)
            energy_tot += energy[idust][P->l[0]][P->l[1]][P->l[2]];

        if (energy_tot > energy_threshold)
            break;
    }
}

/* Propagate a ray through the grid for raytracing. */

void Grid::propagate_ray(Ray *R) {

    int i=0;
    do {
        if (volume[R->l[0]][R->l[1]][R->l[2]] < 
                pi*R->pixel_size*R->pixel_size*R->pixel_size/6.) {
            R->pixel_too_large = true;
            //break;
        }

        double s = next_wall_distance(R);

        double tau_abs = 0;
        double tau_sca = 0;
        double tau_cell = 0;
        double intensity_abs = 0;
        double intensity_sca = 0;
        double intensity_cell = 0;
        for (int idust=0; idust<nspecies; idust++) {
            tau_abs = s*R->current_kext[idust]*(1.-R->current_albedo[idust])*
                dens[idust][R->l[0]][R->l[1]][R->l[2]];
            tau_sca = s*R->current_kext[idust]*R->current_albedo[idust]*
                dens[idust][R->l[0]][R->l[1]][R->l[2]];

            tau_cell += s*R->current_kext[idust]*
                dens[idust][R->l[0]][R->l[1]][R->l[2]];

            intensity_abs += (1.0-exp(-tau_cell)) * (1-R->current_albedo[idust])
                * planck_function(R->nu,temp[idust][R->l[0]][R->l[1]][R->l[2]]);
            intensity_sca += (1.0-exp(-tau_cell)) * R->current_albedo[idust] * 
                scatt[idust][R->l[0]][R->l[1]][R->l[2]][Q->inu];

            intensity_cell += intensity_abs + intensity_sca;
        }

        if (Q->verbose) {
            printf("%2i  %7.5f  %i  %7.4f  %7.4f\n", i, tau_cell, 
                    R->l[0], R->r[0]/au, s*R->n[0]/au);
            printf("%11.1e  %i  %7.4f  %7.4f\n", R->intensity, R->l[1], 
                    R->r[1]/au, s*R->n[1]/au);
            printf("%11.5f  %i  %7.4f  %7.4f\n", R->tau, R->l[2], R->r[2]/au, 
                    s*R->n[2]/au);
        }

        R->intensity += intensity_cell*exp(-R->tau);
        R->tau += tau_cell;

        R->move(s);

        R->l = photon_loc(R);

        i++;
    } while (in_grid(R));
}

/* Propagate a ray through the grid for raytracing. */

void Grid::propagate_ray_from_source(Ray *R) {

    int i=0;
    do {
        double s = next_wall_distance(R);

        double tau_cell = 0;
        for (int idust=0; idust<nspecies; idust++)
            tau_cell += s*R->current_kext[idust]*
                dens[idust][R->l[0]][R->l[1]][R->l[2]];

        if (Q->verbose) {
            printf("%2i  %7.5f  %i  %7.4f  %7.4f\n", i, tau_cell, 
                    R->l[0], R->r[0]/au, s*R->n[0]/au);
            printf("%11.1e  %i  %7.4f  %7.4f\n", R->intensity, R->l[1], 
                    R->r[1]/au, s*R->n[1]/au);
            printf("%11.5f  %i  %7.4f  %7.4f\n", R->tau, R->l[2], R->r[2]/au, 
                    s*R->n[2]/au);
        }

        R->intensity *= exp(-tau_cell);

        R->move(s);

        R->l = photon_loc(R);

        i++;
    } while (in_grid(R));
}

/* Calculate the distance between the photon and the nearest wall. */

double Grid::next_wall_distance(Photon *P) {
    return 0.0;
}

/* Calculate the distance between the photon and the outermost wall. */

double Grid::outer_wall_distance(Photon *P) {
    return 0.0;
}

/* Calculate the smallest absolute distance to the nearest wall. */

double Grid::minimum_wall_distance(Photon *P) {
    return 0.0;
}

/* Calculate the smallest distance across the cell. */

double Grid::smallest_wall_size(Photon *P) {
    return 0.0;
}

/* Calculate the size of the grid. */

double Grid::grid_size() {
    return 0.0;
}

/* Determine which cell the photon is in. */

Vector<int, 3> Grid::photon_loc(Photon *P) {
    return Vector<int, 3>();
}

/* Randomly generate a photon location within a cell. */
 
Vector<double, 3> Grid::random_location_in_cell(int ix, int iy, int iz) {
    return Vector<double, 3>();
}

/* Check whether a photon is on a wall and going parallel to it. */

bool Grid::on_and_parallel_to_wall(Photon *P) {
    return false;
}

/* Check whether a photon is in the boundaries of the grid. */

bool Grid::in_grid(Photon *P) {
    return true;
}

/* Update the temperature in a cell given the number of photons that have 
 * been absorbed in the cell. */

void Grid::update_grid(Vector<int, 3> l) {
    bool not_converged = true;

    for (int idust=0; idust<nspecies; idust++) {
        while (not_converged) {
            double T_old = temp[idust][l[0]][l[1]][l[2]];

            temp[idust][l[0]][l[1]][l[2]]=pow(energy[idust][l[0]][l[1]][l[2]]/
                (4*sigma*dust[idust]->\
                planck_mean_opacity(temp[idust][l[0]][l[1]][l[2]])*
                mass[idust][l[0]][l[1]][l[2]]),0.25);

            // Make sure that there is a minimum temperature that the grid can
            // get to.
            if (temp[idust][l[0]][l[1]][l[2]] < 0.1) 
                temp[idust][l[0]][l[1]][l[2]] = 0.1;

            if ((fabs(T_old-temp[idust][l[0]][l[1]][l[2]])/T_old < 1.0e-2))
                not_converged = false;
        }

        if (Q->use_mrw) {
            rosseland_mean_extinction[idust][l[0]][l[1]][l[2]] = 
                    dust[idust]->rosseland_mean_extinction(
                    temp[idust][l[0]][l[1]][l[2]]);
            planck_mean_opacity[idust][l[0]][l[1]][l[2]] = 
                    dust[idust]->planck_mean_opacity(
                    temp[idust][l[0]][l[1]][l[2]]);
        }
    }
}

void Grid::update_grid() {
    for (int i=0; i<nw1-1; i++)
        for (int j=0; j<nw2-1; j++)
            for (int k=0; k<nw3-1; k++)
                update_grid(Vector<int, 3>(i,j,k));
}

/* Calculate the luminosity of the cell indicated by l. */

double Grid::cell_lum(Vector<int, 3> l) {
    return 4*mass[0][l[0]][l[1]][l[2]]*dust[0]->
        planck_mean_opacity(temp[0][l[0]][l[1]][l[2]])*sigma*
        pow(temp[0][l[0]][l[1]][l[2]],4);
}

double Grid::cell_lum(int idust, int ix, int iy, int iz) {
    return 4*mass[idust][ix][iy][iz]*dust[idust]->
        planck_mean_opacity(temp[idust][ix][iy][iz])*sigma*
        pow(temp[idust][ix][iy][iz],4);
}

double Grid::cell_lum(Vector<int, 3> l, double nu) {
    return 4*pi*mass[0][l[0]][l[1]][l[2]]*
        dust[0]->opacity(nu)*(1. - dust[0]->albdo(nu))*
        planck_function(nu, temp[0][l[0]][l[1]][l[2]]);
}

double Grid::cell_lum(int idust, int ix, int iy, int iz, double nu) {
    return 4*pi*mass[idust][ix][iy][iz]*
        dust[0]->opacity(nu)*(1. - dust[0]->albdo(nu))*
        planck_function(nu, temp[idust][ix][iy][iz]);
}
