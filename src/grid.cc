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
        delete3DArr(rosseland_mean_extinction[i], n1, n2, n3);
        delete3DArr(planck_mean_opacity[i], n1, n2, n3);
        #ifdef _OPENMP
        delete3DLock(lock[i], n1, n2, n3);
        #endif
    }
    for (int i = 0; i < (int) energy.size(); i++) {
        for (int j = 0; j < nspecies; j++)
            delete3DArr(energy[i][j], n1, n2, n3);
        energy[i].clear();
    }
    dens.clear(); temp.clear(); mass.clear(); energy.clear(); dust.clear();
    #ifdef _OPENMP
    lock.clear();
    #endif

    // Deallocate the uses_mrw array.

    delete3DArr(uses_mrw, n1, n2, n3);

    // Make sure the scattering array is deallocated.

    if ((int) scatt.size() > 0)
        deallocate_scattering_array(0);
    
    // Clear the gas.
    
    for (int i = 0; i < ngases; i++) {
        freepymangle(number_dens[i]);
        freepymangle(gas_temp[i]);
        freepymangle(velocity[i]);
    }
    number_dens.clear(); gas_temp.clear(); velocity.clear(); gas.clear();

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
    #ifdef _OPENMP
    omp_lock_t*** __lock = create3DLock(n1, n2, n3);
    #endif

    dens.push_back(__dens);
    temp.push_back(__temp);
    mass.push_back(__mass);
    #ifdef _OPENMP
    lock.push_back(__lock);
    #endif
    if (energy.empty()) { 
        std::vector<double***> ___energy = {__energy};
        energy.push_back(___energy);
    } else
        energy[0].push_back(__energy);

    rosseland_mean_extinction.push_back(__rosseland_mean_extinction);
    planck_mean_opacity.push_back(__planck_mean_opacity);

    // Initialize their values.

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            for (int k = 0; k < n3; k++) {
                if (__dens[i][j][k] < 1.0e-99) __dens[i][j][k] = 1.0e-99;
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

void Grid::add_number_density(py::array_t<double> ___number_dens, 
        py::array_t<double> ___velocity, Gas *G) {

    // Deal with the Python code.

    auto _number_dens_buf = ___number_dens.request();
    auto _velocity_buf = ___velocity.request();

    _number_dens.append(___number_dens);
    _velocity.append(___velocity);

    // Create temperature array.

    py::array_t<double> ___gas_temp = py::array_t<double>(n1*n2*n3);
    ___gas_temp.resize({n1, n2, n3});

    _gas_temp.append(___gas_temp);

    auto _gas_temp_buf = ___gas_temp.request();

    // Now send to C++ useful things.

    double ***__number_dens = pymangle(n1, n2, n3, 
            (double *) _number_dens_buf.ptr);
    double ***__gas_temp = pymangle(n1, n2, n3, 
            (double *) _gas_temp_buf.ptr);
    double ****__velocity = pymangle(n1, n2, n3, 3, 
            (double *) _velocity_buf.ptr);

    number_dens.push_back(__number_dens);
    gas_temp.push_back(__gas_temp);
    velocity.push_back(__velocity);

    // Add the dust to the list of dust classes.

    gas.push_back(G);
    _gas.append(G);
    ngases++;
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

/* Functions to manage the scattering array. */

void Grid::add_scattering_array(py::array_t<double> ___scatt, int nthreads) {
    // Add the array to the Python-connected list of scattering phase
    // functions.

    _scatt.append(___scatt);

    // Get the buffer and create an array useful for C++.

    auto _scatt_buf = ___scatt.request();

    scatt.push_back(pymangle(n1, n2, n3, Q->nnu, (double *) _scatt_buf.ptr));

    // If we are using more than one thread, add additional arrays.

    for (int ithread = 1; ithread < nthreads; ithread++) {
        scatt.push_back(create4DArrValue(n1, n2, n3, Q->nnu, 0.));
    }
}

void Grid::initialize_scattering_array(int nthreads) {
    // Create a 4D scattering array for each of the separate threads.

    for (int ithread = 0; ithread < nthreads; ithread++) {
        scatt.push_back(create4DArrValue(n1, n2, n3, Q->nnu, 0.));
    }
}

void Grid::collapse_scattering_array() {
    for (int ithread = 1; ithread < (int) scatt.size(); ithread++)
        for (int ix=0; ix < n1; ix++)
            for (int iy=0; iy < n2; iy++)
                for (int iz=0; iz < n3; iz++)
                    for (int iq=0; iq < Q->nnu; iq++)
                        scatt[0][ix][iy][iz][iq] += 
                            scatt[ithread][ix][iy][iz][iq];

    deallocate_scattering_array(1);
}

void Grid::deallocate_scattering_array(int start) {
    for (int ithread = start; ithread < (int) scatt.size(); ithread++)
        delete4DArr(scatt[ithread], n1, n2, n3, Q->nnu);
    scatt.erase(scatt.begin()+start, scatt.end());
}

/* Functions to manage the energy array. */

void Grid::add_energy_arrays(int nthreads) {
    // Create a 3D energy array for each of the separate threads.

    for (int ithread = 1; ithread < nthreads; ithread++) {
        std::vector<double***> _energy = {};
        energy.push_back(_energy);
        for (int idust = 0; idust < nspecies; idust++)
            energy[ithread].push_back(create3DArrValue(n1, n2, n3, 0.));
    }
}

void Grid::deallocate_energy_arrays() {
    for (int ithread = 1; ithread < (int) energy.size(); ithread++) {
        for (int idust = 0; idust < nspecies; idust++)
            delete3DArr(energy[ithread][idust], n1, n2, n3);
        energy[ithread].clear();
    }
    energy.erase(energy.begin()+1, energy.end());
}

/* Functions to manage the luminosity array. */

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
    // Determine which dust species should do the absorption.
    idust = 0;
    if (nspecies == 1) {
        idust = 0;
    }
    else {
        double kabs_tot = 0;
        double *kabs_cum = new double[nspecies];
        for (int i=0; i<nspecies; i++) {
            kabs_tot += (1 - P->current_albedo[i])*P->current_kext[i]*
                dens[i][P->l[0]][P->l[1]][P->l[2]];
            kabs_cum[i] = kabs_tot;
        }

        double rand = random_number();
        for (int i=0; i<nspecies; i++) {
            if (rand < kabs_cum[i] / kabs_tot) {
                idust = i;
                break;
            }
        }
        delete[] kabs_cum;
    }

    // Now do the absorption.

    #ifdef _OPENMP
    omp_set_lock(&(lock[idust][P->l[0]][P->l[1]][P->l[2]]));
    #endif
    dust[idust]->absorb(P, temp[idust][P->l[0]][P->l[1]][P->l[2]], Q->bw);
    #ifdef _OPENMP
    omp_unset_lock(&(lock[idust][P->l[0]][P->l[1]][P->l[2]]));
    #endif

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
    // Determine which dust species should do the absorption.
    idust = 0;
    if (nspecies == 1) {
        idust = 0;
    }
    else {
        double kabs_tot = 0;
        double *kabs_cum = new double[nspecies];
        for (int i=0; i<nspecies; i++) {
            kabs_tot += (1 - P->current_albedo[i])*P->current_kext[i]*
                dens[i][P->l[0]][P->l[1]][P->l[2]];
            kabs_cum[i] = kabs_tot;
        }

        double rand = random_number();
        for (int i=0; i<nspecies; i++) {
            if (rand < kabs_cum[i] / kabs_tot) {
                idust = i;
                break;
            }
        }
        delete[] kabs_cum;
    }

    // Now do the absorption.

    #ifdef _OPENMP
    omp_set_lock(&(lock[idust][P->l[0]][P->l[1]][P->l[2]]));
    #endif
    dust[idust]->absorb_mrw(P, temp[idust][P->l[0]][P->l[1]][P->l[2]], Q->bw);
    #ifdef _OPENMP
    omp_unset_lock(&(lock[idust][P->l[0]][P->l[1]][P->l[2]]));
    #endif

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
    // Determine which dust species should do the absorption.
    idust = 0;
    if (nspecies == 1) {
        idust = 0;
    }
    else {
        double ksca_tot = 0;
        double *ksca_cum = new double[nspecies];
        for (int i=0; i<nspecies; i++) {
            ksca_tot += P->current_albedo[i]*P->current_kext[i]*
                dens[i][P->l[0]][P->l[1]][P->l[2]];
            ksca_cum[i] = ksca_tot;
        }

        double rand = random_number();
        for (int i=0; i<nspecies; i++) {
            if (rand < ksca_cum[i] / ksca_tot) {
                idust = i;
                break;
            }
        }
        delete[] ksca_cum;
    }

    // Now do the absorption.

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
        if (nspecies == 1) {
            albedo = P->current_albedo[0];
        }
        else {
            double ksca_tot = 0;
            double kext_tot = 0;
            for (int i=0; i<nspecies; i++) {
                ksca_tot += P->current_albedo[i]*P->current_kext[i]*
                    dens[i][P->l[0]][P->l[1]][P->l[2]];
                kext_tot += P->current_kext[i]*
                    dens[i][P->l[0]][P->l[1]][P->l[2]];
            }
            albedo = ksca_tot / kext_tot;
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
                absorb(P, 0);
                // If we've asked for verbose output, print some info.
                if (Q->verbose) {
                    printf("Absorbing photon at %i  %i  %i\n", P->l[0],
                            P->l[1], P->l[2]);
                    printf("Absorbed in a cell with temperature: %f\n",
                            temp[0][P->l[0]][P->l[1]][P->l[2]]);
                    printf("Re-emitted with direction: %f  %f  %f\n",
                            P->n[0], P->n[1], P->n[2]);
                    printf("Re-emitted with frequency: %e\n", P->nu);
                }
            }
            // Otherwise, scatter the photon.
            else {
                // Now scatter the photon.
                scatter(P, 0);
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
                    absorb_mrw(P, 0);
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
                energy[P->ithread][idust][P->l[0]][P->l[1]][P->l[2]] += 
                    P->energy*s*P->current_kext[idust]*
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
            absorbed_by_source = false;
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

        int i = 0;

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

            // Calculate the relevant optical depths.
            double tau_abs = s * alpha_abs;
            double tau_scat = s * alpha_scat;

            // Calculate the average energy over the course of the cell.
            double average_energy = (1.0 - exp(-tau_abs)) / tau_abs *
                P->energy;
            if (tau_abs < EPSILON)
                average_energy = (1.0 - 0.5*tau_abs) * P->energy;

            // Add some of the energy to the scattering array.
            scatt[P->ithread][P->l[0]][P->l[1]][P->l[2]][Q->inu] += 
                    average_energy * s / (4*pi * 
                    volume[P->l[0]][P->l[1]][P->l[2]]);

            // Absorb some of the photon's energy.
            P->energy *= exp(-tau_abs);

            // Remvove the tau we've used up with this stepl
            tau -= tau_scat;

            // Add to the total absorbed tau.
            total_tau_abs += tau_abs;

            // Move the photon to it's new position.
            P->move(s);

            // If the photon moved to the next cell, update it's location.
            if (s1 < s2) P->l = photon_loc(P);

            i++;
            // Kill the photon if it bounces around too many times...
            if (i > 1000) {
                tau = -1.0;
                printf("!!!!!!! ERROR - Killing photon because it seems to be stuck.\n");
            }
        }

        // If the photon is still in the grid when it reaches it's 
        // destination...
        if (in_grid(P)) {
            P->event_count += 1;
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
        #ifdef _OPENMP
        omp_set_lock(&(lock[idust][P->l[0]][P->l[1]][P->l[2]]));
        #endif
        if (temp[idust][P->l[0]][P->l[1]][P->l[2]] > 3.) {
            for (int ithread = 0; ithread < (int) energy.size(); ithread++)
                energy_threshold += energy[ithread][idust][P->l[0]]
                        [P->l[1]][P->l[2]];
        } else {
            energy_threshold = HUGE_VAL;
            break;
        }
        #ifdef _OPENMP
        omp_unset_lock(&(lock[idust][P->l[0]][P->l[1]][P->l[2]]));
        #endif
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
                energy[P->ithread][idust][P->l[0]][P->l[1]][P->l[2]] += 
                        P->energy*s*planck_mean_opacity[idust][P->l[0]][P->l[1]]
                        [P->l[2]]*dens[idust][P->l[0]][P->l[1]][P->l[2]];

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
                energy[P->ithread][idust][P->l[0]][P->l[1]][P->l[2]] += 
                        P->energy*s*planck_mean_opacity[idust][P->l[0]][P->l[1]]
                        [P->l[2]]*dens[idust][P->l[0]][P->l[1]][P->l[2]];

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
        for (int ithread = 0; ithread < (int) energy.size(); ithread++)
            for (int idust=0; idust<nspecies; idust++)
                energy_tot += energy[ithread][idust][P->l[0]][P->l[1]][P->l[2]];

        if (energy_tot > energy_threshold)
            break;
    }
}

/* Propagate a ray through the grid for raytracing. */

void Grid::propagate_ray(Ray *R) {

    int i=0;
    do {
        double s = next_wall_distance(R);

        if (smallest_wall_size(R) < R->pixel_size)
            R->pixel_too_large = true;

        for (int inu = 0; inu < R->nnu; inu++) {
            double tau_cell = 0;
            double intensity_abs = 0;
            double alpha_ext = 0;
            double alpha_sca = 0;

            for (int idust=0; idust<nspecies; idust++) {
                tau_cell += s*R->current_kext[idust][inu] *
                        dens[idust][R->l[0]][R->l[1]][R->l[2]];

                alpha_ext += R->current_kext[idust][inu] *
                        dens[idust][R->l[0]][R->l[1]][R->l[2]];
                alpha_sca += R->current_kext[idust][inu] *
                        R->current_albedo[idust][inu] *
                        dens[idust][R->l[0]][R->l[1]][R->l[2]];

                intensity_abs += R->current_kext[idust][inu] * 
                        (1 - R->current_albedo[idust][inu]) *
                        dens[idust][R->l[0]][R->l[1]][R->l[2]]
                        * planck_function(R->nu[inu], 
                        temp[idust][R->l[0]][R->l[1]][R->l[2]]);
            }

            double alpha_line = 0;
            double intensity_line = 0;
            for (int igas=0; igas < ngases; igas++) {
                double alpha_this_line = c_l*c_l / (8*pi*gas[igas]->nu[0]*
                        gas[igas]->nu[0]) * gas[igas]->A[0] * 
                        number_dens[igas][R->l[0]][R->l[1]][R->l[2]] * 
                        (exp(h_p * gas[igas]->nu[0] / (k_B * 
                        gas_temp[igas][R->l[0]][R->l[1]][R->l[2]])) - 1) * 
                        line_profile(igas, 0, R->l, R->nu[inu] * (1 - 
                        -R->n.dot(vector_velocity(igas, R)) / c_l));
                        
                alpha_line += alpha_this_line;
                tau_cell += s*alpha_this_line;

                intensity_line += alpha_this_line * planck_function(R->nu[inu],
                        gas_temp[igas][R->l[0]][R->l[1]][R->l[2]]);
            }

            double albedo = alpha_sca / alpha_ext;

            intensity_abs *= (1.0-exp(-tau_cell)) / alpha_ext;
            double intensity_sca = (1.0-exp(-tau_cell)) * albedo * 
                    scatt[0][R->l[0]][R->l[1]][R->l[2]][inu];

            double intensity_cell = intensity_abs + intensity_sca;

            if (Q->verbose) {
                printf("%2i  %7.5f  %i  %7.4f  %7.4f\n", i, tau_cell, 
                        R->l[0], R->r[0]/au, s*R->n[0]/au);
                printf("%11.1e  %i  %7.4f  %7.4f\n", R->intensity[inu], R->l[1],
                        R->r[1]/au, s*R->n[1]/au);
                printf("%11.5f  %i  %7.4f  %7.4f\n", R->tau[inu], R->l[2], 
                        R->r[2]/au, s*R->n[2]/au);
            }

            R->intensity[inu] += intensity_cell*exp(-R->tau[inu]);
            R->tau[inu] += tau_cell;
        }

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

        for (int inu = 0; inu < R->nnu; inu++) {

            double tau_cell = 0;
            for (int idust=0; idust<nspecies; idust++)
                tau_cell += s*R->current_kext[idust][inu]*
                    dens[idust][R->l[0]][R->l[1]][R->l[2]];

            if (Q->verbose) {
                printf("%2i  %7.5f  %i  %7.4f  %7.4f\n", i, tau_cell, 
                        R->l[0], R->r[0]/au, s*R->n[0]/au);
                printf("%11.1e  %i  %7.4f  %7.4f\n", R->intensity[inu], R->l[1],
                        R->r[1]/au, s*R->n[1]/au);
                printf("%11.5f  %i  %7.4f  %7.4f\n", R->tau[inu], R->l[2], 
                        R->r[2]/au, s*R->n[2]/au);
            }

            R->intensity[inu] *= exp(-tau_cell);
        }

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

double Grid::smallest_wall_size(Ray *R) {
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
    for (int idust=0; idust<nspecies; idust++) {
        #ifdef _OPENMP
        omp_set_lock(&(lock[idust][l[0]][l[1]][l[2]]));
        #endif

        bool not_converged = true;

        double total_energy = 0;
        for (int ithread = 0; ithread < (int) energy.size(); ithread++)
            total_energy += energy[ithread][idust][l[0]][l[1]][l[2]];

        while (not_converged) {
            double T_old = temp[idust][l[0]][l[1]][l[2]];

            temp[idust][l[0]][l[1]][l[2]]=pow(total_energy/
                (4*sigma*dust[idust]->
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

        #ifdef _OPENMP
        omp_unset_lock(&(lock[idust][l[0]][l[1]][l[2]]));
        #endif
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
        dust[idust]->opacity(nu)*(1. - dust[idust]->albdo(nu))*
        planck_function(nu, temp[idust][ix][iy][iz]);
}

/* Calculate the line profile of a spectral line. */

Vector<double, 3> Grid::vector_velocity(int igas, Photon *P) {
    return Vector<double, 3>(0., 0., 0.);
}

double Grid::line_profile(int igas, int iline, Vector<int, 3> l, double nu) {
    double gamma_thermal = gas[igas]->nu[iline] / c_l * sqrt(2 * k_B * 
            gas_temp[igas][l[0]][l[1]][l[2]] / (gas[igas]->mu * m_p));

    double profile = exp(-(nu - gas[igas]->nu[iline])*(nu - 
            gas[igas]->nu[iline]) / (gamma_thermal * gamma_thermal)) / 
            (gamma_thermal * sqrt(pi));

    return profile;
}
