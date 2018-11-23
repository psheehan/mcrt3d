#include "mcrt3d.h"

MCRT::MCRT(Grid *_G, Params *_Q) {
    G = _G;

    Q = _Q;
    G->Q = _Q;

    C = new Camera(G, Q);
}

/* Run a Monte Carlo simulation to calculate the temperature throughout the 
 * grid. */

void MCRT::thermal_mc() {
    if (Q->bw)
        mc_iteration();
    else {
        std::vector<double***> told = create4DArr(G->nspecies, G->n1,
                G->n2, G->n3);
        std::vector<double***> treallyold = create4DArr(G->nspecies, G->n1,
                G->n2, G->n3);

        int maxniter = 10;

        equate4DArrs(told, G->temp, G->nspecies, G->n1, G->n2, G->n3);

        int i = 1;
        while (i <= maxniter) {
            printf("Starting iteration # %i \n\n", i);

            equate4DArrs(treallyold, told, G->nspecies, G->n1, G->n2, G->n3);
            equate4DArrs(told, G->temp, G->nspecies, G->n1, G->n2, G->n3);

            mc_iteration();

            G->update_grid();
            set4DArrValue(G->energy, 0.0, G->nspecies, G->n1, G->n2, G->n3);

            if (i > 2)
                if (converged(G->temp, told, treallyold, G->nspecies, G->n1, 
                            G->n2, G->n3))
                    i = maxniter;

            i++;
            printf("\n");
        }

        for (int i=0; i<G->nspecies; i++) {
            delete told[i];
            delete treallyold[i];
        }
    }
}

void MCRT::scattering_mc() {
    // Make sure we've turned the scattering simulation option on.
    bool old_scattering = Q->scattering;
    Q->scattering = true;

    // Set the Grid's scattering array.
    G->initialize_scattering_array();

    // Run the simulation for every frequency bin.
    for (int inu=0; inu<Q->nnu; inu++) {
        printf("inu = %i\n", inu);
        // Set up the right parameters.
        Q->inu = inu;
        Q->nu = Q->scattering_nu[inu];
        Q->dnu = abs(Q->scattering_nu[inu+1] - Q->scattering_nu[inu]);

        G->initialize_luminosity_array(Q->nu);
        mc_iteration();
        G->deallocate_luminosity_array();
    }

    // Reset the scattering simulation to what it was before.
    Q->scattering = old_scattering;
}

void MCRT::mc_iteration() {
    for (int i=0; i<Q->nphot; i++) {
        if (fmod(i+1,Q->nphot/10) == 0) printf("%i\n",i+1);

        Photon *P = G->emit(i);

        if (Q->verbose) {
            printf("Emitting photon # %i\n", i);
            printf("Emitted with direction: %f  %f  %f\n", P->n[0], P->n[1], 
                    P->n[2]);
            printf("Emitted from a cell with temperature: %f\n", 
                    G->temp[0][P->l[0]][P->l[1]][P->l[2]]);
            printf("Emitted with frequency: %e\n", P->nu);
        }

        G->propagate_photon_full(P);

        P->clean();
        delete P;
        if (Q->verbose) printf("Photon has escaped the grid.\n\n");
    }
}

void MCRT::run_image(Image *I) {
    // Set some parameters that are going to be needed.
    Q->scattering_nu = I->nu;
    Q->nnu = I->nnu;

    // Run a scattering simulation.
    scattering_mc();

    // Now, run the image through the camera.
    C->make_image(I);

    // Clean up the appropriate grid parameters.
    G->deallocate_scattering_array();
}

void MCRT::run_spectrum(Spectrum *S) {
    // Set some parameters that are going to be needed.
    Q->scattering_nu = S->nu;
    Q->nnu = S->nnu;

    // Run a scattering simulation.
    scattering_mc();

    // Now, run the image through the camera.
    C->make_spectrum(S);

    // Clean up the appropriate grid parameters.
    G->deallocate_scattering_array();
}
