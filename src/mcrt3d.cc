#include "mcrt3d.h"

MCRT::MCRT(Grid *_G, Params *_Q) {
    printf("Hello\n");
    G = _G;
    printf("This was ok...\n");

    Q = _Q;
    printf("So was this...\n");
    G->Q = _Q;
    printf("If we got here, then something strange is happening...\n");
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
    mc_iteration();
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
