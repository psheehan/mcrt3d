#include "mcrt3d.h"

MCRT::MCRT() {
    Q = new Params();
}

MCRT::~MCRT() {
    delete G; delete C; delete Q;
}

void MCRT::set_cartesian_grid(py::array_t<double> x, py::array_t<double> y,
        py::array_t<double> z) {
    G = new CartesianGrid(x, y, z);
    C = new Camera(G, Q);
}

void MCRT::set_cylindrical_grid(py::array_t<double> r, py::array_t<double> phi,
        py::array_t<double> z) {
    G = new CylindricalGrid(r, phi, z);
    C = new Camera(G, Q);
}

void MCRT::set_spherical_grid(py::array_t<double> r, 
        py::array_t<double> theta, py::array_t<double> phi) {
    G = new SphericalGrid(r, theta, phi);
    C = new Camera(G, Q);
}

/* Run a Monte Carlo simulation to calculate the temperature throughout the 
 * grid. */

void MCRT::thermal_mc(int nphot, bool bw, bool use_mrw, int mrw_gamma, 
        bool verbose) {
    // Make sure parameters are set properly.

    Q->nphot = nphot; 
    Q->bw = bw; 
    Q->use_mrw = use_mrw; 
    Q->mrw_gamma = mrw_gamma;
    Q->verbose = verbose;
    Q->scattering = false;

    // Do the thermal calculation.
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

        delete4DArr(told, G->nspecies, G->n1, G->n2, G->n3);
        delete4DArr(treallyold, G->nspecies, G->n1, G->n2, G->n3);
    }
}

void MCRT::scattering_mc(int nphot, bool verbose) {
    // Make sure parameters are set properly.

    Q->nphot = nphot; 
    Q->use_mrw = false; 
    Q->verbose = verbose;

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

        if (Q->scattering)
            G->propagate_photon_scattering(P);
        else
            G->propagate_photon_full(P);

        delete P;
        if (Q->verbose) printf("Photon has escaped the grid.\n\n");
    }
}

void MCRT::run_image(int nx, int ny, double pixel_size, 
        py::array_t<double> __lam, int nphot, double incl, double pa, 
        double dpc) {

    // Set some parameters that are going to be needed.
    auto _lam_buf = __lam.request();
    if (_lam_buf.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one");
    double *lam = (double *) _lam_buf.ptr;

    Q->nnu = _lam_buf.shape[0];
    Q->scattering_nu = new double[Q->nnu];

    for (int i = 0; i < Q->nnu; i++)
        Q->scattering_nu[i] = c_l / lam[i];

    // Run a scattering simulation.
    scattering_mc(Q->nphot, Q->verbose);

    // Now, run the image through the camera.
    Image *I = C->make_image(nx, ny, pixel_size, __lam, incl, pa, dpc);

    images.append(I);

    // Clean up the appropriate grid parameters.
    G->deallocate_scattering_array();

    Q->nnu = 0;
    delete[] Q->scattering_nu;
}

void MCRT::run_spectrum(py::array_t<double> __lam, int nphot, double incl, 
            double pa, double dpc) {
    // Set some parameters that are going to be needed.

    auto _lam_buf = __lam.request();
    if (_lam_buf.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one");
    double *lam = (double *) _lam_buf.ptr;

    Q-> nnu = _lam_buf.shape[0];
    Q->scattering_nu = new double[Q->nnu];

    for (int i = 0; i < Q->nnu; i++)
        Q->scattering_nu[i] = c_l / lam[i];

    // Run a scattering simulation.
    scattering_mc(Q->nphot, Q->verbose);

    // Now, run the image through the camera.
    Spectrum *S = C->make_spectrum(__lam, incl, pa, dpc);

    spectra.append(S);

    // Clean up the appropriate grid parameters.
    G->deallocate_scattering_array();

    Q->nnu = 0;
    delete[] Q->scattering_nu;
}

PYBIND11_MODULE(mcrt3d, m) {
    py::class_<Dust>(m, "Dust")
        .def(py::init<py::array_t<double>, py::array_t<double>, 
                py::array_t<double>>())
        .def_readonly("lam", &Dust::_lam)
        .def_readonly("nu", &Dust::_nu)
        .def_readonly("kabs", &Dust::_kabs)
        .def_readonly("ksca", &Dust::_ksca)
        .def_readonly("kext", &Dust::_kext)
        .def_readonly("albedo", &Dust::_albedo);

    py::class_<IsotropicDust, Dust>(m, "IsotropicDust")
        .def(py::init<py::array_t<double>, py::array_t<double>, 
                py::array_t<double>>());

    py::class_<Grid>(m, "Grid")
        .def_readonly("volume", &Grid::_volume)
        .def_readonly("density", &Grid::_dens)
        .def_readonly("temperature", &Grid::_temp)
        .def_readonly("dust", &Grid::_dust)
        .def_readonly("scatt", &Grid::_scatt)
        .def("add_density", &Grid::add_density, 
                "Add a density layer to the Grid.")
        .def("add_star", &Grid::add_star, "Add a star to the Grid.", 
                py::arg("x")=0., py::arg("y")=0., py::arg("z")=0., 
                py::arg("mass")=1.989e33, py::arg("radius")=69.634e9, 
                py::arg("temperature")=4000.);

    py::class_<CartesianGrid, Grid>(m, "CartesianGrid")
        .def_readonly("x", &CartesianGrid::x)
        .def_readonly("y", &CartesianGrid::y)
        .def_readonly("z", &CartesianGrid::z);

    py::class_<CylindricalGrid, Grid>(m, "CylindricalGrid")
        .def_readonly("r", &CylindricalGrid::r)
        .def_readonly("phi", &CylindricalGrid::phi)
        .def_readonly("z", &CylindricalGrid::z);

    py::class_<SphericalGrid, Grid>(m, "SphericalGrid")
        .def_readonly("r", &SphericalGrid::r)
        .def_readonly("theta", &SphericalGrid::theta)
        .def_readonly("phi", &SphericalGrid::phi);

    py::class_<MCRT>(m, "MCRT")
        .def(py::init<>())
        .def_readonly("grid", &MCRT::G)
        .def("set_cartesian_grid", &MCRT::set_cartesian_grid,
                "Setup a grid in cartesian coordinates.")
        .def("set_cylindrical_grid", &MCRT::set_cylindrical_grid,
                "Setup a grid in cylindrical coordinates.")
        .def("set_spherical_grid", &MCRT::set_spherical_grid,
                "Setup a grid in spherical coordinates.")
        .def("thermal_mc", &MCRT::thermal_mc, 
                "Calculate the temperature throughout the grid.",
                py::arg("nphot")=1000000, py::arg("bw")=true, 
                py::arg("use_mrw")=false, py::arg("mrw_gamma")=4, 
                py::arg("verbose")=false)
        .def("run_image", &MCRT::run_image, "Generate an image.")
        .def("run_spectrum", &MCRT::run_spectrum, "Generate a spectrum.");
}
