#include "dust.h"

/* Functions to set up the dust. */

Dust::Dust(py::array_t<double> __lam, py::array_t<double> __kabs,
            py::array_t<double> __ksca) {

    random_pool = new Kokkos::Random_XorShift64_Pool<>(/*seed=*/12345);

    auto __lam_buf = __lam.request();
    nlam = __lam_buf.shape[0];

    Kokkos::resize(lam, nlam);
    Kokkos::deep_copy(lam, view_from_array(__lam));
    Kokkos::resize(kabs, nlam);
    Kokkos::deep_copy(kabs, view_from_array(__kabs));
    Kokkos::resize(ksca, nlam);
    Kokkos::deep_copy(ksca, view_from_array(__ksca));

    _lam = array_from_view(lam, 1, {(size_t) nlam});
    _kabs = array_from_view(kabs, 1, {(size_t) nlam});
    _ksca = array_from_view(ksca, 1, {(size_t) nlam});

    // Set up the volume of each cell.

    Kokkos::resize(nu, nlam);
    Kokkos::resize(kext, nlam);
    Kokkos::resize(albedo, nlam);

    _nu = array_from_view(nu, 1, {(size_t) nlam});
    _kext = array_from_view(kext, 1, {(size_t) nlam});
    _albedo = array_from_view(albedo, 1, {(size_t) nlam});

    // Finally, calculate their values.

    for (int i = 0; i < nlam; i++) {
        nu(i) = c_l / lam(i); 
        kext(i) = kabs(i) + ksca(i);
        albedo(i) = ksca(i) / kext(i); 
    }

    // Catch if nu is not increasing.

    if (nu(1) - nu(0) <= 0.)
        throw std::runtime_error("nu must be monotonically increasing.");

    // Finally, make the lookup tables;
    set_lookup_tables();
}

/*Dust::Dust(int _nlam, double *_nu, double *_lam, double *_kabs, 
        double *_ksca, double *_kext, double *_albedo) {
    
    nlam = _nlam;
    nu = _nu;
    lam = _lam;
    kabs = _kabs;
    ksca = _ksca;
    kext = _kext;
    albedo = _albedo;
}*/

Dust::~Dust() {
}

void Dust::set_lookup_tables() {
    // Create the temperature array;
    ntemp = 1000;
    Kokkos::resize(temp, ntemp);

    for (int i = 0; i < ntemp; i++)
        temp(i) = pow(10.,-1.+i*6./(ntemp-1));

    // Calculate the derivatives of kext and albedo.
    Kokkos::resize(dkextdnu, nlam-1);
    Kokkos::deep_copy(dkextdnu, derivative(kext, nu, nlam));
    Kokkos::resize(dalbedodnu, nlam-1);
    Kokkos::deep_copy(dalbedodnu, derivative(albedo, nu, nlam));

    // Calculate the Planck Mean Opacity and its derivative.
    Kokkos::View<double**> tmp_planck("tmp_planck", ntemp, nlam);
    for (int i = 0; i < ntemp; i++)
        for (int j = 0; j < nlam; j++)
            tmp_planck(i,j) = planck_function(nu(j), temp(i)) * kabs(j);

    Kokkos::resize(planck_opacity, ntemp);
    for (int i = 0; i < ntemp; i++)
        planck_opacity(i) = pi / (sigma * pow(temp(i),4)) * 
                integrate(Kokkos::subview(tmp_planck, i, Kokkos::ALL), nu, nlam);

    Kokkos::resize(dplanck_opacity_dT, ntemp-1);
    Kokkos::deep_copy(dplanck_opacity_dT, derivative(planck_opacity, temp, ntemp));

    // Calculate the Planck Mean Opacity and its derivative.
    Kokkos::View<double**> tmp_ross("tmp_ross", ntemp, nlam);
    Kokkos::View<double**> tmp_ross_num("tmp_ross_num", ntemp, nlam);
    for (int i = 0; i < ntemp; i++) {
        for (int j = 0; j < nlam; j++) {
            tmp_ross(i,j) = planck_function_derivative(nu(j), temp(i)) / kext(j);
            tmp_ross_num(i,j) = planck_function_derivative(nu(j), temp(i));
        }
    }

    Kokkos::resize(rosseland_extinction, ntemp);
    for (int i = 0; i < ntemp; i++)
        rosseland_extinction(i) = integrate(Kokkos::subview(tmp_ross_num, i, Kokkos::ALL), nu, nlam) / 
                integrate(Kokkos::subview(tmp_ross, i, Kokkos::ALL), nu, nlam);

    Kokkos::resize(drosseland_extinction_dT, ntemp-1);
    Kokkos::deep_copy(drosseland_extinction_dT, derivative(rosseland_extinction, temp, ntemp));

    // Create the cumulative probability density functions for generating a
    // random nu value, for regular.
    Kokkos::View<double**> tmp("tmp", ntemp, nlam);
    for (int i = 0; i < ntemp; i++)
        for (int j = 0; j < nlam; j++)
            tmp(i,j) = kabs(j) * planck_function(nu(j), temp(i));

    Kokkos::resize(random_nu_CPD, ntemp, nlam);
    for (int i = 0; i < ntemp; i++) {
        auto random_nu_CPD_tmp = cumulative_integrate(Kokkos::subview(tmp, i, Kokkos::ALL), nu, nlam);
        for (int j = 0; j < nlam; j++)
            random_nu_CPD(i,j) = random_nu_CPD_tmp(j);
    }

    Kokkos::resize(drandom_nu_CPD_dT, ntemp-1, nlam);
    Kokkos::deep_copy(drandom_nu_CPD_dT, derivative2D_ax0(random_nu_CPD, temp, ntemp, nlam));

    // Create the cumulative probability density functions for generating a
    // random nu value, for bw.
    for (int i = 0; i < ntemp; i++)
        for (int j = 0; j < nlam; j++)
            tmp(i,j) = kabs(j) * planck_function_derivative(nu(j), temp(i));

    Kokkos::resize(random_nu_CPD_bw, ntemp, nlam);
    for (int i = 0; i < ntemp; i++) {
        auto random_nu_CPD_bw_tmp = cumulative_integrate(Kokkos::subview(tmp, i, Kokkos::ALL), nu, nlam);
        for (int j = 0; j < nlam; j++)
            random_nu_CPD_bw(i,j) = random_nu_CPD_bw_tmp(j);
    }

    Kokkos::resize(drandom_nu_CPD_bw_dT, ntemp-1, nlam);
    Kokkos::deep_copy(drandom_nu_CPD_bw_dT, derivative2D_ax0(random_nu_CPD_bw, temp, ntemp, nlam));
}

/* Scatter a photon isotropically off of dust. */

void Dust::scatter(Photon *P) {
}

/* Absorb and then re-emit a photon from dust. */

void Dust::absorb(Photon *P, double T, bool bw) {
    double cost = -1+2*random_number(random_pool);
    double sint = sqrt(1-pow(cost,2));
    double phi = 2*pi*random_number(random_pool);

    P->n[0] = sint*cos(phi);
    P->n[1] = sint*sin(phi);
    P->n[2] = cost;
    P->invn[0] = 1.0/P->n[0];
    P->invn[1] = 1.0/P->n[1];
    P->invn[2] = 1.0/P->n[2];

    P->nu = random_nu(T,bw);
}

void Dust::absorb_mrw(Photon *P, double T, bool bw) {
    P->nu = random_nu(T,bw);
}

/* Calculate a random frequency for a photon. */

double Dust::random_nu(double T, bool bw) {
    double freq, CPD;

    int i = find_in_arr(T,temp,ntemp);

    double ksi = random_number(random_pool);

    for (int j=1; j < nlam; j++) {
        if (bw)
            CPD = drandom_nu_CPD_bw_dT(i,j) * (T - temp(i)) + 
                random_nu_CPD_bw(i,j);
        else
            CPD = drandom_nu_CPD_dT(i,j) * (T - temp(i)) + 
                random_nu_CPD(i,j);

        if (CPD > ksi) {
            freq = random_number(random_pool) * (nu(j) - nu(j-1)) + nu(j-1);
            break;
        }
    }

    return freq;
}

/* Calculate the opacity of a dust grain at a specific frequency. */

double Dust::opacity(double freq) {
    int l = find_in_arr(freq, nu, nlam);

    double opacity = dkextdnu(l)*(freq-nu(l))+kext(l);

    return opacity;
};

/* Calculate the albedo of a dust grain at a specific frequency. */

double Dust::albdo(double freq) {
    int l = find_in_arr(freq, nu, nlam);

    double albdo = dalbedodnu(l)*(freq-nu(l))+albedo(l);

    return albdo;
}

/* Calculate the Planck Mean Opacity for a dust grain at a given temperature. */

double Dust::planck_mean_opacity(double T) {
    int n = find_in_arr(T,temp,ntemp);

    double planck_mean_opacity = dplanck_opacity_dT(n)*(T-temp(n))+
        planck_opacity(n);

    return planck_mean_opacity;
}

/* Calculate the Rosseland Mean Extinction for a dust grain at a given 
 * temperature. */

double Dust::rosseland_mean_extinction(double T) {
    int n = find_in_arr(T,temp,ntemp);

    double rosseland_mean_extinction = drosseland_extinction_dT(n)*(T-temp(n))+
        rosseland_extinction(n);

    return rosseland_mean_extinction;
}
