#include "dust.h"

/* Functions to set up the dust. */

Dust::Dust(py::array_t<double> __lam, py::array_t<double> __kabs,
            py::array_t<double> __ksca) {

    _lam = __lam; _kabs = __kabs; _ksca = __ksca;

    // Load the array buffers to get the proper setup info.

    auto _lam_buf = __lam.request(); auto _kabs_buf = __kabs.request();
    auto _ksca_buf = __ksca.request();

    if (_lam_buf.ndim != 1 || _kabs_buf.ndim != 1 || _ksca_buf.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one");

    // Now get the correct format.

    nlam = _lam_buf.shape[0];
    lam = (double *) _lam_buf.ptr;
    kabs = (double *) _kabs_buf.ptr;
    ksca = (double *) _ksca_buf.ptr;

    // Set up the volume of each cell.

    _nu = py::array_t<double>(nlam);
    _kext = py::array_t<double>(nlam);
    _albedo = py::array_t<double>(nlam);

    auto _nu_buf = _nu.request(); nu = (double *) _nu_buf.ptr;
    auto _kext_buf = _kext.request(); kext = (double *) _kext_buf.ptr;
    auto _albedo_buf = _albedo.request(); albedo = (double *) _albedo_buf.ptr;

    // Finally, calculate their values.

    for (int i = 0; i < nlam; i++) {
        nu[i] = c_l / lam[i]; 
        kext[i] = kabs[i] + ksca[i];
        albedo[i] = ksca[i] / kext[i];
    }

    // Catch if nu is not increasing.

    if (nu[1] - nu[0] <= 0.)
        throw std::runtime_error("nu must be monotonically increasing.");

    // Finally, make the lookup tables;
    set_lookup_tables();
}

Dust::Dust(int _nlam, double *_nu, double *_lam, double *_kabs, 
        double *_ksca, double *_kext, double *_albedo) {
    
    nlam = _nlam;
    nu = _nu;
    lam = _lam;
    kabs = _kabs;
    ksca = _ksca;
    kext = _kext;
    albedo = _albedo;
}

Dust::~Dust() {
    delete[] temp; delete[] dkextdnu; delete[] dalbedodnu;
    delete[] planck_opacity; delete[] dplanck_opacity_dT;
    delete[] rosseland_extinction; delete[] drosseland_extinction_dT;
    delete2DArr(random_nu_CPD, ntemp, nlam);
    delete2DArr(drandom_nu_CPD_dT, ntemp-1, nlam);
    delete2DArr(random_nu_CPD_bw, ntemp, nlam);
    delete2DArr(drandom_nu_CPD_bw_dT, ntemp-1, nlam);
}

void Dust::set_lookup_tables() {
    // Create the temperature array;
    ntemp = 100;
    temp = new double[ntemp];

    for (int i = 0; i < ntemp; i++)
        temp[i] = pow(10.,-1.+i*6./(ntemp-1));

    // Calculate the derivatives of kext and albedo.
    dkextdnu = derivative(kext, nu, nlam);
    dalbedodnu = derivative(albedo, nu, nlam);

    // Calculate the Planck Mean Opacity and its derivative.
    double **tmp_planck = create2DArr(ntemp, nlam);
    for (int i = 0; i < ntemp; i++)
        for (int j = 0; j < nlam; j++)
            tmp_planck[i][j] = planck_function(nu[j], temp[i]) * kabs[j];

    planck_opacity = new double[ntemp];
    for (int i = 0; i < ntemp; i++)
        planck_opacity[i] = pi / (sigma * pow(temp[i],4)) * 
                integrate(tmp_planck[i], nu, nlam);

    dplanck_opacity_dT = derivative(planck_opacity, temp, ntemp);

    delete2DArr(tmp_planck, ntemp, nlam);

    // Calculate the Planck Mean Opacity and its derivative.
    double **tmp_ross = create2DArr(ntemp, nlam);
    for (int i = 0; i < ntemp; i++)
        for (int j = 0; j < nlam; j++)
            tmp_ross[i][j] = planck_function(nu[j], temp[i]) / kext[j];

    rosseland_extinction = new double[ntemp];
    for (int i = 0; i < ntemp; i++)
        rosseland_extinction[i] = (sigma * pow(temp[i],4) / pi) / 
                integrate(tmp_ross[i], nu, nlam);

    drosseland_extinction_dT = derivative(rosseland_extinction, temp, ntemp);

    delete2DArr(tmp_ross, ntemp, nlam);

    // Create the cumulative probability density functions for generating a
    // random nu value, for regular.
    double **tmp = create2DArr(ntemp, nlam);
    for (int i = 0; i < ntemp; i++)
        for (int j = 0; j < nlam; j++)
            tmp[i][j] = kabs[j] * planck_function(nu[j], temp[i]);

    random_nu_CPD = new double*[ntemp];
    for (int i = 0; i < ntemp; i++)
        random_nu_CPD[i] = cumulative_integrate(tmp[i], nu, nlam);

    drandom_nu_CPD_dT = derivative2D_ax0(random_nu_CPD, temp, ntemp, nlam);

    // Create the cumulative probability density functions for generating a
    // random nu value, for bw.
    for (int i = 0; i < ntemp; i++)
        for (int j = 0; j < nlam; j++)
            tmp[i][j] = kabs[j] * planck_function_derivative(nu[j], temp[i]);

    random_nu_CPD_bw = new double*[ntemp];
    for (int i = 0; i < ntemp; i++)
        random_nu_CPD_bw[i] = cumulative_integrate(tmp[i], nu, nlam);

    drandom_nu_CPD_bw_dT = derivative2D_ax0(random_nu_CPD_bw, temp, ntemp, 
            nlam);

    delete2DArr(tmp, ntemp, nlam);
}

/* Scatter a photon isotropically off of dust. */

void Dust::scatter(Photon *P) {
}

/* Absorb and then re-emit a photon from dust. */

void Dust::absorb(Photon *P, double T, bool bw) {
    double cost = -1+2*random_number();
    double sint = sqrt(1-pow(cost,2));
    double phi = 2*pi*random_number();

    P->n[0] = sint*cos(phi);
    P->n[1] = sint*sin(phi);
    P->n[2] = cost;
    P->invn[0] = 1.0/P->n[0];
    P->invn[1] = 1.0/P->n[1];
    P->invn[2] = 1.0/P->n[2];

    P->nu = random_nu(T,bw);
}

/* Calculate a random frequency for a photon. */

double Dust::random_nu(double T, bool bw) {
    double freq, CPD;

    int i = find_in_arr(T,temp,ntemp);

    double ksi = random_number();

    for (int j=1; j < nlam; j++) {
        if (bw)
            CPD = drandom_nu_CPD_bw_dT[i][j] * (T - temp[i]) + 
                random_nu_CPD_bw[i][j];
        else
            CPD = drandom_nu_CPD_dT[i][j] * (T - temp[i]) + 
                random_nu_CPD[i][j];

        if (CPD > ksi) {
            freq = random_number() * (nu[j] - nu[j-1]) + nu[j-1];
            break;
        }
    }

    return freq;
}

/* Calculate the opacity of a dust grain at a specific frequency. */

double Dust::opacity(double freq) {
    int l = find_in_arr(freq, nu, nlam);

    double opacity = dkextdnu[l]*(freq-nu[l])+kext[l];

    return opacity;
};

/* Calculate the albedo of a dust grain at a specific frequency. */

double Dust::albdo(double freq) {
    int l = find_in_arr(freq, nu, nlam);

    double albdo = dalbedodnu[l]*(freq-nu[l])+albedo[l];

    return albdo;
}

/* Calculate the Planck Mean Opacity for a dust grain at a given temperature. */

double Dust::planck_mean_opacity(double T) {
    int n = find_in_arr(T,temp,ntemp);

    double planck_mean_opacity = dplanck_opacity_dT[n]*(T-temp[n])+
        planck_opacity[n];

    return planck_mean_opacity;
}

/* Calculate the Rosseland Mean Extinction for a dust grain at a given 
 * temperature. */

double Dust::rosseland_mean_extinction(double T) {
    int n = find_in_arr(T,temp,ntemp);

    double rosseland_mean_extinction = drosseland_extinction_dT[n]*(T-temp[n])+
        rosseland_extinction[n];

    return rosseland_mean_extinction;
}
