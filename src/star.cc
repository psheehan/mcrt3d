#include "star.h"

/* Functions to set up the sources. */

Star::Star(double x, double y, double z, double _mass, double _radius, \
        double _temperature) {

    r[0] = x;
    r[1] = y;
    r[2] = z;
    mass = _mass;
    radius = _radius;
    temperature = _temperature;
}

void Star::set_blackbody_spectrum(int _nnu, double *_nu, double *_Bnu, 
        double _luminosity, double *_random_nu_CPD) {

    nnu = _nnu;
    nu = _nu;
    Bnu = _Bnu;
    luminosity = _luminosity;
    random_nu_CPD = _random_nu_CPD;
}

/* Emit a photon from the source. */

Photon *Star::emit(int nphot) {
    Photon *P = new Photon();

    double theta = pi*random_number();
    double phi = 2*pi*random_number();

    P->r[0] = radius*sin(theta)*cos(phi);
    P->r[1] = radius*sin(theta)*sin(phi);
    P->r[2] = radius*cos(theta);

    Vector<double, 3> r_hat, theta_hat, phi_hat;

    r_hat[0] = sin(theta)*cos(phi);
    r_hat[1] = sin(theta)*sin(phi);
    r_hat[2] = cos(theta);
    theta_hat[0] = cos(theta)*cos(phi);
    theta_hat[1] = cos(theta)*sin(phi);
    theta_hat[2] = -sin(theta);
    phi_hat[0] = -sin(phi);
    phi_hat[1] = cos(phi);
    phi_hat[2] = 0;

    double cost = random_number();
    double sint = sqrt(1-pow(cost,2));
    phi = 2*pi*random_number();

    P->n = cost*r_hat + sint*cos(phi)*phi_hat + sint*sin(phi)*theta_hat;

    P->invn[0] = 1.0/P->n[0];
    P->invn[1] = 1.0/P->n[1];
    P->invn[2] = 1.0/P->n[2];
    P->l[0] = -1;
    P->l[1] = -1;
    P->l[2] = -1;

    P->energy = luminosity / nphot;

    P->nu = random_nu();

    return P;
};

Photon *Star::emit(double _nu, double _dnu, int nphot) {
    Photon *P = new Photon();

    double theta = pi*random_number();
    double phi = 2*pi*random_number();

    P->r[0] = radius*sin(theta)*cos(phi);
    P->r[1] = radius*sin(theta)*sin(phi);
    P->r[2] = radius*cos(theta);

    Vector<double, 3> r_hat, theta_hat, phi_hat;

    r_hat[0] = sin(theta)*cos(phi);
    r_hat[1] = sin(theta)*sin(phi);
    r_hat[2] = cos(theta);
    theta_hat[0] = cos(theta)*cos(phi);
    theta_hat[1] = cos(theta)*sin(phi);
    theta_hat[2] = -sin(theta);
    phi_hat[0] = -sin(phi);
    phi_hat[1] = cos(phi);
    phi_hat[2] = 0;

    double cost = random_number();
    double sint = sqrt(1-pow(cost,2));
    phi = 2*pi*random_number();

    P->n = cost*r_hat + sint*cos(phi)*phi_hat + sint*sin(phi)*theta_hat;

    P->invn[0] = 1.0/P->n[0];
    P->invn[1] = 1.0/P->n[1];
    P->invn[2] = 1.0/P->n[2];
    P->l[0] = -1;
    P->l[1] = -1;
    P->l[2] = -1;

    P->nu = _nu;

    P->energy = 4*pi*pi*radius*radius*flux(_nu) / nphot;

    return P;
};

Ray *Star::emit_ray(double _nu, double _dnu, double _pixelsize, \
        Vector<double, 3> _n, int nphot) {
    Ray *R = new Ray();

    double theta = pi*random_number();
    double phi = 2*pi*random_number();

    R->r[0] = radius*sin(theta)*cos(phi);
    R->r[1] = radius*sin(theta)*sin(phi);
    R->r[2] = radius*cos(theta);

    R->n = -_n;

    R->invn[0] = 1.0/R->n[0];
    R->invn[1] = 1.0/R->n[1];
    R->invn[2] = 1.0/R->n[2];
    R->l[0] = -1;
    R->l[1] = -1;
    R->l[2] = -1;

    R->nu = _nu;

    // Now get the total energy of the photon.
    R->intensity = flux(_nu) * pi*radius*radius / (_pixelsize * _pixelsize) /
        nphot;

    R->tau = 0.0;
    R->pixel_size = _pixelsize;
    R->pixel_too_large = false;

    return R;
};

/* Get a random frequency drawn from the spectrum of the source. */

double Star::random_nu() {
    double freq;
    double ksi = random_number();

    for (int i=1; i<nnu; i++) {
        if (random_nu_CPD[i] > ksi) {
            freq = random_number() * (nu[i] - nu[i-1]) + nu[i-1];
            break;
        }
    }

    return freq;
};

double Star::intercept_distance(Photon *P) {
    double s = HUGE_VAL;

    double r = P->r.norm();

    if (!equal(r, radius, EPSILON)) {
        double b = P->r*P->n;
        double c = r*r - radius*radius;
        double d = b*b - c;

        if (d >= 0) {
            double sr1 = -b + sqrt(d);
            if ((sr1 < s) && (sr1 > 0)) s = sr1;
            double sr2 = -b - sqrt(d);
            if ((sr2 < s) && (sr2 > 0)) s = sr2;
        }
    }

    return s;
}

double Star::flux(double freq) {
    int l = find_in_arr(freq, nu, nnu);

    double dBnudnu = (Bnu[l+1] - Bnu[l])/(nu[l+1] - nu[l]);

    double flux = dBnudnu*(freq-nu[l])+Bnu[l];

    return flux;
};
