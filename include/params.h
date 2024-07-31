#ifndef PARAMS_H
#define PARAMS_H

struct Params {
    // Constructor class.

    Params ();

    // General purpose parameters.
    int nphot;
    bool bw;
    bool verbose;
    bool use_mrw;
    double mrw_gamma;

    // Parameters for doing the scattering simulation.
    bool scattering;
    int inu;
    int nnu;
    double nu;
    double dnu;
    Kokkos::View<double*> scattering_nu{"scatering_nu", 0};

    // Parameters for imaging.

    bool raytrace_dust;
    bool raytrace_gas;

    // Functions to set various properties.
    void set_nphot(int _nphot);
    void set_bw(bool _bw);
    void set_scattering(bool _scattering);
    void set_verbose(bool _verbose);
    void set_mrw(bool _use_mrw);
    void set_mrw_gamma(double _mrw_gamma);
};

#endif
