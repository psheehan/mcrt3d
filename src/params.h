#ifndef PARAMS_H
#define PARAMS_H

struct Params {
    int nphot;
    bool bw;
    bool scattering;
    double scattering_nu;
    bool verbose;
    bool use_mrw;
    double mrw_gamma;

    void set_nphot(int _nphot);
    void set_bw(bool _bw);
    void set_scattering(bool _scattering);
    void set_verbose(bool _verbose);
    void set_mrw(bool _use_mrw);
    void set_mrw_gamma(double _mrw_gamma);
};

#endif
