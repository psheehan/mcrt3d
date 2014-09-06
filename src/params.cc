#ifndef PARAMS_CC
#define PARAMS_CC

struct Params {
    int nphot;
    bool bw;
    bool scattering;
    double scattering_nu;
    bool verbose;
    bool use_mrw;
    double mrw_gamma;
};

#endif
