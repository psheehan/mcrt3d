#include "params.h"

Params::Params() {
    nphot = 1000000;
    bw = true;
    verbose = false;
    use_mrw = false;
    mrw_gamma = 4;
    scattering = false;
}

void Params::set_nphot(int _nphot) {
    nphot = _nphot;
}

void Params::set_bw(bool _bw) {
    bw = _bw;
}

void Params::set_scattering(bool _scattering) {
    scattering = _scattering;
}

void Params::set_verbose(bool _verbose) {
    verbose = _verbose;
}

void Params::set_mrw(bool _use_mrw) {
    use_mrw = _use_mrw;
}

void Params::set_mrw_gamma(double _mrw_gamma) {
    mrw_gamma = _mrw_gamma;
}
