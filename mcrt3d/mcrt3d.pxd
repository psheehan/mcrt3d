import cython

from libcpp cimport bool
from libcpp.vector cimport vector

import numpy as np
cimport numpy as np

cdef extern from "../include/params.h":
    cppclass Params:
        int nphot
        bool bw
        bool scattering
        double scattering_nu
        bool verbose
        bool use_mrw
        double mrw_gamma

        void set_nphot(int _nphot)
        void set_bw(bool _bw)
        void set_scattering(bool _scattering)
        void set_verbose(bool _verbose)
        void set_mrw(bool _use_mrw)
        void set_mrw_gamma(double _mrw_gamma)

cdef extern from "../include/dust.h":
    cppclass Dust:
        int nlam
        double *nu
        double *lam
        double *kabs
        double *ksca
        double *kext
        double *albedo
        double *dkextdnu
        double *dalbedodnu

        int ntemp
        double *temp
        double *planck_opacity
        double *dplanck_opacity_dT
        double *rosseland_extinction
        double *drosseland_extinction_dT
        double **random_nu_CPD
        double **random_nu_CPD_bw
        double **drandom_nu_CPD_dT
        double **drandom_nu_CPD_bw_dT

        Dust(int _nlam, double *_nu, double *_lam, double *_kabs, 
                double *_ksca, double *_kext, double *_albedo)
            
        void set_lookup_tables(int _ntemp, double *_temp, 
                double *_planck_opacity, double *_rosseland_extinction, 
                double *_dplanck_opacity_dT, double *_drosseland_extinction_dT,
                double *_dkextdnu, double *dalbedodnu, double *_random_nu_CPD, 
                double *_random_nu_CPD_bw, double *_drandom_nu_CPD_dT, 
                double *_drandom_nu_CPD_bw_dT)

cdef extern from "../include/isotropic_dust.h":
    cppclass IsotropicDust(Dust):
        IsotropicDust(int _nlam, double *_nu, double *_lam, double *_kabs, 
                double *_ksca, double *_kext, double *_albedo)
            

cdef extern from "../include/source.h":
    cppclass Source:
        double *nu
        double *Bnu
        int nnu

cdef extern from "../include/grid.h":
    cppclass Grid:
        int n1
        int n2
        int n3
        int nw1
        int nw2
        int nw3
        double *w1
        double *w2
        double *w3

        #std::vector<double***> dens
        #std::vector<double***> energy
        #std::vector<double***> temp
        #std::vector<double***> mass
        #std::vector<double****> scatt
        vector[double***] dens
        vector[double***] energy
        vector[double***] temp
        vector[double***] mass
        vector[double****] scatt
        double ***volume

        int nspecies
        #std::vector<Dust*> dust
        vector[Dust*] dust

        int nsources
        #std::vector<Source*> sources
        vector[Source*] sources
        double total_lum

        Params *Q

        int ny
        double *y
        double *f
        double *dydf

        Grid(int _n1, int _n2, int _n3, int _nw1, int _nw2, int _nw3, 
                double *_w1, double *_w2, double *_w3, double *_volume)

        void add_density(double *_dens, double *_temp, double *_mass, 
                Dust *D)
        void add_source(Source *S)
        void set_mrw_tables(double *y, double *f, double *dydf, int ny)

cdef extern from "../include/cartesian_grid.h":
    cppclass CartesianGrid(Grid)

cdef extern from "../include/cylindrical_grid.h":
    cppclass CylindricalGrid(Grid)

cdef extern from "../include/spherical_grid.h":
    cppclass SphericalGrid(Grid)

cdef extern from "../include/mcrt3d.h":
    cppclass MCRT:
        Grid *G
        Params *Q

        MCRT(Grid *_G, Params *_Q)

        void thermal_mc()
        void scattering_mc()
        void mc_iteration()
