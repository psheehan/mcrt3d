import cython

import numpy
cimport numpy

from ..mcrt3d cimport Dust, IsotropicDust

cdef class DustObj:
    cdef Dust *obj

    cdef int nlam, ntemp

    cdef numpy.ndarray lam, nu, kabs, ksca, kext, albedo
    cdef numpy.ndarray temp, planck_opacity, dplanck_opacity_dT, \
            rosseland_extinction, drosseland_extinction_dT, dkextdnu, \
            dalbedodnu, random_nu_CPD, random_nu_CPD_bw, drandom_nu_CPD_dT, \
            drandom_nu_CPD_bw_dT
