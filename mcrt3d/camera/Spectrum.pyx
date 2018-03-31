import cython

import numpy
cimport numpy

from ..mcrt3d cimport Spectrum

from ..constants.physics import c

cdef class SpectrumObj:
    def __init__(self, double r, double incl, double pa, \
            numpy.ndarray[double, ndim=1, mode="c"] lam, \
            double pixel_size, int nlam):

        self.r = r
        self.incl = incl
        self.pa = pa
        self.lam = lam
        self.pixel_size = pixel_size
        self.nnu = nlam
        self.nnu = nlam

        cdef numpy.ndarray[double, ndim=1, mode="c"] nu = c / (lam * 1.0e-4)
        cdef numpy.ndarray[double, ndim=1, mode="c"] intensity = \
                numpy.zeros((self.nnu,),dtype=float)

        self.nu = nu
        self.intensity = intensity

        self.obj = new Spectrum(self.r, self.incl, self.pa, &intensity[0], \
                &nu[0], self.pixel_size, self.nnu)

    property intensity:
        def __get__(self):
            return self.intensity

    property lam:
        def __get__(self):
            return self.lam

    property nu:
        def __get__(self):
            return self.nu

    def __del__(self):
        del self.obj

