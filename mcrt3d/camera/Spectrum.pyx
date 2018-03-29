import cython

import numpy
cimport numpy

from ..mcrt3d cimport Spectrum

cdef class SpectrumObj:
    def __init__(self, double r, double incl, double pa, \
            numpy.ndarray[double, ndim=1, mode="c"] nu, \
            double pixel_size, int nnu):

        self.r = r
        self.incl = incl
        self.pa = pa
        self.nu = nu
        self.pixel_size = pixel_size
        self.nnu = nnu

        cdef numpy.ndarray[double, ndim=1, mode="c"] intensity = \
                numpy.zeros((nnu,),dtype=float)

        self.intensity = intensity

        self.obj = new Spectrum(self.r, self.incl, self.pa, &intensity[0], \
                &nu[0], self.pixel_size, self.nnu)

    property intensity:
        def __get__(self):
            return self.intensity

    property nu:
        def __get__(self):
            return self.nu

    def __del__(self):
        del self.obj

