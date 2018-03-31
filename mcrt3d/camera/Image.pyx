import cython

import numpy
cimport numpy

from ..mcrt3d cimport Image

from ..constants.physics import c

cdef class ImageObj:

    def __init__(self, double r, double incl, double pa, \
            numpy.ndarray[double, ndim=1, mode="c"] x, \
            numpy.ndarray[double, ndim=1, mode="c"] y, \
            int nx, int ny, numpy.ndarray[double, ndim=1, mode="c"] lam, \
            double pixel_size, int nlam):

        self.r = r
        self.incl = incl
        self.pa = pa
        self.x = x
        self.y = y
        self.nx = nx
        self.ny = ny
        self.lam = lam
        self.pixel_size = pixel_size
        self.nlam = nlam
        self.nnu = nlam

        cdef numpy.ndarray[double, ndim=1, mode="c"] nu = c / (lam * 1.0e-4)
        cdef numpy.ndarray[double, ndim=3, mode="c"] intensity = \
                numpy.zeros((nx, ny, self.nnu),dtype=float)

        self.nu = nu
        self.intensity = intensity

        self.obj = new Image(self.r, self.incl, self.pa, &x[0], \
                &y[0], &intensity[0,0,0], self.nx, self.ny, &nu[0], \
                self.pixel_size, self.nnu)

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

