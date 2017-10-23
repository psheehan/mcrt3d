import cython

import numpy
cimport numpy

from ..mcrt3d cimport Image

cdef class ImageObj:
    #cdef Image *obj
    #cdef int nx, ny, nnu
    #cdef double r, incl, pa, pixel_size
    #cdef numpy.ndarray x, y, nu, intensity

    def __init__(self, double r, double incl, double pa, \
            numpy.ndarray[double, ndim=1, mode="c"] x, \
            numpy.ndarray[double, ndim=1, mode="c"] y, \
            int nx, int ny, numpy.ndarray[double, ndim=1, mode="c"] nu, \
            double pixel_size, int nnu):

        self.r = r
        self.incl = incl
        self.pa = pa
        self.x = x
        self.y = y
        self.nx = nx
        self.ny = ny
        self.nu = nu
        self.pixel_size = pixel_size
        self.nnu = nnu

        cdef numpy.ndarray[double, ndim=3, mode="c"] intensity = \
                numpy.zeros((nx, ny, nnu),dtype=float)

        self.intensity = intensity

        print("Hello!")
        self.obj = new Image(self.r, self.incl, self.pa, &x[0], \
                &y[0], &intensity[0,0,0], self.nx, self.ny, &nu[0], \
                self.pixel_size, self.nnu)
        print("Goodbye")

    def __del__(self):
        del self.obj

