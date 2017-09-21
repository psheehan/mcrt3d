import cython

import numpy
cimport numpy

from ..mcrt3d cimport Image

cdef class ImageObj:
    cdef Image *obj
    cdef int nx, ny, nnu
    cdef double r, incl, pa, pixel_size
    cdef numpy.ndarray x, y, nu, intensity
