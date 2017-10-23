import cython

import numpy
cimport numpy

from ..mcrt3d cimport Grid

cdef class GridObj:
    cdef Grid *obj

    cdef str coordsystem
    cdef int n1, n2, n3, nw1, nw2, nw3
    cdef numpy.ndarray w1, w2, w3
    cdef int ny
    cdef numpy.ndarray y, f, dydf

    cdef numpy.ndarray xx, yy, zz
    cdef numpy.ndarray r, rho, phi, theta

    cdef numpy.ndarray volume

    cdef list density, mass, temperature, dust, sources
