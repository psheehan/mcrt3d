import cython

import numpy
cimport numpy

from ..mcrt3d cimport Grid

cdef class GridObj:
    cdef Grid *obj

    cdef int n1, n2, n3, nw1, nw2, nw3
    cdef int ny
    cdef numpy.ndarray y, f, dydf
