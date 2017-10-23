import cython

import numpy
cimport numpy

from ..mcrt3d cimport Source, Star

cdef class StarObj:
    cdef Star *obj

    cdef int nnu
    cdef double x, y, z, mass, temperature, radius, luminosity
    cdef numpy.ndarray nu, Bnu, random_nu_CPD
