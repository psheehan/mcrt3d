import cython

import numpy
cimport numpy

from ..mcrt3d cimport Grid
from .Grid cimport GridObj

cdef class SphericalGridObj(GridObj):
    pass
