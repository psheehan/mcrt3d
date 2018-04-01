import cython

import numpy
cimport numpy

from ..mcrt3d cimport Camera, ParamsObj
from ..grid.Grid cimport GridObj
from .Image cimport ImageObj

cdef class CameraObj:
    #cdef Camera *obj
    #cdef GridObj grid
    #cdef ParamsObj params

    cdef int nx, ny, nlam
    cdef double pixel_size
    cdef numpy.ndarray x, y, lam
