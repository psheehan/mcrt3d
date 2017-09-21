import cython

import numpy
cimport numpy

from ..mcrt3d cimport Camera, ParamsObj
from ..grid.Grid cimport GridObj
from .Image cimport ImageObj

cdef class CameraObj:
    cdef Camera *obj
    cdef GridObj grid
    cdef ParamsObj params

    def __init__(self, GridObj G, ParamsObj Q):
        self.grid = G
        self.params = Q

        self.obj = new Camera(G.obj, Q.obj)

    def __del__(self):
        del self.obj

    def make_image(self, ImageObj image):
        self.obj.make_image(image.obj)
