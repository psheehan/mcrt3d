import cython

import numpy
cimport numpy

from ..mcrt3d cimport Camera, ParamsObj
from ..grid.Grid cimport GridObj
from .Image cimport ImageObj

cdef class CameraObj:
    #def __init__(self, GridObj G, ParamsObj Q):
    #    self.grid = G
    #    self.params = Q

    #    self.obj = new Camera(G.obj, Q.obj)

    #def __del__(self):
    #    del self.obj

    #def make_image(self, ImageObj image):
    #    self.obj.make_image(image.obj)

    property nx:
        def __get__(self):
            return self.nx
        def __set__(self, int val):
            self.nx = val

    property ny:
        def __get__(self):
            return self.ny
        def __set__(self, int val):
            self.ny = val

    property nlam:
        def __get__(self):
            return self.nlam
        def __set__(self, int val):
            self.nlam = val

    property pixel_size:
        def __get__(self):
            return self.pixel_size
        def __set__(self, double val):
            self.pixel_size = val

    property lam:
        def __get__(self):
            return self.lam
        def __set__(self, numpy.ndarray[double, ndim=1, mode="c"] val):
            self.lam = val
            self.nlam = len(val)

    property x:
        def __get__(self):
            return self.x
        def __set__(self, numpy.ndarray[double, ndim=1, mode="c"] val):
            self.x = val

    property y:
        def __get__(self):
            return self.y
        def __set__(self, numpy.ndarray[double, ndim=1, mode="c"] val):
            self.y = val
