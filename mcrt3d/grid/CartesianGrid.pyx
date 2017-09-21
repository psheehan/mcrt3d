import cython

import numpy
cimport numpy

from ..mcrt3d cimport CartesianGrid
from .Grid cimport GridObj

cdef class CartesianGridObj(GridObj):
    def __init__(self, numpy.ndarray[double, ndim=1, mode="c"] w1, \
            numpy.ndarray[double, ndim=1, mode="c"] w2, \
            numpy.ndarray[double, ndim=1, mode="c"] w3):

        self.coordsystem = "cartesian"

        self.x = 0.5*(w1[0:w1.size-1] + w1[1:w1.size])
        self.y = 0.5*(w2[0:w2.size-1] + w2[1:w2.size])
        self.z = 0.5*(w3[0:w3.size-1] + w3[1:w3.size])

        self.nw1 = w1.size
        self.nw2 = w2.size
        self.nw3 = w3.size

        self.n1 = w1.size-1
        self.n2 = w2.size-1
        self.n3 = w3.size-1

        self.w1 = w1
        self.w2 = w2
        self.w3 = w3

        cdef numpy.ndarray[double, ndim=3, mode="c"] volume = \
                numpy.zeros((self.n1, self.n2, self.n3), dtype=float)
        for i in range(self.n1):
            for j in range(self.n2):
                for k in range(self.n3):
                    self.volume[i,j,k] = (self.w1[i+1] - self.w1[i])* \
                        (self.w2[j+1] - self.w2[j])*(self.w3[k+1] - self.w3[k])

        self.volume = volume

        self.obj = new CartesianGrid(self.n1, self.n2, self.n3, self.nw1, \
                self.nw2, self.nw3, &w1[0], &w2[0], &w3[0], &volume[0,0,0])

        super(CartesianGridObj, self).__init__()

    def __del__(self):
        del self.obj

