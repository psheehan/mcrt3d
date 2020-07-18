import cython

import numpy
cimport numpy

from ..mcrt3d cimport SphericalGrid
from .Grid cimport GridObj

cdef class SphericalGridObj(GridObj):
    def __init__(self, numpy.ndarray[double, ndim=1, mode="c"] w1, \
            numpy.ndarray[double, ndim=1, mode="c"] w2, \
            numpy.ndarray[double, ndim=1, mode="c"] w3):

        self.coordsystem = "spherical"

        self.r = 0.5*(w1[0:w1.size-1] + w1[1:w1.size])
        self.theta = 0.5*(w2[0:w2.size-1] + w2[1:w2.size])
        self.phi = 0.5*(w3[0:w3.size-1] + w3[1:w3.size])

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
                    volume[i,j,k] = (self.w1[i+1]**3 - self.w1[i]**3)* \
                        (self.w3[k+1] - self.w3[k])* \
                        (numpy.cos(self.w2[j]) - numpy.cos(self.w2[j+1]))/3

        self.volume = volume

        if abs(numpy.cos(self.w2.max())) < 1.0e-6:
            self.mirror_symmetry = True

            self.volume *= 2
        else:
            self.mirror_symmetry = False

        self.obj = new SphericalGrid(self.n1, self.n2, self.n3, self.nw1, \
                self.nw2, self.nw3, &w1[0], &w2[0], &w3[0], &volume[0,0,0], \
                self.mirror_symmetry)

        super(SphericalGridObj, self).__init__()

    def __del__(self):
        del self.obj

