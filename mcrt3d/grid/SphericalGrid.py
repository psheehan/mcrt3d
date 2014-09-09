from ..mcrt3d import lib
from .Grid import Grid
import numpy

class SphericalGrid:
    def __init__(self):
        self.obj = lib.new_SphericalGrid()
        self.coordsystem = "spherical"

        super(SphericalGrid, self).__init__()

    def set_walls(self, w1, w2, w3):
        self.r = 0.5*(w1[0:w1.size-1] + w1[1:w1.size])
        self.theta = 0.5*(w2[0:w2.size-1] + w2[1:w2.size])
        self.phi = 0.5*(w3[0:w3.size-1] + w3[1:w3.size])

        self.w1 = w1
        self.w2 = w2
        self.w3 = w3

        self.volume = numpy.zeros((w1.size-1, w2.size-1, w3.size-1), dtype=float)
        for i in range(self.volume.shape[0]):
            for j in range(self.volume.shape[1]):
                for k in range(self.volume.shape[2]):
                    self.volume[i,j,k] = (self.w1[i+1]**3 - self.w1[i]**3)* \
                        (self.w3[k+1] - self.w3[k])* \
                        (numpy.cos(self.w2[j]) - numpy.cos(self.w2[j+1]))/3

        lib.set_walls(self.obj, w1.size-1, w2.size-1, w3.size-1, \
                w1.size, w2.size, w3.size, w1, w2, w3, self.volume)
