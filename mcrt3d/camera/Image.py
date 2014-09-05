from ..mcrt3d import lib
import numpy

class Image:
    def __init__(self, r, incl, pa, x, y, nx, ny, nu, pixel_size, nnu):
        self.r = r
        self.incl = incl
        self.pa = pa
        self.x = x
        self.y = y
        self.intensity = numpy.zeros((nx, ny, nnu),dtype=float)
        self.nx = nx
        self.ny = ny
        self.nu = nu
        self.pixel_size = pixel_size
        self.nnu = nnu

        self.obj = lib.new_Image(self.r, self.incl, self.pa, self.x, \
                self.y, self.intensity, self.nx, self.ny, self.nu, \
                self.pixel_size, self.nnu)
